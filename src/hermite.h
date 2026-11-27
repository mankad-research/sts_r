// hermite.h -- tabulated log-partition backend (quintic Hermite).
//
// Within one E-step mv, kappa_t and kappa_s are fixed, so each topic's
// normalizer logZ_k is a fixed 1-D function of the scalar alpha_s[k], shared by
// every document (see normalizer.h).  Rather than sweeping all V words on every
// objective, gradient and Hessian evaluation, logZ_k is tabulated once per
// E-step on a uniform grid together with its first two derivatives, and
// interpolated thereafter.
//
// Why quintic rather than cubic.  The gradient needs logZ' and the Hessian
// needs logZ''.  A cubic Hermite fit to (logZ, logZ') at the nodes reproduces
// logZ to O(h^4) but its second derivative is only O(h^2) accurate, which would
// make the Hessian -- and hence sigma.ss and the ELBO -- the weakest link.
// Tabulating logZ'' as well and matching all three at both ends of each
// interval gives a quintic with
//
//     |logZ   - p  | = O(h^6)
//     |logZ'  - p' | = O(h^5)
//     |logZ'' - p''| = O(h^4)
//
// so every quantity the E-step consumes is nodally exact and the Hessian is no
// worse than O(h^4).  The three tabulated quantities come from one pass over V
// per node (the same accumulation ExactNormalizer does), so the table costs
// n_nodes * V * K work and then replaces n_eval * V * K, where n_eval is the
// number of distinct alpha evaluations in the E-step.  The win is therefore
// n_eval / n_nodes, which grows with D: the regime the method is for.
//
// The derivatives of logZ are central moments of kappa_s under beta, so they
// are bounded by powers of max|kappa_s| -- which is what makes the error
// bound above usable in practice rather than merely asymptotic.
//
// Chebyshev interpolation was considered and rejected: the exponent range here
// is wide (alpha_s spans roughly [-13, 6] on the profiling corpus), so a global
// polynomial converges too slowly. A uniform piecewise grid is the right tool.
//
// Outside the tabulated range the exact backend is used, so the result is never
// an extrapolation; out-of-range hits are counted for diagnostics.

#ifndef STS_HERMITE_H
#define STS_HERMITE_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "fpu.h"
#include <atomic>
#include <vector>
#include "normalizer.h"

namespace sts {

class HermiteNormalizer : public Normalizer {
public:
  // a_lo/a_hi bracket the alpha_s values the E-step is expected to visit;
  // spacing is the largest grid step h.
  //
  // Spacing per topic (sts 1.7.1).  The interpolation error grows with
  // S_k = max_v |kappa_s[v,k]|: the j-th derivative of logZ_k is a cumulant of
  // kappa_s under beta, bounded by a multiple of S_k^j, so the error in the
  // value is O((h S_k)^6) and in the j-th derivative O(S_k^j (h S_k)^(6-j)).
  // At a fixed h it is negligible for S_k of order 1-10 and large once one
  // word's kappa_s reaches the hundreds, which happens on small corpora run
  // for many iterations.  Each topic therefore gets its own step
  //     h_k = min(spacing, max_hS / S_k),
  // so h_k * S_k <= max_hS and the relative accuracy is the same for every
  // topic.  max_hS <= 0 disables this (h_k = spacing, as in 1.7.0).
  //
  // A topic whose table would need more than max_nodes nodes is evaluated
  // exactly instead (max_nodes == 0: no limit).  In logZ = "auto" mode the
  // E-step sets max_nodes from the expected number of evaluations, so a
  // topic is tabulated only when that is cheaper than summing over V.
  HermiteNormalizer(const arma::mat& lin_base, const arma::mat& kappa_s,
                    double a_lo, double a_hi, double spacing, int cores = 1,
                    double max_hS = 0.0, double max_nodes = 0.0)
    : exact_(lin_base, kappa_s), K_(lin_base.n_cols),
      a_lo_(a_lo), n_out_(0), n_exact_topic_(0) {
    h_.set_size(K_); n_.set_size(K_); off_.set_size(K_ + 1); a_hi_.set_size(K_);
    use_exact_.assign(K_, 0);
    off_[0] = 0;
    for (arma::uword k = 0; k < K_; ++k) {
      const double S = arma::max(arma::abs(kappa_s.col(k)));
      double hk = spacing;
      if (max_hS > 0.0 && S * spacing > max_hS) hk = max_hS / S;
      const double nk = std::ceil((a_hi - a_lo) / hk) + 1.0;
      h_[k] = hk;
      if (max_nodes > 0.0 && nk > max_nodes) {
        use_exact_[k] = 1;
        n_[k] = 0;
        a_hi_[k] = a_lo_;
      } else {
        n_[k] = std::max<arma::uword>(2, static_cast<arma::uword>(nk));
        a_hi_[k] = a_lo_ + (n_[k] - 1) * hk;
      }
      off_[k + 1] = off_[k] + n_[k];
    }
    f_.set_size(off_[K_]); d1_.set_size(off_[K_]); d2_.set_size(off_[K_]);

    // Every node is an independent exact evaluation, so the table is built
    // in parallel over nodes; each entry is computed identically whatever the
    // thread count, so the table is too.
    TableWorker w(*this);
    RcppParallel::parallelFor(0, off_[K_], w, 8, cores < 1 ? 1 : cores);
  }

  LogZTriple eval(unsigned int k, double a) const {
    if (use_exact_[k]) {                    // topic not tabulated
      n_exact_topic_.fetch_add(1, std::memory_order_relaxed);
      return exact_.eval(k, a);
    }
    if (a < a_lo_ || a > a_hi_[k]) {        // never extrapolate
      n_out_.fetch_add(1, std::memory_order_relaxed);
      return exact_.eval(k, a);
    }

    const double h = h_[k];
    const arma::uword nk = n_[k], o = off_[k];
    double pos = (a - a_lo_) / h;
    arma::uword j = static_cast<arma::uword>(pos);
    if (j >= nk - 1) j = nk - 2;               // right edge lands in last cell
    const double t = pos - static_cast<double>(j);

    const double f0 = f_[o + j],     f1 = f_[o + j + 1];
    const double g0 = d1_[o + j],    g1 = d1_[o + j + 1];
    const double s0 = d2_[o + j],    s1 = d2_[o + j + 1];

    // Quintic Hermite basis on t in [0,1]; see the header comment for the
    // interpolation conditions each one satisfies.
    const double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;

    const double H0 =  1.0 - 10.0 * t3 + 15.0 * t4 -  6.0 * t5;
    const double H1 =  t         -  6.0 * t3 +  8.0 * t4 -  3.0 * t5;
    const double H2 =  0.5 * (t2 -  3.0 * t3 +  3.0 * t4 -        t5);
    const double H3 =         10.0 * t3 - 15.0 * t4 +  6.0 * t5;
    const double H4 =         -4.0 * t3 +  7.0 * t4 -  3.0 * t5;
    const double H5 =  0.5 * (       t3 -  2.0 * t4 +        t5);

    const double dH0 = -30.0 * t2 +  60.0 * t3 - 30.0 * t4;
    const double dH1 = 1.0 - 18.0 * t2 + 32.0 * t3 - 15.0 * t4;
    const double dH2 = 0.5 * (2.0 * t -  9.0 * t2 + 12.0 * t3 -  5.0 * t4);
    const double dH3 =  30.0 * t2 -  60.0 * t3 + 30.0 * t4;
    const double dH4 = -12.0 * t2 +  28.0 * t3 - 15.0 * t4;
    const double dH5 = 0.5 * (3.0 * t2 -  8.0 * t3 +  5.0 * t4);

    const double ddH0 = -60.0 * t + 180.0 * t2 - 120.0 * t3;
    const double ddH1 = -36.0 * t +  96.0 * t2 -  60.0 * t3;
    const double ddH2 = 0.5 * (2.0 - 18.0 * t + 36.0 * t2 - 20.0 * t3);
    const double ddH3 =  60.0 * t - 180.0 * t2 + 120.0 * t3;
    const double ddH4 = -24.0 * t +  84.0 * t2 -  60.0 * t3;
    const double ddH5 = 0.5 * (6.0 * t - 24.0 * t2 + 20.0 * t3);

    const double hg0 = h * g0,     hg1 = h * g1;
    const double hs0 = h * h * s0, hs1 = h * h * s1;

    LogZTriple out;
    out.f  =  f0 * H0  + hg0 * H1  + hs0 * H2
            + f1 * H3  + hg1 * H4  + hs1 * H5;
    out.d1 = (f0 * dH0 + hg0 * dH1 + hs0 * dH2
            + f1 * dH3 + hg1 * dH4 + hs1 * dH5) / h;
    out.d2 = (f0 * ddH0 + hg0 * ddH1 + hs0 * ddH2
            + f1 * ddH3 + hg1 * ddH4 + hs1 * ddH5) / (h * h);
    return out;
  }

  // Total nodes over the tabulated topics.
  arma::uword n_nodes() const { return off_[K_]; }
  double range_lo() const { return a_lo_; }
  // Upper end of the table: the largest over tabulated topics (a_lo if none).
  double range_hi() const {
    double m = a_lo_;
    for (arma::uword k = 0; k < K_; ++k) if (!use_exact_[k] && a_hi_[k] > m) m = a_hi_[k];
    return m;
  }
  long n_outside() const { return n_out_.load(std::memory_order_relaxed); }
  long n_exact_topic_evals() const { return n_exact_topic_.load(std::memory_order_relaxed); }
  // Per-topic step (NA for topics evaluated exactly).
  arma::vec spacing_by_topic() const {
    arma::vec s = h_;
    for (arma::uword k = 0; k < K_; ++k) if (use_exact_[k]) s[k] = NA_REAL;
    return s;
  }
  int n_exact_topics() const {
    int c = 0; for (arma::uword k = 0; k < K_; ++k) c += use_exact_[k]; return c;
  }

private:
  struct TableWorker : public RcppParallel::Worker {
    FpuState fpu_;   // x87 state of the constructing (main) thread; see fpu.h
    HermiteNormalizer& hn;
    explicit TableWorker(HermiteNormalizer& hn_) : hn(hn_) {}
    void operator()(std::size_t i0, std::size_t i1) {
      FpuGuard fpu_guard(fpu_);
      // node i belongs to the topic k with off_[k] <= i < off_[k+1]
      arma::uword k = 0;
      while (k + 1 < hn.K_ && hn.off_[k + 1] <= i0) ++k;
      for (std::size_t i = i0; i < i1; ++i) {
        while (hn.off_[k + 1] <= i) ++k;
        const arma::uword j = i - hn.off_[k];
        const LogZTriple t = hn.exact_.eval(k, hn.a_lo_ + j * hn.h_[k]);
        hn.f_[i] = t.f; hn.d1_[i] = t.d1; hn.d2_[i] = t.d2;
      }
    }
  };

  ExactNormalizer exact_;
  arma::uword K_;
  double a_lo_;
  arma::vec h_, a_hi_;
  arma::uvec n_, off_;
  std::vector<int> use_exact_;
  arma::vec f_, d1_, d2_;
  mutable std::atomic<long> n_out_;
  mutable std::atomic<long> n_exact_topic_;
};

}  // namespace sts

#endif  // STS_HERMITE_H
