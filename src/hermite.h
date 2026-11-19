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
#include <atomic>
#include <vector>
#include "normalizer.h"

namespace sts {

class HermiteNormalizer : public Normalizer {
public:
  // a_lo/a_hi bracket the alpha_s values the E-step is expected to visit;
  // spacing is the grid step h.
  HermiteNormalizer(const arma::mat& lin_base, const arma::mat& kappa_s,
                    double a_lo, double a_hi, double spacing)
    : exact_(lin_base, kappa_s), K_(lin_base.n_cols),
      a_lo_(a_lo), h_(spacing), n_out_(0) {
    // Node count chosen so that the last node is at or beyond a_hi.
    n_nodes_ = static_cast<arma::uword>(std::ceil((a_hi - a_lo) / spacing)) + 1;
    if (n_nodes_ < 2) n_nodes_ = 2;
    a_hi_ = a_lo_ + (n_nodes_ - 1) * h_;

    f_.set_size(n_nodes_, K_);
    d1_.set_size(n_nodes_, K_);
    d2_.set_size(n_nodes_, K_);

    for (arma::uword k = 0; k < K_; ++k) {
      for (arma::uword j = 0; j < n_nodes_; ++j) {
        const LogZTriple t = exact_.eval(k, a_lo_ + j * h_);
        f_(j, k) = t.f; d1_(j, k) = t.d1; d2_(j, k) = t.d2;
      }
    }
  }

  LogZTriple eval(unsigned int k, double a) const {
    if (a < a_lo_ || a > a_hi_) {           // never extrapolate
      n_out_.fetch_add(1, std::memory_order_relaxed);
      return exact_.eval(k, a);
    }

    double pos = (a - a_lo_) / h_;
    arma::uword j = static_cast<arma::uword>(pos);
    if (j >= n_nodes_ - 1) j = n_nodes_ - 2;   // right edge lands in last cell
    const double t = pos - static_cast<double>(j);

    const double f0 = f_(j, k),     f1 = f_(j + 1, k);
    const double g0 = d1_(j, k),    g1 = d1_(j + 1, k);
    const double s0 = d2_(j, k),    s1 = d2_(j + 1, k);

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

    const double hg0 = h_ * g0,      hg1 = h_ * g1;
    const double hs0 = h_ * h_ * s0, hs1 = h_ * h_ * s1;

    LogZTriple out;
    out.f  =  f0 * H0  + hg0 * H1  + hs0 * H2
            + f1 * H3  + hg1 * H4  + hs1 * H5;
    out.d1 = (f0 * dH0 + hg0 * dH1 + hs0 * dH2
            + f1 * dH3 + hg1 * dH4 + hs1 * dH5) / h_;
    out.d2 = (f0 * ddH0 + hg0 * ddH1 + hs0 * ddH2
            + f1 * ddH3 + hg1 * ddH4 + hs1 * ddH5) / (h_ * h_);
    return out;
  }

  arma::uword n_nodes() const { return n_nodes_; }
  double range_lo() const { return a_lo_; }
  double range_hi() const { return a_hi_; }
  long n_outside() const { return n_out_.load(std::memory_order_relaxed); }

private:
  ExactNormalizer exact_;
  arma::uword K_, n_nodes_;
  double a_lo_, a_hi_, h_;
  arma::mat f_, d1_, d2_;
  mutable std::atomic<long> n_out_;
};

}  // namespace sts

#endif  // STS_HERMITE_H
