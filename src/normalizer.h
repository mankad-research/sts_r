// normalizer.h -- per-topic log-partition function and its first two derivatives.
//
// For topic k the normalizer
//
//     logZ_k(a) = log sum_{v=1}^{V} exp( mv_v + kappa_t[v,k] + kappa_s[v,k] * a )
//
// depends on a document only through the scalar a = alpha_s[k].  Within one
// E-step mv, kappa_t and kappa_s are fixed, so logZ_k is a fixed 1-D function
// shared by every document.
//
// The identities that make this worth isolating (verified numerically to
// ~5e-16 against the shipped lgaecpp/esthcpp reductions):
//
//     logZ_k'(a)  = sum_v kappa_s[v,k] * beta[v,k]          = E_beta[kappa_s]
//     logZ_k''(a) = sum_v kappa_s[v,k]^2 * beta[v,k] - (logZ_k')^2
//                                                          = Var_beta[kappa_s]
//
// These are exactly the two full-V reductions that lgaecpp (kap_beta_sum_vec)
// and esthcpp (kap_beta_sum_rows) compute, and logZ_k itself is the only
// full-V quantity the objective needs.  Every other term in the objective,
// gradient and Hessian touches only the document's own rows.  Supplying the
// triple (logZ, logZ', logZ'') per topic therefore reduces the per-document
// E-step cost from O(V*K) to O(n_d*K), where n_d is the document's
// unique-word count.
//
// Two backends implement the same interface:
//   ExactNormalizer   -- one fused pass over V per (topic, a).  Reproduces 1.6
//                        up to floating-point reassociation.
//   HermiteNormalizer -- tabulated once per E-step, quintic Hermite
//                        interpolation (see hermite.h).  Approximate, with a
//                        controllable error bound and an exact fallback
//                        outside the tabulated range.

#ifndef STS_NORMALIZER_H
#define STS_NORMALIZER_H

#include <RcppArmadillo.h>

namespace sts {

// The triple returned for one topic at one value of a.
struct LogZTriple {
  double f;    // logZ_k(a)
  double d1;   // logZ_k'(a)
  double d2;   // logZ_k''(a)
};

class Normalizer {
public:
  virtual ~Normalizer() {}
  // Value and first two derivatives of logZ_k at a.
  virtual LogZTriple eval(unsigned int k, double a) const = 0;
  // Value only; backends may have a cheaper path.
  virtual double value(unsigned int k, double a) const { return eval(k, a).f; }
};

// ---------------------------------------------------------------- exact ----
//
// lin_base is the V x K matrix mv + kappa_t, precomputed once per E-step
// (1.6 rebuilt it on every objective and gradient call).  One pass over V
// accumulates all three quantities with a single exp() per element.
class ExactNormalizer : public Normalizer {
public:
  ExactNormalizer(const arma::mat& lin_base, const arma::mat& kappa_s)
    : lin_base_(lin_base), kappa_s_(kappa_s),
      V_(lin_base.n_rows), K_(lin_base.n_cols) {}

  LogZTriple eval(unsigned int k, double a) const {
    const double* lb = lin_base_.colptr(k);
    const double* ks = kappa_s_.colptr(k);

    // Stabilise exactly as 1.6's lpbdcpp/lgaecpp do: subtract the column max.
    double m = -std::numeric_limits<double>::infinity();
    for (arma::uword v = 0; v < V_; ++v) {
      const double l = lb[v] + ks[v] * a;
      if (l > m) m = l;
    }

    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (arma::uword v = 0; v < V_; ++v) {
      const double e = std::exp(lb[v] + ks[v] * a - m);
      s0 += e;
      s1 += ks[v] * e;
      s2 += ks[v] * ks[v] * e;
    }

    LogZTriple out;
    out.f  = m + std::log(s0);
    out.d1 = s1 / s0;
    out.d2 = s2 / s0 - out.d1 * out.d1;
    return out;
  }

  double value(unsigned int k, double a) const {
    const double* lb = lin_base_.colptr(k);
    const double* ks = kappa_s_.colptr(k);
    double m = -std::numeric_limits<double>::infinity();
    for (arma::uword v = 0; v < V_; ++v) {
      const double l = lb[v] + ks[v] * a;
      if (l > m) m = l;
    }
    double s0 = 0.0;
    for (arma::uword v = 0; v < V_; ++v) s0 += std::exp(lb[v] + ks[v] * a - m);
    return m + std::log(s0);
  }

  arma::uword V() const { return V_; }
  arma::uword K() const { return K_; }

private:
  const arma::mat& lin_base_;
  const arma::mat& kappa_s_;
  arma::uword V_, K_;
};

}  // namespace sts

#endif  // STS_NORMALIZER_H
