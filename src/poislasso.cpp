// poislasso.cpp -- dedicated Poisson-lasso solver for the kappa M-step.
//
// 1.6's opt.kappa() calls glmnet() once per vocabulary word: on poliblog2k
// (V = 2632, K = 10, numGroups = 2) that is 2.22 ms per word, or 5.8 s per
// M-step, for a design that is only 20 x 20 with two nonzeros per row. Almost
// all of it is glmnet's per-call R and S4 overhead rather than arithmetic.
//
// The design has exploitable structure. With
//     X = [x1_agg, alpha_agg],   rows indexed by (group g, topic k)
// row (g,k) holds a 1 in column k and abar[g,k] in column K+k, and zeros
// elsewhere. So the fit separates into K independent two-coefficient Poisson
// lasso problems with numGroups rows each, coupled only through the shared
// lambda path and the joint AIC. Verified on the real design: topic-separable,
// two rows and two coefficients per topic.
//
// What has to be reproduced exactly, since the AIC choice depends on all of
// it (all three determined empirically against glmnet 4.1-8):
//   lambda_max = max_j |sum_i x_ij (mu0_i - y_i)| / n,  mu0 = exp(offset),
//                n = numGroups*K        [matched to 2.4e-16 relative]
//   the grid    = lambda_max * r^t, r = lambda.min.ratio^(1/(nlambda-1)),
//                 built at full length regardless of truncation
//   truncation  = stop after the first lambda with dev.ratio > devmax (0.999)
//
// The per-lambda subproblems are solved to a tighter tolerance than glmnet's
// default thresh, by coordinate descent with a safeguarded 1-D Newton step per
// coordinate. The penalised Poisson objective is strictly convex in the
// identified directions, so its minimiser is unique and the algorithm used to
// reach it does not matter -- unlike the E-step, where the BFGS path had to be
// preserved, here only the optimum needs to agree.

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

namespace {

// Poisson deviance contribution, with the y log(y/mu) term understood as 0
// at y = 0 (its limit).
inline double dev_term(double y, double mu) {
  const double a = (y > 0.0) ? y * std::log(y / mu) : 0.0;
  return 2.0 * (a - (y - mu));
}

// One topic's two-coefficient subproblem at a single lambda.
//
//   minimise  (1/n) sum_g [ -y_g eta_g + exp(eta_g) ] + lambda (|b1| + |b2|)
//   with      eta_g = o_g + b1 + b2 a_g
//
// b1/b2 are updated in place from their warm-started values. The 1/n uses the
// GLOBAL n, which is what couples the topics to a common lambda scale.
// One topic's subproblem. E_g = exp(offset_g) is precomputed once per
// (word, topic) and reused for every lambda on the path, since
//     exp(o_g + b1 + b2 a_g) = E_g * e^{b1} * e^{b2 a_g}.
//
// The intercept coordinate then has a CLOSED FORM for any numGroups:
//     grad1(b1) = (1/n) ( e^{b1} S - Y ),   S = sum_g E_g e^{b2 a_g},
//                                           Y = sum_g y_g
// so stationarity grad1 + lambda*s1 = 0 gives e^{b1} = (Y - n lambda s1)/S
// directly -- no iteration. Only the slope coordinate needs a 1-D Newton
// solve, and it is warm-started along the lambda path. That is what makes
// this cheaper than a generic coordinate-descent lasso.
struct TopicFit {
  const double* y;    // numGroups
  const double* E;    // numGroups, exp(offset)
  const double* a;    // numGroups, abar[, k]
  int G;
  double n;           // global nobs (couples topics to a common lambda scale)
  double inv_n;
  double Y;           // sum_g y_g
  double Ya;          // sum_g a_g y_g
  double SE;          // sum_g E_g          (S at b2 = 0)
  double SaE;         // sum_g a_g E_g      (slope gradient at b2 = 0, per e^{b1})

  // Closed-form intercept update given b2. Returns the new b1.
  double update_b1(double b2, double lambda) const {
    double S = SE;
    if (b2 != 0.0) { S = 0.0; for (int g = 0; g < G; ++g) S += E[g] * std::exp(b2 * a[g]); }
    const double g0 = (S - Y) * inv_n;            // gradient at b1 = 0
    if (std::fabs(g0) <= lambda * (1.0 + 1e-10)) return 0.0;   // see note at lambda_max
    const double s1 = (g0 > 0.0) ? -1.0 : 1.0;
    const double rhs = Y - n * lambda * s1;
    if (!(rhs > 0.0)) return 0.0;                 // mean would be non-positive
    return std::log(rhs / S);
  }

  // Slope update given b1, by safeguarded Newton on the scalar stationarity
  // condition. grad2 is strictly increasing in b2 (the a_g^2 weights are
  // positive), so the root is unique once the sign is fixed.
  double update_b2(double b1, double b2, double lambda) const {
    const double c = std::exp(b1);
    // gradient at b2 = 0 needs no exp(): sum_g a_g E_g is precomputed
    const double g0 = (c * SaE - Ya) * inv_n;
    if (std::fabs(g0) <= lambda * (1.0 + 1e-10)) return 0.0;   // see note at lambda_max
    const double s2 = (g0 > 0.0) ? -1.0 : 1.0;
    double v = b2;
    if (v * s2 <= 0.0) v = 0.0;                   // warm start on wrong side

    for (int it = 0; it < 40; ++it) {
      // gradient and Hessian share one exp() per group
      double sg = 0.0, sh = 0.0;
      for (int g = 0; g < G; ++g) {
        const double e = E[g] * std::exp(v * a[g]);
        sg += a[g] * e; sh += a[g] * a[g] * e;
      }
      const double gv = (c * sg - Ya) * inv_n + lambda * s2;
      if (std::fabs(gv) < 1e-13) break;
      const double hv = c * sh * inv_n;
      if (!(hv > 0.0) || !std::isfinite(hv)) break;
      double step = -gv / hv;
      if (step >  2.0) step =  2.0;               // guard exp() overshoot
      if (step < -2.0) step = -2.0;
      double vn = v + step;
      if (vn * s2 < 0.0) vn = 0.0;
      if (std::fabs(vn - v) < 1e-14) { v = vn; break; }
      v = vn;
    }
    return v;
  }

  // Closed-form solve for numGroups == 2 (a binary initializationVar, the
  // common case). With two rows and two coefficients the stationarity
  // conditions are LINEAR in the fitted means mu_g = E_g exp(b1 + b2 a_g):
  //   both active, signs (s1,s2):  mu1 + mu2         = Y  - n lambda s1
  //                                a1 mu1 + a2 mu2   = Ya - n lambda s2
  // a 2x2 solve for mu, then two logs give (b1, b2). The intercept-only case
  // is one log; the slope-only case is the 1-D Newton in update_b2. The
  // problem is strictly convex when a1 != a2, so exactly one sign pattern
  // satisfies the KKT conditions; the warm-start pattern is tried first.
  // Returns false (caller falls back to solve()) if a1 == a2 or no pattern
  // passes, which can only happen through rounding at a boundary.
  bool solve_g2(double& c1, double& c2, double lambda, double* mu) const {
    const double a1 = a[0], a2 = a[1], da = a2 - a1;
    if (!(std::fabs(da) > 1e-12)) return false;
    const double tolf = 1.0 + 1e-10, nl = n * lambda;
    const double w1 = c1, w2 = c2;                     // warm start

    auto try_zero = [&]() -> bool {
      if (std::fabs((SE - Y) * inv_n) > lambda * tolf) return false;
      if (std::fabs((SaE - Ya) * inv_n) > lambda * tolf) return false;
      c1 = 0.0; c2 = 0.0; mu[0] = E[0]; mu[1] = E[1]; return true;
    };
    auto try_b1 = [&](double s1) -> bool {
      const double rhs = Y - nl * s1;
      if (!(rhs > 0.0)) return false;
      const double e = rhs / SE, b = std::log(e);
      if (!(b * s1 > 0.0)) return false;
      const double m1 = E[0] * e, m2 = E[1] * e;
      if (std::fabs((a1 * m1 + a2 * m2 - Ya) * inv_n) > lambda * tolf) return false;
      c1 = b; c2 = 0.0; mu[0] = m1; mu[1] = m2; return true;
    };
    auto try_both = [&](double s1, double s2) -> bool {
      const double R1 = Y - nl * s1, R2 = Ya - nl * s2;
      const double m1 = (a2 * R1 - R2) / da, m2 = (R2 - a1 * R1) / da;
      if (!(m1 > 0.0 && m2 > 0.0)) return false;
      const double L1 = std::log(m1 / E[0]), L2 = std::log(m2 / E[1]);
      const double b2 = (L2 - L1) / da, b1 = L1 - b2 * a1;
      if (!(b1 * s1 > 0.0 && b2 * s2 > 0.0)) return false;
      c1 = b1; c2 = b2; mu[0] = m1; mu[1] = m2; return true;
    };
    auto try_b2 = [&]() -> bool {
      const double v = update_b2(0.0, w2, lambda);
      if (v == 0.0) return false;
      const double m1 = E[0] * std::exp(v * a1), m2 = E[1] * std::exp(v * a2);
      if (std::fabs((m1 + m2 - Y) * inv_n) > lambda * tolf) return false;
      c1 = 0.0; c2 = v; mu[0] = m1; mu[1] = m2; return true;
    };

    const double sg1 = (w1 >= 0.0) ? 1.0 : -1.0, sg2 = (w2 >= 0.0) ? 1.0 : -1.0;
    // warm-start pattern first
    if (w1 == 0.0 && w2 == 0.0) { if (try_zero()) return true; }
    else if (w2 == 0.0)         { if (try_b1(sg1)) return true; }
    else if (w1 != 0.0)         { if (try_both(sg1, sg2)) return true; }
    else                        { if (try_b2()) return true; }
    // then every pattern
    if (try_zero()) return true;
    if (try_b1(1.0) || try_b1(-1.0)) return true;
    if (try_both(1.0, 1.0) || try_both(1.0, -1.0) ||
        try_both(-1.0, 1.0) || try_both(-1.0, -1.0)) return true;
    if (try_b2()) return true;
    c1 = w1; c2 = w2;                                   // restore for fallback
    return false;
  }

  // Full solve at one lambda, warm-started from (c1, c2).
  //
  // Pure coordinate descent is slow here: the intercept column (1,...,1) and
  // the slope column (a_1,...,a_G) are strongly correlated, so alternating
  // exact coordinate updates zig-zag and converge only linearly. Instead,
  // coordinate updates are used to settle the active set and signs, and once
  // both coefficients are active the stationarity system
  //     (1/n) X_k' (mu - y) + lambda * s = 0
  // is solved jointly by a 2x2 Newton iteration, which converges
  // quadratically. A sign flip during Newton means the active set guessed was
  // wrong, so it falls back to coordinate descent.
  void solve(double& c1, double& c2, double lambda, double tol,
             long& n_sweep, long& n_newton) const {
    // phase 1: coordinate sweeps until the active pattern is stable
    for (int sweep = 0; sweep < 100; ++sweep) {
      const double p1 = c1, p2 = c2;
      c1 = update_b1(c2, lambda);
      c2 = update_b2(c1, c2, lambda);
      ++n_sweep;
      if (std::fabs(c1 - p1) < tol && std::fabs(c2 - p2) < tol) return;
      const bool stable = ((p1 != 0.0) == (c1 != 0.0)) &&
                          ((p2 != 0.0) == (c2 != 0.0)) &&
                          p1 * c1 >= 0.0 && p2 * c2 >= 0.0;
      if (stable && sweep >= 1 && c1 != 0.0 && c2 != 0.0) break;
    }
    if (c1 == 0.0 || c2 == 0.0) return;   // single-active cases end above

    // phase 2: joint Newton with fixed signs
    const double s1 = (c1 > 0.0) ? 1.0 : -1.0;
    const double s2 = (c2 > 0.0) ? 1.0 : -1.0;
    for (int it = 0; it < 50; ++it) {
      ++n_newton;
      const double e1 = std::exp(c1);
      double S0 = 0.0, S1 = 0.0, S2 = 0.0;
      for (int g = 0; g < G; ++g) {
        const double mu = E[g] * e1 * std::exp(c2 * a[g]);
        S0 += mu; S1 += a[g] * mu; S2 += a[g] * a[g] * mu;
      }
      const double r1 = (S0 - Y)  * inv_n + lambda * s1;
      const double r2 = (S1 - Ya) * inv_n + lambda * s2;
      if (std::fabs(r1) < 1e-13 && std::fabs(r2) < 1e-13) return;
      const double h11 = S0 * inv_n, h12 = S1 * inv_n, h22 = S2 * inv_n;
      const double det = h11 * h22 - h12 * h12;
      if (!(det > 0.0) || !std::isfinite(det)) break;
      double d1 = -( h22 * r1 - h12 * r2) / det;
      double d2 = -(-h12 * r1 + h11 * r2) / det;
      const double mx = std::max(std::fabs(d1), std::fabs(d2));
      if (mx > 2.0) { d1 *= 2.0 / mx; d2 *= 2.0 / mx; }
      const double n1 = c1 + d1, n2 = c2 + d2;
      if (n1 * s1 <= 0.0 || n2 * s2 <= 0.0) {
        // active set was wrong: finish with coordinate descent from here
        for (int sweep = 0; sweep < 200; ++sweep) {
          const double p1 = c1, p2 = c2;
          c1 = update_b1(c2, lambda);
          c2 = update_b2(c1, c2, lambda);
          ++n_sweep;
          if (std::fabs(c1 - p1) < tol && std::fabs(c2 - p2) < tol) break;
        }
        return;
      }
      if (mx < 1e-14) { c1 = n1; c2 = n2; return; }
      c1 = n1; c2 = n2;
    }
  }
};

// Solve one topic at one lambda and return its deviance contribution. Uses the
// closed form when numGroups == 2, the iterative solver otherwise or as a
// fallback. Shared by the M-step and the path diagnostic so they agree.
inline double solve_point(const TopicFit& tf, double& c1, double& c2, double lambda,
                          double tol, long& n_sweep, long& n_newton) {
  if (tf.G == 2) {
    double mu[2];
    if (tf.solve_g2(c1, c2, lambda, mu))
      return dev_term(tf.y[0], mu[0]) + dev_term(tf.y[1], mu[1]);
  }
  tf.solve(c1, c2, lambda, tol, n_sweep, n_newton);
  const double ec1 = std::exp(c1);
  double d = 0.0;
  for (int g = 0; g < tf.G; ++g) d += dev_term(tf.y[g], tf.E[g] * ec1 * std::exp(c2 * tf.a[g]));
  return d;
}

}  // namespace

// Replacement for opt.kappa()'s per-word glmnet loop.
//
// Words are independent, so the loop runs in parallel over words with
// RcppParallel. Each word's fit depends only on its own data, every output is
// a disjoint row, and the diagnostic counters are kept per word and summed
// serially, so the result is identical for any number of threads.
namespace {

struct OptKappaWorker : public RcppParallel::Worker {
  const std::vector<arma::mat>& phi;      // numGroups x (V x K)
  const arma::mat& phi_sum;               // G x K
  const arma::mat& abar;                  // G x K
  const arma::vec& mv;
  const arma::mat& kappa_t_old;
  const arma::mat& kappa_s_old;
  const int nlambda;
  const double r, devmax, tol;
  const double damping;                   // weight on the OLD kappa (1.6: 0.5)
  const int select;                       // 0 = AIC (1.6), 1 = last path point,
                                          // 2 = solve at a given lambda per word
  const arma::vec& lambda_fixed;          // used when select == 2
  const int G;
  const arma::uword K;
  arma::mat& kt;
  arma::mat& ks;
  std::vector<double>& npath, &nsweep, &nnewton;
  arma::vec& lambda_sel;                  // lambda at the selected point, per word

  OptKappaWorker(const std::vector<arma::mat>& phi_, const arma::mat& phi_sum_,
                 const arma::mat& abar_, const arma::vec& mv_,
                 const arma::mat& kto, const arma::mat& kso, int nlambda_,
                 double r_, double devmax_, double tol_, double damping_, int select_,
                 const arma::vec& lambda_fixed_, int G_, arma::uword K_,
                 arma::mat& kt_, arma::mat& ks_, std::vector<double>& npath_,
                 std::vector<double>& nsweep_, std::vector<double>& nnewton_,
                 arma::vec& lambda_sel_)
    : phi(phi_), phi_sum(phi_sum_), abar(abar_), mv(mv_), kappa_t_old(kto),
      kappa_s_old(kso), nlambda(nlambda_), r(r_), devmax(devmax_), tol(tol_),
      damping(damping_), select(select_), lambda_fixed(lambda_fixed_), G(G_), K(K_),
      kt(kt_), ks(ks_), npath(npath_), nsweep(nsweep_), nnewton(nnewton_),
      lambda_sel(lambda_sel_) {}

  void operator()(std::size_t v0, std::size_t v1) {
    const int n = G * (int)K;
    const double inv_n = 1.0 / (double)n;
    // per-thread scratch, laid out topic-major (index k*G + g)
    std::vector<double> yv(K * G), Ev(K * G), av(K * G), Yk(K), Yak(K), SEk(K), SaEk(K);
    std::vector<double> b1(K), b2(K), bb1(K), bb2(K);

    for (std::size_t v = v0; v < v1; ++v) {
      long n_sweep = 0, n_newton = 0, np = 0;

      // ---- lambda_max and nulldev over the whole stacked problem ----
      double gmax = 0.0, nulldev = 0.0;
      for (arma::uword k = 0; k < K; ++k) {
        double g_int = 0.0, g_slp = 0.0;
        for (int g = 0; g < G; ++g) {
          const double yy = phi[g](v, k);
          const double mu0 = std::exp(phi_sum(g, k) + mv[v]);
          g_int += (mu0 - yy);
          g_slp += abar(g, k) * (mu0 - yy);
          nulldev += dev_term(yy, mu0);
        }
        gmax = std::max(gmax, std::fabs(g_int));
        gmax = std::max(gmax, std::fabs(g_slp));
      }
      const double lambda_max = gmax * inv_n;

      // ---- per-(word, topic) data, computed once for the whole path ----
      for (arma::uword k = 0; k < K; ++k) {
        double Y = 0.0, Ya = 0.0, SE = 0.0, SaE = 0.0;
        for (int g = 0; g < G; ++g) {
          const arma::uword i = k * G + g;
          yv[i] = phi[g](v, k);
          Ev[i] = std::exp(phi_sum(g, k) + mv[v]);
          av[i] = abar(g, k);
          Y += yv[i]; Ya += av[i] * yv[i];
          SE += Ev[i]; SaE += av[i] * Ev[i];
        }
        Yk[k] = Y; Yak[k] = Ya; SEk[k] = SE; SaEk[k] = SaE;
      }

      // ---- walk the path, tracking AIC ----
      for (arma::uword k = 0; k < K; ++k) { b1[k] = b2[k] = bb1[k] = bb2[k] = 0.0; }
      double best_ic = std::numeric_limits<double>::infinity();

      double lambda = lambda_max;
      double sel_lambda = lambda_max;
      // select == 2: walk the usual grid down to the fixed lambda (for warm
      // starts), then solve exactly there and take that point.
      const double lam_fix = (select == 2) ? lambda_fixed[v] : 0.0;
      bool fixed_done = false;
      for (int t = 0; t < nlambda; ++t, lambda *= r) {
        if (select == 2 && (lambda <= lam_fix || t == nlambda - 1)) {
          lambda = std::max(lam_fix, 0.0); fixed_done = true;
        }
        double dev = 0.0;
        int df = 0;
        // lambda_max is, by definition, the smallest lambda at which every
        // coefficient is zero, and glmnet sets the first path solution to
        // zero explicitly. Do the same rather than re-deriving it: the
        // gradient at zero inside the solver is accumulated differently from
        // lambda_max and the two can differ by an ulp, which on some words
        // let a ~1e-15 coefficient "enter" at the first point and moved the
        // AIC choice to the runner-up.
        if (t == 0 && !(select == 2 && fixed_done && lam_fix < lambda_max)) { dev = nulldev; df = 0; }
        else for (arma::uword k = 0; k < K; ++k) {
          const double* y = &yv[k * G];
          const double* E = &Ev[k * G];
          const double* a = &av[k * G];
          TopicFit tf{ y, E, a, G, (double)n, inv_n, Yk[k], Yak[k], SEk[k], SaEk[k] };
          double c1 = b1[k], c2 = b2[k];                 // warm start along the path
          dev += solve_point(tf, c1, c2, lambda, tol, n_sweep, n_newton);
          b1[k] = c1; b2[k] = c2;
          if (c1 != 0.0) ++df;
          if (c2 != 0.0) ++df;
        }

        const double ic = dev + 2.0 * (double)df;
        // select == 1 keeps the least-penalised point reached (diagnostic:
        // removes the AIC-driven active-set changes between iterations)
        if (select == 1 || (select == 2 && fixed_done) || (select == 0 && ic < best_ic)) {
          best_ic = ic; sel_lambda = lambda;
          for (arma::uword k = 0; k < K; ++k) { bb1[k] = b1[k]; bb2[k] = b2[k]; }
        }
        ++np;
        if (select == 2) { if (fixed_done) break; else continue; }

        // Exact early stop for AIC selection. Along a lasso path the l1 norm of
        // the solution is non-decreasing as lambda falls, and the problem
        // separates by topic, so a topic with any nonzero coefficient keeps at
        // least one at every later point. Deviance is non-negative. Every
        // later point therefore has AIC >= 2 * (number of active topics now),
        // and once that reaches the best AIC so far no later point can be
        // selected (selection takes the first strict minimum). Saves ~17% of
        // path points on poliblog2k with identical selections.
        if (select == 0) {
          int active = 0;
          for (arma::uword k = 0; k < K; ++k) if (b1[k] != 0.0 || b2[k] != 0.0) ++active;
          if (2.0 * (double)active >= best_ic) break;
        }
        // glmnet's truncation: stop after the first lambda whose dev.ratio
        // exceeds devmax. Guard nulldev == 0 (a word absent from every group).
        const double devratio = (nulldev > 0.0) ? 1.0 - dev / nulldev : 1.0;
        if (devratio > devmax) break;
      }

      // ---- damped update; damping = 0.5 reproduces 1.6's (coef + kappa)/2 ----
      for (arma::uword k = 0; k < K; ++k) {
        if (damping == 0.5) {                 // keep 1.6's exact arithmetic
          kt(v, k) = 0.5 * (bb1[k] + kappa_t_old(v, k));
          ks(v, k) = 0.5 * (bb2[k] + kappa_s_old(v, k));
        } else {
          kt(v, k) = (1.0 - damping) * bb1[k] + damping * kappa_t_old(v, k);
          ks(v, k) = (1.0 - damping) * bb2[k] + damping * kappa_s_old(v, k);
        }
      }
      npath[v] = (double)np; nsweep[v] = (double)n_sweep; nnewton[v] = (double)n_newton;
      lambda_sel[v] = sel_lambda;
    }
  }
};

}  // namespace

// phi_list : numGroups matrices, each V x K
// abar     : numGroups x K, the within-group column means of alpha_s
// Returns the DAMPED kappa, matching 1.6's kappa <- (coef + kappa)/2.
// [[Rcpp::export]]
Rcpp::List sts_optkappa_cpp(const Rcpp::List& phi_list,
                            const arma::mat& abar,
                            const arma::vec& mv,
                            const arma::mat& kappa_t_old,
                            const arma::mat& kappa_s_old,
                            int nlambda,
                            double lambda_min_ratio,
                            double devmax,
                            double tol,
                            int cores = 1,
                            double damping = 0.5,
                            int select = 0,
                            Rcpp::Nullable<Rcpp::NumericVector> lambda_fixed = R_NilValue) {
  const int G = phi_list.size();
  std::vector<arma::mat> phi(G);
  for (int g = 0; g < G; ++g) phi[g] = Rcpp::as<arma::mat>(phi_list[g]);
  const arma::uword V = phi[0].n_rows;
  const arma::uword K = phi[0].n_cols;
  if (cores < 1) cores = 1;

  // offset's word-independent part: log(sum_v phi[[g]][v,k]), safelog'd at
  // -1000 exactly as safelog() in R/helperFunctions.R does.
  arma::mat phi_sum(G, K);
  for (int g = 0; g < G; ++g)
    for (arma::uword k = 0; k < K; ++k) {
      const double l = std::log(arma::accu(phi[g].col(k)));
      phi_sum(g, k) = (l < -1000.0 || !std::isfinite(l)) ? -1000.0 : l;
    }
  const double r = std::pow(lambda_min_ratio, 1.0 / (double)(nlambda - 1));

  arma::mat kt(V, K), ks(V, K);
  std::vector<double> npath(V), nsweep(V), nnewton(V);
  arma::vec lam_fix(V, arma::fill::zeros), lam_sel(V, arma::fill::zeros);
  if (select == 2) {
    if (lambda_fixed.isNull()) Rcpp::stop("select = 2 needs lambda_fixed");
    lam_fix = Rcpp::as<arma::vec>(lambda_fixed.get());
    if (lam_fix.n_elem != V) Rcpp::stop("lambda_fixed must have length V");
  }
  OptKappaWorker w(phi, phi_sum, abar, mv, kappa_t_old, kappa_s_old, nlambda,
                   r, devmax, tol, damping, select, lam_fix, G, K, kt, ks,
                   npath, nsweep, nnewton, lam_sel);
  RcppParallel::parallelFor(0, V, w, 16, cores);

  double tp = 0.0, ts = 0.0, tn = 0.0;
  for (arma::uword v = 0; v < V; ++v) { tp += npath[v]; ts += nsweep[v]; tn += nnewton[v]; }
  return Rcpp::List::create(
    Rcpp::Named("kappa_t") = kt,
    Rcpp::Named("kappa_s") = ks,
    Rcpp::Named("lambda") = lam_sel,
    Rcpp::Named("meanPath") = tp / (double)V,
    Rcpp::Named("sweepsPerSolve") = ts / (tp * (double)K),
    Rcpp::Named("newtonPerSolve") = tn / (tp * (double)K));
}

// Single-word diagnostic: the full path the solver walks for one word, in the
// same stacked (group, topic) layout glmnet uses. For tests and for checking
// individual disagreements against glmnet.
//   y, offset : length G*K, ordered (g=1,k=1..K), (g=2,k=1..K), ...
//   abar      : G x K
// [[Rcpp::export]]
Rcpp::List sts_poislasso_path_cpp(const arma::vec& y, const arma::vec& offset,
                                  const arma::mat& abar, int nlambda,
                                  double lambda_min_ratio, double devmax,
                                  double tol) {
  const int G = abar.n_rows;
  const arma::uword K = abar.n_cols;
  const int n = G * (int)K;
  const double inv_n = 1.0 / (double)n;
  std::vector<double> yv(K * G), Ev(K * G), av(K * G), Yk(K), Yak(K), SEk(K), SaEk(K);
  double gmax = 0.0, nulldev = 0.0;
  for (arma::uword k = 0; k < K; ++k) {
    double Y = 0.0, Ya = 0.0, gi = 0.0, gs = 0.0, SE = 0.0, SaE = 0.0;
    for (int g = 0; g < G; ++g) {
      const arma::uword src = (arma::uword)g * K + k, i = k * G + g;
      yv[i] = y[src]; Ev[i] = std::exp(offset[src]); av[i] = abar(g, k);
      Y += yv[i]; Ya += av[i] * yv[i]; SE += Ev[i]; SaE += av[i] * Ev[i];
      gi += Ev[i] - yv[i]; gs += av[i] * (Ev[i] - yv[i]);
      nulldev += dev_term(yv[i], Ev[i]);
    }
    Yk[k] = Y; Yak[k] = Ya; SEk[k] = SE; SaEk[k] = SaE;
    gmax = std::max(gmax, std::max(std::fabs(gi), std::fabs(gs)));
  }
  const double lambda_max = gmax * inv_n;
  const double r = std::pow(lambda_min_ratio, 1.0 / (double)(nlambda - 1));

  std::vector<double> b1(K, 0.0), b2(K, 0.0);
  std::vector<double> lam, dev, dfv, obj;
  arma::mat beta(2 * K, nlambda, arma::fill::zeros);
  long ns = 0, nn = 0;
  double lambda = lambda_max;
  int t = 0;
  for (; t < nlambda; ++t, lambda *= r) {
    double d = 0.0, pen = 0.0; int df = 0;
    if (t == 0) d = nulldev;               // all-zero by construction, as above
    else for (arma::uword k = 0; k < K; ++k) {
      TopicFit tf{ &yv[k * G], &Ev[k * G], &av[k * G], G, (double)n, inv_n, Yk[k], Yak[k], SEk[k], SaEk[k] };
      double c1 = b1[k], c2 = b2[k];
      d += solve_point(tf, c1, c2, lambda, tol, ns, nn);
      b1[k] = c1; b2[k] = c2;
      if (c1 != 0.0) ++df;
      if (c2 != 0.0) ++df;
      pen += std::fabs(c1) + std::fabs(c2);
      beta(k, t) = c1; beta(K + k, t) = c2;
    }
    lam.push_back(lambda); dev.push_back(d); dfv.push_back(df);
    obj.push_back(0.5 * d * inv_n + lambda * pen);   // glmnet's penalised objective
    const double dr = (nulldev > 0.0) ? 1.0 - d / nulldev : 1.0;
    if (dr > devmax) { ++t; break; }
  }
  return Rcpp::List::create(
    Rcpp::Named("lambda") = lam, Rcpp::Named("dev") = dev,
    Rcpp::Named("df") = dfv, Rcpp::Named("obj") = obj,
    Rcpp::Named("nulldev") = nulldev,
    Rcpp::Named("beta") = beta.cols(0, std::min<int>(t, nlambda) - 1));
}
