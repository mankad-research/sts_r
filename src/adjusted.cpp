// adjusted.cpp -- the kappa M-step of Chen and Mankad (2025, Management
// Science 71(7)), Section 4.2: Poisson regression with sample aggregation and
// the Proposition 2 adjustment (kappaEstimation = "adjusted").
//
// For word v and topic k the documents are aggregated within the groups g
// defined by initializationVar. The aggregated count
//     y_g = sum_{d in g} phi_{d,v,k}
// is Poisson with mean
//     sum_{d in g} phibar_{d,k} exp(m_v + kt + ks * alphahat_{d,k}),
// which equals (sum_{d in g} phibar_{d,k}) exp(m_v + kt + ks * abar_g) when the
// group's "equivalent" sentiment-discourse value is (equation 13)
//     abar_g = sum_{d in g} u_d alphahat_{d,k},
//     u_d    = phibar_{d,k} exp(ks alphahat_{d,k}) / sum_{d in g} phibar_{d,k} exp(ks alphahat_{d,k}).
// Since abar_g depends on ks, the paper starts from ks = 0 (abar_g is then the
// phibar-weighted mean of alphahat in the group), fits (kt, ks) by Poisson
// regression, recomputes abar_g from the new ks, and repeats until abar
// converges. That is what this file does, one (v, k) pair at a time: the
// aggregated design is separable by topic (row (g,k) has a 1 in column k and
// abar_g in column K + k) and the weights u_d for topic k involve only ks for
// that topic, so every (v, k) is an independent two-coefficient problem.
//
// The paper's regression is unpenalised. The maximum-likelihood estimate does
// not exist when a word never occurs in some group (y_g = 0 drives the fit to
// -Inf) and is not identified when the abar_g coincide, so a ridge term
// (ridge/2)(kt^2 + ks^2) is added to the log-likelihood as a numerical
// safeguard; at the default ridge = 1e-6 it changes well-determined estimates
// by a relative amount of order ridge / (expected count). Chen and Mankad's
// replication code made the same choice for the same reason (glmnet with
// alpha = 0 and a 1e-6 lambda ratio).
//
// The Proposition 2 iteration is run to convergence of abar (or maxInner
// passes; see the fallback below), and the result does not depend on the
// number of threads.

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

namespace {

// One group's documents for one topic: their alphahat and log phibar.
struct GroupDocs {
  std::vector<double> a;      // alphahat_{d,k}
  std::vector<double> lw;     // log phibar_{d,k}; documents with phibar = 0 dropped
};

// abar_g(ks): the u-weighted mean of alphahat (equation 13), computed with the
// usual max-shift so that large |ks| cannot overflow.
inline double weighted_abar(const GroupDocs& gd, double ks) {
  const std::size_t n = gd.a.size();
  if (n == 0) return 0.0;
  double mx = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < n; ++i) mx = std::max(mx, gd.lw[i] + ks * gd.a[i]);
  double s0 = 0.0, s1 = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double w = std::exp(gd.lw[i] + ks * gd.a[i] - mx);
    s0 += w; s1 += w * gd.a[i];
  }
  return s1 / s0;
}

// Ridge-safeguarded Poisson regression with two coefficients and G rows:
// maximise sum_g [ y_g eta_g - exp(eta_g) ] - ridge/2 (kt^2 + ks^2),
// eta_g = off_g + kt + ks x_g, by Newton's method with backtracking. The
// objective is strictly concave, so the maximiser is unique. (kt, ks) carry
// the warm start in and the solution out; returns the Newton iterations used,
// or -1 if the iteration limit was reached.
int fit_two_coef(const std::vector<double>& y, const std::vector<double>& off,
                 const std::vector<double>& x, double ridge,
                 double& kt, double& ks, int maxit = 200, double tol = 1e-12) {
  const std::size_t G = y.size();
  auto objective = [&](double a, double b) {
    double f = -0.5 * ridge * (a * a + b * b);
    for (std::size_t g = 0; g < G; ++g) {
      if (!std::isfinite(off[g])) continue;   // topic absent from this group
      const double eta = off[g] + a + b * x[g];
      if (eta > 700.0) return -std::numeric_limits<double>::infinity();
      f += y[g] * eta - std::exp(eta);
    }
    return f;
  };
  double f = objective(kt, ks);
  if (!std::isfinite(f)) {               // a warm start that overflows: restart
    kt = 0.0; ks = 0.0; f = objective(kt, ks);
  }
  for (int it = 1; it <= maxit; ++it) {
    double ga = -ridge * kt, gb = -ridge * ks;
    double h11 = ridge, h12 = 0.0, h22 = ridge;   // minus the Hessian
    for (std::size_t g = 0; g < G; ++g) {
      if (!std::isfinite(off[g])) continue;
      const double mu = std::exp(off[g] + kt + ks * x[g]);
      const double r  = y[g] - mu;
      ga += r; gb += r * x[g];
      h11 += mu; h12 += mu * x[g]; h22 += mu * x[g] * x[g];
    }
    const double det = h11 * h22 - h12 * h12;
    double da, db;
    if (det > 0.0 && std::isfinite(det)) {
      da = ( h22 * ga - h12 * gb) / det;
      db = (-h12 * ga + h11 * gb) / det;
    } else {                              // numerically singular: gradient step
      da = ga / h11; db = gb / h22;
    }
    // backtracking on the (concave) objective
    double step = 1.0, fn = objective(kt + da, ks + db);
    int halvings = 0;
    while (!(fn >= f) && halvings < 60) { step *= 0.5; fn = objective(kt + step * da, ks + step * db); ++halvings; }
    if (!(fn >= f)) return it;            // no ascent possible: at the optimum to rounding
    kt += step * da; ks += step * db;
    const double change = std::max(std::fabs(step * da) / (1.0 + std::fabs(kt)),
                                   std::fabs(step * db) / (1.0 + std::fabs(ks)));
    f = fn;
    if (change < tol) return it;
  }
  return -1;
}

struct AdjustedWorker : public RcppParallel::Worker {
  const std::vector<arma::mat>& phi;                    // G matrices, V x K
  const std::vector< std::vector<GroupDocs> >& docs;    // [k][g]
  const arma::mat& logtot;                              // G x K: log sum_v phi_g[v,k]
  const arma::vec& mv;
  const double ridge, tol; const int maxInner;
  arma::mat& kt_out; arma::mat& ks_out;
  arma::imat& inner_out;                                // Proposition 2 iterations, -1 = not converged

  AdjustedWorker(const std::vector<arma::mat>& phi_, const std::vector< std::vector<GroupDocs> >& docs_,
                 const arma::mat& logtot_, const arma::vec& mv_, double ridge_, double tol_, int maxInner_,
                 arma::mat& kt_, arma::mat& ks_, arma::imat& inner_)
    : phi(phi_), docs(docs_), logtot(logtot_), mv(mv_), ridge(ridge_), tol(tol_), maxInner(maxInner_),
      kt_out(kt_), ks_out(ks_), inner_out(inner_) {}

  void operator()(std::size_t v0, std::size_t v1) {
    const int G = (int)phi.size();
    const int K = (int)logtot.n_cols;
    std::vector<double> y(G), off(G), abar(G), abar_new(G);
    for (std::size_t v = v0; v < v1; ++v) {
      for (int k = 0; k < K; ++k) {
        double ysum = 0.0, musum = 0.0;
        for (int g = 0; g < G; ++g) {
          y[g] = phi[g](v, k);
          off[g] = logtot(g, k) + mv[v];
          if (std::isfinite(off[g])) { ysum += y[g]; musum += std::exp(off[g]); }
        }
        // start: ks = 0 (as in the paper), kt at its ks = 0 maximum-likelihood value
        double ks = 0.0;
        double kt = (ysum > 0.0 && musum > 0.0) ? std::log(ysum / musum) : 0.0;
        for (int g = 0; g < G; ++g) abar[g] = weighted_abar(docs[k][g], 0.0);
        int inner = -1;
        double kt1 = 0.0, ks1 = 0.0;
        for (int it = 1; it <= maxInner; ++it) {
          fit_two_coef(y, off, abar, ridge, kt, ks);
          if (it == 1) { kt1 = kt; ks1 = ks; }
          double delta = 0.0, scale = 1.0;
          for (int g = 0; g < G; ++g) {
            abar_new[g] = weighted_abar(docs[k][g], ks);
            delta = std::max(delta, std::fabs(abar_new[g] - abar[g]));
            scale = std::max(scale, std::fabs(abar[g]));
          }
          if (delta <= tol * scale) { inner = it; break; }
          abar.swap(abar_new);
        }
        // No fixed point reached. Equation (13) need not have one: tilting
        // moves each group's weighted mean towards its largest alphahat, so
        // when the groups overlap heavily abar_1(ks) - abar_2(ks) can change
        // sign and the map ks -> fit(abar(ks)) has a pole and no fixed point
        // (seen on synthetic data with G = 2). Fall back to the first pass,
        // the fit at the phibar-weighted means (ks = 0 weights), which is
        // always defined.
        if (inner < 0) { kt = kt1; ks = ks1; }
        kt_out(v, k) = kt; ks_out(v, k) = ks; inner_out(v, k) = inner;
      }
    }
  }
};

}  // namespace

// kappa M-step for kappaEstimation = "adjusted" (Chen and Mankad 2025, Section
// 4.2 and Proposition 2). phi_list: numGroups V x K matrices of aggregated
// expected word-topic counts; alphaS, phiD: D x K sentiment-discourse values
// and per-document expected topic counts; group: 1-based group of each
// document. Returns the new kappa_t and kappa_s (undamped) and, per (v, k),
// the number of Proposition 2 iterations (-1 where maxInner was reached).
// [[Rcpp::export]]
Rcpp::List sts_optkappa_adjusted_cpp(const Rcpp::List& phi_list, const arma::mat& alphaS,
                                     const arma::mat& phiD, const Rcpp::IntegerVector& group,
                                     const arma::vec& mv, double ridge = 1e-6, double tol = 1e-8,
                                     int maxInner = 100, int cores = 1) {
  const int G = phi_list.size();
  const arma::uword D = alphaS.n_rows, K = alphaS.n_cols;
  if (phiD.n_rows != D || phiD.n_cols != K || (arma::uword)group.size() != D)
    Rcpp::stop("sts_optkappa_adjusted_cpp: alphaS, phiD and group disagree in size");
  std::vector<arma::mat> phi(G);
  for (int g = 0; g < G; ++g) phi[g] = Rcpp::as<arma::mat>(phi_list[g]);
  const arma::uword V = phi[0].n_rows;
  if (mv.n_elem != V) Rcpp::stop("sts_optkappa_adjusted_cpp: mv has the wrong length");

  std::vector< std::vector<GroupDocs> > docs(K, std::vector<GroupDocs>(G));
  for (arma::uword d = 0; d < D; ++d) {
    const int g = group[d] - 1;
    if (g < 0 || g >= G) Rcpp::stop("sts_optkappa_adjusted_cpp: group out of range");
    for (arma::uword k = 0; k < K; ++k) {
      if (phiD(d, k) > 0.0) {
        docs[k][g].a.push_back(alphaS(d, k));
        docs[k][g].lw.push_back(std::log(phiD(d, k)));
      }
    }
  }
  arma::mat logtot(G, K);
  for (int g = 0; g < G; ++g)
    for (arma::uword k = 0; k < K; ++k) {
      const double s = arma::accu(phi[g].col(k));
      logtot(g, k) = (s > 0.0) ? std::log(s) : -std::numeric_limits<double>::infinity();
    }

  arma::mat kt(V, K), ks(V, K);
  arma::imat inner(V, K);
  AdjustedWorker w(phi, docs, logtot, mv, ridge, tol, maxInner, kt, ks, inner);
  RcppParallel::parallelFor(0, V, w, 8, cores < 1 ? 1 : cores);
  return Rcpp::List::create(Rcpp::Named("kappa_t") = kt, Rcpp::Named("kappa_s") = ks,
                            Rcpp::Named("inner") = inner);
}
