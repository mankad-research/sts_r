// spectral.cpp -- C++ versions of the two expensive steps of the spectral
// (anchor-word) initialization ported from stm: fastAnchor() and the
// per-word exponentiated-gradient fits inside recoverL2().
//
// Why: at V = 6084 the R versions took 7-20 s (fastAnchor) and 2.5-7 s
// (recoverL2) with peak memory 2.1-2.4 GB, roughly ten dense V x V copies:
// each fastAnchor pass materialises Qbar^2, an outer product and the
// difference. At stm's 10,000-word cap that is ~6-8 GB.
//
// Fidelity: the anchors are chosen by argmax, so a last-bit difference could
// change the initialization. These ports therefore reproduce R's arithmetic
// exactly on a given platform: the matrix-vector products use the same BLAS
// routine R's %*% dispatches to (dgemv), row sums and sums accumulate in long
// double as R's rowSums()/sum() do, and products are formed before they are
// subtracted, in R's order.

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

namespace {

// rowSums(Q^2) with R's accumulation: long double per row, columns in order.
struct RowSqWorker : public RcppParallel::Worker {
  const arma::mat& Q;
  std::vector<double>& out;
  RowSqWorker(const arma::mat& Q_, std::vector<double>& out_) : Q(Q_), out(out_) {}
  void operator()(std::size_t i0, std::size_t i1) {
    const arma::uword V = Q.n_cols;
    for (std::size_t i = i0; i < i1; ++i) {
      long double acc = 0.0L;
      for (arma::uword j = 0; j < V; ++j) {
        const double x = Q(i, j);
        const double sq = x * x;
        acc += sq;
      }
      out[i] = (double)acc;
    }
  }
};

// Qbar <- Qbar - ip %o% q, column by column (disjoint columns per thread).
struct Rank1Worker : public RcppParallel::Worker {
  arma::mat& Q;
  const std::vector<double>& ip;
  const std::vector<double>& q;
  Rank1Worker(arma::mat& Q_, const std::vector<double>& ip_, const std::vector<double>& q_)
    : Q(Q_), ip(ip_), q(q_) {}
  void operator()(std::size_t j0, std::size_t j1) {
    const arma::uword V = Q.n_rows;
    for (std::size_t j = j0; j < j1; ++j) {
      double* col = Q.colptr(j);
      const double qj = q[j];
      for (arma::uword i = 0; i < V; ++i) {
        const double p = ip[i] * qj;   // the outer-product entry, as R forms it
        col[i] = col[i] - p;
      }
    }
  }
};

}  // namespace

// Greedy anchor selection (stabilised Gram-Schmidt), as stm's fastAnchor().
// Works on one internal copy of Qbar; returns 1-based anchor indices.
// [[Rcpp::export]]
Rcpp::IntegerVector sts_fastanchor_cpp(const arma::mat& Qbar_in, int K, int cores = 1) {
  arma::mat Q = Qbar_in;
  const int V = (int)Q.n_rows;
  if (cores < 1) cores = 1;
  std::vector<double> rss(V), ip(V), q(V);
  std::vector<int> basis;
  std::vector<char> is_basis(V, 0);

  RowSqWorker w0(Q, rss);
  RcppParallel::parallelFor(0, V, w0, 64, cores);

  for (int it = 0; it < K; ++it) {
    // which.max: first index of the maximum
    int b = 0; double mx = rss[0];
    for (int i = 1; i < V; ++i) if (rss[i] > mx) { mx = rss[i]; b = i; }
    basis.push_back(b); is_basis[b] = 1;
    const double normalizer = 1.0 / std::sqrt(mx);
    for (int j = 0; j < V; ++j) Q(b, j) = Q(b, j) * normalizer;
    for (int j = 0; j < V; ++j) q[j] = Q(b, j);

    // innerproducts <- Qbar %*% Qbar[b, ]. R dispatches this to BLAS dgemv
    // ("N", alpha 1, beta 0); Armadillo's matrix-vector product makes the
    // same call for a matrix this size, so the result is identical.
    {
      const arma::vec qv(q.data(), V, false, true);
      const arma::vec ipv = Q * qv;
      for (int i = 0; i < V; ++i) ip[i] = ipv[i];
    }
    for (int bb : basis) ip[bb] = 0.0;                   // project[basis, ] <- 0

    Rank1Worker w1(Q, ip, q);
    RcppParallel::parallelFor(0, V, w1, 16, cores);

    RowSqWorker w2(Q, rss);
    RcppParallel::parallelFor(0, V, w2, 64, cores);
    for (int bb : basis) rss[bb] = 0.0;
  }
  Rcpp::IntegerVector out(K);
  for (int i = 0; i < K; ++i) out[i] = basis[i] + 1;
  return out;
}

namespace {

// One word's exponentiated-gradient fit, as stm's expgrad() with L2 loss.
// ytX: that word's row of Qbar %*% t(X) (K); XtX: K x K.
struct ExpgradWorker : public RcppParallel::Worker {
  const arma::mat& ytX;      // V x K
  const arma::mat& XtX;      // K x K
  const std::vector<int>& anchor_pos;   // -1, or position among the anchors
  const double tol; const int max_its;
  arma::mat& W;              // V x K output weights
  std::vector<int>& its_out;

  ExpgradWorker(const arma::mat& ytX_, const arma::mat& XtX_, const std::vector<int>& ap,
                double tol_, int max_its_, arma::mat& W_, std::vector<int>& its_)
    : ytX(ytX_), XtX(XtX_), anchor_pos(ap), tol(tol_), max_its(max_its_), W(W_), its_out(its_) {}

  void operator()(std::size_t i0, std::size_t i1) {
    const int K = (int)XtX.n_rows;
    std::vector<double> alpha(K), grad(K), z(K), yx(K);
    for (std::size_t i = i0; i < i1; ++i) {
      if (anchor_pos[i] >= 0) {
        for (int k = 0; k < K; ++k) W(i, k) = 0.0;
        W(i, anchor_pos[i]) = 1.0;
        its_out[i] = 0;
        continue;
      }
      for (int k = 0; k < K; ++k) { alpha[k] = 1.0 / (double)K; yx[k] = ytX(i, k); }
      const double eta = 50.0;
      double sse_old = std::numeric_limits<double>::infinity();
      int its = 1;
      for (;;) {
        // alpha %*% XtX  (1 x K times K x K: R uses dgemv "T"; reference
        // dgemv accumulates each output from zero in row order)
        for (int j = 0; j < K; ++j) {
          double t = 0.0;
          for (int r = 0; r < K; ++r) { const double p = XtX(r, j) * alpha[r]; t = t + p; }
          z[j] = t;
        }
        long double ss = 0.0L;
        for (int k = 0; k < K; ++k) { grad[k] = yx[k] - z[k]; const double g2 = grad[k] * grad[k]; ss += g2; }
        const double sse = (double)ss;
        const double c = 2.0 * eta;
        double mx = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < K; ++k) { grad[k] = c * grad[k]; if (grad[k] > mx) mx = grad[k]; }
        long double sa = 0.0L;
        for (int k = 0; k < K; ++k) { alpha[k] = alpha[k] * std::exp(grad[k] - mx); sa += alpha[k]; }
        const double s = (double)sa;
        for (int k = 0; k < K; ++k) alpha[k] = alpha[k] / s;
        const bool converged = std::fabs(std::sqrt(sse_old) - std::sqrt(sse)) < tol;
        if (converged || its == max_its) break;
        sse_old = sse;
        ++its;
      }
      for (int k = 0; k < K; ++k)
        W(i, k) = (alpha[k] <= 0.0) ? std::numeric_limits<double>::epsilon() : alpha[k];
      its_out[i] = its;
    }
  }
};

}  // namespace

// The per-word fits of recoverL2(): returns the V x K weight matrix (anchors
// get unit vectors; non-positive weights are set to .Machine$double.eps, as in
// the R code). Words are independent, so this runs in parallel and the result
// does not depend on the thread count.
// [[Rcpp::export]]
Rcpp::List sts_expgrad_rows_cpp(const arma::mat& ytX, const arma::mat& XtX,
                                const Rcpp::IntegerVector& anchor,
                                double tol = 1e-7, int max_its = 500, int cores = 1) {
  const arma::uword V = ytX.n_rows, K = ytX.n_cols;
  std::vector<int> ap(V, -1);
  for (int a = 0; a < anchor.size(); ++a) ap[anchor[a] - 1] = a;
  arma::mat W(V, K);
  std::vector<int> its(V);
  ExpgradWorker w(ytX, XtX, ap, tol, max_its, W, its);
  RcppParallel::parallelFor(0, V, w, 16, cores < 1 ? 1 : cores);
  return Rcpp::List::create(Rcpp::Named("weights") = W, Rcpp::Named("its") = its);
}
