// estep.cpp -- the whole per-document E-step in C++.
//
// 1.6 ran, for each document and each of ~45 BFGS objective/gradient calls:
//   * a full V x K softmax inside lpbdcpp,
//   * a second, independent full V x K softmax inside lgaecpp,
//   * a rebuild of the V x K matrix mv + kappa_t on every one of those calls,
// then a full V x K softmax in esthcpp for the Hessian and another in R to
// accumulate phi -- even though only the document's own n_d rows are ever
// used, and the only genuinely full-V quantities are the three scalars per
// topic (logZ_k, logZ_k', logZ_k'') described in normalizer.h.
//
// This file does instead, per document:
//   * form the document's n_d x K slices of (mv + kappa_t) and kappa_s once,
//   * fuse objective and gradient behind a one-point memo so beta is built
//     once per distinct alpha rather than twice,
//   * reuse the same beta_bar for the Hessian and for the phi accumulation
//     (beta_bar on the document's rows is identically the BZ[cdv_ind, ] that
//     1.6 recomputed over all V rows in R),
//   * drive BFGS with the vmmin port in vmmin.h so the algorithm is unchanged.
//
// Exactness: every change here is an algebraic rearrangement or a removal of
// redundant work.  Results match 1.6 up to floating-point reassociation, not
// bit-for-bit, because sums are accumulated in a different order.

#include <RcppArmadillo.h>
#include "normalizer.h"
#include "vmmin.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace sts;

// --------------------------------------------------------------------------
// One document's variational problem: objective, gradient and Hessian, all
// restricted to the document's rows plus the per-topic logZ triples.
// --------------------------------------------------------------------------
class DocProblem {
public:
  DocProblem(const Normalizer& nrm,
             const arma::mat& lb_doc,     // n_d x K   (mv + kappa_t)[ind, ]
             const arma::mat& ks_doc,     // n_d x K   kappa_s[ind, ]
             const arma::vec& cnt,        // n_d       word counts
             const arma::vec& mu_d,       // 2K-1
             const arma::mat& Sigma_Inv)  // (2K-1)^2
    : nrm_(nrm), lb_(lb_doc), ks_(ks_doc), cnt_(cnt), mu_(mu_d), Si_(Sigma_Inv),
      nd_(lb_doc.n_rows), K_(lb_doc.n_cols), C_(arma::accu(cnt)),
      bd_(lb_doc.n_rows, lb_doc.n_cols), bb_(lb_doc.n_rows, lb_doc.n_cols),
      lZ_(lb_doc.n_cols), lZ1_(lb_doc.n_cols), lZ2_(lb_doc.n_cols),
      etaexp_(lb_doc.n_cols), cached_(2 * lb_doc.n_cols - 1),
      has_cache_(false), obj_(0.0), sum_eta_(0.0) {}

  // vmmin minimises, so fn/gr return the negated objective and gradient.
  double fn(const arma::vec& a) { prepare(a); return -obj_; }

  void gr(const arma::vec& a, arma::vec& g) {
    prepare(a);
    const arma::uword K = K_;
    arma::vec full(2 * K - 1);

    // g1_k = sum_i c_i * beta_bar[i,k] - C * theta_k      (first K-1 kept)
    // g2_k = sum_i c_i * beta_bar[i,k] * (kappa_s[i,k] - logZ'_k)
    for (arma::uword k = 0; k < K; ++k) {
      double s1 = 0.0, s2 = 0.0;
      const double* bb = bb_.colptr(k);
      const double* ks = ks_.colptr(k);
      const double  d1 = lZ1_[k];
      for (arma::uword i = 0; i < nd_; ++i) {
        const double w = cnt_[i] * bb[i];
        s1 += w;
        s2 += w * (ks[i] - d1);
      }
      if (k < K - 1) full[k] = s1 - C_ * etaexp_[k] / sum_eta_;
      full[K - 1 + k] = s2;
    }
    g = -(full - Si_ * (a - mu_));
  }

  // Observed-information matrix, identical in structure to 1.6's esthcpp but
  // built from the document's rows and the logZ'' values.
  arma::mat hessian(const arma::vec& a) {
    prepare(a);
    const arma::uword K = K_;

    arma::vec kbar_col(nd_);
    arma::mat A(nd_, K), Ac(nd_, K), Bk(nd_, K), Bkc(nd_, K);
    arma::vec S1(K), S2(K), hss1(K);

    for (arma::uword k = 0; k < K; ++k) {
      const double* bb = bb_.colptr(k);
      const double* ks = ks_.colptr(k);
      const double d1 = lZ1_[k], d2 = lZ2_[k];
      double s1 = 0.0, s2 = 0.0, sh = 0.0;
      for (arma::uword i = 0; i < nd_; ++i) {
        const double kb = ks[i] - d1;
        const double w  = cnt_[i] * bb[i];
        A(i, k)   = bb[i];
        Ac(i, k)  = w;
        Bk(i, k)  = bb[i] * kb;
        Bkc(i, k) = w * kb;
        s1 += w;
        s2 += w * kb;
        sh += w * (kb * kb - d2);
      }
      S1[k] = s1; S2[k] = s2; hss1[k] = sh;
    }

    arma::vec theta = etaexp_ / sum_eta_;
    // The prevalence block's diagonal is esthcpp's g1 = S1 - C*theta, not S1.
    arma::mat h_pp = arma::diagmat(S1 - C_ * theta) - A.t()  * Ac
                     + C_ * theta * theta.t();
    arma::mat h_ps = arma::diagmat(S2)   - Bk.t() * Ac;
    arma::mat h_ss = arma::diagmat(hss1) - Bk.t() * Bkc;

    arma::mat H(2 * K - 1, 2 * K - 1, arma::fill::zeros);
    H.submat(0, 0, K - 2, K - 2)                 = h_pp.submat(0, 0, K - 2, K - 2);
    H.submat(K - 1, 0, 2 * K - 2, K - 2)         = h_ps.cols(0, K - 2);
    H.submat(0, K - 1, K - 2, 2 * K - 2)         = h_ps.cols(0, K - 2).t();
    H.submat(K - 1, K - 1, 2 * K - 2, 2 * K - 2) = h_ss;

    // The h_pp and h_ss blocks come from cross-products (A.t() * Ac and
    // Bk.t() * Bkc), which are symmetric in exact arithmetic but not
    // bit-symmetric: entry (i,j) and (j,i) are accumulated independently and
    // can differ in the last bits.  1.6 handed such a matrix to R's chol(),
    // which runs LAPACK dpotrf on the UPPER triangle and never inspects the
    // lower one -- so R silently symmetrised it upwards.  arma::chol() does
    // check symmetry and would otherwise reject every document, driving the
    // diagonal fix-up below on all D documents and corrupting sigma.ss and
    // the ELBO.  symmatu() reproduces R's semantics exactly.
    return arma::symmatu(-H + Si_);
  }

  // beta_bar on the document's rows == 1.6's BZ[cdv_ind, ].
  const arma::mat& beta_bar() const { return bb_; }
  double objective() const { return obj_; }

private:
  void prepare(const arma::vec& a) {
    if (has_cache_ && arma::approx_equal(a, cached_, "absdiff", 0.0)) return;

    const arma::uword K = K_;
    for (arma::uword k = 0; k < K; ++k) {
      const LogZTriple t = nrm_.eval(k, a[K - 1 + k]);
      lZ_[k] = t.f; lZ1_[k] = t.d1; lZ2_[k] = t.d2;
    }

    sum_eta_ = 0.0;
    for (arma::uword k = 0; k + 1 < K; ++k) { etaexp_[k] = std::exp(a[k]); sum_eta_ += etaexp_[k]; }
    etaexp_[K - 1] = 1.0; sum_eta_ += 1.0;

    // One pass builds beta on the document's rows, its mixture row sums (the
    // objective's word term) and beta_bar -- the work 1.6 did three times.
    double words = 0.0;
    for (arma::uword i = 0; i < nd_; ++i) {
      double rs = 0.0;
      for (arma::uword k = 0; k < K; ++k) {
        const double b = std::exp(lb_(i, k) + ks_(i, k) * a[K - 1 + k] - lZ_[k]);
        bd_(i, k) = b;
        const double t = etaexp_[k] * b;
        bb_(i, k) = t;
        rs += t;
      }
      words += cnt_[i] * std::log(rs);
      for (arma::uword k = 0; k < K; ++k) bb_(i, k) /= rs;
    }

    arma::vec resid = a - mu_;
    obj_ = -0.5 * arma::as_scalar(resid.t() * Si_ * resid)
           + words - C_ * std::log(sum_eta_);

    cached_ = a;
    has_cache_ = true;
  }

  const Normalizer& nrm_;
  const arma::mat& lb_;
  const arma::mat& ks_;
  const arma::vec& cnt_;
  const arma::vec& mu_;
  const arma::mat& Si_;
  arma::uword nd_, K_;
  double C_;
  arma::mat bd_, bb_;
  arma::vec lZ_, lZ1_, lZ2_, etaexp_, cached_;
  bool has_cache_;
  double obj_, sum_eta_;
};

// --------------------------------------------------------------------------
// Full E-step over all documents.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List sts_estep_cpp(const Rcpp::List& documents,
                         const arma::mat& alpha_init,
                         const arma::mat& kappa_t,
                         const arma::mat& kappa_s,
                         const arma::vec& mv,
                         const arma::mat& mu,
                         const arma::mat& Sigma_Inv,
                         double logdet_sigma,
                         const arma::uvec& group,
                         int numGroups,
                         int maxit) {
  const arma::uword D = documents.size();
  const arma::uword V = mv.n_elem;
  const arma::uword K = kappa_t.n_cols;
  const arma::uword P = 2 * K - 1;

  // Precomputed once per E-step instead of once per objective/gradient call.
  arma::mat lin_base(V, K);
  for (arma::uword k = 0; k < K; ++k) lin_base.col(k) = mv + kappa_t.col(k);

  ExactNormalizer nrm(lin_base, kappa_s);

  arma::mat alpha = alpha_init;
  arma::mat sigma_ss(P, P, arma::fill::zeros);
  arma::mat phiD(D, K, arma::fill::zeros);
  std::vector<arma::mat> phi(numGroups, arma::mat(V, K, arma::fill::zeros));
  double elbo = 0.0;
  long fncount = 0, grcount = 0, nfail = 0, ncholfix = 0;

  for (arma::uword d = 0; d < D; ++d) {
    const arma::umat doc = Rcpp::as<arma::umat>(documents[d]);
    const arma::uvec ind = doc.row(0).t() - 1;              // 0-based vocab rows
    const arma::vec  cnt = arma::conv_to<arma::vec>::from(doc.row(1).t());

    // The document's slices, formed once -- not once per BFGS call.
    arma::mat lb_doc = lin_base.rows(ind);
    arma::mat ks_doc = kappa_s.rows(ind);
    arma::vec mu_d   = mu.col(d);

    DocProblem prob(nrm, lb_doc, ks_doc, cnt, mu_d, Sigma_Inv);

    arma::vec a = alpha.row(d).t();
    VmminResult vr = vmmin(a, prob, maxit);
    alpha.row(d) = a.t();
    fncount += vr.fncount; grcount += vr.grcount; nfail += vr.fail;

    arma::mat H = prob.hessian(a);
    arma::mat R;
    if (!arma::chol(R, H)) {
      // 1.6's fix-up: inflate each diagonal entry to at least its off-diagonal
      // absolute row sum, then retry.
      arma::vec dv = H.diag();
      arma::vec mg = arma::sum(arma::abs(H), 1) - arma::abs(dv);
      for (arma::uword i = 0; i < P; ++i) dv[i] = std::max(dv[i], mg[i]);
      H.diag() = dv;
      ncholfix++;
      if (!arma::chol(R, H)) Rcpp::stop("Cholesky failed after diagonal fix-up");
    }
    arma::mat Ri = arma::inv(arma::trimatu(R));
    sigma_ss += Ri * Ri.t();                                 // == chol2inv(R)
    elbo += prob.objective() + 0.5 * (arma::accu(arma::log(R.diag())) - logdet_sigma);

    // phi accumulation reuses beta_bar: no further softmax over V.
    const arma::mat& bb = prob.beta_bar();
    arma::mat contrib = bb;
    contrib.each_col() %= cnt;
    phi[group[d] - 1].rows(ind) += contrib;
    phiD.row(d) = arma::sum(contrib, 0);
  }

  Rcpp::List phi_out(numGroups);
  for (int g = 0; g < numGroups; ++g) phi_out[g] = phi[g];

  return Rcpp::List::create(
    Rcpp::Named("alpha")    = alpha,
    Rcpp::Named("elbo")     = elbo,
    Rcpp::Named("sigma_ss") = sigma_ss,
    Rcpp::Named("phi")      = phi_out,
    Rcpp::Named("phiD")     = phiD,
    Rcpp::Named("fncount")  = (double)fncount,
    Rcpp::Named("grcount")  = (double)grcount,
    Rcpp::Named("nfail")    = (double)nfail,
    Rcpp::Named("ncholfix") = (double)ncholfix);
}
