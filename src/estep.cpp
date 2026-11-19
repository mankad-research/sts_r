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
#include <RcppParallel.h>
#include <memory>
#include <string>
#include "normalizer.h"
#include "hermite.h"
#include "vmmin.h"
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

using namespace sts;

// Direct access to the two normalizer backends, for the accuracy and
// convergence-rate tests in tests/testthat and for the error analysis.
// mode 0 = exact (one pass over V), mode 1 = quintic Hermite table.
// [[Rcpp::export]]
Rcpp::List sts_logZ_probe(const arma::mat& kappa_t, const arma::mat& kappa_s,
                          const arma::vec& mv, const arma::vec& avals,
                          int mode, double a_lo, double a_hi, double spacing) {
  const arma::uword V = mv.n_elem, K = kappa_t.n_cols, n = avals.n_elem;
  arma::mat lin_base(V, K);
  for (arma::uword k = 0; k < K; ++k) lin_base.col(k) = mv + kappa_t.col(k);

  ExactNormalizer ex(lin_base, kappa_s);
  arma::mat f(n, K), d1(n, K), d2(n, K);
  double nodes = NA_REAL, lo = NA_REAL, hi = NA_REAL, nout = NA_REAL;

  if (mode == 0) {
    for (arma::uword i = 0; i < n; ++i)
      for (arma::uword k = 0; k < K; ++k) {
        const LogZTriple t = ex.eval(k, avals[i]);
        f(i, k) = t.f; d1(i, k) = t.d1; d2(i, k) = t.d2;
      }
  } else {
    HermiteNormalizer hm(lin_base, kappa_s, a_lo, a_hi, spacing);
    for (arma::uword i = 0; i < n; ++i)
      for (arma::uword k = 0; k < K; ++k) {
        const LogZTriple t = hm.eval(k, avals[i]);
        f(i, k) = t.f; d1(i, k) = t.d1; d2(i, k) = t.d2;
      }
    nodes = (double)hm.n_nodes(); lo = hm.range_lo(); hi = hm.range_hi();
    nout = (double)hm.n_outside();
  }

  return Rcpp::List::create(
    Rcpp::Named("f") = f, Rcpp::Named("d1") = d1, Rcpp::Named("d2") = d2,
    Rcpp::Named("n_nodes") = nodes, Rcpp::Named("lo") = lo,
    Rcpp::Named("hi") = hi, Rcpp::Named("n_outside") = nout);
}

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
  // logZ_k at the last prepared alpha (the optimum, once hessian() has run)
  const arma::vec& logZ() const { return lZ_; }

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
// Full E-step over all documents, multithreaded with RcppParallel.
//
// Determinism: results must not depend on the number of threads, so nothing
// is accumulated in a thread-dependent order.
//  * Documents are split into a FIXED set of blocks determined by D alone.
//    Each block is processed serially by one thread, accumulating its own
//    sigma.ss and ELBO terms, and the blocks are summed in block order.
//  * Per-document outputs (alpha, phiD, logZ at the optimum) go to disjoint
//    rows.
//  * phi, where documents collide on shared vocabulary rows, is built in a
//    second pass that is parallel over WORDS (disjoint rows of phi), adding a
//    word's documents in ascending order -- the same order and arithmetic as
//    the serial loop, so phi is bitwise identical to it.
// No R API is touched inside the threads: documents are converted up front,
// and errors are recorded per block and raised after the parallel region.
// --------------------------------------------------------------------------
namespace {

struct BlockOut {
  arma::mat sigma_ss;
  double elbo = 0.0, elbo_obj = 0.0, elbo_lap = 0.0;
  long fncount = 0, grcount = 0, nfail = 0, ncholfix = 0;
  std::string err;
};

struct EstepWorker : public RcppParallel::Worker {
  const std::vector<arma::uvec>& ind;
  const std::vector<arma::vec>&  cnt;
  const Normalizer& nrm;
  const arma::mat& lin_base;
  const arma::mat& kappa_s;
  const arma::mat& mu;
  const arma::mat& Sigma_Inv;
  const double logdet_sigma;
  const int maxit;
  const std::vector<arma::uword>& block_start;
  arma::mat& alpha;
  arma::mat& phiD;
  arma::mat& lZ_opt;
  std::vector<BlockOut>& out;
  const arma::uword P;

  EstepWorker(const std::vector<arma::uvec>& ind_, const std::vector<arma::vec>& cnt_,
              const Normalizer& nrm_, const arma::mat& lin_base_, const arma::mat& kappa_s_,
              const arma::mat& mu_, const arma::mat& Sigma_Inv_, double logdet_sigma_,
              int maxit_, const std::vector<arma::uword>& block_start_,
              arma::mat& alpha_, arma::mat& phiD_, arma::mat& lZ_opt_,
              std::vector<BlockOut>& out_, arma::uword P_)
    : ind(ind_), cnt(cnt_), nrm(nrm_), lin_base(lin_base_), kappa_s(kappa_s_),
      mu(mu_), Sigma_Inv(Sigma_Inv_), logdet_sigma(logdet_sigma_), maxit(maxit_),
      block_start(block_start_), alpha(alpha_), phiD(phiD_), lZ_opt(lZ_opt_),
      out(out_), P(P_) {}

  void operator()(std::size_t b0, std::size_t b1) {
    for (std::size_t b = b0; b < b1; ++b) {
      BlockOut& bo = out[b];
      bo.sigma_ss.zeros(P, P);
      try {
        for (arma::uword d = block_start[b]; d < block_start[b + 1]; ++d) {
          // The document's slices, formed once -- not once per BFGS call.
          arma::mat lb_doc = lin_base.rows(ind[d]);
          arma::mat ks_doc = kappa_s.rows(ind[d]);
          arma::vec mu_d   = mu.col(d);
          DocProblem prob(nrm, lb_doc, ks_doc, cnt[d], mu_d, Sigma_Inv);

          arma::vec a = alpha.row(d).t();
          VmminResult vr = vmmin(a, prob, maxit);
          alpha.row(d) = a.t();
          bo.fncount += vr.fncount; bo.grcount += vr.grcount; bo.nfail += vr.fail;

          arma::mat H = prob.hessian(a);
          arma::mat R;
          if (!arma::chol(R, H)) {
            // 1.6's fix-up: inflate each diagonal entry to at least its
            // off-diagonal absolute row sum, then retry.
            arma::vec dv = H.diag();
            arma::vec mg = arma::sum(arma::abs(H), 1) - arma::abs(dv);
            for (arma::uword i = 0; i < P; ++i) dv[i] = std::max(dv[i], mg[i]);
            H.diag() = dv;
            bo.ncholfix++;
            if (!arma::chol(R, H))
              throw std::runtime_error("Cholesky failed after diagonal fix-up");
          }
          arma::mat Ri;
          if (!arma::inv(Ri, arma::trimatu(R)))
            throw std::runtime_error("inverse of the Cholesky factor failed");
          bo.sigma_ss += Ri * Ri.t();                        // == chol2inv(R)
          const double lap = 0.5 * (arma::accu(arma::log(R.diag())) - logdet_sigma);
          bo.elbo     += prob.objective() + lap;
          bo.elbo_obj += prob.objective();
          bo.elbo_lap += lap;

          // phiD reuses beta_bar; phi is built in the word pass below.
          arma::mat contrib = prob.beta_bar();
          contrib.each_col() %= cnt[d];
          phiD.row(d)   = arma::sum(contrib, 0);
          lZ_opt.row(d) = prob.logZ().t();
        }
      } catch (std::exception& e) {
        bo.err = e.what();
      } catch (...) {
        bo.err = "unknown error in the E-step";
      }
    }
  }
};

// Second pass: phi, parallel over vocabulary words. Replicates the arithmetic
// of DocProblem::prepare() at the optimum exactly, row by row.
struct PhiWorker : public RcppParallel::Worker {
  const std::vector<arma::uword>& wptr;   // CSR: documents containing word v
  const std::vector<arma::uword>& wdoc;
  const std::vector<double>&      wcnt;
  const arma::mat& lin_base;
  const arma::mat& kappa_s;
  const arma::mat& alpha;
  const arma::mat& lZ_opt;
  const arma::uvec& group;                // 1-based
  std::vector<arma::mat>& phi;
  const arma::uword K;

  PhiWorker(const std::vector<arma::uword>& wptr_, const std::vector<arma::uword>& wdoc_,
            const std::vector<double>& wcnt_, const arma::mat& lin_base_,
            const arma::mat& kappa_s_, const arma::mat& alpha_, const arma::mat& lZ_opt_,
            const arma::uvec& group_, std::vector<arma::mat>& phi_, arma::uword K_)
    : wptr(wptr_), wdoc(wdoc_), wcnt(wcnt_), lin_base(lin_base_), kappa_s(kappa_s_),
      alpha(alpha_), lZ_opt(lZ_opt_), group(group_), phi(phi_), K(K_) {}

  void operator()(std::size_t v0, std::size_t v1) {
    std::vector<double> t(K);
    for (std::size_t v = v0; v < v1; ++v) {
      for (arma::uword p = wptr[v]; p < wptr[v + 1]; ++p) {
        const arma::uword d = wdoc[p];
        double rs = 0.0;
        for (arma::uword k = 0; k < K; ++k) {
          const double eta = (k + 1 < K) ? std::exp(alpha(d, k)) : 1.0;
          const double bk  = std::exp(lin_base(v, k) + kappa_s(v, k) * alpha(d, K - 1 + k)
                                      - lZ_opt(d, k));
          t[k] = eta * bk;
          rs += t[k];
        }
        arma::mat& ph = phi[group[d] - 1];
        for (arma::uword k = 0; k < K; ++k) ph(v, k) += (t[k] / rs) * wcnt[p];
      }
    }
  }
};

}  // namespace

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
                         int maxit,
                         int logZMode,
                         double logZSpacing,
                         double logZMargin,
                         int cores = 1) {
  const arma::uword D = documents.size();
  const arma::uword V = mv.n_elem;
  const arma::uword K = kappa_t.n_cols;
  const arma::uword P = 2 * K - 1;
  if (cores < 1) cores = 1;

  // Documents converted on the main thread: no R API inside the workers.
  std::vector<arma::uvec> ind(D);
  std::vector<arma::vec>  cnt(D);
  for (arma::uword d = 0; d < D; ++d) {
    const arma::umat doc = Rcpp::as<arma::umat>(documents[d]);
    ind[d] = doc.row(0).t() - 1;                            // 0-based vocab rows
    cnt[d] = arma::conv_to<arma::vec>::from(doc.row(1).t());
  }

  // Precomputed once per E-step instead of once per objective/gradient call.
  arma::mat lin_base(V, K);
  for (arma::uword k = 0; k < K; ++k) lin_base.col(k) = mv + kappa_t.col(k);

  // Normalizer backend.  logZMode 0 keeps the exact pass over V; mode 1
  // tabulates logZ, logZ' and logZ'' once for this E-step over the range of
  // alpha_s the documents currently occupy, padded by logZMargin, and
  // interpolates.  Anything outside the table still goes to the exact
  // backend, so no value is ever extrapolated.  (The table is built before
  // the parallel region; lookups are read-only.)
  double a_lo = 0.0, a_hi = 0.0;
  if (logZMode == 1) {
    const arma::mat as = alpha_init.cols(K - 1, P - 1);
    a_lo = as.min() - logZMargin;
    a_hi = as.max() + logZMargin;
  }
  ExactNormalizer   nrm_exact(lin_base, kappa_s);
  std::unique_ptr<HermiteNormalizer> nrm_herm;
  if (logZMode == 1) {
    nrm_herm.reset(new HermiteNormalizer(lin_base, kappa_s, a_lo, a_hi, logZSpacing, cores));
  }
  const Normalizer& nrm = (logZMode == 1)
      ? static_cast<const Normalizer&>(*nrm_herm)
      : static_cast<const Normalizer&>(nrm_exact);

  // ---- pass 1: documents, in a fixed block decomposition ----
  const arma::uword nb = std::max<arma::uword>(1, std::min<arma::uword>(D, 256));
  std::vector<arma::uword> block_start(nb + 1);
  for (arma::uword b = 0; b <= nb; ++b) block_start[b] = (b * D) / nb;

  arma::mat alpha = alpha_init;
  arma::mat phiD(D, K, arma::fill::zeros);
  arma::mat lZ_opt(D, K, arma::fill::zeros);
  std::vector<BlockOut> out(nb);

  EstepWorker ew(ind, cnt, nrm, lin_base, kappa_s, mu, Sigma_Inv, logdet_sigma,
                 maxit, block_start, alpha, phiD, lZ_opt, out, P);
  RcppParallel::parallelFor(0, nb, ew, 1, cores);

  arma::mat sigma_ss(P, P, arma::fill::zeros);
  double elbo = 0.0, elbo_obj = 0.0, elbo_lap = 0.0;
  long fncount = 0, grcount = 0, nfail = 0, ncholfix = 0;
  for (arma::uword b = 0; b < nb; ++b) {                    // fixed order
    if (!out[b].err.empty()) Rcpp::stop("E-step: " + out[b].err);
    sigma_ss += out[b].sigma_ss;
    elbo += out[b].elbo; elbo_obj += out[b].elbo_obj; elbo_lap += out[b].elbo_lap;
    fncount += out[b].fncount; grcount += out[b].grcount;
    nfail += out[b].nfail; ncholfix += out[b].ncholfix;
  }

  // ---- pass 2: phi, parallel over words ----
  std::vector<arma::uword> wptr(V + 1, 0);
  for (arma::uword d = 0; d < D; ++d)
    for (arma::uword i = 0; i < ind[d].n_elem; ++i) wptr[ind[d][i] + 1]++;
  for (arma::uword v = 0; v < V; ++v) wptr[v + 1] += wptr[v];
  std::vector<arma::uword> wdoc(wptr[V]), fill(wptr.begin(), wptr.end() - 1);
  std::vector<double> wcnt(wptr[V]);
  for (arma::uword d = 0; d < D; ++d)                       // ascending d per word
    for (arma::uword i = 0; i < ind[d].n_elem; ++i) {
      const arma::uword v = ind[d][i], p = fill[v]++;
      wdoc[p] = d; wcnt[p] = cnt[d][i];
    }

  std::vector<arma::mat> phi(numGroups, arma::mat(V, K, arma::fill::zeros));
  PhiWorker pw(wptr, wdoc, wcnt, lin_base, kappa_s, alpha, lZ_opt, group, phi, K);
  RcppParallel::parallelFor(0, V, pw, 64, cores);

  Rcpp::List phi_out(numGroups);
  for (int g = 0; g < numGroups; ++g) phi_out[g] = phi[g];

  return Rcpp::List::create(
    Rcpp::Named("alpha")    = alpha,
    Rcpp::Named("elbo")     = elbo,
    Rcpp::Named("elboObj")  = elbo_obj,
    Rcpp::Named("elboLap")  = elbo_lap,
    Rcpp::Named("sigma_ss") = sigma_ss,
    Rcpp::Named("phi")      = phi_out,
    Rcpp::Named("phiD")     = phiD,
    Rcpp::Named("fncount")  = (double)fncount,
    Rcpp::Named("grcount")  = (double)grcount,
    Rcpp::Named("nfail")    = (double)nfail,
    Rcpp::Named("ncholfix") = (double)ncholfix,
    Rcpp::Named("logZNodes")   = logZMode == 1 ? (double)nrm_herm->n_nodes() : NA_REAL,
    Rcpp::Named("logZOutside") = logZMode == 1 ? (double)nrm_herm->n_outside() : NA_REAL,
    Rcpp::Named("logZLo")      = logZMode == 1 ? a_lo : NA_REAL,
    Rcpp::Named("logZHi")      = logZMode == 1 ? a_hi : NA_REAL);
}
