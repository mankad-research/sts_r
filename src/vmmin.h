// vmmin.h -- port of R's BFGS minimiser (vmmin, R/src/appl/optim.c) so that
// moving the E-step into C++ does not change the optimisation algorithm.
//
// 1.6 called stats::optim(method = "BFGS", control = list(fnscale = -1,
// maxit = 500)) once per document.  optim() minimises fn/fnscale, so with
// fnscale = -1 it maximises; the port below minimises the negated objective
// and the caller negates the result back.  Defaults mirror optim()'s BFGS
// defaults: reltol = sqrt(.Machine$double.eps), abstol = -Inf, and the
// unit parscale / all-free mask that sts() uses, so the iterate sequence is
// identical to 1.6's up to floating-point reassociation inside the objective.
//
// Reference: R 4.3.3, src/appl/optim.c, function vmmin().  Constants
// stepredn, acctol and reltest keep their R values.

#ifndef STS_VMMIN_H
#define STS_VMMIN_H

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

namespace sts {

static const double vm_stepredn = 0.2;
static const double vm_acctol   = 0.0001;
static const double vm_reltest  = 10.0;

struct VmminResult {
  double value;    // minimum of the *negated* objective found
  int fncount;
  int grcount;
  int fail;        // 0 = converged within maxit, 1 = hit maxit
};

// Problem must provide:
//   double fn(const arma::vec& par)                       // to MINIMISE
//   void   gr(const arma::vec& par, arma::vec& g)         // gradient of fn
template <class Problem>
VmminResult vmmin(arma::vec& b, Problem& prob, int maxit = 500,
                  double abstol = -std::numeric_limits<double>::infinity(),
                  double reltol = 1.490116119384766e-08) {
  const int n = static_cast<int>(b.n_elem);
  VmminResult res;
  res.fncount = res.grcount = 0;
  res.fail = 0;

  if (maxit <= 0) {
    res.value = prob.fn(b);
    return res;
  }

  arma::vec g(n, arma::fill::zeros);
  std::vector<double> t(n), X(n), c(n);
  // B is the lower triangle of the inverse-Hessian approximation, as in R.
  std::vector< std::vector<double> > B(n, std::vector<double>(n, 0.0));

  bool accpoint, enough;
  double f, gradproj, s, steplength, D1, D2;
  int count, funcount, gradcount, iter = 0, ilast;

  f = prob.fn(b);
  if (!std::isfinite(f)) Rcpp::stop("initial value in 'vmmin' is not finite");
  res.value = f;
  funcount = gradcount = 1;
  prob.gr(b, g);
  iter++;
  ilast = gradcount;

  do {
    if (ilast == gradcount) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) B[i][j] = 0.0;
        B[i][i] = 1.0;
      }
    }
    for (int i = 0; i < n; i++) { X[i] = b[i]; c[i] = g[i]; }

    gradproj = 0.0;
    for (int i = 0; i < n; i++) {
      s = 0.0;
      for (int j = 0; j <= i; j++)     s -= B[i][j] * g[j];
      for (int j = i + 1; j < n; j++)  s -= B[j][i] * g[j];
      t[i] = s;
      gradproj += s * g[i];
    }

    if (gradproj < 0.0) {              // search direction is downhill
      steplength = 1.0;
      accpoint = false;
      do {
        count = 0;
        for (int i = 0; i < n; i++) {
          b[i] = X[i] + steplength * t[i];
          if (vm_reltest + X[i] == vm_reltest + b[i]) count++;   // no change
        }
        if (count < n) {
          f = prob.fn(b);
          funcount++;
          accpoint = (std::isfinite(f) &&
                      f <= res.value + gradproj * steplength * vm_acctol);
          if (!accpoint) steplength *= vm_stepredn;
        }
      } while (!(count == n || accpoint));

      enough = (f > abstol) &&
               std::fabs(f - res.value) > reltol * (std::fabs(res.value) + reltol);
      if (!enough) { count = n; res.value = f; }

      if (count < n) {                 // making progress
        res.value = f;
        prob.gr(b, g);
        gradcount++;
        iter++;
        D1 = 0.0;
        for (int i = 0; i < n; i++) {
          t[i] = steplength * t[i];
          c[i] = g[i] - c[i];
          D1 += t[i] * c[i];
        }
        if (D1 > 0) {
          D2 = 0.0;
          for (int i = 0; i < n; i++) {
            s = 0.0;
            for (int j = 0; j <= i; j++)    s += B[i][j] * c[j];
            for (int j = i + 1; j < n; j++) s += B[j][i] * c[j];
            X[i] = s;
            D2 += s * c[i];
          }
          D2 = 1.0 + D2 / D1;
          for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
              B[i][j] += (D2 * t[i] * t[j] - X[i] * t[j] - t[i] * X[j]) / D1;
            }
          }
        } else {
          ilast = gradcount;
        }
      } else {                         // no progress
        if (ilast < gradcount) { count = 0; ilast = gradcount; }
      }
    } else {                           // uphill search
      count = 0;
      if (ilast == gradcount) count = n;
      else ilast = gradcount;
    }

    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * n) ilast = gradcount;   // periodic restart
  } while (count != n || ilast != gradcount);

  res.fail = (iter < maxit) ? 0 : 1;
  res.fncount = funcount;
  res.grcount = gradcount;
  return res;
}

}  // namespace sts

#endif  // STS_VMMIN_H
