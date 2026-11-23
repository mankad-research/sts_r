## The spectral (anchor-word) initialization: the C++ fastAnchor() must choose
## the same anchors as stm's R code, and the C++ recoverL2() fits must match
## their R references (sts 1.7's simplex least squares by default, stm's
## exponentiated gradient under options(sts.recoverL2 = "expgrad")).

spectral_Q <- function() {
  cp <- gad_corpus()
  docs <- sts:::doc.to.ijv(cp$documents)
  mat <- Matrix::sparseMatrix(docs$i, docs$j, x = docs$v)
  wprob <- Matrix::colSums(mat); wprob <- wprob / sum(wprob)
  Q <- sts:::gram(mat); Q <- Q / rowSums(Q)
  list(Q = Q, wprob = wprob)
}

test_that("fastAnchor in C++ selects the same anchors as the R version", {
  s <- spectral_Q()
  for (K in c(3, 10)) {
    a_R <- sts:::fastAnchor_R(s$Q, K)
    expect_identical(sts:::fastAnchor(s$Q, K, 1L), as.integer(a_R))
    expect_identical(sts:::fastAnchor(s$Q, K, 2L), as.integer(a_R))
  }
})

## recoverL2 (sts 1.7 default): simplex-constrained least squares by projected
## gradient. The C++ fits must match the R reference, reach the optimum (KKT
## conditions), and be insensitive to rounding-level changes in the inputs,
## which is what makes the initialization reproducible across platforms.
test_that("recoverL2 in C++ reproduces the R version", {
  s <- spectral_Q()
  for (K in c(3, 5, 10)) {
    a <- sts:::fastAnchor_R(s$Q, K)
    A_R <- sts:::recoverL2_R(s$Q, a, s$wprob)$A
    A_1 <- sts:::recoverL2(s$Q, a, s$wprob, 1L)$A
    expect_equal(A_1, A_R, tolerance = 1e-12)
  }
})

test_that("the simplex least-squares fits converge to the optimum", {
  s <- spectral_Q()
  a <- sts:::fastAnchor_R(s$Q, 10)
  X <- s$Q[a, ]; XtX <- tcrossprod(X); ytX <- s$Q %*% t(X)
  r <- sts:::sts_simplexls_rows_cpp(ytX, XtX, as.integer(a), sts:::simplexLS_L(XtX),
                                    1e-14, 1000L, 1L)
  expect_true(all(r$converged == 1L))
  words <- setdiff(seq_len(nrow(s$Q)), a)
  ## KKT: the gradient equals its minimum on the support and is no smaller off it
  kkt <- vapply(words, function(i) {
    w <- r$weights[i, ]; w[w <= .Machine$double.eps] <- 0
    g <- 2 * (drop(w %*% XtX) - ytX[i, ]); m <- min(g); on <- w > 0
    max(abs(g[on] - m), pmax(m - g[!on], 0))
  }, numeric(1))
  expect_lt(max(kkt), 1e-12)
  ## at least as good as stm's exponentiated gradient on every word
  f <- function(w, i) sum((s$Q[i, ] - drop(w %*% X))^2)
  eg <- sts:::sts_expgrad_rows_cpp(ytX, XtX, as.integer(a), 1e-7, 500L, 1L)$weights
  dobj <- vapply(words, function(i) f(eg[i, ], i) - f(r$weights[i, ], i), numeric(1))
  expect_true(all(dobj > -1e-12))
  ## a 1e-15 relative perturbation of the inputs barely moves the solution
  set.seed(1)
  P <- ytX * (1 + 1e-15 * rnorm(length(ytX)))
  r2 <- sts:::sts_simplexls_rows_cpp(P, XtX, as.integer(a), sts:::simplexLS_L(XtX),
                                     1e-14, 1000L, 1L)
  expect_lt(max(abs(r$weights - r2$weights)), 1e-12)
})

## stm's exponentiated-gradient fit (eta = 50), kept for
## options(sts.recoverL2 = "expgrad"). It usually does not converge on these
## data: the iterates oscillate until max.its, and rounding-level differences
## grow by about two orders of magnitude every few iterations, so after 500
## iterations the C++ and R results agree only where the arithmetic is
## bit-identical. The port is compared with expgrad() over the first 10
## iterations, where rounding differences are still below 1e-10.
test_that("the C++ exponentiated-gradient fit reproduces expgrad()", {
  s <- spectral_Q()
  a <- sts:::fastAnchor_R(s$Q, 5)
  X <- s$Q[a, ]
  XtX <- tcrossprod(X)
  ytX <- s$Q %*% t(X)
  W <- sts:::sts_expgrad_rows_cpp(ytX, XtX, as.integer(a), 1e-7, 10L, 1L)$weights
  words <- setdiff(seq_len(nrow(s$Q)), a)
  W_R <- t(vapply(words, function(i) {
    w <- sts:::expgrad(X, s$Q[i, ], XtX, max.its = 10)$par
    w[w <= 0] <- .Machine$double.eps
    w
  }, numeric(length(a))))
  expect_equal(W[words, ], W_R, tolerance = 1e-9)
  expect_equal(W[a, ], diag(length(a)))
})

test_that("recoverL2 in C++ does not depend on the number of threads", {
  s <- spectral_Q()
  a <- sts:::fastAnchor_R(s$Q, 5)
  for (solver in c("simplex", "expgrad")) {
    old <- options(sts.recoverL2 = solver)
    A_1 <- sts:::recoverL2(s$Q, a, s$wprob, 1L)$A
    A_2 <- sts:::recoverL2(s$Q, a, s$wprob, 2L)$A
    options(old)
    expect_identical(A_1, A_2)
    expect_equal(dim(A_1), c(5L, nrow(s$Q)))
    expect_equal(rowSums(A_1), rep(1, 5), tolerance = 1e-12)
  }
})
