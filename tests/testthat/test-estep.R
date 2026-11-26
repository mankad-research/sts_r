## The C++ E-step (src/estep.cpp): thread invariance, the interpolated
## backend against the exact one, and the optional Newton optimizer.

test_that("the E-step gives identical results on one and two threads", {
  f <- gad_fit(); cp <- gad_corpus()
  r1 <- run_estep(f, cp, logZMode = 1L, cores = 1L)
  r2 <- run_estep(f, cp, logZMode = 1L, cores = 2L)
  expect_identical(r1$alpha, r2$alpha)
  expect_identical(r1$elbo, r2$elbo)
  expect_identical(r1$phi, r2$phi)
  expect_identical(r1$sigma_ss, r2$sigma_ss)
})

test_that("the interpolated log-partition reproduces the exact E-step", {
  f <- gad_fit(); cp <- gad_corpus()
  ex <- run_estep(f, cp, logZMode = 0L)
  hm <- run_estep(f, cp, logZMode = 1L)
  expect_equal(hm$elbo, ex$elbo, tolerance = 1e-10)
  expect_equal(hm$alpha, ex$alpha, tolerance = 1e-7)
  for (g in seq_along(ex$phi)) expect_equal(hm$phi[[g]], ex$phi[[g]], tolerance = 1e-8)
  expect_equal(ex$nfail, 0)
  expect_equal(hm$nfail, 0)
})

test_that("the E-step bound decomposes into its data-fit and Laplace terms", {
  f <- gad_fit(); cp <- gad_corpus()
  r <- run_estep(f, cp)
  expect_equal(r$elbo, r$elboObj + r$elboLap, tolerance = 1e-12)
  ## with maxit = 0 the documents are not re-optimised: alpha is returned as given
  r0 <- run_estep(f, cp, maxit = 0L)
  expect_identical(r0$alpha, f$alpha)
})

test_that("the bound is the Laplace ELBO of Chen and Mankad (2025), Proposition 1", {
  ## Reference built from the exported sts 1.6 functions, which form beta over
  ## the full vocabulary: lpbdcpp() is the document objective f_d and
  ## esthcpp() its negative Hessian. Proposition 1:
  ##   ELBO_d = f_d - (1/2) log det(-H_d) - (1/2) log det(Sigma).
  ## Evaluated at the E-step's own modes, where every negative Hessian is
  ## positive definite (away from a mode the E-step adds a diagonal fix-up
  ## before the Cholesky factorization, which the reference does not).
  f <- gad_fit(); cp <- gad_corpus(); V <- length(cp$vocab)
  r <- run_estep(f, cp)
  expect_equal(r$ncholfix, 0)
  ldS <- as.numeric(determinant(f$sigma, logarithm = TRUE)$modulus)
  ref_obj <- ref_ld <- 0
  for (d in seq_along(cp$documents)) {
    doc <- cp$documents[[d]]; a <- r$alpha[d, ]
    ref_obj <- ref_obj + lpbdcpp(a, f$sigma_inv, f$kappa$kappa_t, f$kappa$kappa_s, f$mu[, d], doc, V, f$mv)
    H <- esthcpp(a, f$kappa$kappa_t, f$kappa$kappa_s, f$sigma_inv, doc, V, f$mv)
    ref_ld <- ref_ld + as.numeric(determinant(H, logarithm = TRUE)$modulus)
  }
  D <- length(cp$documents)
  expect_equal(r$elboObj, ref_obj, tolerance = 1e-10)
  expect_equal(r$halfLogdetNegH, ref_ld / 2, tolerance = 1e-10)
  expect_equal(r$elbo, ref_obj - ref_ld / 2 - D * ldS / 2, tolerance = 1e-10)
  ## boundForm = 1 (options(sts.bound = "1.6")) restores the term of sts
  ## 1.0-1.7.0, (1/2)(sum log diag chol(-H_d) - log det Sigma); nothing else
  ## changes.
  r16 <- run_estep(f, cp, boundForm = 1L)
  expect_equal(r16$elboLap, (ref_ld / 2 - D * ldS) / 2, tolerance = 1e-10)
  expect_identical(r16$elboObj, r$elboObj)
  expect_identical(r16$alpha, r$alpha)
  expect_identical(r16$sigma_ss, r$sigma_ss)
})

test_that("the (mu, Sigma) update never lowers the bound", {
  ## It maximizes the bound given the E-step's means and covariances, and
  ## re-evaluating the Laplace covariance at fixed alpha can only raise it.
  ## (The term used before 1.7.1 fell in this step on most runs.)
  cp <- gad_corpus(); set.seed(2024)
  old <- options(sts.traceBound = TRUE); on.exit(options(old))
  f <- suppressWarnings(sts(~ treatment * pid_rep, ~ noTreatment, cp, K = 3, maxIter = 15,
                            convTol = 0, patience = 100L, verbose = FALSE, cores = 1L))
  bt <- f$boundTrace
  expect_true(all(bt$afterMuSigma - bt$afterE >= -1e-9 * abs(bt$afterE)))
})

test_that("option sts.bound is validated", {
  cp <- gad_corpus()
  old <- options(sts.bound = "wrong"); on.exit(options(old))
  expect_error(sts(~ treatment * pid_rep, ~ noTreatment, cp, K = 3, maxIter = 2, verbose = FALSE),
               "sts.bound")
})

test_that("the Newton E-step reaches the BFGS optimum", {
  f <- gad_fit(); cp <- gad_corpus()
  b <- run_estep(f, cp, logZMode = 0L, optimizer = 0L)
  n <- run_estep(f, cp, logZMode = 0L, optimizer = 1L)
  expect_equal(n$elbo, b$elbo, tolerance = 1e-6)
  ## BFGS stops on a relative change in the objective of 1.5e-8, which leaves
  ## alpha accurate to roughly 1e-4 to 1e-3 in weakly identified directions;
  ## Newton converges further. The per-document objective can also have more
  ## than one local maximum, so allow a small fraction of documents to differ.
  d <- apply(abs(n$alpha - b$alpha), 1, max)
  expect_lt(median(d), 1e-3)
  expect_gt(mean(d < 1e-2), 0.98)
})
