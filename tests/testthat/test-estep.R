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
