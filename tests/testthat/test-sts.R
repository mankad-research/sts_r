## End-to-end behaviour of sts(): output structure, reproducibility, thread
## invariance, the convergence rule and argument handling.

fit_gad <- function(..., K = 3, maxIter = 8, seed = 11) {
  cp <- gad_corpus(); set.seed(seed)
  suppressWarnings(sts(~ treatment * pid_rep, ~ noTreatment, cp, K = K, maxIter = maxIter,
                       verbose = FALSE, ...))
}

test_that("sts() returns a well-formed STS object", {
  f <- gad_fit(); cp <- gad_corpus()
  D <- length(cp$documents); V <- length(cp$vocab); K <- 3
  expect_s3_class(f, "STS")
  expect_equal(dim(f$alpha), c(D, 2 * K - 1))
  expect_equal(dim(f$kappa$kappa_t), c(V, K))
  expect_equal(dim(f$kappa$kappa_s), c(V, K))
  expect_equal(dim(f$sigma), c(2 * K - 1, 2 * K - 1))
  expect_equal(f$sigma %*% f$sigma_inv, diag(2 * K - 1), tolerance = 1e-8)
  expect_identical(f$vocab, cp$vocab)
  expect_true(all(is.finite(f$elbo)))
  expect_equal(nrow(f$diagnostics), length(f$elbo))
  expect_equal(f$diagnostics$elbo, f$elbo)
  expect_equal(f$elbo[f$bestIter], max(f$elbo))
})

test_that("fits are reproducible and identical for any number of threads", {
  a <- fit_gad(cores = 1L)
  b <- fit_gad(cores = 1L)
  c2 <- fit_gad(cores = 2L)
  expect_identical(a$elbo, b$elbo)
  expect_identical(a$alpha, b$alpha)
  expect_identical(a$elbo, c2$elbo)
  expect_identical(a$alpha, c2$alpha)
  expect_identical(a$kappa, c2$kappa)
})

## One iteration isolates the interpolation error. Over several iterations
## AIC selection can tip at a near-tie (the selected lambda moves to the
## neighbouring grid point), which turns a ~1e-12 difference into ~1e-3 in
## kappa, so the longer comparison is loose.
test_that("exact and interpolated log-partition fits agree", {
  e <- fit_gad(logZ = "exact", maxIter = 1)
  i <- fit_gad(logZ = "interpolated", maxIter = 1)
  expect_equal(i$elbo, e$elbo, tolerance = 1e-10)
  expect_equal(i$alpha, e$alpha, tolerance = 1e-5)
  e <- fit_gad(logZ = "exact", maxIter = 4)
  i <- fit_gad(logZ = "interpolated", maxIter = 4)
  expect_equal(i$elbo, e$elbo, tolerance = 1e-4)
  expect_equal(i$alpha, e$alpha, tolerance = 1e-2)
})

test_that("the best iterate is returned and stopping follows allow.neg.change", {
  f <- fit_gad(maxIter = 60, patience = 2L)
  e <- f$elbo
  expect_equal(f$bestIter, which.max(e))
  if (length(e) < 60) {
    ## stopped early: either converged, or `patience` consecutive decreases
    dec <- diff(e) < 0
    converged <- abs(diff(tail(e, 2))) / abs(tail(e, 1)) < 1e-5
    expect_true(converged || all(tail(dec, 2)))
  }
  ## allow.neg.change = TRUE: stop at the first decrease, return the last iterate
  g <- fit_gad(maxIter = 60, allow.neg.change = TRUE)
  d <- diff(g$elbo)
  expect_true(all(head(d, -1) >= 0))
})

test_that("a fit that stops at maxIter warns and is marked unconverged", {
  cp <- gad_corpus(); set.seed(11)
  w <- character()
  f <- withCallingHandlers(
    sts(~ treatment * pid_rep, ~ noTreatment, cp, K = 3, maxIter = 3, verbose = FALSE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("stopped at maxIter = 3 before converging", w)))
  expect_false(f$converged)
  expect_output(print(f), "did not converge")
  g <- fit_gad(maxIter = 60, convTol = 1e-2)
  expect_true(g$converged)
  expect_lt(length(g$elbo), 60)
})

test_that("argument handling", {
  cp <- gad_corpus()
  ## maxIter = 1 also triggers the stopped-at-maxIter warning; muffle only that
  expect_warning(
    withCallingHandlers(
      sts(~ treatment, ~ noTreatment, cp, K = 3, maxIter = 1, verbose = FALSE,
          kappaEstimation = "adjusted"),
      warning = function(cnd) if (grepl("stopped at maxIter", conditionMessage(cnd)))
        invokeRestart("muffleWarning")),
    "removed in sts 1.7")
  expect_error(
    sts(~ treatment, ~ noTreatment, cp, K = 3, maxIter = 1, verbose = FALSE,
        kappaEstimation = "glm"),
    "must be \"lasso\"")
})

test_that("post-estimation functions run on a fit", {
  f <- gad_fit(); cp <- gad_corpus()
  expect_output(print(f))
  expect_no_error(suppressMessages(capture.output(summary(f))))
  expect_s3_class(plot(f, n = 5), "ggplot")
  ho <- make.heldout(cp$documents, cp$vocab, seed = 1)
  ho$meta <- cp$meta
  set.seed(5)
  fh <- suppressWarnings(sts(~ treatment * pid_rep, ~ noTreatment, ho, K = 3, maxIter = 3,
                             verbose = FALSE))
  hl <- heldoutLikelihood(fh, ho$missing)
  expect_true(is.finite(hl$expected.heldout))
  expect_lt(hl$expected.heldout, 0)
})
