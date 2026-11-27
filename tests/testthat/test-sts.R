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

test_that("the best iterate is returned and stopping follows the windowed rule", {
  ## Default rule (1.7.1): stop at iteration t when the best bound so far has
  ## improved by less than convTol (relative) over the last `patience`
  ## iterations; return the iterate with the highest bound.
  win <- function(e, t, P) { b <- cummax(e); (b[t] - b[t - P]) / (0.1 + abs(b[t - P])) }
  f <- fit_gad(maxIter = 60, patience = 5L, convTol = 1e-4)
  e <- f$elbo; n <- length(e)
  expect_equal(f$bestIter, which.max(e))
  expect_equal(f$elbo[f$bestIter], max(f$elbo))
  expect_gt(n, 5)
  if (n < 60) {
    expect_true(f$converged)
    expect_lt(win(e, n, 5), 1e-4)
  }
  ## and the test did not fire at any earlier iteration
  if (n > 6) for (t in 6:(n - 1)) expect_gte(win(e, t, 5), 1e-4)
  ## it cannot fire before iteration patience + 1, and fires there when
  ## convTol is huge
  g <- fit_gad(maxIter = 60, patience = 3L, convTol = 1)
  expect_length(g$elbo, 4)
  expect_true(g$converged)
  ## a single small step no longer stops the fit: 1.7.0 stopped whenever one
  ## step changed the bound by less than convTol
  h <- fit_gad(maxIter = 12, patience = 20L, convTol = 1)
  expect_length(h$elbo, 12)
  expect_false(h$converged)
  ## allow.neg.change = TRUE: stop at the first decrease, return the last iterate
  a <- fit_gad(maxIter = 60, allow.neg.change = TRUE)
  d <- diff(a$elbo)
  expect_true(all(head(d, -1) >= 0))
  expect_equal(a$bestIter, length(a$elbo))
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

test_that("verbose fits print top words every 10 iterations", {
  ## sts 1.7 before this fix stopped with an error at iteration 10 when
  ## verbose = TRUE (the default): the interim printTopWords() call was
  ## given a plain list, which its input check rejects.
  cp <- gad_corpus(); set.seed(11)
  out <- capture.output(f <- suppressWarnings(
    sts(~ treatment * pid_rep, ~ noTreatment, cp, K = 3, maxIter = 10, patience = 100L)))
  expect_s3_class(f, "STS")
  expect_length(f$elbo, 10)
  expect_true(any(grepl("Completing Iteration 10", out)))
  expect_true(any(grepl("^Topic 1 Avg sentiment-discourse:", out)))
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

test_that("internal options: least-squares Gamma and the iteration hook", {
  cp <- gad_corpus()
  old <- options(sts.gammaUpdate = "ols",
                 sts.iterHook = function(iter, st) c(iter, length(st$elbo)))
  on.exit(options(old))
  set.seed(2024)
  f <- sts(~ treatment * pid_rep, ~ noTreatment, cp, K = 3, maxIter = 3,
           allow.neg.change = TRUE, convTol = -1, verbose = FALSE)
  X <- cbind(1, model.matrix(~ treatment * pid_rep, cp$meta)[, -1])
  expect_equal(unname(f$gamma), unname(qr.coef(qr(X), f$alpha)), tolerance = 1e-10)
  expect_equal(f$hookTrace, lapply(1:3, function(i) c(i, i)))
  options(sts.gammaUpdate = "bad")
  expect_error(sts(~ treatment, ~ noTreatment, cp, K = 3, maxIter = 1, verbose = FALSE),
               "sts.gammaUpdate")
})
