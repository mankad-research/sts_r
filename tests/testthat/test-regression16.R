## Regression tests against sts 1.6.
##
## fixtures/ref16_gadarian.rds holds two fits made with sts 1.6 (CRAN) on the
## same corpus and seeds: K = 3 (seed 101) and K = 5 (seed 102), maxIter = 50,
## default settings. sts 1.7 changed two defaults that alter results -- the
## kappa update no longer averages with the previous value, and the best
## iterate is returned -- and its default solvers differ from 1.6's in the last
## digits. The configuration below switches all of that back, so 1.7 must
## reproduce the 1.6 fits up to floating-point reassociation:
##   kappaSolver = "glmnet", logZ = "exact", allow.neg.change = TRUE,
##   estepOptimizer = "bfgs", options(sts.kappaDamping = 0.5).
##
## Exact agreement relies on the same glmnet, BLAS and libm as the machine
## that made the fixtures (Ubuntu 24.04, R 4.3.3, glmnet 4.1-8): AIC selection
## is discontinuous at near-ties, so on another platform a rounding-level
## difference can tip one word's penalty and the fits then drift apart. The
## test therefore runs only when asked for, with STS_REGRESSION_TESTS=true.

test_that("1.7 reproduces sts 1.6 in the 1.6 configuration", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("STS_REGRESSION_TESTS"), "true"),
              "set STS_REGRESSION_TESTS=true to compare with the sts 1.6 fixtures")
  ref <- readRDS(test_path("fixtures", "ref16_gadarian.rds"))
  cp <- gad_corpus()
  old <- options(sts.kappaDamping = 0.5); on.exit(options(old))
  for (nm in names(ref)) {
    r <- ref[[nm]]; K <- ncol(r$kappa_t)
    set.seed(r$seed)
    f <- suppressWarnings(sts(~ treatment * pid_rep, ~ noTreatment, cp, K = K, maxIter = 50,
                              verbose = FALSE, kappaSolver = "glmnet", logZ = "exact",
                              allow.neg.change = TRUE, estepOptimizer = "bfgs", cores = 1L))
    expect_equal(length(f$elbo), length(r$elbo), label = paste(nm, "iterations"))
    expect_equal(f$elbo, r$elbo, tolerance = 1e-8, label = paste(nm, "bound"))
    expect_equal(f$alpha, r$alpha, tolerance = 1e-6, label = paste(nm, "alpha"))
    expect_equal(f$kappa$kappa_s, r$kappa_s, tolerance = 1e-6, label = paste(nm, "kappa_s"))
    expect_equal(unname(f$gamma), r$gamma, tolerance = 1e-6, label = paste(nm, "gamma"))
  }
})
