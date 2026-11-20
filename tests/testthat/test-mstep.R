## The kappa M-step solver (src/poislasso.cpp) against glmnet.
##
## For each word, sts 1.6 fit a lasso-penalised Poisson regression with glmnet
## over a 250-point lambda path and picked the penalty by AIC. The C++ solver
## reproduces glmnet's lambda grid and path truncation and solves each point
## to a tighter tolerance; with glmnet's own tolerance tightened (thresh =
## 1e-12) the two must select the same penalty and agree on the coefficients.
## Tests that compare with glmnet do not run on CRAN, so that a change inside
## glmnet cannot fail this package's checks there.

mstep_state <- function() {
  f <- gad_fit(); cp <- gad_corpus(); g <- gad_groups(cp)
  e <- run_estep(f, cp)
  K <- ncol(f$kappa$kappa_t)
  alphaS <- e$alpha[, 1:K + K - 1]
  abar <- do.call(rbind, lapply(split(alphaS, g$group), function(x) colMeans(matrix(x, ncol = K))))
  list(phi = e$phi, abar = abar, alphaS = alphaS, mv = f$mv, kappa = f$kappa, K = K,
       group = g$group, numGroups = g$numGroups)
}

cpp_mstep <- function(s, damping = 0, cores = 1L)
  sts:::sts_optkappa_cpp(s$phi, s$abar, s$mv, s$kappa$kappa_t, s$kappa$kappa_s,
                         nlambda = 250L, lambda_min_ratio = 0.001, devmax = 0.999,
                         tol = 1e-10, cores = cores, damping = damping)

## one word through glmnet, as 1.6's opt.kappa() set it up, at a tight threshold
glmnet_word <- function(s, v, thresh = 1e-12) {
  K <- s$K; G <- s$numGroups
  y  <- unlist(lapply(s$phi, function(m) m[v, ]))
  off <- log(unlist(lapply(s$phi, colSums))) + s$mv[v]
  x1 <- matrix(rep(diag(K), G), ncol = K, byrow = TRUE)
  xa <- do.call(rbind, lapply(seq_len(G), function(g) diag(s$abar[g, ], K)))
  m <- sts:::.glmnet_fit(cbind(x1, xa), y, family = "poisson", offset = off,
                      standardize = FALSE, intercept = FALSE, lambda.min.ratio = 0.001,
                      nlambda = 250, alpha = 1, maxit = 1e6, thresh = thresh)
  ic <- (1 - m$dev.ratio) * m$nulldev + 2 * m$df
  list(coef = as.numeric(m$beta[, which.min(ic)]), fit = m, y = y, off = off)
}

test_that("the C++ solver agrees with glmnet at a tight convergence threshold", {
  skip_on_cran()                              # depends on glmnet internals
  skip_if_not_installed("glmnet")
  s <- mstep_state()
  r <- cpp_mstep(s)
  set.seed(3); words <- sort(sample(nrow(r$kappa_t), 60))
  d <- vapply(words, function(v) {
    cf <- glmnet_word(s, v)$coef
    max(abs(cf - c(r$kappa_t[v, ], r$kappa_s[v, ])))
  }, numeric(1))
  ## AIC is a discontinuous function of the data: allow the odd near-tie
  expect_gte(mean(d < 1e-6), 0.95)
  expect_lt(median(d), 1e-8)
})

test_that("the solver walks glmnet's lambda path", {
  skip_on_cran()                              # depends on glmnet internals
  skip_if_not_installed("glmnet")
  s <- mstep_state()
  v <- which.max(rowSums(s$phi[[1]]))        # the most frequent word
  gw <- glmnet_word(s, v)
  p <- sts:::sts_poislasso_path_cpp(gw$y, gw$off, s$abar, 250L, 0.001, 0.999, 1e-10)
  n <- length(gw$fit$lambda)
  expect_equal(length(p$lambda), n)
  expect_equal(p$lambda, gw$fit$lambda, tolerance = 1e-12)
  expect_equal(p$nulldev, gw$fit$nulldev, tolerance = 1e-10)
  dev_glmnet <- (1 - gw$fit$dev.ratio) * gw$fit$nulldev
  expect_equal(p$dev, dev_glmnet, tolerance = 1e-6)
  expect_equal(p$df, as.numeric(gw$fit$df))
})

test_that("the M-step gives identical results on one and two threads", {
  s <- mstep_state()
  expect_identical(cpp_mstep(s, cores = 1L), cpp_mstep(s, cores = 2L))
})

test_that("the M-step replaces kappa; damping = 0.5 restores the sts 1.6 half-step", {
  s <- mstep_state()
  r0 <- cpp_mstep(s, damping = 0)
  r5 <- cpp_mstep(s, damping = 0.5)
  expect_equal(r5$kappa_t, (r0$kappa_t + s$kappa$kappa_t) / 2, tolerance = 1e-14)
  expect_equal(r5$kappa_s, (r0$kappa_s + s$kappa$kappa_s) / 2, tolerance = 1e-14)
  expect_identical(r0$lambda, r5$lambda)
  ## opt.kappa() uses no damping unless the internal option asks for it
  k <- sts:::opt.kappa(s$phi, s$alphaS, s$kappa, numGroups = s$numGroups, group = s$group,
                       V = length(s$mv), mv = s$mv)
  expect_identical(unclass(k)[c("kappa_t", "kappa_s")], list(kappa_t = r0$kappa_t, kappa_s = r0$kappa_s))
})

test_that("the glmnet and C++ solvers give the same M-step on most words", {
  skip_on_cran()                              # one glmnet fit per word: slow
  skip_if_not_installed("glmnet")
  s <- mstep_state()
  kc <- sts:::opt.kappa(s$phi, s$alphaS, s$kappa, numGroups = s$numGroups, group = s$group,
                        V = length(s$mv), mv = s$mv, solver = "cpp")
  kg <- sts:::opt.kappa(s$phi, s$alphaS, s$kappa, numGroups = s$numGroups, group = s$group,
                        V = length(s$mv), mv = s$mv, solver = "glmnet")
  ## glmnet's default thresh = 1e-5 leaves some words at a different penalty
  d <- apply(abs(cbind(kc$kappa_t - kg$kappa_t, kc$kappa_s - kg$kappa_s)), 1, max)
  expect_gte(mean(d < 1e-3), 0.9)
})
