## The log-partition backends. Each topic's normaliser is
##   logZ_k(a) = log sum_v exp(m_v + kappa_t[v,k] + kappa_s[v,k] a),
## whose first two derivatives are the mean and variance of kappa_s[, k] under
## the topic-word distribution. The E-step needs only these three numbers per
## topic, which is what makes its cost O(n_d K) instead of O(V K).

logZ_direct <- function(kt, ks, mv, a) {
  K <- ncol(kt); n <- length(a)
  out <- list(f = matrix(0, n, K), d1 = matrix(0, n, K), d2 = matrix(0, n, K))
  for (k in seq_len(K)) for (i in seq_len(n)) {
    l <- mv + kt[, k] + ks[, k] * a[i]; m <- max(l); p <- exp(l - m); s <- sum(p); p <- p / s
    out$f[i, k]  <- m + log(s)
    out$d1[i, k] <- sum(p * ks[, k])
    out$d2[i, k] <- sum(p * ks[, k]^2) - out$d1[i, k]^2
  }
  out
}

test_that("the exact backend computes logZ and its first two derivatives", {
  f <- gad_fit()
  set.seed(1); a <- runif(50, -2.5, 2.5)
  ex  <- sts:::sts_logZ_probe(f$kappa$kappa_t, f$kappa$kappa_s, f$mv, a, 0L, -3, 3, 0.01)
  ref <- logZ_direct(f$kappa$kappa_t, f$kappa$kappa_s, f$mv, a)
  expect_equal(ex$f,  ref$f,  tolerance = 1e-12)
  expect_equal(ex$d1, ref$d1, tolerance = 1e-12)
  expect_equal(ex$d2, ref$d2, tolerance = 1e-12)
})

test_that("quintic Hermite interpolation is accurate and converges at the quintic rates", {
  f <- gad_fit(); kt <- f$kappa$kappa_t; ks <- f$kappa$kappa_s; mv <- f$mv
  set.seed(1); a <- runif(200, -2.5, 2.5)
  ex <- sts:::sts_logZ_probe(kt, ks, mv, a, 0L, -3, 3, 0.01)
  err <- function(h) {
    hm <- sts:::sts_logZ_probe(kt, ks, mv, a, 1L, -3, 3, h)
    c(f = max(abs(hm$f - ex$f)), d1 = max(abs(hm$d1 - ex$d1)), d2 = max(abs(hm$d2 - ex$d2)))
  }
  ## at the default spacing, relative to the theoretical scale: the quintic
  ## Hermite error is O(h^6 f^(6)) for the value, O(h^5 f^(6)) for d1 and
  ## O(h^4 f^(6)) for d2, and the j-th derivative of logZ in a is a cumulant
  ## of kappa_s, bounded by S^j with S = max|kappa_s|. The constants observed
  ## are 1e-5 to 1e-4 on fits with S from 3 to 6.5.
  h <- 0.01; S <- max(abs(ks))
  e <- err(h)
  expect_lt(e[["f"]]  / (h^6 * S^6), 1e-3)
  expect_lt(e[["d1"]] / (h^5 * S^6), 1e-3)
  expect_lt(e[["d2"]] / (h^4 * S^6), 1e-3)
  ## observed orders, from spacings 0.2 -> 0.05 (theory: 6, 5, 4)
  E <- sapply(c(0.2, 0.1, 0.05), err)
  slope <- function(r) unname(coef(lm(log2(E[r, ]) ~ log2(c(0.2, 0.1, 0.05))))[2])
  expect_gt(slope("f"),  5.5)
  expect_gt(slope("d1"), 4.5)
  expect_gt(slope("d2"), 3.5)
})

test_that("values outside the interpolation range fall back to the exact backend", {
  f <- gad_fit(); kt <- f$kappa$kappa_t; ks <- f$kappa$kappa_s; mv <- f$mv
  a <- c(-6, -3.5, 4, 7)
  hm <- sts:::sts_logZ_probe(kt, ks, mv, a, 1L, -3, 3, 0.01)
  ex <- sts:::sts_logZ_probe(kt, ks, mv, a, 0L, -3, 3, 0.01)
  expect_equal(hm$n_outside, length(a) * ncol(kt))
  expect_identical(hm$f, ex$f)
  expect_identical(hm$d1, ex$d1)
  expect_identical(hm$d2, ex$d2)
})

## sts 1.7.1: each topic's grid step is min(spacing, max_hS / S_k). With one
## spacing for every topic (1.7.0) the error grows like (h S_k)^6, and S_k
## reaches the hundreds on small corpora run for many iterations.
big_kappa <- function(S = 400, at = 0.5) {
  f <- gad_fit(); kt <- f$kappa$kappa_t; ks <- f$kappa$kappa_s; mv <- f$mv
  ## one word whose kappa_s is S, placed so that it starts to dominate topic 1
  ## near a = `at`, which is where logZ_1 changes curvature over a width 1/S
  z <- sts:::sts_logZ_probe(kt, ks, mv, at, 0L, -3, 3, 0.01)$f[1, 1]
  ks[1, 1] <- S; kt[1, 1] <- z - mv[1] - S * at
  list(kt = kt, ks = ks, mv = mv, S = S, fit = f)
}

test_that("the tables stay accurate when one word's kappa_s is in the hundreds", {
  b <- big_kappa()
  set.seed(1); a <- c(seq(0.4, 0.6, length.out = 2001), runif(500, -2.5, 2.5))
  ex  <- sts:::sts_logZ_probe(b$kt, b$ks, b$mv, a, 0L, -3, 3, 0.01)
  old <- sts:::sts_logZ_probe(b$kt, b$ks, b$mv, a, 1L, -3, 3, 0.01)
  new <- sts:::sts_logZ_probe(b$kt, b$ks, b$mv, a, 1L, -3, 3, 0.01, max_hS = 0.1)
  Sk <- apply(abs(b$ks), 2, max)
  expect_equal(as.vector(new$spacing), pmin(0.01, 0.1 / Sk))
  expect_equal(as.vector(old$spacing), rep(0.01, ncol(b$ks)))
  ## errors in logZ, logZ' / S and logZ'' / S^2 for the topic with large S
  err <- function(p) c(max(abs(p$f[, 1] - ex$f[, 1])), max(abs(p$d1[, 1] - ex$d1[, 1])) / b$S,
                       max(abs(p$d2[, 1] - ex$d2[, 1])) / b$S^2)
  expect_gt(err(old)[1], 1e-4)       # the 1.7.0 tables
  e <- err(new)
  expect_lt(e[1], 1e-10)
  expect_lt(e[2], 1e-8)
  expect_lt(e[3], 1e-7)
  ## the other topics are unchanged
  expect_identical(new$f[, -1], old$f[, -1])
})

test_that("auto evaluates a topic exactly when its table would cost more", {
  b <- big_kappa(); cp <- gad_corpus(); f <- b$fit
  f$kappa$kappa_t <- b$kt; f$kappa$kappa_s <- b$ks
  ex <- run_estep(f, cp, logZMode = 0L, maxit = 0L)
  au <- run_estep(f, cp, logZMode = 2L, maxit = 0L)
  tb <- run_estep(f, cp, logZMode = 1L, maxit = 0L)
  old <- run_estep(f, cp, logZMode = 1L, maxit = 0L, logZMaxhS = 0)   # 1.7.0
  ## 341 documents: about 6,800 evaluations, so a table of more than 1,705
  ## nodes is not worth building; topic 1 needs (range / 0.00025) nodes
  expect_equal(au$logZExactTopics, 1)
  expect_true(is.na(au$logZSpacingByTopic[1]))
  expect_equal(tb$logZExactTopics, 0)
  expect_equal(au$elbo, ex$elbo, tolerance = 1e-10)
  expect_equal(tb$elbo, ex$elbo, tolerance = 5e-9)
  ## 160x on Linux, 70x on Windows (the forced tables differ by platform at 1e-5)
  expect_gt(abs(old$elbo - ex$elbo), 20 * abs(tb$elbo - ex$elbo))
})
