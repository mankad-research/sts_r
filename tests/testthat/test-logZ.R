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
