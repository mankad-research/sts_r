## The spectral (anchor-word) initialization: C++ ports of stm's fastAnchor()
## and recoverL2() must choose the same anchors and the same initial
## topic-word matrix as the R code they replace.

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

## stm's exponentiated-gradient fit (eta = 50) usually does not converge on
## these data: the iterates oscillate until max.its, and rounding-level
## differences (BLAS, compiler, R version) grow by about two orders of
## magnitude every few iterations. So after 500 iterations the C++ and R
## results agree only where the arithmetic is bit-identical. The C++ port is
## compared with expgrad() over the first 10 iterations instead, where
## rounding differences are still below 1e-10.
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
  A_1 <- sts:::recoverL2(s$Q, a, s$wprob, 1L)$A
  A_2 <- sts:::recoverL2(s$Q, a, s$wprob, 2L)$A
  expect_identical(A_1, A_2)
  expect_equal(dim(A_1), c(5L, nrow(s$Q)))
  expect_equal(rowSums(A_1), rep(1, 5), tolerance = 1e-12)
})
