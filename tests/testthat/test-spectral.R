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

test_that("recoverL2 in C++ reproduces the R version", {
  s <- spectral_Q()
  a <- sts:::fastAnchor_R(s$Q, 5)
  A_R <- sts:::recoverL2_R(s$Q, a, s$wprob)$A
  A_1 <- sts:::recoverL2(s$Q, a, s$wprob, 1L)$A
  A_2 <- sts:::recoverL2(s$Q, a, s$wprob, 2L)$A
  expect_equal(A_1, A_R, tolerance = 1e-12)
  expect_identical(A_1, A_2)
})
