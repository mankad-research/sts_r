## Text preparation, the "stm" initialization and the post-estimation
## helpers that the other files do not exercise.

test_that("textProcessor() and prepDocuments() reproduce the stored corpus", {
  skip_on_cran()          # depends on tm and SnowballC, which can change
  temp <- textProcessor(documents = gadarian$open.ended.response,
                        metadata = gadarian, verbose = FALSE)
  expect_s3_class(temp, "textProcessor")
  expect_output(print(temp), "documents")
  out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
  cp <- gad_corpus()
  expect_identical(out$vocab, cp$vocab)
  expect_identical(lapply(out$documents, unname), lapply(cp$documents, unname))
  ## words in fewer documents than lower.thresh are dropped
  out5 <- prepDocuments(temp$documents, temp$vocab, temp$meta, lower.thresh = 5,
                        verbose = FALSE)
  df <- tabulate(unlist(lapply(out5$documents, function(d) d[1, ])), length(out5$vocab))
  expect_gt(min(df), 5)
  expect_lt(length(out5$vocab), length(out$vocab))
})

test_that("topicSemanticCoherence() follows Mimno et al. (2011)", {
  f <- gad_fit(); cp <- gad_corpus()
  sc <- topicSemanticCoherence(f, cp, M = 5)
  expect_length(sc, 3)
  expect_true(all(is.finite(sc)))
  ## topic 1 directly: sum over word pairs of log((1 + D(v_m, v_l)) / D(v_l)),
  ## l ranked above m among the topic's top words
  beta <- sts:::.sts_beta(f, sts:::.sts_alpha_s(f))
  top <- order(beta[, 1], decreasing = TRUE)[1:5]
  has <- sapply(top, function(v) vapply(cp$documents, function(d) v %in% d[1, ], logical(1)))
  ref <- 0
  for (m in 2:5) for (l in 1:(m - 1))
    ref <- ref + log((1 + sum(has[, m] & has[, l])) / sum(has[, l]))
  expect_equal(sc[1], ref)
  expect_error(topicSemanticCoherence(f, list(documents = cp$documents, vocab = cp$vocab[-1])),
               "vocabulary")
})

test_that("findRepresentativeDocs() and plotRepresentativeDocs()", {
  f <- gad_fit(); cp <- gad_corpus()
  txt <- sprintf("document %d", seq_along(cp$documents))   # the fixture keeps no text
  docs <- findRepresentativeDocs(f, txt, topic = 2, n = 4)
  theta <- sts:::.sts_theta(f)
  expect_identical(docs$index, order(theta[, 2], decreasing = TRUE)[1:4])
  expect_identical(docs$text, txt[docs$index])
  q <- quantile(f$alpha[, 3 + 2 - 1], c(0.25, 0.5, 0.75))
  expect_identical(docs$sentiment_quartile,
                   as.numeric(findInterval(f$alpha[docs$index, 4], q, left.open = TRUE) + 1))
  expect_error(findRepresentativeDocs(f, txt, topic = 4), "topic")
  expect_error(findRepresentativeDocs(f, txt[-1], topic = 1), "elements")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plotRepresentativeDocs(docs, text.cex = 0.7, width = 80))
})

test_that("the stm initialization runs and starts kappa from the STM fit", {
  skip_on_cran()          # fits an STM first
  cp <- gad_corpus()
  set.seed(2024)
  f <- suppressWarnings(sts(~ treatment * pid_rep, ~ noTreatment, cp, K = 3,
                            initialization = "stm", maxIter = 2, verbose = FALSE))
  expect_s3_class(f, "STS")
  expect_length(f$elbo, 2)
  expect_true(all(is.finite(f$elbo)))
  expect_gt(sum(abs(f$kappa$kappa_s)), 0)
})
