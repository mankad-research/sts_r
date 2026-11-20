## Post-estimation functions and S3 methods.

## sts 1.6's heldoutLikelihood(), kept verbatim as the reference
heldout16 <- function(object, missing) {
  mv <- object$mv; kappa <- object$kappa; alpha <- object$alpha
  K <- (1 + ncol(alpha)) / 2
  heldout <- vector(length = length(missing$index))
  expeta <- exp(cbind(alpha[, 1:(K - 1)], 0)); theta <- expeta / rowSums(expeta)
  alpha_s <- alpha[, 1:K + K - 1]
  for (i in 1:length(missing$index)) {
    docid <- missing$index[i]; words <- missing$docs[[i]][1, ]
    full_kappa <- exp(mv + kappa$kappa_t + kappa$kappa_s %*% diag(alpha_s[docid, ]))
    beta <- t(apply(full_kappa, 1, function(m) m / colSums(full_kappa)))
    probs <- theta[docid, ] %*% t(beta)[, words]
    probs <- rep(probs, missing$docs[[i]][2, ])
    heldout[i] <- mean(log(probs))
  }
  mean(heldout, na.rm = TRUE)
}

heldout_fit <- function() {
  if (is.null(.sts_test_cache$heldout)) {
    cp <- gad_corpus()
    ho <- make.heldout(cp$documents, cp$vocab, seed = 1); ho$meta <- cp$meta
    set.seed(5)
    f <- suppressWarnings(sts(~ treatment * pid_rep, ~ noTreatment, ho, K = 3, maxIter = 4, verbose = FALSE))
    .sts_test_cache$heldout <- list(fit = f, ho = ho)
  }
  .sts_test_cache$heldout
}

test_that("heldoutLikelihood() matches sts 1.6 and does not overflow", {
  h <- heldout_fit()
  new <- heldoutLikelihood(h$fit, h$ho$missing)
  expect_equal(new$expected.heldout, heldout16(h$fit, h$ho$missing), tolerance = 1e-12)
  expect_length(new$doc.heldout, length(h$ho$missing$index))
  expect_equal(new$ntokens, vapply(h$ho$missing$docs, function(d) sum(d[2, ]), numeric(1)))
  ## coefficients large enough to overflow exp(): 1.6 gives NaN, 1.7 a finite value
  big <- h$fit
  big$mv <- big$mv + 800
  expect_true(is.nan(suppressWarnings(heldout16(big, h$ho$missing))))
  expect_true(is.finite(heldoutLikelihood(big, h$ho$missing)$expected.heldout))
})

test_that("topic-word summaries use a stable normalisation", {
  f <- gad_fit()
  b <- sts:::.sts_beta(f, sts:::.sts_alpha_s(f))
  expect_equal(colSums(b), rep(1, 3), tolerance = 1e-12)
  ## the same distributions as 1.6's exp-then-normalise computation
  a <- colMeans(f$alpha[, 3:5])
  e <- exp(f$mv + f$kappa$kappa_t + f$kappa$kappa_s %*% diag(a))
  expect_equal(b, sweep(e, 2, colSums(e), "/"), tolerance = 1e-12)
  ## and stay finite where that computation overflows
  big <- f; big$mv <- big$mv + 800
  expect_true(all(is.finite(sts:::.sts_beta(big, a))))
  expect_output(printTopWords(big), "Topic 1 Avg")
  expect_true(all(is.finite(topicExclusivity(big))))
})

test_that("coef() on an STS fit returns the labelled Gamma", {
  f <- gad_fit()
  g <- coef(f)
  expect_identical(unname(g), unname(f$gamma))
  expect_identical(rownames(g), c("(Intercept)", "treatment", "pid_rep", "treatment:pid_rep"))
  expect_identical(colnames(coef(f, "prevalence")), c("Topic 1 prevalence", "Topic 2 prevalence"))
  expect_identical(ncol(coef(f, "sentiment")), 3L)
  expect_output(print(f), "estimates are from iteration")
  expect_true(is.call(f$call))
})

test_that("estimateRegns() returns an STSRegns object with working methods", {
  f <- gad_fit(); cp <- gad_corpus()
  set.seed(9)
  r <- estimateRegns(f, ~ treatment * pid_rep, cp, nsims = 5, nsims_se = 50)
  expect_s3_class(r, "STSRegns")
  expect_length(r, 6)
  expect_identical(names(r)[c(1, 4)], c("Topic 1 prevalence", "Topic 1 sentiment"))
  ## still usable as the plain list of tables sts 1.6 returned
  expect_identical(colnames(r[[1]]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  expect_identical(rownames(r[[1]]), c("Intercept", "treatment", "pid_rep", "treatment:pid_rep"))
  expect_output(printRegnTables(r, topics = 2), "Topic 2 sentiment")
  expect_output(print(r), "Topic 3 prevalence")

  cf <- coef(r)
  expect_equal(dim(cf), c(4, 6))
  expect_equal(cf[, "Topic 2 sentiment"], r[["Topic 2 sentiment"]][, "Estimate"])
  expect_equal(dim(coef(r, topic = 1, type = "sentiment")), c(4, 1))

  V <- vcov(r, topic = 2, type = "prevalence")
  expect_equal(dim(V), c(4, 4))
  expect_equal(sqrt(diag(V)), r[["Topic 2 prevalence"]][, "Std. Error"])
  expect_length(vcov(r), 6)

  ci <- confint(r, topic = 1, type = "sentiment", level = 0.9)
  tab <- r[["Topic 1 sentiment"]]
  q <- qt(0.95, nrow(f$alpha) - 4)
  expect_equal(unname(ci[, 1]), unname(tab[, 1] - q * tab[, 2]))
  expect_identical(colnames(ci), c("5 %", "95 %"))
  expect_equal(dim(confint(r, parm = "treatment", topic = 1, type = "prevalence")), c(1, 2))
  expect_length(confint(r), 6)

  d <- as.data.frame(r)
  expect_equal(nrow(d), 24)
  expect_setequal(unique(d$type), c("prevalence", "sentiment"))
  expect_equal(d$estimate[d$topic == 3 & d$type == "sentiment"], unname(cf[, "Topic 3 sentiment"]))

  ## prevalence draws are topic proportions: intercepts sum to one across topics
  expect_equal(sum(cf["Intercept", 1:3]), 1, tolerance = 0.02)
})

test_that("functions reject objects that are not STS fits", {
  cp <- gad_corpus()
  expect_error(printTopWords(list(a = 1)), "must be a model fitted by sts")
  expect_error(heldoutLikelihood(cp, list()), "must be a model fitted by sts")
  expect_error(estimateRegns(cp, ~ treatment, cp), "must be a model fitted by sts")
  expect_error(topicExclusivity(1:3), "must be a model fitted by sts")
  f <- gad_fit()
  expect_error(estimateRegns(f, ~ treatment, list(meta = cp$meta[1:10, ])), "rows but the model has")
  expect_error(findRepresentativeDocs(f, rep("a", 341), topic = 7), "single topic number")
  expect_error(findRepresentativeDocs(f, c("a", "b"), topic = 1), "has 2 elements")
  expect_error(printRegnTables(list(1, 2, 3)), "returned by estimateRegns")
  expect_error(heldoutLikelihood(f, list(index = 1)), "make.heldout")
})
