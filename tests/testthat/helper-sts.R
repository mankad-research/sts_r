## Shared fixtures for the sts test suite.
##
## gadarian_corpus.rds is the Gadarian and Albertson (2014) corpus as
## prepared by textProcessor() and prepDocuments() (341 documents, 462 words),
## stored so that changes in tm or SnowballC cannot alter the inputs of the
## regression tests. noTreatment = -1 / 1 is the initialization variable.
##
## Thread use is capped at two throughout, as CRAN requires.

gad_corpus <- function() readRDS(test_path("fixtures", "gadarian_corpus.rds"))

## A short fit, cached for the session so the tests that only need some fitted
## state do not each pay for it.
.sts_test_cache <- new.env(parent = emptyenv())
gad_fit <- function(K = 3, maxIter = 5) {
  key <- paste(K, maxIter)
  if (is.null(.sts_test_cache[[key]])) {
    cp <- gad_corpus()
    set.seed(2024)
    .sts_test_cache[[key]] <- suppressWarnings(
      sts(~ treatment * pid_rep, ~ noTreatment, cp, K = K, maxIter = maxIter,
          verbose = FALSE, cores = 1L))
  }
  .sts_test_cache[[key]]
}

## Group indices exactly as sts() builds them from initializationVar.
gad_groups <- function(cp) {
  x <- cp$meta$noTreatment
  vals <- sort(unique(x))
  list(group = as.integer(match(x, vals)), numGroups = length(vals))
}

## Arguments for one E-step (sts_estep_cpp) at the state of a fitted model.
estep_args <- function(fit, cp, logZMode = 0L, cores = 1L, optimizer = 0L, maxit = 500L) {
  g <- gad_groups(cp)
  list(documents = cp$documents, alpha_init = fit$alpha,
       kappa_t = fit$kappa$kappa_t, kappa_s = fit$kappa$kappa_s, mv = fit$mv,
       mu = fit$mu, Sigma_Inv = fit$sigma_inv,
       logdet_sigma = as.numeric(determinant(fit$sigma, logarithm = TRUE)$modulus),
       group = g$group, numGroups = as.integer(g$numGroups), maxit = as.integer(maxit),
       logZMode = as.integer(logZMode), logZSpacing = 0.01, logZMargin = 2,
       cores = as.integer(cores), optimizer = as.integer(optimizer))
}
run_estep <- function(...) do.call(sts:::sts_estep_cpp, estep_args(...))
