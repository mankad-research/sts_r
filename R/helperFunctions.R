##
#A series of fast softmax functions mostly wrappers around matrixStats package functions -- code from the STM package
##
logsoftmax <- function(x) {
  x - lse(x)
}

lse <- function(x) {
 matrixStats::logSumExp(x)
}

row.lse <- function(mat) {
  matrixStats::rowLogSumExps(mat)
}
col.lse <- function(mat) {
  matrixStats::colLogSumExps(mat)
}

softmax <- function(x) {
  exp(x - lse(x))
}

safelog <- function(x, min=-1000) {
  out <- log(x)
  out[which(out< min)] <- min
  out
}

safeexp = function(x, max=1e8) {
  out <- exp(x)
  out[which(out > max)] <- max
  out
}

## ---- Helpers for STS objects (sts 1.7) -------------------------------------

## Stop with an informative message unless `object` is a fitted STS model.
.check_sts <- function(object, arg = "object") {
  if (!inherits(object, "STS"))
    stop(sprintf("'%s' must be a model fitted by sts() (class \"STS\"), not an object of class \"%s\".",
                 arg, paste(class(object), collapse = "/")), call. = FALSE)
  need <- c("alpha", "kappa", "mv", "sigma")
  miss <- setdiff(need, names(object))
  if (length(miss))
    stop(sprintf("'%s' is missing component(s) %s; was it created by sts()?",
                 arg, paste(miss, collapse = ", ")), call. = FALSE)
  invisible(object)
}

## Number of topics.
.sts_K <- function(object) as.integer((ncol(object$alpha) + 1) / 2)

## Topic proportions theta (D x K), computed with a log-sum-exp so that large
## prevalence values cannot overflow.
.sts_theta <- function(object) {
  K <- .sts_K(object)
  eta <- cbind(object$alpha[, seq_len(K - 1), drop = FALSE], 0)
  exp(eta - matrixStats::rowLogSumExps(eta))
}

## Topic-word distributions (V x K) at sentiment-discourse values `a` (one per
## topic): beta[v, k] proportional to exp(m_v + kappa_t[v, k] + kappa_s[v, k] a_k),
## normalized over the vocabulary with a log-sum-exp. log = TRUE returns
## log(beta).
.sts_beta <- function(object, a, log = FALSE) {
  K <- ncol(object$kappa$kappa_t)
  if (length(a) != K) stop("internal: need one sentiment-discourse value per topic")
  L <- object$mv + object$kappa$kappa_t + sweep(object$kappa$kappa_s, 2, a, `*`)
  L <- sweep(L, 2, matrixStats::colLogSumExps(L), `-`)
  if (log) L else exp(L)
}

## Sentiment-discourse values summarised per topic (mean or a quantile).
.sts_alpha_s <- function(object, stat = "mean", p = NULL) {
  K <- .sts_K(object)
  A <- object$alpha[, seq_len(K) + K - 1, drop = FALSE]
  if (stat == "mean") colMeans(A) else apply(A, 2, stats::quantile, p, names = FALSE)
}
