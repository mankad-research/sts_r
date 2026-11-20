#' Summarize a Fitted STS Model
#' 
#' Function to report on the contents of STS objects
#' 
#' \code{print} reports the size of the model and how the fit ended;
#' \code{summary} also runs \code{\link{printTopWords}}.
#' 
#' @aliases summary.STS print.STS
#' @param object,x An STS object.
#' @param \dots Additional arguments affecting the summary
#' @method summary STS
#' @export
summary.STS <- function(object,...) {
  print.STS(object)
  cat("Most likely words:\n")
  printTopWords(object)
  invisible(object)
}

#' @rdname summary.STS
#' @method print STS
#' @export
print.STS <- function(x,...) {
  .check_sts(x, "x")
  cat(sprintf("A topic model with %i topics, %i documents and a %i word dictionary.\n",
              .sts_K(x), nrow(x$alpha), nrow(x$kappa$kappa_s)))
  if (!is.null(x$elbo)) {
    best <- if (is.null(x$bestIter)) length(x$elbo) else x$bestIter
    cat(sprintf("Fitted in %d iterations; the estimates are from iteration %d (approximate ELBO %.6g).\n",
                length(x$elbo), best, x$elbo[best]))
  }
  invisible(x)
}

#' Coefficients of an STS Model
#'
#' The variational-EM estimates of the regression coefficients (Gamma) that
#' link the document-level covariates to topic prevalence and
#' sentiment-discourse.
#'
#' These are point estimates. For standard errors and tests that account for
#' the uncertainty in the latent document-level variables, use
#' \code{\link{estimateRegns}}, whose prevalence tables are on the scale of the
#' topic proportions rather than the log-ratio scale used here.
#'
#' @param object An STS object.
#' @param type \code{"both"} (default), \code{"prevalence"} or
#' \code{"sentiment"}.
#' @param \dots Ignored.
#'
#' @return A matrix with one row per coefficient (intercept first) and one
#' column per outcome: the prevalence of topics 1 to K-1, as log-ratios against
#' topic K, and the sentiment-discourse of topics 1 to K.
#' @importFrom stats coef
#' @export
coef.STS <- function(object, type = c("both", "prevalence", "sentiment"), ...) {
  .check_sts(object)
  type <- match.arg(type)
  K <- .sts_K(object)
  g <- object$gamma
  if (is.null(colnames(g)))   # fits from sts <= 1.6
    colnames(g) <- c(sprintf("Topic %d prevalence", seq_len(K - 1)), sprintf("Topic %d sentiment", seq_len(K)))
  switch(type, both = g,
         prevalence = g[, seq_len(K - 1), drop = FALSE],
         sentiment  = g[, seq_len(K) + K - 1, drop = FALSE])
}
