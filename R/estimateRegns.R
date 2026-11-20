#' Regression Table Estimation
#' 
#' Estimates regression tables for prevalence and sentiment/discourse.
#' 
#' Estimate Gamma coefficients (along with standard errors, p-values, etc.) to  
#' assess how document-level meta-data determine prevalence and sentiment/discourse  
#'  
#' @param object an sts object
#' @param prevalence_sentiment A formula object with no response variable or a 
#' design matrix with the covariates. If a formula, the variables must be 
#' contained in corpus$meta.
#' @param corpus The document term matrix to be modeled in a sparse term count matrix with one row
#' per document and one column per term. The object must be a list of with each element 
#' corresponding to a document. Each document is represented
#' as an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears. This is the format produced by \code{\link{prepDocuments}}. 
#' 
#' @param nsims Number of draws of the document-level variables (the
#' "method of composition"); default 25, as in earlier versions.
#' @param nsims_se Number of draws from each regression's sampling distribution
#' per document-level draw; default 250.
#'
#' @details The regressions account for uncertainty in the latent document-level
#' variables by the method of composition: \code{nsims} times, the prevalence
#' and sentiment-discourse values of every document are drawn from a normal
#' distribution centred at their estimates with covariance \code{object$sigma},
#' each outcome is regressed on the covariates, and \code{nsims_se} draws are
#' taken from each fitted regression's sampling distribution. Estimates and
#' standard errors are the mean and standard deviation of the pooled draws.
#' Results depend on the random number generator; call \code{set.seed()} first
#' for reproducible tables.
#'
#' @return An object of class \code{"STSRegns"}: a list of \eqn{2K}
#' coefficient tables (columns Estimate, Std. Error, t value, Pr(>|t|)), named
#' \code{"Topic k prevalence"} and \code{"Topic k sentiment"}. The first \eqn{K}
#' are for prevalence (the expected topic proportions) and the last \eqn{K} for
#' sentiment-discourse. It has \code{print}, \code{coef}, \code{vcov},
#' \code{confint} and \code{as.data.frame} methods, and can be indexed like
#' the plain list returned by sts 1.6 and earlier.
#' @seealso \code{\link{printRegnTables}}, \code{\link{coef.STSRegns}}
#'
#' @examples \donttest{
#' library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' regns <- estimateRegns(sts_estimate, ~treatment*pid_rep, out)
#' printRegnTables(x = regns)
#' }
#' @export
estimateRegns = function(object, prevalence_sentiment, corpus, nsims = 25, nsims_se = 250) {
  .check_sts(object)
  alpha <- object$alpha
  sigma <- object$sigma

  if (inherits(prevalence_sentiment, "formula")) {
    if (is.null(corpus$meta)) stop("'corpus' has no $meta data frame for the formula.", call. = FALSE)
    X <- model.matrix(prevalence_sentiment, data = corpus$meta)[, -1, drop = FALSE]
  } else {
    X <- as.matrix(prevalence_sentiment)
  }
  D <- nrow(alpha)
  K <- .sts_K(object)
  if (nrow(X) != D)
    stop(sprintf("The covariates have %d rows but the model has %d documents.", nrow(X), D), call. = FALSE)

  ## Regression on a cached QR decomposition: lighter than lm(), same
  ## estimates and covariance (summary.lm's sigma^2 (X'X)^-1).
  qx <- qr(cbind(1, X))
  p  <- qx$rank
  if (p < ncol(qx$qr)) stop("The covariate matrix is rank deficient.", call. = FALSE)
  XtXinv <- chol2inv(qx$qr[seq_len(p), seq_len(p), drop = FALSE])
  rdf <- D - p
  fit_one <- function(y) {
    beta <- qr.coef(qx, y)
    res  <- qr.resid(qx, y)
    list(est = beta, vcov = (sum(res^2) / rdf) * XtXinv)
  }

  ## Draw all documents' latent variables at once: alpha + Z R with
  ## R'R = sigma. (sts <= 1.6 drew document by document with
  ## mvtnorm::rmvnorm; the distribution is the same, the random stream is not.)
  R <- tryCatch(chol(sigma), error = function(e) {
    ev <- eigen(sigma, symmetric = TRUE)
    t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  })
  P <- ncol(alpha)
  storage <- vector("list", 2 * K)
  for (i in seq_len(nsims)) {
    draws <- alpha + matrix(stats::rnorm(D * P), D, P) %*% R
    eta <- cbind(draws[, seq_len(K - 1), drop = FALSE], 0)
    theta <- exp(eta - matrixStats::rowLogSumExps(eta))
    ysims <- cbind(theta, draws[, seq_len(K) + K - 1, drop = FALSE])
    for (k in seq_len(2 * K)) storage[[k]][[i]] <- fit_one(ysims[, k])
  }

  terms <- c("Intercept", colnames(X))
  if (length(terms) != p) terms <- c("Intercept", paste0("X", seq_len(p - 1)))
  tables <- vector("list", 2 * K)
  vcovs  <- vector("list", 2 * K)
  for (k in seq_len(2 * K)) {
    sims <- do.call(rbind, lapply(storage[[k]], function(x) mvtnorm::rmvnorm(nsims_se, x$est, x$vcov)))
    est  <- colMeans(sims)
    V    <- stats::cov(sims)
    se   <- sqrt(diag(V))
    tval <- est / se
    pval <- 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE)
    tab  <- cbind(est, se, tval, pval)
    dimnames(tab) <- list(terms, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    dimnames(V)   <- list(terms, terms)
    tables[[k]] <- tab
    vcovs[[k]]  <- V
  }
  labs <- c(sprintf("Topic %d prevalence", seq_len(K)), sprintf("Topic %d sentiment", seq_len(K)))
  names(tables) <- names(vcovs) <- labs
  structure(tables, class = "STSRegns", vcov = vcovs, df.residual = rdf, K = K,
            nsims = nsims, nsims_se = nsims_se)
}

## Which tables of an STSRegns object a (topic, type) selection refers to.
.regns_index <- function(object, topic = NULL, type = c("both", "prevalence", "sentiment")) {
  type <- match.arg(type)
  K <- attr(object, "K"); if (is.null(K)) K <- length(object) / 2
  if (is.null(topic)) topic <- seq_len(K)
  if (!all(topic %in% seq_len(K)))
    stop(sprintf("'topic' must be between 1 and %d.", K), call. = FALSE)
  switch(type, prevalence = topic, sentiment = topic + K, both = c(topic, topic + K))
}

#' Methods for Estimated Regression Tables
#'
#' Extract coefficients, their covariance matrices and confidence intervals
#' from the regression tables returned by \code{\link{estimateRegns}}, or
#' arrange them in a data frame.
#'
#' Covariances are those of the pooled simulation draws, so they include the
#' uncertainty in the document-level variables; confidence intervals use the
#' t distribution with the residual degrees of freedom of the regressions, as
#' the p-values in the tables do. Covariances between different tables are not
#' estimated.
#'
#' @param object,x An object of class \code{"STSRegns"} from
#' \code{\link{estimateRegns}}.
#' @param topic Topics to return; default all.
#' @param type \code{"both"} (default), \code{"prevalence"} or
#' \code{"sentiment"}.
#' @param parm Coefficients to return (names or positions); default all.
#' @param level Confidence level.
#' @param row.names,optional Ignored; for compatibility with the generic.
#' @param \dots Further arguments, passed to \code{\link{printRegnTables}} by
#' \code{print}; ignored otherwise.
#'
#' @return \code{coef}: a matrix with one row per coefficient and one column per
#' selected table. \code{vcov} and \code{confint}: a matrix when a single table
#' is selected, otherwise a named list of matrices. \code{as.data.frame}: one
#' row per coefficient and table, with columns \code{topic}, \code{type},
#' \code{term}, \code{estimate}, \code{std.error}, \code{statistic} and
#' \code{p.value}.
#'
#' @examples \donttest{
#' library("sts")
#' temp <- textProcessor(documents = gadarian$open.ended.response,
#'                       metadata = gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' set.seed(1)
#' regns <- estimateRegns(sts_estimate, ~ treatment*pid_rep, out)
#' coef(regns)
#' vcov(regns, topic = 1, type = "sentiment")
#' confint(regns, topic = 1, type = "sentiment")
#' head(as.data.frame(regns))
#' }
#' @name STSRegns-methods
#' @aliases coef.STSRegns vcov.STSRegns confint.STSRegns print.STSRegns as.data.frame.STSRegns
NULL

#' @rdname STSRegns-methods
#' @importFrom stats coef
#' @export
coef.STSRegns <- function(object, topic = NULL, type = c("both", "prevalence", "sentiment"), ...) {
  idx <- .regns_index(object, topic, type)
  out <- vapply(unclass(object)[idx], function(t) t[, "Estimate"], numeric(nrow(object[[1]])))
  if (!is.matrix(out)) out <- matrix(out, ncol = length(idx))
  dimnames(out) <- list(rownames(object[[1]]), names(object)[idx])
  out
}

#' @rdname STSRegns-methods
#' @importFrom stats vcov
#' @export
vcov.STSRegns <- function(object, topic = NULL, type = c("both", "prevalence", "sentiment"), ...) {
  V <- attr(object, "vcov")
  if (is.null(V)) stop("This object has no covariance matrices; re-run estimateRegns() with sts >= 1.7.", call. = FALSE)
  idx <- .regns_index(object, topic, type)
  if (length(idx) == 1) V[[idx]] else V[idx]
}

#' @rdname STSRegns-methods
#' @importFrom stats confint
#' @export
confint.STSRegns <- function(object, parm, level = 0.95, topic = NULL,
                             type = c("both", "prevalence", "sentiment"), ...) {
  idx <- .regns_index(object, topic, type)
  df <- attr(object, "df.residual")
  q <- stats::qt(1 - (1 - level) / 2, df)
  pct <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE, scientific = FALSE, digits = 3), "%")
  all_parm <- missing(parm)
  one <- function(tab) {
    ci <- cbind(tab[, "Estimate"] - q * tab[, "Std. Error"], tab[, "Estimate"] + q * tab[, "Std. Error"])
    dimnames(ci) <- list(rownames(tab), pct)
    if (!all_parm) ci <- ci[parm, , drop = FALSE]
    ci
  }
  out <- lapply(unclass(object)[idx], one)
  if (length(idx) == 1) out[[1]] else out
}

#' @rdname STSRegns-methods
#' @method print STSRegns
#' @export
print.STSRegns <- function(x, ...) {
  printRegnTables(x, ...)
}

#' @rdname STSRegns-methods
#' @export
as.data.frame.STSRegns <- function(x, row.names = NULL, optional = FALSE, ...) {
  K <- attr(x, "K"); if (is.null(K)) K <- length(x) / 2
  rows <- lapply(seq_along(x), function(i) {
    tab <- x[[i]]
    data.frame(topic = (i - 1) %% K + 1,
               type = if (i <= K) "prevalence" else "sentiment",
               term = rownames(tab),
               estimate = tab[, 1], std.error = tab[, 2], statistic = tab[, 3], p.value = tab[, 4],
               row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
