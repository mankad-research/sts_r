#' Compute Heldout Log-Likelihood
#' 
#' Compute the heldout log-likelihood of the STS model
#' 
#' @param object an sts object, typically estimated on a corpus created with \code{\link{make.heldout}}
#' @param missing list of which words and documents are in the heldout set
#' 
#' @return A list with \code{expected.heldout}, the average over documents of
#' each document's mean held-out log-likelihood per token; \code{doc.heldout},
#' those per-document values; \code{index}, the documents with held-out words;
#' and \code{ntokens}, the number of held-out tokens in each.
#'
#' @details For a document with topic proportions \eqn{\theta_d} and
#' sentiment-discourse values \eqn{\alpha^{(s)}_d}, each held-out token of word
#' \eqn{v} contributes \eqn{\log \sum_k \theta_{d,k} \beta_{d,k,v}}, where
#' \eqn{\beta_{d,k,\cdot}} is topic \eqn{k}'s word distribution at
#' \eqn{\alpha^{(s)}_{d,k}}. Since sts 1.7 the computation is done on the log
#' scale, so it cannot overflow for large coefficients (sts 1.6 returned
#' \code{NaN} there), and is faster.
#'
#' @examples \donttest{
#' library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' out_ho <- make.heldout(out$documents, out$vocab, seed = 1)
#' out_ho$meta <- out$meta
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out_ho, K = 3, maxIter = 2, verbose = FALSE)
#' heldoutLikelihood(sts_estimate, out_ho$missing)$expected.heldout
#' }
#' @export
heldoutLikelihood <- function (object, missing)
{
  .check_sts(object)
  if (!is.list(missing) || is.null(missing$index) || is.null(missing$docs))
    stop("'missing' must be the $missing component returned by make.heldout().", call. = FALSE)
  D <- nrow(object$alpha)
  if (length(missing$index) && max(missing$index) > D)
    stop("'missing' refers to documents the model does not have; was the model fitted to the make.heldout() corpus?",
         call. = FALSE)
  K <- .sts_K(object)
  logtheta <- log(.sts_theta(object))
  alpha_s  <- object$alpha[, seq_len(K) + K - 1, drop = FALSE]
  base     <- object$mv + object$kappa$kappa_t            # V x K

  n <- length(missing$index)
  heldout <- numeric(n)
  ntokens <- numeric(n)
  for (i in seq_len(n)) {
    docid <- missing$index[i]
    words <- missing$docs[[i]][1, ]
    cnt   <- missing$docs[[i]][2, ]
    ## log beta for this document: normalise each topic over the vocabulary
    L <- base + sweep(object$kappa$kappa_s, 2, alpha_s[docid, ], `*`)
    lse <- matrixStats::colLogSumExps(L)
    ## log sum_k theta_k beta_k(w) for each held-out word
    lw <- sweep(L[words, , drop = FALSE], 2, lse - logtheta[docid, ], `-`)
    lp <- matrixStats::rowLogSumExps(lw)
    heldout[i] <- sum(cnt * lp) / sum(cnt)     # mean over tokens, as in sts 1.6
    ntokens[i] <- sum(cnt)
  }
  list(expected.heldout = mean(heldout, na.rm = TRUE),
       doc.heldout = heldout, index = missing$index, ntokens = ntokens)
}
