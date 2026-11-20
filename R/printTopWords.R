#' Print Top Words that Load Heavily on each Topic
#' 
#' Prints the top words for each document for low, average, and high levels of sentiment-discourse
#' 
#' @param object Model output from sts
#' @param n number of words to print to console for each topic
#' @param lowerPercentile Percentile to calculate a representative negative sentiment document.
#' @param upperPercentile Percentile to calculate a representative positive sentiment document.
#' @examples
#' \donttest{
#' #Examples with the Gadarian Data
#' library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' printTopWords(sts_estimate)
#' }
#' @export
printTopWords = function(object, n = 10, lowerPercentile = 0.05, upperPercentile = 0.95) {
  .check_sts(object)
  K <- .sts_K(object)
  n <- min(n, length(object$vocab))
  ## top words of each topic at the given sentiment-discourse values
  labels <- function(a) {
    beta <- .sts_beta(object, a, log = TRUE)
    apply(beta, 2, function(x) paste(object$vocab[order(x, decreasing = TRUE)[seq_len(n)]], collapse = ", "))
  }
  cat(sprintf("Topic %i Avg sentiment-discourse: %s \n", seq_len(K), labels(.sts_alpha_s(object))), sep = "")
  cat(sprintf("Topic %i Positive sentiment-discourse: %s \n", seq_len(K),
              labels(.sts_alpha_s(object, "quantile", upperPercentile))), sep = "")
  cat(sprintf("Topic %i Negative sentiment-discourse: %s \n", seq_len(K),
              labels(.sts_alpha_s(object, "quantile", lowerPercentile))), sep = "")
  invisible(object)
}
