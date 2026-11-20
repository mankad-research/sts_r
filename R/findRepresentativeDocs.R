#' Find Documents that Load Heavily on a Topic
#' 
#' Extracts documents with the highest prevalence for a given topic
#' 
#' @param object Model output from sts
#' @param corpus_text vector of text documents, usually contained in the output of prepDocuments
#' @param topic a single topic number
#' @param n number of documents to extract
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
#' docs <- findRepresentativeDocs(sts_estimate, out$meta$open.ended.response, topic = 3, n = 4)
#' plotRepresentativeDocs(docs, text.cex = 0.7, width = 100)
#' }
#' @export
findRepresentativeDocs = function(object, corpus_text, topic, n = 3) {
  .check_sts(object)
  K <- .sts_K(object)
  if (length(topic) != 1 || !(topic %in% seq_len(K)))
    stop(sprintf("'topic' must be a single topic number between 1 and %d.", K), call. = FALSE)
  if (length(corpus_text) != nrow(object$alpha))
    stop(sprintf("'corpus_text' has %d elements but the model has %d documents.",
                 length(corpus_text), nrow(object$alpha)), call. = FALSE)
  theta <- .sts_theta(object)
  n <- min(n, nrow(theta))

  ind <- order(theta[,topic], decreasing = TRUE)[seq_len(n)]
  
  # Calculate quartiles
  quartiles <- quantile(object$alpha[,K + topic - 1], probs = c(0.25, 0.5, 0.75))
  
  # Function to determine the quartile of a value
  get_quartile <- function(value, corpus, quartiles) {
    if (value <= quartiles[1]) {
      return(1)
    } else if (value <= quartiles[2]) {
      return(2)
    } else if (value <= quartiles[3]) {
      return(3)
    } else {
      return(4)
    }
  }
  
  titles <- sapply(object$alpha[ind,K + topic - 1], get_quartile, quartiles = quartiles)
  
  data.frame(text = corpus_text[ind], sentiment_quartile = titles, index = ind)
}  
