% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heldoutLikelihood.R
\name{heldoutLikelihood}
\alias{heldoutLikelihood}
\title{Compute Heldout Log-Likelihood}
\usage{
heldoutLikelihood(object, missing)
}
\arguments{
\item{object}{an sts object, typically after applying \code{\link[stm]{make.heldout}}}

\item{missing}{list of which words and documents are in the heldout set}
}
\value{
expected.heldout is the average of the held-out log-likelihood values
for each document.
}
\description{
Compute the heldout log-likelihood of the STS model
}
\examples{
\donttest{
library("tm"); library("stm"); library("sts")
temp<-textProcessor(documents=gadarian$open.ended.response,
metadata=gadarian, verbose = FALSE)
out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
out_ho <- make.heldout(out$documents, out$vocab)
out_ho$meta <- out$meta
## low max iteration number just for testing
sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out_ho, K = 3, maxIter = 2, verbose = FALSE)
heldoutLikelihood(sts_estimate, out_ho$missing)$expected.heldout
}
}
