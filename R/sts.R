#' Variational EM for the Structural Topic and Sentiment-Discourse (STS) Model
#' 
#' 
#' Estimation of the STS Model using variational EM.
#' The function takes sparse representation of a document-term matrix, covariates 
#' for each document, and an integer number of topics and returns fitted model 
#' parameters. See an overview of functions in the package here: 
#' \code{\link{sts-package}}
#' 
#' This is the main function for estimating the Structural Topic and 
#' Sentiment-Discourse (STS) Model. Users provide a corpus of documents and a 
#' number of topics.  Each word in a document comes from exactly one topic and 
#' each document is represented by the proportion of its words that come from 
#' each of the topics. The document-specific content covariates affect how much 
#' (prevalence) and the way in which a topic is discussed (sentiment-discourse). 
#' 
#' @param prevalence_sentiment A formula object with no response variable or a 
#' design matrix with the covariates. The variables must be 
#' contained in corpus$meta.
#' @param initializationVar A formula with a single variable for use in the initialization of latent sentiment. This argument  
#' is usually the key experimental variable (e.g., review rating binary indicator of experiment/control group).
#' @param corpus The document term matrix to be modeled in a sparse term count matrix with one row
#' per document and one column per term. The object must be a list of with each element 
#' corresponding to a document. Each document is represented
#' as an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears. This is the format produced by \code{\link{prepDocuments}} (and used by
#' the stm package). 
#' @param K A positive integer (of size 2 or greater) representing
#' the desired number of topics. 
#' @param maxIter Maximum number of variational EM iterations; default 500 (100
#' before sts 1.7). \code{sts()} warns if it stops here before converging, and
#' the fit's \code{converged} component records whether it converged.
#' @param convTol Convergence tolerance for the variational EM algorithm, a
#' relative change in the bound; default 1e-5. See \code{patience} for how it
#' is used.
#' @param initialization Character argument that allows the user to specify an initialization
#' method. The default choice, \code{"anchor"} to initialize prevalence according to anchor words and 
#' the key experimental covariate identified in argument \code{initializationVar}. One can also use 
#' \code{"stm"}, which uses a fitted STM model (Roberts et al. 2014, 2016) 
#' to initialize coefficients related to prevalence and sentiment-discourse. 
#' The STM model is fit internally (with \code{initializationVar} as the content 
#' covariate), so the stm package does not need to be installed. Its topic 
#' proportions initialize the prevalence variables, and its content-covariate 
#' word effects initialize kappa: for each topic, the sentiment-discourse 
#' coefficients are the slope (and the topic coefficients the intercept) of a 
#' document-weighted least-squares fit of the STM content effects on the centered 
#' levels of \code{initializationVar}, which reproduces the STM topic-word 
#' distributions exactly when \code{initializationVar} has two levels.
#' @param kappaEstimation How kappa is estimated. The only option is
#' \code{"lasso"} (default): for each word, a Poisson regression on the
#' expected word-topic counts aggregated within the levels of
#' \code{initializationVar}, with an L1 penalty whose strength is chosen by AIC
#' along a regularization path. This uses the factorized Poisson approximation
#' to the multinomial of Taddy (2013), Distributed Multinomial Regression; see
#' Chen and Mankad (2025) on the implementation here. The new coefficients
#' replace the old ones at each M-step (sts 1.1 to 1.6 averaged them with the
#' previous values). The \code{"adjusted"} option was removed in sts 1.7; it
#' now gives a warning and uses \code{"lasso"}.
#' @param verbose A logical flag indicating whether information should be
#' printed to the screen.  
#' @param parallelize Deprecated; use \code{cores}. \code{TRUE} is treated
#' as \code{cores = parallel::detectCores() - 1} with a warning.
#' @param stmSeed A prefit STM model object (e.g., from \code{stm::stm()} with 
#' \code{initializationVar} as the content covariate) to initialize the STS model. 
#' Note this is ignored unless initialization = "stm"
#' @param logZ How each topic's log-partition function
#' \eqn{\log Z_k(a) = \log \sum_v \exp(mv_v + \kappa_{t[v,k]} + \kappa_{s[v,k]} a)}
#' is evaluated in the E-step. Within an E-step this is a fixed one-dimensional
#' function of the scalar latent sentiment-discourse value, shared by every
#' document, and it is the only part of the objective, gradient and Hessian that
#' requires a sum over the whole vocabulary. \code{"exact"} sums over all
#' \eqn{V} words on every evaluation, as in sts 1.6. \code{"interpolated"}
#' tabulates \eqn{\log Z_k} together with its first two derivatives once per
#' E-step and uses quintic Hermite interpolation thereafter, which reduces the
#' per-document cost from \eqn{O(VK)} to \eqn{O(n_d K)} for a document with
#' \eqn{n_d} distinct words. Values falling outside the tabulated range are
#' evaluated exactly, so no value is ever extrapolated. The default,
#' \code{"auto"}, decides for each topic at every E-step: it tabulates the
#' topic when that is cheaper than summing over the vocabulary for every
#' document (about 20 evaluations per document, with a factor of 4 to spare)
#' and uses the exact sums otherwise. (sts 1.7.0 decided once, for all
#' topics, at the start.)
#' @param logZSpacing Largest grid step for the interpolated log-partition;
#' default 0.01. Topic \eqn{k} uses the step
#' \eqn{h_k = \min(}\code{logZSpacing}\eqn{, 0.1/S_k)}{h_k = min(logZSpacing, 0.1/S_k)}, where
#' \eqn{S_k = \max_v |\kappa_{s[v,k]}|}{S_k = max |kappa_s[, k]|}: the
#' interpolation error depends on the step through \eqn{h S_k} (it is of order
#' \eqn{(h S_k)^6} in the value and \eqn{S_k^2 (h S_k)^4} in the second
#' derivative), so this keeps it small when a word's coefficient becomes
#' large, as it can on small corpora fitted for many iterations. (sts 1.7.0
#' used \code{logZSpacing} for every topic.) In our tests the error in
#' \eqn{\log Z_k} stayed below about \eqn{10^{-11}} at these steps, for
#' \eqn{S_k} from 1 to 400. Accuracy is not monotone in the step: the
#' interpolated second derivative divides by the square of the step, so very
#' fine grids eventually lose accuracy to rounding.
#' @param logZMargin How far beyond the observed range of latent
#' sentiment-discourse values to tabulate, in the same units; default 2. Larger
#' values trade more grid nodes for fewer exact fallbacks and do not change
#' accuracy.
#' @param allow.neg.change Logical. sts 1.6 and earlier used
#' \code{delta < convTol} as the convergence test, where \code{delta} is the
#' relative change in the bound from one iteration to the next, so any
#' decrease, however large, counted as convergence. The default \code{FALSE}
#' tracks the best bound seen, returns that iterate and stops by the windowed
#' test described under \code{patience}. \code{TRUE} restores the 1.6 rule:
#' stop at the first decrease (or at an increase smaller than \code{convTol})
#' and return the last iterate. (Results from sts 1.6 are not reproduced
#' exactly in any case: 1.7 also changed the kappa update; see
#' \code{kappaEstimation} and NEWS.)
#' @param patience Length of the window for the stopping test, when
#' \code{allow.neg.change = FALSE}; default 20. The fit stops at iteration
#' \eqn{t} when the best bound so far has improved by less than
#' \code{convTol}, relative to its absolute value, over the last
#' \code{patience} iterations, and returns the iterate with the highest
#' bound. The bound is not monotone for this model: the kappa M-step
#' maximizes an l1-penalized Poisson approximation of the bound, with the
#' penalty re-selected by AIC at every iteration, rather than the bound
#' itself, and the bound can fall for one iteration or many. A test on the
#' change over a single iteration then stops wherever one step happens to be
#' small. On the Gadarian data, with three held-out splits and K = 3, 5 and
#' 10, the windowed test with \code{patience = 20} stopped at iterations
#' 64-163, a median 11 below the best bound over 300 iterations; sts 1.7.0's
#' rule (stop when one step changes the bound by less than \code{convTol}
#' or after 5 consecutive decreases) stopped at iterations 6-53 and once
#' returned iteration 1. Changed in sts 1.7.1; see NEWS.
#' @param kappaSolver Solver for the kappa M-step. \code{"cpp"} (default) uses
#' a dedicated Poisson-lasso solver that exploits the separable structure of the
#' aggregated design (one two-coefficient problem per topic, sharing the lambda
#' path and the AIC selection), with a closed-form solution when the
#' initialization variable has two levels. It reproduces glmnet's lambda grid
#' and path truncation and solves each problem to a tighter tolerance than
#' glmnet's default, so on a small fraction of words it selects a different
#' penalty than glmnet did in sts 1.6; in those cases glmnet agrees with it once
#' its convergence threshold is tightened. \code{"glmnet"} uses glmnet as in sts
#' 1.6.
#' @param cores Number of threads for the E-step (parallel over documents), the
#' kappa M-step (parallel over words) and the spectral initialization; default
#' \code{getOption("sts.cores", 1L)}. Results are identical for any number of
#' threads: work is divided into blocks that depend only on the data, and
#' partial sums are combined in a fixed order.
#' @param estepOptimizer Optimizer for each document's variational problem in
#' the E-step. \code{"bfgs"} (default) uses the BFGS algorithm of
#' \code{stats::optim()}, as in sts 1.6 (ported to C++). \code{"newton"} uses
#' Newton's method with the exact Hessian, a capped step, a diagonal shift where
#' the Hessian is not positive definite and a backtracking line search; it needs
#' about 7 objective evaluations per document instead of about 25 and converges
#' to a tighter tolerance. That makes it about three times faster when
#' \code{logZ = "exact"}, where each evaluation sums over the vocabulary, but
#' no faster with the interpolated log-partition, where evaluations are cheap.
#' The per-document objective can have more than one local maximum; the two
#' optimizers reach the same one for over 99.5% of documents in our tests.
#' @return An object of class \code{"STS"}, a list with components
#' 
#' \item{alpha}{Estimated prevalence and sentiment-discourse values for each document and topic} 
#' \item{gamma}{Estimated regression coefficients that determine prevalence and
#' sentiment/discourse for each topic: one row per covariate (with an intercept),
#' and columns for the prevalence of topics 1 to K-1 (relative to topic K) and the
#' sentiment-discourse of topics 1 to K. These are the variational-EM point
#' estimates, also returned by \code{\link{coef.STS}}; for standard errors that
#' account for the uncertainty in the document-level variables, use
#' \code{\link{estimateRegns}}.} 
#' \item{kappa}{Estimated kappa coefficients that determine sentiment-discourse and the topic-word distributions} 
#' \item{sigma_inv}{Inverse of the covariance matrix for the alpha parameters}
#' \item{sigma}{Covariance matrix for the alpha parameters} 
#' \item{elbo}{the approximate ELBO at each iteration of the estimation
#' algorithm: the Laplace variational bound of Chen and Mankad (2025),
#' Proposition 1, summed over documents. sts 1.0 to 1.7.0 reported a different
#' quantity, whose Laplace term had the wrong sign and half the weight on the
#' log-determinant of the negative Hessian; see NEWS for 1.7.1.}
#' \item{converged}{\code{TRUE} if the fit converged, \code{FALSE} if it
#' stopped at \code{maxIter}}
#' \item{mv}{the baseline log-transformed occurrence rate of each word in the corpus}
#' \item{runtime}{Time elapsed in seconds, including initialization (sts 1.6
#' and earlier excluded it)} 
#' \item{vocab}{Vocabulary vector used} 
#' \item{mu}{Mean (fitted) values for alpha based on document-level variables * estimated 
#' Gamma for each document}
#' \item{diagnostics}{A data frame with one row per iteration: the bound
#' (\code{elbo}), its data-fit and Laplace terms (\code{elboObj},
#' \code{elboLap}, where the Laplace term is
#' \eqn{-\frac{1}{2}\log\det(-H_d)-\frac{1}{2}\log\det\Sigma}{-(1/2) log det(-H_d) - (1/2) log det(Sigma)} summed over documents), and the mean relative changes in alpha and kappa
#' (\code{achange}, \code{kchange}). New in sts 1.7.}
#' \item{bestIter}{The iteration whose estimates are returned: the one with
#' the highest bound, or the last with \code{allow.neg.change = TRUE}. New in
#' sts 1.7.}
#' \item{call}{The matched call. New in sts 1.7.}
#' 
#' 
#' @seealso  \code{\link{estimateRegns}}
#' @references 
#' Roberts, M., Stewart, B., Tingley, D., and Airoldi, E. (2013)
#' "The structural topic model and applied social science." In Advances in
#' Neural Information Processing Systems Workshop on Topic Models: Computation,
#' Application, and Evaluation. 
#' 
#' Roberts M., Stewart, B. and Airoldi, E. (2016) "A model of text for
#' experimentation in the social sciences" Journal of the American Statistical
#' Association.
#' 
#' Chen, L. and Mankad, S. (2025) "A Structural Topic and Sentiment-Discourse Model
#' for Text Analysis" Management Science, 71(7), pp.5767-5787.
#' @examples
#' #An example using the Gadarian data.  From raw text to 
#' # fitted model using textProcessor() which leverages the tm Package
#' library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 1, verbose = FALSE)
#' @export
sts = function(prevalence_sentiment, initializationVar, corpus, K, maxIter = 500, convTol = 1e-5, initialization = "anchor", kappaEstimation = "lasso", verbose = TRUE, parallelize = FALSE, stmSeed = NULL, logZ = c("auto", "exact", "interpolated"), logZSpacing = 0.01, logZMargin = 2, allow.neg.change = FALSE, patience = 20L, kappaSolver = c("cpp", "glmnet"), cores = getOption("sts.cores", 1L), estepOptimizer = c("bfgs", "newton")) {
  call <- match.call()
  timer_full <- proc.time()[3]

  # if (!is.null(seed)) {set.seed(seed)}
  if (!(is.character(initialization) && length(initialization) == 1 &&
        initialization %in% c("anchor", "stm"))) {
    stop("initialization must be either \"anchor\" or \"stm\".")
  }
  logZ <- match.arg(logZ)
  if (!(is.numeric(logZSpacing) && length(logZSpacing) == 1 && logZSpacing > 0))
    stop("logZSpacing must be a single positive number.")
  if (!(is.numeric(logZMargin) && length(logZMargin) == 1 && logZMargin >= 0))
    stop("logZMargin must be a single non-negative number.")
  if (!(is.logical(allow.neg.change) && length(allow.neg.change) == 1 &&
        !is.na(allow.neg.change)))
    stop("allow.neg.change must be TRUE or FALSE.")
  if (!(is.numeric(patience) && length(patience) == 1 && patience >= 1))
    stop("patience must be a single integer >= 1.")
  patience <- as.integer(patience)
  kappaSolver <- match.arg(kappaSolver)
  estepOptimizer <- match.arg(estepOptimizer)
  ## Laplace term of the bound.  The default is the bound of Chen and Mankad
  ## (2025), Proposition 1.  Internal option sts.bound = "1.6" restores the
  ## term used by sts 1.0-1.7.0 (wrong sign and half the weight on
  ## log det(-H)); it changes the reported bound and therefore the stopping
  ## point and the returned iterate, and is needed only to reproduce those
  ## releases.
  boundOpt <- getOption("sts.bound", "laplace")
  if (!is.character(boundOpt) || length(boundOpt) != 1L || !boundOpt %in% c("laplace", "1.6"))
    stop("option sts.bound must be \"laplace\" (default) or \"1.6\"", call. = FALSE)
  boundForm <- if (identical(boundOpt, "1.6")) 1L else 0L
  ## Gamma update: stm's variational linear regression (default), or least
  ## squares with the internal option sts.gammaUpdate = "ols" (diagnostics).
  gammaUpdate <- getOption("sts.gammaUpdate", "vb")
  if (!is.character(gammaUpdate) || length(gammaUpdate) != 1L || !gammaUpdate %in% c("vb", "ols"))
    stop('option sts.gammaUpdate must be "vb" (default) or "ols"', call. = FALSE)
  if (!(is.numeric(cores) && length(cores) == 1 && !is.na(cores) && cores >= 1))
    stop("cores must be a single integer >= 1.")
  cores <- as.integer(cores)
  if (isTRUE(parallelize)) {
    ## 1.6 started a PSOCK cluster on detectCores() - 1 workers and ran
    ## foreach() with one task per document and per word; per-task overhead
    ## dominated. 1.7 threads both steps in C++ instead.
    warning("'parallelize' is deprecated; use 'cores' to set the number of ",
            "threads. Using cores = parallel::detectCores() - 1.", call. = FALSE)
    cores <- max(1L, parallel::detectCores() - 1L)
  }
  if(inherits(prevalence_sentiment,"formula")) {
    X <- model.matrix(prevalence_sentiment, data = corpus$meta)[,-1, drop = FALSE]
  } else {
    X <- prevalence_sentiment
  }
  if(inherits(initializationVar,"formula")) {
    X_seed <- model.matrix(initializationVar, data = corpus$meta)[,-1, drop = FALSE]
  } else {
    X_seed <- initializationVar
  }
  


  conv_criteria <- convTol
  topicreportevery <- 10
  if (identical(kappaEstimation, "adjusted")) {
    warning('kappaEstimation = "adjusted" was removed in sts 1.7; using "lasso".', call. = FALSE)
    kappaEstimation <- "lasso"
  }
  if (!identical(kappaEstimation, "lasso"))
    stop('kappaEstimation must be "lasso".', call. = FALSE)
  estimation = kappaEstimation
  nc_X <- ncol(X)
  D <- length(corpus$documents)
  V <- length(corpus$vocab)
  c_d_bar <- unlist(lapply(corpus$documents, function(x) sum(x[2,])))
  wcounts <- tabulate(
    unlist(lapply(corpus$documents, function(x) rep(x[1,], times=x[2,])), use.names=FALSE),
    nbins = V)
  mv <- log(wcounts / sum(wcounts))
  
  if (verbose) {cat("\nStarting Initialization...")}
  
  
  if (verbose) {cat("Beta...")}
  if (initialization == "anchor") {
    ## anchor words (spectral) initialization; see stmInit() in R/stmPort.R
    mod.out <- suppressWarnings(stmInit(corpus$documents, corpus$vocab, K, X = X,
                                        max.em.its = 0, cores = cores))
  } else if (initialization == "stm") {
    ## stm initialization: fit an STM with initializationVar as the content
    ## covariate (internal implementation; the stm package is not required)
    if (is.null(stmSeed))
      mod.out <- suppressWarnings(stmInit(corpus$documents, corpus$vocab, K, X = X,
                                          content = X_seed, max.em.its = 100, cores = cores))
    else
      mod.out <- stmSeed
  }
  
  if (verbose) {cat("Sigma, mu, and alpha...")}
  ## define the other parameters
  Sigma_Inv.est <- diag(1/20, 2*K-1)
  Sigma_Inv.est[1:(K-1),1:(K-1)] <- mod.out$invsigma
  sigma <- solve(Sigma_Inv.est)
  alpha.est <- matrix(0, D, 2*K-1)
  alpha.est[,1:(K-1)] <- mod.out$eta
  alpha.est[,1:K+K-1] <- matrix(X_seed - mean(X_seed), byrow = FALSE)
  # alpha.est[,1:K+K-1] <- matrix(1, byrow = FALSE)
  mu <- opt.mu(lambda=alpha.est, covar=cbind(1,X), enet=0, ic.k=NULL, maxits=1000, method=gammaUpdate)
  Gamma.est <- mu$gamma
  
  ## start sequence to get kappa
  if (verbose) {cat("Kappa...")}
  rho <- 0    
  kappa.est <- list(kappa_t = matrix(0, V, K), kappa_s = matrix(0, V, K))
  
  ## create groups for aggregation (to eventually estimate kappa)
  X_seed_orig <- X_seed
  group <- rep(1, D)
  vals <- sort(unique(X_seed))
  numGroups <- length(vals)
  for (g in 1:numGroups) {
    group[X_seed == vals[g]] <- g
  } 
  
  ## calculate phi
  phi <- vector("list", length = numGroups)
  phiD <- matrix(0, nrow=D, ncol = K)
  for (g in 1:numGroups) {
    phi[[g]] <- matrix(data = 0, nrow = V, ncol = K)
  }  
  BZ <- t(exp(mod.out$beta$logbeta[[1]]))
  BZ <- BZ / rowSums(BZ)
  for (d in 1:D) {
    cdv_ind <- corpus$documents[[d]][1,]
    cdv <- corpus$documents[[d]][2,]
    phi[[group[d]]][cdv_ind,] <- phi[[group[d]]][cdv_ind,] + cdv*BZ[cdv_ind,]
    phiD[d,] <- phiD[d,] + colSums(cdv*BZ[cdv_ind,, drop = FALSE])
  }

  ## estimate kappa 
  if (initialization == "anchor") {
    kappa.est <- opt.kappa(phi = phi, alphaS = alpha.est[,1:K+K-1], kappa = kappa.est, estimation = estimation, numGroups = numGroups, group = group, V = V, mv = mv, solver = kappaSolver, cores = cores)
    kappa.est$kappa_t[which(kappa.est$kappa_t == 0)] <- sign(rnorm(sum(kappa.est$kappa_t == 0)))*1e-4
    kappa.est$kappa_s[which(kappa.est$kappa_s == 0)] <- sign(rnorm(sum(kappa.est$kappa_s == 0)))*1e-4
  } else {
    ## initialize kappa from the content-covariate STM fit (see kappaFromSTM)
    kappa.est <- kappaFromSTM(mod.out, K = K, V = V, vals = vals, group = group)
  }
  
  elbo_vec <- achange_vec <- kchange_vec <- rep(NA, maxIter)
  ## Bound decomposition per iteration: the data-fit term (sum over documents
  ## of the log-posterior at its mode) and the Laplace term.  1.6 computed the
  ## total only; keeping the parts makes a declining bound diagnosable.
  elboObj_vec <- elboLap_vec <- rep(NA, maxIter)
  trace_mu <- trace_kap <- matrix(NA_real_, maxIter, 3)
  kappaLambda <- NULL
  ## best-bound tracking for the windowed stopping test and the returned
  ## iterate (see the convergence block below)
  best_elbo  <- -Inf
  best_iter  <- 1L
  best_state <- NULL
  best_path  <- rep(NA_real_, maxIter)
  ## internal option sts.iterHook: function(iter, state) called after each
  ## iteration with the current estimates as an "STS" object (for example to
  ## record held-out likelihood along a fit); its values are returned in
  ## hookTrace
  iterHook <- getOption("sts.iterHook", NULL)
  if (!is.null(iterHook) && !is.function(iterHook))
    stop("option sts.iterHook must be a function(iter, state)", call. = FALSE)
  hook_trace <- list()
  ntokens <- sum(mod.out$settings$dim$wcounts$x)
  converge <- FALSE
  iter <- 1
  ## ------------------- log-partition backend selection -------------------
  ## Both backends cost O(V*K) per evaluation, so the trade-off is purely the
  ## number of evaluations: the exact path sweeps V once per topic on every
  ## objective evaluation (roughly 20 per document per E-step in practice),
  ## while the table sweeps V once per topic per grid node and then
  ## interpolates.  The table therefore pays off when 20*D comfortably
  ## exceeds the node count -- V and K cancel.  Each topic's grid step is
  ## min(logZSpacing, 0.1 / S_k), where S_k = max|kappa_s[, k]|, which keeps
  ## the interpolation error at the same relative size whatever the size of
  ## kappa_s (sts 1.7.0 used logZSpacing for every topic, and lost accuracy
  ## once a topic's S_k reached the hundreds).  "auto" decides per topic at
  ## every E-step: a topic is tabulated when 20*D exceeds 4 times its node
  ## count and evaluated exactly otherwise.  (1.7.0 decided once, for all
  ## topics, from the initial range of alpha.)
  logZMaxhS <- getOption("sts.logZMaxhS", 0.1)
  if (!(is.numeric(logZMaxhS) && length(logZMaxhS) == 1 && !is.na(logZMaxhS)))
    stop("option sts.logZMaxhS must be a single number", call. = FALSE)
  logZmode_num <- switch(logZ, exact = 0L, interpolated = 1L, auto = 2L)
  if (verbose) {
    alpha_s_rng <- range(alpha.est[, 1:K + K - 1])
    logZnodes_est <- ceiling((diff(alpha_s_rng) + 2 * logZMargin) / logZSpacing) + 1
    cat(switch(logZ,
      exact = "Log-partition: exact.\n",
      interpolated = sprintf("Log-partition: interpolated (h <= %g, ~%d nodes per topic at the start).\n",
                             logZSpacing, logZnodes_est),
      auto = if (20 * D > 4 * logZnodes_est)
        sprintf("Log-partition: auto (tables with h <= %g at the start, ~%d nodes per topic; decided per topic at each E-step).\n",
                logZSpacing, logZnodes_est)
      else "Log-partition: auto (exact at the start: corpus too small for the table to pay off; decided per topic at each E-step).\n"))
  }

  if (verbose) {cat("Done.\nStarting Estimation...\n")}

  while (iter <= maxIter & !converge) {
    t1 <- proc.time()[3]
    alpha.est.old <- alpha.est ## DEBUG
    kappa.est.old <- kappa.est ## DEBUG
    sigma.ss <- diag(0, nrow = 2*K-1)
    elbo <- 0

    ## ----------------------------- E-step -----------------------------
    ## The entire per-document loop runs in C++ (src/estep.cpp): BFGS via a
    ## port of R's vmmin, the Hessian, the Cholesky fix-up and the phi
    ## accumulation.  Relative to 1.6 this forms only the document's own
    ## rows, computes beta once for the objective and gradient instead of
    ## twice, builds mv + kappa_t once per E-step instead of once per
    ## objective call, and reuses beta_bar for phi (1.6 recomputed a full
    ## V x K softmax in R for it).  All of these are algebraic
    ## rearrangements: results match 1.6 up to floating-point
    ## reassociation.
    estep <- sts_estep_cpp(
      documents    = corpus$documents,
      alpha_init   = alpha.est,
      kappa_t      = kappa.est$kappa_t,
      kappa_s      = kappa.est$kappa_s,
      mv           = mv,
      mu           = mu$mu,
      Sigma_Inv    = Sigma_Inv.est,
      ## determinant(logarithm = TRUE) rather than 1.6's log(det(sigma)),
      ## which underflows to -Inf for large K.
      logdet_sigma = as.numeric(determinant(sigma, logarithm = TRUE)$modulus),
      group        = as.integer(group),
      numGroups    = as.integer(numGroups),
      maxit        = 500L,
      logZMode     = logZmode_num,
      logZSpacing  = logZSpacing,
      logZMargin   = logZMargin,
      cores        = cores,
      optimizer    = if (estepOptimizer == "newton") 1L else 0L,
      boundForm    = boundForm,
      logZMaxhS    = logZMaxhS)

    alpha.est <- estep$alpha
    sigma.ss  <- estep$sigma_ss
    elbo      <- estep$elbo
    phi       <- estep$phi
    phiD      <- estep$phiD

    if (verbose) {
      cat(sprintf("Completed E-Step (%d seconds). \n",
                  floor((proc.time() - t1)[3])))
    }
    t1 <- proc.time()
    
    
    # Now do the maximum likelihood updates 
    mu <- opt.mu(lambda=alpha.est, covar=cbind(1,X), enet=0, ic.k=NULL, maxits=1000, method=gammaUpdate)
    Gamma.est <- mu$gamma
    
    sigma <- opt.sigma(nu=sigma.ss, lambda=alpha.est,
                       mu=mu$mu, sigprior=0)
    sigma.chol <- chol(sigma)
    Sigma_Inv.est <- chol2inv(sigma.chol)
    
    ## Diagnostic trace (internal option sts.traceBound): the bound at fixed
    ## alpha after each M-step update, so a decline can be attributed to the
    ## mu/sigma update or to the kappa update. maxit = 0 evaluates without
    ## re-optimising alpha.  The log-partition is evaluated as in the E-step
    ## (same logZ mode and spacing; 1.7.0 always evaluated it exactly here,
    ## so with tables the trace compared exact with interpolated bounds).
    traceB <- isTRUE(getOption("sts.traceBound", FALSE))
    boundAt <- function(kap) {
      r <- sts_estep_cpp(documents = corpus$documents, alpha_init = alpha.est,
             kappa_t = kap$kappa_t, kappa_s = kap$kappa_s, mv = mv, mu = mu$mu,
             Sigma_Inv = Sigma_Inv.est,
             logdet_sigma = as.numeric(determinant(sigma, logarithm = TRUE)$modulus),
             group = as.integer(group), numGroups = as.integer(numGroups),
             maxit = 0L, logZMode = logZmode_num, logZSpacing = logZSpacing,
             logZMargin = logZMargin, cores = cores, boundForm = boundForm,
             logZMaxhS = logZMaxhS)
      c(r$elbo, r$elboObj, r$elboLap)
    }
    if (traceB) trace_mu[iter, ] <- boundAt(kappa.est)

    freezeAfter <- getOption("sts.kappaFreezeAfter", NULL)
    lamFix <- if (!is.null(freezeAfter) && iter > freezeAfter) kappaLambda else NULL
    kappa.est <- opt.kappa(phi = phi, alphaS = alpha.est[,1:K+K-1], kappa = kappa.est, estimation = estimation, numGroups = numGroups, group = group, V = V, mv = mv, solver = kappaSolver, cores = cores, lambda_fixed = lamFix)
    if (is.null(lamFix)) kappaLambda <- attr(kappa.est, "lambda")

    if (traceB) trace_kap[iter, ] <- boundAt(kappa.est)
    
    if(verbose) {
      #M-step message
      timer <- floor((proc.time()-t1)[3])
      msg <- ifelse(timer>1,
                    sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                    "Completed M-Step. \n")
      
      cat(msg)
    }
    
    elbo_vec[iter] <- elbo
    best_path[iter] <- if (iter == 1) elbo else max(best_path[iter - 1], elbo)
    elboObj_vec[iter] <- estep$elboObj
    elboLap_vec[iter] <- estep$elboLap
    if (iter == 1) {
      msg <- sprintf("Completing Iteration %i (approx. per word bound = %.3f) \n", 
                     iter, elbo/ntokens) 
      if(verbose) {
        cat(msg)
      }
    } else {
      new <- elbo_vec[iter]
      old <- elbo_vec[iter-1]
      delta <- (new - old)/(0.1 + abs(old))

      ## Convergence rule.  Up to sts 1.6 this was `delta < conv_criteria`,
      ## which treats ANY decrease in the bound as convergence, however large.
      ## allow.neg.change = TRUE keeps that rule (and returns the last iterate).
      ##
      ## The default rule allows for a bound that is not monotone. The kappa
      ## M-step maximizes an l1-penalized Poisson (DMR) approximation of the
      ## bound with the penalty re-selected by AIC, not the bound itself, so
      ## the bound rises and falls. With the bound of Chen and Mankad (2025,
      ## Prop. 1; see NEWS 1.7.1) and the 1.7 defaults, Gadarian fits (three
      ## held-out splits, K = 3, 5, 10, 300 iterations) decreased in 77-132
      ## of 299 iterations, starting at iteration 2-4, with decreases about as
      ## large as increases. On Yelp (K = 10) the bound peaked at iteration 62
      ## and then fell at every iteration to 110; from the peak on, every kappa
      ## update lowered it.
      ##
      ## A single-step test (|delta| < convTol, or `patience` consecutive
      ## decreases, as in 1.7.0) therefore stops wherever one step happens to
      ## be small: at iterations 6-53 on those Gadarian fits, once returning
      ## iteration 1. The test below is on the best bound so far instead: stop
      ## when it has improved by less than convTol (relative) over the last
      ## `patience` iterations. With patience = 20 the same fits stopped at
      ## iterations 64-163, a median 11 below the 300-iteration maximum,
      ## and Yelp returned its maximum (iteration 62) at iteration 78.
      ## The best-bound state is tracked and returned.
      if (allow.neg.change) {
        converge <- delta < conv_criteria
      } else if (iter > patience) {
        old_best <- best_path[iter - patience]
        converge <- (best_path[iter] - old_best) / (0.1 + abs(old_best)) < conv_criteria
      }
      
      achange <- mean(abs(alpha.est - alpha.est.old)/(0.1 + abs(alpha.est.old)))
      kchange <- mean(abs(c(kappa.est$kappa_t, kappa.est$kappa_s) - c(kappa.est.old$kappa_t, kappa.est.old$kappa_s))/(0.1 + abs(c(kappa.est.old$kappa_t, kappa.est.old$kappa_s))))
      
      achange_vec[iter] <- achange
      kchange_vec[iter] <- kchange
      
      if(verbose) {
        msg <- sprintf("Completing Iteration %i (approx. per word bound = %.3f, relative change = %.3e;\n\t\t\tavg rel change in alpha = %.3e; avg rel change in kappa = %.3e) \n",
                       iter, new/ntokens, delta, achange, kchange)
        cat(msg)
      }
      if (verbose && iter %% topicreportevery==0) {
        ## printTopWords() validates its input (1.7), so give the current
        ## estimates the class of a fitted model
        printTopWords(structure(list(alpha = alpha.est, kappa = kappa.est, sigma = sigma,
                                     mv = mv, vocab = mod.out$vocab), class = "STS"))
      }
    }

    ## Snapshot the best-bound state.  Taken here, after the M-step, so the
    ## snapshot is exactly what the function would have returned had it
    ## stopped at this iteration -- the same point in the cycle 1.6 returns
    ## from, just possibly a different iteration.
    if (!allow.neg.change && elbo > best_elbo) {
      best_elbo <- elbo
      best_iter <- iter
      best_state <- list(alpha = alpha.est, gamma = Gamma.est,
                         kappa = kappa.est, sigma_inv = Sigma_Inv.est,
                         sigma = sigma, mu = mu$mu)
    }
    if (!is.null(iterHook)) {
      hook_trace[[iter]] <- iterHook(iter, structure(list(alpha = alpha.est, gamma = Gamma.est,
        kappa = kappa.est, sigma_inv = Sigma_Inv.est, sigma = sigma, mu = mu$mu, mv = mv,
        vocab = mod.out$vocab, elbo = elbo_vec[seq_len(iter)]), class = "STS"))
    }
    iter <- iter + 1
  }

  elbo_path <- elbo_vec[1:(iter-1)]
  if (!converge) {
    n <- length(elbo_path)
    last_change <- if (n >= 2) (elbo_path[n] - elbo_path[n - 1]) / (0.1 + abs(elbo_path[n - 1])) else NA_real_
    warning(sprintf(paste0("sts() stopped at maxIter = %d before converging ",
                           "(last relative change in the bound %.3e; convTol = %g). ",
                           "Consider a larger maxIter."),
                    maxIter, last_change, convTol), call. = FALSE)
  }
  if (!allow.neg.change && !is.null(best_state) && best_iter < (iter - 1)) {
    if (verbose) {
      cat(sprintf(paste0("Returning the best iterate (%d of %d; bound %.4f ",
                         "vs %.4f at the last iteration).\n"),
                  best_iter, iter - 1, best_elbo, elbo_vec[iter - 1]))
    }
    alpha.est     <- best_state$alpha
    Gamma.est     <- best_state$gamma
    kappa.est     <- best_state$kappa
    Sigma_Inv.est <- best_state$sigma_inv
    sigma         <- best_state$sigma
    mu$mu         <- best_state$mu
  }

  ## label the prevalence/sentiment regression coefficients (1.7)
  dimnames(Gamma.est) <- list(c("(Intercept)", if (is.null(colnames(X))) paste0("X", seq_len(ncol(X))) else colnames(X)),
                              c(sprintf("Topic %d prevalence", seq_len(K - 1)), sprintf("Topic %d sentiment", seq_len(K))))
  final.est <- list(alpha = alpha.est, gamma = Gamma.est, kappa = kappa.est, sigma_inv = Sigma_Inv.est, sigma = sigma, elbo = elbo_path, mv = mv, runtime = proc.time()[3] - timer_full, vocab = mod.out$vocab, mu = mu$mu)
  final.est$converged <- converge
  ## Per-iteration diagnostics (new in 1.7): the bound's two components and
  ## the mean relative change in alpha and kappa.
  final.est$diagnostics <- data.frame(
    iter    = seq_len(iter - 1),
    elbo    = elbo_path,
    elboObj = elboObj_vec[1:(iter - 1)],
    elboLap = elboLap_vec[1:(iter - 1)],
    achange = achange_vec[1:(iter - 1)],
    kchange = kchange_vec[1:(iter - 1)])
  if (isTRUE(getOption("sts.traceBound", FALSE))) {
    final.est$boundTrace <- data.frame(
      iter = seq_len(iter - 1),
      afterE = elbo_path,
      afterMuSigma = trace_mu[1:(iter - 1), 1], afterKappa = trace_kap[1:(iter - 1), 1],
      objAfterMuSigma = trace_mu[1:(iter - 1), 2], objAfterKappa = trace_kap[1:(iter - 1), 2],
      lapAfterMuSigma = trace_mu[1:(iter - 1), 3], lapAfterKappa = trace_kap[1:(iter - 1), 3])
  }
  if (!is.null(iterHook)) final.est$hookTrace <- hook_trace
  final.est$bestIter <- if (allow.neg.change) length(elbo_path) else best_iter
  final.est$call <- call
  class(final.est) <- "STS"

  return(final.est)
}