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
#' @param maxIter A positive integer representing the max number of VEM iterations allowed.
#' @param convTol Convergence tolerance for the variational EM estimation algorithm; Default value = 1e-5.
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
#' @param kappaEstimation A character input specifying how kappa should be estimated. \code{"lasso"} (default) allows for 
#' penalties on the L1 norm.  We estimate a regularization path and then select the optimal
#' shrinkage parameter using AIC. \code{"adjusted"} utilizes the lasso penalty with an adjusted aggregated Poisson regression. 
#' All options use an approximation framework developed in Taddy (2013) called
#' Distributed Multinomial Regression which utilizes a factorized poisson
#' approximation to the multinomial.  See Chen and Mankad (2025) on the implementation here.
#' @param verbose A logical flag indicating whether information should be
#' printed to the screen.  
#' @param parallelize A logical flag indicating whether to parallelize the estimation using all but one CPU cores on your local machine. 
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
#' \code{"auto"}, uses the table when the corpus is large enough for it to pay
#' off and the exact sum otherwise.
#' @param logZSpacing Grid spacing for \code{logZ = "interpolated"}; default
#' 0.01. Accuracy is not monotone in this value: the interpolated second
#' derivative divides by the square of the spacing, so very fine grids
#' eventually lose accuracy to rounding. 0.01 is near the optimum, and 0.02 is
#' slightly faster with essentially the same typical accuracy.
#' @param logZMargin How far beyond the observed range of latent
#' sentiment-discourse values to tabulate, in the same units; default 2. Larger
#' values trade more grid nodes for fewer exact fallbacks and do not change
#' accuracy.
#' @param allow.neg.change Logical. sts 1.6 and earlier used
#' \code{delta < convTol} as the convergence test, so any decrease in the ELBO,
#' however large, counted as convergence; on the Gadarian data this stopped 3 of
#' 4 test fits. The default \code{FALSE} instead tracks the best bound seen and
#' returns that iterate, stopping once the bound has decreased for
#' \code{patience} consecutive iterations or has stopped changing by more than
#' \code{convTol}. This can only improve on the bound the old rule returned from
#' the same stopping point. \code{TRUE} restores the 1.6 behaviour exactly --
#' stop at the first decrease and return the last iterate -- for reproducing
#' results from sts 1.6 or earlier.
#' @param patience Number of consecutive iterations with a decreasing ELBO
#' tolerated before stopping, when \code{allow.neg.change = FALSE}; default 1.
#' Note that the bound is not guaranteed to increase monotonically for this
#' model: on the Gadarian data it rises to a single interior maximum and then
#' declines for as long as the fit is allowed to run, so large values of
#' \code{patience} cost time without improving the returned fit.
#' @return An object of class sts
#' 
#' \item{alpha}{Estimated prevalence and sentiment-discourse values for each document and topic} 
#' \item{gamma}{Estimated regression coefficients that determine prevalence and sentiment/discourse for each topic} 
#' \item{kappa}{Estimated kappa coefficients that determine sentiment-discourse and the topic-word distributions} 
#' \item{sigma_inv}{Inverse of the covariance matrix for the alpha parameters}
#' \item{sigma}{Covariance matrix for the alpha parameters} 
#' \item{elbo}{the ELBO at each iteration of the estimation algorithm}
#' \item{mv}{the baseline log-transformed occurrence rate of each word in the corpus}
#' \item{runtime}{Time elapsed in seconds} 
#' \item{vocab}{Vocabulary vector used} 
#' \item{mu}{Mean (fitted) values for alpha based on document-level variables * estimated 
#' Gamma for each document}
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
sts = function(prevalence_sentiment, initializationVar, corpus, K, maxIter = 100, convTol = 1e-5, initialization = "anchor", kappaEstimation = "lasso", verbose = TRUE, parallelize = FALSE, stmSeed = NULL, logZ = c("auto", "exact", "interpolated"), logZSpacing = 0.01, logZMargin = 2, allow.neg.change = FALSE, patience = 1L) {

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
  
  if (parallelize) {
    ## Only the kappa M-step still uses the foreach backend; the E-step is
    ## threaded inside C++ (see src/estep.cpp).
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
  }


  conv_criteria <- convTol
  topicreportevery <- 10
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
                                        max.em.its = 0))
  } else if (initialization == "stm") {
    ## stm initialization: fit an STM with initializationVar as the content
    ## covariate (internal implementation; the stm package is not required)
    if (is.null(stmSeed))
      mod.out <- suppressWarnings(stmInit(corpus$documents, corpus$vocab, K, X = X,
                                          content = X_seed, max.em.its = 100))
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
  mu <- opt.mu(lambda=alpha.est, covar=cbind(1,X), enet=0, ic.k=NULL, maxits=1000)
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
    kappa.est <- opt.kappa(phi = phi, alphaS = alpha.est[,1:K+K-1], kappa = kappa.est, c_d_bar = c_d_bar, estimation = estimation, numGroups = numGroups, group = group, V = V, mv = mv, phiD =phiD, parallelize = parallelize)
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
  ## best-bound tracking and the consecutive-decrease counter (see the
  ## convergence block below)
  best_elbo  <- -Inf
  best_iter  <- 1L
  best_state <- NULL
  n_decrease <- 0L
  ntokens <- sum(mod.out$settings$dim$wcounts$x)
  converge <- FALSE
  iter <- 1
  ## ------------------- log-partition backend selection -------------------
  ## Both backends cost O(V*K) per evaluation, so the trade-off is purely the
  ## number of evaluations: the exact path sweeps V once per topic on every
  ## objective evaluation (roughly 20 per document per E-step in practice),
  ## while the table sweeps V once per topic per grid node and then
  ## interpolates.  The table therefore pays off when 20*D comfortably
  ## exceeds the node count, which is a statement about D alone -- V and K
  ## cancel.  "auto" requires a 4x margin before switching.
  alpha_s_rng <- range(alpha.est[, 1:K + K - 1])
  logZnodes_est <- ceiling((diff(alpha_s_rng) + 2 * logZMargin) / logZSpacing) + 1
  evals_est <- 20 * D
  logZmode_num <- switch(logZ,
    exact        = 0L,
    interpolated = 1L,
    auto         = if (evals_est > 4 * logZnodes_est) 1L else 0L)
  if (verbose) {
    if (logZmode_num == 1L) {
      cat(sprintf("Log-partition: interpolated (h = %g, ~%d nodes; est. %dx fewer V-sweeps).\n",
                  logZSpacing, logZnodes_est, round(evals_est / logZnodes_est)))
    } else {
      cat(sprintf("Log-partition: exact%s.\n",
                  if (logZ == "auto") " (corpus too small for the table to pay off)" else ""))
    }
  }

  if (verbose) {cat("Done.\nStarting Estimation...\n")}
  timer_full <- proc.time()[3]

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
      logZMargin   = logZMargin)

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
    mu <- opt.mu(lambda=alpha.est, covar=cbind(1,X), enet=0, ic.k=NULL, maxits=1000)
    Gamma.est <- mu$gamma
    
    sigma <- opt.sigma(nu=sigma.ss, lambda=alpha.est,
                       mu=mu$mu, sigprior=0)
    sigma.chol <- chol(sigma)
    Sigma_Inv.est <- chol2inv(sigma.chol)
    
    kappa.est <- opt.kappa(phi = phi, alphaS = alpha.est[,1:K+K-1], kappa = kappa.est, c_d_bar = c_d_bar, estimation = estimation, numGroups = numGroups, group = group, V = V, mv = mv, phiD =phiD, parallelize = parallelize)
    
    if(verbose) {
      #M-step message
      timer <- floor((proc.time()-t1)[3])
      msg <- ifelse(timer>1,
                    sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                    "Completed M-Step. \n")
      
      cat(msg)
    }
    
    elbo_vec[iter] <- elbo
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
      ## which treats ANY decrease in the ELBO as convergence, however large:
      ## a dip of delta = -7e-5 against convTol = 1e-5 stopped the fit.  That
      ## fired in 3 of 4 Gadarian test runs, so reported iteration counts and
      ## final bounds from 1.6 and earlier partly reflect premature stops
      ## rather than convergence.
      ##
      ## Simply allowing the fit to continue past a decrease is NOT the fix.
      ## On the Gadarian data the bound rises to a single interior maximum and
      ## then declines monotonically for as long as it is allowed to run: K=3
      ## with kappaEstimation = "lasso" peaks at iteration 25 and then falls
      ## for 35 consecutive iterations (-227 by iteration 60), and K=3 with
      ## "adjusted" peaks at iteration 8 and falls for 52 (-1554). So the
      ## algorithm does not merely dip; past the peak it diverges, and 1.6's
      ## rule was accidentally stopping one iteration past the best point.
      ##
      ## What is therefore fixed here is which iterate is RETURNED, not the
      ## stopping test. The best-bound state is tracked and returned, so the
      ## result can only improve on what the same stopping point gave before
      ## (+2.09 and +4.40 on the two fits above). Stopping still happens on a
      ## sustained decrease, after `patience` consecutive decreasing
      ## iterations, which keeps the fit out of the declining region.
      ##
      ## allow.neg.change = TRUE restores 1.6 exactly: stop on the first
      ## decrease and return the last iterate, not the best one.
      if (delta < 0) n_decrease <- n_decrease + 1L else n_decrease <- 0L
      if (delta < 0 && abs(delta) > conv_criteria && !allow.neg.change) {
        warning(sprintf(paste0("The ELBO decreased at iteration %d ",
                               "(relative change %.3e). The best iterate so ",
                               "far (iteration %d) will be returned."),
                        iter, delta, which.max(elbo_vec[1:iter])), call. = FALSE)
      }
      converge <- if (allow.neg.change) {
        delta < conv_criteria
      } else {
        abs(delta) < conv_criteria || n_decrease >= patience
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
        printTopWords(list(alpha = alpha.est, kappa = kappa.est, sigma = sigma, mv = mv, vocab = mod.out$vocab))
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
    iter <- iter + 1
  }

  elbo_path <- elbo_vec[1:(iter-1)]
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

  final.est <- list(alpha = alpha.est, gamma = Gamma.est, kappa = kappa.est, sigma_inv = Sigma_Inv.est, sigma = sigma, elbo = elbo_path, mv = mv, runtime = proc.time()[3] - timer_full, vocab = mod.out$vocab, mu = mu$mu)
  ## Per-iteration diagnostics (new in 1.7): the bound's two components and
  ## the mean relative change in alpha and kappa.
  final.est$diagnostics <- data.frame(
    iter    = seq_len(iter - 1),
    elbo    = elbo_path,
    elboObj = elboObj_vec[1:(iter - 1)],
    elboLap = elboLap_vec[1:(iter - 1)],
    achange = achange_vec[1:(iter - 1)],
    kchange = kchange_vec[1:(iter - 1)])
  final.est$bestIter <- if (allow.neg.change) length(elbo_path) else best_iter
  if (parallelize) {
    # Stop the cluster when done
    parallel::stopCluster(cl)
  }
  class(final.est) <- "STS"

  return(final.est)
}