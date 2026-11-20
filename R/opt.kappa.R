# Kappa Estimation
#
# Penalized multinomial regression methods for estimating the kappa coefficients using the
# distributed Poisson approach.  By default, via the main STS function, estimation is
# performed with lasso.


# Internal helper: fit a Poisson-lasso via glmnet with automatic nlambda retry and AIC
# selection.  Returns the selected coefficient vector (length 2*K).
.fit_glmnet_lasso <- function(x, y, offset_agg, penalty.factor, low, upp) {
  mod <- NULL
  inner_count <- 1
  nlambda <- 250
  while (is.null(mod)) {
    mod <- tryCatch(
      glmnet::glmnet(x = x, y = y, family = "poisson",
                     offset = offset_agg, standardize = FALSE,
                     intercept = FALSE, lambda.min.ratio = 0.001,
                     nlambda = nlambda, alpha = 1,
                     lower.limits = low, upper.limits = upp,
                     maxit = 10000, thresh = 1e-05,
                     penalty.factor = penalty.factor),
      warning = function(w) return(NULL),
      error   = function(e) stop(e))
    # if it didn't converge, increase nlambda paths by 50% and retry
    if (is.null(mod)) nlambda <- nlambda + floor(0.5 * nlambda)
    inner_count <- inner_count + 1
    # after 10 failed attempts, accept estimates as-is
    if (inner_count == 10) {
      mod <- suppressWarnings(
        glmnet::glmnet(x = x, y = y, family = "poisson",
                       offset = offset_agg, standardize = FALSE,
                       intercept = FALSE, lambda.min.ratio = 0.001,
                       nlambda = nlambda, alpha = 1,
                       lower.limits = low, upper.limits = upp,
                       maxit = 10000, thresh = 1e-05,
                       penalty.factor = penalty.factor))
    }
  }
  dev     <- (1 - mod$dev.ratio) * mod$nulldev
  ic      <- dev + 2 * mod$df  ## AIC
  as.numeric(mod$beta[, which.min(ic)])
}

opt.kappa <- function(phi, alphaS, kappa, c_d_bar, estimation = estimation, numGroups, group, V, mv, phiD, solver = "cpp", cores = 1L, lambda_fixed = NULL) {
  low <- -Inf
  upp <- Inf
  K   <- ncol(alphaS)

  ## Dedicated C++ solver for estimation = "lasso" (src/poislasso.cpp).
  ## The design separates into K two-coefficient Poisson lasso problems with
  ## numGroups rows each, coupled only through the shared lambda path and the
  ## joint AIC; the solver reproduces glmnet's lambda_max, grid and dev.ratio
  ## truncation exactly and solves each subproblem to a tighter tolerance than
  ## glmnet's default thresh = 1e-5. Where it differs from the glmnet path
  ## (about 2% of words on poliblog2k) glmnet refit at thresh = 1e-12 agrees
  ## with it, so the differences are glmnet's convergence error, not the
  ## solver's. solver = "glmnet" keeps the 1.6 code path below.
  ## kappaEstimation = "adjusted": Chen and Mankad (2025), Section 4.2 and
  ## Proposition 2, as published -- Poisson regression on the aggregated
  ## samples, with each group's sentiment-discourse value set by equation (13)
  ## and recomputed from the new kappa_s until it converges, starting from
  ## kappa_s = 0; no lasso, no AIC selection, no damping. A ridge of 1e-6 is
  ## the only addition, so that the estimate exists when a word is absent from
  ## a group (see src/adjusted.cpp). solver = "glmnet" keeps the 1.6 variant
  ## (lasso + AIC, at most 10 adjustment passes, half-step damping).
  if (estimation == "adjusted" && identical(solver, "cpp")) {
    res <- sts_optkappa_adjusted_cpp(phi_list = phi, alphaS = alphaS, phiD = phiD,
                                     group = as.integer(group), mv = mv,
                                     ridge    = getOption("sts.adjustedRidge", 1e-6),
                                     tol      = getOption("sts.adjustedTol", 1e-8),
                                     maxInner = as.integer(getOption("sts.adjustedMaxInner", 100L)),
                                     cores = cores)
    damp <- getOption("sts.adjustedDamping", 0)   # internal diagnostic; the paper has none
    kappa$kappa_t <- (1 - damp) * res$kappa_t + damp * kappa$kappa_t
    kappa$kappa_s <- (1 - damp) * res$kappa_s + damp * kappa$kappa_s
    attr(kappa, "adjustedInner") <- res$inner      # Proposition 2 passes per (v, k); -1 = not converged
    return(kappa)
  }
  if (estimation == "lasso" && identical(solver, "cpp")) {
    abar <- do.call(rbind, lapply(split(alphaS, group),
                                  function(x) colMeans(matrix(x, ncol = K))))
    res <- sts_optkappa_cpp(phi_list = phi, abar = abar, mv = mv,
                            kappa_t_old = kappa$kappa_t,
                            kappa_s_old = kappa$kappa_s,
                            nlambda = 250L, lambda_min_ratio = 0.001,
                            devmax = 0.999, tol = 1e-10, cores = cores,
                            ## internal diagnostic options (not user-facing):
                            ## weight on the old kappa, and penalty selection
                            ## (0 = AIC as in 1.6, 1 = least-penalised point)
                            damping = getOption("sts.kappaDamping", 0.5),
                            select  = if (is.null(lambda_fixed)) as.integer(getOption("sts.kappaSelect", 0L)) else 2L,
                            lambda_fixed = lambda_fixed)
    kappa$kappa_t <- res$kappa_t
    kappa$kappa_s <- res$kappa_s
    attr(kappa, "lambda") <- res$lambda     # selected penalty per word
    return(kappa)
  }
  penalty.factor <- c(rep(1, K), rep(1, K))

  phi_sum   <- safelog(unlist(lapply(phi, colSums)))  # numGroups*K vector
  x1_agg    <- matrix(rep(diag(K), numGroups), ncol = K, byrow = TRUE)

  alpha_s   <- lapply(split(alphaS, group), matrix, ncol = K)
  alpha_agg <- do.call(rbind, lapply(alpha_s, function(m) diag(colMeans(m))))

  phiD_list <- lapply(split(phiD, group), matrix, ncol = K)
  G         <- length(alpha_s)

  ## glmnet path (kappaSolver = "glmnet", or kappaEstimation = "adjusted").
  ## Runs serially: the foreach backend was removed in 1.7.
  for (v in 1:V) {
    y_agg      <- unlist(lapply(phi, function(m) m[v, ]))
    offset_agg <- phi_sum + mv[v]

    if (estimation == "lasso") {
      coef <- .fit_glmnet_lasso(cbind(x1_agg, alpha_agg), y_agg, offset_agg, penalty.factor, low, upp)

    } else if (estimation == "adjusted") {
      coef           <- rep(0, 2 * K)
      kk             <- 1
      max_inner_iter <- 10
      weighted_alpha <- alpha_s

      while (kk <= max_inner_iter) {
        coef_old <- coef

        for (t in 1:G) {
          weighted_alpha[[t]] <- alpha_s[[t]] * apply(
            safelog(phiD_list[[t]]) + (alpha_s[[t]] %*% diag(coef[(K + 1):(2 * K)])),
            2, softmax)
        }
        alpha_agg_adj <- do.call(rbind, lapply(weighted_alpha, function(x) diag(colSums(x))))

        coef <- .fit_glmnet_lasso(cbind(x1_agg, alpha_agg_adj), y_agg, offset_agg, penalty.factor, low, upp)

        delta_converge <- mean(abs(coef - coef_old) / (0.1 + abs(coef_old)))
        if (delta_converge <= 1e-4) break
        kk <- kk + 1
      }
    }

    kappa$kappa_t[v, ] <- (coef[1:K]              + kappa$kappa_t[v, ]) / 2
    kappa$kappa_s[v, ] <- (coef[(K + 1):(2 * K)]  + kappa$kappa_s[v, ]) / 2
  }
  return(kappa)
}
