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

opt.kappa <- function(phi, alphaS, kappa, c_d_bar, estimation = estimation, numGroups, group, V, mv, phiD, solver = "cpp", cores = 1L) {
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
  if (estimation == "lasso" && identical(solver, "cpp")) {
    abar <- do.call(rbind, lapply(split(alphaS, group),
                                  function(x) colMeans(matrix(x, ncol = K))))
    res <- sts_optkappa_cpp(phi_list = phi, abar = abar, mv = mv,
                            kappa_t_old = kappa$kappa_t,
                            kappa_s_old = kappa$kappa_s,
                            nlambda = 250L, lambda_min_ratio = 0.001,
                            devmax = 0.999, tol = 1e-10, cores = cores)
    kappa$kappa_t <- res$kappa_t
    kappa$kappa_s <- res$kappa_s
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
