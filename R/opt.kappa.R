# Kappa Estimation
#
# Lasso-penalized multinomial regression for the kappa coefficients, via the
# distributed Poisson approach (Taddy 2013) on samples aggregated within the
# levels of initializationVar, with the penalty chosen by AIC. The new
# coefficients replace the old ones; sts 1.1-1.6 averaged them with the
# previous values (internal option sts.kappaDamping = 0.5 restores that, for
# regression tests against 1.6).


# glmnet() with maxit and thresh passed the way the installed glmnet expects.
# glmnet >= 5.0 takes them in `control`; passing them as top-level arguments
# is deprecated there and warns (once per session), which the retry loops in
# .fit_glmnet_lasso() and mnreg() would take for a failed fit. The fits are
# the same either way.
.glmnet_fit <- function(..., maxit, thresh) {
  if ("control" %in% names(formals(glmnet::glmnet)))
    glmnet::glmnet(..., control = list(maxit = maxit, thresh = thresh))
  else
    glmnet::glmnet(..., maxit = maxit, thresh = thresh)
}

# Internal helper: fit a Poisson-lasso via glmnet with automatic nlambda retry and AIC
# selection.  Returns the selected coefficient vector (length 2*K).
.fit_glmnet_lasso <- function(x, y, offset_agg, penalty.factor, low, upp) {
  mod <- NULL
  inner_count <- 1
  nlambda <- 250
  while (is.null(mod)) {
    mod <- tryCatch(
      .glmnet_fit(x = x, y = y, family = "poisson",
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
        .glmnet_fit(x = x, y = y, family = "poisson",
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

opt.kappa <- function(phi, alphaS, kappa, estimation = "lasso", numGroups, group, V, mv, solver = "cpp", cores = 1L, lambda_fixed = NULL) {
  low <- -Inf
  upp <- Inf
  K   <- ncol(alphaS)

  ## Dedicated C++ solver for estimation = "lasso" (src/poislasso.cpp).
  ## The design separates into K two-coefficient Poisson lasso problems with
  ## numGroups rows each, coupled only through the shared lambda path and the
  ## joint AIC; the solver reproduces glmnet's lambda_max, grid and dev.ratio
  ## truncation exactly and solves each subproblem to a tighter tolerance than
  ## the thresh = 1e-5 that sts 1.6 passed to glmnet (glmnet's own default is
  ## 1e-7). Where it differs from the glmnet path
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
                            devmax = 0.999, tol = 1e-10, cores = cores,
                            ## internal diagnostic options (not user-facing):
                            ## weight on the old kappa (sts 1.1-1.6: 0.5, the
                            ## half-step (coef + kappa)/2; 1.7: 0), and penalty selection
                            ## (0 = AIC as in 1.6, 1 = least-penalised point)
                            damping = getOption("sts.kappaDamping", 0),
                            select  = if (is.null(lambda_fixed)) as.integer(getOption("sts.kappaSelect", 0L)) else 2L,
                            lambda_fixed = lambda_fixed)
    kappa$kappa_t <- res$kappa_t
    kappa$kappa_s <- res$kappa_s
    attr(kappa, "lambda") <- res$lambda     # selected penalty per word
    return(kappa)
  }
  penalty.factor <- c(rep(1, K), rep(1, K))
  damping <- getOption("sts.kappaDamping", 0)

  phi_sum   <- safelog(unlist(lapply(phi, colSums)))  # numGroups*K vector
  x1_agg    <- matrix(rep(diag(K), numGroups), ncol = K, byrow = TRUE)

  alpha_s   <- lapply(split(alphaS, group), matrix, ncol = K)
  alpha_agg <- do.call(rbind, lapply(alpha_s, function(m) diag(colMeans(m))))


  ## glmnet path (kappaSolver = "glmnet").
  ## Runs serially: the foreach backend was removed in 1.7.
  for (v in 1:V) {
    y_agg      <- unlist(lapply(phi, function(m) m[v, ]))
    offset_agg <- phi_sum + mv[v]

    coef <- .fit_glmnet_lasso(cbind(x1_agg, alpha_agg), y_agg, offset_agg, penalty.factor, low, upp)

    if (damping == 0.5) {          # sts 1.1-1.6 arithmetic, for regression tests
      kappa$kappa_t[v, ] <- (coef[1:K]              + kappa$kappa_t[v, ]) / 2
      kappa$kappa_s[v, ] <- (coef[(K + 1):(2 * K)]  + kappa$kappa_s[v, ]) / 2
    } else {
      kappa$kappa_t[v, ] <- (1 - damping) * coef[1:K]             + damping * kappa$kappa_t[v, ]
      kappa$kappa_s[v, ] <- (1 - damping) * coef[(K + 1):(2 * K)] + damping * kappa$kappa_s[v, ]
    }
  }
  return(kappa)
}
