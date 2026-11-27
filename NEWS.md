# sts 1.7.1

## Bug fix: the reported bound

* The bound reported in `elbo` (and used by the stopping rule and to choose
  the returned iterate) now matches the approximate ELBO of Chen and Mankad
  (2025, Management Science, Proposition 1): for each document,
  f(alpha_d) - (1/2) log det(-H_d) - (1/2) log det(Sigma), where H_d is the
  Hessian of the document objective at its mode. Every earlier release, and
  the paper's replication code, used f(alpha_d) + (1/4) log det(-H_d) -
  (1/2) log det(Sigma): the wrong sign and half the weight on the
  log-determinant. That quantity is not the Laplace approximation of the
  ELBO: it grows without limit as Sigma shrinks, and it often fell in the
  (Gamma, Sigma) update. The Proposition 1 bound cannot fall in the Sigma
  update. Gamma is updated by stm's variational linear regression, which
  shrinks the coefficients, so the joint update can lower it: it did not in
  the default fits we traced, but did in a diagnostic run with the
  least-penalized kappa, where least squares for Gamma removed every such
  decrease.
* The estimates at any given iteration are unchanged, but the reported
  values, the iteration at which a fit stops and the iterate it returns can
  differ from 1.7.0. `diagnostics$elboLap` is the new Laplace term.
* Internal option `options(sts.bound = "1.6")` restores the old term, for
  reproducing earlier fits. To reproduce an sts 1.6 fit, use
  `kappaSolver = "glmnet"`, `logZ = "exact"`, `allow.neg.change = TRUE`,
  `maxIter = 100` and
  `options(sts.kappaDamping = 0.5, sts.recoverL2 = "expgrad", sts.bound = "1.6")`.

## Stopping rule

* The default stopping test is now on the best bound so far: a fit stops at
  iteration t when the best bound has improved by less than `convTol`
  (relative) over the last `patience` iterations, and it returns the iterate
  with the highest bound. The default `patience` is now 20 (5 in 1.7.0).
  1.7.0 stopped when a single iteration changed the bound by less than
  `convTol`, or after `patience` consecutive decreases.
* Why: the bound is not monotone. The kappa M-step maximizes an l1-penalized
  Poisson approximation of the bound, with the penalty re-selected by AIC,
  rather than the bound itself. With the corrected bound this shows up
  throughout a fit, not only near the end.
  * On the Gadarian data (three held-out splits, K = 3, 5, 10, 300
    iterations) the bound decreased in 78-131 of 299 iterations, starting at
    iteration 2-4.
  * Applied to the corrected bound, 1.7.0's rule would have stopped those
    fits at iterations 6-59, one of them at iteration 6 returning iteration
    1. The new rule stops them at iterations 78-173, a median 11 below the
    best bound over 300 iterations.
  * On the Yelp example of the package paper (K = 10, one held-out split)
    the bound peaks at iteration 61 and then falls at every iteration through
    110, the last one run. 1.7.1 stops at iteration 78 and returns iteration
    61. 1.7.0 stopped at iteration 41 and returned iteration 36.
* The per-iteration warning when the bound decreased is gone; a decrease is
  expected. `diagnostics$elbo` records the bound at every iteration.
* `allow.neg.change = TRUE` still restores 1.6's rule: stop at the first
  decrease and return the last iterate.

## Accuracy of the interpolated log-partition functions

* With `logZ = "interpolated"` or `"auto"`, each topic's table now has its
  own step, min(`logZSpacing`, 0.1 / S_k), where S_k is the largest
  |kappa_s| of topic k. The interpolation error depends on the step only
  through h S_k (it is of order (h S_k)^6 in the value), and 1.7.0 used the
  same step for every topic. That was accurate while S_k was of order 1-10,
  as on the Yelp example (S at most 11.5), but S_k can reach the hundreds
  when a rare word's coefficient grows on a small corpus: in five of nine
  Gadarian fits run for 300 iterations S passed 100, up to 612. At 54
  states along those fits the 1.7.0 tables were off by up to 2.4 in the
  bound (1.4e-4 of its value) and 0.03 in alpha; with the new step, by at
  most 2.2e-7 in the bound and 2.1e-7 in alpha, against exact sums.
* Because the fits are sensitive to rounding-level changes at AIC near-ties,
  the more accurate tables also change fits where S_k stays moderate: on the
  Yelp example S reached 11.5 at iteration 20, where the new step is 0.0087
  instead of 0.01, and the fit now returns iteration 61 instead of 62
  (held-out log-likelihood -6.5109 instead of -6.5107).
* `logZ = "auto"` now decides for each topic at every E-step whether to
  tabulate it (when the table is cheaper than summing over the vocabulary
  for every document, as estimated by 1.7.0's rule) or to use the exact
  sums. 1.7.0 decided once, for all topics, from the initial range of alpha.
  On large corpora nothing changes; on small ones a topic with a large S_k
  is evaluated exactly.
* The internal option `options(sts.logZMaxhS = 0)` restores the 1.7.0
  spacing.
* The internal diagnostic `options(sts.traceBound = TRUE)` now evaluates the
  bound after each M-step update with the same log-partition backend as the
  E-step. 1.7.0 always evaluated it exactly there, so with tables the trace
  compared exact with interpolated bounds.

## `estimateRegns()` uses stm's covariance for the draws

* The document-level variables are now drawn with the average of the
  documents' variational covariances, Sigma - (1/D) sum_d (alpha_d -
  mu_d)(alpha_d - mu_d)', which is the "Global" approximation that stm's
  `estimateEffect()` uses. Up to 1.7.0 they were drawn with `sigma`, the
  prior covariance, which also contains the spread of the documents around
  their fitted means; it overstated the uncertainty, so standard errors
  were too large. With the new draws the standard errors of the covariate
  effects are a median 3-9% smaller on the Gadarian data (K = 3, 5) and
  17-18% smaller on the Yelp example of the package paper (K = 10; range
  4-25%); the estimates change only by simulation noise.
  `options(sts.regnsSigma = "prior")` restores the old draws.

## Package

* `esthcpp()`, `lgaecpp()` and `lpbdcpp()` are no longer exported. They
  were the R-callable C++ routines of sts 1.6's E-step (the Hessian, the
  gradient and the objective of one document), which 1.7 replaced; nothing
  in the package calls them, and the tests keep them as reference
  implementations.
* `diagnostics` has a new column, `ncholfix`: the number of documents whose
  negative Hessian at the E-step's optimum was not positive definite and was
  repaired (by 1.6's diagonal fix-up) before inversion. The internal bound
  trace records the same count for its evaluations.
* New tests for text preparation, the `"stm"` initialization,
  `topicSemanticCoherence()` (against a direct computation of Mimno et al.'s
  measure), `findRepresentativeDocs()` and `plotRepresentativeDocs()`, which
  had none. The benchmark script (`benchmarks/`, not in the package build)
  gains time-to-solution, matched-estimator, ablation and large-K parts.
* The vignette has been rewritten for 1.7.1: the model, the workflow on the
  Gadarian data, convergence and the stopping rule, the estimation settings,
  and how to reproduce earlier releases. Its source is in
  `vignettes/source/`.

# sts 1.7

Version 1.7 makes estimation much faster; changes the kappa update, the
anchor-word initialization and the stopping defaults; and fixes the prevalence
tables of `estimateRegns()`. Results differ from sts 1.6; see "Changes that
affect results" below.

## Speed

On the Yelp example of the package paper (35,961 reviews, 3,041 words,
K = 10) one EM iteration takes 5.4 s on one thread and 2.9 s on two, against
about 517 s in sts 1.6 (96x and 181x).

* The per-document E-step runs in C++ (`src/estep.cpp`). The three sums over
  the vocabulary it needs are each topic's log-partition function and its first
  two derivatives, so the cost per document falls from O(V K) to O(n_d K),
  where n_d is the number of distinct words in the document. BFGS is a port of
  R's own `optim()` code, so the optimizer is unchanged.
* New argument `logZ` (`"auto"`, `"exact"`, `"interpolated"`): the
  log-partition functions can be tabulated once per E-step and interpolated
  with quintic Hermite polynomials (`logZSpacing`, default 0.01; `logZMargin`,
  default 2). At the default spacing the error is far below the effect of
  floating-point reassociation; values outside the table are computed exactly.
  `"auto"` uses interpolation when it is cheaper.
* The kappa M-step uses a dedicated Poisson-lasso solver (`src/poislasso.cpp`)
  instead of one `glmnet()` call per word. It reproduces glmnet's lambda grid,
  both of its path-truncation rules and AIC selection, has a closed form when
  `initializationVar` has two levels, and solves each problem to a tighter
  tolerance than glmnet's default. `kappaSolver = "glmnet"` keeps the old path.
* The spectral initialization (anchor words) runs in C++ and needs a fraction
  of the memory.
* New argument `cores`: the E-step, the M-step and the initialization run in
  parallel with RcppParallel. Results are identical for any number of threads.
  The default is `getOption("sts.cores", 1L)`.
* New argument `estepOptimizer = "newton"`: a Newton E-step with the exact
  Hessian, about 3x faster than BFGS with `logZ = "exact"` (no faster with
  interpolation). BFGS remains the default.

## Changes that affect results

* **Convergence.** sts 1.6 stopped at the first decrease of the approximate
  ELBO and returned the last iterate. The bound is not monotone for this model
  (the lasso penalty is re-selected by AIC at every iteration), so that rule
  often stopped early, one iteration past the best point. sts 1.7 returns the
  iterate with the highest bound and stops after `patience` (default 5)
  consecutive decreases or when the relative change falls below `convTol`.
  `allow.neg.change = TRUE` restores the old stopping rule. New components
  `bestIter` and `diagnostics` (the bound's data-fit and Laplace terms and the
  changes in alpha and kappa at each iteration).
* **Defaults `maxIter = 500` (was 100) and `patience = 5` (was 1).** With the
  faster estimation, fits can run to convergence: on the Gadarian corpus 8 of
  9 fits converged within 351 iterations, while at 100 iterations only 1 had.
  Held-out likelihood was essentially unchanged (the extra iterations mainly
  raise the bound), but the estimates are those of a converged fit rather than
  of wherever the cap fell. `sts()` now warns when it stops at `maxIter`, and
  the fit has a `converged` component. With `patience = 1`, a single dip in the
  (non-monotone) bound could end a fit early; one Gadarian fit stopped at
  iteration 2.
* **Kappa update.** sts 1.1 to 1.6 replaced kappa by the average of the new
  estimate and the previous one at each M-step. This half-step is not part of
  the method of Chen and Mankad (2025) and meant the M-step did not maximize
  the bound; it has been removed. Bounds are higher, and held-out likelihood is
  equal or slightly better on the test corpora.
* **`estimateRegns()` prevalence tables (bug fix).** When simulated topic
  proportions were formed, the reference topic's logit was set to 1 instead of
  0, inflating topic K's proportion and biasing all prevalence estimates. Fixed;
  sentiment-discourse tables were not affected. Prevalence results from earlier
  versions should be recomputed.
* **Anchor-word initialization.** The RecoverL2 step of the spectral
  initialization (`initialization = "anchor"`, the default) finds, for each
  word, the weights on the simplex that best reconstruct its co-occurrence row
  from the anchor words. stm and sts 1.1 to 1.6 solved this by exponentiated
  gradient with a fixed step of 50. That iteration is stable only for steps
  below about K / (largest eigenvalue of the anchor rows' Gram matrix), which
  depends on the corpus. On the Gadarian corpus (about 3.5 for K = 3) 456 of
  the 459 non-anchor words never converged, the iterates oscillated until the
  500-iteration cap, and a rounding-level change in the inputs could change
  the weights completely, so the initialization was neither the
  least-squares solution nor the same on different platforms. On larger
  vocabularies the limit is far above 50 (about 170 on the Yelp corpus of
  the package paper and 550-610 on poliblog5k), and the old solver converged
  on every word, but stopped short of the minimum (median relative excess
  1e-6 to 1e-3). sts 1.7 solves each problem to convergence (projected
  gradient with exact projection onto the simplex), reaching at least as
  good a fit for every word and the same initialization on every platform. Fitted models
  differ from 1.6; on the Gadarian corpus, held-out likelihood was better at
  K = 3, similar at K = 5 and slightly worse at K = 10.
  `options(sts.recoverL2 = "expgrad")` restores the old solver.
* **Kappa M-step path truncation.** glmnet stops a lambda path when the
  deviance ratio exceeds 0.999 or when it has improved by less than a
  relative 1e-4 over the last four points. The C++ solver applied only the
  first rule, so on some words it searched further down the path than glmnet
  and AIC could select a point glmnet never reaches. It now applies both.
* The M-step solver and the log-partition interpolation change results in the
  last digits. Because AIC selection is discontinuous at near-ties, such
  differences (and differences between platforms) occasionally change the
  penalty chosen for a word; fitted topics and regression coefficients are
  essentially unaffected.

To reproduce an sts 1.6 fit (up to floating-point reassociation, and on the
same platform), use `kappaSolver = "glmnet"`, `logZ = "exact"`,
`allow.neg.change = TRUE`, `maxIter = 100` and
`options(sts.kappaDamping = 0.5, sts.recoverL2 = "expgrad")`.

## Removed and deprecated

* `kappaEstimation = "adjusted"` has been removed: it fitted worse than
  `"lasso"` on held-out data in every comparison we ran. It now gives a warning
  and uses `"lasso"`.
* `parallelize` is deprecated in favour of `cores`; `parallelize = TRUE` gives a
  warning and uses `cores = parallel::detectCores() - 1`. The foreach and
  doParallel dependencies have been dropped.

## New methods and functions

* `estimateRegns()` returns an object of class `"STSRegns"` (still the same
  list of tables, now named) with `print()`, `coef()`, `vcov()`, `confint()`
  and `as.data.frame()` methods. It is much faster on large corpora. New
  arguments `nsims` and `nsims_se`. Results depend on the random seed and the
  random stream differs from earlier versions.
* `coef()` for fitted models returns the labelled regression coefficients
  (Gamma); `print()` reports the iterations used and the iterate returned. The
  fit stores the matched call, and `runtime` now includes initialization.

## Bug fixes and robustness

* `heldoutLikelihood()`, `printTopWords()`, `plot()`, `topicExclusivity()`,
  `topicSemanticCoherence()`, `findRepresentativeDocs()` and `esthcpp()` no
  longer overflow for large coefficients (they returned `Inf` or `NaN`), and
  `heldoutLikelihood()` is faster. Results are otherwise unchanged.
* The log-determinant of the covariance matrix no longer underflows to `-Inf`
  for large K.
* Functions that take a fitted model check that they were given one, and
  check their other inputs, with informative errors.
* Compatible with glmnet 5.x, which takes `maxit` and `thresh` in `control`.
  Passing them as top-level arguments now gives a deprecation warning, which
  the `kappaSolver = "glmnet"` path and the STM initialization took for a
  failed fit, refitting the first word of a session with a longer lambda path.

## Package

* A testthat suite covering the numerical components against reference
  implementations and glmnet, thread invariance, and the 1.6 reference fits.
* RcppParallel is a new dependency; GNU make is a system requirement.

# sts 1.6

## Bug fix

* `initialization = "stm"` (including `stmSeed`): the sentiment-discourse
  coefficients `kappa_s` were always initialized to zero, because the STM
  topic coefficients were read where the content-covariate interaction
  coefficients were intended. `kappa_t` and `kappa_s` are now initialized
  from the STM's content-covariate effects (intercept and slope of a
  document-weighted least-squares fit on the centered levels of
  `initializationVar`), which reproduces the STM topic-word distributions
  exactly for a binary `initializationVar`. On the Gadarian example this
  gives a higher final ELBO and faster convergence; fits with
  `initialization = "stm"` will differ from sts 1.5. An informative error is
  now given if `stmSeed` has no content covariate matching
  `initializationVar`.

## Removal of the stm dependency

* The package no longer depends on the stm package, which is scheduled for
  archival on CRAN. Apart from the bug fix above, results are numerically
  identical to sts 1.5 for the same random seed.
  - The spectral (anchor word) initialization used by
    `initialization = "anchor"` (the default) and the STM fit used by
    `initialization = "stm"` are now implemented internally (`R/stmPort.R`,
    `src/stmfuns.cpp`), adapted from stm 1.3.8 (MIT license; the stm authors
    are credited in `Authors@R` and `inst/COPYRIGHTS`).
  - New exported functions `textProcessor()`, `prepDocuments()` and
    `make.heldout()`, ported from stm, so that raw text can be prepared for
    `sts()` and held-out data created without stm.
  - New data set `gadarian` (from stm) used in the examples.
  - A prefit STM model can still be supplied through `stmSeed`.
* `sts()` now gives an informative error for an unknown `initialization`
  value.
* Updated the reference to Chen and Mankad (2025), Management Science 71(7).
* The vignette is now registered through the `R.rsp::asis` engine.
* Added SnowballC (used for stemming in `textProcessor()`) to Imports and
  removed stm, knitr and rmarkdown.

# sts 1.5

* Previous CRAN release.
