## sts17_benchmark.R -- scaling benchmarks for sts 1.7 (and 1.6 for comparison)
##
## Run on your Mac or the NC State HPC. Every fit runs in its own R process, so
## the two package versions never share a session and timings include nothing
## else. Results are appended to <out>/results.csv as they finish; rerunning
## skips rows already there, so the script can be interrupted and resumed.
##
## ---------------------------------------------------------------------------
## Setup (once)
##
##   # 1.6 from CRAN into its own library
##   mkdir -p ~/Rlib16 ~/Rlib17
##   R -e 'install.packages("sts", lib = "~/Rlib16")'
##
##   # 1.7 from the v1.7-speed branch (clone or bundle), optimized build
##   R CMD INSTALL --preclean --library=$HOME/Rlib17 /path/to/sts_r_17
##
##   Do not install 1.7 with devtools::load_all() or pkgload: they compile at
##   -O0 and every timing would be several times too slow.
##
## Usage
##
##   Rscript sts17_benchmark.R --yelp=/path/to/YelpData.csv --out=bench_out \
##       [--lib16=~/Rlib16] [--lib17=~/Rlib17] [--parts=D,K,V,cores,v16,full] \
##       [--maxcores=8] [--quick]
##
##   --parts  which experiments to run (default: D,K,V,cores,v16,full)
##     D       per-iteration time vs number of documents (Yelp subsamples), K = 10
##     K       per-iteration time vs number of topics, D = 10,000
##     V       per-iteration time vs vocabulary size (prepDocuments lower.thresh)
##     cores   full Yelp, K = 10, threads 1, 2, 4, ... up to --maxcores
##     v16     sts 1.6 vs 1.7 on the same subsamples (1.6 is slow: about 15 s
##             per 1,000 documents per iteration at K = 10)
##     full    full Yelp fit to convergence with 1.7 defaults (the JSS model)
##   Added for the IJOC paper (sts 1.7.1):
##     tts     time to a finished fit on the full corpus: 1.7.1 defaults (1 and
##             --maxcores threads, seeds --seeds) and the 1.6 estimator run in
##             1.7.1 (the recipe below), each to its own stopping rule. With
##             --with16full also sts 1.6 itself to its own stop (hours).
##     matched per-iteration time with the estimator held fixed: sts 1.6
##             against 1.7.1 running the 1.6 recipe (same iterates), on the
##             full corpus and the v16 subsamples
##     ablation the 1.6 recipe in 1.7.1, then each 1.7.1 component in turn:
##             C++ kappa solver, interpolated log-partition, threads; full
##             corpus, per-iteration time
##     Kbig    per-iteration time for K = 10, 25, 50, 100 on the full corpus,
##             1 and --maxcores threads
##   --seeds  comma-separated seeds for tts (default 2024,1,2)
##   --quick  smaller grids and fewer iterations, to check the setup first
##
## The 1.6 recipe (reproduces sts 1.6's iterates in 1.7.1; NEWS, 1.7.1):
##   options(sts.kappaDamping = 0.5, sts.recoverL2 = "expgrad", sts.bound = "1.6")
##   sts(..., kappaSolver = "glmnet", logZ = "exact", allow.neg.change = TRUE,
##       maxIter = 100)
##
## On a Slurm cluster, one job per part keeps each within a time limit, e.g.
##   sbatch --cpus-per-task=8 --mem=16G --time=12:00:00 \
##     --wrap "Rscript sts17_benchmark.R --yelp=... --out=bench --parts=cores --maxcores=8"
## (Request as many CPUs as --maxcores. All parts can share one --out folder.)
##
## Per-iteration time is measured as (time for N iterations - time for 1
## iteration) / (N - 1), so initialization is excluded; it is reported
## separately as the time of the 1-iteration run. Convergence checks are
## disabled in those runs (convTol = -1e10) so every run does exactly N iterations.
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default = NULL) {
  hit <- grep(sprintf("^--%s(=|$)", name), args, value = TRUE)
  if (!length(hit)) return(default)
  v <- sub(sprintf("^--%s=?", name), "", hit[1])
  if (v == "") TRUE else v
}
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))

## ============================ worker mode ==================================
## Rscript sts17_benchmark.R --worker --spec=<file.rds>
## Runs one fit and saves a one-row data frame next to the spec.
if (isTRUE(getarg("worker", FALSE))) {
  spec <- readRDS(getarg("spec"))
  suppressPackageStartupMessages(library(sts, lib.loc = spec$lib))
  corpus <- readRDS(spec$corpus)
  if (!is.null(spec$D) && spec$D < length(corpus$documents)) {
    set.seed(spec$subseed)
    keep <- sort(sample(length(corpus$documents), spec$D))
    corpus$documents <- corpus$documents[keep]
    corpus$meta <- corpus$meta[keep, , drop = FALSE]
    ## drop words that no longer occur, re-index
    used <- sort(unique(unlist(lapply(corpus$documents, function(d) d[1, ]))))
    map <- integer(length(corpus$vocab)); map[used] <- seq_along(used)
    corpus$documents <- lapply(corpus$documents, function(d) { d[1, ] <- map[d[1, ]]; d })
    corpus$vocab <- corpus$vocab[used]
  }
  is16 <- spec$version == "1.6"
  if (length(spec$opts)) do.call(options, spec$opts)
  fit_args <- list(prevalence_sentiment = ~ SubPeriod * eliteStatus + stars + meanTemp + sumPrecip + meanWind + city,
                   initializationVar = ~ stars, corpus = corpus, K = spec$K, verbose = FALSE)
  ## never stop early in timing runs: with convTol = 0, sts 1.6 and
  ## allow.neg.change = TRUE stop at the first decrease of the bound
  if (!spec$converge) fit_args$convTol <- -1e10
  if (!is16) {
    fit_args$cores <- spec$cores
    if (!spec$converge && is.null(spec$extra$allow.neg.change)) fit_args$patience <- 1e6L
  }
  if (length(spec$extra)) fit_args[names(spec$extra)] <- spec$extra
  run <- function(maxIter) {
    set.seed(spec$seed)
    gc(reset = TRUE)
    t0 <- proc.time()[["elapsed"]]
    f <- suppressWarnings(do.call(sts, c(fit_args, list(maxIter = maxIter))))
    list(fit = f, secs = proc.time()[["elapsed"]] - t0, rmem = sum(gc()[, 6]))
  }
  peak_rss <- function() {               # whole-process peak, Linux only
    s <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
    h <- grep("^VmHWM", s, value = TRUE)
    if (length(h)) as.numeric(gsub("[^0-9]", "", h)) / 1024 else NA_real_
  }
  r1 <- run(1L)
  if (spec$converge) {
    rN <- run(spec$maxIter); n <- length(rN$fit$elbo)
    per <- (rN$secs - r1$secs) / max(n - 1, 1)
  } else {
    rN <- run(spec$iters); n <- length(rN$fit$elbo)
    per <- (rN$secs - r1$secs) / (n - 1)
  }
  out <- data.frame(part = spec$part, version = spec$version, D = length(corpus$documents),
                    V = length(corpus$vocab), K = spec$K, cores = if (is16) 1L else spec$cores,
                    iters = n, init_plus_1iter_s = r1$secs, total_s = rN$secs, per_iter_s = per,
                    final_elbo = tail(rN$fit$elbo, 1), best_elbo = max(rN$fit$elbo),
                    best_iter = if (is.null(rN$fit$bestIter)) n else rN$fit$bestIter,
                    r_heap_mb = max(r1$rmem, rN$rmem), peak_rss_mb = peak_rss(),
                    host = Sys.info()[["nodename"]], R = paste(R.version$major, R.version$minor, sep = "."),
                    sts = as.character(packageVersion("sts", lib.loc = spec$lib)),
                    when = format(Sys.time(), "%Y-%m-%d %H:%M"), stringsAsFactors = FALSE)
  if (isTRUE(spec$save_fit)) saveRDS(rN$fit, sub("\\.rds$", "_fit.rds", spec$result))
  saveRDS(out, spec$result)
  quit(save = "no")
}

## ============================ driver mode ==================================
yelp_csv <- getarg("yelp")
out_dir  <- path.expand(getarg("out", "bench_out"))
lib16    <- path.expand(getarg("lib16", "~/Rlib16"))
lib17    <- path.expand(getarg("lib17", "~/Rlib17"))
parts    <- strsplit(getarg("parts", "D,K,V,cores,v16,full"), ",")[[1]]
maxcores <- as.integer(getarg("maxcores", parallel::detectCores()))
quick    <- isTRUE(getarg("quick", FALSE))
seeds    <- as.integer(strsplit(getarg("seeds", "2024,1,2"), ",")[[1]])
with16full <- isTRUE(getarg("with16full", FALSE))
recipe16_opts  <- list(sts.kappaDamping = 0.5, sts.recoverL2 = "expgrad", sts.bound = "1.6")
recipe16_extra <- list(kappaSolver = "glmnet", logZ = "exact", allow.neg.change = TRUE)
if (is.null(yelp_csv)) stop("--yelp=/path/to/YelpData.csv is required")
dir.create(file.path(out_dir, "runs"), recursive = TRUE, showWarnings = FALSE)
for (l in c(lib16, lib17))
  if (!file.exists(file.path(l, "sts", "DESCRIPTION"))) stop("sts is not installed in ", l)
v17 <- read.dcf(file.path(lib17, "sts", "DESCRIPTION"), "Version")[1]
v16 <- read.dcf(file.path(lib16, "sts", "DESCRIPTION"), "Version")[1]
cat(sprintf("sts in %s: %s;  sts in %s: %s\n", lib17, v17, lib16, v16))
if (!startsWith(v17, "1.6.9") && !startsWith(v17, "1.7")) stop("--lib17 does not hold sts 1.7")

## ---- prepare the corpora once (as in the JSS replication script) ----------
prep <- function(thresh) {
  f <- file.path(out_dir, sprintf("yelp_corpus_thresh%d.rds", thresh))
  if (file.exists(f)) return(f)
  cat(sprintf("Preparing Yelp corpus (lower.thresh = %d)...\n", thresh))
  code <- sprintf('
    suppressPackageStartupMessages(library(sts, lib.loc = "%s"))
    pf <- "%s"
    if (file.exists(pf)) processed <- readRDS(pf) else {
      yelp <- read.csv("%s")
      processed <- textProcessor(yelp$text, metadata = yelp, verbose = FALSE)
      saveRDS(processed, pf)
    }
    out <- prepDocuments(processed$documents, processed$vocab, processed$meta, lower.thresh = %d, verbose = FALSE)
    out$meta$SubPeriod <- factor(out$meta$SubPeriod, levels = c("Pre-COVID", "Buildup", "Onset", "Stay_at_Home", "New_Normal"))
    out$meta$eliteStatus <- factor(out$meta$eliteStatus, levels = c("Non-Elite", "Elite"))
    saveRDS(out[c("documents", "vocab", "meta")], "%s")',
    lib17, file.path(out_dir, "yelp_processed.rds"), normalizePath(yelp_csv), thresh, f)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)))
  if (status != 0 || !file.exists(f)) stop("corpus preparation failed")
  f
}

results_file <- file.path(out_dir, "results.csv")
done_keys <- if (file.exists(results_file)) read.csv(results_file, stringsAsFactors = FALSE)$key else character()

run_spec <- function(key, ...) {
  if (key %in% done_keys) { cat("[skip]", key, "\n"); return(invisible()) }
  spec <- list(...)
  spec$result <- file.path(out_dir, "runs", paste0(key, ".rds"))
  spec$lib <- if (spec$version == "1.6") lib16 else lib17
  spec_file <- file.path(out_dir, "runs", paste0(key, "_spec.rds"))
  saveRDS(spec, spec_file)
  cat(sprintf("[run ] %-40s ", key)); flush.console()
  env <- c(sprintf("RCPP_PARALLEL_NUM_THREADS=%d", max(1L, spec$cores)), "OMP_NUM_THREADS=1",
           "OPENBLAS_NUM_THREADS=1", "VECLIB_MAXIMUM_THREADS=1", "MKL_NUM_THREADS=1")
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(self), "--worker", paste0("--spec=", shQuote(spec_file))),
                    env = env)
  if (status != 0 || !file.exists(spec$result)) { cat("FAILED\n"); return(invisible()) }
  r <- readRDS(spec$result); r$key <- key
  cat(sprintf("D=%d V=%d K=%d cores=%d: %.2f s/iter (init+1 iter %.1f s)\n",
              r$D, r$V, r$K, r$cores, r$per_iter_s, r$init_plus_1iter_s))
  write.table(r, results_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(results_file), append = file.exists(results_file))
}

iters <- if (quick) 3L else 6L
c30 <- prep(30)                     # the JSS corpus: D = 35,961, V = 3,041
Dall <- length(readRDS(c30)$documents)

if ("D" %in% parts) {
  Ds <- if (quick) c(2000, 5000) else c(2000, 5000, 10000, 20000, Dall)
  for (D in Ds) for (cc in unique(c(1L, min(2L, maxcores))))
    run_spec(sprintf("D_%d_c%d", D, cc), part = "D", version = "1.7", corpus = c30, D = D, subseed = 1,
             K = 10, cores = cc, iters = iters, converge = FALSE, seed = 2024)
}
if ("K" %in% parts) {
  Ks <- if (quick) c(5, 10) else c(5, 10, 20, 40)
  for (K in Ks)
    run_spec(sprintf("K_%d", K), part = "K", version = "1.7", corpus = c30, D = min(10000, Dall), subseed = 1,
             K = K, cores = 1L, iters = iters, converge = FALSE, seed = 2024)
}
if ("V" %in% parts) {
  th <- if (quick) c(100, 30) else c(200, 100, 30, 10, 5)
  for (t in th)
    run_spec(sprintf("V_thresh%d", t), part = "V", version = "1.7", corpus = prep(t), D = NULL,
             K = 10, cores = 1L, iters = iters, converge = FALSE, seed = 2024)
}
if ("cores" %in% parts) {
  cs <- unique(pmin(c(1L, 2L, 4L, 8L, 16L, 32L), maxcores))
  for (cc in cs)
    run_spec(sprintf("cores_%d", cc), part = "cores", version = "1.7", corpus = c30, D = NULL,
             K = 10, cores = cc, iters = iters, converge = FALSE, seed = 2024)
}
if ("v16" %in% parts) {
  Ds <- if (quick) 1000 else c(1000, 2000, 5000, 10000)
  for (D in Ds) for (ver in c("1.6", "1.7"))
    run_spec(sprintf("v16_%s_D%d", ver, D), part = "v16", version = ver, corpus = c30, D = D, subseed = 1,
             K = 10, cores = 1L, iters = 3L, converge = FALSE, seed = 2024)
}
if ("full" %in% parts) {
  run_spec("full_c1", part = "full", version = "1.7", corpus = c30, D = NULL, K = 10, cores = 1L,
           maxIter = if (quick) 20L else 500L, converge = TRUE, seed = 2024, save_fit = TRUE)
  if (maxcores > 1)
    run_spec(sprintf("full_c%d", maxcores), part = "full", version = "1.7", corpus = c30, D = NULL, K = 10,
             cores = maxcores, maxIter = if (quick) 20L else 500L, converge = TRUE, seed = 2024)
}

if ("tts" %in% parts) {
  for (sd in seeds) {
    for (cc in unique(c(1L, maxcores)))
      run_spec(sprintf("tts_171_c%d_s%d", cc, sd), part = "tts", version = "1.7", corpus = c30, D = NULL, K = 10,
               cores = cc, maxIter = if (quick) 20L else 500L, converge = TRUE, seed = sd)
    run_spec(sprintf("tts_recipe16_c1_s%d", sd), part = "tts", version = "1.7", corpus = c30, D = NULL, K = 10,
             cores = 1L, maxIter = if (quick) 10L else 100L, converge = TRUE, seed = sd,
             opts = recipe16_opts, extra = recipe16_extra)
    if (with16full)
      run_spec(sprintf("tts_16_s%d", sd), part = "tts", version = "1.6", corpus = c30, D = NULL, K = 10,
               cores = 1L, maxIter = if (quick) 5L else 100L, converge = TRUE, seed = sd)
  }
}
if ("matched" %in% parts) {
  Ds <- if (quick) 1000 else c(1000, 2000, 5000, 10000, Dall)
  for (D in Ds) {
    if (D == Dall || !("v16" %in% parts))
      run_spec(sprintf("matched_16_D%d", D), part = "matched", version = "1.6", corpus = c30, D = D, subseed = 1,
               K = 10, cores = 1L, iters = 3L, converge = FALSE, seed = 2024)
    run_spec(sprintf("matched_recipe16_D%d", D), part = "matched", version = "1.7", corpus = c30, D = D, subseed = 1,
             K = 10, cores = 1L, iters = if (quick) 3L else 4L, converge = FALSE, seed = 2024,
             opts = recipe16_opts, extra = recipe16_extra)
  }
}
if ("ablation" %in% parts) {
  steps <- list(
    recipe16   = list(opts = recipe16_opts, extra = recipe16_extra, cores = 1L),
    cppkappa   = list(opts = recipe16_opts, extra = list(kappaSolver = "cpp", logZ = "exact", allow.neg.change = TRUE), cores = 1L),
    tables     = list(opts = recipe16_opts, extra = list(kappaSolver = "cpp", logZ = "interpolated", allow.neg.change = TRUE), cores = 1L),
    threads    = list(opts = recipe16_opts, extra = list(kappaSolver = "cpp", logZ = "interpolated", allow.neg.change = TRUE), cores = maxcores),
    defaults171 = list(opts = list(), extra = list(), cores = maxcores))
  for (nm in names(steps)) { st <- steps[[nm]]
    run_spec(sprintf("ablation_%s", nm), part = "ablation", version = "1.7", corpus = c30, D = NULL, K = 10,
             cores = st$cores, iters = if (quick) 3L else 5L, converge = FALSE, seed = 2024, opts = st$opts, extra = st$extra) }
}
if ("Kbig" %in% parts) {
  Ks <- if (quick) c(10, 25) else c(10, 25, 50, 100)
  for (K in Ks) for (cc in unique(c(1L, maxcores)))
    run_spec(sprintf("Kbig_%d_c%d", K, cc), part = "Kbig", version = "1.7", corpus = c30, D = NULL, K = K,
             cores = cc, iters = 3L, converge = FALSE, seed = 2024)
}

if (file.exists(results_file)) {
  res <- read.csv(results_file, stringsAsFactors = FALSE)
  cat("\n==== summary (seconds per EM iteration) ====\n")
  print(res[, c("key", "version", "D", "V", "K", "cores", "iters", "per_iter_s", "init_plus_1iter_s", "peak_rss_mb")],
        row.names = FALSE, digits = 3)
  v <- res[res$part == "v16", ]
  if (nrow(v)) {
    w <- reshape(v[, c("D", "version", "per_iter_s")], idvar = "D", timevar = "version", direction = "wide")
    if (all(c("per_iter_s.1.6", "per_iter_s.1.7") %in% names(w))) {
      w$speedup <- w$per_iter_s.1.6 / w$per_iter_s.1.7
      cat("\n==== 1.6 vs 1.7, one thread ====\n"); print(w, row.names = FALSE, digits = 3)
    }
  }
  cat("\nAll results:", normalizePath(results_file), "\n")
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
