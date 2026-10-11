## 21_bayes_sample1_stages.R
## source: Results/scripts/bayes_12h.R
source("00_setup.R")

## ==========================================================================
## GACR sienaBayes -- staged full-scale run, >= 12 hours, proper weak priors.
## Full model: 1 network + 5 behaviours, 43 evaluation effects, 205 parameters,
## p1 = 19 random (density, recip, gwespFF, linear+quad x5, + 6 rates),
## p2 = 30 fixed (all 15 personality selection terms, all cross-trait controls,
## gender, inPopSqrt/outActSqrt, and the 10 influence terms).
## Priors: priorSigEta = 1 (proper, weak), priorKappa = 1.0.
## Sampling: nSampVarying = 3, nSampConst = 10, nrunMHBatches = 20.
## Runs in chunks; after EACH chunk it glues to the accumulated chain, saves,
## and reports median/max split-Rhat and median ESS for mu and eta.
## ==========================================================================
suppressMessages(library(RSienaTest))
set.seed(2012)
SC <- WORK
WD <- file.path(SC, "bayes12h"); dir.create(WD, showWarnings = FALSE); setwd(WD)
s    <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds"))
ans0 <- readRDS(file.path(OBJ, "gacr_bayes_init.rds"))
E <- s$eff; class(E) <- "data.frame"
p2 <- sum(E$include & !E$basicRate & !E$randomEffects & E$type == "eval")
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
labMu <- paste(E$name[selR], E$shortName[selR])
labEta <- paste(E$name[E$include & !E$basicRate & !E$randomEffects & E$type == "eval"],
                E$shortName[E$include & !E$basicRate & !E$randomEffects & E$type == "eval"])

BUDGET <- 12.5 * 3600      # seconds
NWARM  <- 100
CHUNK  <- 250

splitRhat <- function(x) { n <- length(x); h <- floor(n / 2)
  if (h < 10) return(NA)
  a <- x[1:h]; b <- x[(n - h + 1):n]
  W <- mean(c(var(a), var(b))); if (W <= 0) return(NA)
  gm <- mean(c(mean(a), mean(b)))
  B <- h * sum((c(mean(a), mean(b)) - gm)^2)
  sqrt((((h - 1) / h) * W + B / h) / W) }
rk <- function(x) qnorm((rank(x) - 3/8) / (length(x) + 1/4))
ess <- function(x) { n <- length(x)
  r <- acf(x, lag.max = min(300, n - 2), plot = FALSE)$acf[-1]
  k <- which(r < 0.05)[1]; if (is.na(k)) k <- length(r)
  max(1, n / (1 + 2 * sum(r[1:max(1, k - 1)]))) }

diag1 <- function(M, tag, labs) {
  R <- apply(M, 2, function(c) splitRhat(rk(c)))
  S <- apply(M, 2, ess)
  cat(sprintf("  %-4s n=%4d | Rhat med %.3f max %.3f | <1.01: %2d/%2d | ESS med %4.0f min %4.0f\n",
      tag, nrow(M), median(R, na.rm = TRUE), max(R, na.rm = TRUE),
      sum(R < 1.01, na.rm = TRUE), length(R), median(S), min(S)))
  w <- order(-R)[1:min(3, length(R))]
  cat("       worst Rhat:", paste(sprintf("%s %.3f", labs[w], R[w]), collapse = " | "), "\n")
  invisible(data.frame(par = labs, Rhat = round(R, 4), ESS = round(S)))
}
report <- function(g, stage, nw) {
  mu <- g$ThinPosteriorMu; mu <- mu[!is.na(mu[, 1]), , drop = FALSE]
  et <- g$ThinPosteriorEta; et <- et[!is.na(et[, 1]), , drop = FALSE]
  if (nw > 0 && nrow(mu) > nw) { mu <- mu[-(1:nw), , drop = FALSE]; et <- et[-(1:nw), , drop = FALSE] }
  cat("\n############ CUMULATIVE DIAGNOSTICS after stage", stage, "############\n")
  a <- diag1(mu, "mu", labMu); b <- diag1(et, "eta", labEta)
  write.csv(rbind(cbind(block = "mu", a), cbind(block = "eta", b)),
            file.path(SC, sprintf("bayes12h_diag_stage%d.csv", stage)), row.names = FALSE)
  flush.console()
}

runChunk <- function(prev, nwarm, nmain, stage) {
  algo <- sienaAlgorithmCreate(projname = paste0("b12h_s", stage), seed = 2012)
  args <- list(algo, data = s$grp, effects = s$eff,
               priorMu = s$priorMu, priorSigma = s$priorSig, priorKappa = 1.0,
               priorSigEta = rep(1, p2),
               nmain = nmain, nrunMHBatches = 20,
               nSampVarying = 3, nSampConst = 10, saveFreq = 20,
               nbrNodes = 8, clusterType = "PSOCK", silentstart = TRUE)
  if (is.null(prev)) { args$nprewarm <- 20; args$nwarm <- nwarm; args$prevAns <- ans0 }
  else               { args$prevBayes <- prev }
  do.call(sienaBayes, args)
}

T0 <- proc.time()[3]; stage <- 0; acc <- NULL; prev <- NULL
cat("start:", format(Sys.time()), " budget", BUDGET / 3600, "hours\n"); flush.console()
repeat {
  stage <- stage + 1
  cat("\n===== stage", stage, "=====", format(Sys.time()), "\n"); flush.console()
  g <- try(runChunk(prev, NWARM, CHUNK, stage), silent = TRUE)
  if (inherits(g, "try-error")) { cat("stage", stage, "FAILED:\n"); print(g); break }
  saveRDS(g, file.path(SC, sprintf("bayes12h_stage%d.rds", stage)))
  acc <- if (is.null(acc)) g else {
    gl <- try(glueBayes(acc, g, nwarm2 = 0), silent = TRUE)
    if (inherits(gl, "try-error")) { cat("  glueBayes failed; reporting this chunk only\n"); g } else gl }
  saveRDS(acc, file.path(SC, "bayes12h_accumulated.rds"))
  report(acc, stage, if (stage == 1) NWARM else 0)
  el <- proc.time()[3] - T0
  cat("  elapsed", round(el / 3600, 2), "h of", BUDGET / 3600, "\n"); flush.console()
  prev <- g
  if (el > BUDGET) { cat("\nbudget reached.\n"); break }
}
cat("\nTOTAL:", round((proc.time()[3] - T0) / 3600, 2), "hours |", stage, "stages\n")
cat("DONE", format(Sys.time()), "\n")
