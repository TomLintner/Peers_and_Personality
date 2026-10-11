## 20_bayes_sample1.R
## source: RSiena_build/bayes/bayes_full.R
source("00_setup.R")

## =====================================================================
## GACR multilevel SAOM -- FULL Bayesian estimation
## 27 classrooms, 1 network + 5 Big Five behaviours.
## Random: density, recip, gwespFF + linear/quad for each trait; rates
## always random => p = 19.  Fixed (eta): the other 30 evaluation effects.
## 8 cores, via the patched parRapply re-nesting in RSienaTest 1.2-30.
## seed 2012 throughout.
## nprewarm 50 + nwarm 200 + nmain 1000 = 1250 iterations
## measured 22 s/iteration on 8 cores => ~7.6 h
## saveFreq 20 => PartialBayesResult.RData refreshed roughly every 7 min,
##   and that file is usable as prevBayes to continue if anything stops.
## =====================================================================
suppressMessages(library(RSienaTest))
set.seed(2012)
SC <- WORK
WD <- file.path(SC, "bayes_full"); dir.create(WD, showWarnings = FALSE)
setwd(WD)

s    <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds"))
ans0 <- readRDS(file.path(OBJ, "gacr_bayes_init.rds"))
cat("groups:", length(s$grp), "| random dimension p =", length(s$priorMu), "\n")
cat("start:", format(Sys.time()), "\n"); flush.console()

algo <- sienaAlgorithmCreate(projname = "gacr_bayes_full", seed = 2012)
T0 <- proc.time()
gb <- try(sienaBayes(algo, data = s$grp, effects = s$eff,
                     priorMu = s$priorMu, priorSigma = s$priorSig,
                     priorKappa = 0.1,
                     nwarm = 200, nmain = 1000, nrunMHBatches = 20,
                     saveFreq = 20,
                     prevAns = ans0,
                     nbrNodes = 8, clusterType = "PSOCK",
                     silentstart = TRUE),
          silent = TRUE)
el <- (proc.time() - T0)[3]
cat("\n=== elapsed:", round(el / 3600, 2), "hours ===\n")
cat("finish:", format(Sys.time()), "\n")

if (inherits(gb, "try-error")) {
  cat("*** RUN FAILED ***\n"); print(gb)
  cat("PartialBayesResult.RData in", WD, "can be used as prevBayes.\n")
} else {
  saveRDS(gb, file.path(SC, "gacr_bayes_full.rds"))
  cat("saved gacr_bayes_full.rds\n")

  cat("\n================ POSTERIOR SUMMARY ================\n")
  r <- try(shortBayesResults(gb), silent = TRUE)
  if (!inherits(r, "try-error")) print(r) else cat("shortBayesResults failed\n")

  ## --- convergence: split the main phase and compare halves (Geweke-style) ---
  cat("\n================ CONVERGENCE CHECK ================\n")
  mu <- try(gb$ThinPosteriorMu, silent = TRUE)
  if (!inherits(mu, "try-error") && !is.null(mu)) {
    ok <- !is.na(mu[, 1]); mu <- mu[ok, , drop = FALSE]
    n <- nrow(mu); cat("retained draws:", n, "\n")
    if (n >= 20) {
      h1 <- seq_len(floor(n / 3)); h2 <- (n - floor(n / 3) + 1):n
      m1 <- colMeans(mu[h1, , drop = FALSE]); m2 <- colMeans(mu[h2, , drop = FALSE])
      sd1 <- apply(mu[h1, , drop = FALSE], 2, sd); sd2 <- apply(mu[h2, , drop = FALSE], 2, sd)
      zst <- (m1 - m2) / sqrt(sd1^2 / length(h1) + sd2^2 / length(h2))
      nm <- s$eff; class(nm) <- "data.frame"
      sel <- nm$include & ((!nm$basicRate & nm$randomEffects) | (nm$basicRate & nm$group == 1))
      lab <- paste(nm$name[sel], nm$shortName[sel])
      out <- data.frame(par = lab, first_third = round(m1, 3), last_third = round(m2, 3),
                        z = round(zst, 2), flag = ifelse(abs(zst) > 2, "<<<", ""))
      print(out, row.names = FALSE)
      cat("\n|z| > 2 suggests the chain had not settled for that parameter.\n")
      cat("parameters flagged:", sum(abs(zst) > 2, na.rm = TRUE), "of", length(zst), "\n")
    }
  }

  ## --- MDS of group-wise posterior means: too many random effects? ---
  cat("\n=== MDS diagnostic (dense core + shifting outliers = too many random effects) ===\n")
  p1 <- try({ pdf(file.path(SC, "gacr_bayes_MDS.pdf"), width = 8, height = 8)
              plotPostMeansMDS(gb); dev.off(); "ok" }, silent = TRUE)
  cat("plotPostMeansMDS:", if (inherits(p1, "try-error")) "failed" else "written to gacr_bayes_MDS.pdf", "\n")

  ## --- trace plots ---
  p2 <- try({ pdf(file.path(SC, "gacr_bayes_traces.pdf"), width = 11, height = 8)
              op <- par(mfrow = c(4, 5), mar = c(3, 3, 2, 1))
              for (j in seq_len(ncol(gb$ThinPosteriorMu)))
                plot(gb$ThinPosteriorMu[, j], type = "l", main = paste("mu", j),
                     xlab = "", ylab = "")
              par(op); dev.off(); "ok" }, silent = TRUE)
  cat("trace plots:", if (inherits(p2, "try-error")) "failed" else "written to gacr_bayes_traces.pdf", "\n")
}
cat("\nDONE", format(Sys.time()), "\n")
