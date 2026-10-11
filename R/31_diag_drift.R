## 31_diag_drift.R
## source: Results/scripts/oz_drift.R
source("00_setup.R")

## Is the OZ chain drifting, or is the rising max Rhat just noise?
## For each parameter: OLS slope over the whole chain, expressed in posterior
## SDs per 1000 draws, and the standardised difference between the first and
## second half. Systematic drift makes split-Rhat WORSE as the chain lengthens;
## pure noise does not.
suppressMessages(library(RSienaTest))
s <- readRDS(file.path(OBJ, "oz_bayes_setup.rds")); E <- s$eff; class(E) <- "data.frame"
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
fxi  <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labMu <- E$effectName[selR]; labEta <- E$effectName[fxi]
MU <- NULL; ET <- NULL
for (k in 1:20) {
  f <- file.path(CHN, sprintf("ozbayes_stage%d.rds", k)); if (!file.exists(f)) next
  g <- readRDS(f); ok <- !is.na(g$ThinPosteriorMu[,1]); if (k==1) ok[1:100] <- FALSE
  MU <- rbind(MU, g$ThinPosteriorMu[ok,,drop=FALSE]); ET <- rbind(ET, g$ThinPosteriorEta[ok,,drop=FALSE])
}
cat("draws:", nrow(MU), "\n")
rep1 <- function(M, labs, tag) {
  n <- nrow(M); i <- seq_len(n); h <- floor(n/2)
  out <- t(sapply(seq_len(ncol(M)), function(j) {
    x <- M[,j]; sd0 <- sd(x)
    sl <- coef(lm(x ~ i))[2] * 1000 / sd0                 # SD per 1000 draws
    d  <- (mean(x[(n-h+1):n]) - mean(x[1:h])) / sd0        # half-to-half, in SD
    c(slope = sl, halfdiff = d) }))
  rownames(out) <- labs
  o <- order(-abs(out[,"halfdiff"]))
  cat("\n=== ", tag, ": largest first-half vs second-half shifts ===\n", sep="")
  print(round(out[o[1:10], , drop=FALSE], 3))
  cat("  |halfdiff| > 0.25 SD:", sum(abs(out[,"halfdiff"]) > 0.25), "of", nrow(out), "\n")
  cat("  |halfdiff| > 0.10 SD:", sum(abs(out[,"halfdiff"]) > 0.10), "of", nrow(out), "\n")
  out }
a <- rep1(MU, labMu, "mu"); b <- rep1(ET, labEta, "eta")
write.csv(rbind(data.frame(block="mu", par=labMu, a), data.frame(block="eta", par=labEta, b)),
          "oz_drift.csv", row.names = FALSE)
