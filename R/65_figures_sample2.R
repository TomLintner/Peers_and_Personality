## 65_figures_sample2.R
## source: Results/scripts_v2/oz_plots.R
source("00_setup.R")

## Moment of Change trace plots over the full 6000-draw chain, with the 1200-draw burn-in
## marked so the transient is visible rather than hidden.
suppressMessages(library(RSienaTest))
s <- readRDS(file.path(OBJ, "oz_bayes_setup.rds")); E <- s$eff; class(E) <- "data.frame"
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
fxi  <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labMu <- E$effectName[selR]; labEta <- E$effectName[fxi]
MU <- NULL; ET <- NULL; bnd <- c()
for (k in 1:20) {
  f <- sprintf(file.path(CHN, "ozbayes_stage%d.rds"), k); if (!file.exists(f)) next
  g <- readRDS(f); ok <- !is.na(g$ThinPosteriorMu[,1]); if (k==1) ok[1:100] <- FALSE
  MU <- rbind(MU, g$ThinPosteriorMu[ok,,drop=FALSE])
  ET <- rbind(ET, g$ThinPosteriorEta[ok,,drop=FALSE])
  bnd <- c(bnd, nrow(MU))
}
BURN <- 1200
cat("draws:", nrow(MU), "| burn-in marked at", BURN, "\n")
trace <- function(M, labs, file, ttl) {
  png(file, width = 1900, height = 1500, res = 110)
  op <- par(mfrow = c(ceiling(ncol(M)/5), 5), mar = c(2.2, 2.6, 1.8, 0.6),
            mgp = c(1.5, 0.5, 0), oma = c(0, 0, 2.6, 0), cex.main = 0.95)
  for (j in seq_len(ncol(M))) {
    plot(M[,j], type = "l", col = "grey40", lwd = 0.4, xlab = "", ylab = "", main = labs[j])
    rect(0, par("usr")[3], BURN, par("usr")[4], col = adjustcolor("grey70", .30), border = NA)
    lines(M[,j], col = "grey40", lwd = 0.4)
    abline(v = BURN, col = "purple", lwd = 1.6)
    abline(h = mean(M[(BURN+1):nrow(M), j]), col = "firebrick", lwd = 1.4)
    lines(seq_len(nrow(M)), filter(M[,j], rep(1/101, 101), sides = 2),
          col = "orange", lwd = 1.3)
  }
  par(op); mtext(ttl, outer = TRUE, line = 0.8, cex = 1.1, font = 2)
  dev.off() }
trace(MU, labMu, "figures/moment_of_change_trace_mu.png",
  sprintf("Moment of Change sienaBayes: population means mu, %d draws (grey = 1200-draw burn-in discarded; red = retained mean; orange = 101-draw MA)", nrow(MU)))
trace(ET, labEta, "figures/moment_of_change_trace_eta.png",
  sprintf("Moment of Change sienaBayes: fixed effects eta, %d draws (grey = 1200-draw burn-in discarded)", nrow(ET)))
cat("written the Moment of Change trace figures\n")
