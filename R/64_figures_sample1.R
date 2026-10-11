## 64_figures_sample1.R
## source: Results/scripts_v2/gacr_plots.R
source("00_setup.R")

## Defending Victims of Bullying Bayesian diagnostics on the full 2650-draw chain:
##  (1) trace plots for mu and eta
##  (2) Snijders' MDS of group-wise posterior means (too many random effects?)
suppressMessages(library(RSienaTest))
SC <- getwd()
s <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds")); E <- s$eff; class(E) <- "data.frame"
G <- readRDS(file.path(OBJ, "gacr_siena.rds"))
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
labMu  <- E$effectName[selR]
fxi <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labEta <- paste(E$effectName[fxi], ifelse(E$name[fxi]=="net","","(beh)"))

MU <- NULL; ET <- NULL; TP <- NULL; bnd <- c()
for (k in 1:9) {
  f <- sprintf(file.path(CHN, "bayes12h_stage%d.rds"), k); if (!file.exists(f)) next
  g <- readRDS(f); ok <- !is.na(g$ThinPosteriorMu[, 1]); if (k == 1) ok[1:100] <- FALSE
  MU <- rbind(MU, g$ThinPosteriorMu[ok, , drop = FALSE])
  ET <- rbind(ET, g$ThinPosteriorEta[ok, , drop = FALSE])
  tp <- g$ThinParameters
  if (!is.null(tp)) TP <- if (is.null(TP)) tp[ok, , , drop = FALSE] else
      array(c(aperm(TP, c(2,3,1)), aperm(tp[ok,,,drop=FALSE], c(2,3,1))),
            dim = c(dim(TP)[2], dim(TP)[3], dim(TP)[1] + sum(ok)))
  bnd <- c(bnd, nrow(MU))
}
if (length(dim(TP)) == 3 && dim(TP)[1] != nrow(MU)) TP <- aperm(TP, c(3,1,2))
cat("draws:", nrow(MU), "| ThinParameters dim:", paste(dim(TP), collapse=" x "), "\n")

trace <- function(M, labs, file, ttl) {
  png(file, width = 1900, height = 1500, res = 110)
  op <- par(mfrow = c(ceiling(ncol(M)/5), 5), mar = c(2.2, 2.6, 1.8, 0.6),
            mgp = c(1.5, 0.5, 0), oma = c(0, 0, 2.4, 0), cex.main = 0.95)
  for (j in seq_len(ncol(M))) {
    plot(M[, j], type = "l", col = "grey35", lwd = 0.5, xlab = "", ylab = "",
         main = labs[j])
    abline(v = head(bnd, -1), col = "steelblue", lty = 3)
    abline(h = mean(M[, j]), col = "firebrick", lwd = 1.4)
    lines(seq_len(nrow(M)), filter(M[, j], rep(1/51, 51), sides = 2),
          col = "orange", lwd = 1.3)
  }
  par(op); mtext(ttl, outer = TRUE, line = 0.6, cex = 1.15, font = 2)
  dev.off() }
trace(MU, labMu, "figures/defending_victims_trace_mu.png",
      sprintf("Defending Victims of Bullying sienaBayes: trace of population means mu (%d draws; red = mean, orange = 51-draw MA, blue = stage joins)", nrow(MU)))
trace(ET, labEta, "figures/defending_victims_trace_eta.png",
      sprintf("Defending Victims of Bullying sienaBayes: trace of fixed effects eta (%d draws)", nrow(ET)))

q("no")
## ---- skipped, see gacr_mds2.R ----
pm <- apply(TP, c(2, 3), mean)                  # groups x random params
z <- scale(pm); z[is.na(z)] <- 0
d <- dist(z); mds <- cmdscale(d, k = 2)
cen <- colMeans(mds); rad <- sqrt(rowSums(sweep(mds, 2, cen)^2))
png("figures/defending_victims_mds.png", width = 1150, height = 1000, res = 130)
par(mar = c(4, 4, 3.4, 1))
plot(mds, type = "n", xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main = "Group-wise posterior means, MDS\n(dense core + detached outliers = too many random effects)")
sz <- G$desc$n
symbols(mds[,1], mds[,2], circles = rep(diff(range(mds))/45, nrow(mds)),
        inches = FALSE, add = TRUE, fg = "grey60",
        bg = adjustcolor("steelblue", 0.25))
text(mds, labels = G$classes, cex = 0.6)
points(cen[1], cen[2], pch = 3, col = "firebrick", cex = 1.6, lwd = 2)
for (m in c(1, 2, 3) * median(rad))
  symbols(cen[1], cen[2], circles = m, inches = FALSE, add = TRUE,
          fg = adjustcolor("firebrick", 0.5), lty = 2)
legend("topleft", bty = "n", cex = 0.7,
       legend = c("classroom", "centroid", "1x / 2x / 3x median radius"),
       pch = c(21, 3, NA), lty = c(NA, NA, 2),
       col = c("grey60", "firebrick", "firebrick"),
       pt.bg = c(adjustcolor("steelblue", .25), NA, NA))
dev.off()
o <- order(-rad)
cat("\n=== distance from MDS centroid ===\n")
print(data.frame(class = G$classes[o], n = G$desc$n[o], dist = round(rad[o], 2)),
      row.names = FALSE)
cat("\nmedian", round(median(rad),2), "| max", round(max(rad),2),
    "| max/median", round(max(rad)/median(rad),2),
    "| beyond 3x median:", sum(rad > 3*median(rad)), "\n")
cat("\nwritten the Defending Victims of Bullying trace and MDS figures\n")
