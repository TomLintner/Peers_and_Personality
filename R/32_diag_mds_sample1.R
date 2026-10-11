## 32_diag_mds_sample1.R
## source: Results/scripts_v2/gacr_mds2.R
source("00_setup.R")

## Snijders' MDS diagnostic on the full 2650-draw chain.
## Accumulate the group x parameter posterior MEAN across stages by summing
## (avoids building a 2650 x 27 x 49 array) and restrict to the 19 RANDOM
## parameters -- the fixed eta components are identical across groups and
## would contribute nothing but noise.
suppressMessages(library(RSienaTest))
s <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds")); E <- s$eff; class(E) <- "data.frame"
G <- readRDS(file.path(OBJ, "gacr_siena.rds"))
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
randNames <- E$effectName[selR]

SUM <- NULL; N <- 0
for (k in 1:9) {
  f <- sprintf(file.path(CHN, "bayes12h_stage%d.rds"), k); if (!file.exists(f)) next
  g <- readRDS(f); tp <- g$ThinParameters
  ok <- !is.na(tp[, 1, 1]); if (k == 1) ok[1:100] <- FALSE
  tp <- tp[ok, , , drop = FALSE]
  if (is.null(SUM)) { SUM <- apply(tp, c(2,3), sum); dn <- dimnames(g$ThinParameters)[[3]] }
  else SUM <- SUM + apply(tp, c(2,3), sum)
  N <- N + sum(ok)
}
pm <- SUM / N
colnames(pm) <- dn
cat("draws used:", N, "| group x param matrix:", paste(dim(pm), collapse=" x "), "\n")
keep <- grepl("rate", dn) | dn == "outdegree (density)" | dn == "reciprocity" | grepl("GWESP", dn) | grepl("linear shape|quadratic shape", dn)
cat("random parameters matched:", sum(keep), "of", length(randNames), "\n")
cat("  ", paste(dn[keep], collapse=" | "), "\n")
pmr <- pm[, keep, drop = FALSE]

z <- scale(pmr); z[is.na(z)] <- 0
mds <- cmdscale(dist(z), k = 2)
cen <- colMeans(mds); rad <- sqrt(rowSums(sweep(mds, 2, cen)^2))
cls <- G$classes; nn <- G$desc$n

png("figures/defending_victims_mds.png", width = 1200, height = 1050, res = 130)
par(mar = c(4, 4, 3.8, 1))
lim <- range(mds) * 1.12
plot(mds, type = "n", xlim = lim, ylim = lim, asp = 1,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main = "Defending Victims of Bullying: group-wise posterior means, MDS\ndense core + detached outliers = too many random effects")
for (m in c(1,2,3) * median(rad))
  symbols(cen[1], cen[2], circles = m, inches = FALSE, add = TRUE,
          fg = adjustcolor("firebrick", 0.45), lty = 2)
points(mds, pch = 21, cex = 1.9 + 2.2*(nn-min(nn))/diff(range(nn)),
       bg = adjustcolor("steelblue", 0.30), col = "grey45")
text(mds, labels = cls, cex = 0.5)
points(cen[1], cen[2], pch = 3, col = "firebrick", cex = 1.7, lwd = 2)
legend("topleft", bty = "n", cex = 0.68,
       legend = c("classroom (size ~ n actors)", "centroid", "1x / 2x / 3x median radius"),
       pch = c(21, 3, NA), lty = c(NA, NA, 2),
       col = c("grey45", "firebrick", "firebrick"),
       pt.bg = c(adjustcolor("steelblue", .30), NA, NA))
dev.off()

o <- order(-rad)
cat("\n=== distance from MDS centroid (27 classrooms) ===\n")
print(data.frame(class = cls[o], n = nn[o], dist = round(rad[o], 2)), row.names = FALSE)
cat("\nmedian", round(median(rad),2), "| max", round(max(rad),2),
    "| max/median", round(max(rad)/median(rad),2),
    "| beyond 3x median:", sum(rad > 3*median(rad)), "\n")
write.csv(data.frame(class = cls, n = nn, mds1 = mds[,1], mds2 = mds[,2],
                     dist = rad), "tables/mds_coords_defending_victims.csv", row.names = FALSE)
