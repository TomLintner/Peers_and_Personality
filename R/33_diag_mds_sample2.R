## 33_diag_mds_sample2.R
## source: Results/scripts_v2/oz_mds.R
source("00_setup.R")

## Snijders MDS diagnostic for Moment of Change: 18 classrooms, 19 random effects.
suppressMessages(library(RSienaTest))
s <- readRDS(file.path(OBJ, "oz_bayes_setup.rds"))
cls <- names(s$grp)
nn  <- sapply(s$grp, function(d) attr(d$depvars[[1]], "netdims")[1])
SUM <- NULL; N <- 0; dn <- NULL
for (k in 1:20) {
  f <- sprintf(file.path(CHN, "ozbayes_stage%d.rds"), k); if (!file.exists(f)) next
  g <- readRDS(f); tp <- g$ThinParameters
  ok <- !is.na(tp[,1,1]); if (k == 1) ok[1:100] <- FALSE
  tp <- tp[ok,,,drop=FALSE]
  if (is.null(SUM)) { SUM <- apply(tp, c(2,3), sum); dn <- dimnames(g$ThinParameters)[[3]] }
  else SUM <- SUM + apply(tp, c(2,3), sum)
  N <- N + sum(ok)
}
pm <- SUM / N; colnames(pm) <- dn
cat("draws:", N, "| group x param:", paste(dim(pm), collapse=" x "), "\n")
keep <- grepl("rate", dn) | dn == "outdegree (density)" | dn == "reciprocity" |
        grepl("GWESP", dn) | grepl("linear shape|quadratic shape", dn)
cat("random params matched:", sum(keep), "\n")
pmr <- pm[, keep, drop = FALSE]
z <- scale(pmr); z[is.na(z)] <- 0
mds <- cmdscale(dist(z), k = 2)
cen <- colMeans(mds); rad <- sqrt(rowSums(sweep(mds, 2, cen)^2))
png("figures/moment_of_change_mds.png", width = 1200, height = 1050, res = 130)
par(mar = c(4,4,3.8,1)); lim <- range(mds)*1.12
plot(mds, type="n", xlim=lim, ylim=lim, asp=1, xlab="MDS dimension 1", ylab="MDS dimension 2",
     main=sprintf("Moment of Change: group-wise posterior means, MDS (%d draws)\ndense core + detached outliers = too many random effects", N))
for (m in c(1,2,3)*median(rad))
  symbols(cen[1], cen[2], circles=m, inches=FALSE, add=TRUE, fg=adjustcolor("firebrick",.45), lty=2)
points(mds, pch=21, cex=1.9+2.2*(nn-min(nn))/max(1,diff(range(nn))),
       bg=adjustcolor("darkorange",.30), col="grey45")
text(mds, labels=cls, cex=0.5)
points(cen[1], cen[2], pch=3, col="firebrick", cex=1.7, lwd=2)
legend("topleft", bty="n", cex=0.68,
       legend=c("classroom (size ~ n actors)","centroid","1x / 2x / 3x median radius"),
       pch=c(21,3,NA), lty=c(NA,NA,2), col=c("grey45","firebrick","firebrick"),
       pt.bg=c(adjustcolor("darkorange",.30),NA,NA))
dev.off()
o <- order(-rad)
cat("\n=== distance from MDS centroid (18 classrooms) ===\n")
print(data.frame(class=cls[o], n=nn[o], dist=round(rad[o],2)), row.names=FALSE)
cat("\nmedian", round(median(rad),2), "| max", round(max(rad),2),
    "| max/median", round(max(rad)/median(rad),2),
    "| beyond 3x median:", sum(rad > 3*median(rad)), "\n")
