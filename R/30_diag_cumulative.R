## 30_diag_cumulative.R
## source: Results/scripts/cum_diag.R
source("00_setup.R")

## Cumulative diagnostics built by stacking the per-stage draw matrices directly.
## glueBayes cannot re-glue an already-glued object ("argument is of length zero"),
## so the in-script accumulation froze at 650 draws. Stacking the matrices is
## equivalent for Rhat/ESS purposes: prevBayes continues the SAME chain state,
## so the concatenation is one continuous chain.
suppressMessages(library(RSienaTest))
SC <- WORK
s <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds"))
E <- s$eff; class(E) <- "data.frame"
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
labMu <- paste(E$name[selR], E$shortName[selR])
fxi <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labEta <- paste(E$name[fxi], E$shortName[fxi])

NWARM1 <- 100
stack <- function() {
  MU <- NULL; ET <- NULL; used <- c()
  for (k in 1:9) {
    f <- file.path(CHN, sprintf("bayes12h_stage%d.rds", k))
    if (!file.exists(f)) next
    g <- readRDS(f)
    mu <- g$ThinPosteriorMu; mu <- mu[!is.na(mu[, 1]), , drop = FALSE]
    et <- g$ThinPosteriorEta; et <- et[!is.na(et[, 1]), , drop = FALSE]
    if (k == 1) { mu <- mu[-(1:NWARM1), , drop = FALSE]; et <- et[-(1:NWARM1), , drop = FALSE] }
    MU <- rbind(MU, mu); ET <- rbind(ET, et); used <- c(used, sprintf("s%d:%d", k, nrow(mu)))
  }
  cat("stages stacked:", paste(used, collapse = "  "), "\n")
  list(mu = MU, et = ET)
}
splitRhat <- function(x) { n <- length(x); h <- floor(n / 2)
  if (h < 10) return(NA)
  a <- x[1:h]; b <- x[(n - h + 1):n]
  W <- mean(c(var(a), var(b))); if (W <= 0) return(NA)
  gm <- mean(c(mean(a), mean(b)))
  B <- h * sum((c(mean(a), mean(b)) - gm)^2)
  sqrt((((h - 1) / h) * W + B / h) / W) }
rk  <- function(x) qnorm((rank(x) - 3/8) / (length(x) + 1/4))
ess <- function(x) { n <- length(x)
  r <- acf(x, lag.max = min(500, n - 2), plot = FALSE)$acf[-1]
  k <- which(r < 0.05)[1]; if (is.na(k)) k <- length(r)
  max(1, n / (1 + 2 * sum(r[1:max(1, k - 1)]))) }

d <- stack()
out <- list()
for (blk in c("mu", "et")) {
  M <- d[[blk]]; labs <- if (blk == "mu") labMu else labEta
  R <- apply(M, 2, function(c) splitRhat(rk(c))); S <- apply(M, 2, ess)
  cat(sprintf("\n%-4s n=%4d | Rhat med %.4f max %.4f | <1.01: %2d/%2d | ESS med %5.0f min %4.0f\n",
      toupper(blk), nrow(M), median(R, na.rm = TRUE), max(R, na.rm = TRUE),
      sum(R < 1.01, na.rm = TRUE), length(R), median(S), min(S)))
  o <- order(-R)[1:5]
  cat("  worst five:\n")
  for (i in o) cat(sprintf("    %-26s Rhat %.3f  ESS %4.0f\n", labs[i], R[i], S[i]))
  out[[blk]] <- data.frame(block = blk, par = labs, Rhat = round(R, 4), ESS = round(S))
}
write.csv(do.call(rbind, out), file.path(SC, "bayes_cumulative_diag.csv"), row.names = FALSE)
cat("\nwritten bayes_cumulative_diag.csv\n")
