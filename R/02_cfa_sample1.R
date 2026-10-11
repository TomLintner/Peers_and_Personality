## 02_cfa_sample1.R
## source: RSiena_build/scripts/gacr_cfa2.R
source("00_setup.R")

suppressMessages(library(lavaan))
SC <- WORK
d <- readRDS(file.path(OBJ, "gacr.rds"))

num <- function(x) suppressWarnings(as.numeric(ifelse(x == "", NA, x)))
rnd <- function(x) { x[!is.na(x) & x == 0] <- NA; floor(x + 0.5) }
A <- sapply(paste0("abfi", 1:44), function(v) rnd(num(d[[v]])))
B <- sapply(paste0("bbfi", 1:44), function(v) rnd(num(d[[v]])))
colnames(A) <- colnames(B) <- paste0("i", 1:44)
key <- list(E = c(1, -6, 11, 16, -21, 26, -31, 36),
            A = c(-2, 7, -12, 17, 22, -27, 32, -37, 42),
            C = c(3, -8, 13, -18, -23, 28, 33, 38, -43),
            N = c(4, -9, 14, 19, -24, 29, -34, 39),
            O = c(5, 10, 15, 20, 25, 30, -35, 40, -41, 44))
recode <- function(M) { for (tr in names(key)) for (k in key[[tr]]) if (k < 0) M[, abs(k)] <- 6 - M[, abs(k)]; M }
A <- recode(A); B <- recode(B)

serp <- function(items, cnt) { a <- integer(length(items))
  for (i in seq_along(items)) { p <- which.min(cnt); a[i] <- p; cnt[p] <- cnt[p] + 1 }; list(a = a, cnt = cnt) }
pd <- list()
for (tr in names(key)) {
  k <- key[[tr]]; idx <- abs(k); fwd <- k > 0; X <- A[, idx, drop = FALSE]
  ir <- sapply(seq_along(idx), function(j) cor(X[, j], rowMeans(X[, -j, drop = FALSE], na.rm = TRUE), use = "pairwise"))
  itF <- idx[fwd][order(ir[fwd], decreasing = TRUE)]; itR <- idx[!fwd][order(ir[!fwd], decreasing = TRUE)]
  aF <- serp(itF, c(0, 0, 0)); aR <- serp(itR, aF$cnt)
  pd[[tr]] <- lapply(1:3, function(p) c(itF[aF$a == p], itR[aR$a == p]))
}
mk <- function(M, tag) { o <- list()
  for (tr in names(key)) for (p in 1:3) { it <- pd[[tr]][[p]]; X <- M[, it, drop = FALSE]
    v <- rowMeans(X, na.rm = TRUE); v[rowSums(!is.na(X)) < ceiling(length(it) / 2)] <- NA
    o[[paste0(tr, "p", p, "_", tag)]] <- v }
  as.data.frame(o) }
P <- cbind(mk(A, "1"), mk(B, "2")); P$clid <- d$clid
P <- P[rowSums(!is.na(P[, grep("p[123]_", names(P))])) > 0, ]
cat("rows in measurement model:", nrow(P), "\n")

cat("\n=== observed parcel retest correlations (T1-T2) ===\n")
for (tr in names(key)) cat(sprintf("  %s: %s\n", tr, paste(sprintf("%.2f",
  sapply(1:3, function(p) cor(P[[paste0(tr,"p",p,"_1")]], P[[paste0(tr,"p",p,"_2")]], use = "pairwise"))), collapse = "  ")))

## marker-variable identification: parcel 1 loading fixed to 1 in both waves
syn <- function(tr, lv) {
  v1 <- paste0(tr, "p", 1:3, "_1"); v2 <- paste0(tr, "p", 1:3, "_2")
  L1 <- if (lv == "configural") c("", "la2*", "la3*") else c("", "l2*", "l3*")
  L2 <- c("", "l2*", "l3*")
  s <- c(paste0("f1 =~ ", paste0(L1, v1, collapse = " + ")),
         paste0("f2 =~ ", paste0(L2, v2, collapse = " + ")),
         paste0(v1, " ~~ ", v2))
  if (lv == "scalar") s <- c(s, paste0(v1, " ~ t", 1:3, "*1"), paste0(v2, " ~ t", 1:3, "*1"),
                             "f1 ~ 0*1", "f2 ~ 1")
  s
}
cols <- c("chisq", "df", "cfi", "tli", "rmsea", "srmr")
res <- list(); fits <- list()
for (tr in names(key)) for (lv in c("configural", "metric", "scalar")) {
  f <- try(cfa(paste(syn(tr, lv), collapse = "\n"), data = P, missing = "fiml", estimator = "MLR"), silent = TRUE)
  fits[[paste(tr, lv)]] <- f
  m <- if (inherits(f, "try-error") || !lavInspect(f, "converged")) rep(NA, 6) else
    fitMeasures(f, c("chisq.scaled", "df.scaled", "cfi.robust", "tli.robust", "rmsea.robust", "srmr"))
  res[[length(res) + 1]] <- data.frame(trait = tr, level = lv, setNames(as.list(round(as.numeric(m), 3)), cols))
}
R <- do.call(rbind, res); R$dCFI <- NA; R$dRMSEA <- NA
for (i in seq_len(nrow(R))) if (R$level[i] != "configural") {
  R$dCFI[i] <- round(R$cfi[i] - R$cfi[i - 1], 3); R$dRMSEA[i] <- round(R$rmsea[i] - R$rmsea[i - 1], 3) }
cat("\n=== GACR longitudinal invariance, per trait (MLR, FIML, n =", nrow(P), ") ===\n")
print(R, row.names = FALSE)
cat("\ninvariance supported if dCFI >= -0.010 and dRMSEA <= 0.015\n")

cat("\n=== latent mean change at T2 (T1 fixed at 0), scalar models ===\n")
for (tr in names(key)) {
  f <- fits[[paste(tr, "scalar")]]
  pe <- parameterEstimates(f); r <- pe[pe$lhs == "f2" & pe$op == "~1", ]
  vv <- pe[pe$lhs == "f1" & pe$op == "~~" & pe$rhs == "f1", "est"]
  cat(sprintf("  %s  mean_T2 = %+.3f (SE %.3f, p = %.3f)  = %+.2f SD_T1\n",
      tr, r$est, r$se, r$pvalue, r$est / sqrt(vv)))
}
saveRDS(list(P = P, fits = fits, R = R, pd = pd), file.path(OBJ, "gacr_cfa2.rds"))
