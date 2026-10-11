## 01_bins_sample1.R
## source: RSiena_build/scripts/gacr_bins2.R
source("00_setup.R")

suppressMessages(library(lavaan))
SC <- WORK
z <- readRDS(file.path(OBJ, "gacr_cfa2.rds")); P <- z$P; fits <- z$fits
d <- readRDS(file.path(OBJ, "gacr.rds")); tr5 <- c("E", "A", "C", "N", "O")

cat("=== score test on the scalar equality constraints (A, C) ===\n")
for (tr in c("A", "C")) {
  st <- lavTestScore(fits[[paste(tr, "scalar")]])
  cat("\n", tr, " overall:", sprintf("chisq=%.2f df=%.0f p=%.3f", st$test$X2, st$test$df, st$test$p.value), "\n")
  print(st$uni, row.names = FALSE)
}

getsc <- function(f, tr) {
  s <- lavPredict(f, method = "regression"); ci <- lavInspect(f, "case.idx")
  o <- matrix(NA_real_, nrow(P), 2, dimnames = list(NULL, c("f1", "f2")))
  o[ci, ] <- as.matrix(s)[, c("f1", "f2")]
  # a wave with no observed parcel data must stay missing
  o[rowSums(!is.na(P[, paste0(tr, "p", 1:3, "_1")])) == 0, 1] <- NA
  o[rowSums(!is.na(P[, paste0(tr, "p", 1:3, "_2")])) == 0, 2] <- NA
  o
}
FS <- lapply(setNames(tr5, tr5), function(tr) getsc(fits[[paste(tr, "scalar")]], tr))

cat("\n=== factor scores after masking waves with no data ===\n")
for (tr in tr5) {
  f <- fits[[paste(tr, "scalar")]]
  lc <- lavInspect(f, "cor.lv")["f1", "f2"]
  s1 <- FS[[tr]][, 1]; s2 <- FS[[tr]][, 2]; ok <- !is.na(s1) & !is.na(s2)
  cat(sprintf("  %s  n_T1=%3d n_T2=%3d both=%3d | latent r(T1,T2)=%.2f | score r (complete)=%.2f | SD %.2f/%.2f\n",
      tr, sum(!is.na(s1)), sum(!is.na(s2)), sum(ok), lc, cor(s1[ok], s2[ok]),
      sd(s1, na.rm = TRUE), sd(s2, na.rm = TRUE)))
}

cat("\n=== 5-category behaviour variables (one set of pooled quintile cuts, both waves) ===\n")
BIN <- list(); CUT <- list()
for (tr in tr5) {
  s1 <- FS[[tr]][, 1]; s2 <- FS[[tr]][, 2]
  cuts <- as.numeric(quantile(c(s1, s2), c(.2, .4, .6, .8), na.rm = TRUE)); CUT[[tr]] <- cuts
  b1 <- as.integer(cut(s1, c(-Inf, cuts, Inf))); b2 <- as.integer(cut(s2, c(-Inf, cuts, Inf)))
  BIN[[tr]] <- cbind(w1 = b1, w2 = b2); ok <- !is.na(b1) & !is.na(b2)
  cat(sprintf("\n%s  cuts: %s\n", tr, paste(sprintf("%+.3f", cuts), collapse = "  ")))
  cat("   T1 per bin:", sprintf("%3d", table(factor(b1, 1:5))),
      "\n   T2 per bin:", sprintf("%3d", table(factor(b2, 1:5))), "\n")
  cat(sprintf("   both waves n=%d | unchanged %.0f%% | move 1 %.0f%% | move >=2 %.0f%% (max %d) | r=%.2f\n",
      sum(ok), 100 * mean(b1[ok] == b2[ok]), 100 * mean(abs(b1[ok] - b2[ok]) == 1),
      100 * mean(abs(b1[ok] - b2[ok]) >= 2), max(abs(b1[ok] - b2[ok])), cor(b1[ok], b2[ok])))
}

ridx <- as.integer(rownames(P))
out <- data.frame(aid = d$aid[ridx], bid = d$bid[ridx], clid = d$clid[ridx],
                  part = d$part[ridx], sex = d$sex[ridx], stringsAsFactors = FALSE)
for (tr in tr5) { out[[paste0(tr, "_w1")]] <- BIN[[tr]][, 1]; out[[paste0(tr, "_w2")]] <- BIN[[tr]][, 2]
  out[[paste0(tr, "_sc1")]] <- round(FS[[tr]][, 1], 4); out[[paste0(tr, "_sc2")]] <- round(FS[[tr]][, 2], 4) }
write.csv(out, file.path(SC, "gacr_bfi_bins.csv"), row.names = FALSE)
saveRDS(list(scores = FS, bins = BIN, cuts = CUT, ids = out[, 1:5]), file.path(OBJ, "gacr_bins.rds"))

mat <- file.path(RAW, "matrices + matrix format/Matrices without duplicities")
w2 <- unique(substr(sub("-output-.*$", "", list.files(mat, pattern = "-bn1[.]csv$")), 1, 4))
sel <- sprintf("%04d", as.integer(out$clid)) %in% w2
cat("\n=== coverage inside the 27 two-wave matrix classes (", sum(sel), "students ) ===\n")
for (tr in tr5) cat(sprintf("  %s  T1=%3d  T2=%3d  both=%3d\n", tr,
  sum(!is.na(out[[paste0(tr, "_w1")]][sel])), sum(!is.na(out[[paste0(tr, "_w2")]][sel])),
  sum(!is.na(out[[paste0(tr, "_w1")]][sel]) & !is.na(out[[paste0(tr, "_w2")]][sel]))))
cat("\nwritten gacr_bfi_bins.csv:", nrow(out), "rows\n")
