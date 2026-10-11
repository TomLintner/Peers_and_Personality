## 05_bins_sample2_refine.R
## source: RSiena_build/scripts/oz_bins2.R
source("00_setup.R")

SC <- WORK
z <- readRDS(file.path(OBJ, "oz_bins.rds")); D <- z$D
trait <- list(E = c("E2", "E8"), A = c("A1", "A5"), C = c("C1", "C5"),
              N = c("N2", "N8"), O = c("O4", "O8", "O9"))
rnd <- function(x) ifelse(is.na(x), NA_integer_, as.integer(pmin(5, pmax(1, floor(x + 0.5)))))
sc <- function(v, w) { X <- D[paste0(v, "_", w)]; m <- rowMeans(X, na.rm = TRUE)
  m[rowSums(!is.na(X)) == 0] <- NA; m }

stats <- function(b1, b2, lab) {
  ok <- !is.na(b1) & !is.na(b2)
  cat(sprintf("  %-18s T1:%s  T2:%s  n=%d unchanged %2.0f%% move1 %2.0f%% move>=2 %2.0f%% (max %d) r=%.2f\n",
      lab, paste(sprintf("%3d", table(factor(b1, 1:5))), collapse = ""),
      paste(sprintf("%3d", table(factor(b2, 1:5))), collapse = ""),
      sum(ok), 100 * mean(b1[ok] == b2[ok]), 100 * mean(abs(b1[ok] - b2[ok]) == 1),
      100 * mean(abs(b1[ok] - b2[ok]) >= 2), max(abs(b1[ok] - b2[ok])), cor(b1[ok], b2[ok])))
}
cat("=== OZ: rounded composite vs best single item (bins are integers 1-5, fixed cuts) ===\n")
for (tr in names(trait)) {
  v <- trait[[tr]]
  stats(rnd(sc(v, 1)), rnd(sc(v, 2)), paste0(tr, " all ", length(v), " items"))
  for (one in v) stats(rnd(sc(one, 1)), rnd(sc(one, 2)), paste0(tr, " ", one, " alone"))
  cat("\n")
}

## final choice: composite for N, O, C; single item E8 for E; single item A1(rev) for A
chosen <- list(E = "E8", A = "A1", C = c("C1", "C5"), N = c("N2", "N8"), O = c("O4", "O8", "O9"))
cat("=== chosen OZ behaviour variables ===\n")
out <- D[, c("id", "idClassroom", "sexBirthFemale", "genderFemale")]
for (tr in names(chosen)) {
  v <- chosen[[tr]]
  b1 <- rnd(sc(v, 1)); b2 <- rnd(sc(v, 2))
  stats(b1, b2, paste0(tr, " = ", paste(v, collapse = "+")))
  out[[paste0(tr, "_w1")]] <- b1; out[[paste0(tr, "_w2")]] <- b2
  out[[paste0(tr, "_sc1")]] <- round(sc(v, 1), 4); out[[paste0(tr, "_sc2")]] <- round(sc(v, 2), 4)
}
write.csv(out, file.path(SC, "oz_bfi_bins.csv"), row.names = FALSE)
cat("\nwritten oz_bfi_bins.csv:", nrow(out), "rows\n")

cat("\n=== per-classroom coverage (chosen variables, N as example) ===\n")
tb <- table(out$idClassroom, !is.na(out$N_w1) & !is.na(out$N_w2))
print(cbind(n = rowSums(tb), bothWavesBFI = tb[, "TRUE"]))
