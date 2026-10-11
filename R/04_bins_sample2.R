## 04_bins_sample2.R
## source: RSiena_build/scripts/oz_bins.R
source("00_setup.R")

suppressMessages(library(lavaan))
setwd(RAW)
SC <- WORK
num <- function(x) suppressWarnings(as.numeric(ifelse(x == "" | x == "-oth-", NA, x)))
b <- read.csv("OZ_baseline_survey_data.csv", check.names = FALSE, colClasses = "character")
f <- read.csv("OZ_followup_survey_data.csv", check.names = FALSE, colClasses = "character")   # all 436 rows: the 5 mislabelled are real follow-ups

trait <- list(E = c("E2", "E8"), A = c("A1", "A5"), C = c("C1", "C5"),
              N = c("N2", "N8"), O = c("O4", "O8", "O9"))
revit <- c("E2", "N2", "O9", "C5", "A1")          # empirically established keying
items <- unlist(trait, use.names = FALSE)

ids <- union(b$idStudent_final, f$idStudent_final)
D <- data.frame(id = ids, stringsAsFactors = FALSE)
ib <- match(ids, b$idStudent_final); jf <- match(ids, f$idStudent_final)
D$idClassroom <- ifelse(!is.na(ib), b$idClassroom[ib], f$idClassroom[jf])
D$sexBirthFemale <- ifelse(!is.na(ib), num(b$sexBirthFemale)[ib], num(f$sexBirthFemale)[jf])
gb <- num(b$genderFemale); gf <- num(f$genderFemale)                     # '-oth-' -> NA by num()
D$genderFemale <- ifelse(!is.na(ib) & !is.na(gb[ib]), gb[ib], gf[jf])
for (v in items) {
  x1 <- num(b[[paste0("bfi_bfi", v)]])[ib]; x2 <- num(f[[paste0("bfi_bfi", v)]])[jf]
  if (v %in% revit) { x1 <- 6 - x1; x2 <- 6 - x2 }
  D[[paste0(v, "_1")]] <- x1; D[[paste0(v, "_2")]] <- x2
}
cat("OZ persons:", nrow(D), " classrooms:", length(unique(D$idClassroom)), "\n")
cat("with any BFI T1:", sum(rowSums(!is.na(D[paste0(items, "_1")])) > 0),
    " T2:", sum(rowSums(!is.na(D[paste0(items, "_2")])) > 0), "\n")
cat("gender: sexBirthFemale missing", sum(is.na(D$sexBirthFemale)),
    "| genderFemale missing (incl. '-oth-' as NA)", sum(is.na(D$genderFemale)), "\n")

syn <- function(tr, lv) {
  v <- trait[[tr]]; k <- length(v)
  v1 <- paste0(v, "_1"); v2 <- paste0(v, "_2")
  Lc <- if (k == 3) c("", "l2*", "l3*") else c("", "l2*")
  L1 <- if (lv == "configural") sub("^l", "la", Lc) else Lc
  s <- c(paste0("f1 =~ ", paste0(L1, v1, collapse = " + ")),
         paste0("f2 =~ ", paste0(Lc, v2, collapse = " + ")))
  if (k == 3) s <- c(s, paste0(v1, " ~~ ", v2))      # residual covariances only identified with 3 indicators
  if (lv == "scalar") s <- c(s, paste0(v1, " ~ t", 1:k, "*1"), paste0(v2, " ~ t", 1:k, "*1"),
                             "f1 ~ 0*1", "f2 ~ 1")
  s
}
cols <- c("chisq", "df", "cfi", "rmsea", "srmr")
res <- list(); fits <- list()
for (tr in names(trait)) {
  lvs <- if (length(trait[[tr]]) == 3) c("configural", "metric", "scalar") else c("metric", "scalar")
  for (lv in lvs) {
    fit <- try(cfa(paste(syn(tr, lv), collapse = "\n"), data = D, missing = "fiml", estimator = "MLR"), silent = TRUE)
    fits[[paste(tr, lv)]] <- fit
    ok <- !inherits(fit, "try-error") && lavInspect(fit, "converged")
    m <- if (!ok) rep(NA, 5) else fitMeasures(fit, c("chisq.scaled", "df.scaled", "cfi.robust", "rmsea.robust", "srmr"))
    res[[length(res) + 1]] <- data.frame(trait = tr, k = length(trait[[tr]]), level = lv,
                                         setNames(as.list(round(as.numeric(m), 3)), cols))
  }
}
R <- do.call(rbind, res)
cat("\n=== OZ longitudinal CFA per trait (2-item traits: no residual covariances, not identified) ===\n")
print(R, row.names = FALSE)

getsc <- function(fit, tr) {
  v <- trait[[tr]]
  s <- lavPredict(fit, method = "regression"); ci <- lavInspect(fit, "case.idx")
  o <- matrix(NA_real_, nrow(D), 2, dimnames = list(NULL, c("f1", "f2")))
  o[ci, ] <- as.matrix(s)[, c("f1", "f2")]
  o[rowSums(!is.na(D[paste0(v, "_1")])) == 0, 1] <- NA
  o[rowSums(!is.na(D[paste0(v, "_2")])) == 0, 2] <- NA
  o
}
cat("\n=== factor scores vs simple item means (do the scores add anything?) ===\n")
FS <- list()
for (tr in names(trait)) {
  fit <- fits[[paste(tr, "scalar")]]; FS[[tr]] <- getsc(fit, tr)
  v <- trait[[tr]]
  m1 <- rowMeans(D[paste0(v, "_1")], na.rm = TRUE); m2 <- rowMeans(D[paste0(v, "_2")], na.rm = TRUE)
  lc <- lavInspect(fit, "cor.lv")["f1", "f2"]
  cat(sprintf("  %s  r(score, item mean) T1=%.3f T2=%.3f | latent r(T1,T2)=%.2f | distinct score values=%d\n",
      tr, cor(FS[[tr]][, 1], m1, use = "pairwise"), cor(FS[[tr]][, 2], m2, use = "pairwise"),
      lc, length(unique(round(c(FS[[tr]]), 6)[!is.na(c(FS[[tr]]))]))))
}

cat("\n=== OZ 5-category behaviour variables (pooled quintile cuts) ===\n")
BIN <- list(); CUT <- list()
for (tr in names(trait)) {
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
out <- D[, c("id", "idClassroom", "sexBirthFemale", "genderFemale")]
for (tr in names(trait)) { out[[paste0(tr, "_w1")]] <- BIN[[tr]][, 1]; out[[paste0(tr, "_w2")]] <- BIN[[tr]][, 2]
  out[[paste0(tr, "_sc1")]] <- round(FS[[tr]][, 1], 4); out[[paste0(tr, "_sc2")]] <- round(FS[[tr]][, 2], 4) }
write.csv(out, file.path(SC, "oz_bfi_bins.csv"), row.names = FALSE)
saveRDS(list(scores = FS, bins = BIN, cuts = CUT, D = D, R = R), file.path(OBJ, "oz_bins.rds"))
cat("\nwritten oz_bfi_bins.csv:", nrow(out), "rows\n")
