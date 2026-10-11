## 03_build_sample1.R
## source: RSiena_build/scripts/build_gacr.R
source("00_setup.R")

suppressMessages(library(RSiena))
SC  <- WORK
MD  <- file.path(RAW, "matrices + matrix format/Matrices without duplicities")
tr5 <- c("E", "A", "C", "N", "O")

bf <- read.csv(file.path(SC, "gacr_bfi_bins.csv"), colClasses = "character")
for (v in grep("_w[12]$", names(bf), value = TRUE)) bf[[v]] <- as.integer(bf[[v]])
for (v in grep("_sc[12]$", names(bf), value = TRUE)) bf[[v]] <- as.numeric(bf[[v]])
bf$sex <- as.integer(bf$sex)
byA <- setNames(seq_len(nrow(bf)), bf$aid); byB <- setNames(seq_len(nrow(bf)), bf$bid)

## sex comes from the full survey file (complete for all 910, incl. actors with no BFI)
gg <- readRDS(file.path(OBJ, "gacr.rds"))
sexByA <- setNames(as.integer(gg$sex), trimws(gg$aid))
sexByB <- setNames(as.integer(gg$sex), trimws(gg$bid))

classes <- sort(unique(substr(sub("-output-.*$", "", list.files(MD, pattern = "-bn1[.]csv$")), 1, 4)))
cat("two-wave classes:", length(classes), "\n")

readnet <- function(cl, wave) {
  it <- if (wave == 1) "an1" else "bn1"
  x <- read.csv(file.path(MD, paste0(cl, sprintf("%02d", wave), "-output-", it, ".csv")),
                check.names = FALSE, row.names = 1)
  m <- as.matrix(x); stopifnot(identical(rownames(m), colnames(m)))
  pn <- substr(sprintf("%08s", rownames(m)), 7, 8)          # person number within class
  o  <- order(pn); m <- m[o, o, drop = FALSE]
  list(m = m, ids = rownames(m), pn = pn[o])
}

DATA <- list(); DESC <- list()
for (cl in classes) {
  w1 <- readnet(cl, 1); w2 <- readnet(cl, 2)
  stopifnot(identical(w1$pn, w2$pn))                        # same actors, same order
  n <- length(w1$pn)

  clean <- function(m) {
    nonmem <- apply(m, 2, function(cc) any(cc == 10, na.rm = TRUE))   # not a class member this wave
    m[m %in% c(-1, -9)] <- NA                               # row-wise non-response -> missing
    if (any(nonmem)) { m[nonmem, ] <- 10; m[, nonmem] <- 10 }         # RSiena structural zeros
    diag(m) <- 0                                            # self-nominations disregarded
    list(m = m, nonmem = nonmem)
  }
  c1 <- clean(w1$m); c2 <- clean(w2$m); M1 <- c1$m; M2 <- c2$m
  stopifnot(all(M1 %in% c(0, 1, 10, NA)), all(M2 %in% c(0, 1, 10, NA)))
  stopifnot(all(diag(M1) == 0), all(diag(M2) == 0))
  for (j in which(c1$nonmem)) stopifnot(all(M1[j, -j] == 10), all(M1[-j, j] == 10))
  for (j in which(c2$nonmem)) stopifnot(all(M2[j, -j] == 10), all(M2[-j, j] == 10))
  net <- sienaDependent(array(c(M1, M2), dim = c(n, n, 2)))

  i1 <- byA[w1$ids]; i2 <- byB[w2$ids]                      # rows in the BFI file (NA if absent)
  pick <- function(idx, col) { v <- rep(NA_real_, n); ok <- !is.na(idx)
                               v[ok] <- as.numeric(bf[[col]])[idx[ok]]; v }
  sex <- unname(sexByA[w1$ids]); miss <- is.na(sex)
  sex[miss] <- unname(sexByB[w2$ids])[miss]

  beh <- lapply(setNames(tr5, tr5), function(t)
    cbind(as.integer(pick(i1, paste0(t, "_w1"))), as.integer(pick(i2, paste0(t, "_w2")))))
  # control covariates: T1 factor score, falling back to T2 when T1 is missing
  t1sc <- lapply(setNames(tr5, tr5), function(t) {
    v <- pick(i1, paste0(t, "_sc1")); m <- is.na(v)
    v[m] <- pick(i2, paste0(t, "_sc2"))[m]; v })

  DATA[[cl]] <- list(net = net, beh = beh, sex = sex, t1sc = t1sc, n = n, ids1 = w1$ids, ids2 = w2$ids)

  tieR <- function(m) { m[m == 10] <- NA; m }
  a <- tieR(M1); b <- tieR(M2); both <- !is.na(a) & !is.na(b)
  jac <- sum(a[both] == 1 & b[both] == 1) / sum(a[both] == 1 | b[both] == 1)
  DESC[[cl]] <- data.frame(class = cl, n = n,
    densT1 = round(mean(a[row(a) != col(a)] == 1, na.rm = TRUE), 3),
    densT2 = round(mean(b[row(b) != col(b)] == 1, na.rm = TRUE), 3),
    missT1 = round(mean(is.na(a[row(a) != col(a)])), 2),
    missT2 = round(mean(is.na(b[row(b) != col(b)])), 2),
    jaccard = round(jac, 2), nonMemT1 = sum(c1$nonmem), nonMemT2 = sum(c2$nonmem),
    behBoth = sum(!is.na(beh$N[, 1]) & !is.na(beh$N[, 2])), stringsAsFactors = FALSE)
}
D <- do.call(rbind, DESC)
cat("\n=== GACR per-class descriptives (like-most network) ===\n"); print(D, row.names = FALSE)
cat("\ntotals: actors", sum(D$n), "| mean density T1", round(mean(D$densT1), 3),
    "T2", round(mean(D$densT2), 3), "| mean Jaccard", round(mean(D$jaccard), 2),
    "| min Jaccard", min(D$jaccard), "\n")
cat("classes with Jaccard < 0.30:", sum(D$jaccard < 0.30), " < 0.20:", sum(D$jaccard < 0.20), "\n")

## behaviour sanity per class
cat("\n=== behaviour variability per class (n distinct values, both waves pooled) ===\n")
for (t in tr5) {
  nv <- sapply(DATA, function(z) length(unique(na.omit(c(z$beh[[t]])))))
  nm <- sapply(DATA, function(z) sum(is.na(z$beh[[t]][, 1]) & is.na(z$beh[[t]][, 2])))
  cat(sprintf("  %s  min distinct=%d  classes with <3 distinct=%d  actors missing both waves: %d of %d\n",
              t, min(nv), sum(nv < 3), sum(nm), sum(D$n)))
}

## ---- one sienaGroup per trait ----
groups <- list()
for (t in tr5) {
  objs <- list()
  for (cl in classes) {
    z <- DATA[[cl]]
    args <- list(net = z$net,
                 beh = sienaDependent(z$beh[[t]], type = "behavior", allowOnly = FALSE),
                 sex = coCovar(as.numeric(z$sex)))
    for (o in setdiff(tr5, t)) args[[paste0("cov", o)]] <- coCovar(as.numeric(z$t1sc[[o]]))
    objs[[cl]] <- do.call(sienaDataCreate, args)
  }
  groups[[t]] <- sienaGroupCreate(objs)
  cat("\n--- sienaGroup for", t, "---\n")
  print(c(groups = length(objs), actors = sum(sapply(objs, function(o) o$observations * 0 + 1))))
}
saveRDS(list(groups = groups, desc = D, data = DATA, classes = classes),
        file.path(OBJ, "gacr_siena.rds"))
cat("\nsaved gacr_siena.rds\n")
cat("\n=== example: print of the E group ===\n")
print(groups[["E"]])
