## 06_build_sample2.R
## source: RSiena_build/scripts/build_oz.R
source("00_setup.R")

suppressMessages(library(RSiena))
setwd(RAW)
SC <- WORK
tr5 <- c("E", "A", "C", "N", "O")

b <- read.csv("OZ_baseline_survey_data.csv", check.names = FALSE, colClasses = "character")
f <- read.csv("OZ_followup_survey_data.csv", check.names = FALSE, colClasses = "character")
bf <- read.csv(file.path(SC, "oz_bfi_bins.csv"), colClasses = "character")
for (v in grep("_w[12]$", names(bf), value = TRUE)) bf[[v]] <- as.integer(bf[[v]])
for (v in grep("_sc[12]$", names(bf), value = TRUE)) bf[[v]] <- as.numeric(bf[[v]])
bf$sexBirthFemale <- as.numeric(bf$sexBirthFemale); bf$genderFemale <- as.numeric(bf$genderFemale)
rowOf <- setNames(seq_len(nrow(bf)), bf$id)

ros <- function(dat) sub("^peerReportFriends_", "", grep("^peerReportFriends_", names(dat), value = TRUE))
r1 <- ros(b); r2 <- ros(f)
classes <- sort(unique(c(b$idClassroom, f$idClassroom)))

DATA <- list(); DESC <- list()
for (cl in classes) {
  ro1 <- r1[substr(r1, 1, 4) == cl]; ro2 <- r2[substr(r2, 1, 4) == cl]
  rp1 <- b$idStudent_final[b$idClassroom == cl]; rp2 <- f$idStudent_final[f$idClassroom == cl]
  nodes <- sort(union(union(ro1, ro2), union(rp1, rp2))); n <- length(nodes)

  build <- function(dat, roster, resp) {
    M <- matrix(NA_real_, n, n, dimnames = list(nodes, nodes))
    ix <- match(resp, dat$idStudent_final)
    V  <- suppressWarnings(matrix(as.numeric(as.matrix(dat[ix, paste0("peerReportFriends_", roster)])),
                                 nrow = length(ix)))
    M[match(resp, nodes), match(roster, nodes)] <- V
    inros <- nodes %in% roster; present <- inros | nodes %in% resp
    M[, !inros] <- 10                 # not on this wave's roster -> nobody could nominate them
    M[!present, ] <- 10               # not in the class at this wave at all
    diag(M) <- 0                      # self-nominations disregarded
    list(M = M, inros = inros, present = present)
  }
  c1 <- build(b, ro1, rp1); c2 <- build(f, ro2, rp2); M1 <- c1$M; M2 <- c2$M
  stopifnot(all(M1 %in% c(0, 1, 10, NA)), all(M2 %in% c(0, 1, 10, NA)))
  stopifnot(all(diag(M1) == 0), all(diag(M2) == 0))

  net <- sienaDependent(array(c(M1, M2), dim = c(n, n, 2)))
  ii <- rowOf[nodes]
  pick <- function(col) { v <- rep(NA_real_, n); ok <- !is.na(ii); v[ok] <- bf[[col]][ii[ok]]; v }
  beh  <- lapply(setNames(tr5, tr5), function(t) cbind(as.integer(pick(paste0(t, "_w1"))),
                                                       as.integer(pick(paste0(t, "_w2")))))
  t1sc <- lapply(setNames(tr5, tr5), function(t) pick(paste0(t, "_sc1")))
  DATA[[cl]] <- list(net = net, beh = beh, t1sc = t1sc, n = n, nodes = nodes,
                     genderF = pick("genderFemale"), sexF = pick("sexBirthFemale"))

  tR <- function(m) { m[m == 10] <- NA; m }
  a <- tR(M1); bb <- tR(M2); both <- !is.na(a) & !is.na(bb)
  DESC[[cl]] <- data.frame(class = cl, n = n,
    densT1 = round(mean(a[row(a) != col(a)] == 1, na.rm = TRUE), 3),
    densT2 = round(mean(bb[row(bb) != col(bb)] == 1, na.rm = TRUE), 3),
    missT1 = round(mean(is.na(a[row(a) != col(a)])), 2),
    missT2 = round(mean(is.na(bb[row(bb) != col(bb)])), 2),
    jaccard = round(sum(a[both] == 1 & bb[both] == 1) / sum(a[both] == 1 | bb[both] == 1), 2),
    offRosT1 = sum(!c1$inros), offRosT2 = sum(!c2$inros),
    noGender = sum(is.na(DATA[[cl]]$genderF)),
    behBoth = sum(!is.na(beh$N[, 1]) & !is.na(beh$N[, 2])), stringsAsFactors = FALSE)
}
D <- do.call(rbind, DESC)
cat("=== OZ per-class descriptives (friendship network) ===\n"); print(D, row.names = FALSE)
cat("\ntotals: actors", sum(D$n), "| mean density T1", round(mean(D$densT1), 3),
    "T2", round(mean(D$densT2), 3), "| mean Jaccard", round(mean(D$jaccard), 2),
    "| min", min(D$jaccard), "\n")
cat("classes with Jaccard < 0.30:", sum(D$jaccard < 0.30), " < 0.20:", sum(D$jaccard < 0.20), "\n")
cat("actors with no gender:", sum(D$noGender), "of", sum(D$n), "\n")

cat("\n=== behaviour: category occupancy pooled over classes ===\n")
for (t in tr5) {
  v <- unlist(lapply(DATA, function(z) c(z$beh[[t]])))
  cat(sprintf("  %s  %s   missing both waves: %d\n", t,
      paste(sprintf("%d:%d", 1:5, table(factor(v, 1:5))), collapse = "  "),
      sum(unlist(lapply(DATA, function(z) is.na(z$beh[[t]][, 1]) & is.na(z$beh[[t]][, 2]))))))
  nv <- sapply(DATA, function(z) length(unique(na.omit(c(z$beh[[t]])))))
  if (any(nv < 3)) cat("     classes with <3 distinct values:", names(nv)[nv < 3], "\n")
}

groups <- list()
for (t in tr5) {
  objs <- list()
  for (cl in classes) {
    z <- DATA[[cl]]
    args <- list(net = z$net,
                 beh = sienaDependent(z$beh[[t]], type = "behavior", allowOnly = FALSE),
                 genderF = coCovar(z$genderF), sexF = coCovar(z$sexF))
    for (o in setdiff(tr5, t)) args[[paste0("cov", o)]] <- coCovar(z$t1sc[[o]])
    objs[[cl]] <- do.call(sienaDataCreate, args)
  }
  groups[[t]] <- sienaGroupCreate(objs)
}
saveRDS(list(groups = groups, desc = D, data = DATA, classes = classes), file.path(OBJ, "oz_siena.rds"))
cat("\nsaved oz_siena.rds -", length(classes), "groups per trait\n")
print(groups[["N"]])
