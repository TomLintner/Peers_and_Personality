## 40_assortativity.R
## source: Results/scripts/homophily.R
source("00_setup.R")

## Wave-1 personality homophily, per classroom and per dataset.
##
## Measure choice. The trait variables are 5-category ordinal scores, which is
## how RSiena treats them, so a categorical index (E-I, Coleman, Freeman) would
## throw away the ordering and require an arbitrary dichotomisation. The right
## family is network autocorrelation of a scalar attribute:
##   PRIMARY  Newman's assortativity coefficient r -- the Pearson correlation
##            between sender and receiver trait values taken over the observed
##            directed ties. Bounded [-1, 1], 0 = random mixing, directly
##            comparable across classrooms of different size and density.
##   CHECK    Moran's I with the adjacency matrix as weights (the classic
##            network autocorrelation statistic).
##   MODEL    simGap: mean RSiena similarity over tied dyads minus the mean over
##            all dyads. This is exactly the quantity the simX effect acts on,
##            so it links the descriptive result to the SAOM.
## Significance comes from permuting trait values over actors within each
## classroom (1,000 permutations), which holds the network fixed and destroys
## only the attribute-tie association.
suppressMessages(library(RSienaTest))
set.seed(2012)
NPERM <- 1000
tr <- c("E", "A", "C", "N", "O")

metrics <- function(A, v) {
  ## A: binary adjacency (wave 1), NA where unobserved. v: trait, NA allowed.
  ok <- !is.na(v)
  if (sum(ok) < 5) return(NULL)
  A <- A[ok, ok, drop = FALSE]; v <- v[ok]
  n <- length(v)
  if (sd(v) == 0) return(NULL)
  diag(A) <- NA
  E <- which(A == 1, arr.ind = TRUE)
  m <- nrow(E)
  if (m < 5) return(NULL)
  ## Newman assortativity: correlation of (sender, receiver) values over edges
  r <- suppressWarnings(cor(v[E[, 1]], v[E[, 2]]))
  ## Moran's I with adjacency weights
  W <- A; W[is.na(W)] <- 0
  d <- v - mean(v); S0 <- sum(W)
  I <- if (S0 > 0) (n / S0) * sum(W * outer(d, d)) / sum(d^2) else NA
  ## RSiena similarity: 1 - |vi - vj| / range, over ties vs over all dyads
  rng <- diff(range(v)); if (rng == 0) rng <- 1
  S <- 1 - abs(outer(v, v, "-")) / rng
  diag(S) <- NA
  simTie <- mean(S[which(A == 1)], na.rm = TRUE)
  simAll <- mean(S, na.rm = TRUE)
  list(r = r, I = I, simTie = simTie, simAll = simAll,
       simGap = simTie - simAll, n = n, m = m, dens = m / (n * (n - 1)))
}

perm_null <- function(A, v) {
  ok <- !is.na(v); A <- A[ok, ok, drop = FALSE]; v <- v[ok]
  diag(A) <- NA
  E <- which(A == 1, arr.ind = TRUE)
  replicate(NPERM, { p <- sample(v); suppressWarnings(cor(p[E[, 1]], p[E[, 2]])) })
}

run <- function(file, tag) {
  g <- readRDS(file)
  out <- list()
  for (cl in names(g)) {
    z <- g[[cl]]
    N <- as.matrix(z$depvars$net[, , 1])
    A <- ifelse(N %in% c(1, 11), 1, ifelse(N %in% c(0, 10), 0, NA))
    dim(A) <- dim(N)
    for (t in tr) {
      v <- as.numeric(z$depvars[[t]][, , 1])
      mm <- metrics(A, v); if (is.null(mm)) next
      nul <- perm_null(A, v)
      z_sc <- (mm$r - mean(nul, na.rm = TRUE)) / sd(nul, na.rm = TRUE)
      p <- (1 + sum(abs(nul - mean(nul, na.rm = TRUE)) >=
                    abs(mm$r - mean(nul, na.rm = TRUE)), na.rm = TRUE)) / (NPERM + 1)
      out[[length(out) + 1]] <- data.frame(
        dataset = tag, class = cl, trait = t, n = mm$n, ties = mm$m,
        density = round(mm$dens, 4),
        r_assort = round(mm$r, 4), moran_I = round(mm$I, 4),
        sim_tie = round(mm$simTie, 4), sim_all = round(mm$simAll, 4),
        sim_gap = round(mm$simGap, 4),
        null_mean = round(mean(nul, na.rm = TRUE), 4),
        null_sd = round(sd(nul, na.rm = TRUE), 4),
        r_excess = round(mm$r - mean(nul, na.rm = TRUE), 4),
        z_perm = round(z_sc, 3), p_perm = round(p, 4),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

G <- run("gacr_bayes_group_fixed.rds", "GACR")
Z <- run("oz_bayes_group_fixed.rds",   "OZ")
H <- rbind(G, Z)
H$trait <- factor(H$trait, levels = tr)
write.csv(H, "homophily_classroom.csv", row.names = FALSE)
cat("classroom-level rows:", nrow(H), "\n")

## ---- dataset-level pooling ----------------------------------------------
boot_ci <- function(r, w, B = 4000) {
  k <- length(r); if (k < 3) return(c(NA, NA))
  bs <- replicate(B, { i <- sample(k, k, TRUE)
                       sum(r[i] * w[i]) / sum(w[i]) })
  unname(quantile(bs, c(.025, .975)))
}
pool <- do.call(rbind, lapply(split(H, list(H$dataset, H$trait), drop = TRUE), function(d) {
  d <- d[!is.na(d$r_excess), ]
  w <- d$ties
  ## Pool the permutation-CENTRED assortativity. The raw edgewise correlation
  ## carries a structural negative bias in small directed networks (the null
  ## mean is about -1/(n-1), not 0), so r - E[r_null] is the unbiased quantity.
  rw <- sum(d$r_excess * w) / sum(w)
  ci <- boot_ci(d$r_excess, w)
  stouffer <- sum(d$z_perm, na.rm = TRUE) / sqrt(sum(!is.na(d$z_perm)))
  data.frame(dataset = d$dataset[1], trait = d$trait[1], classes = nrow(d),
             actors = sum(d$n), ties = sum(d$ties),
             r_excess_pooled = round(rw, 4),
             CI_low = round(ci[1], 4), CI_high = round(ci[2], 4),
             r_raw_pooled = round(sum(d$r_assort * w) / sum(w), 4),
             null_mean_pooled = round(sum(d$null_mean * w) / sum(w), 4),
             r_sd_between = round(sd(d$r_excess), 4),
             moran_pooled = round(sum(d$moran_I * w) / sum(w), 4),
             simgap_pooled = round(sum(d$sim_gap * w) / sum(w), 4),
             z_stouffer = round(stouffer, 3),
             p_stouffer = signif(2 * pnorm(-abs(stouffer)), 3),
             classes_p05 = sum(d$p_perm < .05),
             stringsAsFactors = FALSE)
}))
pool <- pool[order(pool$dataset, pool$trait), ]
write.csv(pool, "homophily_dataset.csv", row.names = FALSE)
cat("\n=== dataset-level wave-1 homophily (permutation-centred assortativity) ===\n")
print(pool[, c("dataset","trait","classes","ties","r_excess_pooled","CI_low","CI_high",
               "r_raw_pooled","null_mean_pooled","r_sd_between","simgap_pooled",
               "z_stouffer","p_stouffer","classes_p05")],
      row.names = FALSE)
