## 60_posterior_table_sample1.R
## source: Results/scripts/full_table.R
source("00_setup.R")

## Full posterior table: mean, SD, 95% credible interval, probability of
## direction (PD), split-Rhat and ESS, for every parameter estimated so far.
## Chain = stages 1..k stacked (stage 1's 100 warming draws discarded).
## For the random effects, also the posterior mean between-classroom SD,
## i.e. sqrt of the diagonal of Sigma.
suppressMessages(library(RSienaTest))
SC <- WORK
s <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds"))
E <- s$eff; class(E) <- "data.frame"
selR  <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
labMu <- paste(E$name[selR], E$effectName[selR])
fxi   <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labEta<- paste(E$name[fxi], E$effectName[fxi])

MU <- NULL; ET <- NULL; SGl <- list(); used <- c()
for (k in 1:9) {
  f <- file.path(CHN, sprintf("bayes12h_stage%d.rds", k)); if (!file.exists(f)) next
  g <- readRDS(f)
  ok <- !is.na(g$ThinPosteriorMu[, 1])
  if (k == 1) ok[1:100] <- FALSE
  MU <- rbind(MU, g$ThinPosteriorMu[ok, , drop = FALSE])
  ET <- rbind(ET, g$ThinPosteriorEta[ok, , drop = FALSE])
  sg <- g$ThinPosteriorSigma
  if (!is.null(sg)) SGl[[length(SGl) + 1]] <-
      t(apply(sg[ok, , , drop = FALSE], 1, function(m) sqrt(diag(m))))
  used <- c(used, sprintf("s%d:%d", k, sum(ok)))
}
SD_between <- if (length(SGl)) colMeans(do.call(rbind, SGl)) else rep(NA, ncol(MU))

splitRhat <- function(x) { n <- length(x); h <- floor(n / 2)
  if (h < 10) return(NA); a <- x[1:h]; b <- x[(n - h + 1):n]
  W <- mean(c(var(a), var(b))); if (W <= 0) return(NA)
  gm <- mean(c(mean(a), mean(b))); B <- h * sum((c(mean(a), mean(b)) - gm)^2)
  sqrt((((h - 1) / h) * W + B / h) / W) }
rk  <- function(x) qnorm((rank(x) - 3/8) / (length(x) + 1/4))
ess <- function(x) { n <- length(x)
  r <- acf(x, lag.max = min(500, n - 2), plot = FALSE)$acf[-1]
  k <- which(r < 0.05)[1]; if (is.na(k)) k <- length(r)
  max(1, n / (1 + 2 * sum(r[1:max(1, k - 1)]))) }

mk <- function(M, labs, blk, sdb = NULL) {
  d <- data.frame(block = blk, parameter = labs,
    M    = round(colMeans(M), 3),
    SD   = round(apply(M, 2, sd), 3),
    lo95 = round(apply(M, 2, quantile, .025), 3),
    hi95 = round(apply(M, 2, quantile, .975), 3),
    PD   = round(pmax(colMeans(M > 0), colMeans(M < 0)), 3),
    Rhat = round(apply(M, 2, function(c) splitRhat(rk(c))), 3),
    ESS  = round(apply(M, 2, ess)))
  if (!is.null(sdb)) d <- cbind(d[, 1:6], SD_btw = round(sdb, 3), d[, 7:9])
  d }

tabMu  <- mk(MU, labMu,  "mu",  SD_between)
tabEta <- mk(ET, labEta, "eta")
cat("chain: stacked", paste(used, collapse = " "), "->", nrow(MU), "post-warm draws\n")
cat("M/SD = posterior mean and SD of the parameter; SD_btw = posterior mean\n",
    "between-classroom SD (sqrt diag Sigma); PD = posterior prob. of direction.\n")
cat("\n============ RANDOM EFFECTS: population means mu ============\n")
print(tabMu[, -1], row.names = FALSE)
cat("\n============ FIXED EFFECTS eta (constant across classrooms) ============\n")
print(tabEta[, -1], row.names = FALSE)
write.csv(rbind(tabMu[, c("block","parameter","M","SD","lo95","hi95","PD","Rhat","ESS")],
                tabEta),
          file.path(SC, "bayes_full_table.csv"), row.names = FALSE)
write.csv(tabMu, file.path(SC, "bayes_full_table_mu.csv"), row.names = FALSE)
cat("\nwritten bayes_full_table.csv / _mu.csv\n")
