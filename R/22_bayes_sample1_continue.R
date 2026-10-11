## 22_bayes_sample1_continue.R
## source: Results/scripts/bayes_cont.R
source("00_setup.R")

## ==========================================================================
## Continue the GACR sienaBayes chain from stage 1.
## BUDGET IS NOW IN ITERATIONS, not wall clock -- immune to the machine
## sleeping, which is what truncated the previous run (13.1 h clock but only
## ~4.5 h of compute).
## 1200 further main iterations at the measured ~40 s/iteration ~ 13 h COMPUTE,
## in 4 chunks of 300. Each chunk: save, glue, report cumulative Rhat + ESS.
## Model and priors unchanged: 43 evaluation effects, p1 = 19, p2 = 30,
## priorSigEta = 1, priorKappa = 1.0, nSampVarying = 3, nSampConst = 10.
## Behaviour side is linear + quad + avSim + recipDeg per trait; no effFrom.
## ==========================================================================
suppressMessages(library(RSienaTest))
set.seed(2012)
SC <- WORK
WD <- file.path(SC, "bayes12h"); dir.create(WD, showWarnings = FALSE); setwd(WD)
s <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds"))
E <- s$eff; class(E) <- "data.frame"
p2 <- sum(E$include & !E$basicRate & !E$randomEffects & E$type == "eval")
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
labMu <- paste(E$name[selR], E$shortName[selR])
fxi <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labEta <- paste(E$name[fxi], E$shortName[fxi])

CHUNK <- 300; NSTAGE <- 4          # 1200 iterations ~ 13 h of compute

splitRhat <- function(x) { n <- length(x); h <- floor(n / 2)
  if (h < 10) return(NA)
  a <- x[1:h]; b <- x[(n - h + 1):n]
  W <- mean(c(var(a), var(b))); if (W <= 0) return(NA)
  gm <- mean(c(mean(a), mean(b)))
  B <- h * sum((c(mean(a), mean(b)) - gm)^2)
  sqrt((((h - 1) / h) * W + B / h) / W) }
rk  <- function(x) qnorm((rank(x) - 3/8) / (length(x) + 1/4))
ess <- function(x) { n <- length(x)
  r <- acf(x, lag.max = min(400, n - 2), plot = FALSE)$acf[-1]
  k <- which(r < 0.05)[1]; if (is.na(k)) k <- length(r)
  max(1, n / (1 + 2 * sum(r[1:max(1, k - 1)]))) }
diag1 <- function(M, tag, labs) {
  R <- apply(M, 2, function(c) splitRhat(rk(c))); S <- apply(M, 2, ess)
  cat(sprintf("  %-4s n=%4d | Rhat med %.3f max %.3f | <1.01: %2d/%2d | ESS med %4.0f min %4.0f\n",
      tag, nrow(M), median(R, na.rm = TRUE), max(R, na.rm = TRUE),
      sum(R < 1.01, na.rm = TRUE), length(R), median(S), min(S)))
  w <- order(-R)[1:3]
  cat("       worst:", paste(sprintf("%s %.3f", labs[w], R[w]), collapse = " | "), "\n")
  data.frame(par = labs, Rhat = round(R, 4), ESS = round(S)) }
report <- function(g, stage, nw) {
  mu <- g$ThinPosteriorMu; mu <- mu[!is.na(mu[, 1]), , drop = FALSE]
  et <- g$ThinPosteriorEta; et <- et[!is.na(et[, 1]), , drop = FALSE]
  if (nw > 0 && nrow(mu) > nw) { mu <- mu[-(1:nw), , drop = FALSE]; et <- et[-(1:nw), , drop = FALSE] }
  cat("\n######## CUMULATIVE after stage", stage, "########\n")
  a <- diag1(mu, "mu", labMu); b <- diag1(et, "eta", labEta)
  write.csv(rbind(cbind(block = "mu", a), cbind(block = "eta", b)),
            file.path(SC, sprintf("bayes12h_diag_stage%d.csv", stage)), row.names = FALSE)
  flush.console() }

acc  <- readRDS(file.path(SC, "bayes12h_stage1.rds"))
prev <- acc
cat("resuming from stage 1:", sum(!is.na(acc$ThinPosteriorMu[, 1])), "draws\n")
cat("budget:", NSTAGE * CHUNK, "further iterations (~",
    round(NSTAGE * CHUNK * 40 / 3600, 1), "h of compute )\n")
cat("start:", format(Sys.time()), "\n\n"); flush.console()

totIter <- 0; totComp <- 0
for (k in 2:(NSTAGE + 1)) {
  cat("===== stage", k, "=====", format(Sys.time()), "\n"); flush.console()
  t0 <- proc.time()[3]
  algo <- sienaAlgorithmCreate(projname = paste0("b12h_s", k), seed = 2012)
  g <- try(sienaBayes(algo, data = s$grp, effects = s$eff,
                      priorMu = s$priorMu, priorSigma = s$priorSig, priorKappa = 1.0,
                      priorSigEta = rep(1, p2),
                      nmain = CHUNK, nrunMHBatches = 20,
                      nSampVarying = 3, nSampConst = 10, saveFreq = 20,
                      prevBayes = prev, nbrNodes = 8, clusterType = "PSOCK",
                      silentstart = TRUE), silent = TRUE)
  el <- proc.time()[3] - t0
  if (inherits(g, "try-error")) { cat("stage", k, "FAILED:\n"); print(g); break }
  saveRDS(g, file.path(CHN, sprintf("bayes12h_stage%d.rds", k)))
  totIter <- totIter + CHUNK
  cat(sprintf("  stage %d done: %d iterations, %.2f h wall, %.1f s/iteration\n",
              k, CHUNK, el / 3600, el / CHUNK)); flush.console()
  gl <- try(glueBayes(acc, g, nwarm2 = 0), silent = TRUE)
  if (inherits(gl, "try-error")) cat("  glueBayes failed; keeping previous accumulation\n")
  else { acc <- gl; saveRDS(acc, file.path(SC, "bayes12h_accumulated.rds")) }
  report(acc, k, 100)
  cat(sprintf("  cumulative NEW iterations: %d of %d\n\n", totIter, NSTAGE * CHUNK))
  prev <- g
}
cat("\nFINISHED", format(Sys.time()), "| new iterations:", totIter, "\n")
