## 10_mom_sample1.R
## source: RSiena_build/scripts/gacr_conv.R
source("00_setup.R")

## =====================================================================
## GACR multigroup co-evolution model: 1 network + 5 Big Five behaviours
## Reproducible: set.seed(2012) in R and seed = 2012 in every algorithm object.
## Convergence strategy: tconv.max has a floor of about sqrt(p / n3) at the
## solution. With p = 205, n3 = 3000 gives ~0.26 -- below 0.25 is unreachable.
## So: locate the solution at n3 = 3000, then raise n3 to shrink the floor.
## =====================================================================
suppressMessages(library(RSiena))
set.seed(2012)

SC <- WORK
setwd(SC)
G <- readRDS(file.path(OBJ, "gacr_siena.rds"))
tr5 <- c("E", "A", "C", "N", "O")

objs <- lapply(G$data, function(z) {
  args <- list(net = z$net)
  for (t in tr5) args[[t]] <- sienaDependent(z$beh[[t]], type = "behavior", allowOnly = FALSE)
  args$sex <- coCovar(as.numeric(z$sex))
  do.call(sienaDataCreate, args)
})
grp <- sienaGroupCreate(objs)
saveRDS(grp, file.path(OBJ, "gacr_coev_group.rds"))

eff <- getEffects(grp)
eff <- includeEffects(eff, gwespFF, inPopSqrt, outActSqrt)
eff <- includeEffects(eff, altX,  interaction1 = "sex")
eff <- includeEffects(eff, egoX,  interaction1 = "sex")
eff <- includeEffects(eff, sameX, interaction1 = "sex")
for (t in tr5) {
  eff <- includeEffects(eff, egoX, interaction1 = t)
  eff <- includeEffects(eff, altX, interaction1 = t)
  eff <- includeEffects(eff, simX, interaction1 = t)
}
for (t in tr5) {
  eff <- includeEffects(eff, avSim,    name = t, interaction1 = "net")
  eff <- includeEffects(eff, recipDeg, name = t, interaction1 = "net")
}
E <- eff; class(E) <- "data.frame"; inc <- E[E$include, ]
p <- nrow(inc)
cat("groups:", length(objs), "| actors:", sum(sapply(G$data, function(z) z$n)), "\n")
cat("parameters:", p, "(", sum(inc$type == "eval"), "eval +", sum(inc$type == "rate"), "rate )\n")
cat("predicted convergence floor sqrt(p/n3):",
    paste(sprintf("n3=%d -> %.3f", c(3000, 10000, 20000, 40000),
                  sqrt(p / c(3000, 10000, 20000, 40000))), collapse = " | "), "\n\n")
flush.console()

step <- function(prev, n3, nsub, tag) {
  set.seed(2012)
  alg <- sienaAlgorithmCreate(projname = paste0("gacr_conv_", tag), seed = 2012,
                              n3 = n3, nsub = nsub)
  t0 <- proc.time()
  a <- try(siena07(alg, data = grp, effects = eff, prevAns = prev, thetaBound = 50,
                   useCluster = TRUE, nbrNodes = 8, batch = TRUE), silent = TRUE)
  if (inherits(a, "try-error")) { cat("[", tag, "] FAILED:\n"); print(a); return(prev) }
  cat(sprintf("[%s] n3=%-6d nsub=%d -> tconv.max = %.4f | max|t_conv| = %.4f | above 0.10: %d/%d | %.1f min\n",
              tag, n3, nsub, a$tconv.max, max(abs(a$tconv)), sum(abs(a$tconv) > 0.10),
              length(a$tconv), (proc.time() - t0)[3] / 60))
  flush.console()
  saveRDS(a, file.path(SC, paste0("gacr_conv_", tag, ".rds")))
  a
}

T0 <- proc.time()
a <- step(NULL, 3000,  4, "s1")     # locate the solution
if (a$tconv.max > 0.25) a <- step(a, 10000, 2, "s2")
if (a$tconv.max > 0.25) a <- step(a, 20000, 2, "s3")
if (a$tconv.max > 0.25) a <- step(a, 40000, 2, "s4")
cat("\nTOTAL:", round((proc.time() - T0)[3] / 60, 1), "min\n")
saveRDS(a, file.path(MOM, "gacr_final_ans.rds"))

cat("\n=== FINAL CONVERGENCE ===\n")
cat("overall max convergence ratio:", round(a$tconv.max, 4),
    if (a$tconv.max < 0.25) "  <-- CONVERGED (< 0.25)" else "  <-- still above 0.25", "\n")
cat("max individual |t_conv|:", round(max(abs(a$tconv)), 4),
    "| above 0.10:", sum(abs(a$tconv) > 0.10), "of", length(a$tconv), "\n")

est <- a$theta; se <- sqrt(diag(a$covtheta)); tv <- est / se; pv <- 2 * pnorm(-abs(tv))
pt <- a$effects$name; nm <- a$effects$effectName; ty <- a$effects$type
tab <- data.frame(part = pt, type = ty, effect = nm, est = round(est, 3), se = round(se, 3),
                  t = round(tv, 2), p = format.pval(pv, digits = 2, eps = 1e-4),
                  sig = ifelse(pv < .001, "***", ifelse(pv < .01, "**", ifelse(pv < .05, "*", ""))))
cat("\n================ NETWORK EVALUATION FUNCTION ================\n")
print(tab[pt == "net" & ty == "eval", c("effect", "est", "se", "t", "p", "sig")], row.names = FALSE)
for (t in tr5) { cat("\n--", t, "--\n")
  print(tab[pt == t & ty == "eval", c("effect", "est", "se", "t", "p", "sig")], row.names = FALSE) }
cat("\n=== rate parameters ===\n")
for (q in c("net", tr5)) { r <- est[pt == q & ty == "rate"]
  cat(sprintf("  %-4s n=%2d median %6.2f  range %5.2f - %6.2f\n", q, length(r), median(r), min(r), max(r))) }

cat("\n=== joint tests ===\n")
wald <- function(ii) { W <- drop(t(est[ii]) %*% solve(a$covtheta[ii, ii]) %*% est[ii])
  sprintf("chi2 = %6.2f  df = %2d  p = %s", W, length(ii),
          format.pval(pchisq(W, length(ii), lower.tail = FALSE), digits = 3, eps = 1e-4)) }
for (b in c("similarity", "alter", "ego")) {
  ii <- which(nm %in% paste(tr5, b) & pt == "net")
  cat(sprintf("  network %-11s: %s  est: %s\n", b, wald(ii),
              paste(sprintf("%+.3f", est[ii]), collapse = " "))) }
for (b in c("average similarity", "reciprocated degree")) {
  ii <- which(nm %in% paste(tr5, b))
  cat(sprintf("  behaviour %-20s: %s  est: %s\n", b, wald(ii),
              paste(sprintf("%+.3f", est[ii]), collapse = " "))) }
