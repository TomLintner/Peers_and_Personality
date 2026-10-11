## 12_gof_phase3.R
## source: RSiena_build/scripts/gof_phase3.R
source("00_setup.R")

## =====================================================================
## Goodness of fit for the published method-of-moments specification.
## Answers coauthor comment 46 ("sienaGOF, at least for MoM").
##
## The converged fits already exist and are NOT re-estimated:
##   output/gacr_final_ans.rds  205 parameters, 27 groups, tconv.max 0.2456
##   output/oz_final_ans.rds    151 parameters, 18 groups, tconv.max 0.2210
## They were produced by gacr_conv.R / oz_conv.R with RSiena, set.seed(2012),
## thetaBound = 50, and n3 escalated to 40000 to push tconv.max under its
## sqrt(p/n3) floor. Both carry returnDeps = FALSE and a sims list of NULLs,
## so sienaGOF refuses them. All that is needed is one Phase-3-only pass at
## the converged estimates with returnDeps = TRUE: nsub = 0, prevAns = the
## saved fit. No Phase 1, no Phase 2, no re-estimation.
##
## Usage:  Rscript gof_phase3.R s1     (Sample 1, GACR, 27 classes)
##         Rscript gof_phase3.R s2     (Sample 2, OZ,   18 classes)
## =====================================================================
suppressMessages(library(RSiena))
set.seed(2012)

ARG <- commandArgs(trailingOnly = TRUE)
S   <- if (length(ARG)) ARG[1] else "s1"
B   <- WORK
O   <- MOM
G   <- file.path(OUT, "gof")
dir.create(G, showWarnings = FALSE)
setwd(G)

CFG <- list(
  s1 = list(fit = "gacr_final_ans.rds", grp = "gacr_coev_group.rds",
            lab = "Sample 1 (Defending Victims of Bullying, 27 classes)", nodes = 14),
  s2 = list(fit = "oz_final_ans.rds",   grp = "oz_coev_group.rds",
            lab = "Sample 2 (Moment of Change, 18 classes)",            nodes = 14))[[S]]
LOG <- file.path(G, paste0("gof_", S, "_log.txt"))
tr5 <- c("E", "A", "C", "N", "O")

say <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(..., collapse = ""))
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE); flush.console() }
cat("", file = LOG)
say("=== ", CFG$lab, " ===")

## ---------------------------------------------------------------- data and specification
grp <- readRDS(file.path(O, CFG$grp))
say("group object: ", length(grp), " groups")

## rebuilt verbatim from gacr_conv.R / oz_conv.R -- the two scripts use the same specification
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

ans <- readRDS(file.path(O, CFG$fit))
say("loaded ", CFG$fit, ": ", length(ans$theta), " parameters, tconv.max ", round(ans$tconv.max, 4))

## hard stop if the rebuilt specification does not match the one that was fitted
E <- eff; class(E) <- "data.frame"; inc <- E[E$include, ]
if (nrow(inc) != length(ans$theta))
  stop(sprintf("specification mismatch: rebuilt %d parameters, fit has %d",
               nrow(inc), length(ans$theta)))
## Compare the substantive identity of each effect, not its printed name. The display name renders
## the effect's parm at include time, so a freshly built object shows "GWESP I -> K -> J (#)" where
## the stored fit has "GWESP I -> K -> J (69)" -- the same gwespFF effect with the same parm.
key <- function(d) paste(d$name, d$shortName, d$type, d$interaction1, d$interaction2, d$parm, sep = "|")
k_new <- key(inc); k_fit <- key(ans$effects)
if (!identical(k_new, k_fit)) {
  d <- which(k_new != k_fit)
  stop(sprintf("specification mismatch at %d positions, first: rebuilt '%s' vs fit '%s'",
               length(d), k_new[d[1]], k_fit[d[1]])) }
say("specification verified: ", nrow(inc),
    " parameters match the saved fit on dependent variable, effect, type, interaction and parm")
## and report the printed-name differences, which should be cosmetic only
nd <- which(as.character(inc$effectName) != as.character(ans$effects$effectName))
if (length(nd)) for (i in nd)
  say("  display name differs (same effect): '", inc$effectName[i], "' vs '",
      ans$effects$effectName[i], "'")

## ---------------------------------------------------------------- Phase 3 only
f_dep <- file.path(G, paste0("mom_", S, "_deps.rds"))
if (file.exists(f_dep)) {
  ansd <- readRDS(f_dep); say("Phase 3 skipped, ", basename(f_dep), " already exists")
} else {
  alg <- sienaAlgorithmCreate(projname = paste0("gof_", S), seed = 2012, n3 = 1000, nsub = 0)
  say("Phase 3 only at the converged estimates, n3 = 1000, returnDeps = TRUE, ",
      CFG$nodes, " nodes")
  t0 <- proc.time()
  ansd <- siena07(alg, data = grp, effects = eff, prevAns = ans, thetaBound = 50,
                  returnDeps = TRUE, useCluster = TRUE, nbrNodes = CFG$nodes, batch = TRUE)
  say("  done in ", round((proc.time() - t0)[3] / 60, 1), " min | returnDeps = ",
      isTRUE(ansd$returnDeps), " | sims = ", sum(!sapply(ansd$sims, is.null)))
  saveRDS(ansd, f_dep)
  say("saved ", basename(f_dep), " (",
      round(as.numeric(object.size(ansd)) / 1024^2), " MB in memory)")
}
rm(ans); gc()

## ---------------------------------------------------------------- sienaGOF
gnames <- names(grp)
net_stats <- list(
  "outdegree distribution" = function(...) OutdegreeDistribution(..., levls = 0:8, cumulative = FALSE),
  "indegree distribution"  = function(...) IndegreeDistribution(...,  levls = 0:8, cumulative = FALSE),
  "triad census"           = function(...) TriadCensus(...))

res <- list(); rows <- list()
f_gof <- file.path(G, paste0("gof_", S, ".rds"))
say("sienaGOF over ", length(gnames), " groups x (3 network + 5 behaviour) statistics")

## A sienaGOF object is a list of sienaGofTest elements; with join = TRUE there is exactly one.
## Its fields are p, ObservedTestStat (the joint Mahalanobis distance), SimulatedTestStat,
## Simulations, Observations, InvCovSimStats and Rank. There is no "joinedTests" attribute.
grab <- function(o) {
  t1 <- o[[1]]
  c(p = as.numeric(t1$p), MHD = as.numeric(t1$ObservedTestStat), rank = as.numeric(t1$Rank)) }

for (gi in seq_along(gnames)) {
  g <- gnames[gi]
  for (nm in names(net_stats)) {
    o <- tryCatch(sienaGOF(ansd, net_stats[[nm]], varName = "net", groupName = g,
                           join = TRUE, verbose = FALSE),
                  error = function(e) conditionMessage(e))
    if (inherits(o, "sienaGOF")) {
      v <- grab(o); res[[paste(g, nm)]] <- o
      rows[[paste(g, nm)]] <- data.frame(group = g, variable = "net", statistic = nm,
                                         p = round(unname(v["p"]), 4), MHD = round(unname(v["MHD"]), 2), rank = unname(v["rank"]),
                                         row.names = NULL)
    } else {
      rows[[paste(g, nm)]] <- data.frame(group = g, variable = "net", statistic = nm,
                                         p = NA_real_, MHD = NA_real_, rank = NA_real_, row.names = NULL)
      say("  ! ", g, " ", nm, ": ", o) } }
  for (bv in tr5) {
    o <- tryCatch(sienaGOF(ansd, BehaviorDistribution, varName = bv, groupName = g,
                           join = TRUE, levls = 1:5, verbose = FALSE),
                  error = function(e) conditionMessage(e))
    if (inherits(o, "sienaGOF")) {
      v <- grab(o); res[[paste(g, bv)]] <- o
      rows[[paste(g, bv)]] <- data.frame(group = g, variable = bv, statistic = "behaviour distribution",
                                         p = round(unname(v["p"]), 4), MHD = round(unname(v["MHD"]), 2), rank = unname(v["rank"]),
                                         row.names = NULL)
    } else {
      rows[[paste(g, bv)]] <- data.frame(group = g, variable = bv, statistic = "behaviour distribution",
                                         p = NA_real_, MHD = NA_real_, rank = NA_real_, row.names = NULL)
      say("  ! ", g, " ", bv, ": ", o) } }
  if (gi %% 3 == 0 || gi == length(gnames)) {
    say("  ", gi, "/", length(gnames), " groups")
    saveRDS(res, f_gof)
    write.csv(do.call(rbind, rows), paste0("gof_", S, "_table.csv"), row.names = FALSE) } }

saveRDS(res, f_gof)
tab <- do.call(rbind, rows)
write.csv(tab, paste0("gof_", S, "_table.csv"), row.names = FALSE)

## ---------------------------------------------------------------- summary
say("=== fit by statistic, across groups ===")
for (st in unique(tab$statistic)) {
  z <- tab$p[tab$statistic == st & !is.na(tab$p)]
  say(sprintf("  %-24s n=%2d  median p = %.3f  p<.05 in %2d  p<.01 in %2d",
              st, length(z), median(z), sum(z < .05), sum(z < .01))) }
say("=== behaviour fit by trait ===")
for (bv in tr5) {
  z <- tab$p[tab$variable == bv & !is.na(tab$p)]
  if (length(z)) say(sprintf("  %-2s  n=%2d  median p = %.3f  p<.05 in %2d",
                             bv, length(z), median(z), sum(z < .05))) }
say("wrote gof_", S, "_table.csv and gof_", S, ".rds in ", G)
say("DONE")
