## 51_meta_two_study.R
## source: RSiena_build/scripts/meta08.R
source("00_setup.R")

suppressMessages(library(RSiena))
set.seed(2012)
SC <- WORK
setwd(SC)
g <- readRDS(file.path(MOM, "gacr_final_ans.rds"))
z <- readRDS(file.path(MOM, "oz_final_ans.rds"))
cat("GACR: tconv.max", round(g$tconv.max, 4), "| params", length(g$theta), "\n")
cat("OZ  : tconv.max", round(z$tconv.max, 4), "| params", length(z$theta), "\n\n")

meta <- try(siena08(g, z, projname = "meta_gacr_oz", bound = 10, alpha = 0.05), silent = TRUE)
if (inherits(meta, "try-error")) { cat("siena08 FAILED:\n"); print(meta) } else {
  saveRDS(meta, file.path(SC, "meta_gacr_oz.rds"))
  cat("siena08 completed. components:", length(meta), "\n")
  nm <- names(meta)
  cat("effects analysed:", sum(!grepl("^Rate|rate", nm)), "of", length(nm), "\n")
}

## --- manual two-study random-effects comparison, effect by effect ---
## (with 2 studies siena08's own tests are thin; this makes the contrast explicit)
gi <- which(g$effects$type == "eval"); zi <- which(z$effects$type == "eval")
gk <- paste(g$effects$name[gi], g$effects$effectName[gi])
zk <- paste(z$effects$name[zi], z$effects$effectName[zi])
common <- intersect(gk, zk)
cat("\ncommon evaluation effects:", length(common), "\n")
a <- gi[match(common, gk)]; b <- zi[match(common, zk)]
eg <- g$theta[a]; sg <- sqrt(diag(g$covtheta))[a]
ez <- z$theta[b]; sz <- sqrt(diag(z$covtheta))[b]
# fixed-effect pooled estimate + Cochran Q heterogeneity test (1 df)
w1 <- 1 / sg^2; w2 <- 1 / sz^2
pool <- (w1 * eg + w2 * ez) / (w1 + w2); sep <- sqrt(1 / (w1 + w2))
Q <- w1 * (eg - pool)^2 + w2 * (ez - pool)^2
out <- data.frame(effect = common,
                  GACR = round(eg, 3), se_G = round(sg, 3),
                  OZ = round(ez, 3), se_O = round(sz, 3),
                  pooled = round(pool, 3), se_p = round(sep, 3),
                  z_pool = round(pool / sep, 2),
                  p_pool = round(2 * pnorm(-abs(pool / sep)), 4),
                  Q = round(Q, 2), p_Q = round(pchisq(Q, 1, lower.tail = FALSE), 4))
out$verdict <- ifelse(out$p_Q < .05, "HETEROGENEOUS",
               ifelse(out$p_pool < .05, "consistent+sig", "consistent+ns"))
cat("\n================ TWO-STUDY META-ANALYSIS ================\n")
cat("(pooled = inverse-variance fixed effect; Q = Cochran heterogeneity, 1 df)\n\n")
print(out, row.names = FALSE)
write.csv(out, file.path(SC, "meta_two_study.csv"), row.names = FALSE)

cat("\n=== summary ===\n")
cat("effects with significant heterogeneity (p_Q < .05):", sum(out$p_Q < .05), "of", nrow(out), "\n")
cat("  ", paste(out$effect[out$p_Q < .05], collapse = " | "), "\n")
cat("effects consistent across studies AND significant when pooled:",
    sum(out$verdict == "consistent+sig"), "\n")
cat("  ", paste(out$effect[out$verdict == "consistent+sig"], collapse = " | "), "\n")
