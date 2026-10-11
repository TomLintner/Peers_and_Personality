## 50_heterogeneity_27classes.R
## source: RSiena_build/meta27/meta27_final.R
source("00_setup.R")

## Two-stage meta-analysis across the 27 GACR classrooms: final report.
## mu.ml       = pooled estimate across classrooms
## sigma.ml    = BETWEEN-CLASSROOM standard deviation of the effect
## Qstat       = Cochran heterogeneity statistic, df = n1 - 1
suppressMessages(library(RSiena))
SC <- WORK
r <- readRDS(file.path(OBJ, "gacr_meta27.rds"))

tab <- function(trait) {
  m <- r[[trait]]$meta
  nm <- names(m)[!grepl("Rate", names(m))]
  do.call(rbind, lapply(nm, function(k) {
    x <- m[[k]]
    if (!is.list(x) || is.null(x[["mu.ml"]])) return(NULL)
    ph <- pchisq(x$Qstat, x$n1 - 1, lower.tail = FALSE)
    data.frame(trait = trait,
               effect = sub("^eval : ", "", k),
               classes = x$n1,
               pooled = round(x$mu.ml, 3),
               se = round(x$mu.ml.se, 3),
               z = round(x$mu.ml / x$mu.ml.se, 2),
               sig = ifelse(abs(x$mu.ml / x$mu.ml.se) > 1.96, "*", ""),
               sd_between = round(x$sigma.ml, 3),
               sd_CI = paste0("[", round(x$sigma.confint[1], 2), ",",
                                   round(x$sigma.confint[2], 2), "]"),
               Q = round(x$Qstat, 1),
               p_het = round(ph, 4),
               het = ifelse(ph < .05, "VARIES", ""))
  }))
}
all <- do.call(rbind, lapply(names(r), tab))

cat("================ POOLED ESTIMATES AND BETWEEN-CLASSROOM VARIATION ================\n")
cat("(27 GACR classrooms fitted separately, combined with siena08)\n\n")
cat("--- structural effects, as estimated alongside each trait ---\n")
st <- all[all$effect %in% c("outdegree (density)", "reciprocity",
                            "GWESP I -> K -> J (69)", "same sex"), ]
print(st[, c("trait", "effect", "classes", "pooled", "se", "z", "sig",
             "sd_between", "sd_CI", "Q", "p_het", "het")], row.names = FALSE)

cat("\n--- personality effects ---\n")
pe <- all[grepl("^trait", all$effect), ]
print(pe[, c("trait", "effect", "classes", "pooled", "se", "z", "sig",
             "sd_between", "sd_CI", "Q", "p_het", "het")], row.names = FALSE)

cat("\n================ SUMMARY ================\n")
cat("effects showing significant between-classroom variation (p < .05):",
    sum(all$p_het < .05), "of", nrow(all), "\n")
v <- all[all$p_het < .05, c("trait", "effect", "sd_between", "p_het")]
if (nrow(v)) print(v, row.names = FALSE)
cat("\npersonality effects with a significant POOLED estimate:\n")
s <- pe[pe$sig == "*", c("trait", "effect", "pooled", "se", "z")]
if (nrow(s)) print(s, row.names = FALSE) else cat("  none\n")
write.csv(all, file.path(SC, "gacr_meta27_table.csv"), row.names = FALSE)
cat("\nwritten gacr_meta27_table.csv\n")
