## 13_gof_summarise.R
## source: RSiena_build/scripts/gof_summarize.R
source("00_setup.R")

## Pulls the two GOF tables together, and checks one thing that bears on how the triad-census
## misfit should be read: sienaGOF tests the CONSTANT-parameter multigroup fit, whereas the paper
## reports a multilevel model in which the structural parameters vary across classrooms. If the
## between-classroom SDs on density, reciprocity and GWESP are large, part of the misfit is the
## homogeneity assumption rather than a missing effect.
B <- RAW
G <- file.path(B, "RSiena_build", "gof")
S <- file.path(B, "EJP_submission", "supplemental")

lab <- c(s1 = "Sample 1", s2 = "Sample 2")
all <- do.call(rbind, lapply(c("s1", "s2"), function(s) {
  t <- read.csv(file.path(G, sprintf("gof_%s_table.csv", s)), stringsAsFactors = FALSE)
  t$sample <- lab[s]; t }))

cat("=== uniformity check: under a correct specification p should be uniform on [0,1] ===\n")
cat("   so expect median p near .50 and about 5 per cent below .05\n\n")
agg <- function(d) data.frame(
  n = nrow(d), median_p = round(median(d$p, na.rm = TRUE), 3),
  prop_lt_05 = round(mean(d$p < .05, na.rm = TRUE), 3),
  prop_lt_01 = round(mean(d$p < .01, na.rm = TRUE), 3),
  median_MHD = round(median(d$MHD, na.rm = TRUE), 1))
out <- do.call(rbind, lapply(split(all, list(all$sample, all$statistic), drop = TRUE), function(d)
  cbind(sample = d$sample[1], statistic = d$statistic[1], agg(d))))
print(out[order(out$statistic, out$sample), ], row.names = FALSE)

cat("\n=== triad census, classroom by classroom ===\n")
for (s in unique(all$sample)) {
  x <- all[all$sample == s & all$statistic == "triad census", ]
  x <- x[order(x$p), ]
  cat("\n--", s, "-- worst-fitting first\n")
  print(data.frame(class = x$group, p = x$p, MHD = x$MHD, rank = x$rank), row.names = FALSE) }

cat("\n=== between-classroom SDs of the structural parameters in the reported Bayesian model ===\n")
for (f in c("bayes_defending_victims.csv", "bayes_moment_of_change.csv")) {
  p <- read.csv(file.path(S, f), stringsAsFactors = FALSE)
  k <- p[p$parameter %in% c("outdegree (density)", "reciprocity", "GWESP I -> K -> J (#)",
                            "indegree-popularity (sqrt)", "outdegree-activity (sqrt)"), ]
  cat("\n--", sub("bayes_", "", sub(".csv", "", f)), "\n")
  print(data.frame(parameter = k$parameter, M = k$M, SD_within = k$SD,
                   SD_between = k$SD_between), row.names = FALSE) }

write.csv(all, file.path(S, "gof_classroom_table.csv"), row.names = FALSE)
write.csv(out[order(out$statistic, out$sample), ], file.path(S, "gof_summary.csv"), row.names = FALSE)
cat("\nwrote supplemental/gof_classroom_table.csv and supplemental/gof_summary.csv\n")
