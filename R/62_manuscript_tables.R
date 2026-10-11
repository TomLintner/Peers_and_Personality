## 62_manuscript_tables.R
## source: Results/make_tables.R
source("00_setup.R")

## Build APA-style tables for the paper draft.
setwd(OUT)
T <- "tables"

## ---------- TABLE 1: sample descriptives ----------
g <- readRDS(file.path(OBJ, "gacr_bayes_group_fixed.rds"))
o <- readRDS(file.path(OBJ, "oz_bayes_group_fixed.rds"))
sz <- function(G) sapply(G, function(d) attr(d$depvars[[1]], "netdims")[1])
dj <- function(G) sapply(G, function(d) { A <- d$depvars[[1]]
  a <- A[,,1]; b <- A[,,2]
  a <- ifelse(a %in% c(1,11),1,ifelse(a %in% c(0,10),0,NA))
  b <- ifelse(b %in% c(1,11),1,ifelse(b %in% c(0,10),0,NA))
  ok <- !is.na(a)&!is.na(b)
  n11<-sum(a==1&b==1&ok); n10<-sum(a==1&b==0&ok); n01<-sum(a==0&b==1&ok)
  c(d1=mean(a,na.rm=TRUE), d2=mean(b,na.rm=TRUE),
    jac=if(n11+n10+n01==0) NA else n11/(n11+n10+n01)) })
G1 <- dj(g); G2 <- dj(o)
tab1 <- data.frame(
  Characteristic = c("Classrooms", "Pupils", "Class size, range", "Interval between waves",
                     "Personality instrument", "Sociometric question",
                     "Network density, Wave 1", "Network density, Wave 2",
                     "Jaccard index, mean", "Jaccard index, range"),
  Sample_1_Defending_victims_of_bullying = c(length(g), sum(sz(g)),
    paste(range(sz(g)), collapse = " to "), "6 months",
    "BFI-44 (44 items)", "Liked most",
    sprintf("%.3f", mean(G1["d1",])), sprintf("%.3f", mean(G1["d2",])),
    sprintf("%.2f", mean(G1["jac",], na.rm=TRUE)),
    sprintf("%.2f to %.2f", min(G1["jac",],na.rm=TRUE), max(G1["jac",],na.rm=TRUE))),
  Sample_2_Moment_of_change = c(length(o), sum(sz(o)),
    paste(range(sz(o)), collapse = " to "), "8 weeks",
    "BFI-2 short form (11 items)", "Friendship",
    sprintf("%.3f", mean(G2["d1",])), sprintf("%.3f", mean(G2["d2",])),
    sprintf("%.2f", mean(G2["jac",], na.rm=TRUE)),
    sprintf("%.2f to %.2f", min(G2["jac",],na.rm=TRUE), max(G2["jac",],na.rm=TRUE))),
  stringsAsFactors = FALSE)
write.csv(tab1, file.path(T,"paper_table1_samples.csv"), row.names=FALSE)
cat("=== TABLE 1 ===\n"); print(tab1, row.names=FALSE)

## ---------- TABLE 2: wave-1 homophily ----------
h <- read.csv(file.path(T,"homophily_dataset.csv"), stringsAsFactors=FALSE)
trn <- c(E="Extraversion", A="Agreeableness", C="Conscientiousness",
         N="Neuroticism", O="Openness")
h$Trait <- trn[h$trait]
h$Sample <- ifelse(h$dataset=="Defending Victims of Bullying", "Sample 1 (Defending Victims of Bullying)", "Sample 2 (Moment of Change)")
tab2 <- data.frame(Sample=h$Sample, Trait=h$Trait, k=h$classes,
  Excess_r = sprintf("%.3f", h$r_excess_pooled),
  CI = sprintf("[%.3f, %.3f]", h$CI_low, h$CI_high),
  SD_between = sprintf("%.3f", h$r_sd_between),
  z = sprintf("%.2f", h$z_stouffer),
  p = ifelse(h$p_stouffer < .001, "< .001", sub("^0","",sprintf("%.3f", h$p_stouffer))),
  stringsAsFactors=FALSE)
tab2 <- tab2[order(tab2$Sample, match(tab2$Trait, trn)), ]
write.csv(tab2, file.path(T,"paper_table2_homophily.csv"), row.names=FALSE)
cat("\n=== TABLE 2 ===\n"); print(tab2, row.names=FALSE)

## ---------- TABLE 3: posterior estimates ----------
cmp <- read.csv(file.path(T,"sample_comparison.csv"), stringsAsFactors=FALSE)
lab <- c("constant net rate (period 1)"="Rate of network change",
  "outdegree (density)"="Outdegree (density)", "reciprocity"="Reciprocity",
  "GWESP I -> K -> J (#)"="Transitive closure (GWESP)",
  "indegree-popularity (sqrt)"="Indegree popularity", "outdegree-activity (sqrt)"="Outdegree activity",
  "sex alter"="Gender alter", "sex ego"="Gender ego", "same sex"="Same gender")
for (t in names(trn)) {
  lab[paste(t,"alter")] <- paste(trn[t], "alter")
  lab[paste(t,"ego")] <- paste(trn[t], "ego")
  lab[paste(t,"similarity")] <- paste(trn[t], "similarity")
  lab[paste("rate",t,"(period 1)")] <- paste("Rate of change,", trn[t])
  lab[paste(t,"linear shape")] <- paste(trn[t], "linear shape")
  lab[paste(t,"quadratic shape")] <- paste(trn[t], "quadratic shape")
  lab[paste(t,"average similarity")] <- paste(trn[t], "average similarity")
  lab[paste(t,"reciprocated degree")] <- paste(trn[t], "reciprocated degree")
}
cmp$Parameter <- ifelse(cmp$parameter %in% names(lab), lab[cmp$parameter], cmp$parameter)
grp <- function(p) {
  if (grepl("Rate", p)) "Rate parameters"
  else if (p %in% c("Outdegree (density)","Reciprocity","Transitive closure (GWESP)",
                    "Indegree popularity","Outdegree activity","Gender alter","Gender ego","Same gender")) "Network structure and gender"
  else if (grepl("alter$|ego$|similarity$", p) & !grepl("average", p)) "Personality selection"
  else if (grepl("average similarity", p)) "Peer influence"
  else if (grepl("reciprocated degree", p)) "Reciprocated degree on trait change"
  else "Trait shape" }
cmp$Block <- sapply(cmp$Parameter, grp)
ord <- c("Network structure and gender","Personality selection","Peer influence",
         "Trait shape","Reciprocated degree on trait change","Rate parameters")
cmp$Block <- factor(cmp$Block, levels=ord)
cmp <- cmp[order(cmp$Block), ]
tab3 <- data.frame(Block=as.character(cmp$Block), Parameter=cmp$Parameter,
  Defending_victims_of_bullying=trimws(cmp$Defending_victims_of_bullying), Defending_victims_of_bullying_SDbtw=ifelse(is.na(cmp$Defending_victims_of_bullying_SDbtw),"",sprintf("%.2f",cmp$Defending_victims_of_bullying_SDbtw)),
  Moment_of_change=trimws(cmp$Moment_of_change), Moment_of_change_SDbtw=ifelse(is.na(cmp$Moment_of_change_SDbtw),"",sprintf("%.2f",cmp$Moment_of_change_SDbtw)),
  stringsAsFactors=FALSE)
write.csv(tab3, file.path(T,"paper_table3_estimates.csv"), row.names=FALSE)
cat("\n=== TABLE 3 (first 20 rows) ===\n"); print(head(tab3, 20), row.names=FALSE)
cat("\nrows:", nrow(tab3), "\n")
cat("\nwritten paper_table1_samples.csv, paper_table2_homophily.csv, paper_table3_estimates.csv\n")
