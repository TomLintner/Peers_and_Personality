## 70_supplemental_descriptives.R
## source: EJP_submission/trait_descriptives.R
source("00_setup.R")

## Supplemental Table S1: descriptive statistics, retest reliability and zero-order correlations of the
## five ordinal trait variables that entered the models, per sample. Required by the EJP Open Science
## Policy (effects statement) and the Evaluation Criteria (reporting standards).
suppressMessages({library(RSienaTest); library(officer); library(flextable)})
R <- OUT
setwd(OUT)
TR <- c(E = "Extraversion", A = "Agreeableness", C = "Conscientiousness", N = "Neuroticism", O = "Openness to experience")
grab <- function(f) { G <- readRDS(file.path(R, "objects", f)); tr <- c("E","A","C","N","O")
  out <- lapply(tr, function(b) do.call(rbind, lapply(G$data, function(g) g$beh[[b]]))); names(out) <- tr; out }
summ <- function(L, tag) {
  k <- names(L); W1 <- sapply(L, function(m) m[, 1]); W2 <- sapply(L, function(m) m[, 2])
  d <- data.frame(Trait = unname(TR[k]),
    M1 = round(colMeans(W1, na.rm = TRUE), 2), SD1 = round(apply(W1, 2, sd, na.rm = TRUE), 2),
    M2 = round(colMeans(W2, na.rm = TRUE), 2), SD2 = round(apply(W2, 2, sd, na.rm = TRUE), 2),
    rtt = round(sapply(k, function(b) cor(L[[b]][, 1], L[[b]][, 2], use = "complete.obs")), 2),
    n1 = colSums(!is.na(W1)), n2 = colSums(!is.na(W2)), row.names = NULL)
  cr <- round(cor(W1, use = "pairwise.complete.obs"), 2); dimnames(cr) <- list(k, k)
  cat("\n=== ", tag, " ===\n", sep = ""); print(d, row.names = FALSE); cat("\nwave-1 correlations:\n"); print(cr)
  list(d = d, cr = cr) }
s1 <- summ(grab("gacr_siena.rds"), "Sample 1 (early adolescent, 27 classes)")
s2 <- summ(grab("oz_siena.rds"), "Sample 2 (middle adolescent, 18 classes)")

## ---- write the supplemental table as a Word document ----
FONT <- "Times New Roman"
apa <- function(ft) { ft <- font(ft, part = "all", fontname = FONT); ft <- fontsize(ft, part = "all", size = 10)
  ft <- border_remove(ft); ft <- hline_top(ft, part = "header", border = fp_border(color = "black", width = 1))
  ft <- hline_bottom(ft, part = "header", border = fp_border(color = "black", width = 1))
  ft <- hline_bottom(ft, part = "body", border = fp_border(color = "black", width = 1))
  ft <- padding(ft, part = "all", padding.top = 1.5, padding.bottom = 1.5); autofit(ft) }
mk <- function(x, lab) { d <- x$d
  names(d) <- c("Trait", "M wave 1", "SD wave 1", "M wave 2", "SD wave 2", "Retest r", "n wave 1", "n wave 2")
  ft <- apa(flextable(d)); ft <- bold(ft, part = "header"); ft }
mkcor <- function(x) { cr <- x$cr; cr[upper.tri(cr, diag = TRUE)] <- NA
  d <- data.frame(Trait = unname(TR[rownames(cr)]), cr, check.names = FALSE, row.names = NULL)
  names(d)[-1] <- unname(TR[colnames(cr)]); d[] <- lapply(d, function(z) ifelse(is.na(z), "", as.character(z)))
  ft <- apa(flextable(d)); bold(ft, part = "header") }
B <- fp_text(font.family = FONT, font.size = 12); BB <- update(B, bold = TRUE); BI <- update(B, italic = TRUE)
NT <- fp_text(font.family = FONT, font.size = 10)
P <- fp_par(line_spacing = 1.5)
doc <- read_docx()
doc <- body_add_fpar(doc, fpar(ftext("Supplemental Material", BB), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("Peers and Personality: Friendship Selection Without Peer Influence on the Big Five in Two Longitudinal Classroom Network Samples", BI), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("", B), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("Table S1", BB), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("Descriptive statistics and retest reliability of the ordinal trait variables entering the models", BI), fp_p = P))
doc <- body_add_flextable(doc, mk(s1, "S1")); doc <- body_add_fpar(doc, fpar(ftext("Sample 1 (early adolescent sample, 27 classes).", NT), fp_p = P))
doc <- body_add_flextable(doc, mk(s2, "S2")); doc <- body_add_fpar(doc, fpar(ftext("Sample 2 (middle adolescent sample, 18 classes).", NT), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("Note. Values describe the five-category ordinal variables that entered the stochastic actor-oriented models, not the underlying scale scores. Retest r is the Pearson correlation between the two waves among students observed at both.", NT), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("", B), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("Table S2", BB), fp_p = P))
doc <- body_add_fpar(doc, fpar(ftext("Zero-order correlations among the trait variables at the first wave", BI), fp_p = P))
doc <- body_add_flextable(doc, mkcor(s1)); doc <- body_add_fpar(doc, fpar(ftext("Sample 1.", NT), fp_p = P))
doc <- body_add_flextable(doc, mkcor(s2)); doc <- body_add_fpar(doc, fpar(ftext("Sample 2.", NT), fp_p = P))
print(doc, target = "supplemental/Supplemental_Material_Tables_S1_S2.docx")
cat("\nwritten supplemental/Supplemental_Material_Tables_S1_S2.docx\n")
