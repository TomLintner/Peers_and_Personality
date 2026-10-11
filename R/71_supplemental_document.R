## 71_supplemental_document.R
## source: EJP_submission/make_supplemental.R
source("00_setup.R")

## Assembles the whole Supplemental Material into one captioned Word document:
##   Table S1  descriptive statistics and retest reliability of the ordinal trait variables
##   Table S2  zero-order correlations among the trait variables at the first wave
##   Table S3  full posterior summary, Sample 1
##   Table S4  full posterior summary, Sample 2
##   Table S5  classroom-level wave-1 assortativity, Sample 1
##   Table S6  classroom-level wave-1 assortativity, Sample 2
##   Figures S1-S6  convergence diagnostics (MDS of the group-level deviations; trace plots)
suppressMessages({library(RSienaTest); library(officer); library(flextable)})
BASE <- RAW
setwd(file.path(BASE, "EJP_submission"))
S <- "supplemental"
TR <- c(E = "Extraversion", A = "Agreeableness", C = "Conscientiousness",
        N = "Neuroticism", O = "Openness to experience")
S1LAB <- "Defending Victims of Bullying"   # Sample 1, 27 classes, early adolescent, six months
S2LAB <- "Moment of Change"                # Sample 2, 18 classes, middle adolescent, eight weeks

## ---------------------------------------------------------------- data
grab <- function(f) { G <- readRDS(file.path(BASE, "Results", "objects", f)); tr <- c("E","A","C","N","O")
  out <- lapply(tr, function(b) do.call(rbind, lapply(G$data, function(g) g$beh[[b]]))); names(out) <- tr; out }
summ <- function(L) {
  k <- names(L); W1 <- sapply(L, function(m) m[, 1]); W2 <- sapply(L, function(m) m[, 2])
  d <- data.frame(Trait = unname(TR[k]),
    M1 = round(colMeans(W1, na.rm = TRUE), 2), SD1 = round(apply(W1, 2, sd, na.rm = TRUE), 2),
    M2 = round(colMeans(W2, na.rm = TRUE), 2), SD2 = round(apply(W2, 2, sd, na.rm = TRUE), 2),
    rtt = round(sapply(k, function(b) cor(L[[b]][, 1], L[[b]][, 2], use = "complete.obs")), 2),
    n1 = colSums(!is.na(W1)), n2 = colSums(!is.na(W2)), row.names = NULL)
  cr <- round(cor(W1, use = "pairwise.complete.obs"), 2); dimnames(cr) <- list(k, k)
  list(d = d, cr = cr) }
s1 <- summ(grab("gacr_siena.rds")); s2 <- summ(grab("oz_siena.rds"))

post <- function(lab) {
  f <- if (lab == S1LAB) "bayes_defending_victims.csv" else "bayes_moment_of_change.csv"
  x <- read.csv(file.path(S, f), stringsAsFactors = FALSE)
  x$block <- ifelse(grepl("random", x$block), "Varying across classrooms", "Constant across classrooms")
  d <- data.frame(Block = x$block, Parameter = x$parameter,
    M = sprintf("%.3f", x$M), SD = sprintf("%.3f", x$SD),
    CI = sprintf("[%.3f, %.3f]", x$CI_low, x$CI_high),
    SDb = ifelse(is.na(x$SD_between), "", sprintf("%.3f", x$SD_between)),
    PD = sprintf("%.3f", x$PD), Rhat = sprintf("%.3f", x$Rhat), ESS = format(round(x$ESS)),
    row.names = NULL)
  names(d) <- c("Block", "Parameter", "M", "SD", "95% CI", "SD between classrooms",
                "PD", "R-hat", "ESS")
  d }
p1 <- post(S1LAB); p2 <- post(S2LAB)

hom <- read.csv(file.path(S, "homophily_classroom.csv"), stringsAsFactors = FALSE)
assort <- function(lab) {
  x <- hom[hom$dataset == lab, ]
  x <- x[order(x$class, match(x$trait, c("E","A","C","N","O"))), ]
  d <- data.frame(Class = x$class, Trait = unname(TR[x$trait]), n = x$n, Ties = x$ties,
    Density = sprintf("%.3f", x$density), r = sprintf("%.3f", x$r_assort),
    Null = sprintf("%.3f (%.3f)", x$null_mean, x$null_sd),
    z = sprintf("%.2f", x$z_perm), p = sprintf("%.3f", x$p_perm), row.names = NULL)
  names(d) <- c("Class", "Trait", "n", "Ties", "Density", "r", "Null M (SD)", "z", "p")
  d }
a1 <- assort(S1LAB); a2 <- assort(S2LAB)

## ---------------------------------------------------------------- styling
FONT <- "Times New Roman"
B  <- fp_text(font.family = FONT, font.size = 12)
BB <- update(B, bold = TRUE); BI <- update(B, italic = TRUE)
NT <- fp_text(font.family = FONT, font.size = 10)
P  <- fp_par(text.align = "left", line_spacing = 1.5)
P1 <- fp_par(text.align = "left", line_spacing = 1.15)
apa <- function(ft, size = 9) {
  ft <- font(ft, part = "all", fontname = FONT); ft <- fontsize(ft, part = "all", size = size)
  ft <- border_remove(ft)
  ft <- hline_top(ft, part = "header", border = fp_border(color = "black", width = 1))
  ft <- hline_bottom(ft, part = "header", border = fp_border(color = "black", width = 1))
  ft <- hline_bottom(ft, part = "body", border = fp_border(color = "black", width = 1))
  ft <- padding(ft, part = "all", padding.top = 1.5, padding.bottom = 1.5)
  ft <- valign(ft, valign = "top", part = "all"); bold(ft, part = "header") }
## autofit() sizes columns to their content and happily overruns the page, so every table is then
## laid out to exactly the printable width of its section
fit <- function(ft) set_table_properties(ft, layout = "autofit", width = 1)
tbl <- function(d, size = 9) fit(apa(autofit(flextable(d)), size))
tblw <- function(d, size, j, wj) fit(width(apa(autofit(flextable(d)), size), j = j, width = wj))
hd  <- function(doc, t) body_add_fpar(doc, fpar(ftext(t, BB), fp_p = P))
it  <- function(doc, t) body_add_fpar(doc, fpar(ftext(t, BI), fp_p = P))
nn  <- function(doc, t) body_add_fpar(doc, fpar(ftext(t, NT), fp_p = P1))
sp  <- function(doc)    body_add_fpar(doc, fpar(ftext("", B), fp_p = P1))

## image sizing from the PNG header, so the aspect ratio is preserved
imgdim <- function(f, w = 6.3) {
  h <- tryCatch({ d <- dim(png::readPNG(f)); w * d[1] / d[2] }, error = function(e) w * 0.72)
  c(w = w, h = min(h, 8.6)) }
figs <- function(doc, file, num, cap, note = "") {
  z <- imgdim(file.path(S, file))
  doc <- hd(doc, paste0("Figure S", num))
  doc <- it(doc, cap)
  doc <- body_add_img(doc, file.path(S, file), width = z["w"], height = z["h"])
  if (nzchar(note)) doc <- nn(doc, note)
  sp(doc) }

mk <- function(x) { d <- x$d
  names(d) <- c("Trait", "M wave 1", "SD wave 1", "M wave 2", "SD wave 2", "Retest r", "n wave 1", "n wave 2")
  tbl(d, 10) }
mkcor <- function(x) { cr <- x$cr; cr[upper.tri(cr, diag = TRUE)] <- NA
  d <- data.frame(Trait = unname(TR[rownames(cr)]), cr, check.names = FALSE, row.names = NULL)
  names(d)[-1] <- unname(TR[colnames(cr)])
  d[] <- lapply(d, function(z) ifelse(is.na(z), "", as.character(z)))
  tbl(d, 10) }

PORT <- prop_section(page_size = page_size(width = 8.5, height = 11, orient = "portrait"),
  page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1, header = .5, footer = .5, gutter = 0))
LAND <- prop_section(page_size = page_size(width = 11, height = 8.5, orient = "landscape"),
  page_margins = page_mar(top = .8, bottom = .8, left = .8, right = .8, header = .4, footer = .4, gutter = 0))

## ---------------------------------------------------------------- document
## the default (last) section must be US Letter portrait too, otherwise officer leaves the
## trailing block on its own A4 default and the figures print on a different page size
d <- read_docx(); d <- body_set_default_section(d, PORT)
d <- hd(d, "Supplemental Material")
d <- it(d, paste("Peers and Personality: Friendship Selection Without Peer Influence on the Big Five",
                 "in Two Longitudinal Classroom Network Samples"))
d <- sp(d)
d <- nn(d, paste("Sample 1 is the early adolescent sample: 27 classes followed over six months.",
                 "Sample 2 is the middle adolescent sample: 18 classes followed over eight weeks.",
                 "Tables S1 and S2 describe the trait variables as they entered the models.",
                 "Tables S3 and S4 give the complete posterior summaries, including every parameter",
                 "omitted from the tables in the main text. Tables S5 and S6 give the classroom-level",
                 "assortativity values that Figure 1 of the main text summarizes. Figures S1 to S6",
                 "are the convergence diagnostics."))
d <- sp(d)

## Table S1
d <- hd(d, "Table S1")
d <- it(d, "Descriptive statistics and retest reliability of the ordinal trait variables entering the models")
d <- body_add_flextable(d, mk(s1)); d <- nn(d, "Sample 1 (early adolescent sample, 27 classes).")
d <- sp(d)
d <- body_add_flextable(d, mk(s2)); d <- nn(d, "Sample 2 (middle adolescent sample, 18 classes).")
d <- nn(d, paste("Note. Values describe the five-category ordinal variables that entered the stochastic",
                 "actor-oriented models, not the underlying scale scores. Retest r is the Pearson",
                 "correlation between the two waves among students observed at both waves."))
d <- sp(d)

## Table S2
d <- hd(d, "Table S2")
d <- it(d, "Zero-order correlations among the trait variables at the first wave")
d <- body_add_flextable(d, mkcor(s1)); d <- nn(d, "Sample 1.")
d <- sp(d)
d <- body_add_flextable(d, mkcor(s2)); d <- nn(d, "Sample 2.")
d <- nn(d, paste("Note. Pearson correlations among the five-category ordinal variables at the first wave,",
                 "pairwise complete."))

## Tables S3 and S4, landscape
d <- body_end_block_section(d, block_section(PORT))
d <- hd(d, "Table S3")
d <- it(d, "Complete posterior summary, Sample 1 (early adolescent sample, 27 classes)")
d <- body_add_flextable(d, tblw(p1, 8, 2, 2.6))
d <- nn(d, paste("Note. Posterior means, standard deviations and 95% credible intervals from the Bayesian",
                 "multilevel stochastic actor-oriented model. Parameters in the first block vary across",
                 "classrooms and are summarized by the population mean and the between-classroom standard",
                 "deviation; parameters in the second block are constant across classrooms. PD is the",
                 "probability of direction, the posterior probability that the parameter has the sign of its",
                 "mean. R-hat is the rank-normalized split R-hat and ESS the effective sample size."))
d <- sp(d)
d <- hd(d, "Table S4")
d <- it(d, "Complete posterior summary, Sample 2 (middle adolescent sample, 18 classes)")
d <- body_add_flextable(d, tblw(p2, 8, 2, 2.6))
d <- nn(d, "Note. As Table S3.")

## Tables S5 and S6
d <- hd(d, "Table S5")
d <- it(d, "Classroom-level wave-1 assortativity by trait, Sample 1")
d <- body_add_flextable(d, tbl(a1, 7.5))
d <- nn(d, paste("Note. r is Newman's assortativity coefficient for the ordinal trait variable over the",
                 "wave-1 friendship network of each class. Null M (SD) is the mean and standard deviation of",
                 "r across 1,000 within-class permutations of the trait values, which hold the network fixed.",
                 "z is the standardized deviation of the observed coefficient from that null distribution and",
                 "p the corresponding two-sided permutation p value. Classes are identified by the",
                 "project's internal codes."))
d <- sp(d)
d <- hd(d, "Table S6")
d <- it(d, "Classroom-level wave-1 assortativity by trait, Sample 2")
d <- body_add_flextable(d, tbl(a2, 7.5))
d <- nn(d, "Note. As Table S5. Class sizes differ slightly between traits because of item-level missingness.")
d <- body_end_block_section(d, block_section(LAND))

## Figures
d <- figs(d, "defending_victims_mds.png", 1,
  "Multidimensional scaling of the classroom-level parameter deviations, Sample 1",
  paste("Note. Each point is one classroom, positioned by the posterior means of its deviations from the",
        "population parameters. The display is diagnostic: it shows whether any classroom is extreme enough",
        "to drive the population estimates."))
d <- figs(d, "moment_of_change_mds.png", 2,
  "Multidimensional scaling of the classroom-level parameter deviations, Sample 2",
  "Note. As Figure S1.")
d <- figs(d, "defending_victims_trace_mu.png", 3,
  "Trace plots of the population means, Sample 1",
  paste("Note. Post-burn-in draws for every parameter that varies across classrooms. The chains are",
        "shown after the burn-in reported in the method."))
d <- figs(d, "defending_victims_trace_eta.png", 4,
  "Trace plots of the parameters constant across classrooms, Sample 1", "Note. As Figure S3.")
d <- figs(d, "moment_of_change_trace_mu.png", 5,
  "Trace plots of the population means, Sample 2", "Note. As Figure S3.")
d <- figs(d, "moment_of_change_trace_eta.png", 6,
  "Trace plots of the parameters constant across classrooms, Sample 2", "Note. As Figure S3.")

print(d, target = file.path(S, "Supplemental_Material.docx"))
cat("written supplemental/Supplemental_Material.docx\n")
cat("S3 rows:", nrow(p1), "| S4 rows:", nrow(p2),
    "| S5 rows:", nrow(a1), "| S6 rows:", nrow(a2), "\n")
