## 41_assortativity_tables.R
## source: Results/scripts_v2/homophily_out.R
source("00_setup.R")

## Wave-1 homophily: figure (same design as the posterior ridgelines) + Excel.
suppressMessages({library(ggplot2); library(ggridges); library(svglite); library(openxlsx)})
set.seed(2012)
COL <- c("Defending Victims of Bullying" = "#0072B2", "Moment of Change" = "#D55E00")
trn <- c(E = "Extraversion", A = "Agreeableness", C = "Conscientiousness",
         N = "Neuroticism",  O = "Openness")

H <- read.csv("tables/homophily_classroom.csv", stringsAsFactors = FALSE)
P <- read.csv("tables/homophily_dataset.csv",   stringsAsFactors = FALSE)
recode <- function(x) { m <- c(GACR = "Defending Victims of Bullying", OZ = "Moment of Change"); ifelse(x %in% names(m), m[x], x) }
H$dataset <- unname(recode(H$dataset)); P$dataset <- unname(recode(P$dataset))
H$dataset <- factor(H$dataset, levels = c("Defending Victims of Bullying", "Moment of Change"))
H$lab <- factor(trn[H$trait], levels = rev(trn))

## pooled estimate with 80% and 95% bootstrap CIs over classrooms
boot <- function(d, B = 4000) {
  w <- d$ties; r <- d$r_excess; k <- nrow(d)
  bs <- replicate(B, { i <- sample(k, k, TRUE); sum(r[i] * w[i]) / sum(w[i]) })
  c(M = sum(r * w) / sum(w), setNames(quantile(bs, c(.025, .10, .90, .975)),
                                      c("lo", "lo80", "hi80", "hi")))
}
S <- do.call(rbind, lapply(split(H, list(H$dataset, H$trait), drop = TRUE), function(d)
       data.frame(dataset = d$dataset[1], trait = d$trait[1], t(boot(d)))))
S$lab <- factor(trn[S$trait], levels = levels(H$lab))
S$y <- as.numeric(S$lab)
S$ynudge <- S$y - ifelse(S$dataset == "Defending Victims of Bullying", 0.16, 0.30)

p <- ggplot(H, aes(x = r_excess, y = lab, fill = dataset, colour = dataset)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey45", linewidth = .4) +
  geom_density_ridges(alpha = .42, scale = .80, rel_min_height = .01,
                      linewidth = .35, colour = "white",
                      jittered_points = TRUE, point_size = 1.15, point_alpha = .55,
                      point_shape = 16,
                      position = position_points_jitter(width = 0, height = .045)) +
  geom_linerange(data = S, inherit.aes = FALSE,
                 aes(xmin = lo, xmax = hi, y = ynudge, colour = dataset), linewidth = .55) +
  geom_linerange(data = S, inherit.aes = FALSE,
                 aes(xmin = lo80, xmax = hi80, y = ynudge, colour = dataset), linewidth = 1.5) +
  geom_point(data = S, inherit.aes = FALSE,
             aes(x = M, y = ynudge, colour = dataset),
             shape = 21, size = 2.1, fill = "white", stroke = .9) +
  scale_fill_manual(values = COL, name = NULL) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(title = "Wave-1 personality homophily in classroom friendship networks",
       subtitle = paste0(
         "Newman assortativity of each trait over observed ties, centred on a within-classroom permutation null\n",
         "(1,000 permutations). Positive = friends resemble each other more than chance. Ridges and dots show the\n",
         "individual classrooms, 27 in Defending Victims of Bullying and 18 in Moment of Change.\n",
         "The interval below each ridge is the tie-weighted pooled estimate."),
       x = "excess assortativity  (r observed - r expected under random mixing)", y = NULL,
       caption = paste0("Point = pooled estimate; thick bar = 80% CI; thin bar = 95% CI, ",
                        "bootstrapped over classrooms. Dots are individual classrooms.")) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = .3),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey94", linewidth = .3),
        axis.text.y  = element_text(size = 10, colour = "grey15"),
        axis.text.x  = element_text(size = 9,  colour = "grey30"),
        plot.title    = element_text(face = "bold", size = 13.5, colour = "grey10"),
        plot.subtitle = element_text(size = 9.5, colour = "grey40", margin = margin(b = 10)),
        plot.caption  = element_text(size = 8, colour = "grey50", hjust = 0),
        legend.position = "top", legend.justification = "left",
        legend.key.height = unit(10, "pt"),
        plot.margin = margin(12, 16, 10, 10)) +
  coord_cartesian(clip = "off")
ggsave("figures/fig_homophily.png", p, width = 9, height = 5.6, dpi = 400, bg = "white")
ggsave("figures/fig_homophily.svg", p, width = 9, height = 5.6, device = svglite, bg = "white")
cat("written fig_homophily.png / .svg\n")

## ---------------- Excel ----------------
P$trait_name <- trn[P$trait]
P <- P[, c("dataset","trait","trait_name","classes","actors","ties",
           "r_excess_pooled","CI_low","CI_high","r_raw_pooled","null_mean_pooled",
           "r_sd_between","moran_pooled","simgap_pooled",
           "z_stouffer","p_stouffer","classes_p05")]
H2 <- H[, c("dataset","class","trait","n","ties","density","r_assort","null_mean",
            "null_sd","r_excess","z_perm","p_perm","moran_I","sim_tie","sim_all","sim_gap")]

wb <- createWorkbook(); modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
hdr <- createStyle(textDecoration="bold", fgFill="#1F3864", fontColour="white",
                   halign="center", valign="center", border="TopBottomLeftRight", wrapText=TRUE)
ttl <- createStyle(textDecoration="bold", fontSize=13, fontColour="#1F3864")
sub <- createStyle(fontSize=9, fontColour="#555555")
n4  <- createStyle(numFmt="0.0000"); n3 <- createStyle(numFmt="0.000"); i0 <- createStyle(numFmt="0")
sig <- createStyle(fgFill="#E2EFDA", textDecoration="bold")
wrapL <- createStyle(halign="left", valign="top", wrapText=TRUE)

addWorksheet(wb, "README", tabColour="#C00000")
rd <- data.frame(Item=c("Question","Data","Measure (primary)","Why not E-I or Coleman",
  "Permutation centring","Significance","Measure (check)","Measure (model-linked)",
  "Pooling","Between-classroom SD","Missing data","Structural zeros","Seed"),
  Detail=c(
  "How much do adolescents at wave 1 nominate friends with similar Big Five scores?",
  "Wave-1 networks and the 5-category trait variables used in the SAOM. Defending Victims of Bullying 27 classrooms / 632 pupils; Moment of Change 18 classrooms / 543 pupils.",
  "Newman's assortativity coefficient: the Pearson correlation between sender and receiver trait values taken over the observed directed ties.",
  "Those indices are built for categorical attributes. The traits are 5-category ordinal, and dichotomising them would discard the ordering that the SAOM itself uses.",
  "In small directed networks the edgewise correlation is biased downward - the permutation null averages about -0.05 to -0.06 here, not 0. The reported r_excess is the observed r minus the mean of its own within-classroom permutation null, which removes that bias. r_raw_pooled shows the uncentred value for transparency.",
  "1,000 within-classroom permutations of the trait over actors, holding the network fixed. z_perm is the standardised excess; classroom p-values are combined across classrooms with Stouffer's method.",
  "Moran's I with the adjacency matrix as weights - the classic network autocorrelation statistic. Same ordering of traits.",
  "sim_gap: mean RSiena similarity (1 - |xi - xj| / range) over tied dyads minus the mean over all dyads. This is the quantity the simX effect in the SAOM operates on.",
  "Classroom estimates pooled with weights proportional to the number of ties; 95% CI bootstrapped over classrooms (4,000 resamples).",
  "r_sd_between is the raw SD of the classroom-level estimates - homophily itself varies a lot from classroom to classroom.",
  "Actors with a missing trait score at wave 1 are dropped from that trait's calculation. Classrooms with fewer than 5 usable actors or ties, or with no trait variance, are skipped.",
  "Network codes 10 (structural zero) and 11 (structural one) are read as absent and present respectively; NA stays missing.",
  "set.seed(2012)"), stringsAsFactors=FALSE)
writeData(wb,"README","Wave-1 personality homophily", startRow=1); addStyle(wb,"README",ttl,1,1)
writeData(wb,"README",paste("Generated",format(Sys.Date())),startRow=2); addStyle(wb,"README",sub,2,1)
writeData(wb,"README",rd,startRow=4,headerStyle=hdr)
addStyle(wb,"README",createStyle(textDecoration="bold",valign="top"),5:(4+nrow(rd)),1,gridExpand=TRUE)
addStyle(wb,"README",wrapL,5:(4+nrow(rd)),2,gridExpand=TRUE)
setColWidths(wb,"README",1,26); setColWidths(wb,"README",2,115)

addWorksheet(wb, "Dataset_level", tabColour="#1F3864")
writeData(wb,"Dataset_level","Pooled wave-1 homophily by dataset and trait",startRow=1)
addStyle(wb,"Dataset_level",ttl,1,1)
writeData(wb,"Dataset_level",
  "r_excess_pooled is the primary measure. Green = 95% CI excludes zero.",startRow=2)
addStyle(wb,"Dataset_level",sub,2,1)
writeData(wb,"Dataset_level",P,startRow=4,headerStyle=hdr,borders="surrounding")
rows <- 5:(4+nrow(P))
for(cc in c(7:14)) addStyle(wb,"Dataset_level",n4,rows,cc,gridExpand=TRUE,stack=TRUE)
addStyle(wb,"Dataset_level",n3,rows,15,gridExpand=TRUE,stack=TRUE)
for(i in which(P$CI_low > 0)) addStyle(wb,"Dataset_level",sig,4+i,7,stack=TRUE)
freezePane(wb,"Dataset_level",firstActiveRow=5,firstActiveCol=4)
setColWidths(wb,"Dataset_level",1:ncol(P),"auto")

addWorksheet(wb, "Classroom_level", tabColour="#1F3864")
writeData(wb,"Classroom_level","Wave-1 homophily, every classroom x trait",startRow=1)
addStyle(wb,"Classroom_level",ttl,1,1)
writeData(wb,"Classroom_level",
  "r_excess = observed assortativity minus its permutation-null mean. Green = permutation p < .05.",startRow=2)
addStyle(wb,"Classroom_level",sub,2,1)
writeData(wb,"Classroom_level",H2,startRow=4,headerStyle=hdr,borders="surrounding")
rows <- 5:(4+nrow(H2))
for(cc in c(6:12,13:16)) addStyle(wb,"Classroom_level",n4,rows,cc,gridExpand=TRUE,stack=TRUE)
for(i in which(H2$p_perm < .05)) addStyle(wb,"Classroom_level",sig,4+i,10,stack=TRUE)
freezePane(wb,"Classroom_level",firstActiveRow=5,firstActiveCol=4)
setColWidths(wb,"Classroom_level",1:ncol(H2),"auto")

saveWorkbook(wb,"tables/Homophily_wave1.xlsx",overwrite=TRUE)
cat("written Homophily_wave1.xlsx | sheets:", paste(names(wb),collapse=" | "), "\n")
