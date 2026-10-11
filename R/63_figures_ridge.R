## 63_figures_ridge.R
## source: Results/scripts_v2/ridge_plots.R
source("00_setup.R")

## Publication-ready posterior ridge plots, Defending Victims of Bullying vs Moment of Change overlaid.
## Draws are the burn-in-pruned chains: Defending Victims of Bullying 2650-900 = 1750; Moment of Change 6000-1200 = 4800.
## Every figure is written as both a 400 dpi PNG and a vector SVG.
suppressMessages({library(RSienaTest); library(ggplot2); library(ggridges); library(svglite)})

COL <- c("Defending Victims of Bullying" = "#0072B2", "Moment of Change" = "#D55E00")     # Okabe-Ito, colourblind-safe
CAP <- paste0("Posterior densities from the burn-in-pruned chains: Defending Victims of Bullying 1,750 draws,\n",
              "Moment of Change 4,800 draws. Point = posterior mean; thick bar = 80% CrI; thin bar = 95% CrI.")

## ---- extract draws -------------------------------------------------------
draws <- function(setup, pat, kmax, burn, tag) {
  s <- readRDS(setup); E <- s$eff; class(E) <- "data.frame"
  selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
  fxi  <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
  MU <- NULL; ET <- NULL; SG <- NULL
  for (k in 1:kmax) {
    f <- sprintf(pat, k); if (!file.exists(f)) next
    g <- readRDS(f); ok <- !is.na(g$ThinPosteriorMu[, 1]); if (k == 1) ok[1:100] <- FALSE
    MU <- rbind(MU, g$ThinPosteriorMu[ok, , drop = FALSE])
    ET <- rbind(ET, g$ThinPosteriorEta[ok, , drop = FALSE])
    sg <- g$ThinPosteriorSigma
    if (!is.null(sg))
      SG <- rbind(SG, t(apply(sg[ok, , , drop = FALSE], 1, function(m) sqrt(diag(m)))))
  }
  MU <- MU[-(1:burn), , drop = FALSE]; ET <- ET[-(1:burn), , drop = FALSE]
  SG <- SG[-(1:burn), , drop = FALSE]
  colnames(MU) <- colnames(SG) <- E$effectName[selR]; colnames(ET) <- E$effectName[fxi]
  bind <- function(M, blk) data.frame(dataset = tag, block = blk,
      param = rep(colnames(M), each = nrow(M)), value = as.vector(M),
      stringsAsFactors = FALSE)
  rbind(bind(MU, "mu"), bind(ET, "eta"), bind(SG, "sigma"))
}
ALL <- rbind(draws(file.path(OBJ, "gacr_bayes_setup2.rds"), file.path(CHN, "bayes12h_stage%d.rds"),  9,  900, "Defending Victims of Bullying"),
             draws(file.path(OBJ, "oz_bayes_setup.rds"),    file.path(CHN, "ozbayes_stage%d.rds"),  20, 1200, "Moment of Change"))
ALL$dataset <- factor(ALL$dataset, levels = c("Defending Victims of Bullying", "Moment of Change"))
D  <- ALL[ALL$block != "sigma", ]
SD <- ALL[ALL$block == "sigma", ]

summ <- function(d) {
  s <- aggregate(value ~ dataset + param, d, function(x)
        c(M = mean(x), lo = unname(quantile(x, .025)), hi = unname(quantile(x, .975)),
          lo80 = unname(quantile(x, .10)), hi80 = unname(quantile(x, .90))))
  data.frame(s[, 1:2], as.data.frame(s$value))
}

## ---- plot builder --------------------------------------------------------
ridge <- function(dat, params, labels, title, subtitle, xlab, base,
                  w = 9, h = 6, facet = NULL, facetlab = NULL, freescale = "free_x") {
  d <- dat[dat$param %in% params, ]
  d$lab <- factor(labels[d$param], levels = rev(unique(labels[params])))
  S <- summ(d); S$lab <- factor(labels[S$param], levels = levels(d$lab))
  S$y <- as.numeric(S$lab)
  if (!is.null(facet)) {
    d$fc <- factor(facetlab[facet[d$param]], levels = facetlab)
    S$fc <- factor(facetlab[facet[S$param]], levels = facetlab)
    ## With scales = "free" ggplot drops unused y levels per panel and
    ## re-indexes them, so the interval rows must be numbered within panel.
    if (freescale == "free")
      S <- do.call(rbind, lapply(split(S, S$fc), function(z) {
             z$y <- as.numeric(droplevels(z$lab)); z }))
  }
  S$ynudge <- S$y - ifelse(S$dataset == "Defending Victims of Bullying", 0.16, 0.30)

  p <- ggplot(d, aes(x = value, y = lab, fill = dataset, colour = dataset)) +
    geom_vline(xintercept = 0, linetype = "22", colour = "grey45", linewidth = .4) +
    geom_density_ridges(alpha = .48, scale = .88, rel_min_height = .004,
                        linewidth = .35, colour = "white", panel_scaling = TRUE) +
    geom_linerange(data = S, inherit.aes = FALSE,
                   aes(xmin = lo, xmax = hi, y = ynudge, colour = dataset),
                   linewidth = .55, alpha = .95) +
    geom_linerange(data = S, inherit.aes = FALSE,
                   aes(xmin = lo80, xmax = hi80, y = ynudge, colour = dataset),
                   linewidth = 1.5, alpha = .95) +
    geom_point(data = S, inherit.aes = FALSE,
               aes(x = M, y = ynudge, colour = dataset),
               shape = 21, size = 2.1, fill = "white", stroke = .9) +
    scale_fill_manual(values = COL, name = NULL) +
    scale_colour_manual(values = COL, name = NULL) +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL, caption = CAP) +
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
          strip.text = element_text(face = "bold", size = 10.5, colour = "grey20"),
          plot.margin = margin(12, 16, 10, 10)) +
    coord_cartesian(clip = "off")
  if (!is.null(facet)) p <- p + facet_wrap(~fc, scales = freescale, nrow = 1)
  ggsave(paste0("figures/", base, ".png"), p, width = w, height = h, dpi = 400, bg = "white")
  ggsave(paste0("figures/", base, ".svg"), p, width = w, height = h, device = svglite, bg = "white")
  cat("written", base, ".png / .svg\n"); invisible(p)
}

tr  <- c("E", "A", "C", "N", "O")
trn <- c(E = "Extraversion", A = "Agreeableness", C = "Conscientiousness",
         N = "Neuroticism",  O = "Openness")

## 1. personality -> ties ---------------------------------------------------
selp <- c(paste(tr, "alter"), paste(tr, "ego"), paste(tr, "similarity"))
ridge(D, selp, setNames(rep(trn[tr], 3), selp),
      "Effects of personality on tie formation",
      "Fixed effects (eta) in the network function. Alter = being nominated; ego = nominating; similarity = homophily.",
      "posterior value (log-odds)", "fig_selection", w = 10, h = 5.4,
      facet = setNames(rep(c("alter", "ego", "similarity"), each = 5), selp),
      facetlab = c(alter = "Alter (being nominated)",
                   ego = "Ego (nominating others)",
                   similarity = "Similarity (homophily)"))

## 2. peers -> personality --------------------------------------------------
infp <- c(paste(tr, "average similarity"), paste(tr, "reciprocated degree"))
ridge(D, infp, setNames(rep(trn[tr], 2), infp),
      "Effects of the network on personality traits",
      "Fixed effects (eta) in the behaviour function: assimilation to friends, and the effect of mutual ties on the trait.",
      "posterior value", "fig_influence", w = 8.5, h = 4.8,
      facet = setNames(rep(c("avSim", "recipDeg"), each = 5), infp),
      facetlab = c(avSim = "Peer influence (average similarity)",
                   recipDeg = "Number of mutual ties (reciprocated degree)"))

## 3. structure -------------------------------------------------------------
strp <- c("outdegree (density)", "reciprocity", "GWESP I -> K -> J (#)",
          "indegree-popularity (sqrt)", "outdegree-activity (sqrt)",
          "same sex", "sex alter", "sex ego")
ridge(D, strp,
      setNames(c("Density", "Reciprocity", "Transitive closure (GWESP)",
                 "Indegree-popularity", "Outdegree-activity",
                 "Gender homophily", "Gender (alter)", "Gender (ego)"), strp),
      "Structural and gender effects on tie formation",
      "Density, reciprocity and closure vary across classrooms (random effects);\nthe degree and gender terms are held constant across classrooms (fixed effects).",
      "posterior value (log-odds)", "fig_structure", w = 9.5, h = 5.6,
      facet = setNames(c(rep("endo", 3), rep("cov", 5)), strp),
      facetlab = c(endo = "Endogenous network structure (random)",
                   cov  = "Degree and gender (fixed)"))

## 4. baseline trait dynamics ----------------------------------------------
behp <- c(paste(tr, "linear shape"), paste(tr, "quadratic shape"))
ridge(D, behp, setNames(rep(trn[tr], 2), behp),
      "Baseline dynamics of the personality traits",
      "Random effects (population means) in the behaviour function. Negative quadratic = scores regress toward the classroom mean.",
      "posterior value", "fig_behaviour", w = 8.5, h = 4.8,
      facet = setNames(rep(c("lin", "quad"), each = 5), behp),
      facetlab = c(lin = "Linear shape (overall drift)",
                   quad = "Quadratic shape (self-correction)"))

## 5. between-classroom heterogeneity --------------------------------------
sdp <- c(paste("rate", tr, "(period 1)"),
         "outdegree (density)", "reciprocity", "GWESP I -> K -> J (#)",
         paste(tr, "linear shape"), paste(tr, "quadratic shape"))
sdl <- setNames(c(paste(trn[tr], "rate"),
                  "Density", "Reciprocity", "Transitive closure (GWESP)",
                  paste(trn[tr], "linear"), paste(trn[tr], "quadratic")), sdp)
sdf <- setNames(c(rep("rate", 5), rep("net", 3), rep("shape", 10)), sdp)
ridge(SD, sdp, sdl,
      "How much do classrooms differ? Between-classroom standard deviations",
      paste0("Posterior of the square root of the diagonal of Sigma: the SD of each random effect across classrooms.\n",
             "Larger values mean the effect is less transferable from one classroom to another. The network rate parameter\n",
             "is omitted for scale (between-classroom SD 2.54 in Defending Victims of Bullying,\n",
             "5.40 in Moment of Change)."),
      "between-classroom SD", "fig_heterogeneity", w = 11, h = 6.4,
      facet = sdf,
      facetlab = c(rate = "Rate parameters", net = "Network structure",
                   shape = "Trait shape"),
      freescale = "free")
cat("done\n")
