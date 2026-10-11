## 24_bayes_setup_sample2.R
## source: Results/scripts/oz_bayes_setup.R
source("00_setup.R")

## OZ sienaBayes setup -- identical specification and priors to GACR.
## 43 evaluation effects, random = density, recip, gwespFF, linear+quad x5,
## plus the 6 basic rates -> p1 = 19.  Priors as GACR: priorMu 0 except
## density -2 / recip +1; priorSigma diag 0.1 except density 0.6, recip 0.3,
## gwespFF 0.3.
suppressMessages(library(RSienaTest))
set.seed(2012)
SC <- getwd()
grp <- readRDS(file.path(OBJ, "oz_coev_group.rds"))
cat("OZ group object class:", paste(class(grp), collapse="/"), "\n")
cat("number of groups:", length(grp), "\n")
ns <- sapply(grp, function(d) d$nodeSets[[1]][[1]])
cat("group sizes:", paste(ns, collapse=" "), " | total actors:", sum(ns), "\n")
cat("names:", paste(names(grp), collapse=" "), "\n\n")

tr5 <- c("E","A","C","N","O")
eff <- getEffects(grp)
eff <- includeEffects(eff, gwespFF, inPopSqrt, outActSqrt)
eff <- includeEffects(eff, altX,  interaction1 = "sex")
eff <- includeEffects(eff, egoX,  interaction1 = "sex")
eff <- includeEffects(eff, sameX, interaction1 = "sex")
for (t in tr5) {
  eff <- includeEffects(eff, egoX, interaction1 = t)
  eff <- includeEffects(eff, altX, interaction1 = t)
  eff <- includeEffects(eff, simX, interaction1 = t)
  eff <- includeEffects(eff, avSim,    name = t, interaction1 = "net")
  eff <- includeEffects(eff, recipDeg, name = t, interaction1 = "net")
}
eff <- setEffect(eff, density, random = TRUE)
eff <- setEffect(eff, recip,   random = TRUE)
eff <- setEffect(eff, gwespFF, random = TRUE)
for (t in tr5) {
  eff <- setEffect(eff, linear, name = t, random = TRUE)
  eff <- setEffect(eff, quad,   name = t, random = TRUE)
}
E <- eff; class(E) <- "data.frame"
sel <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
reffs <- E[sel, c("name","shortName")]
p <- nrow(reffs)
p2 <- sum(E$include & !E$basicRate & !E$randomEffects & E$type == "eval")
cat("evaluation effects:", sum(E$include & E$type=="eval"),
    "| p1 (random) =", p, "| p2 (fixed) =", p2, "\n")
cat("random set:", paste(paste(reffs$name, reffs$shortName), collapse=" | "), "\n\n")

priorMu <- rep(0, p); priorSig <- diag(0.1, p)
priorMu[reffs$shortName == "density"] <- -2
priorMu[reffs$shortName == "recip"]   <-  1
diag(priorSig)[reffs$shortName == "density"] <- 0.6
diag(priorSig)[reffs$shortName == "recip"]   <- 0.3
diag(priorSig)[reffs$shortName == "gwespFF"] <- 0.3

cat("=== rough MoM init (nsub=2, n3=1000) for prevAns ===\n"); flush.console()
t0 <- proc.time()
algo0 <- sienaAlgorithmCreate(projname = "oz_bayes_init", seed = 2012,
                              nsub = 2, n3 = 1000, cond = FALSE)
ans0 <- try(siena07(algo0, data = grp, effects = eff, batch = TRUE,
                    useCluster = TRUE, nbrNodes = 8), silent = TRUE)
if (inherits(ans0, "try-error")) { cat("MoM init FAILED:\n"); print(ans0); ans0 <- NULL
} else { cat("MoM init done in", round((proc.time()-t0)[3]/60,1), "min | tconv.max",
             round(ans0$tconv.max,3), "\n"); saveRDS(ans0, "oz_bayes_init.rds") }
saveRDS(list(grp = grp, eff = eff, priorMu = priorMu, priorSig = priorSig),
        "oz_bayes_setup.rds")
cat("\nwritten oz_bayes_setup.rds\n")
