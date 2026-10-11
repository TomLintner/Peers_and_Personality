## 25_bayes_rebuild_sample2.R
## source: Results/scripts/oz_rebuild.R
source("00_setup.R")

## Rebuild the OZ group for likelihood-based estimation:
##  (a) forbidden changes (structural zero at T1 -> tie at T2) set to NA
##  (b) allowOnly = FALSE on the network (as in the GACR Bayes rebuild)
##  (c) covariate named "sex" from sexBirthFemale, matching GACR effect labels
## Classes 0503 and 0602 stay excluded (degenerate: 1 and 0 changers).
suppressMessages(library(RSienaTest))
set.seed(2012)
tr5 <- c("E","A","C","N","O")
O <- readRDS(file.path(OBJ, "oz_siena.rds"))
keep <- names(readRDS(file.path(OBJ, "oz_coev_group.rds")))
cat("classes:", length(keep), "\n\n=== forbidden changes ===\n")
tot <- 0
objs <- lapply(keep, function(cl) {
  z <- O$data[[cl]]
  n1 <- z$net[,,1]; n2 <- z$net[,,2]
  imp <- (n1==10 & n2==1) | (n1==11 & n2==0); imp[is.na(imp)] <- FALSE
  if (any(imp)) { cat(sprintf("  %s: %d cell(s) -> NA\n", cl, sum(imp))); tot <<- tot + sum(imp)
                  n1[imp] <- NA }
  args <- list(net = sienaDependent(array(c(n1, n2), dim = c(dim(n1), 2)), allowOnly = FALSE))
  for (t in tr5) args[[t]] <- sienaDependent(z$beh[[t]], type = "behavior", allowOnly = FALSE)
  args$sex <- coCovar(z$sexF)
  do.call(sienaDataCreate, args)
})
names(objs) <- keep
cat("total cells set to NA:", tot, "\n")
grp <- sienaGroupCreate(objs)
cat("element class:", paste(class(grp[[1]]), collapse="/"), "| groups:", length(grp), "\n")
ic <- try(unlist(lapply(grp, RSienaTest:::checkImpossibleChanges)), silent = TRUE)
if (inherits(ic, "try-error")) { cat("checkImpossibleChanges ERRORED\n"); print(ic) } else {
  cat("impossible changes per group:", sum(ic), "total\n")
  if (sum(ic) > 0) print(ic[ic > 0]) }
saveRDS(grp, "oz_bayes_group_fixed.rds")

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
for (t in tr5) { eff <- setEffect(eff, linear, name = t, random = TRUE)
                 eff <- setEffect(eff, quad,   name = t, random = TRUE) }
E <- eff; class(E) <- "data.frame"
sel <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
reffs <- E[sel, c("name","shortName")]; p <- nrow(reffs)
p2 <- sum(E$include & !E$basicRate & !E$randomEffects & E$type=="eval")
cat("\neval effects:", sum(E$include & E$type=="eval"), "| p1:", p, "| p2:", p2, "\n")
priorMu <- rep(0, p); priorSig <- diag(0.1, p)
priorMu[reffs$shortName=="density"] <- -2; priorMu[reffs$shortName=="recip"] <- 1
diag(priorSig)[reffs$shortName=="density"] <- 0.6
diag(priorSig)[reffs$shortName=="recip"]   <- 0.3
diag(priorSig)[reffs$shortName=="gwespFF"] <- 0.3
saveRDS(list(grp=grp, eff=eff, priorMu=priorMu, priorSig=priorSig), "oz_bayes_setup.rds")
cat("written oz_bayes_setup.rds\n")
