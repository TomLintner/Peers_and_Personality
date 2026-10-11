## 26_bayes_sample2.R
## source: Results/scripts/oz_all.R
source("00_setup.R")

## OZ: MoM init on the repaired group, then staged sienaBayes.
## Settings identical to the GACR chain (priorKappa 1.0, priorSigEta 1,
## nSampVarying 3, nSampConst 10, nrunMHBatches 20, 8 cores).
suppressMessages(library(RSienaTest))
set.seed(2012)
SC <- WORK
WD <- file.path(SC, "ozbayes"); dir.create(WD, showWarnings = FALSE); setwd(WD)
s <- readRDS(file.path(OBJ, "oz_bayes_setup.rds"))
E <- s$eff; class(E) <- "data.frame"
p2 <- sum(E$include & !E$basicRate & !E$randomEffects & E$type == "eval")
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
labMu  <- paste(E$name[selR], E$shortName[selR])
fxi <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labEta <- paste(E$name[fxi], E$shortName[fxi])
cat("groups:", length(s$grp), "| p1:", sum(selR), "| p2:", p2, "\n"); flush.console()

IF <- file.path(SC, "oz_bayes_init2.rds")
if (file.exists(IF)) { ans0 <- readRDS(IF); cat("init: reusing oz_bayes_init2.rds\n")
} else {
  cat("\n=== MoM init on repaired group (nsub=2, n3=1000) ===\n"); flush.console()
  t0 <- proc.time()[3]
  a0 <- sienaAlgorithmCreate(projname = "oz_init2", seed = 2012, nsub = 2, n3 = 1000, cond = FALSE)
  ans0 <- try(siena07(a0, data = s$grp, effects = s$eff, batch = TRUE,
                      useCluster = TRUE, nbrNodes = 8), silent = TRUE)
  if (inherits(ans0, "try-error")) { cat("init FAILED\n"); print(ans0); ans0 <- NULL
  } else { cat("init done", round((proc.time()[3]-t0)/60,1), "min | tconv.max",
               round(ans0$tconv.max,3), "\n"); saveRDS(ans0, IF) }
  flush.console() }

CHUNK <- 300; NSTAGE <- 5; NWARM1 <- 100
splitRhat <- function(x){n<-length(x);h<-floor(n/2);if(h<10)return(NA)
  a<-x[1:h];b<-x[(n-h+1):n];W<-mean(c(var(a),var(b)));if(W<=0)return(NA)
  gm<-mean(c(mean(a),mean(b)));B<-h*sum((c(mean(a),mean(b))-gm)^2)
  sqrt((((h-1)/h)*W+B/h)/W)}
rk<-function(x)qnorm((rank(x)-3/8)/(length(x)+1/4))
ess<-function(x){n<-length(x);r<-acf(x,lag.max=min(500,n-2),plot=FALSE)$acf[-1]
  k<-which(r<0.05)[1];if(is.na(k))k<-length(r);max(1,n/(1+2*sum(r[1:max(1,k-1)])))}
cumreport <- function(upto) {
  MU<-NULL;ET<-NULL;used<-c()
  for (k in 1:upto) { f<-file.path(SC,sprintf("ozbayes_stage%d.rds",k)); if(!file.exists(f))next
    g<-readRDS(f); ok<-!is.na(g$ThinPosteriorMu[,1]); if(k==1)ok[1:NWARM1]<-FALSE
    MU<-rbind(MU,g$ThinPosteriorMu[ok,,drop=FALSE]); ET<-rbind(ET,g$ThinPosteriorEta[ok,,drop=FALSE])
    used<-c(used,sprintf("s%d:%d",k,sum(ok))) }
  cat("\n######## OZ CUMULATIVE through stage",upto,"########\n  stacked:",paste(used,collapse=" "),"\n")
  out<-list()
  for (b in c("mu","eta")) {
    M<-if(b=="mu")MU else ET; labs<-if(b=="mu")labMu else labEta
    R<-apply(M,2,function(c)splitRhat(rk(c))); S<-apply(M,2,ess)
    cat(sprintf("  %-3s n=%4d | Rhat med %.4f max %.4f | <1.01: %2d/%2d | ESS med %4.0f min %4.0f\n",
        b,nrow(M),median(R,na.rm=TRUE),max(R,na.rm=TRUE),sum(R<1.01,na.rm=TRUE),length(R),median(S),min(S)))
    w<-order(-R)[1:3]
    cat("      worst:",paste(sprintf("%s %.3f (ESS %.0f)",labs[w],R[w],S[w]),collapse=" | "),"\n")
    out[[b]]<-data.frame(block=b,par=labs,Rhat=round(R,4),ESS=round(S)) }
  write.csv(do.call(rbind,out),file.path(SC,sprintf("ozbayes_diag_stage%d.csv",upto)),row.names=FALSE)
  flush.console() }

prev <- NULL
cat("\nchain start:", format(Sys.time()), "| budget", NSTAGE*CHUNK, "iterations\n"); flush.console()
for (k in 1:NSTAGE) {
  cat("\n===== OZ stage", k, "=====", format(Sys.time()), "\n"); flush.console()
  t0 <- proc.time()[3]
  algo <- sienaAlgorithmCreate(projname = paste0("oz_b", k), seed = 2012)
  args <- list(algo, data = s$grp, effects = s$eff,
               priorMu = s$priorMu, priorSigma = s$priorSig, priorKappa = 1.0,
               priorSigEta = rep(1, p2), nmain = CHUNK, nrunMHBatches = 20,
               nSampVarying = 3, nSampConst = 10, saveFreq = 20,
               nbrNodes = 8, clusterType = "PSOCK", silentstart = TRUE)
  if (k == 1) { args$nprewarm <- 20; args$nwarm <- NWARM1; args$prevAns <- ans0 }
  else        { args$prevBayes <- prev }
  g <- try(do.call(sienaBayes, args), silent = FALSE)
  el <- proc.time()[3] - t0
  if (inherits(g, "try-error")) { cat("OZ stage", k, "FAILED\n"); print(g); break }
  saveRDS(g, file.path(CHN, sprintf("ozbayes_stage%d.rds", k)))
  cat(sprintf("  OZ stage %d done: %d iterations, %.2f h wall, %.1f s/iteration\n",
              k, CHUNK, el/3600, el/CHUNK)); flush.console()
  cumreport(k); prev <- g
}
cat("\nOZ FINISHED", format(Sys.time()), "\n")
