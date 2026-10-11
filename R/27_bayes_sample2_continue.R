## 27_bayes_sample2_continue.R
## source: Results/scripts/oz_cont.R
source("00_setup.R")

## OZ continuation: resume from the highest existing stage and keep adding
## 300-iteration stages until EVERY parameter (mu and eta) has split-Rhat
## below 1.01, or stage KMAX is reached.  Settings identical to stages 1-5.
suppressMessages(library(RSienaTest))
set.seed(2012)
SC <- WORK
WD <- file.path(SC, "ozbayes"); dir.create(WD, showWarnings = FALSE); setwd(WD)
s <- readRDS(file.path(OBJ, "oz_bayes_setup.rds"))
E <- s$eff; class(E) <- "data.frame"
p2   <- sum(E$include & !E$basicRate & !E$randomEffects & E$type == "eval")
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
labMu  <- paste(E$name[selR], E$shortName[selR])
fxi <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labEta <- paste(E$name[fxi], E$shortName[fxi])

CHUNK <- 300; KMAX <- 11; NWARM1 <- 100; TARGET <- 1.01
splitRhat <- function(x){n<-length(x);h<-floor(n/2);if(h<10)return(NA)
  a<-x[1:h];b<-x[(n-h+1):n];W<-mean(c(var(a),var(b)));if(W<=0)return(NA)
  gm<-mean(c(mean(a),mean(b)));B<-h*sum((c(mean(a),mean(b))-gm)^2);sqrt((((h-1)/h)*W+B/h)/W)}
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
  out<-list(); worst<-0
  for (b in c("mu","eta")) {
    M<-if(b=="mu")MU else ET; labs<-if(b=="mu")labMu else labEta
    R<-apply(M,2,function(c)splitRhat(rk(c))); S<-apply(M,2,ess)
    worst <- max(worst, max(R, na.rm=TRUE))
    cat(sprintf("  %-3s n=%4d | Rhat med %.4f max %.4f | <1.01: %2d/%2d | ESS med %4.0f min %4.0f\n",
        b,nrow(M),median(R,na.rm=TRUE),max(R,na.rm=TRUE),sum(R<TARGET,na.rm=TRUE),length(R),median(S),min(S)))
    ab <- which(R >= TARGET)
    if (length(ab)) cat("      still >= 1.01:",
        paste(sprintf("%s %.3f(ESS %.0f)", labs[ab], R[ab], S[ab]), collapse=" | "), "\n")
    out[[b]]<-data.frame(block=b,par=labs,Rhat=round(R,4),ESS=round(S)) }
  write.csv(do.call(rbind,out),file.path(SC,sprintf("ozbayes_diag_stage%d.csv",upto)),row.names=FALSE)
  flush.console(); worst }

done <- which(sapply(1:KMAX, function(k) file.exists(file.path(CHN, sprintf("ozbayes_stage%d.rds", k)))))
k0 <- max(done)
cat("resuming after stage", k0, "| target: all Rhat <", TARGET, "| cap: stage", KMAX, "\n")
w <- cumreport(k0)
if (w < TARGET) { cat("\nAlready converged at stage", k0, "- nothing to do.\n"); quit("no") }
prev <- readRDS(file.path(SC, sprintf("ozbayes_stage%d.rds", k0)))
cat("start:", format(Sys.time()), "\n"); flush.console()
for (k in (k0+1):KMAX) {
  cat("\n===== OZ stage", k, "=====", format(Sys.time()), "\n"); flush.console()
  t0 <- proc.time()[3]
  algo <- sienaAlgorithmCreate(projname = paste0("oz_b", k), seed = 2012)
  g <- try(sienaBayes(algo, data = s$grp, effects = s$eff,
                      priorMu = s$priorMu, priorSigma = s$priorSig, priorKappa = 1.0,
                      priorSigEta = rep(1, p2), nmain = CHUNK, nrunMHBatches = 20,
                      nSampVarying = 3, nSampConst = 10, saveFreq = 20,
                      prevBayes = prev, nbrNodes = 8, clusterType = "PSOCK",
                      silentstart = TRUE), silent = FALSE)
  el <- proc.time()[3] - t0
  if (inherits(g, "try-error")) { cat("OZ stage", k, "FAILED\n"); print(g); break }
  saveRDS(g, file.path(CHN, sprintf("ozbayes_stage%d.rds", k)))
  cat(sprintf("  OZ stage %d done: %d iterations, %.2f h wall, %.1f s/iteration\n",
              k, CHUNK, el/3600, el/CHUNK)); flush.console()
  w <- cumreport(k); prev <- g
  if (w < TARGET) { cat("\n*** ALL Rhat <", TARGET, "at stage", k, "- stopping. ***\n"); break }
  if (k == KMAX) cat("\n*** reached stage cap", KMAX, "with max Rhat", round(w,4), "***\n")
}
cat("\nOZ CONTINUATION FINISHED", format(Sys.time()), "\n")
