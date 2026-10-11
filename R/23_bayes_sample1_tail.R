## 23_bayes_sample1_tail.R
## source: Results/scripts/gacr_tail.R
source("00_setup.R")

## Same drift + tail diagnostics that exposed OZ's 1200-draw burn-in,
## applied to the GACR chain (9 stages, 2650 post-warm draws).
suppressMessages(library(RSienaTest))
s <- readRDS(file.path(OBJ, "gacr_bayes_setup2.rds")); E <- s$eff; class(E) <- "data.frame"
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
fxi  <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labMu <- E$effectName[selR]; labEta <- E$effectName[fxi]
MU <- NULL; ET <- NULL
for (k in 1:9) {
  f <- file.path(CHN, sprintf("bayes12h_stage%d.rds", k)); if (!file.exists(f)) next
  g <- readRDS(f); ok <- !is.na(g$ThinPosteriorMu[,1]); if (k==1) ok[1:100] <- FALSE
  MU <- rbind(MU, g$ThinPosteriorMu[ok,,drop=FALSE]); ET <- rbind(ET, g$ThinPosteriorEta[ok,,drop=FALSE])
}
N <- nrow(MU); cat("GACR full chain:", N, "draws\n")
splitRhat <- function(x){n<-length(x);h<-floor(n/2);if(h<10)return(NA)
  a<-x[1:h];b<-x[(n-h+1):n];W<-mean(c(var(a),var(b)));if(W<=0)return(NA)
  gm<-mean(c(mean(a),mean(b)));B<-h*sum((c(mean(a),mean(b))-gm)^2);sqrt((((h-1)/h)*W+B/h)/W)}
rk<-function(x)qnorm((rank(x)-3/8)/(length(x)+1/4))
ess<-function(x){n<-length(x);r<-acf(x,lag.max=min(500,n-2),plot=FALSE)$acf[-1]
  k<-which(r<0.05)[1];if(is.na(k))k<-length(r);max(1,n/(1+2*sum(r[1:max(1,k-1)])))}
drift<-function(M){n<-nrow(M);h<-floor(n/2)
  apply(M,2,function(x)(mean(x[(n-h+1):n])-mean(x[1:h]))/sd(x))}
cat(sprintf("\n%-10s %6s | %-26s | %-26s\n","keep","draws","MU Rhat max/<1.01/ESSmed","ETA Rhat max/<1.01/ESSmed"))
cat(strrep("-",76),"\n")
for (bn in c(0, 300, 600, 900, 1325)) {
  A <- MU[(bn+1):N,,drop=FALSE]; B <- ET[(bn+1):N,,drop=FALSE]
  RA<-apply(A,2,function(c)splitRhat(rk(c))); RB<-apply(B,2,function(c)splitRhat(rk(c)))
  DA<-drift(A); DB<-drift(B)
  cat(sprintf("burn %5d %6d | %.3f  %2d/19  %5.0f       | %.3f  %2d/30  %5.0f\n",
      bn, nrow(A), max(RA,na.rm=TRUE), sum(RA<1.01,na.rm=TRUE), median(apply(A,2,ess)),
      max(RB,na.rm=TRUE), sum(RB<1.01,na.rm=TRUE), median(apply(B,2,ess))))
  cat(sprintf("             drift >0.25SD: mu %2d/19 eta %2d/30 | >0.10SD: mu %2d/19 eta %2d/30\n",
      sum(abs(DA)>.25), sum(abs(DB)>.25), sum(abs(DA)>.10), sum(abs(DB)>.10)))
}
D <- drift(ET); o <- order(-abs(D))[1:8]
cat("\n=== GACR eta: largest half-to-half shifts (full chain) ===\n")
print(data.frame(par=labEta[o], drift_SD=round(D[o],3)), row.names=FALSE)
