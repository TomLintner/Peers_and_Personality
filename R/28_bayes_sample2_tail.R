## 28_bayes_sample2_tail.R
## source: Results/scripts/oz_tail.R
source("00_setup.R")

## Is a TAIL of the OZ chain internally clean, even though the whole chain drifts?
## For each candidate burn-in fraction, recompute on the retained portion only:
##   - split-Rhat (rank-normalised) and ESS
##   - half-to-half drift in posterior-SD units
## If a tail is stationary, its drift counts collapse and its Rhat clears.
suppressMessages(library(RSienaTest))
s <- readRDS(file.path(OBJ, "oz_bayes_setup.rds")); E <- s$eff; class(E) <- "data.frame"
selR <- E$include & ((!E$basicRate & E$randomEffects) | (E$basicRate & E$group == 1))
fxi  <- E$include & !E$basicRate & !E$randomEffects & E$type == "eval"
labMu <- E$effectName[selR]; labEta <- E$effectName[fxi]
MU <- NULL; ET <- NULL
for (k in 1:20) {
  f <- file.path(CHN, sprintf("ozbayes_stage%d.rds", k)); if (!file.exists(f)) next
  g <- readRDS(f); ok <- !is.na(g$ThinPosteriorMu[,1]); if (k==1) ok[1:100] <- FALSE
  MU <- rbind(MU, g$ThinPosteriorMu[ok,,drop=FALSE]); ET <- rbind(ET, g$ThinPosteriorEta[ok,,drop=FALSE])
}
N <- nrow(MU); cat("full chain:", N, "draws\n")
splitRhat <- function(x){n<-length(x);h<-floor(n/2);if(h<10)return(NA)
  a<-x[1:h];b<-x[(n-h+1):n];W<-mean(c(var(a),var(b)));if(W<=0)return(NA)
  gm<-mean(c(mean(a),mean(b)));B<-h*sum((c(mean(a),mean(b))-gm)^2);sqrt((((h-1)/h)*W+B/h)/W)}
rk<-function(x)qnorm((rank(x)-3/8)/(length(x)+1/4))
ess<-function(x){n<-length(x);r<-acf(x,lag.max=min(500,n-2),plot=FALSE)$acf[-1]
  k<-which(r<0.05)[1];if(is.na(k))k<-length(r);max(1,n/(1+2*sum(r[1:max(1,k-1)])))}
drift <- function(M){n<-nrow(M);h<-floor(n/2)
  apply(M,2,function(x) (mean(x[(n-h+1):n])-mean(x[1:h]))/sd(x))}

cat(sprintf("\n%-8s %6s | %-28s | %-28s\n", "keep", "draws",
            "MU  Rhat max /<1.01/ ESSmed", "ETA Rhat max /<1.01/ ESSmed"))
cat(strrep("-", 78), "\n")
res <- list()
BURN <- 1200
for (frac in c(1.00, 0.75, (N-BURN)/N, 0.50, 0.40, 0.25)) {
  i0 <- N - floor(N*frac) + 1
  A <- MU[i0:N,,drop=FALSE]; B <- ET[i0:N,,drop=FALSE]
  RA<-apply(A,2,function(c)splitRhat(rk(c))); RB<-apply(B,2,function(c)splitRhat(rk(c)))
  SA<-apply(A,2,ess); SB<-apply(B,2,ess); DA<-drift(A); DB<-drift(B)
  cat(sprintf("last %3.0f%% %6d | %.3f  %2d/19  %5.0f          | %.3f  %2d/30  %5.0f\n",
      100*frac, nrow(A), max(RA,na.rm=TRUE), sum(RA<1.01,na.rm=TRUE), median(SA),
      max(RB,na.rm=TRUE), sum(RB<1.01,na.rm=TRUE), median(SB)))
  cat(sprintf("           drift >0.25SD: mu %2d/19  eta %2d/30 | >0.10SD: mu %2d/19  eta %2d/30\n",
      sum(abs(DA)>.25), sum(abs(DB)>.25), sum(abs(DA)>.10), sum(abs(DB)>.10)))
  res[[as.character(frac)]] <- list(RB=RB,SB=SB,DB=DB,n=nrow(A))
}
## worst eta parameters in the last third
r <- res[[as.character(1/3)]]
o <- order(-r$RB)[1:8]
cat("\n=== last third: worst eta parameters ===\n")
print(data.frame(par=labEta[o], Rhat=round(r$RB[o],3), ESS=round(r$SB[o]),
                 drift_SD=round(r$DB[o],3)), row.names=FALSE)
