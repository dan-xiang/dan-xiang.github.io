#R code for fitting
normsq <- function(v) sum(v^2)
norm <- function(v) sqrt(sum(v^2))
zeta0 <- function(y, alpha=1, df1=1, df2=Inf, L){
  # zeta_{s, alpha}(y) for scalar y and inverse power measure with index 0 < alpha < 2
  # scalar version only
  if(missing(L)){if(df2 > 20) L <- 10 + 3*y^2/2 else L <- 10 + (18+log(2*df2/df1)) * y^2/df2}
  if(df2==Inf && abs(y) > 37.6) return(.Machine$double.xmax);
  r <- 1:ceiling(L);
  if((ysq <- normsq(y)) == 0) return(0.0);
  coef1 <- log(abs(2*r-2-alpha)) - log(2*r) + log(ysq) - log(df1+2*r-2);
  if(df2 == Inf) return(sum(exp(cumsum(coef1))));
  coef2  <- log(df1+df2+2*r-2) - log(df2 + ysq);
  return(sum(exp(cumsum(coef1 + coef2))));
}
zeta <- function(y, alpha=1, df1=1, df2=Inf){
  # vector version of zeta0
  for(i in 1:length(y)) y[i] <- zeta0(y[i], alpha, df1, df2)
  return(y)
}

mixfit0 <- function(zeta, rho=0.05, maxcycle=10, tol=1.0e-5){
  # fits mixture model with given zeta-vector
  zeta <- zeta - 1
  if(sum(zeta) <= 0) return(list(rho=0.0, cycle=0, llik=0, activity=rep(0, length(zeta))))
  for(cycle in 1:maxcycle){ 
    z <- zeta / (1 + rho * zeta)      # vector of derivatives
    delta <- sum(z) / sum(z^2) / rho
    rho <- min(rho * exp(delta), 2*rho)
    if(abs(delta) < tol) break
  }
  llik <- sum(log(1 + rho*zeta))
  list(rho=rho, cycle=cycle, llik=llik, activity=rho * (zeta+1)/(1 + rho*zeta))
}

mixfit <- function(y, alpha, rho=0.05, df2=Inf, maxcycle=10, tol=1.0e-3, div0=1.09, div1=2.5){
  # fits mixture model with general inverse-power measure
  trace <- matrix(0, maxcycle, 4);  colnames(trace) <- c("alpha", "llik", "dalpha", "rho")
  if(!missing(alpha)) return(mixfit0(zeta(y, alpha, df2=df2), rho=rho))
  alpha <- 1;  dalpha <- 0.2;
  fit0 <- mixfit0(zeta(y, alpha, df2=df2), rho=rho, maxcycle=3, tol=1.0e-2)
  trace[1,] <- c(alpha, fit0$llik, dalpha, fit0$rho)
  for(cycle in 2:maxcycle){
    fit1 <- mixfit0(zeta(y, alpha+dalpha, df2=df2), rho=fit0$rho, maxcycle=3)
    trace[cycle,] <- c(alpha+dalpha, fit1$llik, dalpha, fit1$rho)
    if(fit1$llik > fit0$llik){
      alpha <- alpha + dalpha
      fit0 <- fit1
      dalpha <- dalpha / div0
    } else {if(cycle == 2) dalpha <- -dalpha / div0 else dalpha <- -dalpha / div1 }
    if(abs(dalpha) < tol) break()
  }
  alpha <- trace[which.max(trace[1:cycle, 2]), 1]
  if(cycle > 4){
    coef <- lm(llik~alpha+I(alpha^2), data = data.frame(trace)[(cycle-3):cycle,])$coef
    alpha <- -coef[2]/(2*coef[3])
  }
  fit1 <- mixfit0(zeta(y, alpha, df2=df2), rho=fit1$rho)
  list(rho=fit1$rho, cycle=cycle, llik=fit1$llik, activity=fit1$activity, alpha=alpha, trace=trace[1:cycle,])
}

# Example computation with scores on 6 d.f. for hivdata z scores

hivdat <- log(hivdata)
treated <- rowMeans(hivdat[,5:8])
control <- rowMeans(hivdat[,1:4])
sd_treat <- apply(hivdat[,5:8],1,sd)
sd_control <- apply(hivdat[,1:4],1,sd)
pool_var <- (sd_treat^2*3+sd_control^2*3)/6*(1/4+1/4)
t_scores <- (treated-control)/sqrt(pool_var)
z_scores <- qnorm(pt(t_scores,df=6))
hist(z_scores,prob=F,breaks=70)

# computations for tables 2 and 3

# fit d
inds <- order(abs(z_scores))
fit0 <- mixfit(z_scores, df2=Inf)
lfdr0 <- 1-fit0$activity
fit1 <- mixfit(t_scores, df2=6)
lfdr1 <- 1-fit1$activity
rev(tail(z_scores[inds])) # print top 6 z scores
rev(tail(t_scores[inds])) # print top 6 t-ratios
# produce values for table 2
as.numeric(formatC(rev(tail(lfdr0[inds])), format = "e", digits = 3))*10^4
as.numeric(formatC(rev(tail(lfdr1[inds])), format = "e", digits = 3))*10^4
formatC(rev(tail(lfdr0[inds])/tail(lfdr1[inds])), format = "e", digits = 2)


# repeat with d=1 fixed
inds <- order(abs(z_scores))
fit0 <- mixfit(z_scores, alpha=1, df2=Inf)
lfdr0 <- 1-fit0$activity
fit1 <- mixfit(t_scores, alpha=1, df2=6)
lfdr1 <- 1-fit1$activity
rev(tail(z_scores[inds]))
rev(tail(t_scores[inds]))
as.numeric(formatC(rev(tail(lfdr0[inds])), format = "e", digits = 3))*10^4
as.numeric(formatC(rev(tail(lfdr1[inds])), format = "e", digits = 3))*10^4
rev(tail(lfdr0[inds])/tail(lfdr1[inds]))



# compare with BH set
p_values <- 2*pnorm(abs(z_scores),lower=F)
q <- 0.1
R_bh <- sum(p.adjust(p_values,method="BH")<q)
BH_set <- which(p.adjust(p_values,method="BH")<q)
mean(lfdr0[BH_set])
mean(lfdr1[BH_set])









