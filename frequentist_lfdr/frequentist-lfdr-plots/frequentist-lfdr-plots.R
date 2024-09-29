###### R code for reproducing plots in frequentist lfdr manuscript

# load packages
install.packages("ggplot2")
library("ggplot2"); theme_set(theme_classic(base_size = 17))
install.packages("latex2exp")
library(latex2exp)

# load data
dat <- read.csv("mhhb_nma_data_corrected.csv")
dat$tstat <- with(dat, cohens_d / sqrt(variance_d))
dat$pvalue <- with(dat, pt(-tstat, df=n_comparison-2))

trunc.pvals <- with(dat, 40*pvalue[pvalue < 0.025 & tstat > 0])
pi0.hat <- 2*mean(trunc.pvals > 0.5)
sum(p.adjust(trunc.pvals, method="BH") < 0.1 / pi0.hat)
R_bh <- sum(round(p.adjust(trunc.pvals, method="BH"), 2)<0.1/pi0.hat)
bh_thresh <- sort(trunc.pvals)[R_bh]
mbreaks <- c(seq(0,bh_thresh,length=9),seq(0.316,1,length=20)) # so that no bins are empty

df <- data.frame(pvals=trunc.pvals)
df$unscaled_pvals <- dat$pvalue[dat$pvalue<0.025]
df$group = ifelse(df$pvals < 0.5, 'A', 'B')
df$Fm <- c(1:length(trunc.pvals))
df$ord_pvals <- sort(df$pvals)

########### Figure 1: P-value histogram with bFDR estimate

df$group = ifelse((df$pvals < 0.203) | (df$pvals > bh_thresh), 'A', 'B')
pi0prop <- pi0.hat # put pi0hat at the correct relative height on the count hist
count.pi0.hat <- 9*pi0prop
ggplot(df, aes(x=pvals,fill=as.factor(group))) + 
  geom_histogram(breaks=mbreaks) + 
  scale_fill_manual(values = c('A' = 'gray','B' = 'red'))+
  labs(title="         Histogram with BH(q) threshold, q = 0.1",x="right-tailed p-value (selection-adjusted)", y = "count") +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())+guides(fill = "none")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),expand = c(0,0))+
  scale_y_continuous(expand = c(0,-7))+
  guides(fill = "none")+
  geom_vline(xintercept=sort(df$pvals)[R_bh],lwd=0.5,colour="red")+
  annotate("text", x=bh_thresh, y=-7,
           label=TeX(r'($\tau_{q}^{BH}$)'), parse=TRUE, col="red",size=4.5)+
  annotate("text", x=0.5, y=30, 
           label=TeX(r'($\widehat{FDP}(\[0.20,0.27\])=0.32$)'),
           parse=TRUE, col="red",size=6)+
  coord_cartesian(xlim=c(0,1),clip = "off")


########### Figure 2: Broken line plot

k <- 115
ggplot(df, aes(x=Fm,y=sort(pvals))) + geom_point() + 
  labs(x="k",y=TeX(r'($p_{(k)}$)'),title="            Secant line FDP estimator")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_text(hjust=0.45),
        axis.title.y = element_text(hjust=0.6))+guides(fill = "none")+
  scale_x_continuous(breaks = seq(0, nrow(df), by = 60),expand=c(0,0))+
  coord_cartesian(ylim=c(0,0.32),xlim=c(0,230),clip = "off")+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),expand=c(0,0))+
  geom_segment(aes(x=0,y=0,
                   xend=k,yend=df$ord_pvals[k]),color="red")+
  geom_segment(aes(x=k,y=df$ord_pvals[k],
                   xend=R_bh,yend=df$ord_pvals[R_bh]),color="red")+
  geom_segment(aes(x=R_bh,y=df$ord_pvals[R_bh],xend=R_bh,yend=0),color="red",linetype="dashed")+
  geom_segment(aes(x=0,y=df$ord_pvals[R_bh],xend=R_bh,yend=df$ord_pvals[R_bh]),color="red",linetype="dashed")+
  annotate("text", x=90, y=0.22, 
           label=TeX(r'($\widehat{FDP}(\[0.01,0.27\])= \hat{\pi}_0 m \times \frac{0.26}{87}=0.22$)'),
           parse=TRUE,color="red",size=5)+
  annotate("text", x=128, y=0.14, 
           label=TeX(r'($Slope = \frac{0.26}{87}$)'),
           parse=TRUE,color="red",size=5)+
  annotate("text", x=50, y=0.035,
           label=TeX(r'($Slope = \frac{0.01}{115}$)'),
           parse=TRUE,color="red",size=5)+
  annotate("text", x=-9, y=df$ord_pvals[R_bh],
           label=TeX(r'($0.27$)'), parse=TRUE,color="red",size=5)+
  annotate("text", x=-9, y=df$ord_pvals[k],
           label=TeX(r'($0.01$)'), parse=TRUE,color="red",size=5)+
  annotate("text", x=k, y=-0.015,
           label=TeX(r'($115$)'), parse=TRUE,color="red",size=5)+
  annotate("text", x=R_bh, y=-.015, 
           label=TeX(r'($202$)'), parse=TRUE,color="red",size=5)



########### Figure 3: Scatterplot of clfdr versus mlfdr with m=1000 tests

install.packages("ramify")
library(ramify)

# sinkhorn algorithm from https://rdrr.io/github/aaamini/nett/src/R/sinkhorn.R

l2_norm = function(x) sqrt(sum(x^2))

sinkhorn_knopp = function(A, sums = rep(1, nrow(A)),
                          niter = 100, tol = 1e-8, sym = FALSE, verb = FALSE) {
  
  delta = Inf
  r = c = rep(1, nrow(A))
  converged = FALSE
  t = 1
  while( t <= niter && !converged) {
    r = sums / (A %*% c)
    cnew = sums / (t(A) %*% r)
    
    if (sym) cnew = (cnew + r)/2
    
    delta = l2_norm(cnew-c)
    if (delta < tol) converged = TRUE
    if (verb) nett::printf("err = %3.5e\n", delta)
    c = cnew
    t = t+1
  }
  list(r = as.numeric(r), c = as.numeric(c))
}

# sinkhorn approximation of matrix permanent (returns log permanent)
sinkhorn_perm <- function(B) {
  n <- dim(B)[1]
  h <- rep(1,n)
  I <- diag(1,n)
  J <- 1/(n)*t(t(rep(1,n))) %*% t(rep(1,n))
  t <- n/(n-1)
  
  out = sinkhorn_knopp(B,sums=h)
  r <- out$r
  c <- out$c
  A <- diag(r) %*% B %*% diag(c)
  ata <- t(A) %*% A
  decomp <- eigen(I+t^2*J-t^2*ata)
  approx <- sum(log(c(1:n)))-1/2*sum(log(abs(decomp$values)))
  return(approx-n*log(n)-sum(log(r*c)))
}

rds_file <- "scatterplot_data.rds"
if (file.exists(rds_file)) {scatterplot_data <- readRDS(rds_file)} else {
set.seed(12)
m <- 1000
a <- 1/4
b <- 1
m0 <- 800
m1 <- m-m0
pi0 <- m0/m
H <- c(0,1,rep(0,m0-1),rep(1,m1-1))
pvals <- rep(0,m)
K <- 6 # number of realizations
comps <- matrix(NA,nrow=m,ncol=K)
simps <- matrix(NA,nrow=m,ncol=K)
for (i in 1:K) {
  print(i)
  pvals[as.logical(H)] <- rbeta(m1,a,b)
  pvals[as.logical(1-H)] <- runif(m0)
  ord_pvals <- sort(pvals)
  dens_mat <- matrix(1,nrow=m,ncol=m)
  for (j in 1:m1) {
    dens_mat[as.logical(H),][j,] <- dbeta(pvals,a,b)
  }
  comp_lfdr <- rep(0,m)
  simp_lfdr <- pi0/(pi0+(1-pi0)*dbeta(ord_pvals,a,b))
  for (k in 1:m) {
    print(k)
    ind <- order(pvals)[k]
    X_ratio <- exp(sinkhorn_perm(dens_mat[-2,-ind])-sinkhorn_perm(dens_mat[-1,-ind]))
    comp_lfdr[k] <- pi0/(pi0+(1-pi0)*dbeta(ord_pvals[k],a,b)*X_ratio)
  }
  comps[,i] <- comp_lfdr
  simps[,i] <- simp_lfdr
}

tmp <- matrix(NA,nrow=m,ncol=K)
for (i in 1:K) {
  tmp[,i] <- i
}
simple_lfdr <- flatten(simps,across="columns")
compound_lfdr <- flatten(comps,across="columns")
realization <- flatten(tmp,across="columns")
df <- data.frame(simple_lfdr = simple_lfdr,
                 compound_lfdr = compound_lfdr,
                 realization = realization)
df$realization <- as.factor(df$realization)

saveRDS(df, "scatterplot_data.rds") # save the simulation data (takes a day to produce)

}

df <- scatterplot_data

ggplot(df, aes(x = simple_lfdr, y = compound_lfdr, color = realization)) +
  geom_point(size = .25) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(x = "marginal lfdr", y = "compound lfdr",title="   Compound versus marginal lfdr, m=1000")+
  xlim(0, 1)+
  ylim(0, 1)+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size=4)))


########### Figure 4: Small experiment to illustrate clfdr with m=6 tests

library(combinat)

mat_perm <- function(M) {
  n <- dim(M)[1]
  x <- permn(c(1:n))
  perm_sum <- 0
  for (i in 1:length(x)) {
    sig <- x[[i]]
    prod <- 1
    for (j in 1:n) {
      prod <- prod*M[j,sig[j]]
    }
    perm_sum <- perm_sum+prod
  }
  return(perm_sum)
}

set.seed(123)
m <- 6
m0 <- 4
m1 <- m-m0
pi0 <- m0/m
a=1/4
b=1
H <- c(0,1,rep(0,m0-1),rep(1,m1-1))
clfdr_fun <- function(p,a,b) {
  comp_lfdr <- rep(0,m) 
  dens_mat <- matrix(1,nrow=m,ncol=m)
  for (j in 1:m1) {
    dens_mat[as.logical(H),][j,] <- dbeta(p,a,b)
  }
  for (k in 1:m) {
    print(k)
    X_ratio <- mat_perm(dens_mat[-2,-k])/mat_perm(dens_mat[-1,-k])
    comp_lfdr[k] <- pi0/(pi0+(1-pi0)*dbeta(p[k],a,b)*X_ratio)
  }
  return(comp_lfdr)
}
pvals1 <- c(runif(1),rbeta(1,a,b),runif(m0-1),rbeta(m1-1,a,b))
clfdr1 <- clfdr_fun(pvals1,a,b)
pvals2 <- c(runif(1),rbeta(1,a,b),runif(m0-1),rbeta(m1-1,a,b))
clfdr2 <- clfdr_fun(pvals2,a,b)
pvals3 <- c(runif(1),rbeta(1,a,b),runif(m0-1),rbeta(m1-1,a,b))
clfdr3 <- clfdr_fun(pvals3,a,b)

df1 <- data.frame(x = pvals1, y = clfdr1, realization = "1")
df2 <- data.frame(x = pvals2, y = clfdr2, realization = "2")
df3 <- data.frame(x = pvals3, y = clfdr3, realization = "3")

# Combine data frames
combined_df <- rbind(df1, df2, df3)

t <- seq(0,1,length=100)
mlfdr <- pi0/(pi0+(1-pi0)*dbeta(t,a,b))
df_mlfdr <- data.frame(x=t,y=mlfdr)

# Plot clfdr and mlfdr
ggplot() +
  geom_point(data=combined_df, aes(x = x, y = y, color = realization)) +
  geom_line(data=combined_df, aes(x = x, y = y, color = realization)) +
  geom_line(data=df_mlfdr, aes(x=x,y=y,color="mlfdr")) +
  theme_minimal()+
  scale_y_continuous(limits = c(0, 1))+
  labs(y="local fdr",
       x="p-value",
       title="    Three realizations of compound lfdr, m=6") +
  theme(legend.text = element_text(size = 18),
        legend.position = c(0.80, 0.35),
        legend.title = element_text(size=18),
        axis.text = element_text(size=18),
        axis.title = element_text(size=18),
        plot.title = element_text(size=18))


########### Figure 5: Larger experiment to show clfdr and mlfdr \approx local FDP

rds_file <- "calibration_plot_data.rds"
if (file.exists(rds_file)) {dat <- readRDS(rds_file)} else {
set.seed(12)
a <- 1/8
b <- 1
m0 <- 320
m1 <- 80
m <- m0+m1
pi0 <- m0/m
H <- c(0,1,rep(0,m0-1),rep(1,m1-1))
pvals <- rep(0,m)
N <- 10^3 # 10^4
comps <- rep(0,N)
simps <- rep(0,N)

pvals[as.logical(H)] <- rbeta(m1,a,b)
pvals[as.logical(1-H)] <- runif(m0)
ord_pvals <- sort(pvals)
dens_mat <- matrix(1,nrow=m,ncol=m)
for (j in 1:m1) {
  dens_mat[as.logical(H),][j,] <- dbeta(pvals,a,b)
}
comp_lfdr <- rep(0,m)
simp_lfdr <- pi0/(pi0+(1-pi0)*dbeta(pvals,a,b))
for (k in 1:m) {
  print(k)
  X_ratio <- exp(sinkhorn_perm(dens_mat[-2,-k])-sinkhorn_perm(dens_mat[-1,-k]))
  comp_lfdr[k] <- pi0/(pi0+(1-pi0)*dbeta(pvals[k],a,b)*X_ratio)
}
lfdp <- rep(NA,m)
for (i in 16:(m-15)) {
  grid_n <- 31
  inds <- order(pvals)
  left <- i - 15
  right <- i + 15
  lfdp[i] <- sum(1-H[inds][left:right])/grid_n
}
dat <- data.frame(lfdp = lfdp[rank(pvals)],
                  clfdr = comp_lfdr,
                  mlfdr = simp_lfdr,
                  pvals=pvals)
saveRDS(dat, "calibration_plot_data.rds") # save the simulation data (takes a day to produce)
}

ggplot(data = dat, aes(x = pvals)) +
  geom_point(aes(x = pvals, y = clfdr, color = "clfdr"),size=1) +
  geom_point(aes(x = pvals, y = lfdp, color = "local FDP"),size=1) + 
  geom_line(aes(x = pvals, y = mlfdr, color = "mlfdr"),size=.75) +
  labs(x = "p-values", y = "l-values",
       title="Comparison of clfdr, mlfdr, and local FDP") +
  scale_color_manual(values = c("clfdr" = "red", "local FDP" = "blue", "mlfdr" = "black")) +
  scale_shape_manual(values = c("clfdr" = 2, "lfdp" = 3)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.80, 0.35),
        legend.text = element_text(size = 20))


########### Figure 6: Hockey stick plot


set.seed(12)
m <- 500
a <- 0.05
b <- 1
m0 <- 250
m1 <- m-m0
pi0 <- m0/m
H <- c(0,1,rep(0,m0-1),rep(1,m1-1))
pvals <- rep(0,m)
pvals[as.logical(H)] <- rbeta(m1,a,b)
pvals[as.logical(1-H)] <- runif(m0)

df <- data.frame(ord_pvals=sort(pvals),rank=c(1:m),H=as.factor(H[order(pvals)]))
ggplot(df,aes(x=rank,y=ord_pvals,color=H))+
  geom_point()+
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
  labs(x="Rank k",y=TeX(r'(Order statistic $p_{(k)}$)'),
       title="   Comparison of bFDR and FDR")+
  scale_color_manual(values = c("1" = "blue", "0" = "red"))+
  theme_minimal()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "none")+
    annotate("text", x=120, y=0.75,
             label=TeX(r'(FDR(0,0.1) \approx 9\%)'), parse=TRUE,size=6)+
    annotate("text", x=120, y=0.6,
             label=TeX(r'(bFDR(0,0.1) \approx 60\%)'), parse=TRUE,size=6)

###### compute boundary FDR
# df$num_H <- as.numeric(df$H)-1
# sum(df$ord_pvals<0.1)
# sum(1-df$num_H[230:244])/(244-230+1)
###### compute FDR
#sum(df$ord_pvals<0.1)
#sum(1-df$num_H[1:244])/244

