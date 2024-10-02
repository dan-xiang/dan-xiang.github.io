###### R code for reproducing plots in frequentist lfdr manuscript

# load packages
install.packages("ggplot2")
library("ggplot2"); theme_set(theme_classic(base_size = 17))
install.packages("latex2exp")
library(latex2exp)

library(devtools)
devtools::install_github("twbattaglia/MicrobeDS")
library(phyloseq)
library(MicrobeDS)
library(fdrtool)



############ Figure 1: Strimmer's lfdr estimate on microbiome data

################# stratified p-value calculation (9/27/2024)

rds_file <- "micro-biome-strat-pvals.rds"

if (file.exists(rds_file)) {perm.pvals <- readRDS(rds_file)} else {
  
  data('qa10394')
  # the sample data is available from the data source
  con <- file("10394_20180418-105453.txt","r")
  qa_whole_txt <- readLines(con)
  sample_data_qa10394 <- data.frame(do.call(rbind, strsplit(qa_whole_txt[2:length(qa_whole_txt)], '\t')))
  colnames(sample_data_qa10394) <- strsplit(qa_whole_txt[1], '\t')[[1]]
  close(con)
  
  
  names(sample_data_qa10394)
  
  OTU <- otu_table(qa10394)
  
  colSums(OTU)
  
  sample_data_qa10394$sample_name
  
  homo.OTU4 <- OTU[,grep("H4.rep.*fresh",colnames(OTU))]
  homo.OTU2 <- OTU[,grep("H2.rep.*fresh",colnames(OTU))]
  colSums(homo.OTU4)
  
  colnames(OTU)
  
  
  which.nonzero <- which(rowSums(homo.OTU4) > 0 & rowSums(homo.OTU2) > 0)
  which.lots <- which(rowSums(homo.OTU2 > 0) + rowSums(homo.OTU4 > 0) > 5)
  
  grp1 <- OTU[,grep("H[123789].95etoh.rep",colnames(OTU))]
  grp2 <- OTU[,grep("H[123789].rep.*fresh",colnames(OTU))]
  #grp1 <- OTU[,grep("H[123456789].95etoh.rep",colnames(OTU))]
  #grp2 <- OTU[,grep("H[123456789].rep.*fresh",colnames(OTU))]
  which.lots <- which(rowSums(grp1 > 0) + rowSums(grp2 > 0) > 10)
  
  
  ## stratified permutation test (from ChatGPT)
  stratified_permutation_test <- function(X, groups, strata, B) {
    n <- nrow(X)
    p <- ncol(X)
    unique_strata <- unique(strata)
    R <- matrix(NA, n, p)
    
    # Compute ranks within each stratum
    for (s in unique_strata) {
      idx <- which(strata == s)
      X_s <- X[idx, , drop = FALSE]
      ranks_s <- apply(X_s, 2, rank)
      R[idx, ] <- ranks_s
    }
    
    # Encode groups as +1 and -1
    groups <- as.factor(groups)
    group_levels <- levels(groups)
    if (length(group_levels) != 2) {
      stop("The 'groups' variable must have exactly two levels.")
    }
    groups_coded <- ifelse(groups == group_levels[1], 1, -1)
    
    # Compute observed test statistics
    T_obs <- colSums(R * groups_coded)
    
    # Initialize matrix to store permutation test statistics
    T_perm <- matrix(NA, B, p)
    
    # Perform permutations
    for (b in 1:B) {
      permuted_g <- rep(NA, n)
      for (s in unique_strata) {
        idx <- which(strata == s)
        permuted_g[idx] <- sample(groups_coded[idx])
      }
      T_perm[b, ] <- colSums(R * permuted_g)
    }
    
    # Compute p-values
    p_values <- (colSums(abs(T_perm) >= matrix(abs(T_obs), nrow = B, ncol = p, byrow = TRUE)) + 1) / (B + 1)
    
    # Return a list containing the observed test statistics and p-values
    return(list(T_obs = T_obs, p_values = p_values))
  }
  
  dim(grp1)
  
  counts.table <- t(cbind(grp1, grp2)[which.lots,])
  rel.counts.table <- t(apply(counts.table, MARGIN = 1, function(x) x/sum(x)))
  grp.labels <- rep(1:2, c(ncol(grp1), ncol(grp2)))
  strata <- substr(rownames(rel.counts.table), 7, 9)
  
  strat.test <- stratified_permutation_test(rel.counts.table, groups = grp.labels, strata = as.factor(strata), B = 10^4)
  plot(ecdf(strat.test$p_values))
  
  fdrmod <- fdrtool(strat.test$p_values, statistic="pvalue")
  
  perm.pvals <- strat.test$p_values
  saveRDS(perm.pvals, "micro-biome-strat-pvals.rds")
  
}


lfdr_strim <- fdrtool(perm.pvals, statistic="pvalue")
lfdr_est <- lfdr_strim$lfdr
Fdr_est <- lfdr_strim$qval
pi0_est <- lfdr_strim$param[3]
fhat <- pi0_est / lfdr_est

df1 <- data.frame(perm.pvals = perm.pvals, fhat = fhat-pi0_est)
df2 <- data.frame(perm.pvals = sort(perm.pvals), lfdr_est = sort(lfdr_est), Fdr_est=sort(Fdr_est))

# First plot: Histogram of perm.pvals with null and alternative components
ggplot(data = df1, aes(x = perm.pvals)) +
  geom_histogram(aes(y = ..density..), 
                 breaks = seq(0, 1, length.out = 51),  # 50 bins from 0 to 1
                 fill = "gray", 
                 color = "black") +
  scale_y_continuous(breaks = seq(0, 15, by = 2)) + 
  scale_x_continuous(limits = c(0, 1),
                     breaks=seq(0,1,by=0.2)) +  # Ensure the x-axis runs from 0 to 1
  geom_line(aes(x = sort(perm.pvals), 
                y = sort(fhat, decreasing = TRUE),
                color="Alternative component"),
            size = .65) +
  coord_cartesian(ylim = c(0, 15)) + 
  geom_line(aes(x=perm.pvals,
                y=pi0_est,
                color="Null component"),
            linetype="dashed",size=0.65)+
  labs(x = "p-value", y = "Density") +
  scale_color_manual(name = "Legend", values = c("Alternative component" = "blue","Null component"="red")) +
  scale_linetype_manual(values = c("Alternative component" = "solid", 
                                   "Null component" = "dashed")) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),  # Increase axis label size
    axis.text = element_text(size = 14),   # Increase numeric labels size
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title size and center it
    legend.position = c(0.7, 0.8),  # Place legend inside the plot
    legend.background = element_blank(),  # Remove legend background
    legend.title = element_blank(),  # Remove the "Legend" title
    legend.text = element_text(size = 16),  # Increase font size of legend text
    legend.spacing = unit(1, "cm"),  # Increase vertical space between legend items
    legend.key.height = unit(1, "cm"),  # Increase height of the legend keys
    axis.title.y = element_text(margin = margin(r = 10))
  )

# Second plot: lfdr and tail fdr versus p-values
ggplot(data = df2, aes(x = perm.pvals)) +
  geom_line(aes(y=lfdr_est,color = "lfdr")) +
  geom_line(aes(y=Fdr_est,color = "tail Fdr"),linetype="dashed",size=0.65) +
  geom_hline(yintercept = 0, color = "black",size=0.25) +
  geom_vline(xintercept = 0, color = "black",size=0.25) +
  labs(#title = "lfdr estimate",
    x = "p-value", y = "local and tail fdr") +
  scale_y_continuous(breaks = seq(0, 1, by = .2)) + 
  #ylim(0,7.5)+
  scale_x_continuous(breaks = seq(0,1,by=.2)) +  # Ensure the x-axis runs from 0 to 1
  scale_color_manual(name = "Legend", values = c("lfdr" = "black","tail Fdr"="black")) +
  scale_linetype_manual(values = c("lfdr" = "solid", 
                                   "tail Fdr" = "dashed")) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),  # Increase axis label size
    axis.text = element_text(size = 14),   # Increase numeric labels size
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title size and center it
    legend.position = c(0.7, 0.7),  # Place legend inside the plot
    legend.background = element_blank(),  # Remove legend background
    legend.title = element_blank(),  # Remove the "Legend" title
    legend.text = element_text(size = 16),  # Increase font size of legend text
    legend.spacing = unit(1, "cm"),  # Increase vertical space between legend items
    legend.key.height = unit(1, "cm")  # Increase height of the legend keys
  )


############### Figure 2: GGM example: histogram plot

d <- 80
n <- d * 10

set.seed(12)
m <- choose(d, 2)
m1 <- ceiling(m * .05)
non.null.vals <- -(1 + runif(m1)) / sqrt(n) * 2

Theta <- diag(d)
Theta[upper.tri(Theta)] <- sample(c(non.null.vals, rep(0,m-m1))) * 2
Theta <- (Theta + t(Theta))/2

round(Theta,2)[1:10,1:10]
Sigma <- solve(Theta)

Z <- matrix(rnorm(n * d), nrow=n)
X <- Z %*% chol(Sigma)

tstats <- matrix(NA, d, d)
pvals <- matrix(NA, d, d)
for(j in 1:d) {
  mod <- lm(X[,j] ~ 0 + X[,-j])
  tvals <- summary(mod)$coef[,3]
  tstats[-j,j] <- tvals
  pvals[-j,j] <- summary(mod)$coef[,4]
}

vec.tstats <- tstats[upper.tri(tstats)]

bin.width <- 0.2

# Create the histogram with ggplot
df <- data.frame(x=vec.tstats)

# Modify the histogram plot with larger axis labels and custom y-axis ticks
ggplot(data = df, aes(x = x)) +
  geom_histogram(binwidth = 0.2,
                 fill = "gray",
                 color = "black") +
  xlim(1.8, 4.2) +
  ylim(0, dt(2, df = n - d) * m * bin.width * 1.8) +
  labs(title = "Histogram of t-statistics", 
       x = "t-statistic", 
       y = "Count") +
  stat_function(fun = function(x) dt(x, df = n - d) * (m-m1) * bin.width, 
                aes(color = "Expected null count"), size = 0.5) +
  geom_rug(sides = "b", color = "black") +
  scale_y_continuous(breaks = seq(0, dt(2, df = n - d) * m * bin.width * 1.8, by = 10)) +  # Set y-axis breaks by 10
  scale_color_manual(name = "Legend", values = c("Expected null count" = "red")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),  # Increase axis label size
    axis.text = element_text(size = 14),   # Increase numeric labels size
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title size and center it
    legend.position = c(0.7, 0.8),  # Place legend inside the plot (x, y between 0 and 1)
    legend.background = element_blank(),  # Remove legend background
    legend.title = element_blank(),  # Remove the "Legend" title
    legend.text = element_text(size = 16) 
  )


################## Figure 3: GGM example, calibration plot


rds_file <- "ggm_d=80_many_iter.rds"

if (file.exists(rds_file)) {data <- readRDS(rds_file)} else {

d <- 80
n <- d * 10

set.seed(123)
m <- choose(d, 2)
pi0 <- 0.95
m1 <- ceiling(m * (1-pi0))
non.null.vals <- -(1 + runif(m1)) / sqrt(n) * 2
Theta <- diag(d)
Theta[upper.tri(Theta)] <- sample(c(non.null.vals, rep(0,m-m1))) * 2
Theta <- (Theta + t(Theta))/2

round(Theta,2)[1:10,1:10]
Sigma <- solve(Theta)

fdp <- function(x,y,H,vals) {
  R <- sum((vals > x) & (vals < y))
  V <- sum((vals > x) & (vals < y) & (H==0))
  if (R>0) return(V/R) else 0
}

many <- 10^4
ngrid <- 41
xgrid <- seq(0,1,length=ngrid)
calib_p <- rep(NA,ngrid-1)
calib_q <- rep(NA,ngrid-1)
calib_l <- rep(NA,ngrid-1)
pcount_R <- rep(0,ngrid-1)
pcount_V <- rep(0,ngrid-1)
qcount_R <- rep(0,ngrid-1)
qcount_V <- rep(0,ngrid-1)
lcount_R <- rep(0,ngrid-1)
lcount_V <- rep(0,ngrid-1)

for (k in 1:many) {
  print(k)
  
  Z <- matrix(rnorm(n * d), nrow=n)
  X <- Z %*% chol(Sigma)
  
  tstats <- matrix(NA, d, d)
  pvals <- matrix(NA, d, d)
  for(j in 1:d) {
    #print(j)
    mod <- lm(X[,j] ~ 0 + X[,-j])
    tvals <- summary(mod)$coef[,3]
    tstats[-j,j] <- tvals
    pvals[-j,j] <- summary(mod)$coef[,4]
  }
  
  vec.tstats <- tstats[upper.tri(tstats)]
  
  H <- (Theta[upper.tri(Theta)] != 0)+0
  vec.pvals <- pvals[upper.tri(pvals)]
  
  qvals <- p.adjust(vec.pvals,method="BH")
  tmp <- fdrtool(vec.pvals,statistic="pvalue") # strimmer
  lvals <- tmp$lfdr
  for (i in 1:(ngrid-1)) {
    R <- sum((vec.pvals > xgrid[i]) & (vec.pvals < xgrid[i+1]))
    pcount_R[i] <- pcount_R[i] + R
    V <- fdp(xgrid[i],xgrid[i+1],H,vec.pvals)*R
    pcount_V[i] <- pcount_V[i] + V
    
    R <- sum((qvals > xgrid[i]) & (qvals < xgrid[i+1]))
    qcount_R[i] <- qcount_R[i] + R
    V <- fdp(xgrid[i],xgrid[i+1],H,qvals)*R
    qcount_V[i] <- qcount_V[i] + V
    
    R <- sum((lvals > xgrid[i]) & (lvals < xgrid[i+1]))
    lcount_R[i] <- lcount_R[i] + R
    V <- fdp(xgrid[i],xgrid[i+1],H,lvals)*R
    lcount_V[i] <- lcount_V[i] + V
  }
}

for (i in 1:(ngrid-1)) {
  calib_p[i] <- pcount_V[i]/pcount_R[i]
  calib_q[i] <- qcount_V[i]/qcount_R[i]
  calib_l[i] <- lcount_V[i]/lcount_R[i]
}

centering <- (xgrid[1:(ngrid-1)] + xgrid[2:ngrid])/2
data <- data.frame(
  #x = rep(xgrid[1:(ngrid-1)], 3),
  x = rep(centering, 3),
  value = c(calib_p, calib_q, calib_l),
  Statistic = factor(rep(c("p-value", "q-value", "l-value"), 
                         each = length(xgrid[1:(ngrid-1)])),
                     levels = c("p-value", "q-value", "l-value"))
)

saveRDS(data, "ggm_d=80_many_iter.rds") # save the simulation data (takes a while to produce)

}


# Creating the ggplot
p <- ggplot(data, aes(x = x, y = value, color = Statistic)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_color_manual(
    name = "Statistic",
    values = c("p-value" = "black", 
               "q-value" = "red", 
               "l-value" = "blue"),
    labels = c(
      TeX("p-value"), 
      TeX("q-value"), 
      TeX("lfdr estimate")
    )
  ) +
  labs(y = "Fraction of nulls", x = "Reported statistic") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme(
    aspect.ratio = 1,  # Makes the plot square
    legend.position = c(0.8, 0.2),  # Positioning the legend inside the plot area
    legend.background = element_blank(),  # Remove the box around the legend
    legend.key = element_blank(),  # Remove the background of the legend keys
    legend.text.align = 0  # Aligns the legend text to the left
  )

# Display the plot
print(p)



# load data for Nudge example
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

########### Figure 4: P-value histogram with bFDR estimate

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


########### Figure 5: Broken line plot

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



########### Figure 7: Scatterplot of clfdr versus mlfdr with m=1000 tests

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


########### Figure 8: Small experiment to illustrate clfdr with m=6 tests

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


########### Larger experiment to show clfdr and mlfdr \approx local FDP (old manuscript)

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


