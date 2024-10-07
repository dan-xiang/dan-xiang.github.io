# local knockoffs+ simulations

library(knockoff)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Iso)
library(doParallel)



####################### Figure 1: ballot theorem formula for bFDR of local knockoffs(+)

rds_file <- "bFDR_ballot.rds"

if (file.exists(rds_file)) {df_long <- readRDS(rds_file)} else {
  

coefs <- function(q,m) {
  C0 <- 1 # initial condition
  if (m == 0) {
    return(C0)
  } else {
    C <- rep(1,m)
    for (b in 1:m) {
      num <- choose(b,c(0:b))
      den <- choose(floor(q*b)+b-1,c(0:b))
      C[b] <- -sum((num/den)[1:b]*c(C0,C)[1:b])*den[b+1]
    }
    return(c(C0,C))
  }
}

# test case 
coefs(1,10) # should return 1, followed by -1's

# compute the number of binary strings where (# votes for A) > q * (# votes for B)
Nq <- function(a,b,q) {
  if ((b == 0) & (q == Inf)) {
    return(0)
  } else if (a > q*b) {
    C <- coefs(q,b)
    ratio <- choose(b,c(0:b))/choose(a+b-1,c(0:b))
    return(a/(a+b)*sum(C*ratio)*choose(a+b,a))
  } else {
    return(0)
  }
}

# compute the number of binary strings where (# votes for A) > q * [(# votes for B) + 1]
Nq_plus <- function(a,b,q) {
  if ((b == 0) & (q == Inf)) {
    return(0)
  } else if (a > q*(b+1)) {
    C <- coefs(q,b)
    ratio <- choose(b,c(0:b))/choose(a+b-1,c(0:b))
    return(a/(a+b)*sum(C*ratio)*choose(a+b,a))
  } else {
    return(0)
  }
}

# r is the rank of the last rejection
# q is the tuning parameter for the local knockoffs procedure
# p is the number of hypotheses (number of betas)
g <- function(r,q,p) {
  sum1 <- 0
  for (y in 0:r) {
    sum1 <- sum1 + Nq(y,r-y,1/q)
  }
  sum2 <- 0
  if (p > r) {
    if (q!=1) {
      for (z in 0:(p-r)) {
        sum2 <- sum2 + Nq(z,p-r-z,q)
      }
    } else if (q==1) {
      for (z in 0:(p-r)) {
        sum2 <- sum2 + Nq(z+1,p-r+1-(z+1),q)
      }
    }
    
  } else if (p==r) {
    sum2 <- 1
  }
  
  return(sum1*sum2)
}

# test cases
g(1,1,3) # should be 2
g(1,1,4) # should be 3
g(2,1,4) # should be 2
g(1,0,10) # should be 0 (because we take the smallest argmin)
g(3,1,3) # should be 2


# ballot theorem formula for local knockoffs

ballot_loc_kn <- function(alpha,p) {
  out <- 0
  for (r in 1:p) {
    out <- out + g(r,alpha,p)
  }
  return(out*2^(-p))
}


# ballot theorem formula for local knockoffs+

# r is an integer
# q is the tuning parameter for the local knockoffs procedure
# p is the number of hypotheses (number of betas)
# a is the marginal success probability of the (binary) p-value
ballot_pred <- function(r,q,p,a) {
  sum1 <- 0
  for (y in 0:(r-1)) {
    bin_prob <- log(1/2)+(r-1-y)*log(a) + y*log(1-a)
    sum1 <- sum1 + Nq_plus(y+1,r-(y+1),1/q)*exp(bin_prob)
  }
  sum2 <- 0
  if (p > r) {
    if (q!=1) {
      for (z in 0:(p-r)) {
        bin_prob <- log(dbinom(z,size=p-r,prob=a))-log(choose(p-r,z))
        sum2 <- sum2 + Nq(z,p-r-z,q)*exp(bin_prob)
      }
    } else if (q==1) {
      for (z in 0:(p-r)) {
        bin_prob <- log(dbinom(z,size=p-r,prob=a))-log(choose(p-r,z))
        sum2 <- sum2 + Nq(z+1,p-r-z,q)*exp(bin_prob)
      }
    }
  } else if (p==r) {sum2 <- 1}
  return(sum1*sum2)
}



p <- 64
a <- 1/2
m <- 200
alphas <- seq(0,1,length=m)
bFDR_lockn <- rep(NA,m)
bFDR_lockn_plus <- rep(NA,m)
j <- 0
for (alpha in alphas) {
  j <- j+1
  print(j)
  
  # predict bFDR using ballot theorem
  out <- 0
  for (r in 1:p) {
    out <- out + ballot_pred(r,alpha,p,a) 
  }
  bFDR_lockn_plus[j] <- out
  bFDR_lockn[j] <- ballot_loc_kn(alpha,p)
}

# Create the data frame
df <- data.frame(bFDR_lockn = bFDR_lockn,
                 bFDR_lockn_plus = bFDR_lockn_plus,
                 alpha = alphas)

# Reshape the data to long format
df_long <- df %>%
  pivot_longer(cols = c(bFDR_lockn,bFDR_lockn_plus), names_to = "variable", values_to = "value")

saveRDS(df_long, "bFDR_ballot.rds")

}


ggplot(df_long, aes(x = alpha, y = value, color = variable)) +
  geom_line() +  # Increase the line size for better visibility
  geom_vline(xintercept = 0, linetype = "solid", color = "black",size=0.35) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black",size=0.35) + 
  geom_abline(slope = 1, intercept = 0, linetype="dashed") +
  ylim(0, 1) +
  labs(title = "Global null",
       x = "Tuning parameter q",
       y = "bFDR(q)") +
  scale_linetype_manual(values = c("bFDR_lockn" = "solid", "bFDR_lockn_plus" = "solid")) +
  scale_color_manual(values = c("bFDR_lockn" = "black", "bFDR_lockn_plus" = "red"),
                     labels = c("bFDR_lockn" = "Local knockoffs", 
                                "bFDR_lockn_plus" = "Local knockoffs+")) +
  guides(linetype = "none") +  # Remove linetype from legend
  theme_minimal() +
  theme(
    #legend.position = c(0.8, 0.3),  # Place legend in bottom right
    legend.position = "none",  # remove legend
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 14),   # Increase x-axis numbers size
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 20),  # Increase plot title size
    legend.text = element_text(size = 16),  # Increase legend text size
    legend.title = element_blank(),
    legend.spacing = unit(1, "cm"),
    legend.key.height = unit(1, "cm")
  )





######### Generate a plot of bFDR(q) versus q for trim_kn+ vs loc_kn+


rds_file <- "trimmed-lockn-bFDR.rds"

if (file.exists(rds_file)) {df_long <- readRDS(rds_file)} else {

set.seed(1)
# Problem parameters
#n=150
#p=64
#s=16
#n = 50          # number of observations
#p = 16           # number of variables
#s = 4            # number of variables with nonzero coefficients
n = 500          # number of observations
p = 100           # number of variables
s = 10            # number of variables with nonzero coefficients
amplitude = 4.5   # signal amplitude (for noise level = 1)
#amplitude = 1.75   # signal amplitude (for noise level = 1)
#amplitude = 10^3   # signal amplitude (for noise level = 1)

# Generate the variables from a multivariate normal distribution
rho = 0 # 0.6
Sigma = toeplitz(rho^(0:(p-1)))
N <- 10^4
m <- 20
qs <- seq(1/m,1,length=m)
bFDR_local_kn <- rep(NA,m)
se_local_kn <- rep(NA,m)
bFDR_trimmed_kn <- rep(NA,m)
se_trimmed_kn <- rep(NA,m)
count <- 0
for (q in qs) {
  count <- count + 1
  print(count)
  bFDR <- rep(NA,N)
  lfdr_bdry <- rep(NA,N)
  for (k in 1:N) {
    #print(k)
    X = matrix(rnorm(n*p),n)  %*% chol(Sigma)
    # Generate the response from a linear model
    nonzero = sample(p, s)
    beta = amplitude * (1:p %in% nonzero) / sqrt(n)
    y.sample = function(X) X %*% beta + rnorm(n)
    y = y.sample(X)
    alpha=q
    result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=alpha)
    #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_lambdasmax,fdr=alpha)
    #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic =  stat.lasso_lambdadiff,fdr=alpha)
    #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdadiff,fdr=alpha)
    #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_coefdiff,fdr=alpha)
    W <- result$statistic
    js <- 1:p
    inds <- order(abs(result$statistic))
    # compute boundary FDR for local knock-offs
    b <- (1-sign(rev(W[inds]))) / 2
    Bk <- cumsum(b)
    Zk <- 1/(1+alpha) + Bk - c(1:p)*alpha/(1+alpha)
    R <- if (sum(Zk<0)>0) which.min(Zk) else 0
    bFDR[k] <- if (R>0) 1-(rev(beta[inds])[R] != 0) else 0
    n_pos <- sum(W>0)
    r <- sum(1-b[1:R])
    W_bdry <- if (r > 0) W[W>0][order(W[W>0])][n_pos-r+1] else Inf
    bdry_ind <- if (W_bdry < Inf) which(rev(W[inds])==W_bdry)[1] else 0
    lfdr_bdry[k] <- if (bdry_ind > 0) {
      # if bdry_ind == p, then set bdry_int <- bdry_ind - 1
      while(pava(c(1,b[min(p,(bdry_ind+1)):p]))[1] > alpha/(1+alpha) || b[bdry_ind]==1) {
        bdry_ind <- bdry_ind - 1
        if (bdry_ind == 0) break
      }
      # check the null status of the last 0 in the rejection set
      rej_inds <- which(rev(b[inds])[1:bdry_ind]==0)
      last_rej <- if (length(rej_inds) > 0) max(rej_inds) else 0
      if (last_rej > 0) 1-(rev(beta[inds])[last_rej] != 0) else 0
    } else 0
  }
  bFDR_local_kn[count] <- mean(bFDR)
  se_local_kn[count] <- sd(bFDR)/sqrt(N)
  bFDR_trimmed_kn[count] <- mean(lfdr_bdry)
  se_trimmed_kn[count] <- sd(lfdr_bdry)/sqrt(N)
}

bFDR_trimmed_kn
se_trimmed_kn
df <- data.frame(bFDR_local_kn = bFDR_local_kn,
                 bFDR_trimmed_kn = bFDR_trimmed_kn,
                 q = qs)

# Reshape the data to long format
df_long <- df %>%
  pivot_longer(cols = c(bFDR_local_kn,bFDR_trimmed_kn), names_to = "variable", values_to = "value")

saveRDS(df_long, "trimmed-lockn-bFDR.rds")

}

# Plot using ggplot2
ggplot(df_long, aes(x = q, y = value, color = variable, linetype = variable)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black",size=0.35) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black",size=0.35) + 
  geom_point(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value), size = 1) +  # Plot bFDR_pred as a line
  geom_point(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value), size= 1) +  # Plot bFDR_true as points
  geom_line(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value, color = variable), size = 0.5, show.legend = FALSE) +  # connect dots
  geom_line(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value, color = variable), size = 0.5,linetype = "solid", show.legend = FALSE) +  # connnect dots
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Add the dotted diagonal line
  labs(title = "Non-exchangeable",
       x = "q",
       y = "bFDR(q)",
       color = "Variable") +
  scale_color_manual(values = c("bFDR_local_kn" = "red", "bFDR_trimmed_kn" = "blue"),
                     labels = c("bFDR_local_kn" = "Local knockoffs+", 
                                "bFDR_trimmed_kn" = "Trimmed knockoffs+")) +
  theme_minimal()+
  ylim(0, 1) + 
  theme(
    plot.title = element_text(size = 20),  # Increase title size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),   # Increase x-axis numbers size
    axis.text.y = element_text(size = 14),
    #legend.text = element_text(size = 16),   # Increase legend text size
    legend.title = element_blank(), # Increase y-axis label size
    #legend.position = c(0.75, 0.3),
    legend.position = "none", # remove legend
    legend.spacing = unit(1, "cm"),
    legend.key.height = unit(1, "cm")
  )





############# Figure 2: Plot of W-statistics among rejected hypotheses


rds_file <- "knockoff-W-stats.rds"

if (file.exists(rds_file)) {df <- readRDS(rds_file)} else {

set.seed(123)
# Problem parameters
n = 2000          # number of observations
p = 800           # number of variables
k = 80            # number of variables with nonzero coefficients
#amplitude = 8   # signal amplitude (for noise level = 1)
amplitude = 4.5   # signal amplitude (for noise level = 1)

# Generate the variables from a multivariate normal distribution
mu = rep(0,p)
rho = 0.25
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

# Generate the response from a linear model
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero) / sqrt(n)
tmp <- rep("null",p)
tmp[nonzero] <- "non-null"
truth_value <- as.factor(tmp)
fdp = function(selected) sum(selected == 0) / max(1, length(selected))

y.sample = function(X) X %*% beta + rnorm(n)
y = y.sample(X)

alpha=0.2
ell=alpha/(1+alpha)
result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=alpha)
W <- result$statistic
js <- 1:p
inds <- order(abs(result$statistic))
R <- length(result$selected)

par(mfrow=c(1,1))
n_pos <- sum(W>0)
W_bdry <- W[W>0][order(W[W>0])][n_pos-R+1]
lastR <- (p-R+1):p
big_inds <- which(abs(W)>=W_bdry)
R2 <- length(big_inds)
inds2 <- order(abs(W[big_inds]))
J <- length(big_inds)
W_stats <- rev(W[big_inds][inds2])
truth_stats <- rev(truth_value[big_inds][inds2])

# Assuming W_stats, truth_stats, and J are defined, create a data frame for ggplot
df <- data.frame(
  index = 1:J,
  W = W_stats,
  truth = factor(truth_stats)  # Convert truth_stats to a factor for coloring
)

saveRDS(df, "knockoff-W-stats.rds")
}

# Create the plot in ggplot2
ggplot(data = df, aes(x = index, y = W, color = truth)) +
  geom_point(size=1.3) +
  labs(title = "W-statistics ordered by absolute value",
       x = "Rank",
       y = "W") +
  scale_color_manual(values = c("null" = "red", "non-null" = "black"),  # Reverse colors: 1 = red, 2 = blue
                     labels = c("Non-null", "Null")) + 
  geom_hline(yintercept = 0, color = "black",size=.25) +
  geom_vline(xintercept = 0, color = "black",size=.25) +
  scale_x_continuous(breaks = seq(0,90,by=10)) + 
  scale_y_continuous(breaks = seq(-450,450,by=150)) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),  # Increase axis label size
    axis.text = element_text(size = 14),   # Increase numeric labels size
    axis.title.x = element_text(margin = margin(t = 15)),
    plot.title = element_text(size = 18, hjust = 0.5),  # Increase title size and center it
    legend.position = c(0.25, 0.53),  # Place legend inside the plot
    legend.background = element_blank(),  # Remove legend background
    legend.title = element_blank(),  # Remove the "Legend" title
    legend.text = element_text(size = 16),  # Increase font size of legend text
    legend.spacing = unit(1, "cm"),  # Increase vertical space between legend items
    legend.key.height = unit(1, "cm")  # Increase height of the legend keys
  )






#################### Figure 3: Support line characterization of local knockoffs


rds_file <- "support-line-picture.rds"

if (file.exists(rds_file)) {data <- readRDS(rds_file)} else {

set.seed(9)
n = 1500          # number of observations
p = 350           # number of variables
k = 35            # number of variables with nonzero coefficients
amplitude = 4.5   # signal amplitude (for noise level = 1)

# Generate the variables from a multivariate normal distribution
mu = rep(0,p)
rho = 0.25
Sigma = toeplitz(rho^(0:(p-1)))
X = matrix(rnorm(n*p),n) %*% chol(Sigma)

# Generate the response from a linear model
nonzero = sample(p, k)
beta = amplitude * (1:p %in% nonzero) / sqrt(n)
fdp = function(selected) sum(selected == 0) / max(1, length(selected))

N <- 1000
bFDR <- rep(NA,N)
y.sample = function(X) X %*% beta + rnorm(n)
y = y.sample(X)

alpha=0.3
result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=alpha)
W <- result$statistic
js <- 1:p
inds <- order(abs(result$statistic))
R <- length(result$selected)
raw_lvals <- (rev(W[inds])<0)+0
iso_lfdr <- pava(raw_lvals,decreasing=F)
m<-100
b <- (1-sign(rev(W[inds])))/2
Bk <- cumsum(b)
Zk <- Bk-c(1:p)*alpha/(1+alpha)
trimR <- if (sum(Zk< -1/(1+alpha))>0) which.min(Zk) else 0

gcm <- c(0,cumsum(iso_lfdr)[2:length(iso_lfdr)])

# Assuming Bk is your vector of data
# Create a data frame for the first M elements
M <- 70
data <- data.frame(Index = 1:M, 
                   Value = Bk[1:M], 
                   GCM <- gcm[1:M],
                   secant <- alpha/(1+alpha)*c(1:M),
                   tangent <- alpha/(1+alpha)*(c(1:M)-trimR)+Bk[trimR])

saveRDS(data, "support-line-picture.rds")
}

# Plot the data using ggplot2
ggplot(data, aes(x = Index)) +
  geom_point(aes(y=Value)) + # Use geom_line() for line plots
  geom_vline(xintercept = 0, linetype = "solid", color = "black",size=0.25) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black",size=0.25) + 
  geom_line(aes(y = GCM), color = "black") +
  geom_line(aes(y = secant), color = "red") +
  geom_line(aes(y = tangent), color = "blue") +
  theme_minimal() + # Optional: Sets a minimalistic theme
  labs(x = "k", y = "Bk")+
  theme(
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    axis.text.x = element_text(size = 14),   # Increase x-axis numbers size
    axis.text.y = element_text(size = 14),   # Increase y-axis numbers size
    plot.title = element_text(size = 18)     # Increase plot title size
  )






############# Figure 5: bFDR curves for correlated features, and different W-statistic

rds_file <- "trimmed-lockn-bFDR-LCD.rds"

if (file.exists(rds_file)) {df_long <- readRDS(rds_file)} else {
  
  set.seed(1)
  # Problem parameters
  n = 150
  p = 64
  s = 16
  amplitude = 4.5   # signal amplitude
  
  # Generate the variables from a multivariate normal distribution
  rho = 0 # 0.6
  Sigma = toeplitz(rho^(0:(p-1)))
  N <- 10^4
  m <- 20
  qs <- seq(1/m,1,length=m)
  bFDR_local_kn <- rep(NA,m)
  se_local_kn <- rep(NA,m)
  bFDR_trimmed_kn <- rep(NA,m)
  se_trimmed_kn <- rep(NA,m)
  count <- 0
  for (q in qs) {
    count <- count + 1
    print(count)
    bFDR <- rep(NA,N)
    lfdr_bdry <- rep(NA,N)
    for (k in 1:N) {
      #print(k)
      X = matrix(rnorm(n*p),n)  %*% chol(Sigma)
      # Generate the response from a linear model
      nonzero = sample(p, s)
      beta = amplitude * (1:p %in% nonzero) / sqrt(n)
      y.sample = function(X) X %*% beta + rnorm(n)
      y = y.sample(X)
      alpha=q
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_lambdasmax,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic =  stat.lasso_lambdadiff,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdadiff,fdr=alpha)
      result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_coefdiff,fdr=alpha)
      W <- result$statistic
      js <- 1:p
      inds <- order(abs(result$statistic))
      # compute boundary FDR for local knock-offs
      b <- (1-sign(rev(W[inds]))) / 2
      Bk <- cumsum(b)
      Zk <- 1/(1+alpha) + Bk - c(1:p)*alpha/(1+alpha)
      R <- if (sum(Zk<0)>0) which.min(Zk) else 0
      bFDR[k] <- if (R>0) 1-(rev(beta[inds])[R] != 0) else 0
      n_pos <- sum(W>0)
      r <- sum(1-b[1:R])
      W_bdry <- if (r > 0) W[W>0][order(W[W>0])][n_pos-r+1] else Inf
      bdry_ind <- if (W_bdry < Inf) which(rev(W[inds])==W_bdry)[1] else 0
      lfdr_bdry[k] <- if (bdry_ind > 0) {
        # if bdry_ind == p, then set bdry_int <- bdry_ind - 1
        while(pava(c(1,b[min(p,(bdry_ind+1)):p]))[1] > alpha/(1+alpha) || b[bdry_ind]==1) {
          bdry_ind <- bdry_ind - 1
          if (bdry_ind == 0) break
        }
        # check the null status of the last 0 in the rejection set
        rej_inds <- which(rev(b[inds])[1:bdry_ind]==0)
        last_rej <- if (length(rej_inds) > 0) max(rej_inds) else 0
        if (last_rej > 0) 1-(rev(beta[inds])[last_rej] != 0) else 0
      } else 0
    }
    bFDR_local_kn[count] <- mean(bFDR)
    se_local_kn[count] <- sd(bFDR)/sqrt(N)
    bFDR_trimmed_kn[count] <- mean(lfdr_bdry)
    se_trimmed_kn[count] <- sd(lfdr_bdry)/sqrt(N)
  }
  
  bFDR_trimmed_kn
  se_trimmed_kn
  df <- data.frame(bFDR_local_kn = bFDR_local_kn,
                   bFDR_trimmed_kn = bFDR_trimmed_kn,
                   q = qs)
  
  # Reshape the data to long format
  df_long <- df %>%
    pivot_longer(cols = c(bFDR_local_kn,bFDR_trimmed_kn), names_to = "variable", values_to = "value")
  
  saveRDS(df_long, "trimmed-lockn-bFDR-LCD.rds")
  
}

# Plot using ggplot2
ggplot(df_long, aes(x = q, y = value, color = variable, linetype = variable)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black",size=0.35) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black",size=0.35) + 
  geom_point(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value), size = 1) +  # Plot bFDR_pred as a line
  geom_point(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value), size= 1) +  # Plot bFDR_true as points
  geom_line(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value, color = variable), size = 0.5, show.legend = FALSE) +  # connect dots
  geom_line(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value, color = variable), size = 0.5,linetype = "solid", show.legend = FALSE) +  # connnect dots
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Add the dotted diagonal line
  labs(title = "Lasso coefficient difference",
       x = "q",
       y = "bFDR(q)",
       color = "Variable") +
  scale_color_manual(values = c("bFDR_local_kn" = "red", "bFDR_trimmed_kn" = "blue"),
                     labels = c("bFDR_local_kn" = "Local knockoffs+", 
                                "bFDR_trimmed_kn" = "Trimmed knockoffs+")) +
  theme_minimal()+
  ylim(0, 1) + 
  theme(
    plot.title = element_text(size = 20),  # Increase title size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),   # Increase x-axis numbers size
    axis.text.y = element_text(size = 14),
    #legend.text = element_text(size = 16),   # Increase legend text size
    legend.title = element_blank(), # Increase y-axis label size
    #legend.position = c(0.75, 0.3),
    legend.position = "none", # remove legend
    legend.spacing = unit(1, "cm"),
    legend.key.height = unit(1, "cm")
  )







rds_file <- "trimmed-lockn-bFDR-cor.rds"

if (file.exists(rds_file)) {df_long <- readRDS(rds_file)} else {
  
  set.seed(1)
  # Problem parameters
  n = 150
  p = 64
  s = 16
  amplitude = 4.5   # signal amplitude
  
  # Generate the variables from a multivariate normal distribution
  rho = 0.6
  Sigma = toeplitz(rho^(0:(p-1)))
  N <- 10^4
  m <- 20
  qs <- seq(1/m,1,length=m)
  bFDR_local_kn <- rep(NA,m)
  se_local_kn <- rep(NA,m)
  bFDR_trimmed_kn <- rep(NA,m)
  se_trimmed_kn <- rep(NA,m)
  count <- 0
  for (q in qs) {
    count <- count + 1
    print(count)
    bFDR <- rep(NA,N)
    lfdr_bdry <- rep(NA,N)
    for (k in 1:N) {
      #print(k)
      X = matrix(rnorm(n*p),n)  %*% chol(Sigma)
      # Generate the response from a linear model
      nonzero = sample(p, s)
      beta = amplitude * (1:p %in% nonzero) / sqrt(n)
      y.sample = function(X) X %*% beta + rnorm(n)
      y = y.sample(X)
      alpha=q
      result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_lambdasmax,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic =  stat.lasso_lambdadiff,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdadiff,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_coefdiff,fdr=alpha)
      W <- result$statistic
      js <- 1:p
      inds <- order(abs(result$statistic))
      # compute boundary FDR for local knock-offs
      b <- (1-sign(rev(W[inds]))) / 2
      Bk <- cumsum(b)
      Zk <- 1/(1+alpha) + Bk - c(1:p)*alpha/(1+alpha)
      R <- if (sum(Zk<0)>0) which.min(Zk) else 0
      bFDR[k] <- if (R>0) 1-(rev(beta[inds])[R] != 0) else 0
      n_pos <- sum(W>0)
      r <- sum(1-b[1:R])
      W_bdry <- if (r > 0) W[W>0][order(W[W>0])][n_pos-r+1] else Inf
      bdry_ind <- if (W_bdry < Inf) which(rev(W[inds])==W_bdry)[1] else 0
      lfdr_bdry[k] <- if (bdry_ind > 0) {
        # if bdry_ind == p, then set bdry_int <- bdry_ind - 1
        while(pava(c(1,b[min(p,(bdry_ind+1)):p]))[1] > alpha/(1+alpha) || b[bdry_ind]==1) {
          bdry_ind <- bdry_ind - 1
          if (bdry_ind == 0) break
        }
        # check the null status of the last 0 in the rejection set
        rej_inds <- which(rev(b[inds])[1:bdry_ind]==0)
        last_rej <- if (length(rej_inds) > 0) max(rej_inds) else 0
        if (last_rej > 0) 1-(rev(beta[inds])[last_rej] != 0) else 0
      } else 0
    }
    bFDR_local_kn[count] <- mean(bFDR)
    se_local_kn[count] <- sd(bFDR)/sqrt(N)
    bFDR_trimmed_kn[count] <- mean(lfdr_bdry)
    se_trimmed_kn[count] <- sd(lfdr_bdry)/sqrt(N)
  }
  
  bFDR_trimmed_kn
  se_trimmed_kn
  df <- data.frame(bFDR_local_kn = bFDR_local_kn,
                   bFDR_trimmed_kn = bFDR_trimmed_kn,
                   q = qs)
  
  # Reshape the data to long format
  df_long <- df %>%
    pivot_longer(cols = c(bFDR_local_kn,bFDR_trimmed_kn), names_to = "variable", values_to = "value")
  
  saveRDS(df_long, "trimmed-lockn-bFDR-cor.rds")
  
}

# Plot using ggplot2
ggplot(df_long, aes(x = q, y = value, color = variable, linetype = variable)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black",size=0.35) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black",size=0.35) + 
  geom_point(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value), size = 1) +  # Plot bFDR_pred as a line
  geom_point(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value), size= 1) +  # Plot bFDR_true as points
  geom_line(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value, color = variable), size = 0.5, show.legend = FALSE) +  # connect dots
  geom_line(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value, color = variable), size = 0.5,linetype = "solid", show.legend = FALSE) +  # connnect dots
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Add the dotted diagonal line
  labs(title = "Correlated case",
       x = "q",
       y = "bFDR(q)",
       color = "Variable") +
  scale_color_manual(values = c("bFDR_local_kn" = "red", "bFDR_trimmed_kn" = "blue"),
                     labels = c("bFDR_local_kn" = "Local knockoffs+", 
                                "bFDR_trimmed_kn" = "Trimmed knockoffs+")) +
  theme_minimal()+
  ylim(0, 1) + 
  theme(
    plot.title = element_text(size = 20),  # Increase title size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),   # Increase x-axis numbers size
    axis.text.y = element_text(size = 14),
    #legend.text = element_text(size = 16),   # Increase legend text size
    legend.title = element_blank(), # Increase y-axis label size
    #legend.position = c(0.75, 0.3),
    legend.position = "none", # remove legend
    legend.spacing = unit(1, "cm"),
    legend.key.height = unit(1, "cm")
  )






############### iid entries of X


rds_file <- "trimmed-lockn-bFDR-ind.rds"

if (file.exists(rds_file)) {df_long <- readRDS(rds_file)} else {
  
  set.seed(1)
  # Problem parameters
  n = 150
  p = 64
  s = 16
  amplitude = 4.5   # signal amplitude
  
  # Generate the variables from a multivariate normal distribution
  rho = 0 # 0.6
  Sigma = toeplitz(rho^(0:(p-1)))
  N <- 10^4
  m <- 20
  qs <- seq(1/m,1,length=m)
  bFDR_local_kn <- rep(NA,m)
  se_local_kn <- rep(NA,m)
  bFDR_trimmed_kn <- rep(NA,m)
  se_trimmed_kn <- rep(NA,m)
  count <- 0
  for (q in qs) {
    count <- count + 1
    print(count)
    bFDR <- rep(NA,N)
    lfdr_bdry <- rep(NA,N)
    for (k in 1:N) {
      #print(k)
      X = matrix(rnorm(n*p),n)  %*% chol(Sigma)
      # Generate the response from a linear model
      nonzero = sample(p, s)
      beta = amplitude * (1:p %in% nonzero) / sqrt(n)
      y.sample = function(X) X %*% beta + rnorm(n)
      y = y.sample(X)
      alpha=q
      result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_lambdasmax,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic =  stat.lasso_lambdadiff,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdadiff,fdr=alpha)
      #result = knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.lasso_coefdiff,fdr=alpha)
      W <- result$statistic
      js <- 1:p
      inds <- order(abs(result$statistic))
      # compute boundary FDR for local knock-offs
      b <- (1-sign(rev(W[inds]))) / 2
      Bk <- cumsum(b)
      Zk <- 1/(1+alpha) + Bk - c(1:p)*alpha/(1+alpha)
      R <- if (sum(Zk<0)>0) which.min(Zk) else 0
      bFDR[k] <- if (R>0) 1-(rev(beta[inds])[R] != 0) else 0
      n_pos <- sum(W>0)
      r <- sum(1-b[1:R])
      W_bdry <- if (r > 0) W[W>0][order(W[W>0])][n_pos-r+1] else Inf
      bdry_ind <- if (W_bdry < Inf) which(rev(W[inds])==W_bdry)[1] else 0
      lfdr_bdry[k] <- if (bdry_ind > 0) {
        # if bdry_ind == p, then set bdry_int <- bdry_ind - 1
        while(pava(c(1,b[min(p,(bdry_ind+1)):p]))[1] > alpha/(1+alpha) || b[bdry_ind]==1) {
          bdry_ind <- bdry_ind - 1
          if (bdry_ind == 0) break
        }
        # check the null status of the last 0 in the rejection set
        rej_inds <- which(rev(b[inds])[1:bdry_ind]==0)
        last_rej <- if (length(rej_inds) > 0) max(rej_inds) else 0
        if (last_rej > 0) 1-(rev(beta[inds])[last_rej] != 0) else 0
      } else 0
    }
    bFDR_local_kn[count] <- mean(bFDR)
    se_local_kn[count] <- sd(bFDR)/sqrt(N)
    bFDR_trimmed_kn[count] <- mean(lfdr_bdry)
    se_trimmed_kn[count] <- sd(lfdr_bdry)/sqrt(N)
  }
  
  bFDR_trimmed_kn
  se_trimmed_kn
  df <- data.frame(bFDR_local_kn = bFDR_local_kn,
                   bFDR_trimmed_kn = bFDR_trimmed_kn,
                   q = qs)
  
  # Reshape the data to long format
  df_long <- df %>%
    pivot_longer(cols = c(bFDR_local_kn,bFDR_trimmed_kn), names_to = "variable", values_to = "value")
  
  saveRDS(df_long, "trimmed-lockn-bFDR-ind.rds")
  
}

# Plot using ggplot2
ggplot(df_long, aes(x = q, y = value, color = variable, linetype = variable)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black",size=0.35) +  
  geom_hline(yintercept = 0, linetype = "solid", color = "black",size=0.35) + 
  geom_point(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value), size = 1) +  # Plot bFDR_pred as a line
  geom_point(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value), size= 1) +  # Plot bFDR_true as points
  geom_line(data = subset(df_long, variable == "bFDR_local_kn"), aes(x = q, y = value, color = variable), size = 0.5, show.legend = FALSE) +  # connect dots
  geom_line(data = subset(df_long, variable == "bFDR_trimmed_kn"), aes(x = q, y = value, color = variable), size = 0.5,linetype = "solid", show.legend = FALSE) +  # connnect dots
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Add the dotted diagonal line
  labs(title = "Independent case",
       x = "q",
       y = "bFDR(q)",
       color = "Variable") +
  scale_color_manual(values = c("bFDR_local_kn" = "red", "bFDR_trimmed_kn" = "blue"),
                     labels = c("bFDR_local_kn" = "Local knockoffs+", 
                                "bFDR_trimmed_kn" = "Trimmed knockoffs+")) +
  theme_minimal()+
  ylim(0, 1) + 
  theme(
    plot.title = element_text(size = 20),  # Increase title size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),   # Increase x-axis numbers size
    axis.text.y = element_text(size = 14),
    #legend.text = element_text(size = 16),   # Increase legend text size
    legend.title = element_blank(), # Increase y-axis label size
    #legend.position = c(0.75, 0.3),
    legend.position = "none", # remove legend
    legend.spacing = unit(1, "cm"),
    legend.key.height = unit(1, "cm")
  )


















