# =============================================================================
# Master script: Local FDR estimation with an unknown symmetric null
# Reproduces all figures in the paper (paper figure numbering)
# =============================================================================
#
# REQUIRED FILES (place in working directory before running):
#   stacked_W_lsm.rds                 -- 1000 x 1000 LSM W-statistics (Figs 7, 9)
#   W-and-betas-list-largep.rds       -- list: [[1]]=beta, [[2]]=W matrix (Fig 6)
#   bootstrap-knockoffs.rds           -- precomputed bootstrap knockoff matrices
#   bootstrap-knockoffs-OLS-indep.rds -- OLS bootstrap, independent X (appendix)
#   bootstrap-knockoffs-OLS-cor.rds   -- OLS bootstrap, correlated X (appendix)
#   bdhi8.csv                         -- proteomics data (Fig 10)
#
# NOTE ON SLOW SECTIONS:
#   Sections marked [SLOW] take minutes to hours. They save .rds files so
#   subsequent runs can skip the computation and load from disk.
#
# NOTE ON METHODS:
#   Figure 6 uses the inline 7th-degree polynomial fit (computationally
#   tractable on ~1M W-statistics, where the constrained spline method would
#   be too slow).
#   Figure 10 uses compute.clar.poly() (older natural cubic spline, df=3).
#   All other figures use compute.clar() (current constrained-spline method
#   with adaptive knot placement and left-boundary zero-derivative constraint).
#
# FIGURE NUMBERING (paper order):
#   Figure  1 -- HIV barplot, three drugs (APV/IDV/NFV)
#   Figure  2 -- NFV illustration: histogram + isotonic + logistic regression estimate
#   Figure  5 -- Gaussian two-groups: comparison of clar estimators
#   Figure  6 -- Knockoffs LSM histogram + logistic regression + true clar/lfdr
#   Figure  7 -- Knockoffs single realization clar estimate
#   Figure  8 -- Gaussian: bootstrap and true SEs
#   Figure  9 -- Knockoffs: bootstrap and true SEs
#   Figure 10 -- Proteomics (BDHI-8) clar estimate
#   Figure 11 -- Pooled HIV clar with bootstrap SE bands
#
# Figures 3 and 4 are theoretical illustrations; no R code required.
# Appendix figures (per-drug clar, ecdf diagnostics, delta method, pairs
# bootstrap) appear at the end.
# =============================================================================


# =============================================================================
# SECTION 0: Libraries and shared functions
# =============================================================================

library(knockoff)
library(Iso)
library(splines)
library(glmnet)

# --- Current constrained spline method (used in Figs 2, 5, 7, 8, 9, 11) ---

compute.clar <- function(w, thresh = c(0.9, 0.7, 0.3, 0.1),
                         min.pts.betw.knots = 30) {
  df <- data.frame(x = abs(w),
                   y = 1 - (sign(w) + 1) / 2)
  out <- pava(df$y[order(df$x)], decreasing = TRUE)
  inds <- order(df$x)
  w.inds <- rep(NA, length(thresh))
  for (i in 1:length(thresh)) {
    t <- thresh[i]
    w.inds[i] <- min(which(out / (1 - out) < t))
  }
  w.inds <- sort(w.inds)
  while (min(diff(w.inds)) < min.pts.betw.knots && length(w.inds) > 2) {
    w.inds <- w.inds[-(which.min(diff(w.inds)) + 1)]
  }
  pava.knots <- unique(df$x[inds][w.inds])
  diff1 <- pava.knots[2] - pava.knots[1]
  diff2 <- pava.knots[length(pava.knots)] - pava.knots[length(pava.knots) - 1]
  left_bdry  <- max(pava.knots[1] - diff1, pava.knots[1] / 2)
  right_bdry <- pava.knots[length(pava.knots)] + diff2
  boundary_knots <- c(left_bdry, right_bdry)
  df$cubic <- ns(df$x, knots = pava.knots, Boundary.knots = boundary_knots)
  X <- df$cubic
  left_boundary <- boundary_knots[1]
  x_left1 <- left_boundary - 0.1
  x_left2 <- left_boundary - 0.2
  X_left1 <- predict(df$cubic, x_left1)
  X_left2 <- predict(df$cubic, x_left2)
  basis_derivs <- (X_left1 - X_left2) / (x_left1 - x_left2)
  p <- ncol(X)
  d <- basis_derivs
  X_constrained <- X[, -p] - X[, p] %*% t(d[-p] / d[p])
  glmod <- glm(df$y ~ X_constrained, family = binomial())
  p.hat <- glmod$fitted.values
  odds  <- p.hat / (1 - p.hat)
  x.inds <- order(df$x)
  list(x = df$x[x.inds], y = odds[x.inds],
       gren.y = out / (1 - out),
       raw.y = df$y, knots = pava.knots,
       boundary_knots = boundary_knots, glmod = glmod,
       knot.inds = w.inds,
       pava.odds.at.knots = (out / (1 - out))[w.inds])
}

# --- Older natural cubic spline method (df=3, used in Fig 10) ---

compute.clar.poly <- function(w) {
  df <- data.frame(x = abs(w), y = 1 - (sign(w) + 1) / 2)
  df$cubic <- ns(df$x, df = 3)
  glmod <- glm(y ~ cubic, family = binomial, data = df)
  p.hat <- glmod$fitted.values
  odds  <- p.hat / (1 - p.hat)
  x.inds <- order(df$x)
  list(x = df$x[x.inds], y = odds[x.inds], glmod = glmod,
       p.hat = p.hat, clar.hat = odds)
}

# --- KDE and isotonic helpers ---

kde.clar <- function(w) {
  dens <- density(w, kernel = "epanechnikov")
  kde  <- approxfun(dens$x, dens$y)
  max.w <- max(w)
  x <- seq(0, max.w, length = 10^4)
  clar.kde <- sapply(x, function(xi) kde(-xi) / kde(xi))
  data.frame(x = x, y = clar.kde)
}

isotonic.clar <- function(w) {
  df  <- data.frame(x = abs(w), y = 1 - (sign(w) + 1) / 2)
  out <- pava(df$y[order(df$x)], decreasing = TRUE)
  data.frame(x = df$x[order(df$x)], y = out / (1 - out))
}

true.clar <- function(w, pi0, mu, sigma) {
  f     <- pi0 * dnorm(abs(w), sd = sigma) + (1 - pi0) * dnorm(abs(w), mean = mu, sd = sigma)
  f.neg <- pi0 * dnorm(abs(w), sd = sigma) + (1 - pi0) * dnorm(-abs(w), mean = mu, sd = sigma)
  inds  <- order(abs(w))
  data.frame(x = abs(w)[inds], y = (f.neg / f)[inds])
}

# --- Delta method standard errors ---

delta.method <- function(glmod) {
  p.hat    <- glmod$fitted.values
  clar.hat <- p.hat / (1 - p.hat)
  Sig2.hat <- vcov(glmod)
  X        <- model.matrix(glmod)
  nab.g    <- as.matrix(clar.hat * X)
  n <- nrow(X)
  se <- sapply(1:n, function(i) sqrt(nab.g[i, ] %*% Sig2.hat %*% nab.g[i, ]))
  list(clar.hat = clar.hat,
       lower.bd = clar.hat - 2 * se,
       upper.bd = clar.hat + 2 * se)
}

# --- Shared graphical parameters ---

par_std <- function() {
  par(mar = c(3.6, 3.6, 0.5, 0.5), mgp = c(2.2, 0.6, 0), las = 0)
}
std_cex.lab    <- 1.25
std_cex.legend <- 1.22
std_lwd        <- 1.25
std_lwd.dashed <- 1
std_cex_points <- 1

# =============================================================================
# HIV DATA LOADING (used by Figures 1, 2, 11, and appendix per-drug plots)
# =============================================================================
# Loads protease-inhibitor (PI) drug class data from the Stanford HIVDB
# archive, runs the knockoff filter for each drug, and applies per-drug
# outlier removal. The resulting list `results.w` is used downstream.
# =============================================================================

drug_class <- 'PI'
base_url <- 'https://web.archive.org/web/20170701000000/http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
gene_url <- paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep = '/')
tsm_url  <- paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep = '/')

gene_df <- read.delim(gene_url, na.strings = c('NA', ''), stringsAsFactors = FALSE)
tsm_df  <- read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) <- c('Position', 'Mutations')
tsm_df$Position <- as.numeric(tsm_df$Position)

# --- Helpers ---

grepl_rows <- function(pattern, df) {
  cell_matches <- apply(df, c(1, 2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}

flatten_matrix <- function(M, sep = '.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep = sep)))
  x
}

get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

knockoff_threshold <- function(w, q) {
  t.vals <- sort(unique(abs(w[w != 0])))
  for (t in t.vals) {
    if ((1 + sum(w <= -t)) / max(1, sum(w >= t)) <= q) return(t)
  }
  return(Inf)
}

knockoff_fun <- function(X, y, q = 0.20) {
  y <- log(as.numeric(as.character(unlist(y))))
  keep <- !is.na(y); X <- X[keep, ]; y <- y[keep]
  X <- X[, colSums(X) >= 3]
  X <- X[, colSums(abs(cor(X) - 1) < 1e-4) == 1]
  knock.gen <- function(x) create.fixed(x, method = 'equi')
  result <- knockoff.filter(X, y, fdr = q, knockoffs = knock.gen,
                            statistic = stat.glmnet_lambdasmax)
  w <- result$statistic; names(w) <- colnames(X)
  return(w)
}

# --- Construct X (mutation indicator matrix) and Y (log-fold-resistance) ---

pos_start  <- which(names(gene_df) == 'P1')
pos_cols   <- seq.int(pos_start, ncol(gene_df))
valid_rows <- grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[, pos_cols])
gene_df    <- gene_df[valid_rows, ]
muts <- c(LETTERS, 'i', 'd')
X_hiv <- outer(muts, as.matrix(gene_df[, pos_cols]), Vectorize(grepl))
X_hiv <- aperm(X_hiv, c(2, 3, 1))
dimnames(X_hiv)[[3]] <- muts
X_hiv <- t(apply(X_hiv, 1, flatten_matrix))
mode(X_hiv) <- 'numeric'
X_hiv <- X_hiv[, colSums(X_hiv) != 0]
Y_hiv <- gene_df[, 4:(pos_start - 1)]

# --- Run knockoffs for all drugs ---

set.seed(1)
fdr <- 0.20
results.w <- lapply(Y_hiv, function(y) knockoff_fun(X_hiv, y, fdr))

# Per-drug outlier removal
results.w$RTV <- results.w$RTV[results.w$RTV > -300]
results.w$ATV <- results.w$ATV[results.w$ATV > -150]
results.w$NFV <- results.w$NFV[results.w$NFV > -200]

all.drug.names <- c('APV', 'ATV', 'IDV', 'LPV', 'NFV', 'RTV', 'SQV')

# =============================================================================
# FIGURE 1: HIV barplot — aggregate (left) + thirds (right) for APV/IDV/NFV
# =============================================================================

# --- compute_thirds: split rejection set into thirds by W-statistic strength ---

compute_thirds <- function(w, tsm.positions, q) {
  thresh   <- knockoff_threshold(w, q)
  rej.inds <- which(w >= thresh)
  rej.inds <- rej.inds[order(w[rej.inds], decreasing = TRUE)]
  R        <- length(rej.inds)
  if (R < 3) stop("Fewer than 3 rejections")
  remainder <- R %% 3; n_base <- floor(R / 3)
  
  is.false <- function(inds) {
    positions <- get_position(names(w)[inds])
    !(positions %in% tsm.positions)
  }
  
  if (remainder == 0) {
    idx1 <- rej.inds[1:n_base]
    idx2 <- rej.inds[(n_base + 1):(2 * n_base)]
    idx3 <- rej.inds[(2 * n_base + 1):R]
    fd1 <- sum(is.false(idx1)); fd2 <- sum(is.false(idx2)); fd3 <- sum(is.false(idx3))
    td1 <- n_base - fd1; td2 <- n_base - fd2; td3 <- n_base - fd3
    inds1 <- idx1; inds2 <- idx2; inds3 <- idx3
  } else if (remainder == 1) {
    idx1_core <- rej.inds[1:n_base]; b1 <- rej.inds[n_base + 1]
    idx2_core <- rej.inds[(n_base + 2):(2 * n_base)]; b2 <- rej.inds[2 * n_base + 1]
    idx3_core <- rej.inds[(2 * n_base + 2):R]
    b1f <- is.false(b1); b2f <- is.false(b2)
    fd1 <- sum(is.false(idx1_core)) + (1/3) * b1f
    td1 <- n_base - sum(is.false(idx1_core)) + (1/3) * (1 - b1f)
    fd2 <- sum(is.false(idx2_core)) + (2/3) * b1f + (2/3) * b2f
    td2 <- length(idx2_core) - sum(is.false(idx2_core)) + (2/3) * (1 - b1f) + (2/3) * (1 - b2f)
    fd3 <- sum(is.false(idx3_core)) + (1/3) * b2f
    td3 <- length(idx3_core) - sum(is.false(idx3_core)) + (1/3) * (1 - b2f)
    inds1 <- idx1_core; inds2 <- c(idx2_core, b1, b2); inds3 <- idx3_core
  } else {
    idx1_core <- rej.inds[1:n_base]; b1 <- rej.inds[n_base + 1]
    idx2_core <- rej.inds[(n_base + 2):(2 * n_base + 1)]; b2 <- rej.inds[2 * n_base + 2]
    idx3_core <- rej.inds[(2 * n_base + 3):R]
    b1f <- is.false(b1); b2f <- is.false(b2)
    fd1 <- sum(is.false(idx1_core)) + (2/3) * b1f
    td1 <- n_base - sum(is.false(idx1_core)) + (2/3) * (1 - b1f)
    fd2 <- sum(is.false(idx2_core)) + (1/3) * b1f + (1/3) * b2f
    td2 <- length(idx2_core) - sum(is.false(idx2_core)) + (1/3) * (1 - b1f) + (1/3) * (1 - b2f)
    fd3 <- sum(is.false(idx3_core)) + (2/3) * b2f
    td3 <- length(idx3_core) - sum(is.false(idx3_core)) + (2/3) * (1 - b2f)
    inds1 <- c(idx1_core, b1); inds2 <- idx2_core; inds3 <- c(idx3_core, b2)
  }
  
  make_third <- function(inds, fd, td, w.inds)
    list(inds = inds, fdp = fd / (td + fd), w.range = range(w[w.inds]), n = td + fd)
  
  list(w.threshold = thresh, rej.inds = rej.inds,
       thirds = list(
         make_third(inds1, fd1, td1, rej.inds[1:n_base]),
         make_third(inds2, fd2, td2,
                    rej.inds[(n_base + 1):(2 * n_base + min(remainder, 1) + (remainder == 2))]),
         make_third(inds3, fd3, td3,
                    rej.inds[(2 * n_base + remainder - (remainder > 0) + 1):R])
       ))
}

thirds.all <- lapply(all.drug.names, function(drug)
  compute_thirds(results.w[[drug]], tsm.positions = tsm_df$Position, q = fdr))
names(thirds.all) <- all.drug.names

# --- Aggregate counts and sign-based FDP estimates for selected drugs ---

selected.drugs <- c('APV', 'IDV', 'NFV')
n.drugs <- length(selected.drugs)

agg.counts <- matrix(0, nrow = 2, ncol = n.drugs)
agg.fdp    <- numeric(n.drugs)
agg.sign   <- numeric(n.drugs)
for (i in seq_along(selected.drugs)) {
  drug <- selected.drugs[i]
  w    <- results.w[[drug]]
  thresh <- knockoff_threshold(w, fdr)
  rej    <- w[w >= thresh]
  positions <- get_position(names(rej))
  n.false <- sum(!(positions %in% tsm_df$Position))
  n.total <- length(rej)
  agg.counts[1, i] <- n.total - n.false
  agg.counts[2, i] <- n.false
  agg.fdp[i]  <- n.false / max(1, n.total)
  agg.sign[i] <- sum(w <= -thresh) / max(1, sum(w >= thresh))
}

# Thirds counts and sign-based FDP per third
sel.fdp    <- numeric(n.drugs * 3)
sel.sign   <- numeric(n.drugs * 3)
sel.counts <- matrix(0, nrow = 2, ncol = n.drugs * 3)
for (i in seq_along(selected.drugs)) {
  drug <- selected.drugs[i]
  w    <- results.w[[drug]]
  thirds <- thirds.all[[drug]]$thirds
  for (j in 1:3) {
    col.idx <- 3 * (i - 1) + j
    third   <- thirds[[j]]
    fd <- third$fdp * third$n; td <- third$n - fd
    sel.counts[1, col.idx] <- td
    sel.counts[2, col.idx] <- fd
    sel.fdp[col.idx] <- third$fdp
    w_lo <- third$w.range[1]; w_hi <- third$w.range[2]
    neg_count <- sum(w >= -w_hi & w <= -w_lo)
    pos_count <- sum(w >=  w_lo & w <=  w_hi)
    sel.sign[col.idx] <- neg_count / max(1, pos_count)
  }
}

# --- Plot ---

pdf("fig1.pdf", width = 14, height = 6.5)
layout(matrix(1:2, 1, 2), widths = c(1, 1))

# Left panel: aggregate
par(mar = c(4.5, 4, 3, 1))
agg.n <- colSums(agg.counts)
offset.agg <- max(agg.n) * 0.04
bar.w.agg  <- 3
half.w.agg <- bar.w.agg / 2
spacing.agg <- rep(0.35, n.drugs)
bp.agg <- barplot(agg.counts,
                  col = c('navy', 'orange'),
                  width = bar.w.agg,
                  ylim = c(0, max(agg.n) + offset.agg * 6),
                  ylab = 'Number of Discoveries',
                  space = spacing.agg, xaxt = 'n',
                  cex.lab = 1.3, cex.axis = 1.2)
mtext('Knockoffs aggregate (tail) analysis', side = 3, line = 0.5, cex = 1.8, col = 'black')
axis(1, at = bp.agg, labels = selected.drugs,
     las = 1, cex.axis = 1.4, tick = FALSE, line = 0.3)
text(x = bp.agg, y = agg.n + offset.agg * 1.5,
     labels = sprintf('%.0f%%', agg.fdp * 100),
     cex = 1.4, col = 'darkorange')
text(x = bp.agg, y = agg.n + offset.agg * 3.5,
     labels = sprintf('%.0f%%', agg.sign * 100),
     cex = 1.4, col = 'forestgreen')
segments(x0 = bp.agg - half.w.agg, x1 = bp.agg + half.w.agg,
         y0 = agg.n * (1 - agg.sign), y1 = agg.n * (1 - agg.sign),
         col = 'forestgreen', lwd = 4)

# Right panel: thirds
par(mar = c(4.5, 1.5, 3, 1))
n.values <- colSums(sel.counts)
offset.thirds <- max(n.values) * 0.04

bar.w.thirds <- 1.5
within.gap  <- 0.1
between.gap <- 0.6
spacing.thirds <- rep(c(between.gap, within.gap, within.gap), n.drugs)

bp <- barplot(sel.counts,
              col = c('navy', 'orange'),
              width = bar.w.thirds,
              ylim = c(0, max(agg.n) + offset.agg * 6),
              ylab = '',
              space = spacing.thirds, xaxt = 'n',
              cex.lab = 1.3, cex.axis = 1.2)
mtext('Local FDP by thirds', side = 3, line = 0.5, cex = 1.8, col = 'black')
bar.positions <- sapply(seq_along(selected.drugs),
                        function(i) mean(bp[(3 * i - 2):(3 * i)]))
axis(1, at = bar.positions, labels = selected.drugs,
     las = 1, cex.axis = 1.4, tick = FALSE, line = 0.3)
segments(x0 = bp - 0.6, x1 = bp + 0.6,
         y0 = n.values * (1 - sel.sign), y1 = n.values * (1 - sel.sign),
         col = 'forestgreen', lwd = 4)
text(x = bp, y = n.values + offset.agg * 1.5,
     labels = sprintf('%.0f%%', sel.fdp * 100),
     cex = 1.4, col = 'darkorange')
text(x = bp, y = n.values + offset.agg * 3.5,
     labels = sprintf('%.0f%%', sel.sign * 100),
     cex = 1.4, col = 'forestgreen')
legend('topright',
       legend = c('True Discoveries', 'False Discoveries',
                  'Sign-based FDP estimate'),
       pch = c(15, 15, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 4),
       col = c('navy', 'orange', 'forestgreen'),
       pt.cex = 1.8, cex = 1.55, bty = 'n',
       inset = c(0.02, 0.1))

dev.off()
cat("Saved: fig1.pdf\n")

# =============================================================================
# FIGURE 2: NFV two-panel illustration
# Left: histogram with blue mirror overlay.
# Right: clar curve (logistic regression estimate) with isotonic step function,
#        histogram dots, TSM-based FDP, and knot locations.
# =============================================================================

w <- results.w$NFV   # uses outlier-filtered NFV W-statistics from HIV data section

# --- Common bin grid (shared by histogram and dot estimates) ---
bin.width  <- 20
w.nz       <- w[w != 0]
max.abs    <- max(abs(w.nz))
pos.breaks <- seq(0, ceiling(max.abs / bin.width) * bin.width, by = bin.width)

# Two histograms with identical breaks — drive both picture and dots
h.pos <- hist( w.nz[w.nz > 0], breaks = pos.breaks, plot = FALSE)
h.neg <- hist(-w.nz[w.nz < 0], breaks = pos.breaks, plot = FALSE)
count1   <- h.pos$counts
count0   <- h.neg$counts
grid     <- h.pos$mids
hist_est <- count0 / pmax(1, count1)

# TSM-based FDP on the same grid
H <- rep(NA, length(names(w.nz)))
for (i in seq_along(names(w.nz))) {
  H[i] <- if (length(setdiff(get_position(names(w.nz)[i]), tsm_df$Position)) > 0) 0 else 1
}
tsm_fdps <- numeric(length(grid))
for (k in seq_along(grid)) {
  inds <- which(w.nz > pos.breaks[k] & w.nz <= pos.breaks[k + 1])
  tsm_fdps[k] <- sum(H[inds] == 0) / max(1, length(inds))
}

clar <- compute.clar(w)

# --- Plot ---

pdf("fig2.pdf", width = 12, height = 4.5)
layout(matrix(1:2, 1, 2), widths = c(1, 1.1))
par(mar = c(4.5, 4.5, 2.5, 1))

# Left panel: bin-aligned histogram with mirror overlay
neg.breaks.flipped <- -rev(pos.breaks)
full.breaks <- c(neg.breaks.flipped[-length(neg.breaks.flipped)], pos.breaks)
h.full <- hist(w.nz, breaks = full.breaks, plot = FALSE)
plot(h.full, xlim = c(-100, 250),
     col = 'grey80', border = 'grey50',
     main = '', xlab = 'W-statistic (NFV)', ylab = 'Frequency',
     cex.lab = 1.3, cex.axis = 1.2)
rect(h.neg$breaks[-length(h.neg$breaks)], 0,
     h.neg$breaks[-1], h.neg$counts,
     col = NA, border = 'blue', lwd = 2)

# Right panel: clar curve with isotonic step, histogram dots, TSM triangles
plot(grid, hist_est, type = 'n',
     xlim = c(0, 250), ylim = c(0, 1.5),
     xlab = 'W-statistic (absolute value)', ylab = 'clar(w)',
     main = '', cex.lab = 1.3, cex.axis = 1.2)
lines(clar$x, clar$gren.y, type = 's', lwd = 2, col = 11)
lines(clar$x, clar$y, lwd = 2.5, col = 'black')
points(grid, tsm_fdps, pch = 17, cex = 1.4, col = 'darkorange')
points(grid, hist_est, pch = 16, cex = 0.9, col = 'blue')
abline(h = 0, lwd = 0.5); abline(v = 0, lwd = 0.5)
abline(v = clar$knots,          lty = 2, col = 'grey50', lwd = 1.5)
abline(v = clar$boundary_knots, lty = 3, col = 'grey30', lwd = 1.2)
legend('topright',
       c('logistic regression estimate', 'isotonic regression estimate',
         'histogram estimate', 'TSM-based FDP'),
       col = c('black', 11, 'blue', 'darkorange'),
       lty = c(1, 1, NA, NA), lwd = c(2.5, 2, NA, NA),
       pch = c(NA, NA, 16, 17), cex = 1.0, bty = 'n')

layout(1)
dev.off()
cat("Saved: fig2.pdf\n")
# =============================================================================
# FIGURE 5: iid Gaussian two-groups model — clar comparison
# Uses constrained spline method (compute.clar).
# =============================================================================

set.seed(123)
pi0 <- 0.9; mu_gauss <- 2.5; sigma_gauss <- 1
m   <- 10^4
H   <- rbinom(m, 1, prob = 1 - pi0)
m0  <- m - sum(H); m1 <- sum(H)
w_gauss <- rep(NA, m)
w_gauss[H == 0] <- rnorm(m0, mean = 0,        sd = sigma_gauss)
w_gauss[H == 1] <- rnorm(m1, mean = mu_gauss, sd = sigma_gauss)
clar.est_gauss  <- compute.clar(w_gauss, c(0.9, 0.7, 0.3, 0.1))
clar.est.kde    <- kde.clar(w_gauss)
clar.est.iso    <- isotonic.clar(w_gauss)
truth_gauss     <- true.clar(w_gauss, pi0, mu_gauss, sigma_gauss)

# --- Plot: two panels in one PDF ---
pdf("fig5.pdf", width = 12, height = 4.5)
layout(matrix(1:2, 1, 2), widths = c(1, 1))
par(mar = c(4.5, 4.5, 2.5, 1))

# Left panel: linear scale
plot(truth_gauss$x, truth_gauss$y, col = "black", type = "l", lwd = std_lwd,
     xlim = c(0, 5), ylim = c(0, 1.1),
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "clar estimate", cex.lab = std_cex.lab)
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
lines(clar.est.iso$x,    clar.est.iso$y,    col = 11,       lwd = std_lwd)
lines(clar.est.kde$x,    clar.est.kde$y,    col = "magenta", lwd = std_lwd)
lines(clar.est_gauss$x,  clar.est_gauss$y,  col = "blue",    lwd = std_lwd)

# Right panel: log scale
plot(truth_gauss$x, log(truth_gauss$y), col = "black", type = "l", lwd = std_lwd,
     xlim = c(0, 5),
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "log of clar estimate", cex.lab = std_cex.lab)
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
lines(clar.est.iso$x,    log(clar.est.iso$y),   col = 11,       lwd = std_lwd)
lines(clar.est.kde$x,    log(clar.est.kde$y),   col = "magenta", lwd = std_lwd)
lines(clar.est_gauss$x,  log(clar.est_gauss$y), col = "blue",    lwd = std_lwd)
legend("bottomleft",
       c("isotonic regression estimate", "f(-w)/f(w)", "kde", "logistic regression estimate"),
       col = c(11, "black", "magenta", "blue"), lty = 1, lwd = std_lwd,
       cex = std_cex.legend/1.1, bty = "n", y.intersp = 1.25,
       inset = c(0.05, 0.02))

layout(1)
dev.off()
cat("Saved: fig5.pdf\n")

# =============================================================================
# FIGURE 6: Numerical illustration — knockoffs histogram + logistic/lfdr/clar overlay
# Uses inline 7th-degree polynomial fit (tractable on ~1M W-statistics). Requires W-and-betas-list-largep.rds and stacked_W_lsm.rds.
# =============================================================================

# Simulation setup (same seed as original)
set.seed(1)
n <- 3000; p <- 1000; k <- 100; amplitude <- 4.5
rho <- 0.25
Sigma <- toeplitz(rho^(0:(p - 1)))
X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
nonzero <- sample(p, k)
beta <- amplitude * (1:p %in% nonzero) / sqrt(n)
y.sample <- function(X) X %*% beta + rnorm(n)
y <- y.sample(X)
result <- knockoff.filter(X, y, knockoffs = create.fixed,
                          statistic = stat.glmnet_lambdasmax)
w_fig1 <- result$statistic
X_k_fig1 <- result$Xk

# Left panel: histogram with null/non-null split
# W-and-betas-list-largep.rds stores list[[1]] = beta, list[[2]] = W matrix (N runs x p)
if (!file.exists("W-and-betas-list-largep.rds")) {
  stop("W-and-betas-list-largep.rds not found. Place it in the working directory to proceed.")
}
out.tmp   <- readRDS("W-and-betas-list-largep.rds")
beta_fig1 <- out.tmp[[1]]
W_fig1    <- out.tmp[[2]]
nulls     <- c(W_fig1[, which(beta_fig1 == 0)])
non.nulls <- c(W_fig1[, which(beta_fig1 != 0)])

pdf("fig6.pdf", width = 12, height = 4.5)
layout(matrix(1:2, 1, 2), widths = c(1, 1.1))
par(mar = c(4.5, 4.5, 2.5, 1))

# Left panel: histogram with null/non-null split
hist(non.nulls, breaks = 200, col = "cyan", xlim = c(-400, 500),
     ylim = c(0, 31000), main = "", xlab = "W-statistic", cex.lab = std_cex.lab)
hist(nulls, breaks = 200, col = "blue", add = TRUE)
legend("topleft", c("null", "non-null"), col = c("blue", "cyan"),
       lwd = 3, lty = 1, cex = std_cex.legend, bty = "n")

# Right panel: logistic regression clar + true clar + lfdr
# Both the logistic fit AND the true lfdr/clar all come from W-and-betas-list-largep.rds[[2]]
w_all_fig1 <- c(W_fig1)   # flatten N_runs x p matrix column-major
df_fig1 <- data.frame(y = 1 - (sign(w_all_fig1) + 1) / 2, x = abs(w_all_fig1))
df_fig1$x2 <- df_fig1$x^2; df_fig1$x3 <- df_fig1$x^3; df_fig1$x4 <- df_fig1$x^4
df_fig1$x5 <- df_fig1$x^5; df_fig1$x6 <- df_fig1$x^6; df_fig1$x7 <- df_fig1$x^7
glmod_fig1    <- glm(y ~ x + x2 + x3 + x4 + x5 + x6 + x7, family = binomial, data = df_fig1)
p.hat_fig1    <- glmod_fig1$fitted.values
clar.hat_fig1 <- p.hat_fig1 / (1 - p.hat_fig1)
x.inds_fig1   <- order(df_fig1$x)

# True clar and lfdr: computed from the same W_fig1 matrix row by row
fdp <- function(selected) sum(beta_fig1[selected] == 0) / max(1, length(selected))
nsr <- function(selected, W) sum(W[selected] < 0) / max(1, length(selected))
N_runs_fig1 <- nrow(W_fig1)
K    <- 100
ws_fig1 <- seq(0, 1.2 * max(abs(w_all_fig1)), length = K + 1)
lfdr.lsm     <- matrix(NA, nrow = N_runs_fig1, ncol = K)
clar.pop.lsm <- matrix(NA, nrow = N_runs_fig1, ncol = K)
for (j in 1:N_runs_fig1) {
  w.lsm <- W_fig1[j, ]
  for (i in 1:K) {
    inds_bin <- which((w.lsm < ws_fig1[i + 1]) & (w.lsm > ws_fig1[i]))
    lfdr.lsm[j, i]     <- fdp(inds_bin)
    inds_abs <- which((abs(w.lsm) < ws_fig1[i + 1]) & (abs(w.lsm) > ws_fig1[i]))
    clar.pop.lsm[j, i] <- nsr(inds_abs, w.lsm)
  }
}
mean.lfdr.lsm <- colMeans(lfdr.lsm)
mean.clar.lsm <- colMeans(clar.pop.lsm)
clar.pop_fig1 <- mean.clar.lsm / (1 - mean.clar.lsm)
xaxis_fig1    <- (ws_fig1[2:(K + 1)] + ws_fig1[1:K]) / 2

plot(df_fig1$x[x.inds_fig1], clar.hat_fig1[x.inds_fig1],
     type = "l", col = "blue", lwd = std_lwd,
     ylim = c(0, 1.2), xlim = c(0, 500),
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "lfdr quantity", cex.lab = std_cex.lab)
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
lines(xaxis_fig1, mean.lfdr.lsm, col = "black",     lwd = std_lwd)
lines(xaxis_fig1, clar.pop_fig1, col = "darkorange", lwd = std_lwd, lty = "dashed")
legend("topright", c("logistic regression", "lfdr",
                     "clar"),
       col = c("blue", "black", "darkorange"),
       lty = c(1, 1, 2), lwd = std_lwd,
       cex = std_cex.legend, bty = "n", y.intersp = 1.25)

layout(1)
dev.off()

# =============================================================================
# FIGURE 7: Knockoffs W-statistics — clar estimate with true clar and lfdr
# Uses constrained spline method (compute.clar). Requires stacked_W_lsm.rds.
# [SLOW: 1000-run averaging loop, loads from stacked_W_lsm.rds]
# =============================================================================

set.seed(1)
n <- 3000; p <- 1000; k <- 100; amplitude <- 4.5
rho <- 0.25
Sigma <- toeplitz(rho^(0:(p - 1)))
X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
nonzero <- sample(p, k)
beta <- amplitude * (1:p %in% nonzero) / sqrt(n)
y.sample <- function(X) X %*% beta + rnorm(n)
y <- y.sample(X)
result.lsm <- knockoff.filter(X, y, knockoffs = create.fixed,
                              statistic = stat.glmnet_lambdasmax)
w_kf   <- result.lsm$statistic
X_k_kf <- result.lsm$Xk

clar.est_kf <- compute.clar(w_kf)

# True clar and lfdr from 1000 runs (load from disk, or regenerate if missing)
if (file.exists("stacked_W_lsm.rds")) {
  W_lsm <- readRDS("stacked_W_lsm.rds")
} else {
  message("stacked_W_lsm.rds not found — running 1000 knockoff iterations. This will take ~2 hours.")
  N_runs <- 1000
  W_lsm  <- matrix(NA, nrow = N_runs, ncol = p)
  for (j in 1:N_runs) {
    if (j %% 100 == 0) message("  iteration ", j, " of ", N_runs)
    W_lsm[j, ] <- stat.glmnet_lambdasmax(X, X_k_kf, y.sample(X))
  }
  saveRDS(W_lsm, "stacked_W_lsm.rds")
}
K <- 100
ws_kf <- seq(0, 1.2 * max(abs(w_kf)), length = K + 1)
N_runs <- nrow(W_lsm)
lfdr.lsm_kf    <- matrix(NA, nrow = N_runs, ncol = K)
clar.pop.lsm_kf <- matrix(NA, nrow = N_runs, ncol = K)
fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
nsr <- function(selected, W) sum(W[selected] < 0) / max(1, length(selected))
for (j in 1:N_runs) {
  w.lsm <- W_lsm[j, ]
  for (i in 1:K) {
    inds_bin <- which((w.lsm < ws_kf[i + 1]) & (w.lsm > ws_kf[i]))
    lfdr.lsm_kf[j, i]     <- fdp(inds_bin)
    inds_abs <- which((abs(w.lsm) < ws_kf[i + 1]) & (abs(w.lsm) > ws_kf[i]))
    clar.pop.lsm_kf[j, i] <- nsr(inds_abs, w.lsm)
  }
}
mean.lfdr_kf <- colMeans(lfdr.lsm_kf)
mean.clar_kf <- colMeans(clar.pop.lsm_kf)
clar.pop_kf  <- mean.clar_kf / (1 - mean.clar_kf)
xaxis_kf     <- (ws_kf[2:(K + 1)] + ws_kf[1:K]) / 2

pdf("fig7.pdf", width = 6, height = 4)
par_std()
plot(clar.est_kf$x, clar.est_kf$y, type = "l", col = "blue", lwd = std_lwd,
     ylim = c(-0.1, 1.3),
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "clar estimate", cex.lab = std_cex.lab)
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
lines(xaxis_kf, mean.lfdr_kf, col = "black",     lwd = std_lwd)
lines(xaxis_kf, clar.pop_kf,  col = "darkorange", lwd = std_lwd, lty = "dashed")
legend("topright", c("logistic regression estimate", "true lfdr (1000 runs average)",
                     "true clar (1000 runs average)"),
       col = c("blue", "black", "darkorange"),
       lty = c(1, 1, 2), lwd = std_lwd,
       cex = std_cex.legend/1.2, bty = "n", y.intersp = 1.25)
dev.off()
# =============================================================================
# FIGURE 8: iid Gaussian — bootstrap standard errors
# Uses constrained spline method (compute.clar).
# [SLOW: N=200 bootstrap iterations]
# =============================================================================

# Uses w_gauss and clar.est_gauss from Figure 3 section above
N_boot <- 200; K_boot <- 50
upper.rng_gauss <- 1.2 * max(abs(w_gauss))
ws_boot <- seq(0, upper.rng_gauss, length = K_boot)
new.clars_gauss <- matrix(NA, nrow = N_boot, ncol = K_boot)
for (j in 1:N_boot) {
  w.boot <- sample(w_gauss, replace = TRUE)
  clar.hat.boot <- compute.clar(w.boot)
  for (i in 1:K_boot) {
    ind <- which.min(abs(ws_boot[i] - clar.hat.boot$x))
    new.clars_gauss[j, i] <- clar.hat.boot$y[ind]
  }
}
boot.se_gauss <- apply(new.clars_gauss, 2, sd)

# Original estimate on a grid
og.clar_gauss <- sapply(ws_boot, function(wi) {
  clar.est_gauss$y[which.min(abs(wi - clar.est_gauss$x))]
})

# True SE from N=200 fresh draws [SLOW]
N_true <- 200
true.clars_gauss <- matrix(NA, nrow = N_true, ncol = K_boot)
for (j in 1:N_true) {
  H_new <- rbinom(m, 1, prob = 1 - pi0)
  w_new <- ifelse(H_new == 0,
                  rnorm(m, 0, sigma_gauss),
                  rnorm(m, mu_gauss, sigma_gauss))
  clar.hat.new <- compute.clar(w_new)
  for (i in 1:K_boot) {
    ind <- which.min(abs(ws_boot[i] - clar.hat.new$x))
    true.clars_gauss[j, i] <- clar.hat.new$y[ind]
  }
}
true.se_gauss <- apply(true.clars_gauss, 2, sd)

# --- Plot: two panels in one PDF ---
pdf("fig8.pdf", width = 12, height = 4.5)
layout(matrix(1:2, 1, 2), widths = c(1, 1))
par(mar = c(4.5, 4.5, 2.5, 1))

# Left: bootstrap SE
plot(clar.est_gauss$x, clar.est_gauss$y, type = "l", col = "blue", lwd = std_lwd,
     xlim = c(0, 5), ylim = c(0, 1.2),
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "clar estimate", cex.lab = std_cex.lab)
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
points(ws_boot, og.clar_gauss, pch = 16, col = "blue", cex = std_cex_points)
arrows(ws_boot, og.clar_gauss - 2 * boot.se_gauss,
       ws_boot, og.clar_gauss + 2 * boot.se_gauss,
       angle = 90, code = 3, length = 0.05)
legend("topright", c("logistic regression estimate", "+/- 2 bootstrap SE"),
       col = c("blue", "black"), lty = 1, lwd = std_lwd,
       cex = std_cex.legend/1.1, bty = "n", y.intersp = 1.25,
       inset = c(0.02, 0.02))

# Right: true SE
plot(ws_boot, og.clar_gauss, type = "l", col = "blue", lwd = std_lwd,
     xlim = c(0, 5), ylim = c(0, 1.2),
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "clar estimate", cex.lab = std_cex.lab)
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
points(ws_boot, og.clar_gauss, pch = 16, col = "blue", cex = std_cex_points)
arrows(ws_boot, og.clar_gauss - 2 * true.se_gauss,
       ws_boot, og.clar_gauss + 2 * true.se_gauss,
       angle = 90, code = 3, length = 0.05)
legend("topright", c("logistic regression estimate", "+/- 2 SE"),
       col = c("blue", "black"), lty = 1, lwd = std_lwd,
       cex = std_cex.legend/1.1, bty = "n", y.intersp = 1.25,
       inset = c(0.02, 0.02))

layout(1)
dev.off()
cat("Saved: fig8.pdf\n")

# FIGURE 9: Knockoffs — parametric bootstrap standard errors
# Uses constrained spline method (compute.clar). Requires stacked_W_lsm.rds.
# [SLOW: parametric bootstrap]
# =============================================================================
# Uses w_kf, X, y, X_k_kf, beta from Figure 7 section above

# Parametric bootstrap setup (Lasso + OLS)
cv.lasso <- cv.glmnet(X, y)
beta.hat.Lasso <- as.matrix(coef(cv.lasso))
X.ols <- cbind(rep(1, n), X)[, which(beta.hat.Lasso != 0)]
lmod_ols <- lm(y ~ X.ols)
beta.hat.OLS <- lmod_ols$coefficients[-2]
new.beta <- beta.hat.Lasso
new.beta[beta.hat.Lasso != 0] <- beta.hat.OLS
y.sample.param.boot <- function(X) cbind(rep(1, n), X) %*% new.beta + rnorm(n)

# Common grid anchored to observed W-statistics — used by both panels
B_pb <- 20; K_pb <- 50
upper.rng_common <- 1.2 * max(abs(w_kf))
ws_grid <- seq(0, upper.rng_common, length = K_pb)

# Parametric bootstrap draws
w.boots_kf <- matrix(NA, nrow = B_pb, ncol = p)
for (i in 1:B_pb) {
  new.y <- y.sample.param.boot(X)
  w.boots_kf[i, ] <- stat.glmnet_lambdasmax(X, X_k_kf, new.y)
}
new.clars_kf <- matrix(NA, nrow = B_pb, ncol = K_pb)
for (j in 1:B_pb) {
  clar.boot <- compute.clar(w.boots_kf[j, ])
  for (i in 1:K_pb) {
    new.clars_kf[j, i] <- clar.boot$y[which.min(abs(ws_grid[i] - clar.boot$x))]
  }
}
boot.se_kf <- apply(new.clars_kf, 2, sd)

# og.clar: original single realization evaluated on ws_grid
og.clar_kf <- sapply(ws_grid, function(wi) {
  clar.est_kf$y[which.min(abs(wi - clar.est_kf$x))]
})

# True SE from stacked_W_lsm.rds [SLOW]
if (file.exists("stacked_W_lsm.rds")) {
  W_lsm <- readRDS("stacked_W_lsm.rds")
} else {
  message("stacked_W_lsm.rds not found — running 1000 knockoff iterations. This will take ~2 hours.")
  N_runs <- 1000
  W_lsm  <- matrix(NA, nrow = N_runs, ncol = p)
  for (j in 1:N_runs) {
    if (j %% 100 == 0) message("  iteration ", j, " of ", N_runs)
    W_lsm[j, ] <- stat.glmnet_lambdasmax(X, X_k_kf, y.sample(X))
  }
  saveRDS(W_lsm, "stacked_W_lsm.rds")
}
true.clars_kf <- matrix(NA, nrow = nrow(W_lsm), ncol = K_pb)
for (j in 1:nrow(W_lsm)) {
  tmp <- compute.clar(W_lsm[j, ])
  for (i in 1:K_pb) {
    true.clars_kf[j, i] <- tmp$y[which.min(abs(ws_grid[i] - tmp$x))]
  }
}
true.se_kf   <- apply(true.clars_kf, 2, sd)
true.mean_kf <- apply(true.clars_kf, 2, mean)

# --- Plot: two panels in one PDF ---
pdf("fig9.pdf", width = 12, height = 4.5)
layout(matrix(1:2, 1, 2), widths = c(1, 1))
par(mar = c(4.5, 4.5, 2.5, 1))

# Left: parametric bootstrap SE
plot(clar.est_kf$x, clar.est_kf$y, type = "l", col = "blue", lwd = std_lwd,
     ylim = c(-0.1, 1.3), xlim = c(0, upper.rng_common), cex.lab = std_cex.lab,
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "clar estimate")
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
points(ws_grid, og.clar_kf, pch = 16, col = "blue", cex = std_cex_points)
arrows(ws_grid, og.clar_kf - 2 * boot.se_kf, ws_grid, og.clar_kf + 2 * boot.se_kf,
       angle = 90, code = 3, length = 0.05)
legend("topright", c("logistic regression estimate", "+/- 2 bootstrap SE"),
       col = c("blue", "black"), lty = 1, lwd = std_lwd,
       cex = std_cex.legend/1.1, bty = "n", y.intersp = 1.25,
       inset = c(0.02, 0.02))

# Right: true SE
plot(ws_grid, og.clar_kf, type = "l", col = "blue", lwd = std_lwd,
     ylim = c(-0.1, 1.3), xlim = c(0, upper.rng_common), cex.lab = std_cex.lab,
     xlab = expression("absolute value of " * W * "-statistic"),
     ylab = "clar estimate")
points(ws_grid, og.clar_kf, pch = 16, col = "blue", cex = std_cex_points)
abline(h = 0, lty = "dashed", lwd = std_lwd.dashed)
abline(v = 0, lty = "dashed", lwd = std_lwd.dashed)
arrows(ws_grid, og.clar_kf - 2 * true.se_kf, ws_grid, og.clar_kf + 2 * true.se_kf,
       angle = 90, code = 3, length = 0.05)
legend("topright", c("logistic regression estimate", "+/- 2 SE"),
       col = c("blue", "black"), lty = 1, lwd = std_lwd,
       cex = std_cex.legend/1.1, bty = "n", y.intersp = 1.25,
       inset = c(0.02, 0.02))

layout(1)
dev.off()
cat("Saved: fig9.pdf\n")


# =============================================================================
# FIGURE 10: Proteomics — clar estimation for BDHI-8 data
# Uses older natural cubic spline (compute.clar.poly). Requires bdhi8.csv.
# =============================================================================

protein.dat <- read.csv("bdhi8.csv")
w_prot <- protein.dat$z
w_prot <- w_prot[w_prot != 0]   # remove exact zeros
w_prot <- w_prot[w_prot > -5]   # remove outliers

# --- Bin grid (shared by both panels) ---
bin.width  <- 0.2
max.abs    <- max(abs(w_prot))
pos.breaks <- seq(0, ceiling(max.abs / bin.width) * bin.width, by = bin.width)
neg.breaks <- -rev(pos.breaks)
full.breaks <- c(neg.breaks[-length(neg.breaks)], pos.breaks)

# Use w.nz for everything to match the estimate computation
w_prot_nz <- w_prot   # already zero-stripped above
h.pos <- hist( w_prot_nz[w_prot_nz > 0], breaks = pos.breaks, plot = FALSE)
h.neg <- hist(-w_prot_nz[w_prot_nz < 0], breaks = pos.breaks, plot = FALSE)

# Build h.full by combining the two — guarantees boundary consistency
h.full <- list(
  breaks  = full.breaks,
  counts  = c(rev(h.neg$counts), h.pos$counts),
  mids    = (full.breaks[-1] + full.breaks[-length(full.breaks)]) / 2,
  xname   = "w_prot",
  equidist = TRUE
)
class(h.full) <- "histogram"
count1_prot   <- h.pos$counts
count0_prot   <- h.neg$counts
grid_prot     <- h.pos$mids
hist_est_prot <- count0_prot / pmax(1, count1_prot)

# --- Main fit (polynomial method) ---
clar_prot <- compute.clar.poly(w_prot)

# --- Bootstrap CI ---
N_boot_prot    <- 200
K_eval         <- 200
upper.rng_prot <- 1.2 * max(abs(w_prot))
ws_prot        <- seq(0, upper.rng_prot, length = K_eval)

cat("Running bootstrap (B =", N_boot_prot, ")...\n")
boot.clars_prot <- matrix(NA, nrow = N_boot_prot, ncol = K_eval)
for (j in 1:N_boot_prot) {
  if (j %% 50 == 0) cat("  replicate", j, "of", N_boot_prot, "\n")
  w.boot   <- sample(w_prot, replace = TRUE)
  clar.boot <- tryCatch(compute.clar.poly(w.boot), error = function(e) NULL)
  if (!is.null(clar.boot)) {
    boot.clars_prot[j, ] <- approx(clar.boot$x, clar.boot$y,
                                   xout = ws_prot, rule = 1)$y
  }
}
# boot.lower_prot <- apply(boot.clars_prot, 2, quantile, 0.025, na.rm = TRUE)
# boot.upper_prot <- apply(boot.clars_prot, 2, quantile, 0.975, na.rm = TRUE)
# cat("Bootstrap done.\n")
boot.se_prot    <- apply(boot.clars_prot, 2, sd, na.rm = TRUE)
og.clar_prot    <- approx(clar_prot$x, clar_prot$y, xout = ws_prot, rule = 1)$y
boot.lower_prot <- og.clar_prot - 2 * boot.se_prot
boot.upper_prot <- og.clar_prot + 2 * boot.se_prot
cat("Bootstrap done.\n")

# --- Plot ---
pdf("fig10.pdf", width = 11, height = 4.5)
layout(matrix(1:2, 1, 2), widths = c(1, 1.1))
par(mar = c(4.5, 4.5, 2.5, 1))

# Left panel: histogram with blue mirror overlay
plot(h.full, xlim = c(-3,5),
     col = 'grey80', border = 'grey50',
     main = '', xlab = 'log competition ratio', ylab = '',
     cex.lab = 1.3, cex.axis = 1.2)
rect(h.neg$breaks[-length(h.neg$breaks)], 0,
     h.neg$breaks[-1], h.neg$counts,
     col = NA, border = 'blue', lwd = 1.5)

# Right panel: clar curve with bootstrap CI
plot(clar_prot$x, clar_prot$y, type = "l", lwd = 2, ylim = c(0, 1.1),
     xlim = c(0, upper.rng_prot/1.2),
     xlab = "log ratio (absolute value)", ylab = "clar estimate",
     cex.lab = 1.3, cex.axis = 1.2)
abline(h = 0, lwd = 0.5); abline(v = 0, lwd = 0.5)
points(grid_prot, hist_est_prot, pch = 16, cex = 0.75, col = "blue")
lines(ws_prot, boot.lower_prot, lty = "dashed")
lines(ws_prot, boot.upper_prot, lty = "dashed")
legend("topright",
       c("logistic regression estimate", "histogram estimate", 
         expression(phantom(x) %+-% '2 bootstrap SE')),
       col = c("black", "blue", "black"),
       lty = c(1, NA, 2), lwd = c(2, NA, 1), y.intersp = 1.25,
       pch = c(NA, 16, NA), cex = 1.2, bty = "n")

layout(1)
dev.off()
cat("Saved: fig10.pdf\n")

# figure out clar.hat(w=2)
ind <- min(which(clar_prot$x >= 2))
clar.at.2 <- clar_prot$y[ind]
clar.at.2   # 0.5027

ind2 <- min(which(ws_prot >= 2))
se.at.2 <- boot.se_prot[ind2]

CI.left  <- clar.at.2 - 2 * se.at.2
CI.right <- clar.at.2 + 2 * se.at.2
c(CI.left, CI.right)

# =============================================================================
# FIGURE 11: Pooled HIV — clar estimate with bootstrap SE bands
# Uses constrained spline method (compute.clar).
# Requires the HIV data loaded earlier (results.w, tsm_df, all.drug.names).
# =============================================================================

# --- Rescale and pool ---
maxes <- sapply(all.drug.names, function(drug) max(abs(results.w[[drug]])))
w_pooled <- unlist(lapply(all.drug.names, function(drug)
  results.w[[drug]] / maxes[drug]))
clar_pooled <- compute.clar(w_pooled)

# --- Build pooled W with TSM labels ---
w_pooled_list <- lapply(all.drug.names, function(drug) {
  w <- results.w[[drug]]
  positions <- get_position(names(w))
  is_true <- positions %in% tsm_df$Position
  data.frame(w_scaled = w / maxes[drug], is_true = is_true)
})
w_pooled_df <- do.call(rbind, w_pooled_list)

# --- Bin grid (used by both panels) ---
w_pooled_nz <- w_pooled[w_pooled != 0]
bin.width   <- 0.025
max.abs     <- max(abs(w_pooled_nz))
pos.breaks  <- seq(0, ceiling(max.abs / bin.width) * bin.width, by = bin.width)
neg.breaks  <- -rev(pos.breaks)
full.breaks <- c(neg.breaks[-length(neg.breaks)], pos.breaks)

h.full <- hist(w_pooled_nz, breaks = full.breaks, plot = FALSE)
h.neg  <- hist(-w_pooled_nz[w_pooled_nz < 0], breaks = pos.breaks, plot = FALSE)
h.pos  <- hist( w_pooled_nz[w_pooled_nz > 0], breaks = pos.breaks, plot = FALSE)

count1   <- h.pos$counts
count0   <- h.neg$counts
grid.mid <- h.pos$mids
hist.est <- count0 / pmax(1, count1)

# TSM-based FDP on the same bins
tsm_fdps <- numeric(length(grid.mid))
for (k in seq_along(grid.mid)) {
  inds <- which(w_pooled_df$w_scaled >  pos.breaks[k] &
                  w_pooled_df$w_scaled <= pos.breaks[k + 1])
  tsm_fdps[k] <- sum(!w_pooled_df$is_true[inds]) / max(1, length(inds))
}

# --- Bootstrap SE ---
B <- 200
w.eval <- seq(0, 0.5, length.out = 50)
og.clar.grid <- sapply(w.eval, function(wi)
  clar_pooled$y[which.min(abs(wi - clar_pooled$x))])

cat("Running bootstrap (B =", B, ")...\n")
boot.clar <- matrix(NA, nrow = B, ncol = length(w.eval))
for (b in 1:B) {
  if (b %% 50 == 0) cat("  replicate", b, "of", B, "\n")
  w.boot   <- sample(w_pooled, length(w_pooled), replace = TRUE)
  clar.boot <- tryCatch(compute.clar(w.boot), error = function(e) NULL)
  if (!is.null(clar.boot)) {
    boot.clar[b, ] <- sapply(w.eval, function(wi)
      clar.boot$y[which.min(abs(wi - clar.boot$x))])
  }
}
boot.se <- apply(boot.clar, 2, sd, na.rm = TRUE)
cat("Bootstrap done.\n")

# --- Plot ---
pdf("fig11.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))

# Left panel: pooled histogram with mirror overlay
plot(h.full, xlim = c(-0.3, 0.5),
     col = 'grey80', border = 'grey50',
     main = '',
     xlab = 'W / max(|W|)', ylab = '',
     cex.lab = 1.25, cex.axis = 1.2, cex.main = 1.3,
     xaxt = 'n')
axis(1, at = round(seq(-0.3, 0.5, by = 0.1), 1), cex.axis = 1.2)
rect(h.neg$breaks[-length(h.neg$breaks)], 0,
     h.neg$breaks[-1], h.neg$counts,
     col = NA, border = 'blue', lwd = 1.5)

# Right panel: clar curve with SE bands
plot(clar_pooled$x, clar_pooled$y, type = 'l', lwd = 2,
     xlab = '|W| / max(|W|)', ylab = 'clar(w)',
     xlim = c(0, 0.5),
     ylim = c(0, max(1.1, max(og.clar.grid + 2 * boot.se, na.rm = TRUE))),
     cex.lab = 1.25, cex.axis = 1.2, cex.main = 1.3)
abline(h = 0, lwd = 0.5); abline(v = 0, lwd = 0.5)
lines(w.eval, og.clar.grid + 2 * boot.se, lty = 2, lwd = 1.5, col = 'grey40')
lines(w.eval, og.clar.grid - 2 * boot.se, lty = 2, lwd = 1.5, col = 'grey40')
points(grid.mid, tsm_fdps, pch = 17, cex = 1.2, col = 'darkorange')
points(grid.mid, hist.est, pch = 16, cex = 0.75, col = 'blue')
legend('topright',
       c('logistic regression estimate', 'histogram estimate', 'TSM-based FDP',
         expression(phantom(x) %+-% '2 bootstrap SE')),
       col = c('black', 'blue', 'darkorange', 'grey40'),
       lty = c(1, NA, NA, 2), lwd = c(2, 1, 1, 1.5),
       pch = c(NA, 16, 17, NA), cex = 1.1, bty = 'n')

par(mfrow = c(1, 1))
dev.off()
cat("Saved: fig11.pdf\n")

# --- Diagnostic: triangles above the curve ---
clar_at_triangles <- approx(clar_pooled$x, clar_pooled$y,
                            xout = grid.mid, rule = 2)$y
n_above <- sum(tsm_fdps > clar_at_triangles)
n_total <- length(tsm_fdps)
cat("Triangles above the curve:", n_above, "of", n_total, "\n")

above_idx <- which(tsm_fdps > clar_at_triangles)
print(data.frame(
  bin_center = grid.mid[above_idx],
  tsm_fdp    = tsm_fdps[above_idx],
  clar_est   = clar_at_triangles[above_idx],
  excess     = tsm_fdps[above_idx] - clar_at_triangles[above_idx]
))
