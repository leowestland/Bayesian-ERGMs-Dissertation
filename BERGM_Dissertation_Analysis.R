# =============================================================================
# MSc dissertation: consolidated R analysis
# =============================================================================
#
# This file consolidates:
#   1. Florentine marriage network analysis (Chapter 4)
#   2. Florentine frequentist/Bayesian GOF analysis (Chapter 5)
#   3. Islamic State in Europe network analysis (Chapter 6)
#
# The final analytical workflow is preserved below in execution order.
# Superseded duplicate figures, abandoned trial specifications, and commented
# insertion instructions have been removed. Package loading is centralised
# here, and Chapter 6 uses ISE_DATA_DIR rather than setwd().
#
# To use a different IS-E data directory, either edit ISE_DATA_DIR below or set
# it before sourcing this file, for example:
#
#   ISE_DATA_DIR <- "/path/to/Nov13-master"
#   source("Consolidated_Dissertation_Analysis.R")
#
# Chapter 6 outputs are written to ISE_OUTPUT_DIR (the current working
# directory by default).

required_packages <- c(
  "igraph", "ergm", "Bergm", "sna", "brms", "posterior", "coda",
  "ggraph", "ggplot2", "ggforce"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# Preserve the original package attachment order. Several igraph calls in the
# IS-E section are namespace-qualified because statnet/sna mask common names.
invisible(lapply(required_packages, library, character.only = TRUE))

if (!exists("ISE_DATA_DIR", inherits = FALSE)) {
  ISE_DATA_DIR <- Sys.getenv(
    "ISE_DATA_DIR",
    unset = "/Users/leowestland/Desktop/Stats MSc/MSc dissertation/Nov13-master"
  )
}

if (!exists("ISE_OUTPUT_DIR", inherits = FALSE)) {
  ISE_OUTPUT_DIR <- getwd()
}

ise_data_path <- function(...) file.path(ISE_DATA_DIR, ...)
ise_output_path <- function(...) file.path(ISE_OUTPUT_DIR, ...)

# =============================================================================
# PART I: FLORENTINE MARRIAGE NETWORK (CHAPTER 4)
# =============================================================================

# =========================================================================
# Florentine marriage network: ERGM and BERGM analysis
# Reproducible version: one master seed, per-block derived seeds
# =========================================================================

# =========================================================================
# Reproducibility
# =========================================================================
#
# Everything stochastic in this script traces back to MASTER_SEED.
#
# block_seed() turns a block's *name* into a seed by hashing it together
# with the master seed. Two consequences worth knowing:
#
#   1. Each block can be re-run on its own and reproduce exactly, without
#      having to re-run everything above it first.
#   2. Adding a new block later does not shift the seeds of existing
#      blocks, so previously generated figures stay identical.
#
# The three RNG streams in play:
#   - R's global stream: rnorm, rbinom, sample, ergm/Bergm MCMC
#   - Stan's own stream: the brm() fit, seeded explicitly below
#   - control.ergm(seed=): NOT used anywhere, deliberately. Passing it
#     overrides the global stream and severs the chain from MASTER_SEED.

MASTER_SEED <- 15092026

block_seed <- function(label, master = MASTER_SEED) {
  h <- 5381
  for (ch in utf8ToInt(as.character(label))) {
    h <- (h * 33 + ch) %% 2147483647
  }
  as.integer((h + master * 7919) %% 2147483647)
}

set.seed(MASTER_SEED)

# =========================================================================
# 0. Shared density estimator
# =========================================================================
#
# Every mode in the write-up comes from this one function. If the table and
# the figures call density() separately with different arguments, the mode
# line can end up off the plotted peak.
#
# lower / upper set the evaluation grid, NOT the kernel. R's density() uses
# a Gaussian kernel, so a little mass still leaks past a hard boundary and
# the density is slightly under-estimated there. For theta_W (posterior
# ~2.6 SD above 0) and theta_GWD (~2.8 SD inside +/-1.9) this is negligible,
# but say so in a footnote rather than leaving it implicit.

post_density <- function(draws, lower = NULL, upper = NULL,
                         n = 10000, bw = "nrd0") {
  density(
    draws,
    from = if (is.null(lower)) min(draws) else lower,
    to   = if (is.null(upper)) max(draws) else upper,
    n    = n,
    bw   = bw
  )
}

post_mode <- function(draws, ...) {
  d <- post_density(draws, ...)
  d$x[which.max(d$y)]
}


# For the write-up: record MASTER_SEED and sessionInfo() together.
# A seed pins the randomness, not the package versions.
# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# =========================================================================
# Data
# =========================================================================

data(florentine)

n_nodes <- network.size(flomarriage)

flomarriage %v% "wealth_c" <-
  flomarriage %v% "wealth" - mean(flomarriage %v% "wealth")

flomarriage %v% "wealth_z" <-
  (flomarriage %v% "wealth_c") / sd(flomarriage %v% "wealth_c")

# Adjacency matrix, used by several blocks below
Y <- as.matrix.network.adjacency(flomarriage)

# ========================================
# Plot of the Florentine Network
# ========================================

par(mfrow = c(1, 1))

# Attributes
wealth <- flomarriage %v% "wealth"

deg <- sna::degree(flomarriage, gmode = "graph")
btw <- sna::betweenness(flomarriage, gmode = "graph")

# Rescale wealth for node sizes
node_size <- 1.5 + 4.5 *
  (wealth - min(wealth)) /
  (max(wealth) - min(wealth))

# Normalise betweenness to [0,1]
btw_scaled <- (btw - min(btw)) /
  (max(btw) - min(btw))

# Yellow -> orange -> red -> dark red
btw_pal <- colorRampPalette(
  c("lightyellow", "orange", "orangered", "darkred")
)

btw_cols <- btw_pal(100)

node_col <- btw_cols[
  1 + round(btw_scaled * 99)
]

# Plot

plot(
  flomarriage,
  displaylabels = TRUE,
  label = network.vertex.names(flomarriage),
  
  label.cex = 0.8,
  
  vertex.cex = node_size,
  vertex.col = node_col,
  
  edge.col = "grey40",
  edge.lwd = 1.5,
  
  main = "Florentine Marriage Network"
)

# ========================================
# Edge only ERGM
# ========================================

# Dyad-independent, so no MCMC, but seeded anyway for safety
set.seed(block_seed("ergm_edges"))

summary(flomarriage ~ edges)

flomodel.01 <- ergm(flomarriage ~ edges)

summary(flomodel.01)

# ========================================
# PRIOR DISTRIBUTIONS on edge model
# ========================================

N_dyads  <- 120       # possible edges
m_edges  <- 20        # observed edges

# Priors
priors <- list(
  "Beta(1,1)" = c(1, 1),
  "Beta(2,2)" = c(2, 2),
  "Beta(2,8)" = c(2, 8)
)

prior_cols <- c("black", "blue", "red")
prior_ltys <- c(1, 2, 3)

# Values of p
p_grid <- seq(0, 1, length.out = 1000)

# Two-panel figure
par(mfrow = c(1, 2))
par(mar = c(5.1, 4.1, 4.1, 2.1))

# 1. PRIOR DISTRIBUTIONS

prior_density <- lapply(priors, function(ab) {
  dbeta(p_grid, ab[1], ab[2])
})

plot(
  p_grid,
  prior_density[[1]],
  type = "l",
  lwd = 2,
  col = prior_cols[1],
  lty = prior_ltys[1],
  ylim = c(0, max(unlist(prior_density))),
  xlim = c(-0.1, 1.1),
  xlab = expression(p),
  ylab = "Density",
  main = "Prior distributions"
)

for (i in 2:length(priors)) {
  lines(
    p_grid,
    prior_density[[i]],
    lwd = 2,
    col = prior_cols[i],
    lty = prior_ltys[i]
  )
}

legend(
  "topright",
  legend = names(priors),
  col = prior_cols,
  lty = prior_ltys,
  lwd = 2,
  bty = "n"
)

# 2. POSTERIOR DISTRIBUTIONS

posterior_density <- lapply(priors, function(ab) {
  
  a_post <- ab[1] + m_edges
  b_post <- ab[2] + N_dyads - m_edges
  
  dbeta(p_grid, a_post, b_post)
})

plot(
  p_grid,
  posterior_density[[1]],
  type = "l",
  lwd = 2,
  col = prior_cols[1],
  lty = prior_ltys[1],
  ylim = c(0, max(unlist(posterior_density))),
  xlab = expression(p),
  ylab = "Density",
  main = "Posterior distributions"
)

for (i in 2:length(priors)) {
  lines(
    p_grid,
    posterior_density[[i]],
    lwd = 2,
    col = prior_cols[i],
    lty = prior_ltys[i]
  )
}

# MLE line
segments(
  x0 = m_edges / N_dyads,
  y0 = 0,
  x1 = m_edges / N_dyads,
  y1 = 12,
  col = "black",
  lwd = 2,
  lty = 4
)

legend(
  "topright",
  legend = c(
    "Beta(21,101)",
    "Beta(22,102)",
    "Beta(22,108)",
    "MLE = 1/6"
  ),
  col = c(prior_cols, "black"),
  lty = c(prior_ltys, 4),
  lwd = 2,
  bty = "n"
)

# ========================================
# Bar chart of Florentine Families Wealth
# ========================================

families <- network.vertex.names(flomarriage)

ord <- order(flomarriage %v% "wealth")

par(mfrow = c(1, 1))
par(mar = c(5, 8, 2, 1))

barplot(
  (flomarriage %v% "wealth")[ord],
  names.arg = network.vertex.names(flomarriage)[ord],
  horiz = TRUE,
  las = 1,
  xlab = "Wealth (thousands of lira)",
  main = "Wealth of Florentine Families",
  border = NA
)

# ========================================
# Total wealth vs wealth homophily
# ========================================

wealth_z <- flomarriage %v% "wealth_z"

dyads <- t(combn(seq_along(wealth_z), 2))

wealth_activity <-
  wealth_z[dyads[, 1]] + wealth_z[dyads[, 2]]

wealth_difference <-
  abs(wealth_z[dyads[, 1]] - wealth_z[dyads[, 2]])

cor(
  wealth_activity,
  wealth_difference,
  method = "pearson"
)

r <- cor(wealth_activity, wealth_difference)
VIF <- 1 / (1 - r^2)

c(correlation = r, VIF = VIF)

fit <- ergm(
  flomarriage ~
    edges +
    nodecov("wealth_z") +
    absdiff("wealth_z")
)

cov2cor(vcov(fit))

# ========================================
# Edge and standardised wealth ERGM
# ========================================

set.seed(block_seed("ergm_wealth"))

summary(flomarriage ~ edges + nodecov("wealth_z"))

flomodel.11 <- ergm(flomarriage ~ edges + nodecov("wealth_z"))

summary(flomodel.11)

# ========================================
# Edge and standardised wealth: dyadic Bayesian logistic model
# ========================================

# Standardised wealth
wealth_z <- flomarriage %v% "wealth_z"

# Build dyad-level dataset
dyads <- data.frame(
  i = rep(seq_len(n_nodes), each = n_nodes),
  j = rep(seq_len(n_nodes), times = n_nodes)
)

# Keep only unique undirected dyads
dyads <- subset(dyads, i < j)

# Response: whether a marriage edge exists
dyads$edge <- mapply(
  function(i, j) Y[i, j],
  dyads$i,
  dyads$j
)

# Wealth covariate corresponding to nodecov("wealth_z")
dyads$wealth_sum <- mapply(
  function(i, j) wealth_z[i] + wealth_z[j],
  dyads$i,
  dyads$j
)

# Stan has its own RNG, independent of R's, so the seed must be passed
# explicitly here. set.seed() at the top of the script does nothing for it.
flomodel.12 <- brm(
  edge ~ wealth_sum,
  data = dyads,
  family = bernoulli(link = "logit"),
  
  prior = c(
    prior(normal(-1, 1.25), class = "Intercept"),
    prior(normal(0, 1), class = "b", coef = "wealth_sum")
  ),
  
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = block_seed("brms_dyadic"),
  
  backend = "cmdstanr"
)

summary(flomodel.12)
plot(flomodel.12)

# ========================================
# Edge and standardised wealth prior posterior plots
# ========================================

# Deterministic block (KDE + optim), so no block_seed() needed.
# Requires post_density() and post_mode() to be defined above this point.

brms_draws <- as_draws_df(flomodel.12)

# Posterior draws
brms_theta_E <- brms_draws$b_Intercept
brms_theta_W <- brms_draws$b_wealth_sum

# Frequentist MLEs
mle_dyadic   <- coef(flomodel.11)
mle_E_dyadic <- mle_dyadic["edges"]
mle_W_dyadic <- mle_dyadic["nodecov.wealth_z"]

# Marginal posterior densities from the shared estimator, so these modes
# are computed the same way as those in Table 4.7 and Figure 4.11
post_E <- post_density(brms_theta_E)
post_W <- post_density(brms_theta_W)

mode_E_dyadic <- post_E$x[which.max(post_E$y)]
mode_W_dyadic <- post_W$x[which.max(post_W$y)]

# Exact joint MAP as a check on the KDE modes. The model is
# dyad-independent, so the log posterior is tractable up to a constant.
log_post_dyadic <- function(theta) {
  eta <- theta[1] + theta[2] * dyads$wealth_sum
  sum(
    dyads$edge * plogis(eta, log.p = TRUE) +
      (1 - dyads$edge) * plogis(-eta, log.p = TRUE)
  ) +
    dnorm(theta[1], mean = -1, sd = 1.25, log = TRUE) +
    dnorm(theta[2], mean =  0, sd = 1,    log = TRUE)
}

map_fit <- optim(
  par     = c(0, 0),
  fn      = function(th) -log_post_dyadic(th),
  method  = "BFGS",
  control = list(reltol = 1e-12)
)

stopifnot(map_fit$convergence == 0)

dyadic_point_estimates <- rbind(
  MLE           = unname(c(mle_E_dyadic, mle_W_dyadic)),
  marginal_mode = c(mode_E_dyadic, mode_W_dyadic),
  joint_MAP     = map_fit$par
)
colnames(dyadic_point_estimates) <- c("theta_E", "theta_W")
round(dyadic_point_estimates, 4)

par(
  mfrow = c(1, 2),
  mar = c(5, 4, 4, 2) + 0.1
)

# Left plot: theta_E
x_E <- seq(-6, 4, length.out = 1000)
prior_E <- dnorm(x_E, mean = -1, sd = 1.25)

plot(
  post_E,
  col = "red",
  lwd = 2,
  lty = 1,
  xlim = c(-4, 4),
  ylim = c(0, 3),
  xlab = expression(theta[E]),
  ylab = "Density",
  main = bquote("Prior and posterior of " ~ theta[E]),
  font.main = 2
)

lines(x_E, prior_E, col = "red", lwd = 2, lty = 2)

# MAP
segments(
  x0 = mode_E_dyadic, y0 = 0,
  x1 = mode_E_dyadic, y1 = max(post_E$y),
  col = "red", lwd = 2, lty = 3
)

# MLE
post_height_E <- approx(post_E$x, post_E$y, xout = mle_E_dyadic, rule = 2)$y

segments(
  x0 = mle_E_dyadic, y0 = 0,
  x1 = mle_E_dyadic, y1 = post_height_E,
  col = "black", lwd = 2, lty = 4
)

legend(
  "topleft",
  legend = c(
    sprintf("MLE \u2248 %.3f", mle_E_dyadic),
    sprintf("MAP \u2248 %.3f", mode_E_dyadic),
    "Posterior",
    "Prior"
  ),
  col = c("black", "red", "red", "red"),
  lty = c(4, 3, 1, 2),
  lwd = 2,
  bty = "n"
)

# Right plot: theta_W
x_W <- seq(-4, 4, length.out = 1000)
prior_W <- dnorm(x_W, mean = 0, sd = 1)

plot(
  post_W,
  col = "blue",
  lwd = 2,
  lty = 1,
  xlim = c(-4, 4),
  ylim = c(0, 3),
  xlab = expression(theta[W]),
  ylab = "",
  main = bquote("Prior and posterior of " ~ theta[W])
)

lines(x_W, prior_W, col = "blue", lwd = 2, lty = 2)

# MAP
segments(
  x0 = mode_W_dyadic, y0 = 0,
  x1 = mode_W_dyadic, y1 = max(post_W$y),
  col = "blue", lwd = 2, lty = 3
)

# MLE
post_height_W <- approx(post_W$x, post_W$y, xout = mle_W_dyadic, rule = 2)$y

segments(
  x0 = mle_W_dyadic, y0 = 0,
  x1 = mle_W_dyadic, y1 = post_height_W,
  col = "black", lwd = 2, lty = 4
)

legend(
  "topleft",
  legend = c(
    paste0("MLE \u2248 ", round(mle_W_dyadic, 3)),
    paste0("MAP \u2248 ", round(mode_W_dyadic, 3)),
    "Posterior",
    "Prior"
  ),
  col = c("black", "blue", "blue", "blue"),
  lty = c(4, 3, 1, 2),
  lwd = 2,
  bty = "n"
)


# -------------------------
# Prior predictive: dyadic logistic model
# -------------------------

set.seed(block_seed("prior_pred_dyadic"))

n_sim <- 100000

# Draw parameter values from priors
theta_E_prior <- rnorm(n_sim, mean = -1, sd = 1.25)
theta_W_prior <- rnorm(n_sim, mean = 0, sd = 1)

# Store simulated edge counts
prior_pred_edges_dyadic <- numeric(n_sim)

for (s in 1:n_sim) {
  
  # Edge probability for every dyad
  dyad_p <- plogis(
    theta_E_prior[s] +
      theta_W_prior[s] * dyads$wealth_sum
  )
  
  # Simulate each dyad and count total edges
  prior_pred_edges_dyadic[s] <- sum(
    rbinom(
      n = nrow(dyads),
      size = 1,
      prob = dyad_p
    )
  )
}

par(mfrow = c(1, 1))

hist(
  prior_pred_edges_dyadic,
  breaks = seq(-0.5, 124.5, by = 1),
  probability = TRUE,
  main = "Prior predictive edge count ",
  xlab = "Number of edges",
  ylab = "Proportion of simulated graphs",
  border = "white"
)

abline(
  v = 20,
  col = "black",
  lwd = 2,
  lty = 4
)

legend(
  "topright",
  inset = c(0, -0.1),
  legend = "Observed = 20",
  col = "black",
  lty = 4,
  lwd = 2,
  bty = "n",
  xpd = TRUE
)

# ========================================
# Florentine marriage network empirical degree distribution
# ========================================

n_edges <- network.edgecount(flomarriage)

p_hat <- n_edges / choose(n_nodes, 2)

degree_values <- 0:(n_nodes - 1)

expected_prop <- dbinom(
  degree_values,
  size = n_nodes - 1,
  prob = p_hat
)

observed_degree <- rowSums(Y)

# Include all possible degrees, even those not observed
observed_prop <- as.numeric(
  prop.table(
    table(
      factor(observed_degree, levels = degree_values)
    )
  )
)

plot_data <- rbind(
  expected_prop,
  observed_prop
)

par(
  mar = c(5.5, 5.5, 2, 2),
  mgp = c(0, 1, 0),
  family = "sans",
  las = 1
)

bar_positions <- barplot(
  plot_data,
  beside = TRUE,
  names.arg = degree_values,
  col = c("grey45", "grey75"),
  border = NA,
  ylim = c(0, 0.40),
  axes = FALSE,
  axisnames = FALSE,
  space = c(0.10, 1),
  main = "Degree Distribution of the Florentine Marriage Network",
  cex.main = 1.4,
  font.main = 2
)

# Degree labels
axis(
  side = 1,
  at = colMeans(bar_positions),
  labels = degree_values,
  tick = TRUE,
  line = 0.5,
  cex.axis = 1
)

# Proportion axis
axis(
  side = 2,
  at = seq(0, 0.4, by = 0.1),
  labels = sprintf("%.1f", seq(0, 0.4, by = 0.1)),
  las = 2,
  cex.axis = 1
)

# Axis labels
mtext(
  "Degree",
  side = 1,
  line = 4,
  cex = 1.4
)

mtext(
  "Proportion of families",
  side = 2,
  line = 4,
  cex = 1.4,
  las = 0
)

legend(
  "topright",
  legend = c(
    "Erdős–Rényi expected",
    "Florentine observed"
  ),
  fill = c("grey45", "grey75"),
  border = NA,
  bty = "n",
  cex = 1.25
)

box(bty = "l")

# -------------------------
# Edges and 2-star
# -------------------------

set.seed(block_seed("ergm_2star"))

summary(flomarriage ~ edges + kstar(2))

flomodel.21 <- ergm(flomarriage ~ edges + kstar(2))

summary(flomodel.21)

# Print the numerical MCMC diagnostics without the default plots
mcmc.diagnostics(
  flomodel.21,
  center = FALSE,
  which = "texts"
)

# Extract the MCMC sample
mcmc_sample <- flomodel.21$sample

# Accommodate models stored as an mcmc.list
if (inherits(mcmc_sample, "mcmc.list")) {
  mcmc_sample <- mcmc_sample[[1]]
}

sample_matrix <- as.matrix(mcmc_sample)
iterations    <- as.numeric(time(mcmc_sample))

# Observed statistics: edges and 2-stars
observed_stats <- summary(
  flomarriage ~ edges + kstar(2)
)

# Undo the centring
sample_matrix <- sweep(
  sample_matrix,
  MARGIN = 2,
  STATS = observed_stats,
  FUN = "+"
)

# Labels used in the plots
plot_labels <- c("edge count", "2-star count")

# Two rows: trace plot and density plot for each statistic
par(
  mfrow = c(2, 2),
  mar = c(3.2, 3.2, 2.1, 1),
  mgp = c(2, 0.7, 0),
  cex.axis = 1.3,
  cex.lab = 1.3,
  cex.main = 1.3
)

trace_ylim <- list(
  c(5, 35),    # Edge-count trace
  c(0, 130)    # 2-star-count trace
)

trace_ticks <- list(
  c(10, 20, 30),
  c(0, 50, 100)
)

hist_breaks <- list(
  seq(5, 35, by = 2),     # Edge-count bar width = 2
  seq(0, 140, by = 10)    # 2-star bar width = 10
)

x_limits <- list(
  c(0, 35),               # Edge-count x-axis range
  c(0, 150)               # 2-star x-axis range
)

x_ticks <- list(
  seq(0, 35, by = 5),     # Edge-count tick marks
  seq(0, 150, by = 30)    # 2-star tick marks
)

density_ylim <- list(
  c(0, 0.1),    # Top-right: edge count
  c(0, 0.025)   # Bottom-right: 2-star count
)

density_ticks <- list(
  seq(0, 0.08, by = 0.04),
  seq(0, 0.025, by = 0.01)
)

for (j in seq_len(ncol(sample_matrix))) {
  
  # Trace plots
  plot(
    iterations,
    sample_matrix[, j],
    type = "l",
    xaxt = "n",
    yaxt = "n",
    ylim = trace_ylim[[j]],
    xlab = if (j == 1) "" else "Iterations",
    ylab = "",
    main = paste("Trace of", plot_labels[j]),
    lwd = 1
  )
  
  axis(
    side = 2,
    at = trace_ticks[[j]],
    las = 1,             # horizontal numbers
    cex.axis = 1.3
  )
  
  axis(
    side = 1,
    at = c(50000, 150000, 250000),
    labels = c("50000", "150000", "250000"),
    cex.axis = 1.3
  )
  
  # Smoothed trend line
  lines(
    lowess(iterations, sample_matrix[, j]),
    lwd = 1.5
  )
  
  # Density histogram
  hist(
    sample_matrix[, j],
    breaks = hist_breaks[[j]],
    probability = TRUE,
    col = "grey80",
    border = "black",
    xlim = x_limits[[j]],
    ylim = density_ylim[[j]],
    xaxt = "n",
    yaxt = "n",
    xlab = if (j == 2) "Count" else "",
    ylab = "",
    main = paste("Density of", plot_labels[j])
  )
  
  # Custom x-axis
  axis(
    side = 1,
    at = x_ticks[[j]],
    cex.axis = 1.3
  )
  
  # Custom density y-axis
  axis(
    side = 2,
    at = density_ticks[[j]],
    labels = format(
      density_ticks[[j]],
      nsmall = if (j == 1) 2 else 3
    ),
    las = 1,
    cex.axis = 1.3
  )
}

# -------------------------
# Edges and GWD, alpha free
# -------------------------

set.seed(block_seed("ergm_gwd_free"))

# If this struggles to converge, try:
#   control.ergm(MCMC.burnin = 50000, MCMC.interval = 1024,
#                MCMC.samplesize = 5000, MCMLE.maxit = 50)
# but do NOT add seed = here, it would break the master-seed chain.
flomodel.31 <- ergm(
  flomarriage ~ edges + gwdegree(fixed = FALSE)
)

summary(flomodel.31)

#mcmc.diagnostics(flomodel.31, center = FALSE, which = "texts")

# -------------------------
# Alpha grid search
# -------------------------

set.seed(block_seed("gwd_grid"))

# Candidate decay values
alpha_grid <- seq(0.1, 2, by = 0.1)

# Store fitted models
gwd_fits <- vector("list", length(alpha_grid))
names(gwd_fits) <- as.character(alpha_grid)

# Fit one model for each fixed alpha.
# The loop is seeded as a unit rather than per-iteration, so the whole
# grid reproduces together.
for (i in seq_along(alpha_grid)) {
  
  alpha <- alpha_grid[i]
  
  gwd_fits[[i]] <- tryCatch(
    ergm(
      flomarriage ~ edges + gwdegree(alpha, fixed = TRUE)
    ),
    error = function(e) NULL
  )
}

grid_results <- do.call(
  rbind,
  lapply(seq_along(alpha_grid), function(i) {
    
    fit <- gwd_fits[[i]]
    
    if (is.null(fit)) {
      return(
        data.frame(
          alpha = alpha_grid[i],
          theta_edges = NA,
          theta_gwd = NA,
          logLik = NA,
          AIC = NA,
          BIC = NA
        )
      )
    }
    
    data.frame(
      alpha = alpha_grid[i],
      theta_edges = coef(fit)[1],
      theta_gwd = coef(fit)[2],
      logLik = as.numeric(logLik(fit)),
      AIC = AIC(fit),
      BIC = BIC(fit)
    )
  })
)

grid_results

valid_results <- grid_results[complete.cases(grid_results), ]

best_row <- valid_results[
  which.min(valid_results$AIC),
]

best_alpha <- best_row$alpha
best_alpha

best_fit <- gwd_fits[[as.character(best_alpha)]]

summary(best_fit)
mcmc.diagnostics(best_fit, center = FALSE)

# Reset the plotting device safely. A bare dev.off() errors when only the
# null device is open, which is what happens under Rscript.
if (dev.cur() > 1) dev.off()
par(mfrow = c(1, 1))

plot(
  valid_results$alpha,
  valid_results$AIC,
  type = "b",
  pch = 19,
  xlab = expression(alpha),
  ylab = "AIC",
  main = expression(paste("AIC profile for fixed ", alpha))
)

abline(
  v = best_alpha,
  col = "red",
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = bquote("Minimum at " ~ alpha == .(best_alpha)),
  col = "red",
  lty = 2,
  lwd = 2,
  bty = "n"
)



# ========================================
# BERGM: edges, wealth and GWD
# ========================================

# Shared settings, kept modest so each fit takes a minute or two.

BI  <- 500
MI  <- 2000
AUX <- 1000
NC  <- 3


form.41 <- flomarriage ~ edges + nodecov("wealth_z") +
  gwdegree(0.5, fixed = TRUE)

# Prior hyperparameters, defined once so the log prior, the prior
# predictive and the prior/posterior figures cannot drift apart.
PRIOR_E   <- list(mean = -1.5, sd = 0.3)
PRIOR_W   <- list(shape = 2, rate = 4)
PRIOR_GWD <- list(min = -1.9, max = 1.9)

log_prior <- function(theta) {
  dnorm(theta[1], mean = PRIOR_E$mean, sd = PRIOR_E$sd, log = TRUE) +
    dgamma(theta[2], shape = PRIOR_W$shape, rate = PRIOR_W$rate, log = TRUE) +
    dunif(theta[3], min = PRIOR_GWD$min, max = PRIOR_GWD$max, log = TRUE)
}

set.seed(block_seed("bergm_41"))

flomodel.41 <- bergm(
  form.41,
  log.prior = log_prior,
  burn.in = BI,
  main.iters = MI,
  aux.iters = AUX,
  nchains = NC
)

summary(flomodel.41)

coda::effectiveSize(flomodel.41$Theta)

# Posterior draws, extracted once and reused throughout
bergm_draws <- as.matrix(flomodel.41$Theta)

theta_E_post   <- bergm_draws[, 1]
theta_W_post   <- bergm_draws[, 2]
theta_GWD_post <- bergm_draws[, 3]

par(
  family = "sans",
  cex.axis = 1.25,
  cex.lab  = 1.25,
  cex.main = 1.35,
  font.lab = 1,
  font.main = 2,
  las = 1,
  mgp = c(2.6, 0.8, 0),
  tcl = -0.3,
  lwd = 1.2
)

# plot(flomodel.41, center = FALSE)

plot_bergm_custom <- function(
    x,
    main_title,
    parameter_labels = expression(
      theta[E],
      theta[W],
      theta[GWD]
    )
) {
  
  par(
    mfrow = c(x$dim, 3),
    oma = c(0, 0, 3, 0),
    mar = c(4, 2.5, 0.55, 1)
  )
  
  for (i in seq_len(x$dim)) {
    
    # Posterior density
    plot(
      density(x$Theta[, i]),
      main = "",
      axes = FALSE,
      xlab = parameter_labels[i],
      ylab = "",
      lwd = 2
    )
    
    axis(1)
    axis(2)
    
    # Trace
    coda::traceplot(
      x$Theta[, i],
      type = "l",
      xlab = "Iterations",
      ylab = ""
    )
    
    # Autocorrelation
    coda::autocorr.plot(
      x$Theta[, i],
      auto.layout = FALSE
    )
  }
  
  mtext(
    main_title,
    side = 3,
    outer = TRUE,
    line = 1,
    cex = 1.4,
    font = 2
  )
}

# plot_bergm_custom(
#   flomodel.41,
#   main_title = "MCMC diagnostics for the edge, wealth and GWD model"
# )

plot_bergm_chains <- function(fit, main.iters) {
  
  draws <- as.matrix(fit$Theta)
  
  # Recover the number of chains
  nchains <- nrow(draws) / main.iters
  
  if (nchains != round(nchains)) {
    stop("main.iters does not match the stored draws.")
  }
  
  nchains <- as.integer(round(nchains))
  nparams <- ncol(draws)
  
  # iteration x chain x parameter
  chain_array <- array(
    draws,
    dim = c(main.iters, nchains, nparams)
  )
  
  chain_colours <- hcl.colors(nchains, "Dark 3")
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(
    mfrow = c(3, 3),
    oma = c(0, 0, 3, 0),
    mar = c(4, 2.5, 0.55, 1)
  )
  
  parameter_labels <- expression(
    theta[E] ~ "(edges)",
    theta[W] ~ "(wealth)",
    theta[plain(GWD)] ~ "(GWD, decay = 0.5)"
  )
  
  for (j in seq_len(nparams)) {
    
    # Posterior density
    plot(
      density(draws[, j]),
      main = "",
      xlab = parameter_labels[j],
      ylab = "Density",
      lwd = 2
    )
    
    # Separate trace plots
    matplot(
      seq_len(main.iters),
      chain_array[, , j],
      type = "l",
      lty = 1,
      col = chain_colours,
      xlab = "Iterations",
      ylab = ""
    )
    
    # Chain-specific ACFs
    acf_values <- sapply(seq_len(nchains), function(k) {
      as.vector(
        acf(
          chain_array[, k, j],
          lag.max = 40,
          plot = FALSE
        )$acf
      )
    })
    
    matplot(
      0:40,
      acf_values,
      type = "l",
      lty = 1,
      col = chain_colours,
      ylim = c(-1, 1),
      xlab = "Lag",
      ylab = "Autocorrelation"
    )
    
    abline(h = 0, col = "grey60")
    
    # Legend in the top-right plot
    if (j == 1) {
      legend(
        "bottomleft",
        legend = paste("Chain", seq_len(nchains)),
        col = chain_colours,
        lty = 1,
        lwd = 1,
        cex = 0.7,
        bty = "n"
      )
    }
  }
  
  mtext(
    "MCMC diagnostics for the edge, wealth and GWD model",
    outer = TRUE,
    line = 1,
    cex = 1.3,
    font = 2
  )
}

plot_bergm_chains(
  flomodel.41,
  main.iters = MI
)

  # still negative, but weaker once the prior cuts the ridge

# -------------------------
# Alpha grid search: edges + wealth + GWD (form.41)
# -------------------------

set.seed(block_seed("gwd_grid_41"))

alpha_grid_41 <- seq(0.1, 2, by = 0.1)

gwd_fits_41 <- vector("list", length(alpha_grid_41))
names(gwd_fits_41) <- as.character(alpha_grid_41)

n_warn_41 <- integer(length(alpha_grid_41))

for (i in seq_along(alpha_grid_41)) {
  
  # Put the numeric alpha into the formula itself, so later calls on a
  # stored fit (gof, simulate) can't pick up a stale loop variable
  f_i <- eval(bquote(
    flomarriage ~ edges + nodecov("wealth_z") +
      gwdegree(.(alpha_grid_41[i]), fixed = TRUE)
  ))
  
  # Count warnings, since MCMLE non-convergence warns rather than errors
  warn_i <- 0
  
  gwd_fits_41[[i]] <- withCallingHandlers(
    tryCatch(ergm(f_i), error = function(e) NULL),
    warning = function(w) {
      warn_i <<- warn_i + 1
      invokeRestart("muffleWarning")
    }
  )
  
  n_warn_41[i] <- warn_i
  message("alpha = ", alpha_grid_41[i], " done")
}

extract_row_41 <- function(fit, alpha, n_warn) {
  
  if (is.null(fit)) {
    return(data.frame(
      alpha = alpha, theta_E = NA, theta_W = NA, se_W = NA,
      theta_GWD = NA, logLik = NA, AIC = NA, BIC = NA, warnings = n_warn
    ))
  }
  
  est <- coef(fit)
  se  <- sqrt(diag(vcov(fit)))
  gwd <- grep("gwdeg", names(est))
  
  data.frame(
    alpha     = alpha,
    theta_E   = unname(est["edges"]),
    theta_W   = unname(est["nodecov.wealth_z"]),
    se_W      = unname(se["nodecov.wealth_z"]),
    theta_GWD = unname(est[gwd]),
    logLik    = as.numeric(logLik(fit)),
    AIC       = AIC(fit),
    BIC       = BIC(fit),
    warnings  = n_warn
  )
}

grid_results_41 <- do.call(
  rbind,
  Map(extract_row_41, gwd_fits_41, alpha_grid_41, n_warn_41)
)
rownames(grid_results_41) <- NULL

round(grid_results_41, 3)

valid_41 <- grid_results_41[complete.cases(grid_results_41), ]

best_alpha_41 <- valid_41$alpha[which.min(valid_41$AIC)]
best_alpha_41

best_fit_41 <- gwd_fits_41[[as.character(best_alpha_41)]]

summary(best_fit_41)
mcmc.diagnostics(best_fit_41, center = FALSE, which = "texts")

# -------------------------
# AIC profile and wealth-effect sensitivity
# -------------------------

if (dev.cur() > 1) dev.off()

par(
  mfrow = c(1, 2),
  mar = c(5.1, 4.6, 4.1, 1.5)
)

# Left: AIC profile
plot(
  valid_41$alpha,
  valid_41$AIC,
  type = "b",
  pch = 19,
  xlab = expression(alpha),
  ylab = "AIC",
  main = expression(paste("AIC profile for fixed ", alpha))
)

abline(v = best_alpha_41, col = "red", lty = 2, lwd = 2)
abline(v = 0.5, col = "grey50", lty = 3, lwd = 2)

legend(
  "topright",
  legend = as.expression(list(
    bquote("Minimum at" ~ alpha == .(best_alpha_41)),
    bquote("Used in BERGM:" ~ alpha == 0.5)
  )),
  col = c("red", "grey50"),
  lty = c(2, 3),
  lwd = 2,
  bty = "n"
)

# Right: wealth MLE +/- 2 SE across alpha
lo_W <- valid_41$theta_W - 2 * valid_41$se_W
hi_W <- valid_41$theta_W + 2 * valid_41$se_W

plot(
  valid_41$alpha,
  valid_41$theta_W,
  type = "b",
  pch = 19,
  ylim = range(c(lo_W, hi_W, 0)),
  xlab = expression(alpha),
  ylab = expression(hat(theta)[W]),
  main = expression(paste("Wealth effect across ", alpha))
)

arrows(
  valid_41$alpha, lo_W,
  valid_41$alpha, hi_W,
  angle = 90, code = 3, length = 0.03
)

abline(h = 0, col = "grey60")
abline(v = best_alpha_41, col = "red", lty = 2, lwd = 2)



# ========================================
# Prior predictive edge-count simulation (ERGM)
# ========================================

set.seed(block_seed("prior_pred_ergm"))

S <- 6000  # use more for the final analysis

# Draw parameters from the priors used in log_prior above.
# NB: this previously drew theta_GWD from Uniform(-2, 2) while the model
# assumed Uniform(-1, 1). Both now read from PRIOR_GWD.
theta_prior <- cbind(
  theta_E   = rnorm(S, mean = PRIOR_E$mean, sd = PRIOR_E$sd),
  theta_W   = rgamma(S, shape = PRIOR_W$shape, rate = PRIOR_W$rate),
  theta_GWD = runif(S, min = PRIOR_GWD$min, max = PRIOR_GWD$max)
)

prior_pred_edges_ergm <- rep(NA_real_, S)

for (s in seq_len(S)) {
  
  simulated_stats <- tryCatch(
    simulate(
      form.41,
      coef = unname(theta_prior[s, ]),
      nsim = 1,
      output = "stats",
      control = control.simulate.formula(
        MCMC.burnin = 5000,
        MCMC.interval = 1
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(simulated_stats)) {
    prior_pred_edges_ergm[s] <-
      as.matrix(simulated_stats)[1, "edges"]
  }
  
  if (s %% 100 == 0) {
    message("Completed ", s, " of ", S)
  }
}

n_failed <- sum(is.na(prior_pred_edges_ergm))
n_failed

valid_prior_edges <- prior_pred_edges_ergm[!is.na(prior_pred_edges_ergm)]

par(
  mfrow = c(1, 1),
  mar = c(5.1, 5.1, 4.1, 2.1),
  mgp = c(3.5, 1, 0)
)

hist(
  valid_prior_edges,
  breaks = seq(-0.5, 120.5, by = 1),
  probability = TRUE,
  col = "grey80",
  border = "white",
  xlim = c(0, 60),
  main = "Prior predictive edge count",
  xlab = "Number of edges",
  ylab = "Proportion of simulated graphs"
)

abline(
  v = 20,
  col = "black",
  lwd = 2,
  lty = 4
)

legend(
  "topright",
  legend = "Observed = 20",
  col = "black",
  lty = 4,
  lwd = 2,
  bty = "n"
)

# Empty or complete graphs
mean(valid_prior_edges %in% c(0, 120))

# Nearly empty or nearly complete graphs
mean(valid_prior_edges <= 2 | valid_prior_edges >= 118)

summary(valid_prior_edges)

# ========================================
# Frequentist fit with the identical specification
# ========================================

set.seed(block_seed("ergm_mle_41"))

flomodel.41_mle <- ergm(form.41)

summary(flomodel.41_mle)

mle_values <- coef(flomodel.41_mle)

mle_E   <- mle_values[1]
mle_W   <- mle_values[2]
mle_GWD <- mle_values[3]

mle_values

# ========================================
# Posterior predictive edge-count simulation
# ========================================

set.seed(block_seed("post_pred_ergm"))

S_post <- min(6000, nrow(bergm_draws))

# Select posterior parameter draws
draw_index <- sample(
  seq_len(nrow(bergm_draws)),
  size = S_post,
  replace = FALSE
)

theta_post <- bergm_draws[draw_index, , drop = FALSE]

posterior_pred_edges <- rep(NA_real_, S_post)

for (s in seq_len(S_post)) {
  
  simulated_stats <- tryCatch(
    simulate(
      form.41,
      coef = unname(theta_post[s, ]),
      nsim = 1,
      output = "stats",
      control = control.simulate.formula(
        MCMC.burnin = 5000,
        MCMC.interval = 1
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(simulated_stats)) {
    posterior_pred_edges[s] <-
      as.matrix(simulated_stats)[1, "edges"]
  }
  
  if (s %% 100 == 0) {
    message("Completed ", s, " of ", S_post)
  }
}

sum(is.na(posterior_pred_edges))

valid_post_edges <- posterior_pred_edges[!is.na(posterior_pred_edges)]

par(
  mfrow = c(1, 1),
  mar = c(5.1, 5.1, 4.1, 2.1),
  mgp = c(3.5, 1, 0)
)

hist(
  valid_post_edges,
  breaks = seq(-0.5, 120.5, by = 1),
  probability = TRUE,
  col = "grey70",
  border = "white",
  xlim = c(0, 60),
  main = "Posterior predictive edge count",
  xlab = "Number of edges",
  ylab = "Proportion of simulated graphs"
)

abline(
  v = 20,
  col = "black",
  lwd = 2,
  lty = 4
)

legend(
  "topright",
  legend = "Observed = 20",
  col = "black",
  lty = 4,
  lwd = 2,
  bty = "n"
)

# ========================================
# Prior and posterior predictive overlay
# ========================================

# Identical breaks for both distributions
edge_breaks <- seq(-0.5, 120.5, by = 1)

prior_hist <- hist(
  valid_prior_edges,
  breaks = edge_breaks,
  probability = TRUE,
  plot = FALSE
)

post_hist <- hist(
  valid_post_edges,
  breaks = edge_breaks,
  probability = TRUE,
  plot = FALSE
)

ymax <- 1.1 * max(
  prior_hist$density,
  post_hist$density,
  na.rm = TRUE
)

col_prior <- adjustcolor("grey50", alpha.f = 0.45)
col_post  <- adjustcolor("steelblue", alpha.f = 0.40)

par(
  mfrow = c(1, 1),
  mar = c(5.1, 5.1, 4.1, 2.1),
  mgp = c(3.5, 1, 0)
)

# Prior predictive
plot(
  prior_hist,
  freq = FALSE,
  col = col_prior,
  border = "white",
  xlim = c(0, 60),
  ylim = c(0, ymax),
  main = "Prior and posterior predictive edge counts",
  xlab = "Number of edges",
  ylab = "Proportion of simulated graphs"
)

# Posterior predictive
plot(
  post_hist,
  freq = FALSE,
  col = col_post,
  border = "white",
  add = TRUE
)

# Observed edge count
abline(
  v = 20,
  col = "black",
  lwd = 2,
  lty = 4
)

legend(
  "topright",
  legend = c(
    "Prior predictive",
    "Posterior predictive",
    "Observed = 20"
  ),
  col = c("grey50", "steelblue", "black"),
  pch = c(15, 15, NA),
  pt.cex = 2,
  lty = c(NA, NA, 4),
  lwd = c(NA, NA, 2),
  bty = "n"
)

# =========================================================================
# 2. Table 4.7: posterior summaries with a mode column
# =========================================================================

posterior_summary <- function(draws, lower = NULL, upper = NULL) {
  c(
    mean   = mean(draws),
    median = median(draws),
    mode   = post_mode(draws, lower = lower, upper = upper),
    sd     = sd(draws),
    q2.5   = unname(quantile(draws, 0.025)),
    q97.5  = unname(quantile(draws, 0.975))
  )
}

table_4_7 <- rbind(
  theta_E   = posterior_summary(theta_E_post),
  theta_W   = posterior_summary(theta_W_post,   lower = 0),
  theta_GWD = posterior_summary(theta_GWD_post,
                                lower = PRIOR_GWD$min,
                                upper = PRIOR_GWD$max)
)

round(table_4_7, 3)

# Prior modes, for the comparison made in the surrounding text
prior_modes <- c(
  theta_E   = PRIOR_E$mean,                              # Normal: mode = mean
  theta_W   = (PRIOR_W$shape - 1) / PRIOR_W$rate,        # Gamma(2,4) -> 0.25
  theta_GWD = NA                                          # Uniform: no unique mode
)
prior_modes


# =========================================================================
# 3. How stable is the mode?
# =========================================================================
#
# The mean and median are order statistics of the draws; the mode is not.
# Two cheap checks, both worth a sentence in the text. Expect theta_GWD to
# be the unstable one: its posterior is nearly flat near the peak, so small
# changes in bandwidth or chain move the argmax a long way.

# --- 3a. Bandwidth sensitivity ------------------------------------------

mode_sensitivity <- function(draws, lower = NULL, upper = NULL) {
  sapply(c("nrd0", "nrd", "SJ"), function(b) {
    tryCatch(
      suppressWarnings(
        post_mode(draws, lower = lower, upper = upper, bw = b)
      ),
      error = function(e) NA_real_
    )
  })
}

rbind(
  theta_E   = mode_sensitivity(theta_E_post),
  theta_W   = mode_sensitivity(theta_W_post,   lower = 0),
  theta_GWD = mode_sensitivity(theta_GWD_post,
                               lower = PRIOR_GWD$min,
                               upper = PRIOR_GWD$max)
)

# --- 3b. Chain-to-chain variation ---------------------------------------
#
# Bergm stores the draws chain-blocked (iteration varies fastest), which is
# the same layout plot_bergm_chains() already relies on.

mode_by_chain <- function(draws, main.iters = MI,
                          lower = NULL, upper = NULL) {
  
  nchains <- length(draws) / main.iters
  
  if (nchains != round(nchains)) {
    stop("main.iters does not match the stored draws.")
  }
  
  chain_id <- rep(seq_len(nchains), each = main.iters)
  
  tapply(draws, chain_id, function(z) {
    post_mode(z, lower = lower, upper = upper)
  })
}

rbind(
  theta_E   = mode_by_chain(theta_E_post),
  theta_W   = mode_by_chain(theta_W_post,   lower = 0),
  theta_GWD = mode_by_chain(theta_GWD_post,
                            lower = PRIOR_GWD$min,
                            upper = PRIOR_GWD$max)
)


# =========================================================================
# 4. Figure 4.11: revised plotting function with a mode line
# =========================================================================
#
# Changes from the version in the main script:
#   - the density is built by post_density() on [lower, upper], so the mode
#     drawn here is identical to the mode in table_4_7;
#   - xlim now controls the axis only, not the density grid;
#   - the mode is drawn as a dotted line in the parameter colour and the
#     value is returned invisibly.

plot_prior_posterior <- function(
    draws,
    prior_function,
    mle_value,
    xlim,
    ylim,
    colour,
    x_label,
    plot_title,
    show_ylabel = TRUE,
    lower = NULL,
    upper = NULL,
    show_mode = TRUE
) {
  
  # Ensure the MLE is inside the plotting range
  xlim <- range(c(xlim, mle_value))
  
  x <- seq(xlim[1], xlim[2], length.out = 1000)
  
  posterior_density <- post_density(draws, lower = lower, upper = upper)
  
  mode_value <- posterior_density$x[which.max(posterior_density$y)]
  
  prior_density <- prior_function(x)
  
  plot(
    posterior_density$x,
    posterior_density$y,
    type = "l",
    col  = colour,
    lwd  = 2,
    lty  = 1,
    xlim = xlim,
    ylim = ylim,
    xlab = x_label,
    ylab = if (show_ylabel) "Density" else "",
    main = plot_title,
    font.main = 2
  )
  
  # Prior
  lines(
    x,
    prior_density,
    col = colour,
    lwd = 2,
    lty = 2
  )
  
  # Posterior mode
  if (show_mode) {
    segments(
      x0 = mode_value, y0 = 0,
      x1 = mode_value, y1 = max(posterior_density$y),
      col = colour,
      lwd = 2,
      lty = 3
    )
  }
  
  # Posterior density at the MLE
  mle_height <- approx(
    posterior_density$x,
    posterior_density$y,
    xout = mle_value,
    rule = 3
  )$y
  
  # MLE reference line
  segments(
    x0 = mle_value, y0 = 0,
    x1 = mle_value, y1 = mle_height,
    col = "black",
    lwd = 2,
    lty = 4
  )
  
  legend(
    "topright",
    legend = c(
      paste0("MLE \u2248 ",  round(mle_value,  3)),
      if (show_mode) paste0("MAP \u2248 ", round(mode_value, 3)),
      "Posterior",
      "Prior"
    ),
    col = c("black", if (show_mode) colour, colour, colour),
    lty = c(4, if (show_mode) 3, 1, 2),
    lwd = 2,
    bty = "n"
  )
  
  invisible(mode_value)
}


# --- The figure itself, capturing the modes as it draws them ------------

par(
  mfrow = c(1, 3),
  oma = c(0, 0, 3, 0),
  mar = c(4, 4, 1, 2) + 0.1,
  mgp = c(3, 1, 0)
)

fig_modes <- numeric(3)
names(fig_modes) <- c("theta_E", "theta_W", "theta_GWD")

fig_modes["theta_E"] <- plot_prior_posterior(
  draws = theta_E_post,
  prior_function = function(x) dnorm(x, mean = PRIOR_E$mean, sd = PRIOR_E$sd),
  mle_value = mle_E,
  xlim = c(-3, 0),
  ylim = c(0, 2.1),
  colour = "red",
  x_label = expression(theta[E]),
  plot_title = "",
  show_ylabel = TRUE
)

fig_modes["theta_W"] <- plot_prior_posterior(
  draws = theta_W_post,
  prior_function = function(x) dgamma(x, shape = PRIOR_W$shape, rate = PRIOR_W$rate),
  mle_value = mle_W,
  xlim = c(0, 2.5),
  ylim = c(0, 3),
  colour = "blue",
  x_label = expression(theta[W]),
  plot_title = "",
  show_ylabel = FALSE,
  lower = 0
)

fig_modes["theta_GWD"] <- plot_prior_posterior(
  draws = theta_GWD_post,
  prior_function = function(x) dunif(x, min = PRIOR_GWD$min, max = PRIOR_GWD$max),
  mle_value = mle_GWD,
  xlim = c(-2, 2),
  ylim = c(0, 1),
  colour = "darkgreen",
  x_label = expression(theta[GWD]),
  plot_title = "",
  show_ylabel = FALSE,
  lower = PRIOR_GWD$min,
  upper = PRIOR_GWD$max
)

mtext(
  "Edge, wealth and GWD model priors and posteriors",
  side = 3,
  outer = TRUE,
  line = 1,
  cex = 1.5,
  font = 2
)

# Consistency check: these must match the mode column of table_4_7 exactly.
cbind(figure = fig_modes, table = table_4_7[, "mode"])





# =============================================================================
# PART II: FLORENTINE GOODNESS OF FIT (CHAPTER 5)
# =============================================================================

# =========================================================================
# Frequentist vs Bayesian goodness of fit for flomodel.41
#   flomarriage ~ edges + nodecov("wealth_z") + gwdegree(0.5, fixed = TRUE)
#
# Four panels: degree, edgewise shared partners, dyadwise shared partners,
# geodesic distance. Each bin gets two dodged boxes, MLE and Bayesian,
# with the median as the box line and the mean as an open diamond.
#
# Uses these objects created in Part I:
#   flomarriage, form.41, flomodel.41, flomodel.41_mle, bergm_draws,
#   n_nodes, mle_values, PRIOR_E, PRIOR_W, PRIOR_GWD, block_seed()
# =========================================================================

# =========================================================================
# 0. Design of the comparison
# =========================================================================
#
# Both procedures simulate networks and compare the observed statistics to
# the simulated distribution. The only thing that differs is where theta
# comes from:
#
#   Frequentist  theta = theta_hat, held fixed and reused S times
#                -> spread reflects graph variability at one theta
#
#   Bayesian     theta ~ p(theta | y), a fresh draw for every network
#                -> spread reflects graph variability AND uncertainty about
#                   theta, filtered through the prior
#
# Everything downstream is identical: same simulation engine, same
# burn-in, same number of draws, same statistics. Do not substitute gof()
# for the frequentist arm, because its default simulation controls differ
# from those used in the posterior predictive block, and then part of any
# difference in spread is an artefact of the settings rather than of the
# paradigm.
#
# On the four panels:
#   degree    partly in-model via GWD, so read it as a calibration check
#   esp       closure among families that actually married
#   dsp       closure opportunity across all pairs, tied or not
#   distance  global connectivity
#
# esp and dsp together are the informative pair. The model has a degree
# term but no gwesp, so neither closure statistic was fitted. If esp is
# understated while dsp is about right, the model generates the
# opportunity for triads but not the tendency to close them, which is the
# standard argument for adding gwesp.

# =========================================================================
# 1. Are the prior and the likelihood describing the same model?
# =========================================================================
#
#   theta_W   ~ Gamma(2, 4)          support (0, Inf)     HARD
#   theta_GWD ~ Uniform(-1.9, 1.9)   support [-1.9, 1.9]  HARD
#
# If theta_hat falls outside either, the posterior cannot reach the MLE and
# the two arms are fitting different models.

prior_support <- data.frame(
  parameter = names(mle_values),
  mle       = as.numeric(mle_values),
  lower     = c(-Inf, 0, PRIOR_GWD$min),
  upper     = c(Inf, Inf, PRIOR_GWD$max),
  row.names = NULL
)

prior_support$inside_support <- with(prior_support, mle > lower & mle < upper)

prior_support$prior_quantile <- c(
  pnorm(mle_values[1], PRIOR_E$mean, PRIOR_E$sd),
  pgamma(mle_values[2], PRIOR_W$shape, PRIOR_W$rate),
  punif(mle_values[3], PRIOR_GWD$min, PRIOR_GWD$max)
)

prior_support

# A prior quantile near 0 or 1 means the prior is pulling against the
# likelihood, so the Bayesian boxes will be shifted as well as rescaled.

# =========================================================================
# 2. Observed statistics
# =========================================================================

SP_RANGE <- 0:5    # shared-partner bins used for both esp and dsp

obs_esp <- summary(flomarriage ~ esp(SP_RANGE))
obs_dsp <- summary(flomarriage ~ dsp(SP_RANGE))
obs_nsp <- summary(flomarriage ~ nsp(SP_RANGE))

obs_scalar <- summary(
  flomarriage ~ edges + triangle + absdiff("wealth_z") + gwesp(0.5, fixed = TRUE)
)

obs_degree <- tabulate(
  sna::degree(flomarriage, gmode = "graph") + 1,
  nbins = n_nodes
)

obs_geodist <- {
  g <- sna::geodist(flomarriage, inf.replace = Inf)$gdist
  g <- g[upper.tri(g)]
  c(tabulate(g[is.finite(g)], nbins = n_nodes - 1), sum(!is.finite(g)))
}

# =========================================================================
# 3. Shared simulation engine
# =========================================================================

S_arm      <- 500      # networks per arm
SIM_BURNIN <- 20000    # identical for both arms

simulate_gof_arm <- function(theta_mat, label, burnin = SIM_BURNIN) {
  
  S <- nrow(theta_mat)
  
  out <- list(
    theta   = theta_mat,
    degree  = matrix(NA_real_, S, n_nodes),
    esp     = matrix(NA_real_, S, length(obs_esp)),
    dsp     = matrix(NA_real_, S, length(obs_dsp)),
    nsp     = matrix(NA_real_, S, length(obs_nsp)),
    geodist = matrix(NA_real_, S, n_nodes),
    scalar  = matrix(NA_real_, S, length(obs_scalar),
                     dimnames = list(NULL, names(obs_scalar)))
  )
  
  for (s in seq_len(S)) {
    
    net <- tryCatch(
      simulate(
        form.41,
        coef = unname(theta_mat[s, ]),
        nsim = 1,
        output = "network",
        control = control.simulate.formula(
          MCMC.burnin = burnin,
          MCMC.interval = 1
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(net)) next
    
    d <- sna::degree(net, gmode = "graph")
    
    out$degree[s, ] <- tabulate(d + 1, nbins = n_nodes)
    out$esp[s, ]    <- summary(net ~ esp(SP_RANGE))
    out$dsp[s, ]    <- summary(net ~ dsp(SP_RANGE))
    out$nsp[s, ]    <- summary(net ~ nsp(SP_RANGE))
    
    out$scalar[s, ] <- summary(
      net ~ edges + triangle + absdiff("wealth_z") + gwesp(0.5, fixed = TRUE)
    )
    
    g <- sna::geodist(net, inf.replace = Inf)$gdist
    g <- g[upper.tri(g)]
    
    out$geodist[s, ] <- c(
      tabulate(g[is.finite(g)], nbins = n_nodes - 1),
      sum(!is.finite(g))
    )
    
    if (s %% 100 == 0) message(label, ": ", s, " of ", S)
  }
  
  out
}

# =========================================================================
# 4. The two arms
# =========================================================================

set.seed(block_seed("gof_compare_freq"))

theta_freq <- matrix(
  rep(as.numeric(mle_values), each = S_arm),
  nrow = S_arm
)

arm_freq <- simulate_gof_arm(theta_freq, "frequentist (plug-in MLE)")

set.seed(block_seed("gof_compare_bayes"))

bayes_index <- sample(
  seq_len(nrow(bergm_draws)),
  size = min(S_arm, nrow(bergm_draws)),
  replace = FALSE
)

arm_bayes <- simulate_gof_arm(
  bergm_draws[bayes_index, , drop = FALSE],
  "Bayesian posterior"
)

arms <- list(
  "Frequentist (MLE)" = arm_freq,
  "Bayesian (BERGM)"  = arm_bayes
)

arm_cols <- c("darkred", "steelblue")

# =========================================================================
# 5. Four-panel box-and-whisker figure
# =========================================================================
#
# Boxes are drawn by hand rather than via bxp(), so the whiskers can be
# set to the 2.5% and 97.5% quantiles instead of Tukey's 1.5 x IQR rule.
# Two reasons:
#
#   - it keeps continuity with the 95% intervals used everywhere else in
#     the chapter, so the tables and the figure describe the same thing
#   - with 500 draws and occasional near-complete graphs, Tukey whiskers
#     plus outlier points would fill the upper bins with dots that say
#     nothing beyond what the degeneracy table already reports
#
# Box  = interquartile range
# Line = median
# Diamond = mean, plotted so the skew is visible rather than argued about

box_stats <- function(v, whisker_probs = c(0.025, 0.975)) {
  v <- v[!is.na(v)]
  c(
    quantile(v, whisker_probs[1], names = FALSE),
    quantile(v, 0.25, names = FALSE),
    median(v),
    quantile(v, 0.75, names = FALSE),
    quantile(v, whisker_probs[2], names = FALSE)
  )
}

draw_box <- function(x, st, width, fill, border) {
  # st = c(lower whisker, Q1, median, Q3, upper whisker)
  segments(x, st[1], x, st[2], col = border, lwd = 1.4)
  segments(x, st[4], x, st[5], col = border, lwd = 1.4)
  segments(x - width / 4, st[1], x + width / 4, st[1], col = border, lwd = 1.4)
  segments(x - width / 4, st[5], x + width / 4, st[5], col = border, lwd = 1.4)
  rect(x - width / 2, st[2], x + width / 2, st[4],
       col = fill, border = border, lwd = 1.4)
  segments(x - width / 2, st[3], x + width / 2, st[3], col = border, lwd = 2.6)
}

plot_gof_boxes <- function(
    component,
    obs_vec,
    keep,
    x_labels,
    main,
    xlab,
    ylab,
    box_width = 0.30,
    offsets = c(-0.18, 0.18),
    connect = seq_along(keep),
    show_legend = FALSE,
    legend_pos = "topleft"
) {
  
  k   <- length(keep)
  pos <- seq_len(k)
  
  stats_list <- lapply(arms, function(a) {
    apply(a[[component]][, keep, drop = FALSE], 2, box_stats)
  })
  
  means_list <- lapply(arms, function(a) {
    colMeans(a[[component]][, keep, drop = FALSE], na.rm = TRUE)
  })
  
  ymax <- 1.05 * max(
    unlist(stats_list), unlist(means_list), obs_vec, na.rm = TRUE
  )
  
  plot(
    NA,
    xlim = c(0.5, k + 0.5),
    ylim = c(0, ymax),
    xaxt = "n",
    xlab = xlab, ylab = ylab, main = main, font.main = 2
  )
  
  axis(side = 1, at = pos, labels = x_labels)
  
  for (a in seq_along(arms)) {
    
    st <- stats_list[[a]]
    
    for (j in seq_len(k)) {
      draw_box(
        x = pos[j] + offsets[a],
        st = st[, j],
        width = box_width,
        fill = adjustcolor(arm_cols[a], alpha.f = 0.30),
        border = arm_cols[a]
      )
    }
    
    points(
      pos + offsets[a], means_list[[a]],
      pch = 23, bg = "white", col = arm_cols[a], cex = 0.85, lwd = 1.6
    )
  }
  
  # Observed
  # Observed
  if (length(connect) > 1) {
    lines(pos[connect], obs_vec[connect], col = "black", lwd = 1.6, lty = 2)
  }
  points(pos, obs_vec, pch = 19, cex = 1.0)
  span <- max(abs(offsets)) + box_width / 2
  segments(pos - span, obs_vec, pos + span, obs_vec, col = "black", lwd = 2.2, lty = 3)

  
  if (show_legend) {
    legend(
      legend_pos,
      legend = c(
        "Observed",
        names(arms),
        "Median (box line)",
        "Mean"
      ),
      col = c("black", arm_cols, "grey20", "grey20"),
      pt.bg = c("black", adjustcolor(arm_cols, alpha.f = 0.30), NA, "white"),
      pch = c(19, 22, 22, NA, 23),
      lty = c(2, NA, NA, 1, NA),
      lwd = c(1.6, NA, NA, 2.6, 1.6),
      pt.cex = c(1.0, 1.6, 1.6, NA, 0.85),
      bty = "n", cex = 0.7
    )
  }
}

par(
  mfrow = c(1, 1),
  oma = c(0, 0, 0, 0),
  mar = c(4.2, 3.8, 2.2, 1.2),
  mgp = c(2.6, 0.8, 0)
)

#par(mfrow = c(2, 2), oma = c(0, 0, 3.2, 0), mar = c(4.2, 4.4, 2.2, 1.2), mgp = c(2.6, 0.8, 0))

plot_gof_boxes(
  "degree", obs_degree[1:9], keep = 1:9, x_labels = 0:8,
  main = "Degree Distribution GOF",
  xlab = "Degree", ylab = "Number of families",
  show_legend = TRUE
)

plot_gof_boxes(
  "esp", obs_esp, keep = seq_along(obs_esp), x_labels = SP_RANGE,
  main = "Edgewise Shared Partners Distribution GOF",
  xlab = "Shared partners", ylab = "Number of edges",
  show_legend = TRUE
)

plot_gof_boxes(
  "dsp", obs_dsp, keep = seq_along(obs_dsp), x_labels = SP_RANGE,
  main = "Dyadwise Shared Partners Distribution GOF",
  xlab = "Shared partners", ylab = "Number of dyads",
  show_legend = TRUE
)

plot_gof_boxes(
  "geodist", obs_geodist[c(1:6, n_nodes)], keep = c(1:6, n_nodes),
  x_labels = c(1:6, "NR"),
  connect = 1:6,
  main = "Geodesic Distance Distribution GOF",
  xlab = "Distance", ylab = "Number of dyads",
  show_legend = TRUE
)

plot_gof_boxes(
  "nsp", obs_nsp, keep = seq_along(obs_nsp), x_labels = SP_RANGE,
  main = "Non-edgewise shared partners GOF",
  xlab = "Shared partners", ylab = "Number of non-tied dyads",
  show_legend = TRUE
)

mtext(
  "Degree Distribution GOF",
  side = 3, outer = TRUE, line = 1, cex = 1.3, font = 2
)

par(mfrow = c(1, 1))

# The dsp panel is dominated by the zero bin, since most pairs in a
# 16-node graph with 20 edges share no partners. If that squashes the
# rest, redraw with keep = 2:6 and x_labels = 1:5 as a companion panel and
# say in the caption that the zero bin was dropped for legibility.

# =========================================================================
# 6. Where do the mean and the median disagree?
# =========================================================================
#
# A large gap between the diamond and the box line means that bin's
# predictive distribution is strongly right-skewed. Worth reporting: it is
# the reason the median is the honest central summary here, and it also
# flags which bins the degeneracy table in section 10 is speaking about.

skew_check <- function(component, keep, bin_labels) {
  do.call(rbind, lapply(names(arms), function(nm) {
    m <- arms[[nm]][[component]][, keep, drop = FALSE]
    data.frame(
      arm    = nm,
      bin    = bin_labels,
      mean   = round(colMeans(m, na.rm = TRUE), 2),
      median = round(apply(m, 2, median, na.rm = TRUE), 2),
      gap    = round(colMeans(m, na.rm = TRUE) -
                       apply(m, 2, median, na.rm = TRUE), 2),
      row.names = NULL
    )
  }))
}

skew_check("degree",  1:9, paste("degree", 0:8))
skew_check("esp",     seq_along(obs_esp), paste("esp", SP_RANGE))
skew_check("dsp",     seq_along(obs_dsp), paste("dsp", SP_RANGE))
skew_check("geodist", c(1:6, n_nodes), c(paste("dist", 1:6), "unreachable"))

# =========================================================================
# 7. Spread of the boxes
# =========================================================================

interval_widths <- function(mat, keep, probs = c(0.025, 0.975)) {
  apply(mat[, keep, drop = FALSE], 2, function(v) {
    diff(quantile(v, probs = probs, na.rm = TRUE))
  })
}

width_table <- function(component, keep, bin_labels, probs = c(0.025, 0.975)) {
  
  w <- sapply(arms, function(a) interval_widths(a[[component]], keep, probs))
  
  out <- data.frame(bin = bin_labels, round(w, 2), row.names = NULL)
  names(out) <- c("bin", names(arms))
  
  out$ratio <- round(out[["Bayesian (BERGM)"]] / out[["Frequentist (MLE)"]], 2)
  out
}

# Whisker-to-whisker spread
width_table("degree",  1:9, paste("degree", 0:8))
width_table("esp",     seq_along(obs_esp), paste("esp", SP_RANGE))
width_table("dsp",     seq_along(obs_dsp), paste("dsp", SP_RANGE))
width_table("geodist", c(1:6, n_nodes), c(paste("dist", 1:6), "unreachable"))

# Box height (IQR) only, which is what the eye actually compares
width_table("esp", seq_along(obs_esp), paste("esp", SP_RANGE),
            probs = c(0.25, 0.75))

# ratio > 1  the Bayesian box is taller: the plug-in understates
#            predictive uncertainty by ignoring uncertainty in theta
# ratio ~ 1  parameter uncertainty is small relative to graph variability
#            at this network size
# ratio < 1  the prior is tightening the box. With Gamma(2,4) on theta_W
#            and Uniform(-1.9, 1.9) on theta_GWD this is a real
#            possibility, and it needs defending rather than reporting.

# =========================================================================
# 8. p-values, and why neither is calibrated
# =========================================================================

pval_table <- do.call(rbind, lapply(names(arms), function(nm) {
  
  sim <- arms[[nm]]$scalar
  
  data.frame(
    arm        = nm,
    statistic  = colnames(sim),
    observed   = as.numeric(obs_scalar),
    sim_mean   = colMeans(sim, na.rm = TRUE),
    sim_median = apply(sim, 2, median, na.rm = TRUE),
    sim_lower  = apply(sim, 2, quantile, probs = 0.025, na.rm = TRUE),
    sim_upper  = apply(sim, 2, quantile, probs = 0.975, na.rm = TRUE),
    p_upper    = sapply(seq_len(ncol(sim)), function(j) {
      mean(sim[, j] >= obs_scalar[j], na.rm = TRUE)
    }),
    row.names = NULL
  )
}))

pval_table$p_two_sided <- 2 * pmin(pval_table$p_upper, 1 - pval_table$p_upper)
pval_table$inside_95   <- with(pval_table, observed >= sim_lower & observed <= sim_upper)

pval_table[order(pval_table$statistic, pval_table$arm), ]

# Miscalibrated in opposite directions:
#
#   Frequentist bootstrap p: treats theta_hat as known, so the reference
#     distribution is too narrow and p is pushed towards 0 or 1. Declares
#     misfit too readily.
#
#   Posterior predictive p: the data are used to fit and then to check, so
#     it is conservative and pulled towards 0.5. Declares misfit too
#     reluctantly.

# =========================================================================
# 9. Approximate variance decomposition
# =========================================================================
#
#   Var(T(y)) = E_theta[ Var(T | theta) ] + Var_theta( E[T | theta] )
#                \___ within ___/           \___ between ___/
#
# The frequentist arm estimates Var(T | theta_hat), standing in for the
# within term; the Bayesian arm estimates the total. The gap approximates
# the between term.
#
# Approximate only: Var(T | theta) is not constant in theta, so the
# frequentist arm evaluates the within term at a point rather than
# averaging over the posterior. between_approx can come out negative,
# which is the decomposition breaking down rather than a coding error.

var_decomposition <- data.frame(
  statistic     = colnames(arm_freq$scalar),
  within_at_mle = apply(arm_freq$scalar,  2, var, na.rm = TRUE),
  total_bayes   = apply(arm_bayes$scalar, 2, var, na.rm = TRUE),
  row.names = NULL
)

var_decomposition$between_approx <-
  var_decomposition$total_bayes - var_decomposition$within_at_mle

var_decomposition$prop_from_theta <-
  var_decomposition$between_approx / var_decomposition$total_bayes

# round() on the whole frame fails because of the character column
var_decomposition[, -1] <- round(var_decomposition[, -1], 3)
var_decomposition

# =========================================================================
# 10. Degeneracy rates
# =========================================================================
#
# The whiskers stop at the 2.5% and 97.5% quantiles, so the extremes are
# deliberately not drawn. This table is where they are accounted for.

degeneracy_table <- do.call(rbind, lapply(names(arms), function(nm) {
  
  e <- arms[[nm]]$scalar[, "edges"]
  
  data.frame(
    arm             = nm,
    failed_sims     = mean(is.na(e)),
    exactly_extreme = mean(e %in% c(0, 120), na.rm = TRUE),
    near_extreme    = mean(e <= 2 | e >= 118, na.rm = TRUE),
    median_edges    = median(e, na.rm = TRUE),
    mean_edges      = mean(e, na.rm = TRUE),
    row.names = NULL
  )
}))

degeneracy_table







# =============================================================================
# PART III: IS-E NETWORK (CHAPTER 6)
# =============================================================================

## ===========================================================================
## IS-E network: import, attributes, figures
## Data: Gutfraind & Genkin (2017), https://github.com/sashagutfraind/Nov13
##
## Run top to bottom. Every dependency now comes before its use.
## igraph calls are namespace-qualified because sna is ahead of igraph on
## the search path once statnet is loaded.
## ===========================================================================

## ===========================================================================
## 1. Load graphs and attributes
## ===========================================================================

g   <- read_graph(ise_data_path("ise_terrorNetwork.gml"),         format = "gml")
g_e <- read_graph(ise_data_path("ise_terrorNetworkExtended.gml"), format = "gml")
g_l <- read_graph(ise_data_path("ise_terrorNetworkLimited.gml"),  format = "gml")

## latin1: the source CSVs contain single high bytes for the accented names
## (Mostefai, Belkaid). Reading them correctly removes the need to patch
## individual rows by index.
people <- read.csv(ise_data_path("ise_personDetails.csv"), stringsAsFactors = FALSE,
                   fileEncoding = "latin1")

idx <- match(igraph::V(g)$name, people$p.name)
if (any(is.na(idx))) {
  stop("unmatched vertex names: ",
       paste(igraph::V(g)$name[is.na(idx)], collapse = "; "))
}

igraph::V(g)$age         <- people$p.age[idx]
igraph::V(g)$citizenship <- people$p.citizenship[idx]
igraph::V(g)$status      <- people$p.status[idx]
igraph::V(g)$role        <- people$p.role[idx]

## Tidy up
igraph::V(g)$status[igraph::V(g)$status %in% "arreted"]     <- "arrested"  # source typo
igraph::V(g)$status[is.na(igraph::V(g)$status)]           <- "unknown"
igraph::V(g)$role[is.na(igraph::V(g)$role)]               <- "none"
igraph::V(g)$citizenship[is.na(igraph::V(g)$citizenship)] <- "unknown"

table(igraph::V(g)$status,      useNA = "ifany")
table(igraph::V(g)$citizenship, useNA = "ifany")
table(igraph::V(g)$role,        useNA = "ifany")

## In the database but not the graph (status == "free", no tie generated)
setdiff(people$p.name, igraph::V(g)$name)


## ===========================================================================
## 2. Attack involvement
## ===========================================================================

## The one-mode projection discarded the person-event incidences, so go back
## to the raw relationship dump to recover them.
rels <- read.csv(ise_data_path("ise_relationships.csv"), stringsAsFactors = FALSE,
                 fileEncoding = "latin1")
ties <- read.csv(ise_data_path("ise_networkTieDetails.csv"), stringsAsFactors = FALSE,
                 fileEncoding = "latin1")

paris_sites <- c("Bataclan", "Stade de France", "La Belle Equipe",
                 "La Bonne Biere", "La Casa Nostra",
                 "Le Carillon Bar and Le Petit Cambodge", "Comptoir Voltaire")
brussels_sites <- c("Brussels Airport Zaventem", "Maalbeck Metro")

miss <- setdiff(c(paris_sites, brussels_sites), c(rels$n1.name, rels$n2.name))
if (length(miss)) {
  warning("attack sites not matched: ", paste(miss, collapse = "; "))
}

## INVOLVED_IN only: PRESENT_IN would sweep in anyone who visited a site
attackers_at <- function(sites) {
  r <- rels[rels$type.r. == "INVOLVED_IN", ]
  who <- unique(c(r$n1.name[r$n2.name %in% sites],
                  r$n2.name[r$n1.name %in% sites]))
  intersect(who, igraph::V(g)$name)
}

paris    <- attackers_at(paris_sites)
brussels <- attackers_at(brussels_sites)

nmv <- igraph::V(g)$name
igraph::V(g)$attack <- factor(
  ifelse(nmv %in% paris & nmv %in% brussels, "Paris & Brussels",
         ifelse(nmv %in% paris,                     "Paris",
                ifelse(nmv %in% brussels,                  "Brussels", "Neither"))),
  levels = c("Paris", "Brussels", "Paris & Brussels", "Neither")
)

print(table(igraph::V(g)$attack))

## Everything downstream needs these
stopifnot(!is.null(igraph::V(g)$attack), !is.null(igraph::V(g)$status))


## ===========================================================================
## 3. Shared helpers
## ===========================================================================

## Placeholders that no wrapping rule can tidy
short_names <- c(
  "Unknown accomplice of Salah Abdeslam" = "S. Abdeslam accomplice",
  "Salzburg Fourth Suspect"              = "Salzburg suspect 4",
  "Dusseldorf Suspect"                   = "Dusseldorf suspect",
  "Giessen Suspect"                      = "Giessen suspect",
  "Ukraine Suspect"                      = "Ukraine suspect",
  "SalimBenghalem"                       = "Salim Benghalem",
  "August Recruit"                       = "August recruit"
)

## Split a name across exactly two balanced lines
two_lines <- function(x) {
  vapply(x, function(s) {
    if (is.na(s) || !grepl(" ", s)) return(s)
    w  <- strsplit(s, " ")[[1]]
    cl <- cumsum(nchar(w) + 1) - 1
    k  <- max(1, min(which.min(abs(cl - nchar(s) / 2)), length(w) - 1))
    paste(paste(w[seq_len(k)],  collapse = " "),
          paste(w[-seq_len(k)], collapse = " "), sep = "\n")
  }, character(1), USE.NAMES = FALSE)
}

## Rescale to [0, 1]
rs <- function(v) if (diff(range(v)) == 0) rep(0.5, length(v)) else
  (v - min(v)) / diff(range(v))

shorten <- function(x) {
  x[x %in% names(short_names)] <- short_names[x[x %in% names(short_names)]]
  x
}

iso <- igraph::degree(g) == 0


## ===========================================================================
## 4. Figure 1 — community structure
## ===========================================================================

gcon <- igraph::induced_subgraph(g, !iso)

set.seed(1)
cl <- igraph::cluster_louvain(gcon)

cat("modularity:", round(igraph::modularity(cl), 3), "\n")
print(igraph::sizes(cl))

## Names are keyed off a member rather than a community number, so they
## survive a change of seed. If a cluster splits or merges you get
## "Group N" instead of a silently mislabelled hull.
anchors <- c(
  "Abdelhamid Abaaoud" = "Abaaoud's contacts",
  "Salah Abdeslam"     = "Logistics & exfiltration",
  "Adel Haddadi"       = "Balkan route arrivals",
  "Jawad Benhattal"    = "Euro 2016 plot",
  "Reda Kriket"        = "Kriket plot (Argenteuil)",
  "Ahmet Dahmani"      = "Syria exfiltration",
  "Oussama Atar"       = "ISIS command"
)

lab_for <- function(k) {
  who <- igraph::V(gcon)$name[igraph::membership(cl) == k]
  hit <- anchors[names(anchors) %in% who]
  if (length(hit)) unname(hit[1]) else paste("Group", k)
}

min_size <- 3
memb <- igraph::membership(cl)
keep <- names(igraph::sizes(cl))[igraph::sizes(cl) >= min_size]

igraph::V(g)$comm <- NA_character_
igraph::V(g)$comm[!iso] <- ifelse(as.character(memb) %in% keep,
                                  vapply(memb, lab_for, character(1)),
                                  NA_character_)

## --- layout: connected part, plus a row of isolates beneath -----------------

set.seed(29)
L <- matrix(NA_real_, igraph::vcount(g), 2)
Lm <- igraph::layout_with_fr(gcon)
L[!iso, 1] <- rs(Lm[, 1])
L[!iso, 2] <- rs(Lm[, 2])
if (any(iso)) {
  L[iso, 1] <- seq(0.02, 0.98, length.out = sum(iso))
  L[iso, 2] <- -0.16
}

lay <- ggraph::create_layout(g, layout = "manual", x = L[, 1], y = L[, 2])
lay$deg <- igraph::degree(g)
lay$iso <- iso

nm <- shorten(igraph::V(g)$name)

lay$attacker <- as.character(igraph::V(g)$attack) != "Neither"
lay$died     <- lay$attacker & igraph::V(g)$status == "dead"
lay$attack_shape <- factor(
  ifelse(lay$attacker, as.character(igraph::V(g)$attack), "Not an attacker"),
  levels = c("Paris", "Brussels", "Paris & Brussels", "Not an attacker")
)

## Individuals to label. Must match V(g)$name exactly.
label_these <- c(
  "Abdelhamid Abaaoud",
  "Oussama Atar",
  "Salah Abdeslam",
  "Ibrahim Abdeslam",
  "Reda Kriket",
  "Abu Bakr al-Baghdadi",
  "Abu Muhammad al-Adnani",
  "Moustapha Benhattal",
  "Ahmet Dahmani",
  "Ahmed Almuhamed",
  "Bilal Hadfi",
  "Samy Amimour"
)

unmatched <- setdiff(label_these, igraph::V(g)$name)
if (length(unmatched)) {
  warning("label_these not found in the graph: ",
          paste(unmatched, collapse = "; "))
}

lay$lab <- ifelse(igraph::V(g)$name %in% label_these & !lay$iso,
                  two_lines(nm), NA_character_)

xr      <- range(lay$x)
pad     <- diff(xr)
left_x  <- xr[1] - pad * 0.05
right_x <- xr[2] + pad * 0.05
mid     <- stats::median(lay$x[!is.na(lay$lab)])

lay$lab_l   <- ifelse(!is.na(lay$lab) & lay$x <  mid, lay$lab, NA_character_)
lay$lab_r   <- ifelse(!is.na(lay$lab) & lay$x >= mid, lay$lab, NA_character_)
lay$lab_iso <- ifelse(lay$iso, two_lines(nm), NA_character_)

hull_data <- as.data.frame(lay)[!is.na(lay$comm), ]

p_clusters <- ggraph(lay) +
  geom_mark_hull(
    data = hull_data,
    aes(x = x, y = y, group = comm, label = comm, fill = comm),
    concavity = 4, expand = unit(3.5, "mm"), radius = unit(3.5, "mm"),
    alpha = 0.13, colour = NA,
    label.fontsize = 9, label.fill = alpha("white", 0.75),
    label.colour = "grey20", label.buffer = unit(4, "mm"),
    con.type = "elbow", con.colour = "grey55", con.size = 0.3
  ) +
  geom_edge_link(colour = "grey75", width = 0.4, alpha = 0.8) +
  geom_node_point(aes(size = deg, fill = comm, shape = attack_shape),
                  colour = "grey15", stroke = 0.7) +
  geom_node_point(aes(size = deg, filter = died),
                  shape = 4, colour = "grey10", stroke = 1) +
  geom_node_text(aes(label = lab_l),
                 repel = TRUE, xlim = c(NA, left_x), hjust = 1,
                 direction = "y", size = 3.2, lineheight = 0.85,
                 max.overlaps = Inf, force = 12, force_pull = 0.05,
                 max.iter = 30000, box.padding = 0.5, point.padding = 0.3,
                 min.segment.length = 0, segment.colour = "grey35",
                 segment.size = 0.4, colour = "grey15",
                 family = "sans", na.rm = TRUE) +
  geom_node_text(aes(label = lab_r),
                 repel = TRUE, xlim = c(right_x, NA), hjust = 0,
                 direction = "y", size = 3.2, lineheight = 0.85,
                 max.overlaps = Inf, force = 12, force_pull = 0.05,
                 max.iter = 30000, box.padding = 0.5, point.padding = 0.3,
                 min.segment.length = 0, segment.colour = "grey35",
                 segment.size = 0.4, colour = "grey15",
                 family = "sans", na.rm = TRUE) +
  geom_node_text(aes(label = lab_iso),
                 size = 2.4, lineheight = 0.85, vjust = 1.9,
                 colour = "grey30", family = "sans", na.rm = TRUE) +
  scale_size_continuous(range = c(2.5, 14), guide = "none") +
  scale_fill_discrete(na.value = "grey82", guide = "none") +
  scale_shape_manual(values = c("Paris" = 24, "Brussels" = 22,
                                "Paris & Brussels" = 23,
                                "Not an attacker" = 21),
                     drop = TRUE, name = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 4, fill = "grey60"))) +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.22))) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.08))) +
  coord_cartesian(clip = "off") +
  theme_graph(base_family = "sans", base_size = 10) +
  theme(plot.margin = margin(2, 2, 2, 2),
        legend.position = c(0.005, 0.005),
        legend.justification = c(0, 0),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 20, margin = margin(b = 1)),
        plot.subtitle = element_text(size = 10, colour = "grey35")) +
  labs(
    title = "The Islamic State network in Europe - community structure",
    subtitle = sprintf(
      "Louvain communities of %d or more (modularity %.2f).\nKey individuals named; a cross marks attackers who are dead. Node size = degree. %d unconnected actors below.",
      min_size, igraph::modularity(cl), sum(iso))
  )

print(p_clusters)
ggsave(ise_output_path("ise_clusters.pdf"), p_clusters, width = 12, height = 10,
       device = cairo_pdf)

###########################


## Membership listing
for (k in sort(unique(igraph::membership(cl)))) {
  cat("\n---", lab_for(k), "( n =", sum(igraph::membership(cl) == k), ")---\n")
  cat(paste(sort(igraph::V(gcon)$name[igraph::membership(cl) == k]),
            collapse = "\n"), "\n")
}
cat("\n--- unconnected (", sum(iso), ")---\n")
cat(paste(sort(igraph::V(g)$name[iso]), collapse = "\n"), "\n")


## ===========================================================================
## 5. Descriptives
## ===========================================================================

deg <- igraph::degree(g, mode = "all")

c(mean_degree     = mean(deg),
  median_degree   = median(deg),
  maximum_degree  = max(deg),
  isolated_actors = sum(deg == 0),
  density         = igraph::edge_density(g),
  transitivity    = igraph::transitivity(g))

degree_distribution_df <- data.frame(
  Degree     = 0:max(deg),
  Frequency  = tabulate(deg + 1, nbins = max(deg) + 1)
)
degree_distribution_df$Proportion <-
  degree_distribution_df$Frequency / igraph::vcount(g)

degree_distribution_df

plot(degree_distribution_df$Degree, degree_distribution_df$Frequency,
     type = "h", lwd = 4, lend = 1, col = "steelblue",
     xlab = "Degree", ylab = "Number of actors",
     main = "Degree distribution of the IS-E network", xaxt = "n")
axis(side = 1, at = seq(0, max(deg), by = 5))
points(degree_distribution_df$Degree, degree_distribution_df$Frequency,
       pch = 16, col = "steelblue")

## --- degree by member, horizontal bars --------------------------------------

member_degree <- sort(setNames(deg, igraph::V(g)$name), decreasing = TRUE)
plot_degree   <- rev(member_degree)

png(ise_output_path("IS-E_member_degrees.png"), width = 3000, height = 4800, res = 300)
par(mar = c(5, 15, 3, 3), mgp = c(3, 1, 0))
bar_positions <- barplot(
  plot_degree, names.arg = names(plot_degree), horiz = TRUE, las = 1,
  col = "steelblue", border = NA, cex.names = 0.55,
  xlim = c(0, max(plot_degree) + 5),
  xlab = "Degree", main = "Degree of IS-E network members"
)
text(x = plot_degree + 0.5, y = bar_positions, labels = plot_degree,
     pos = 4, cex = 0.55)
dev.off()


## ===========================================================================
## 6. How much of the network is coded rather than observed?
## ===========================================================================

logical_cols <- sapply(ties, is.logical)
colSums(ties[, logical_cols])

reported  <- sum(ties$linkedToAttacker | ties$linkedToWanted)
cat("dyads resting on a reported person-to-person link:",
    reported, "of", nrow(ties), "\n")


## ===========================================================================
## 7. Convert the igraph projections for ERGM estimation
## ===========================================================================

## The GML files share the same actors, but only the Standard graph was
## assigned the cleaned actor attributes above. Copy those attributes by name
## so nodematch("citizenship") is defined identically in every projection.
copy_vertex_attributes <- function(target, source,
                                   attributes = c("age", "citizenship",
                                                  "status", "role")) {
  idx <- match(igraph::V(target)$name, igraph::V(source)$name)

  if (anyNA(idx)) {
    stop(
      "Unmatched vertices while copying attributes: ",
      paste(igraph::V(target)$name[is.na(idx)], collapse = "; ")
    )
  }

  for (attribute in attributes) {
    target <- igraph::set_vertex_attr(
      target,
      name = attribute,
      value = igraph::vertex_attr(source, attribute)[idx]
    )
  }

  target
}

g_l <- copy_vertex_attributes(g_l, g)
g_e <- copy_vertex_attributes(g_e, g)

## Convert an igraph object to the statnet network class required by ergm.
igraph_to_statnet <- function(graph,
                              attributes = c("age", "citizenship",
                                             "status", "role")) {
  adjacency <- as.matrix(
    igraph::as_adjacency_matrix(graph, sparse = FALSE)
  )

  net <- network::network(
    adjacency,
    matrix.type = "adjacency",
    directed = igraph::is_directed(graph),
    loops = FALSE
  )

  network::set.vertex.attribute(
    net,
    "vertex.names",
    igraph::V(graph)$name
  )

  for (attribute in attributes) {
    network::set.vertex.attribute(
      net,
      attribute,
      igraph::vertex_attr(graph, attribute)
    )
  }

  net
}

net_Limited  <- igraph_to_statnet(g_l)
net_Standard <- igraph_to_statnet(g)
net_Extended <- igraph_to_statnet(g_e)

stopifnot(
  network::network.size(net_Limited) == igraph::vcount(g_l),
  network::network.edgecount(net_Limited) == igraph::ecount(g_l),
  network::network.size(net_Standard) == igraph::vcount(g),
  network::network.edgecount(net_Standard) == igraph::ecount(g),
  network::network.size(net_Extended) == igraph::vcount(g_e),
  network::network.edgecount(net_Extended) == igraph::ecount(g_e)
)



## ===========================================================================
## 8. M6: edges + isolates + citizenship homophily + GWESP(0.25)
## ===========================================================================
##
##   Pr(Y = y | theta, x) ∝ exp{ theta_E    s_E(y)
##                             + theta_iso  s_iso(y)
##                             + theta_C    s_C(y, x)
##                             + theta_GW   s_GWESP(y; 0.25) }

M6 <- "edges + isolates + nodematch('citizenship') + gwesp(0.25, fixed = TRUE)"

f_M6 <- function(net_name) {
  f <- stats::as.formula(paste(net_name, "~", M6))
  environment(f) <- globalenv()      # so simulate() and gof() resolve the LHS later
  f
}

ctrl_M6 <- control.ergm(
  MCMC.burnin     = 20000,
  MCMC.interval   = 2000,
  MCMC.samplesize = 5000,
  MCMLE.maxit     = 40,
  seed            = 20260909
)

fit_L <- ergm(f_M6("net_Limited"),  control = ctrl_M6)   # prior source
fit_S <- ergm(f_M6("net_Standard"), control = ctrl_M6)   # analysis network

summary(fit_L)
summary(fit_S)


## ---------------------------------------------------------------------------
## Fixed-decay sensitivity analysis
## ---------------------------------------------------------------------------

alpha_grid_ise <- c(0.1, 0.2, 0.25, 0.3, 0.4)
alpha_sensitivity <- sapply(alpha_grid_ise, function(a) {
  f <- as.formula(sprintf(
    "net_Standard ~ edges + isolates + nodematch('citizenship') + gwesp(%f, fixed = TRUE)", a))
  environment(f) <- globalenv()
  fit <- try(ergm(f, control = ctrl_M6), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(AIC = NA, BIC = NA, theta_gw = NA))
  cf <- coef(fit)
  c(AIC = AIC(fit), BIC = BIC(fit), theta_gw = unname(cf[length(cf)]))
})
colnames(alpha_sensitivity) <- alpha_grid_ise
round(alpha_sensitivity, 2)


form_M6_alpha_040 <- as.formula(
  "net_Standard ~ edges + isolates + nodematch('citizenship') + gwesp(0.4, fixed = TRUE)"
)
environment(form_M6_alpha_040) <- globalenv()
fit_S_alpha_040 <- ergm(form_M6_alpha_040, control = ctrl_M6)

for (fit in list(`0.25` = fit_S, `0.40` = fit_S_alpha_040)) {
  s <- simulate(fit, nsim = 200, output = "stats",
                control = control.simulate.ergm(MCMC.burnin = 20000,
                                                MCMC.interval = 2000))
  e <- as.numeric(s[, "edges"])
  cat(sprintf("median %5.0f | sd %6.1f | empty %4.1f%% | >50%% dens %4.1f%% | obs 197 at q %.2f\n",
              median(e), sd(e), 100*mean(e == 0), 100*mean(e > 1785), mean(e <= 197)))
}


# =========================================================================
# IS-E covert network: ERGM and BERGM analysis
# Reproducible version: one master seed, per-block derived seeds
#
# Data: Gutfraind & Genkin (2017), https://github.com/sashagutfraind/Nov13
# =========================================================================

library(igraph)
library(ergm)
library(Bergm)
library(intergraph)
library(coda)

# igraph and sna both export degree(), betweenness() and closeness(), with
# incompatible input types. Bergm loads sna, so whichever is attached last
# wins. Bind the igraph versions explicitly rather than relying on the
# attach order; the failure mode otherwise is an error about
# as.edgelist.sna, which appears nowhere in this script.

degree               <- igraph::degree
betweenness          <- igraph::betweenness
closeness            <- igraph::closeness
transitivity         <- igraph::transitivity
components           <- igraph::components
edge_density         <- igraph::edge_density
count_triangles      <- igraph::count_triangles
assortativity_degree <- igraph::assortativity_degree
mean_distance        <- igraph::mean_distance
diameter             <- igraph::diameter
induced_subgraph     <- igraph::induced_subgraph
as_edgelist          <- igraph::as_edgelist
sample_gnm           <- igraph::sample_gnm
layout_with_fr       <- igraph::layout_with_fr
norm_coords          <- igraph::norm_coords
read_graph           <- igraph::read_graph
permute              <- igraph::permute
vcount               <- igraph::vcount
ecount               <- igraph::ecount
V                    <- igraph::V


# =========================================================================
# Reproducibility
# =========================================================================
#
# As in Chapter 4: everything stochastic traces back to MASTER_SEED, and
# block_seed() turns a block's name into a seed so that blocks can be re-run
# independently without shifting one another's streams.
#
# control.ergm(seed=) is NOT used, deliberately. Passing it overrides the
# global stream and severs the chain from MASTER_SEED.

MASTER_SEED <- 14092026

block_seed <- function(label, master = MASTER_SEED) {
  h <- 5381
  for (ch in utf8ToInt(as.character(label))) {
    h <- (h * 33 + ch) %% 2147483647
  }
  as.integer((h + master * 7919) %% 2147483647)
}

set.seed(MASTER_SEED)

setwd("/Users/leowestland/Desktop/Stats MSc/MSc dissertation/Nov13-master")

OUT <- "ch6"
dir.create(OUT, showWarnings = FALSE)


# =========================================================================
# Data
# =========================================================================
#
# Three one-mode projections of the same multimodal database, nested by
# construction: N_limited subset N_standard subset N_extended. They differ in
# which relationship types are treated as generating a tie, and therefore in
# their tolerance for false-positive ties.

g   <- read_graph("ise_terrorNetwork.gml",         format = "gml")  # standard
g_e <- read_graph("ise_terrorNetworkExtended.gml", format = "gml")
g_l <- read_graph("ise_terrorNetworkLimited.gml",  format = "gml")

# The projections share a vertex set but not necessarily an ordering, and
# everything below compares them vertex by vertex.

ref <- V(g)$name
stopifnot(setequal(V(g_l)$name, ref), setequal(V(g_e)$name, ref))
g_l <- permute(g_l, match(V(g_l)$name, ref))
g_e <- permute(g_e, match(V(g_e)$name, ref))

people <- read.csv("ise_personDetails.csv", stringsAsFactors = FALSE,
                   fileEncoding = "latin1")
rels   <- read.csv("ise_relationships.csv", stringsAsFactors = FALSE,
                   fileEncoding = "latin1")

# Citizenship is unrecorded for nine actors, in two forms: NA for two and the
# literal string "unknown" for seven. Both must be caught. Each is given a
# distinct placeholder so that two actors of unrecorded citizenship are not
# counted as a match, since nothing in the data says they share one.

MISSING_CITIZENSHIP <- c("", "unknown", "na", "n/a", "-", "?")

attach_attrs <- function(x) {
  i <- match(V(x)$name, people$p.name)
  stopifnot(!any(is.na(i)))
  
  V(x)$age         <- people$p.age[i]
  V(x)$citizenship <- people$p.citizenship[i]
  V(x)$status      <- people$p.status[i]
  V(x)$role        <- people$p.role[i]
  
  V(x)$status[V(x)$status %in% "arreted"] <- "arrested"   # source typo
  V(x)$status[is.na(V(x)$status)]         <- "unknown"
  V(x)$role[is.na(V(x)$role)]             <- "none"
  
  miss <- is.na(V(x)$citizenship) |
    trimws(tolower(V(x)$citizenship)) %in% MISSING_CITIZENSHIP
  V(x)$citizenship[miss] <- paste0("unknown_", V(x)$name[miss])
  
  x
}

g   <- attach_attrs(g)
g_l <- attach_attrs(g_l)
g_e <- attach_attrs(g_e)

# ergm needs `network` objects, and resolves a formula's left-hand side
# lazily, so these must be visible globally.

net   <- intergraph::asNetwork(g)     # standard: the analysis network
net_l <- intergraph::asNetwork(g_l)
net_e <- intergraph::asNetwork(g_e)

n_nodes <- vcount(g)
n_dyads <- choose(n_nodes, 2)

# Guards. These halt rather than let a later block run on the wrong data;
# the citizenship masking in particular is easy to get silently wrong.

stopifnot(sum(grepl("^unknown_", V(g)$citizenship)) == 9)
stopifnot(summary(net   ~ edges)                    == 197,
          summary(net   ~ isolates)                 ==   8,
          summary(net   ~ nodematch("citizenship")) ==  61,
          summary(net_l ~ edges)                    == 163)

for (nm in c("net_l", "net", "net_e")) {
  x <- get(nm)
  cat(sprintf("%-6s edges %3d | isolates %2d | matched %2d\n", nm,
              summary(x ~ edges), summary(x ~ isolates),
              summary(x ~ nodematch("citizenship"))))
}


# ========================================
# Plot of the IS-E network
# ========================================
#
# Coloured by citizenship, sized by square root of degree so the degree-40
# hub does not swallow the figure. Ties between actors of the same
# citizenship are drawn darker, so the figure shows the homophily that the
# nodematch term measures.

set.seed(block_seed("ise_layout"))

cit <- V(g)$citizenship
cit[grepl("^unknown_", cit)] <- "unknown"

cit_tab  <- sort(table(cit), decreasing = TRUE)
cit_keep <- setdiff(names(cit_tab)[cit_tab >= 3], "unknown")

cit_grp <- factor(
  ifelse(cit %in% cit_keep, cit,
         ifelse(cit == "unknown", "unknown", "other")),
  levels = c(cit_keep, "other", "unknown")
)

cit_pal <- setNames(
  c(hcl.colors(length(cit_keep), "Dark 3"), "grey72", "grey88"),
  levels(cit_grp)
)

el      <- as_edgelist(g, names = FALSE)
matched <- cit[el[, 1]] == cit[el[, 2]] & cit[el[, 1]] != "unknown"

lay <- norm_coords(layout_with_fr(g, niter = 5000), -1, 1, -1, 1)

# Fruchterman-Reingold pushes the eight isolates to the frame edge, which is
# honest but unreadable. Park them in a row underneath instead.
iso <- which(degree(g) == 0)
lay[iso, 1] <- seq(-0.85, 0.85, length.out = length(iso))
lay[iso, 2] <- -1.28

while (!is.null(dev.list())) dev.off()
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))

plot(
  g,
  layout  = lay,
  rescale = FALSE,
  xlim    = c(-1.05, 1.05),
  ylim    = c(-1.4, 1.05),
  
  vertex.color       = cit_pal[as.character(cit_grp)],
  vertex.frame.color = "white",
  vertex.size        = 3.2 + 2.4 * sqrt(degree(g)),
  
  vertex.label       = ifelse(degree(g) >= 10, V(g)$name, NA),
  vertex.label.cex   = 0.55,
  vertex.label.color = "black",
  vertex.label.dist  = 0.9,
  
  edge.color = ifelse(matched, "grey25", "grey82"),
  edge.width = ifelse(matched, 1.5, 0.7),
  
  main = "IS-E network, standard projection"
)

legend("topleft", legend = levels(cit_grp), pt.bg = cit_pal, pch = 21,
       col = "white", pt.cex = 1.5, cex = 0.72, bty = "n",
       title = "Citizenship")

legend("topright", legend = c("shared citizenship", "different"),
       lwd = c(1.5, 0.7), col = c("grey25", "grey82"),
       cex = 0.72, bty = "n")

text(0, -1.38, "isolated actors", cex = 0.62, col = "grey35")


# ========================================
# Nesting of the three projections
# ========================================
#
# Verifies the claim that the projections are nested, and measures how much
# the limited and standard networks have in common. That last number is the
# honest measure of the data reuse discussed in Section 6.4.

edge_key <- function(x) {
  e <- as_edgelist(x, names = TRUE)
  paste(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2]), sep = " -- ")
}

E_l <- edge_key(g_l)
E_s <- edge_key(g)
E_e <- edge_key(g_e)

cat("limited not in standard (want 0):", length(setdiff(E_l, E_s)), "\n")
cat("standard not in extended (want 0):", length(setdiff(E_s, E_e)), "\n")
cat("added by site ties     (S \\ L):", length(setdiff(E_s, E_l)), "\n")
cat("added by locality ties (E \\ S):", length(setdiff(E_e, E_s)), "\n")
cat("limited as a fraction of standard:", round(length(E_l) / length(E_s), 3), "\n")


# ========================================
# Structural summaries
# ========================================

nets <- list(Limited = g_l, Standard = g, Extended = g_e)

summarise_net <- function(x) {
  d   <- degree(x)
  p   <- edge_density(x)
  cmp <- components(x)
  gc  <- induced_subgraph(x, cmp$membership == which.max(cmp$csize))
  
  c(Nodes            = vcount(x),
    Edges            = ecount(x),
    Density          = p,
    Isolates         = sum(d == 0),
    MeanDegree       = mean(d),
    MaxDegree        = max(d),
    Triangles        = sum(count_triangles(x)) / 3,
    Clustering       = transitivity(x, type = "global"),
    Assortativity    = assortativity_degree(x, directed = FALSE),
    Components       = cmp$no,
    LargestComponent = max(cmp$csize),
    MeanGeodesic     = mean_distance(gc, directed = FALSE),
    Diameter         = diameter(gc, directed = FALSE))
}

comp_tab <- t(vapply(nets, summarise_net, numeric(13)))
round(comp_tab, 3)


# ========================================
# Erdos-Renyi benchmark
# ========================================
#
# Every statistic above rises with density, so the raw table cannot separate
# "denser" from "differently shaped". Benchmarking each projection against
# G(n, m) on the same n and m removes the density effect.

set.seed(block_seed("er_benchmark"))

er_compare <- function(x, label, nsim = 1000) {
  n <- vcount(x)
  m <- ecount(x)
  
  obs <- c(Triangles     = sum(count_triangles(x)) / 3,
           Clustering    = transitivity(x, type = "global"),
           DegreeSD      = sd(degree(x)),
           Isolates      = sum(degree(x) == 0),
           Assortativity = assortativity_degree(x, directed = FALSE))
  
  sims <- replicate(nsim, {
    r <- sample_gnm(n, m, directed = FALSE)
    c(sum(count_triangles(r)) / 3,
      transitivity(r, type = "global"),
      sd(degree(r)),
      sum(degree(r) == 0),
      assortativity_degree(r, directed = FALSE))
  })
  
  data.frame(network  = label,
             stat     = names(obs),
             observed = obs,
             er_mean  = rowMeans(sims),
             er_sd    = apply(sims, 1, sd),
             z        = (obs - rowMeans(sims)) / apply(sims, 1, sd),
             row.names = NULL)
}

er_tab <- do.call(rbind, Map(er_compare, nets, names(nets)))
print(er_tab, digits = 3, row.names = FALSE)

# Analytic expected isolate count under G(n, p): E[s_iso] = n (1-p)^(n-1).
# Quote this in the text rather than the simulated value, since it is a
# closed form the reader can check.

for (nm in names(nets)) {
  x <- nets[[nm]]
  p <- edge_density(x)
  cat(sprintf("%-9s p = %.4f  expected isolates %.2f  observed %d\n",
              nm, p, vcount(x) * (1 - p)^(vcount(x) - 1), sum(degree(x) == 0)))
}


# ========================================
# Degree and ESP distributions
# ========================================

nws     <- list(Limited = net_l, Standard = net, Extended = net_e)
max_deg <- max(vapply(nets, function(x) max(degree(x)), numeric(1)))

deg_tab <- vapply(nws, function(x) summary(x ~ degree(0:max_deg)),
                  numeric(max_deg + 1))
esp_tab <- vapply(nws, function(x) summary(x ~ esp(0:10)), numeric(11))

deg_tab[1:16, ]
esp_tab

# The spike at esp7 is the nine-actor leadership clique: C(9,2) = 36 ties,
# each with exactly seven shared partners, present in all three projections.
# The projection rule turns any shared event into a clique, so part of the
# observed clustering is generated by the coding rather than observed.

igraph::clique_num(g)
lapply(igraph::largest_cliques(g), function(v) V(g)$name[v])


# ========================================
# Degree distribution against Erdos-Renyi
# ========================================
#
# The Chapter 6 analogue of Figure 4.6. Three layers: observed counts, the
# exact Binomial(n-1, p) expectation, and a 95% band from simulated G(n, m)
# graphs. The band matters because the point expectation alone does not show
# what sampling variation would produce.

set.seed(block_seed("degree_vs_er"))

KMAX <- 15
p    <- edge_density(g)

sim_deg <- replicate(1000, {
  tabulate(degree(sample_gnm(n_nodes, ecount(g), directed = FALSE)) + 1,
           nbins = KMAX + 1)
})

env      <- apply(sim_deg, 1, quantile, c(0.025, 0.975))
obs_deg  <- deg_tab[1:(KMAX + 1), "Standard"]
exp_deg  <- n_nodes * dbinom(0:KMAX, n_nodes - 1, p)
tail_deg <- sort(degree(g)[degree(g) > KMAX])

par(mfrow = c(1, 1), mar = c(5.1, 5.1, 4.1, 2.1), mgp = c(3.5, 1, 0))

plot(0:KMAX, obs_deg, type = "n",
     xlab = "Degree", ylab = "Number of actors",
     main = "Degree distribution, standard projection",
     ylim = c(0, max(obs_deg, env) * 1.05))

polygon(c(0:KMAX, rev(0:KMAX)), c(env[1, ], rev(env[2, ])),
        col = "grey88", border = NA)

points(0:KMAX, obs_deg, type = "h", lwd = 6, col = "grey45", lend = 1)

lines(0:KMAX, exp_deg, lty = 2)
points(0:KMAX, exp_deg, pch = 1, cex = 0.8)

legend("topright", bty = "n", cex = 0.8,
       legend = c("Observed",
                  expression(paste("Expected under ", G(n, p))),
                  "95% band, G(n, m)"),
       lwd = c(6, 1, 6), lty = c(1, 2, 1),
       col = c("grey45", "black", "grey88"))

mtext(sprintf("%d actors above degree %d (%s)",
              length(tail_deg), KMAX, paste(tail_deg, collapse = ", ")),
      side = 3, line = -1.2, cex = 0.72, col = "grey35")


# ========================================
# Citizenship homophily, before any model
# ========================================
#
# The nodematch count means little on its own. What matters is the comparison
# with the proportion of dyads that would share a citizenship under random
# mixing, which is what the coefficient is measuring against.

chance_match <- sum(choose(as.numeric(table(V(g)$citizenship)), 2)) / n_dyads
obs_match    <- summary(net ~ nodematch("citizenship")) / summary(net ~ edges)

cat(sprintf("matched dyads under random mixing: %.4f\n", chance_match))
cat(sprintf("matched ties observed:             %.4f\n", obs_match))
cat(sprintf("odds ratio: %.3f\n",
            (obs_match / (1 - obs_match)) / (chance_match / (1 - chance_match))))

sort(table(V(g)$citizenship[!grepl("^unknown_", V(g)$citizenship)]),
     decreasing = TRUE)


# =========================================================================
# Model M6
# =========================================================================
#
#   Pr(Y = y | theta, x) propto exp{ theta_E   s_E(y)
#                                  + theta_iso s_iso(y)
#                                  + theta_C   s_C(y, x)
#                                  + theta_GW  s_GWESP(y; 0.25) }
#
# Each term answers a feature of the descriptive analysis: the edge count
# controls density; the isolate count reproduces the empty tail of the degree
# distribution, which is far heavier than any homogeneous model can generate;
# nodematch captures citizenship homophily, citizenship being the only
# recorded attribute plausibly exogenous to tie formation; and GWESP captures
# clustering while discounting each additional shared partner.
#
# The decay is fixed rather than estimated jointly, both to avoid the
# identification difficulties of a curved family and because the Bayesian
# estimation below requires a linear exponential family.

ALPHA <- 0.25

M6 <- sprintf(
  "edges + isolates + nodematch('citizenship') + gwesp(%s, fixed = TRUE)",
  ALPHA
)

f_M6 <- function(net_name) {
  f <- stats::as.formula(paste(net_name, "~", M6))
  environment(f) <- globalenv()
  f
}

CTRL <- control.ergm(
  MCMC.burnin     = 20000,
  MCMC.interval   = 2000,
  MCMC.samplesize = 5000,
  MCMLE.maxit     = 40
)

TERMS <- c("edges", "isolates", "nodematch.citizenship", "gwesp.fixed.0.25")
OBS   <- summary(f_M6("net"))
OBS


# ========================================
# Frequentist fit
# ========================================

set.seed(block_seed("ergm_m6_limited"))
fit_l <- ergm(f_M6("net_l"), control = CTRL)

set.seed(block_seed("ergm_m6_standard"))
fit <- ergm(f_M6("net"), control = CTRL)

summary(fit_l)
summary(fit)

mle_values <- coef(fit)
mle_values

coef_tab <- data.frame(
  term  = names(coef(fit_l)),
  est_L = coef(fit_l),
  se_L  = sqrt(diag(vcov(fit_l))),
  est_S = coef(fit),
  se_S  = sqrt(diag(vcov(fit))),
  row.names = NULL
)

print(coef_tab, digits = 3, row.names = FALSE)

# The two fits are on different networks, so their information criteria are
# not comparable. AIC is used below only within a projection.

c(AIC_limited = AIC(fit_l), AIC_standard = AIC(fit),
  BIC_limited = BIC(fit_l), BIC_standard = BIC(fit))

# A specification omitting the clustering term, for the homophily comparison
# in the text: the citizenship coefficient falls once closure is controlled.

set.seed(block_seed("ergm_no_gwesp"))
fit_nogw <- ergm(net ~ edges + nodematch("citizenship"))
coef(fit_nogw)
exp(coef(fit_nogw)["nodematch.citizenship"])


# ========================================
# Is the frequentist fit sound?
# ========================================
#
# MCMLE reporting convergence says the estimating equations are satisfied. It
# says nothing about whether the fitted model concentrates on empty or
# complete graphs. Simulate and look. At the MLE the simulated statistic
# means must also equal the observed values.

set.seed(block_seed("sim_check"))

sim_check <- function(object, obs, label, nsim = 200) {
  s <- simulate(object, nsim = nsim, output = "stats",
                control = control.simulate.ergm(MCMC.burnin   = 20000,
                                                MCMC.interval = 2000))
  e <- as.numeric(s[, "edges"])
  cat(sprintf("%-9s mean %5.0f | median %5.0f | sd %5.1f | empty %4.1f%% | >50%% dens %4.1f%% | obs %3d at q %.2f\n",
              label, mean(e), median(e), sd(e),
              100 * mean(e == 0), 100 * mean(e > 0.5 * n_dyads),
              obs, mean(e <= obs)))
  invisible(s)
}

s_l <- sim_check(fit_l, 163, "Limited")
s_s <- sim_check(fit,   197, "Standard")

round(colMeans(as.matrix(s_l)), 1)   # compare with summary(f_M6("net_l"))
round(colMeans(as.matrix(s_s)), 1)   # compare with OBS

pdf(file.path(OUT, "m6_mcmc_diagnostics.pdf"), width = 8, height = 10)
mcmc.diagnostics(fit_l)
mcmc.diagnostics(fit)
dev.off()


# ========================================
# Selecting the GWESP decay
# ========================================
#
# AIC alone cannot settle this. As alpha grows, 1-(1-e^-a)^k -> k e^-a, so the
# statistic tends towards 3 * (triangle count), the unstable term the
# geometric weighting exists to replace. Fit therefore improves right up to
# the point of degeneracy, and each value must also be checked by simulation.
#
# Slow: roughly 40 minutes.

set.seed(block_seed("gwesp_alpha_grid"))

ALPHAS <- seq(0.05, 0.50, by = 0.05)

grid_one <- function(a, nsim = 300) {
  f <- stats::as.formula(sprintf(
    "net ~ edges + isolates + nodematch('citizenship') + gwesp(%f, fixed = TRUE)", a))
  environment(f) <- globalenv()
  
  cat("alpha =", a, "\n")
  fit_a <- try(ergm(f, control = CTRL), silent = TRUE)
  
  if (inherits(fit_a, "try-error")) {
    return(data.frame(alpha = a, converged = FALSE, AIC = NA, BIC = NA,
                      theta_E = NA, theta_gw = NA, sim_mean = NA,
                      sim_sd = NA, empty = NA, dense = NA, obs_q = NA))
  }
  
  s  <- simulate(fit_a, nsim = nsim, output = "stats",
                 control = control.simulate.ergm(MCMC.burnin   = 20000,
                                                 MCMC.interval = 2000))
  e  <- as.numeric(s[, "edges"])
  cf <- coef(fit_a)
  
  data.frame(alpha      = a,
             converged  = TRUE,
             AIC        = AIC(fit_a),
             BIC        = BIC(fit_a),
             theta_E    = unname(cf["edges"]),
             theta_gw   = unname(cf[length(cf)]),
             sim_mean   = mean(e),
             sim_sd     = sd(e),
             empty      = 100 * mean(e == 0),
             dense      = 100 * mean(e > 0.5 * n_dyads),
             obs_q      = mean(e <= 197))
}

alpha_grid <- do.call(rbind, lapply(ALPHAS, grid_one))
print(alpha_grid, digits = 4, row.names = FALSE)

setNames(round(alpha_grid$sim_sd, 1), alpha_grid$alpha)
setNames(round(alpha_grid$obs_q,  2), alpha_grid$alpha)


# ========================================
# GWESP decay figure
# ========================================
#
# Stacked panels sharing the x axis, rather than twin y axes on one panel: an
# overlay produces a crossing point whose position depends only on the choice
# of scales, and it lands directly over the chosen value.

ok_a   <- alpha_grid$converged & !is.na(alpha_grid$AIC)
fail_a <- alpha_grid$alpha[!ok_a]
xr_a   <- range(alpha_grid$alpha)

par(mfrow = c(2, 1))

par(mar = c(0.4, 4.6, 2.6, 1.2))
plot(alpha_grid$alpha[ok_a], alpha_grid$AIC[ok_a], type = "b", pch = 16,
     xlim = xr_a, xaxt = "n", xlab = "", ylab = "AIC",
     main = "Selecting the GWESP decay")
abline(v = ALPHA, col = "red", lty = 2, lwd = 2)
if (length(fail_a)) {
  abline(v = min(fail_a), col = "grey55", lty = 3, lwd = 2)
  text(min(fail_a) - 0.018 * diff(par("usr")[1:2]),
       par("usr")[4] - 0.12 * diff(par("usr")[3:4]),
       "no convergence", srt = 90, adj = 1, cex = 0.72, col = "grey40")
}

par(mar = c(4.4, 4.6, 0.6, 1.2))
plot(alpha_grid$alpha[ok_a], alpha_grid$sim_sd[ok_a], type = "b", pch = 16,
     xlim = xr_a, xlab = expression(alpha), ylab = "S.d. of edge count")
abline(v = ALPHA, col = "red", lty = 2, lwd = 2)
if (length(fail_a)) abline(v = min(fail_a), col = "grey55", lty = 3, lwd = 2)


# =========================================================================
# Elicited priors
# =========================================================================
#
# Built from what is known about covert networks, independently of the IS-E
# data: sparse, agnostic about isolates, sceptical about citizenship as a
# proxy for membership, confident that triadic closure is present.
#
# Note that theta_E is NOT the logit density. It is the conditional log-odds
# of a tie between two actors sharing no partners, which is far lower, since
# GWESP absorbs most of the tie-formation propensity. The prior is centred on
# the midpoint of the two frequentist estimates rather than on any density.
#
# The isolate and closure priors were revised after the prior predictive
# check below; the diagnostic is retained because it is the evidence for the
# revision.

dlaplace <- function(x, mu = 0, b = 1, log = FALSE) {
  d <- -base::log(2 * b) - abs(x - mu) / b
  if (log) d else exp(d)
}

rlaplace <- function(n, mu = 0, b = 1) {
  mu + b * (rexp(n) - rexp(n))
}

PRIOR_E   <- list(mean = mean(c(coef(fit_l)["edges"], coef(fit)["edges"])),
                  sd = 0.5)                      # -6.696
PRIOR_ISO <- list(min = -5, max = 5)
PRIOR_C   <- list(mu = 0, b = 0.5)
PRIOR_GW  <- list(shape = 4, rate = 1)

PRIOR_DENS <- list(
  edges                 = function(x) dnorm(x, PRIOR_E$mean, PRIOR_E$sd),
  isolates              = function(x) dunif(x, PRIOR_ISO$min, PRIOR_ISO$max),
  nodematch.citizenship = function(x) dlaplace(x, PRIOR_C$mu, PRIOR_C$b),
  gwesp.fixed.0.25      = function(x) dgamma(x, PRIOR_GW$shape, PRIOR_GW$rate)
)

log_prior <- function(th) {
  dnorm(th[1], PRIOR_E$mean, PRIOR_E$sd, log = TRUE) +
    dunif(th[2], PRIOR_ISO$min, PRIOR_ISO$max, log = TRUE) +
    dlaplace(th[3], PRIOR_C$mu, PRIOR_C$b, log = TRUE) +
    dgamma(th[4], PRIOR_GW$shape, PRIOR_GW$rate, log = TRUE)
}

r_prior <- function(K) {
  cbind(edges                 = rnorm(K, PRIOR_E$mean, PRIOR_E$sd),
        isolates              = runif(K, PRIOR_ISO$min, PRIOR_ISO$max),
        nodematch.citizenship = rlaplace(K, PRIOR_C$mu, PRIOR_C$b),
        gwesp.fixed.0.25      = rgamma(K, PRIOR_GW$shape, PRIOR_GW$rate))
}

PRIOR_MODE <- c(edges                 = PRIOR_E$mean,
                isolates              = NA,
                nodematch.citizenship = PRIOR_C$mu,
                gwesp.fixed.0.25      = (PRIOR_GW$shape - 1) / PRIOR_GW$rate)

PRIOR_TEXT <- c(edges                 = sprintf("N(%.2f, %.2f^2)", PRIOR_E$mean, PRIOR_E$sd),
                isolates              = "Uniform(-3, 3)",
                nodematch.citizenship = "Laplace(0, 0.5)",
                gwesp.fixed.0.25      = "Gamma(9, 2.5)")

# Summaries for the write-up, computed rather than retyped.
c(mean = PRIOR_GW$shape / PRIOR_GW$rate,
  mode = (PRIOR_GW$shape - 1) / PRIOR_GW$rate,
  sd   = sqrt(PRIOR_GW$shape) / PRIOR_GW$rate,
  q025 = qgamma(0.025, PRIOR_GW$shape, PRIOR_GW$rate),
  q975 = qgamma(0.975, PRIOR_GW$shape, PRIOR_GW$rate))

c(q025 = PRIOR_C$b * log(0.05), q975 = -PRIOR_C$b * log(0.05))

set.seed(block_seed("prior_marginals"))
round(apply(r_prior(50000), 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975)), 3)


# ========================================
# Prior predictive edge-count simulation
# ========================================
#
# Draw theta from the prior, simulate a network from M6 at that theta, and
# record its statistics. A usable prior implies networks that bracket the
# observed one: the observed statistics should sit inside the prior
# predictive distribution rather than out in a tail.

set.seed(block_seed("prior_pred_ergm"))

S <- 6000

theta_prior <- r_prior(S)

prior_pred_stats <- matrix(NA_real_, S, 4, dimnames = list(NULL, TERMS))

for (s in seq_len(S)) {
  
  simulated_stats <- tryCatch(
    simulate(
      f_M6("net"),
      coef    = unname(theta_prior[s, ]),
      nsim    = 1,
      output  = "stats",
      control = control.simulate.formula(
        MCMC.burnin   = 5000,
        MCMC.interval = 1
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(simulated_stats)) {
    prior_pred_stats[s, ] <- as.numeric(simulated_stats)
  }
  
  if (s %% 100 == 0) {
    message("Completed ", s, " of ", S)
  }
}

n_failed <- sum(is.na(prior_pred_stats[, "edges"]))
n_failed

valid_prior_edges <- prior_pred_stats[!is.na(prior_pred_stats[, "edges"]),
                                      "edges"]

for (k in TERMS) {
  v <- prior_pred_stats[, k]
  v <- v[!is.na(v)]
  cat(sprintf("%-22s median %8.1f  [%8.1f, %8.1f]   observed %7.1f at q %.2f\n",
              k, median(v), quantile(v, .025), quantile(v, .975),
              OBS[k], mean(v <= OBS[k])))
}

par(mfrow = c(1, 1), mar = c(5.1, 5.1, 4.1, 2.1), mgp = c(3.5, 1, 0))

hist(
  valid_prior_edges,
  breaks      = seq(-5, n_dyads + 10, by = 10),
  probability = TRUE,
  col         = "grey80",
  border      = "white",
  xlim        = c(0, 1000),
  main        = "Prior predictive edge count",
  xlab        = "Number of edges",
  ylab        = "Proportion of simulated graphs"
)

abline(v = 197, col = "black", lwd = 2, lty = 4)

legend("topright", legend = "Observed = 197",
       col = "black", lty = 4, lwd = 2, bty = "n")

# Empty or complete graphs
mean(valid_prior_edges %in% c(0, n_dyads))

# Nearly empty or nearly complete graphs
mean(valid_prior_edges <= 2 | valid_prior_edges >= n_dyads - 2)

summary(valid_prior_edges)


# ========================================
# Where the degenerate draws come from
# ========================================
#
# Diagnostic, not tuning. Regressing the empty-graph indicator on the sampled
# parameter vector identifies which regions of the prior are responsible. On
# the initial specification, with Uniform(-5, 5) on theta_iso and Gamma(4, 1)
# on theta_GWESP, this showed the two acting jointly: with the edge parameter
# near -6.7 a tie between actors sharing no partners has conditional
# probability below 0.002, so the network is sustained only by triadic
# closure, and draws combining weak closure with a positive isolate parameter
# collapse to the empty graph, which the isolate term then reinforces.

empty <- !is.na(prior_pred_stats[, "edges"]) & prior_pred_stats[, "edges"] == 0
mean(empty)

if (sum(empty) >= 5) {
  round(rbind(empty    = colMeans(theta_prior[empty, ]),
              nonempty = colMeans(theta_prior[!empty, ])), 3)
  
  summary(glm(empty ~ ., family = binomial,
              data = data.frame(theta_prior, empty = empty)))
}


# =========================================================================
# BERGM fit
# =========================================================================
#
# The exchange algorithm needs only log pi(theta') - log pi(theta), so a
# non-normal prior requires no change to the sampler beyond replacing that
# one function. See Annex E.

set.seed(block_seed("bergm_m6"))

MI <- 20000

bergm_fit <- bergm(
  f_M6("net"),
  log.prior  = log_prior,        # EDIT to match your modified interface
  burn.in    = 500,
  main.iters = MI,
  aux.iters  = 3000,
  nchains    = 4,
  gamma      = 0.4
)

bergm_draws <- as.matrix(bergm_fit$Theta)
colnames(bergm_draws) <- TERMS

bergm_fit$AR                     # target roughly 0.2 to 0.45
nrow(bergm_draws) / MI           # chains

theta_E_post   <- bergm_draws[, "edges"]
theta_ISO_post <- bergm_draws[, "isolates"]
theta_C_post   <- bergm_draws[, "nodematch.citizenship"]
theta_GW_post  <- bergm_draws[, "gwesp.fixed.0.25"]

effectiveSize(as.mcmc(bergm_draws))

chains <- lapply(seq_len(4), function(i) {
  as.mcmc(bergm_draws[((i - 1) * MI + 1):(i * MI), , drop = FALSE])
})
effectiveSize(as.mcmc.list(chains))
gelman.diag(as.mcmc.list(chains))

cor(bergm_draws)

par(mfrow = c(2, 2), mar = c(4, 4.2, 2.4, 1))
for (k in TERMS) {
  acf(bergm_draws[1:MI, k], lag.max = 200, main = k)
}


# post_density / post_mode are defined in the posterior-summaries block
# below; if this runs first, define them here instead.

post_density <- function(draws, lower = NULL, upper = NULL,
                         n = 10000, bw = "nrd0") {
  density(draws,
          from = if (is.null(lower)) min(draws) else lower,
          to   = if (is.null(upper)) max(draws) else upper,
          n    = n, bw = bw)
}

post_mode <- function(draws, ...) {
  d <- post_density(draws, ...)
  d$x[which.max(d$y)]
}

# theta_GWESP has positive support under its Gamma prior, so its density
# grid starts at zero. The others are unbounded.
LOWER <- list(edges = NULL, isolates = NULL,
              nodematch.citizenship = NULL, gwesp.fixed.0.25 = 0)

bergm_summary <- data.frame(
  term   = TERMS,
  mean   = colMeans(bergm_draws),
  median = apply(bergm_draws, 2, median),
  mode   = sapply(TERMS, function(k) post_mode(bergm_draws[, k], LOWER[[k]])),
  sd     = apply(bergm_draws, 2, sd),
  q2.5   = apply(bergm_draws, 2, quantile, 0.025),
  q97.5  = apply(bergm_draws, 2, quantile, 0.975),
  ESS    = effectiveSize(as.mcmc(bergm_draws)),
  MLE    = as.numeric(mle_values),
  row.names = NULL
)

bergm_summary$MCSE      <- bergm_summary$sd / sqrt(bergm_summary$ESS)
bergm_summary$MLE_in_CI <- with(bergm_summary, MLE > q2.5 & MLE < q97.5)

round(bergm_summary[, -1], 4)
cat("acceptance rate:", round(bergm_fit$AR, 3), "\n")


# ========================================
# Chain diagnostics
# ========================================

par(mfrow = c(4, 2), mar = c(4, 4.2, 2.2, 1), mgp = c(2.6, 0.8, 0))

nch <- nrow(bergm_draws) / MI

for (k in TERMS) {
  matplot(matrix(bergm_draws[, k], nrow = MI), type = "l", lty = 1,
          col = adjustcolor(1:nch, 0.6),
          xlab = "Iteration", ylab = k, main = paste(k, "trace"))
  plot(density(bergm_draws[, k]), main = paste(k, "density"),
       xlab = k, lwd = 2)
}


# ========================================
# Posterior summaries, with modes
# ========================================
#
# Every mode in the write-up comes from this one estimator. If the table and
# the figure call density() separately with different arguments, the mode
# line can end up off the plotted peak.
#
# lower / upper set the evaluation grid, NOT the kernel. R's density() uses a
# Gaussian kernel, so a little mass leaks past a hard boundary and the density
# is slightly under-estimated there. This matters only for theta_GWESP, whose
# Gamma prior gives it positive support; the posterior is many standard
# deviations from zero, so the effect is negligible, but say so in a footnote.

post_density <- function(draws, lower = NULL, upper = NULL,
                         n = 10000, bw = "nrd0") {
  density(draws,
          from = if (is.null(lower)) min(draws) else lower,
          to   = if (is.null(upper)) max(draws) else upper,
          n    = n,
          bw   = bw)
}

post_mode <- function(draws, ...) {
  d <- post_density(draws, ...)
  d$x[which.max(d$y)]
}

posterior_summary <- function(draws, lower = NULL, upper = NULL) {
  c(mean   = mean(draws),
    median = median(draws),
    mode   = post_mode(draws, lower = lower, upper = upper),
    sd     = sd(draws),
    q2.5   = unname(quantile(draws, 0.025)),
    q97.5  = unname(quantile(draws, 0.975)))
}

post_tab <- rbind(
  theta_E   = posterior_summary(theta_E_post),
  theta_ISO = posterior_summary(theta_ISO_post),
  theta_C   = posterior_summary(theta_C_post),
  theta_GW  = posterior_summary(theta_GW_post, lower = 0)
)

round(post_tab, 3)

# Prior mode, posterior mode and MLE side by side: how far the data moved
# each parameter away from what was believed beforehand.
round(cbind(prior_mode = PRIOR_MODE,
            post_mode  = post_tab[, "mode"],
            MLE        = mle_values,
            shift      = post_tab[, "mode"] - PRIOR_MODE), 3)

# Does the 95% credible interval contain the MLE?
mle_values > post_tab[, "q2.5"] & mle_values < post_tab[, "q97.5"]


# ========================================
# How stable is the mode?
# ========================================
#
# The mean and median are order statistics of the draws; the mode is not.
# Two cheap checks, both worth a sentence in the text.

mode_sensitivity <- function(draws, lower = NULL, upper = NULL) {
  sapply(c("nrd0", "nrd", "SJ"), function(b) {
    tryCatch(
      suppressWarnings(post_mode(draws, lower = lower, upper = upper, bw = b)),
      error = function(e) NA_real_
    )
  })
}

rbind(
  theta_E   = mode_sensitivity(theta_E_post),
  theta_ISO = mode_sensitivity(theta_ISO_post),
  theta_C   = mode_sensitivity(theta_C_post),
  theta_GW  = mode_sensitivity(theta_GW_post, lower = 0)
)

mode_by_chain <- function(draws, main.iters = MI, lower = NULL, upper = NULL) {
  
  nchains <- length(draws) / main.iters
  
  if (nchains != round(nchains)) {
    stop("main.iters does not match the stored draws.")
  }
  
  chain_id <- rep(seq_len(nchains), each = main.iters)
  
  tapply(draws, chain_id, function(z) {
    post_mode(z, lower = lower, upper = upper)
  })
}

rbind(
  theta_E   = mode_by_chain(theta_E_post),
  theta_ISO = mode_by_chain(theta_ISO_post),
  theta_C   = mode_by_chain(theta_C_post),
  theta_GW  = mode_by_chain(theta_GW_post, lower = 0)
)


# ========================================
# Prior and posterior figure
# ========================================

plot_prior_posterior <- function(
    draws,
    prior_function,
    mle_value,
    xlim,
    ylim = NULL,
    colour,
    x_label,
    show_ylabel = TRUE,
    lower = NULL,
    upper = NULL,
    show_mode = TRUE
) {
  
  # Ensure the MLE is inside the plotting range
  xlim <- range(c(xlim, mle_value))
  
  x <- seq(xlim[1], xlim[2], length.out = 1000)
  
  posterior_density <- post_density(draws, lower = lower, upper = upper)
  
  mode_value <- posterior_density$x[which.max(posterior_density$y)]
  
  prior_density <- prior_function(x)
  
  if (is.null(ylim)) {
    ylim <- c(0, 1.08 * max(posterior_density$y, prior_density, na.rm = TRUE))
  }
  
  plot(
    posterior_density$x,
    posterior_density$y,
    type = "l",
    col  = colour,
    lwd  = 2,
    lty  = 1,
    xlim = xlim,
    ylim = ylim,
    xlab = x_label,
    ylab = if (show_ylabel) "Density" else "",
    font.main = 2
  )
  
  # Prior
  lines(x, prior_density, col = colour, lwd = 2, lty = 2)
  
  # Posterior mode
  if (show_mode) {
    segments(
      x0 = mode_value, y0 = 0,
      x1 = mode_value, y1 = max(posterior_density$y),
      col = colour, lwd = 2, lty = 3
    )
  }
  
  # Posterior density at the MLE
  mle_height <- approx(
    posterior_density$x,
    posterior_density$y,
    xout = mle_value,
    rule = 3
  )$y
  
  # MLE reference line
  segments(
    x0 = mle_value, y0 = 0,
    x1 = mle_value, y1 = mle_height,
    col = "black", lwd = 2, lty = 4
  )
  
  legend(
    "topright",
    legend = c(
      paste0("MLE \u2248 ", round(mle_value, 3)),
      if (show_mode) paste0("MAP \u2248 ", round(mode_value, 3)),
      "Posterior",
      "Prior"
    ),
    col = c("black", if (show_mode) colour, colour, colour),
    lty = c(4, if (show_mode) 3, 1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.8
  )
  
  invisible(mode_value)
}

par(
  mfrow = c(2, 2),
  oma = c(0, 0, 3, 0),
  mar = c(4, 4, 0.4, 1) + 0.1,
  mgp = c(2.6, 1, 0)
)

fig_modes <- numeric(4)
names(fig_modes) <- TERMS

fig_modes["edges"] <- plot_prior_posterior(
  draws = theta_E_post,
  prior_function = PRIOR_DENS$edges,
  mle_value = unname(mle_values["edges"]),
  xlim = range(c(theta_E_post, PRIOR_E$mean + c(-3, 3) * PRIOR_E$sd)),
  colour = "red",
  x_label = expression(theta[E]),
  show_ylabel = TRUE
)

fig_modes["isolates"] <- plot_prior_posterior(
  draws = theta_ISO_post,
  prior_function = PRIOR_DENS$isolates,
  mle_value = unname(mle_values["isolates"]),
  xlim = range(c(theta_ISO_post, -5.1, 5.1)),
  colour = "darkgreen",
  x_label = expression(theta[iso]),
  show_ylabel = FALSE
)

fig_modes["nodematch.citizenship"] <- plot_prior_posterior(
  draws = theta_C_post,
  prior_function = PRIOR_DENS$nodematch.citizenship,
  mle_value = unname(mle_values["nodematch.citizenship"]),
  xlim = range(c(theta_C_post, -2, 2)),
  colour = "purple",
  x_label = expression(theta[C]),
  show_ylabel = TRUE
)

fig_modes["gwesp.fixed.0.25"] <- plot_prior_posterior(
  draws = theta_GW_post,
  prior_function = PRIOR_DENS$gwesp.fixed.0.25,
  mle_value = unname(mle_values["gwesp.fixed.0.25"]),
  xlim = range(c(theta_GW_post, 0, 8)),
  colour = "blue",
  x_label = expression(theta[GWESP]),
  show_ylabel = FALSE,
  lower = 0
)

mtext(
  "Priors and posteriors, IS-E standard projection",
  side = 3, outer = TRUE, line = 1, cex = 1.3, font = 2
)

# Consistency check: these must match the mode column of post_tab exactly.
cbind(figure = fig_modes, table = post_tab[, "mode"])


# ========================================
# Posterior predictive edge-count simulation
# ========================================

set.seed(block_seed("post_pred_ergm"))

S_post <- min(6000, nrow(bergm_draws))

draw_index <- sample(seq_len(nrow(bergm_draws)), size = S_post, replace = FALSE)

theta_post <- bergm_draws[draw_index, , drop = FALSE]

posterior_pred_edges <- rep(NA_real_, S_post)

for (s in seq_len(S_post)) {
  
  simulated_stats <- tryCatch(
    simulate(
      f_M6("net"),
      coef    = unname(theta_post[s, ]),
      nsim    = 1,
      output  = "stats",
      control = control.simulate.formula(
        MCMC.burnin   = 5000,
        MCMC.interval = 1
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(simulated_stats)) {
    posterior_pred_edges[s] <- as.matrix(simulated_stats)[1, "edges"]
  }
  
  if (s %% 100 == 0) {
    message("Completed ", s, " of ", S_post)
  }
}

sum(is.na(posterior_pred_edges))

valid_post_edges <- posterior_pred_edges[!is.na(posterior_pred_edges)]

cat(sprintf("posterior predictive: median %.0f  [%.0f, %.0f]  observed 197 at q %.2f\n",
            median(valid_post_edges),
            quantile(valid_post_edges, .025),
            quantile(valid_post_edges, .975),
            mean(valid_post_edges <= 197)))

par(mfrow = c(1, 1), mar = c(5.1, 5.1, 4.1, 2.1), mgp = c(3.5, 1, 0))

hist(
  valid_post_edges,
  breaks      = seq(-2.5, n_dyads + 5, by = 5),
  probability = TRUE,
  col         = "grey70",
  border      = "white",
  xlim        = c(0, max(400, quantile(valid_post_edges, 0.995))),
  main        = "Posterior predictive edge count",
  xlab        = "Number of edges",
  ylab        = "Proportion of simulated graphs"
)

abline(v = 197, col = "black", lwd = 2, lty = 4)

legend("topright", legend = "Observed = 197",
       col = "black", lty = 4, lwd = 2, bty = "n")


# ========================================
# Prior and posterior predictive overlay
# ========================================
#
# The figure that shows what the data contributed: the prior predictive is
# broad and includes empty graphs, the posterior predictive concentrates on
# the observed density.

edge_breaks <- seq(-5, n_dyads + 10, by = 10)

prior_hist <- hist(valid_prior_edges, breaks = edge_breaks,
                   probability = TRUE, plot = FALSE)

post_hist <- hist(valid_post_edges, breaks = edge_breaks,
                  probability = TRUE, plot = FALSE)

ymax <- 1.1 * max(prior_hist$density, post_hist$density, na.rm = TRUE)

col_prior <- adjustcolor("grey50", alpha.f = 0.45)
col_post  <- adjustcolor("steelblue", alpha.f = 0.40)

par(mfrow = c(1, 1), mar = c(5.1, 5.1, 4.1, 2.1), mgp = c(3.5, 1, 0))

plot(
  prior_hist,
  freq   = FALSE,
  col    = col_prior,
  border = "white",
  xlim   = c(0, 800),
  ylim   = c(0, ymax),
  main   = "Prior and posterior predictive edge counts",
  xlab   = "Number of edges",
  ylab   = "Proportion of simulated graphs"
)

plot(post_hist, freq = FALSE, col = col_post, border = "white", add = TRUE)

abline(v = 197, col = "black", lwd = 2, lty = 4)

legend(
  "topright",
  legend = c("Prior predictive", "Posterior predictive", "Observed = 197"),
  col    = c("grey50", "steelblue", "black"),
  pch    = c(15, 15, NA),
  pt.cex = 2,
  lty    = c(NA, NA, 4),
  lwd    = c(NA, NA, 2),
  bty    = "n"
)



# =========================================================================
# Frequentist vs Bayesian goodness of fit for M6
#   net ~ edges + isolates + nodematch("citizenship") + gwesp(0.25, fixed = TRUE)
#
# Five panels: degree, edgewise shared partners, dyadwise shared partners,
# non-edgewise shared partners, geodesic distance. Each bin gets two dodged
# boxes, MLE and Bayesian, with the median as the box line and the mean as
# an open diamond.
#
# Assumes these already exist from chapter6_ise_network.R:
#   net, g, f_M6, fit, bergm_draws, n_nodes, n_dyads, mle_values, TERMS,
#   PRIOR_E, PRIOR_ISO, PRIOR_C, PRIOR_GW, dlaplace(), block_seed()
# =========================================================================

library(sna)

# sna is attached last here deliberately, so sna::degree and sna::geodist
# are the ones used below. The igraph bindings made in the main script are
# not needed in this file.

# =========================================================================
# 0. Design of the comparison
# =========================================================================
#
# Both procedures simulate networks and compare the observed statistics to
# the simulated distribution. The only thing that differs is where theta
# comes from:
#
#   Frequentist  theta = theta_hat, held fixed and reused S times
#                -> spread reflects graph variability at one theta
#
#   Bayesian     theta ~ p(theta | y), a fresh draw for every network
#                -> spread reflects graph variability AND uncertainty about
#                   theta, filtered through the prior
#
# Everything downstream is identical: same simulation engine, same burn-in,
# same number of draws, same statistics. Do not substitute gof() for the
# frequentist arm, because its default simulation controls differ from those
# used in the posterior predictive block, and then part of any difference in
# spread is an artefact of the settings rather than of the paradigm.
#
# On the five panels, and which are in-model:
#   degree    only the zero bin is fitted, via the isolates term; the rest
#             of the distribution is a genuine out-of-model check, and the
#             degree-40 hub is where the model will struggle
#   esp       IN MODEL via gwesp, so read it as a calibration check
#   dsp       closure opportunity across all pairs, tied or not. NOT fitted
#   nsp       closure among pairs that are not tied. NOT fitted
#   distance  global connectivity. NOT fitted
#
# dsp and nsp are the informative pair here, and the comparison is the
# mirror image of Chapter 4. There the model had a degree term but no
# closure term; here it has closure but almost no degree structure. If esp
# is reproduced while dsp is overstated, the model is generating too much
# closure opportunity and then closing the right fraction of it, which is
# what a single gwesp term with no gwdsp companion would do.
#
# Note also that 36 of the 197 observed ties lie in a single nine-actor
# clique, every one of them with exactly seven shared partners. That block
# is what the esp panel will be dominated by at the upper bins.

# =========================================================================
# 1. Are the prior and the likelihood describing the same model?
# =========================================================================
#
#   theta_iso   ~ Uniform(-5, 5)      support [-5, 5]     HARD
#   theta_GWESP ~ Gamma(9, 2.5)       support (0, Inf)    HARD
#
# If theta_hat falls outside either, the posterior cannot reach the MLE and
# the two arms are fitting different models.

plaplace <- function(q, mu = 0, b = 1) {
  ifelse(q < mu,
         0.5 * exp((q - mu) / b),
         1 - 0.5 * exp(-(q - mu) / b))
}

prior_support <- data.frame(
  parameter = names(mle_values),
  mle       = as.numeric(mle_values),
  lower     = c(-Inf, PRIOR_ISO$min, -Inf, 0),
  upper     = c( Inf, PRIOR_ISO$max,  Inf, Inf),
  row.names = NULL
)

prior_support$inside_support <- with(prior_support, mle > lower & mle < upper)

prior_support$prior_quantile <- c(
  pnorm(mle_values[1], PRIOR_E$mean, PRIOR_E$sd),
  punif(mle_values[2], PRIOR_ISO$min, PRIOR_ISO$max),
  plaplace(mle_values[3], PRIOR_C$mu, PRIOR_C$b),
  pgamma(mle_values[4], PRIOR_GW$shape, PRIOR_GW$rate)
)

prior_support

# A prior quantile near 0 or 1 means the prior is pulling against the
# likelihood, so the Bayesian boxes will be shifted as well as rescaled.
# theta_E is the one to watch: the prior is centred on the midpoint of the
# two projections' estimates, not on the standard-network MLE.

# =========================================================================
# 2. Observed statistics
# =========================================================================

SP_RANGE  <- 0:10    # shared-partner bins, wide enough to reach the clique
DEG_RANGE <- 0:15    # three actors lie above this; see the tail note below
DIST_KEEP <- 1:6     # diameter is 6 in all three projections

obs_esp <- summary(net ~ esp(SP_RANGE))
obs_dsp <- summary(net ~ dsp(SP_RANGE))
obs_nsp <- summary(net ~ nsp(SP_RANGE))

obs_scalar <- summary(
  net ~ edges + isolates + nodematch("citizenship") +
    gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE) + triangle
)

obs_degree <- tabulate(
  sna::degree(net, gmode = "graph") + 1,
  nbins = n_nodes
)

obs_geodist <- {
  gd <- sna::geodist(net, inf.replace = Inf)$gdist
  gd <- gd[upper.tri(gd)]
  c(tabulate(gd[is.finite(gd)], nbins = n_nodes - 1), sum(!is.finite(gd)))
}

# The degree tail: 17, 24 and 40 in the observed network. No ERGM with only
# an isolates term will reproduce a degree-40 hub, so the upper bins are
# excluded from the panel and reported separately rather than stretching
# the axis to 40 for one point.
sort(sna::degree(net, gmode = "graph")[sna::degree(net, gmode = "graph") > 15])

# =========================================================================
# 3. Shared simulation engine
# =========================================================================
#
# Slower than Chapter 4 by a wide margin: 3570 dyads rather than 120, and
# the networks are returned rather than just their statistics. Budget about
# half an hour per arm.

S_arm      <- 5000      # networks per arm
SIM_BURNIN <- 20000    # identical for both arms

simulate_gof_arm <- function(theta_mat, label, burnin = SIM_BURNIN) {
  
  S <- nrow(theta_mat)
  
  out <- list(
    theta   = theta_mat,
    degree  = matrix(NA_real_, S, n_nodes),
    esp     = matrix(NA_real_, S, length(obs_esp)),
    dsp     = matrix(NA_real_, S, length(obs_dsp)),
    nsp     = matrix(NA_real_, S, length(obs_nsp)),
    geodist = matrix(NA_real_, S, n_nodes),
    scalar  = matrix(NA_real_, S, length(obs_scalar),
                     dimnames = list(NULL, names(obs_scalar)))
  )
  
  for (s in seq_len(S)) {
    
    sim_net <- tryCatch(
      simulate(
        f_M6("net"),
        coef = unname(theta_mat[s, ]),
        nsim = 1,
        output = "network",
        control = control.simulate.formula(
          MCMC.burnin = burnin,
          MCMC.interval = 1
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(sim_net)) next
    
    d <- sna::degree(sim_net, gmode = "graph")
    
    out$degree[s, ] <- tabulate(d + 1, nbins = n_nodes)
    out$esp[s, ]    <- summary(sim_net ~ esp(SP_RANGE))
    out$dsp[s, ]    <- summary(sim_net ~ dsp(SP_RANGE))
    out$nsp[s, ]    <- summary(sim_net ~ nsp(SP_RANGE))
    
    out$scalar[s, ] <- summary(
      sim_net ~ edges + isolates + nodematch("citizenship") +
        gwesp(0.25, fixed = TRUE) + gwdsp(0.25, fixed = TRUE) + triangle
    )
    
    gd <- sna::geodist(sim_net, inf.replace = Inf)$gdist
    gd <- gd[upper.tri(gd)]
    
    out$geodist[s, ] <- c(
      tabulate(gd[is.finite(gd)], nbins = n_nodes - 1),
      sum(!is.finite(gd))
    )
    
    if (s %% 50 == 0) message(label, ": ", s, " of ", S)
  }
  
  out
}

# =========================================================================
# 4. The two arms
# =========================================================================

set.seed(block_seed("gof_ise_freq"))

theta_freq <- matrix(
  rep(as.numeric(mle_values), each = S_arm),
  nrow = S_arm
)

arm_freq <- simulate_gof_arm(theta_freq, "frequentist (plug-in MLE)")

set.seed(block_seed("gof_ise_bayes"))

bayes_index <- sample(
  seq_len(nrow(bergm_draws)),
  size = min(S_arm, nrow(bergm_draws)),
  replace = FALSE
)

arm_bayes <- simulate_gof_arm(
  bergm_draws[bayes_index, , drop = FALSE],
  "Bayesian posterior"
)

arms <- list(
  "Frequentist (MLE)" = arm_freq,
  "Bayesian (BERGM)"  = arm_bayes
)

arm_cols <- c("darkred", "steelblue")

# =========================================================================
# 5. Box-and-whisker figures
# =========================================================================
#
# Boxes are drawn by hand rather than via bxp(), so the whiskers can be set
# to the 2.5% and 97.5% quantiles instead of Tukey's 1.5 x IQR rule. Two
# reasons: it keeps continuity with the 95% intervals used elsewhere, and
# with 500 draws Tukey whiskers plus outlier points would fill the upper
# bins with dots that say nothing beyond what the degeneracy table reports.
#
# Box  = interquartile range
# Line = median
# Diamond = mean, plotted so the skew is visible rather than argued about

box_stats <- function(v, whisker_probs = c(0.025, 0.975)) {
  v <- v[!is.na(v)]
  c(
    quantile(v, whisker_probs[1], names = FALSE),
    quantile(v, 0.25, names = FALSE),
    median(v),
    quantile(v, 0.75, names = FALSE),
    quantile(v, whisker_probs[2], names = FALSE)
  )
}

draw_box <- function(x, st, width, fill, border) {
  # st = c(lower whisker, Q1, median, Q3, upper whisker)
  segments(x, st[1], x, st[2], col = border, lwd = 1.4)
  segments(x, st[4], x, st[5], col = border, lwd = 1.4)
  segments(x - width / 4, st[1], x + width / 4, st[1], col = border, lwd = 1.4)
  segments(x - width / 4, st[5], x + width / 4, st[5], col = border, lwd = 1.4)
  rect(x - width / 2, st[2], x + width / 2, st[4],
       col = fill, border = border, lwd = 1.4)
  segments(x - width / 2, st[3], x + width / 2, st[3], col = border, lwd = 2.6)
}

plot_gof_boxes <- function(
    component,
    obs_vec,
    keep,
    x_labels,
    main,
    xlab,
    ylab,
    box_width = 0.30,
    offsets = c(-0.18, 0.18),
    connect = seq_along(keep),
    show_legend = FALSE,
    legend_pos = "topleft",
    log_y = FALSE
) {
  
  k   <- length(keep)
  pos <- seq_len(k)
  
  stats_list <- lapply(arms, function(a) {
    apply(a[[component]][, keep, drop = FALSE], 2, box_stats)
  })
  
  means_list <- lapply(arms, function(a) {
    colMeans(a[[component]][, keep, drop = FALSE], na.rm = TRUE)
  })
  
  ymax <- 1.05 * max(
    unlist(stats_list), unlist(means_list), obs_vec, na.rm = TRUE
  )
  
  plot(
    NA,
    xlim = c(0.5, k + 0.5),
    ylim = c(0, ymax),
    xaxt = "n",
    xlab = xlab, ylab = ylab, main = main, font.main = 2
  )
  
  axis(side = 1, at = pos, labels = x_labels)
  
  for (a in seq_along(arms)) {
    
    st <- stats_list[[a]]
    
    for (j in seq_len(k)) {
      draw_box(
        x = pos[j] + offsets[a],
        st = st[, j],
        width = box_width,
        fill = adjustcolor(arm_cols[a], alpha.f = 0.30),
        border = arm_cols[a]
      )
    }
    
    points(
      pos + offsets[a], means_list[[a]],
      pch = 23, bg = "white", col = arm_cols[a], cex = 0.85, lwd = 1.6
    )
  }
  
  # Observed
  if (length(connect) > 1) {
    lines(pos[connect], obs_vec[connect], col = "black", lwd = 1.6, lty = 2)
  }
  points(pos, obs_vec, pch = 19, cex = 1.0)
  span <- max(abs(offsets)) + box_width / 2
  segments(pos - span, obs_vec, pos + span, obs_vec,
           col = "black", lwd = 2.2, lty = 3)
  
  if (show_legend) {
    legend(
      legend_pos,
      legend = c("Observed", names(arms), "Median (box line)", "Mean"),
      col   = c("black", arm_cols, "grey20", "grey20"),
      pt.bg = c("black", adjustcolor(arm_cols, alpha.f = 0.30), NA, "white"),
      pch   = c(19, 22, 22, NA, 23),
      lty   = c(2, NA, NA, 1, NA),
      lwd   = c(1.6, NA, NA, 2.6, 1.6),
      pt.cex = c(1.0, 1.6, 1.6, NA, 0.85),
      bty = "n", cex = 0.7
    )
  }
}

par(
  mfrow = c(1, 1),
  oma = c(0, 0, 0, 0),
  mar = c(4.2, 4.4, 2.2, 1.2),
  mgp = c(2.6, 0.8, 0)
)

plot_gof_boxes(
  "degree", obs_degree[DEG_RANGE + 1], keep = DEG_RANGE + 1,
  x_labels = DEG_RANGE,
  main = "Degree Distribution GOF",
  xlab = "Degree", ylab = "Number of actors",
  show_legend = TRUE
)

plot_gof_boxes(
  "esp", obs_esp, keep = seq_along(obs_esp), x_labels = SP_RANGE,
  main = "Edgewise Shared Partners Distribution GOF",
  xlab = "Shared partners", ylab = "Number of edges",
  show_legend = TRUE
)

plot_gof_boxes(
  "dsp", obs_dsp, keep = seq_along(obs_dsp), x_labels = SP_RANGE,
  main = "Dyadwise Shared Partners Distribution GOF",
  xlab = "Shared partners", ylab = "Number of dyads",
  show_legend = TRUE
)

plot_gof_boxes(
  "nsp", obs_nsp, keep = seq_along(obs_nsp), x_labels = SP_RANGE,
  main = "Non-edgewise Shared Partners GOF",
  xlab = "Shared partners", ylab = "Number of non-tied dyads",
  show_legend = TRUE
)

plot_gof_boxes(
  "geodist", obs_geodist[c(DIST_KEEP, n_nodes)],
  keep = c(DIST_KEEP, n_nodes),
  x_labels = c(DIST_KEEP, "NR"),
  connect = seq_along(DIST_KEEP),
  main = "Geodesic Distance Distribution GOF",
  xlab = "Distance", ylab = "Number of dyads",
  show_legend = TRUE
)

# The dsp and nsp panels are dominated by the zero bin: with 3570 dyads and
# 197 edges most pairs share no partners at all. If that squashes the rest,
# redraw with keep = 2:11 and say in the caption that the zero bin was
# dropped for legibility.

plot_gof_boxes(
  "dsp", obs_dsp[-1], keep = 2:length(obs_dsp), x_labels = SP_RANGE[-1],
  main = "Dyadwise Shared Partners GOF (zero bin omitted)",
  xlab = "Shared partners", ylab = "Number of dyads",
  show_legend = TRUE
)

# =========================================================================
# 6. Where do the mean and the median disagree?
# =========================================================================
#
# A large gap between the diamond and the box line means that bin's
# predictive distribution is strongly right-skewed. Worth reporting: it is
# the reason the median is the honest central summary here, and it also
# flags which bins the degeneracy table in section 10 is speaking about.

skew_check <- function(component, keep, bin_labels) {
  do.call(rbind, lapply(names(arms), function(nm) {
    m <- arms[[nm]][[component]][, keep, drop = FALSE]
    data.frame(
      arm    = nm,
      bin    = bin_labels,
      mean   = round(colMeans(m, na.rm = TRUE), 2),
      median = round(apply(m, 2, median, na.rm = TRUE), 2),
      gap    = round(colMeans(m, na.rm = TRUE) -
                       apply(m, 2, median, na.rm = TRUE), 2),
      row.names = NULL
    )
  }))
}

skew_check("degree",  DEG_RANGE + 1, paste("degree", DEG_RANGE))
skew_check("esp",     seq_along(obs_esp), paste("esp", SP_RANGE))
skew_check("dsp",     seq_along(obs_dsp), paste("dsp", SP_RANGE))
skew_check("geodist", c(DIST_KEEP, n_nodes),
           c(paste("dist", DIST_KEEP), "unreachable"))

# =========================================================================
# 7. Spread of the boxes
# =========================================================================

interval_widths <- function(mat, keep, probs = c(0.025, 0.975)) {
  apply(mat[, keep, drop = FALSE], 2, function(v) {
    diff(quantile(v, probs = probs, na.rm = TRUE))
  })
}

width_table <- function(component, keep, bin_labels, probs = c(0.025, 0.975)) {
  
  w <- sapply(arms, function(a) interval_widths(a[[component]], keep, probs))
  
  out <- data.frame(bin = bin_labels, round(w, 2), row.names = NULL)
  names(out) <- c("bin", names(arms))
  
  out$ratio <- round(out[["Bayesian (BERGM)"]] / out[["Frequentist (MLE)"]], 2)
  out
}

width_table("degree",  DEG_RANGE + 1, paste("degree", DEG_RANGE))
width_table("esp",     seq_along(obs_esp), paste("esp", SP_RANGE))
width_table("dsp",     seq_along(obs_dsp), paste("dsp", SP_RANGE))
width_table("geodist", c(DIST_KEEP, n_nodes),
            c(paste("dist", DIST_KEEP), "unreachable"))

# Box height (IQR) only, which is what the eye actually compares
width_table("esp", seq_along(obs_esp), paste("esp", SP_RANGE),
            probs = c(0.25, 0.75))

# ratio > 1  the Bayesian box is taller: the plug-in understates predictive
#            uncertainty by ignoring uncertainty in theta
# ratio ~ 1  parameter uncertainty is small relative to graph variability
#            at this network size
# ratio < 1  the prior is tightening the box. With Uniform(-5, 5) on
#            theta_iso this is unlikely; with Gamma(9, 2.5) on theta_GWESP
#            and N(-6.70, 0.5^2) on theta_E it is a real possibility, and it
#            needs defending rather than reporting. Compare against the
#            prior_quantile column from section 1.

# =========================================================================
# 8. p-values, and why neither is calibrated
# =========================================================================

pval_table <- do.call(rbind, lapply(names(arms), function(nm) {
  
  sim <- arms[[nm]]$scalar
  
  data.frame(
    arm        = nm,
    statistic  = colnames(sim),
    observed   = as.numeric(obs_scalar),
    sim_mean   = colMeans(sim, na.rm = TRUE),
    sim_median = apply(sim, 2, median, na.rm = TRUE),
    sim_lower  = apply(sim, 2, quantile, probs = 0.025, na.rm = TRUE),
    sim_upper  = apply(sim, 2, quantile, probs = 0.975, na.rm = TRUE),
    p_upper    = sapply(seq_len(ncol(sim)), function(j) {
      mean(sim[, j] >= obs_scalar[j], na.rm = TRUE)
    }),
    row.names = NULL
  )
}))

pval_table$p_two_sided <- 2 * pmin(pval_table$p_upper, 1 - pval_table$p_upper)
pval_table$inside_95   <- with(pval_table,
                               observed >= sim_lower & observed <= sim_upper)

pval_table[order(pval_table$statistic, pval_table$arm), ]

# The four in-model statistics are matched by construction in the
# frequentist arm, so their p-values are not tests of anything. Read
# gwdsp and triangle instead: both are out of model, and both speak to
# whether a single closure term is enough.
#
# Miscalibrated in opposite directions:
#
#   Frequentist bootstrap p: treats theta_hat as known, so the reference
#     distribution is too narrow and p is pushed towards 0 or 1. Declares
#     misfit too readily.
#
#   Posterior predictive p: the data are used to fit and then to check, so
#     it is conservative and pulled towards 0.5. Declares misfit too
#     reluctantly.

# =========================================================================
# 9. Approximate variance decomposition
# =========================================================================
#
#   Var(T(y)) = E_theta[ Var(T | theta) ] + Var_theta( E[T | theta] )
#                \___ within ___/           \___ between ___/
#
# The frequentist arm estimates Var(T | theta_hat), standing in for the
# within term; the Bayesian arm estimates the total. The gap approximates
# the between term.
#
# Approximate only: Var(T | theta) is not constant in theta, so the
# frequentist arm evaluates the within term at a point rather than
# averaging over the posterior. between_approx can come out negative, which
# is the decomposition breaking down rather than a coding error.

var_decomposition <- data.frame(
  statistic     = colnames(arm_freq$scalar),
  within_at_mle = apply(arm_freq$scalar,  2, var, na.rm = TRUE),
  total_bayes   = apply(arm_bayes$scalar, 2, var, na.rm = TRUE),
  row.names = NULL
)

var_decomposition$between_approx <-
  var_decomposition$total_bayes - var_decomposition$within_at_mle

var_decomposition$prop_from_theta <-
  var_decomposition$between_approx / var_decomposition$total_bayes

var_decomposition[, -1] <- round(var_decomposition[, -1], 3)
var_decomposition

# With 3570 dyads rather than 120, expect prop_from_theta to be smaller
# than in Chapter 4: graph variability at fixed theta grows with the number
# of dyads, while posterior uncertainty about theta shrinks. That contrast
# between the two case studies is worth a sentence in Chapter 7.

# =========================================================================
# 10. Degeneracy rates
# =========================================================================
#
# The whiskers stop at the 2.5% and 97.5% quantiles, so the extremes are
# deliberately not drawn. This table is where they are accounted for.

degeneracy_table <- do.call(rbind, lapply(names(arms), function(nm) {
  
  e <- arms[[nm]]$scalar[, "edges"]
  
  data.frame(
    arm             = nm,
    failed_sims     = mean(is.na(e)),
    exactly_extreme = mean(e %in% c(0, n_dyads), na.rm = TRUE),
    near_extreme    = mean(e <= 2 | e >= n_dyads - 2, na.rm = TRUE),
    median_edges    = median(e, na.rm = TRUE),
    mean_edges      = mean(e, na.rm = TRUE),
    row.names = NULL
  )
}))

degeneracy_table

# =========================================================================
# 11. The hub and the clique
# =========================================================================
#
# Two features of this network that no term in M6 can reproduce, and which
# the panels above will show as specific failures rather than general
# misfit. Report them as such.

# (a) Can either arm generate a degree-40 actor?
sapply(arms, function(a) {
  max_deg <- apply(a$degree, 1, function(r) max(which(r > 0)) - 1)
  c(median_max_degree = median(max_deg, na.rm = TRUE),
    q975_max_degree   = quantile(max_deg, 0.975, na.rm = TRUE, names = FALSE),
    observed          = max(sna::degree(net, gmode = "graph")),
    prop_reaching_40  = mean(max_deg >= 40, na.rm = TRUE))
})

# (b) The esp7 bin: 37 observed ties, 36 of them the leadership clique.
sapply(arms, function(a) {
  v <- a$esp[, which(SP_RANGE == 7)]
  c(median = median(v, na.rm = TRUE),
    q025   = quantile(v, .025, na.rm = TRUE, names = FALSE),
    q975   = quantile(v, .975, na.rm = TRUE, names = FALSE),
    observed = obs_esp[which(SP_RANGE == 7)],
    p_upper  = mean(v >= obs_esp[which(SP_RANGE == 7)], na.rm = TRUE))
})

# (c) Isolates are in model, so both arms should match 8 on average. If the
# frequentist arm does and the Bayesian arm does not, the Uniform(-5, 5)
# prior is doing something and it needs explaining.
sapply(arms, function(a) {
  v <- a$scalar[, "isolates"]
  c(mean = mean(v, na.rm = TRUE), median = median(v, na.rm = TRUE),
    q025 = quantile(v, .025, na.rm = TRUE, names = FALSE),
    q975 = quantile(v, .975, na.rm = TRUE, names = FALSE),
    observed = obs_scalar["isolates"])
})

save(arms, obs_esp, obs_dsp, obs_nsp, obs_degree, obs_geodist, obs_scalar,
     pval_table, var_decomposition, degeneracy_table, prior_support,
     file = file.path(OUT, "chapter6_gof.RData"))




# =========================================================================
# Typical networks under M6: observed, two at the MLE, two at the MAP
#
# The GOF box plots say whether the fitted statistics are matched. This
# figure says whether the resulting graphs look like the observed one, which
# is a different and more visible question.
#
# Two draws from each parameter vector rather than four from one, so the
# comparison is between paradigms and not only between the model and the
# data. The MLE and the MAP differ most in theta_E and theta_GWESP, and
# those move in compensating directions, so the two pairs may look very
# similar: that itself is worth reporting.
#
# Assumes: net, f_M6, fit, mle_values, post_tab, block_seed(), n_nodes
# =========================================================================

library(intergraph)

set.seed(block_seed("typical_grap5h5s10"))

# Parameter vectors. post_tab is the posterior summary table from the BERGM
# block; the mode column is the MAP.
theta_mle <- as.numeric(mle_values)
theta_map <- as.numeric(post_tab[, "mode"])

round(rbind(MLE = theta_mle, MAP = theta_map), 3)

SIM_BURNIN <- 20000

draw_networks <- function(theta, n = 2) {
  simulate(
    f_M6("net"),
    coef    = unname(theta),
    nsim    = n,
    output  = "network",
    control = control.simulate.formula(
      MCMC.burnin   = SIM_BURNIN,
      MCMC.interval = SIM_BURNIN
    )
  )
}

sims_mle <- draw_networks(theta_mle, 2)
sims_map <- draw_networks(theta_map, 2)

# simulate() returns a bare network when nsim = 1 and a list otherwise.
as_list <- function(x) if (inherits(x, "network")) list(x) else as.list(x)

panels <- c(
  list(observed = intergraph::asIgraph(net)),
  lapply(as_list(sims_mle), intergraph::asIgraph),
  lapply(as_list(sims_map), intergraph::asIgraph)
)

labels <- c("Observed", "MLE draw 1", "MLE draw 2",
            "MAP draw 1", "MAP draw 2")

border <- c("black", "darkred", "darkred", "steelblue", "steelblue")


# =========================================================================
# Node sizing
# =========================================================================
#
# Square root of degree, because area scales with the square of the radius
# and plain degree would make the degree-40 hub swamp everything else.
#
# The scale is fixed across panels using the observed maximum, so a node of
# a given size means the same degree in every panel. Without this the
# simulated graphs are silently rescaled and the missing hub disappears,
# which is the one thing the figure exists to show.

MAXDEG <- max(igraph::degree(panels[[1]]))

node_size <- function(gi) 2.5 + 7 * sqrt(igraph::degree(gi) / MAXDEG)

panel_stats <- function(gi) {
  sprintf("%d edges, %d isolates, C = %.2f, max deg %d",
          igraph::ecount(gi),
          sum(igraph::degree(gi) == 0),
          igraph::transitivity(gi, type = "global"),
          max(igraph::degree(gi)))
}


# =========================================================================
# Figure
# =========================================================================

layout(matrix(c(1, 1, 2, 3,
                1, 1, 4, 5), nrow = 2, byrow = TRUE))

par(mar = c(0, 0, 3.6, 0))

for (i in seq_along(panels)) {
  
  gi <- panels[[i]]
  
  set.seed(block_seed(paste0("layout_", i)))
  set.seed(95)
  
  plot(
    gi,
    layout             = igraph::layout_with_fr(gi, niter = 3000),
    vertex.size        = node_size(gi)+3,
    vertex.label       = NA,
    vertex.color       = if (i == 1) "grey40" else
      adjustcolor(border[i], alpha.f = 0.75),
    vertex.frame.color = "white",
    edge.color         = "grey78",
    edge.width         = 0.8
  )
  
  title(main = labels[i], line = 1.2, cex.main = if (i == 1) 2.95 else 2.95,
        col.main = border[i], font.main = 2)
  mtext(panel_stats(gi), side = 3, line = 0.05,
        cex = if (i == 1) 0.9 else 0.9, col = "grey30")
  
  box(col = "grey85")
}

layout(1)


# =========================================================================
# The numbers behind the figure
# =========================================================================
#
# Node sizes are hard to read precisely, so report the maximum degree
# alongside. This is the feature neither parameter vector can reproduce: M6
# has no term for degree concentration, so the closure it fits is spread
# evenly rather than concentrated around one broker.

set.seed(block_seed("typical_graphs_check10"))

more_mle <- draw_networks(theta_mle, 100)
more_map <- draw_networks(theta_map, 100)

max_deg <- function(sims) sapply(as_list(sims), function(x)
  max(sna::degree(x, gmode = "graph")))

rbind(
  MLE = c(median = median(max_deg(more_mle)),
          q975   = quantile(max_deg(more_mle), .975, names = FALSE),
          prop_reaching_40 = mean(max_deg(more_mle) >= 40)),
  MAP = c(median = median(max_deg(more_map)),
          q975   = quantile(max_deg(more_map), .975, names = FALSE),
          prop_reaching_40 = mean(max_deg(more_map) >= 40))
)

max(sna::degree(net, gmode = "graph"))   # 40


# joint MAP estimate.

library(ks)

H <- ks::Hpi(bergm_draws)
d <- ks::kde(bergm_draws, H = H, eval.points = bergm_draws)
joint_map <- bergm_draws[which.max(d$estimate), ]

round(rbind(joint_MAP = joint_map,
            marginal  = post_tab[, "mode"],
            MLE       = mle_values), 3)
# =========================================================================
# Session record
# =========================================================================
#
# A seed pins the randomness, not the package versions. Record both.

# writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

save(g, g_l, g_e, net, net_l, net_e,
     comp_tab, er_tab, deg_tab, esp_tab, alpha_grid,
     fit, fit_l, coef_tab, mle_values,
     bergm_fit, bergm_draws, post_tab,
     valid_prior_edges, valid_post_edges,
     file = file.path(OUT, "chapter6.RData"))

