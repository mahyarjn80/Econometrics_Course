# --- Load Required Libraries ---
install.packages("devtools")
devtools::install_github("bmelly/discreteQ")

library(discreteQ)
library(xtable)
library(parallel)

# --- Load Data and Initialize Environment ---
rm(list = ls())  # Clear workspace
options(warn = -1)  # Suppress warnings

# Load dataset (update path if needed)
load("../../Data/cps2015.rdata")
attach(data)

# -------------------------------
# Summary Statistics Preparation
# -------------------------------
# Variables used in descriptive stats
vars <- c("lnw", "female", "married", "widowed", "separated", "divorced", 
          "nevermarried", "lhs", "hsg", "sc", "cg", "ad", 
          "ne", "mw", "so", "we", "exp1")

options(digits = 2)

# Weighted means for entire sample, women only, and men only
dstats <- cbind(
  sapply(data[, vars], weighted.mean, weight),
  sapply(data[female == 1, vars], weighted.mean, weight[female == 1]),
  sapply(data[female == 0, vars], weighted.mean, weight[female == 0])
)

# Output in LaTeX-style table (can be exported or rendered in Rmd)
xtable(dstats)

# ----------------------------------
# Distribution Regression Parameters
# ----------------------------------

nind  <- 39  # Number of distribution regressions (one per threshold)
taus  <- c(3:(nind - 2)) / (nind + 1)  # Intermediate quantiles for analysis
alpha <- 0.10  # Confidence level for 90% CI

# Specification of covariates for distribution regression
formula <- lnw ~ widowed + divorced + separated + nevermarried +
  (lhs + hsg + cg + ad) * (exp1 + exp2 + exp3 + exp4) +
  mw + so + we

X <- model.matrix(formula, data = data)  # Design matrix

# Grid of threshold values (quantiles) to run the DR on
ys <- quantile(lnw, probs = c(0.02, c(1:nind) / (nind + 1), 0.98))

R <- 200  # Bootstrap replications
set.seed(1)

# Set up a parallel computing cluster for faster estimation
my_cl <- makePSOCKcluster(24, setup_timeout = 0.5)

# -------------------------------------
# Distribution Regression Estimation
# -------------------------------------
# Running discreteQ() performs a series of binary regressions for wage <= y_s
# Then estimates counterfactual distributions and quantile decompositions

dr.fit <- discreteQ(
  y = lnw,
  d = 1 - female,  # Define "treatment" as being male
  x = X,
  w = weight,
  decomposition = TRUE,  # Enables Oaxaca-like decomposition
  q.range = range(taus),
  method = "logit",  # Logit link for distribution regression
  bsrep = R,
  alpha = alpha,
  ys = ys,
  cl = my_cl
)

# ----------------------------------
# Visualizing Wage Distributions
# ----------------------------------

par(mfrow = c(1, 2))

# Histogram of log hourly wages for women
hist(
  lnw[female == 1],
  breaks = unique(lnw[female == 1]),
  col = "tomato",
  main = "Women: Log Hourly Wages",
  xlab = "Log Wage",
  ylab = "Frequency"
)

# Histogram for men
hist(
  lnw[female == 0],
  breaks = unique(lnw[female == 0]),
  col = "steelblue",
  main = "Men: Log Hourly Wages",
  xlab = "Log Wage",
  ylab = "Frequency"
)

par(mfrow = c(1, 1))

# -------------------------------------------
# Plot Observed & Counterfactual CDFs (F(y))
# -------------------------------------------

plot(
  range(ys), c(0, 1),
  type = "n",
  xlab = "Log Hourly Wage",
  ylab = "Cumulative Probability",
  main = "Observed vs Counterfactual Distributions (90% CI)"
)

# Observed CDF for Women with CI
polygon(c(ys, rev(ys)), c(dr.fit$ub.F0(ys), rev(dr.fit$lb.F0(ys))),
        col = "lightblue", border = NA)
lines(ys, dr.fit$F0(ys), col = "darkgreen", lwd = 2)

# Observed CDF for Men with CI
polygon(c(ys, rev(ys)), c(dr.fit$ub.F1(ys), rev(dr.fit$lb.F1(ys))),
        col = "lightgreen", border = NA)
lines(ys, dr.fit$F1(ys), col = "blue", lwd = 2)

# Counterfactual CDF (what men would earn with women's covariates)
polygon(c(ys, rev(ys)), c(dr.fit$ub.Fc(ys), rev(dr.fit$lb.Fc(ys))),
        col = "grey90", border = NA)
lines(ys, dr.fit$Fc(ys), col = "grey30", lwd = 2)

legend("bottomright", bty = "n", lwd = 2,
       legend = c("Women (Obs)", "Men (Obs)", "Men (CF w/ Women's X)"),
       col = c("darkgreen", "blue", "grey30"))

# -------------------------------------------
# Plot Observed & Counterfactual Quantiles
# -------------------------------------------

plot(
  range(taus), range(ys),
  type = "n",
  xlab = "Quantile Index",
  ylab = "Log Hourly Wage",
  main = "Observed vs Counterfactual Quantiles (90% CI)"
)

# Quantile curves
polygon(c(taus, rev(taus)), c(dr.fit$ub.Q0(taus), rev(dr.fit$lb.Q0(taus))),
        col = "lightblue", border = NA)
lines(taus, dr.fit$Q0(taus), col = "darkgreen", lwd = 2)

polygon(c(taus, rev(taus)), c(dr.fit$ub.Q1(taus), rev(dr.fit$lb.Q1(taus))),
        col = "lightgreen", border = NA)
lines(taus, dr.fit$Q1(taus), col = "blue", lwd = 2)

polygon(c(taus, rev(taus)), c(dr.fit$ub.Qc(taus), rev(dr.fit$lb.Qc(taus))),
        col = "grey90", border = NA)
lines(taus, dr.fit$Qc(taus), col = "grey30", lwd = 2)

legend("topleft", bty = "n", lwd = 2,
       legend = c("Women (Obs)", "Men (Obs)", "Men (CF w/ Women's X)"),
       col = c("darkgreen", "blue", "grey30"))

# -----------------------------------
# Plot Decomposition of Wage Gap
# -----------------------------------

par(mfrow = c(3, 1))

plot(
  dr.fit, which = "observed", support = "continuous",
  main = "Total Gender Wage Gap (90% CI)", ylim = c(-0.15, 0.5)
)
abline(h = 0, col = "grey70")

plot(
  dr.fit, which = "unexplained", support = "continuous",
  main = "Unexplained (Discrimination?) Component", ylim = c(-0.15, 0.5)
)
abline(h = 0, col = "grey70")

plot(
  dr.fit, which = "composition", support = "continuous",
  main = "Composition (Covariate) Effect", ylim = c(-0.15, 0.5)
)
abline(h = 0, col = "grey70")



