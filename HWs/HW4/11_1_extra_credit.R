library(foreign)
library(xtable)
library(plm)
library(gmm)
library(readstata13)

# Load the democracy dataset
data <- read.dta13("../../Data/democracy-balanced-l4.dta")
data <- pdata.frame(data, index = c("id", "year"))

# Define the fixed effects model formula
form.fe <- lgdp ~ dem + lag(lgdp, 1:4) - 1

# Estimate the model on the original data to get true parameters
fe.fit <- plm(form.fe, data, model = "within", effect = "twoways", index = c("id", "year"))
coefs.fe <- coef(fe.fit)  # True coefficients
a_k <- fixef(fe.fit, effect = "individual")  # Country fixed effects
b_t <- fixef(fe.fit, effect = "time")  # Year fixed effects
res <- residuals(fe.fit)
sigma2_hat <- mean(res^2)  # Estimated error variance

# Extract dimensions and data matrices
N <- length(levels(data$id))  # 147 countries
T <- length(levels(data$year))  # 23 years (1987-2009)
years <- levels(data$year)
ids <- levels(data$id)

lgdp_data <- matrix(data$lgdp, nrow = N, ncol = T, byrow = FALSE)
dem_data <- matrix(data$dem, nrow = N, ncol = T, byrow = FALSE)

# Align fixed effects with ids and years
a_k_vec <- a_k[match(ids, names(a_k))]
b_t_vec <- b_t[match(years, names(b_t))]

# Function to simulate synthetic data
simulate_data <- function() {
  epsilon <- matrix(rnorm(N * T, 0, sqrt(sigma2_hat)), nrow = N, ncol = T)
  lgdp_sim <- lgdp_data  # Initialize with original data for t=1:4
  for (t in 5:T) {
    lag1 <- lgdp_sim[, t - 1]
    lag2 <- lgdp_sim[, t - 2]
    lag3 <- lgdp_sim[, t - 3]
    lag4 <- lgdp_sim[, t - 4]
    lgdp_sim[, t] <- coefs.fe[1] * dem_data[, t] +
      coefs.fe[2] * lag1 +
      coefs.fe[3] * lag2 +
      coefs.fe[4] * lag3 +
      coefs.fe[5] * lag4 +
      a_k_vec + b_t_vec[t] +
      epsilon[, t]
  }
  data_sim <- data
  data_sim$lgdp <- as.vector(lgdp_sim)
  return(data_sim)
}

# Analytical bias correction function (from provided code)
abc <- function(data, form, lags, N) {
  data$l1lgdp <- lag(data$lgdp, 1)
  data$l2lgdp <- lag(data$lgdp, 2)
  data$l3lgdp <- lag(data$lgdp, 3)
  data$l4lgdp <- lag(data$lgdp, 4)
  fit <- lm(form, data, x = TRUE, na.action = na.omit)
  res <- fit$residuals
  jac <- solve(t(fit$x) %*% fit$x / length(res))[2:6, 2:6]
  indexes <- c(1:length(res))
  bscore <- rep(0, 5)
  T_eff <- length(res) / N
  for (i in 1:lags) {
    indexes <- indexes[-c(1 + c(0:(N - 1)) * T_eff)]
    lindexes <- indexes - i
    bscore <- bscore + t(fit$x[indexes, 2:6]) %*% res[lindexes] / length(indexes)
  }
  bias <- -jac %*% bscore
  return(as.vector(bias / T_eff))
}

# Define the ABC formula
form.abc <- lgdp ~ dem + l1lgdp + l2lgdp + l3lgdp + l4lgdp + factor(year) + factor(id)

# Function to compute estimates on simulated data
compute_estimates <- function(data_sim) {
  # FE estimate
  fe.fit_sim <- plm(form.fe, data_sim, model = "within", effect = "twoways", index = c("id", "year"))
  coefs_fe_sim <- coef(fe.fit_sim)
  
  # JBC estimate
  fe.fit1 <- plm(form.fe, data_sim, subset = (as.double(year) <= 14), model = "within", effect = "twoways", index = c("id", "year"))
  fe.fit2 <- plm(form.fe, data_sim, subset = (as.double(year) >= 14), model = "within", effect = "twoways", index = c("id", "year"))
  coefs_jbc <- 19 * coefs_fe_sim / 9 - 10 * (coef(fe.fit1) + coef(fe.fit2)) / 18
  
  # ABC4 estimate
  bias_l4 <- abc(data_sim, form.abc, lags = 4, N = N)
  coefs_abc4 <- coefs_fe_sim - bias_l4
  
  return(list(coefs_fe_sim = coefs_fe_sim, coefs_jbc = coefs_jbc, coefs_abc4 = coefs_abc4))
}

# Run Monte Carlo simulation with 500 replications
R <- 500
results <- replicate(R, {
  data_sim <- simulate_data()
  compute_estimates(data_sim)
}, simplify = FALSE)

# Extract estimates
coefs_fe_all <- t(sapply(results, function(x) x$coefs_fe_sim))
coefs_jbc_all <- t(sapply(results, function(x) x$coefs_jbc))
coefs_abc4_all <- t(sapply(results, function(x) x$coefs_abc4))

# Compute bias, standard deviation, and RMSE
true_coefs <- coefs.fe

bias_fe <- colMeans(coefs_fe_all) - true_coefs
sd_fe <- apply(coefs_fe_all, 2, sd)
rmse_fe <- sqrt(colMeans((coefs_fe_all - matrix(true_coefs, nrow = R, ncol = length(true_coefs), byrow = TRUE))^2))

bias_jbc <- colMeans(coefs_jbc_all) - true_coefs
sd_jbc <- apply(coefs_jbc_all, 2, sd)
rmse_jbc <- sqrt(colMeans((coefs_jbc_all - matrix(true_coefs, nrow = R, ncol = length(true_coefs), byrow = TRUE))^2))

bias_abc4 <- colMeans(coefs_abc4_all) - true_coefs
sd_abc4 <- apply(coefs_abc4_all, 2, sd)
rmse_abc4 <- sqrt(colMeans((coefs_abc4_all - matrix(true_coefs, nrow = R, ncol = length(true_coefs), byrow = TRUE))^2))



cat("Monte Carlo Results for the Coefficient of 'dem':\n")
cat("FE Bias:", bias_fe[1], "\n")
cat("FE SD:", sd_fe[1], "\n")
cat("FE RMSE:", rmse_fe[1], "\n")
cat("JBC Bias:", bias_jbc[1], "\n")
cat("JBC SD:", sd_jbc[1], "\n")
cat("JBC RMSE:", rmse_jbc[1], "\n")
cat("ABC4 Bias:", bias_abc4[1], "\n")
cat("ABC4 SD:", sd_abc4[1], "\n")
cat("ABC4 RMSE:", rmse_abc4[1], "\n")