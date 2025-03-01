# Function for robust standard errors
robust_se <- function(model) {
  require(sandwich)
  vcovHC(model, type = "HC3")
}

# Load libraries
library(foreign)  # For read.table
library(AER)      # For ivreg
library(lmtest)   # For coeftest

# Load data
my_data <- read.table("~/Downloads/acemoglu_col1.txt", header = TRUE)
alpha <- 0.05  # Significance level

# OLS estimation
ols_mod <- lm(GDP ~ Exprop + Latitude, data = my_data)
ols_est <- coef(ols_mod)["Exprop"]
ols_se <- coeftest(ols_mod, vcov = robust_se)["Exprop", "Std. Error"]
ols_ci_lower <- ols_est + qnorm(alpha / 2) * ols_se
ols_ci_upper <- ols_est + qnorm(1 - alpha / 2) * ols_se

# First-stage regression
fs_mod <- lm(Exprop ~ log(Mort) + Latitude, data = my_data)
fs_est <- coef(fs_mod)["log(Mort)"]
fs_se <- coeftest(fs_mod, vcov = robust_se)["log(Mort)", "Std. Error"]
fs_fstat <- (fs_est / fs_se)^2
cat("First-stage F-statistic:", round(fs_fstat, 1), "\n")

# TSLS estimation
tsls_mod <- ivreg(GDP ~ Exprop + Latitude | log(Mort) + Latitude, data = my_data)
tsls_est <- coef(tsls_mod)["Exprop"]
tsls_se <- coeftest(tsls_mod, vcov = robust_se)["Exprop", "Std. Error"]
tsls_ci_lower <- tsls_est + qnorm(alpha / 2) * tsls_se
tsls_ci_upper <- tsls_est + qnorm(1 - alpha / 2) * tsls_se

# Weak-instrument-robust CI
beta_grid <- seq(tsls_est - 4 * tsls_se, tsls_est + 6 * tsls_se, length.out = 300)
s_stats <- numeric(length(beta_grid))

# Partial out Latitude
y_res <- residuals(lm(GDP ~ Latitude, data = my_data))
d_res <- residuals(lm(Exprop ~ Latitude, data = my_data))
z_res <- residuals(lm(log(Mort) ~ Latitude, data = my_data))
n_obs <- nrow(my_data)

for (i in seq_along(beta_grid)) {
  resid_term <- y_res - beta_grid[i] * d_res
  s_stats[i] <- n_obs * (mean(resid_term * z_res))^2 / var(resid_term * z_res)
}

crit_val <- qchisq(1 - alpha, 1)
robust_ci_lower <- min(beta_grid[s_stats < crit_val])
robust_ci_upper <- max(beta_grid[s_stats < crit_val])

# Build results table
results <- data.frame(
  Method = c("OLS", "TSLS", "Robust"),
  Estimate = c(ols_est, tsls_est, (robust_ci_lower + robust_ci_upper) / 2),
  Std_Error = c(ols_se, tsls_se, (robust_ci_upper - robust_ci_lower) / (2 * qnorm(1 - alpha / 2))),
  Lower_CI = c(ols_ci_lower, tsls_ci_lower, robust_ci_lower),
  Upper_CI = c(ols_ci_upper, tsls_ci_upper, robust_ci_upper)
)

# Format numeric columns
results[, 2:5] <- round(as.matrix(results[, 2:5]), 3)

# Print results
print(results)
