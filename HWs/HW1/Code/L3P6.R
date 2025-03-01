# Define function for robust standard errors
robust_se <- function(model) {
  require(sandwich)
  vcovHC(model, type = "HC")
}

# Load required libraries
library(foreign)  # For read.dta
library(AER)      # For ivreg
library(lmtest)   # For coeftest

# Set working directory and load data
load(file = "~/Downloads/QOB_Census80_Cohort_30_39.Rdata")
alpha_level <- 0.05  # Significance level

dataset <- QOB_Census80_Cohort_30_39

# Partial out control variables
resid_y <- residuals(lm(LWKLYWGE ~ factor(YOB) + factor(MARRIED) + factor(RACE) + 
                          factor(SMSA) + factor(REGION), data = dataset))
resid_d <- residuals(lm(EDUC ~ factor(YOB) + factor(MARRIED) + factor(RACE) + 
                          factor(SMSA) + factor(REGION), data = dataset))
resid_z <- residuals(lm(I(QOB == 4) ~ factor(YOB) + factor(MARRIED) + factor(RACE) + 
                          factor(SMSA) + factor(REGION), data = dataset))

# OLS estimation
ols_model <- lm(resid_y ~ resid_z)
ols_beta <- coef(ols_model)["resid_z"]
ols_se <- coeftest(ols_model, vcov = robust_se)["resid_z", "Std. Error"]
ols_lower <- ols_beta + qnorm(alpha_level / 2) * ols_se
ols_upper <- ols_beta + qnorm(1 - alpha_level / 2) * ols_se

# First-stage regression
fs_model <- lm(resid_d ~ resid_z)
fs_beta <- coef(fs_model)["resid_z"]
fs_se <- coeftest(fs_model, vcov = robust_se)["resid_z", "Std. Error"]
fs_f_stat <- (fs_beta / fs_se)^2
cat("First-stage F-statistic:", round(fs_f_stat, 5), "\n")

# TSLS estimation
tsls_model <- ivreg(resid_y ~ resid_d | resid_z)
tsls_beta <- coef(tsls_model)["resid_d"]
tsls_se <- coeftest(tsls_model, vcov = robust_se)["resid_d", "Std. Error"]
tsls_lower <- tsls_beta + qnorm(alpha_level / 2) * tsls_se
tsls_upper <- tsls_beta + qnorm(1 - alpha_level / 2) * tsls_se

# Weak-instrument-robust confidence interval
beta_range <- tsls_beta + seq(-100, 100, length.out = 201) / 1500
s_values <- numeric(length(beta_range))
n_size <- length(resid_y)

for (idx in seq_along(beta_range)) {
  temp_resid <- resid_y - beta_range[idx] * resid_d
  s_values[idx] <- n_size * (mean(temp_resid * resid_z))^2 / var(temp_resid * resid_z)
}

crit_threshold <- qchisq(1 - alpha_level, 1)
robust_lower <- min(beta_range[s_values < crit_threshold])
robust_upper <- max(beta_range[s_values < crit_threshold])

# Construct results table
results_table <- data.frame(
  Method = c("OLS", "TSLS", "Robust"),
  Estimate = c(ols_beta, tsls_beta, (robust_lower + robust_upper) / 2),
  Std_Error = c(ols_se, tsls_se, (robust_upper - robust_lower) / (2 * qnorm(1 - alpha_level / 2))),
  Lower_CI = c(ols_lower, tsls_lower, robust_lower),
  Upper_CI = c(ols_upper, tsls_upper, robust_upper)
)

# Format numeric columns to 5 decimal places
results_table[, 2:5] <- round(as.matrix(results_table[, 2:5]), 5)

# Display results
print(results_table)