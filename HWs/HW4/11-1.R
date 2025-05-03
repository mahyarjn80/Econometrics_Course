# Load necessary libraries
library(foreign)
library(xtable)
library(plm)
library(gmm)
library(readstata13)
library(boot)

# Load and set up panel dataset
data_panel <- read.dta13("../../Data/democracy-balanced-l4.dta")
data_panel <- pdata.frame(data_panel, index = c("id", "year"))

# Generate descriptive statistics
options(digits = 2)
summary_stats <- matrix(NA, nrow = 3, ncol = 4,
                        dimnames = list(c("Democracy", "Log GDP", "Observations"),
                                        c("Average", "Std. Dev.", "Dem = 1", "Dem = 0")))
summary_stats[1:2, 1] <- colMeans(data_panel[, c("dem", "lgdp")])
summary_stats[1:2, 2] <- apply(data_panel[, c("dem", "lgdp")], 2, sd)
summary_stats[1:2, 3] <- colMeans(data_panel[data_panel$dem == 1, c("dem", "lgdp")])
summary_stats[1:2, 4] <- colMeans(data_panel[data_panel$dem == 0, c("dem", "lgdp")])
summary_stats[3, ] <- c(nrow(data_panel), nrow(data_panel), 
                        sum(data_panel$dem == 1), sum(data_panel$dem == 0))
xtable(summary_stats)

# Fixed Effects Estimation
model_formula <- lgdp ~ dem + lag(lgdp, 1:4) - 1
fe_est <- plm(model_formula, data = data_panel, model = "within",
              effect = "twoways", index = c("id", "year"))
fe_params <- coef(fe_est)
cluster_vcov <- vcovHC(fe_est, cluster = "group")
fe_cse <- sqrt(diag(cluster_vcov))
lr_effect <- fe_params[1] / (1 - sum(fe_params[2:5]))
lr_jac <- c(1, rep(lr_effect, 4)) / (1 - sum(fe_params[2:5]))
lr_cse <- sqrt(t(lr_jac) %*% cluster_vcov[1:5, 1:5] %*% lr_jac)

# Jackknife Bias Correction
split_year <- 14
fe_split_early <- plm(model_formula, data = data_panel, subset = as.numeric(year) <= split_year,
                      model = "within", effect = "twoways", index = c("id", "year"))
fe_split_late <- plm(model_formula, data = data_panel, subset = as.numeric(year) >= split_year,
                     model = "within", effect = "twoways", index = c("id", "year"))
jbc_params <- 19 * fe_params / 9 - 10 * (coef(fe_split_early) + coef(fe_split_late)) / 18
lr_split_early <- coef(fe_split_early)[1] / (1 - sum(coef(fe_split_early)[2:5]))
lr_split_late <- coef(fe_split_late)[1] / (1 - sum(coef(fe_split_late)[2:5]))
jbc_lr_effect <- 19 * lr_effect / 9 - 10 * (lr_split_early + lr_split_late) / 18

# Bootstrap Standard Errors
n_boot <- 500

# Bootstrap data generation
generate_boot_data <- function(data, params) {
  n_units <- length(unique(data$id))
  n_periods <- length(unique(data$year))
  boot_ids <- kronecker(sample.int(n_units, n_units, replace = TRUE), rep(1, n_periods))
  indices <- (boot_ids - 1) * n_periods + rep(seq_len(n_periods), n_units)
  boot_data <- data[indices, ]
  boot_data$id <- kronecker(seq_len(n_units), rep(1, n_periods))
  boot_data$year <- rep(1987:2009, n_units)
  boot_data <- as.data.frame(boot_data)
  boot_data <- pdata.frame(boot_data, index = c("id", "year"))
  return(boot_data)
}

# Bootstrap statistics calculation
compute_boot_stats <- function(data, formula) {
  fe_fit <- plm(formula, data, model = "within", effect = "twoways", index = c("id", "year"))
  fe_coeffs <- coef(fe_fit)
  lr_fe <- fe_coeffs[1] / (1 - sum(fe_coeffs[2:5]))
  
  split1 <- plm(formula, data, subset = as.numeric(year) <= 14, model = "within",
                effect = "twoways", index = c("id", "year"))
  split2 <- plm(formula, data, subset = as.numeric(year) >= 14, model = "within",
                effect = "twoways", index = c("id", "year"))
  jbc_coeffs <- 19 * fe_coeffs / 9 - 10 * (coef(split1) + coef(split2)) / 18
  lr_split1 <- coef(split1)[1] / (1 - sum(coef(split1)[2:5]))
  lr_split2 <- coef(split2)[1] / (1 - sum(coef(split2)[2:5]))
  jbc_lr <- 19 * lr_fe / 9 - 10 * (lr_split1 + lr_split2) / 18
  
  return(c(fe_coeffs, jbc_coeffs, lr_fe, jbc_lr))
}

# Run bootstrap
boot_results <- boot(data = data_panel, statistic = compute_boot_stats, 
                     sim = "parametric", ran.gen = generate_boot_data, 
                     mle = 0, formula = model_formula, 
                     parallel = "multicore", ncpus = 4, R = n_boot)

# Compute robust bootstrap standard errors
robust_se <- function(x) {
  (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)) /
    (qnorm(0.75) - qnorm(0.25))
}

boot_matrix <- matrix(as.numeric(boot_results$t), nrow = nrow(boot_results$t))
bse_fe <- apply(boot_matrix[, 1:5], 2, robust_se)
bse_jbc <- apply(boot_matrix[, 6:10], 2, robust_se)
bse_lr_fe <- robust_se(boot_matrix[, 11])
bse_lr_jbc <- robust_se(boot_matrix[, 12])

# Build results table for FE and FE-JBC
results_table <- matrix(NA, nrow = 18, ncol = 2,
                        dimnames = list(c("Democracy", "Clustered SE", "Bootstrap SE",
                                          "Lag 1 GDP", "Clustered SE Lag 1", "Bootstrap SE Lag 1",
                                          "Lag 2 GDP", "Clustered SE Lag 2", "Bootstrap SE Lag 2",
                                          "Lag 3 GDP", "Clustered SE Lag 3", "Bootstrap SE Lag 3",
                                          "Lag 4 GDP", "Clustered SE Lag 4", "Bootstrap SE Lag 4",
                                          "Long-Run Effect", "Clustered SE Long-Run", "Bootstrap SE Long-Run"),
                                        c("FE", "FE-JBC")))

# Fill FE results
results_table[c(1, 4, 7, 10, 13), 1] <- fe_params
results_table[c(2, 5, 8, 11, 14), 1] <- fe_cse
results_table[c(3, 6, 9, 12, 15), 1] <- bse_fe
results_table[16, 1] <- lr_effect
results_table[17, 1] <- lr_cse
results_table[18, 1] <- bse_lr_fe

# Fill FE-JBC results
results_table[c(1, 4, 7, 10, 13), 2] <- jbc_params
results_table[c(3, 6, 9, 12, 15), 2] <- bse_jbc
results_table[16, 2] <- jbc_lr_effect
results_table[18, 2] <- bse_lr_jbc

# Scale democracy and long-run effect rows
results_table[c(1, 2, 3, 16, 17, 18), ] <- 100 * results_table[c(1, 2, 3, 16, 17, 18), ]

# Display results table
xtable(results_table, digits = 2)