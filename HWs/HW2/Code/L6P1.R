### BS for Hansen-Singleton

library(AER)
library(sandwich)
library(gmm)
library(quantmod)
library(boot)

# Load and process data (unchanged from original)
econ_dataset <- read.csv("../../../Data/ccapm-long.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
attach(econ_dataset)
consumption_per_capita <- PCEND / (PPCEND * CNP16OV)
inflation_yearly <- PPCEND / Lag(PPCEND, k = 12) - 1
bond_real_yield <- (1 + GS1 / 100 - inflation_yearly)^(1 / 12)
stock_real_yield <- (SP500 / Lag(SP500)) * (Lag(PPCEND) / PPCEND)
consumption_change <- consumption_per_capita / Lag(consumption_per_capita)
returns_data <- na.omit(cbind(consumption_change, stock_real_yield, bond_real_yield))
colnames(returns_data) <- c("Cons_Growth", "Stock_Yield", "Bond_Yield")


# Moment conditions function (unchanged)
compute_moments <- function(params, input_data) {
  endogenous <- input_data[, 1:3]
  instr <- input_data[, -(1:3)]
  beta <- params[1]
  gamma <- params[2]
  res_stock <- beta * endogenous[, 2] * endogenous[, 1]^(-gamma) - 1
  res_bond <- beta * endogenous[, 3] * endogenous[, 1]^(-gamma) - 1
  moments_stock <- res_stock * instr
  moments_bond <- res_bond * instr
  return(cbind(moments_stock, moments_bond))
}

# Prepare data (unchanged)
observed_vars <- na.omit(cbind(consumption_change, bond_real_yield, stock_real_yield))
lag_count <- 1
instr_basic <- matrix(1, nrow = nrow(observed_vars), ncol = 1)
model_data_basic <- na.omit(cbind(observed_vars, instr_basic))
instr_extended <- cbind(
  1,
  Lag(observed_vars[, 2], lag_count),
  Lag(observed_vars[, 3], lag_count)
)
model_data_extended <- na.omit(cbind(observed_vars, instr_extended))


initial_params <- c(0.99, 1)
# gmm <- gmm(compute_moments, model_data_extended, t0 = initial_params, method = "BFGS", type = "iter", vcov = 'iid')


# Load required libraries
library(gmm)
library(boot)

# Function for Empirical Bootstrap
gmm_bs_statistic <- function(data, indices) {
  # Subset data using bootstrap indices
  boot_data <- data[indices, ]

  # Fit GMM model to the bootstrap sample and extract coefficients
  gmm_coef <- gmm(compute_moments, boot_data, t0 = c(0.99, 0.22), type = "iterative", vcov = "iid")$coef

  # Return the estimated coefficients
  return(gmm_coef)
}

# Set seed for reproducibility
set.seed(1)

# Perform bootstrap with 400 iterations
result_gmm_bs <- boot(data = model_data_basic, statistic = gmm_bs_statistic, R = 500)

# Calculate standard errors from bootstrap results
ses <- round(sqrt(apply(result_gmm_bs$t, 2, var)), digits = 4)


ci_param1 <- quantile(result_gmm_bs$t[, 1], c(0.05, 0.95))
ci_param2 <- quantile(result_gmm_bs$t[, 2], c(0.05, 0.95))
legend_text_param1 <- paste("SE =", round(ses[1], 4), "\nCI (90%) = [", round(ci_param1[1], 4), ", ", round(ci_param1[2], 4), "]")
legend_text_param2 <- paste("SE =", round(ses[2], 4), "\nCI (90%) = [", round(ci_param2[1], 4), ", ", round(ci_param2[2], 4), "]")

# Open PDF device
pdf("bootstrap_histograms.pdf", width = 8, height = 6)

# Set layout for two plots (2 rows, 1 column)
par(mfrow = c(2, 1))

# Histogram for Parameter 1
hist(result_gmm_bs$t[, 1], breaks = 23, col = "lightblue",
     main = "Histogram for BS Draws (Parameter beta)", xlab = "Values of Statistic")
legend("topright", legend = legend_text_param1, bty = "n")

# Histogram for Parameter 2
hist(result_gmm_bs$t[, 2], breaks = 23, col = "lightblue",
     main = "Histogram for BS Draws (Parameter alpha)", xlab = "Values of Statistic")
legend("topright", legend = legend_text_param2, bty = "n")

# Close PDF device
dev.off()


# 
# ## Load required libraries
# library(gmm)    # For GMM and CUE estimation
# library(boot)   # For bootstrap
# library(xtable) # For table output
# 
# R <- 500
# # Load and process data (unchanged from original)
# econ_dataset <- read.csv("../../../Data/ccapm-long.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
# attach(econ_dataset)
# consumption_per_capita <- PCEND / (PPCEND * CNP16OV)
# inflation_yearly <- PPCEND / Lag(PPCEND, k = 12) - 1
# bond_real_yield <- (1 + GS1 / 100 - inflation_yearly)^(1 / 12)
# stock_real_yield <- (SP500 / Lag(SP500)) * (Lag(PPCEND) / PPCEND)
# consumption_change <- consumption_per_capita / Lag(consumption_per_capita)
# returns_data <- na.omit(cbind(consumption_change, stock_real_yield, bond_real_yield))
# colnames(returns_data) <- c("Cons_Growth", "Stock_Yield", "Bond_Yield")
# 
# 
# 
# compute_moments <- function(params, input_data) {
#   endogenous <- input_data[, 1:3]
#   instr <- input_data[, -(1:3)]
#   beta <- params[1]
#   gamma <- params[2]
#   res_stock <- beta * endogenous[, 2] * endogenous[, 1]^(-gamma) - 1
#   res_bond <- beta * endogenous[, 3] * endogenous[, 1]^(-gamma) - 1
#   moments_stock <- res_stock * instr
#   moments_bond <- res_bond * instr
#   return(cbind(moments_stock, moments_bond))
# }
# 
# # Function to run GMM or CUE estimation with specified vcov
# run_estimation <- function(data, method, vcov_opt, t0) {
#   if (method == "GMM") {
#     fit <- gmm(compute_moments, data, t0 = t0, type = "iter", vcov = vcov_opt)
#   } else if (method == "CUE") {
#     fit <- gmm(compute_moments, data, t0 = t0, type = "cue", vcov = vcov_opt)
#   }
#   summary_fit <- summary(fit)
#   se <- summary_fit$coefficients[, "Std. Error"]
#   return(se)
# }
# 
# # Bootstrap function to compute standard errors
# bootstrap_se <- function(data, method, t0, R = 400) {
#   boot_stat <- if (method == "GMM") {
#     function(d, i) coef(gmm(compute_moments, d[i, ], t0 = t0, type = "iter", vcov = "iid"))
#   } else {
#     function(d, i) coef(gmm(compute_moments, d[i, ], t0 = t0, type = "cue", vcov = "iid"))
#   }
#   boot_result <- boot(data, boot_stat, R = R)
#   se <- apply(boot_result$t, 2, sd)
#   return(se)
# }
# 
# 
# 
# # Prepare data (unchanged)
# observed_vars <- na.omit(cbind(consumption_change, bond_real_yield, stock_real_yield))
# lag_count <- 1
# instr_basic <- matrix(1, nrow = nrow(observed_vars), ncol = 1)
# model_data_basic <- na.omit(cbind(observed_vars, instr_basic))
# instr_extended <- cbind(
#   1,
#   Lag(observed_vars[, 2], lag_count),
#   Lag(observed_vars[, 3], lag_count)
# )
# model_data_extended <- na.omit(cbind(observed_vars, instr_extended))
# 
# 
# # Define sample sizes for each instrument set
# # For basic instruments
# M_basic <- nrow(model_data_basic)
# sample_sizes_basic <- list(
#   "full" = model_data_basic
# )
# 
# # For extended instruments
# M_extended <- nrow(model_data_extended)
# sample_sizes_extended <- list(
#   "full" = model_data_extended
# )
# 
# # Define methods and vcov options
# methods <- c("GMM")
# vcov_options <- c("iid", "HAC", "MDS")
# 
# # Initialize results storage
# results <- list()
# 
# # Loop over sample sizes
# for (size_name in c("full")) {
#   results[[size_name]] <- list()
#   # Loop over instrument sets
#   for (instr_name in c("basic", "extended")) {
#     # Select the appropriate sample size dataset based on instrument set
#     if (instr_name == "basic") {
#       data <- sample_sizes_basic[[size_name]]
#     } else {
#       data <- sample_sizes_extended[[size_name]]
#     }
#     results[[size_name]][[instr_name]] <- list()
#     # Loop over estimation methods
#     for (method in methods) {
#       # Initial parameter values (customize as needed)
#       t0 <- if (method == "CUE") c(0.99, 0.5) else c(0.99, 1)
#       # Compute bootstrap standard errors
#       boot_se <- bootstrap_se(data, method, t0)
#       results[[size_name]][[instr_name]][[method]]$boot_se <- boot_se
#       # Compute vcov-based standard errors
#       for (vcov_opt in vcov_options) {
#         se <- run_estimation(data, method, vcov_opt, t0)
#         results[[size_name]][[instr_name]][[method]][[vcov_opt]] <- se
#       }
#     }
#   }
# }
# 
# # Generate comparison tables
# for (size_name in names(results)) {
#   for (vcov_opt in vcov_options) {
#     table_data <- data.frame()
#     for (instr_name in c("basic", "extended")) {
#       for (method in methods) {
#         boot_se <- results[[size_name]][[instr_name]][[method]]$boot_se
#         vcov_se <- results[[size_name]][[instr_name]][[method]][[vcov_opt]]
#         diff <- boot_se - vcov_se
#         # Create rows for each parameter (assuming 2 parameters: beta, gamma)
#         row <- data.frame(
#           Method = method,
#           Instruments = instr_name,
#           Parameter = c("beta", "gamma"),
#           Boot_SE = boot_se,
#           vcov_SE = vcov_se,
#           Diff = diff
#         )
#         table_data <- rbind(table_data, row)
#       }
#     }
#     # Print table for each sample size and vcov option
#     cat(sprintf("\nTable for Sample Size: %s, vcov: %s\n", size_name, vcov_opt))
#     print(xtable(table_data,
#                  caption = sprintf("Comparison for Sample Size: %s, vcov: %s", size_name, vcov_opt),
#                  digits = 4))
#   }
# }