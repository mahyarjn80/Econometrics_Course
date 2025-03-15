# Install and load packages
library(AER)
library(sandwich)
library(gmm)
library(quantmod)
library(xtable)

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

# Define vcov options and estimators
vcov_options <- c("iid", "HAC", "MDS")
estimators <- c("Basic_GMM", "Basic_CUE", "Extended_GMM", "Extended_CUE")

# Function to run GMM estimation
run_gmm_estimation <- function(estimator, vcov_opt, model_data, initial_params) {
  if (estimator == "Basic_GMM" || estimator == "Extended_GMM") {
    gmm_result <- gmm(compute_moments, model_data, t0 = initial_params, method = "BFGS", type = "iter", vcov = vcov_opt)
  } else {
    gmm_result <- gmm(compute_moments, model_data, t0 = initial_params, type = "cue", vcov = vcov_opt)
  }
  if (vcov_opt == "CL") {
    cluster_var <- 1:nrow(model_data) # Assume time-based clustering
    gmm_result <- update(gmm_result, vcov = "CL", cluster = cluster_var)
  }
  summary_result <- summary(gmm_result)
  estimates <- summary_result$coef[, 1]
  std_errors <- summary_result$coef[, 2]
  j_stat <- summary_result$stest$test[1]
  p_value <- summary_result$stest$test[2]
  return(list(estimates = estimates, std_errors = std_errors, j_stat = j_stat, p_value = p_value))
}

# Store results
results <- list()
for (est in estimators) {
  model_data <- if (grepl("Basic", est)) model_data_basic else model_data_extended
  initial_params <- if (est == "Basic_CUE") c(0.99, 0.5) else c(0.99, 1)
  results[[est]] <- list()
  for (vcov_opt in vcov_options) {
    results[[est]][[vcov_opt]] <- run_gmm_estimation(est, vcov_opt, model_data, initial_params)
  }
}

# Organize results into tables
create_results_df <- function(est_results) {
  df <- data.frame(
    Vcov = vcov_options,
    Beta = sapply(est_results, function(x) x$estimates[1]),
    Beta_SE = sapply(est_results, function(x) x$std_errors[1]),
    Gamma = sapply(est_results, function(x) x$estimates[2]),
    Gamma_SE = sapply(est_results, function(x) x$std_errors[2]),
    J_Stat = sapply(est_results, function(x) x$j_stat),
    P_Value = sapply(est_results, function(x) x$p_value)
  )
  return(df)
}

results_dfs <- lapply(results, create_results_df)

# Generate LaTeX tables
for (est in names(results_dfs)) {
  caption <- paste("Estimation Results for", est, "Across Vcov Methods")
  label <- paste("tab:", gsub("_", "", tolower(est)), sep = "")
  print(xtable(results_dfs[[est]], caption = caption, label = label, digits = 4), 
        include.rownames = FALSE, file = paste0(est, "_table.tex"))
}

