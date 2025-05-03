# Load required packages
library(foreign)
library(xtable)
library(plm)
library(gmm)
library(readstata13)

# Load and prepare panel dataset
input_data <- read.dta13("../../Data/democracy-balanced-l4.dta")
input_data <- pdata.frame(input_data, index = c("id", "year"))

# Compute summary statistics
options(digits = 2)
stats_matrix <- matrix(NA, nrow = 3, ncol = 4,
                       dimnames = list(c("Democracy", "Log GDP", "Sample Size"),
                                       c("Mean", "Std. Dev.", "Dem = 1", "Dem = 0")))
stats_matrix[1:2, 1] <- colMeans(input_data[, c("dem", "lgdp")])
stats_matrix[1:2, 2] <- apply(input_data[, c("dem", "lgdp")], 2, sd)
stats_matrix[1:2, 3] <- colMeans(input_data[input_data$dem == 1, c("dem", "lgdp")])
stats_matrix[1:2, 4] <- colMeans(input_data[input_data$dem == 0, c("dem", "lgdp")])
stats_matrix[3, ] <- c(nrow(input_data), nrow(input_data), 
                       sum(input_data$dem == 1), sum(input_data$dem == 0))
xtable(stats_matrix)

# Anderson-Hsiao GMM Estimation
ah_formula <- lgdp ~ dem + lag(lgdp, 1:4) | lag(lgdp, 2) + lag(dem, 1)
ah_model <- pgmm(ah_formula, data = input_data, model = "twosteps", effect = "twoways", robust = TRUE)
ah_coeffs <- coef(ah_model)
ah_vcov <- vcovHC(ah_model, cluster = "group")
ah_cse <- sqrt(diag(ah_vcov))
ah_jtest <- sargan(ah_model)$statistic
ah_jdof <- sargan(ah_model)$parameter
ah_lr <- ah_coeffs[1] / (1 - sum(ah_coeffs[2:5]))
ah_jac_lr <- c(1, rep(ah_lr, 4)) / (1 - sum(ah_coeffs[2:5]))
ah_lr_cse <- sqrt(t(ah_jac_lr) %*% ah_vcov[1:5, 1:5] %*% ah_jac_lr)

# Anderson-Hsiao Split-Sample Bias Correction (1 partition)
set.seed(888)
num_splits_1 <- 1
num_units <- length(unique(input_data$id))
ah_split_coeffs <- rep(0, length(ah_coeffs))
ah_split_lr <- 0
for (s in 1:num_splits_1) {
  split_sample <- sample(num_units, ceiling(num_units / 2), replace = FALSE)
  ah_model_1 <- pgmm(ah_formula, data = input_data[as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways", robust = TRUE)
  ah_model_2 <- pgmm(ah_formula, data = input_data[!as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways", robust = TRUE)
  lr_1 <- coef(ah_model_1)[1] / (1 - sum(coef(ah_model_1)[2:5]))
  lr_2 <- coef(ah_model_2)[1] / (1 - sum(coef(ah_model_2)[2:5]))
  ah_split_coeffs <- ah_split_coeffs + ((coef(ah_model_1) + coef(ah_model_2)) / 2) / num_splits_1
  ah_split_lr <- ah_split_lr + ((lr_1 + lr_2) / 2) / num_splits_1
}
ah_jbc_coeffs <- 2 * ah_coeffs - ah_split_coeffs
ah_jbc_lr <- 2 * ah_lr - ah_split_lr

# Anderson-Hsiao Split-Sample Bias Correction (5 partitions)
num_splits_5 <- 5
ah_split_coeffs_5 <- rep(0, length(ah_coeffs))
ah_split_lr_5 <- 0
for (s in 1:num_splits_5) {
  split_sample <- sample(num_units, ceiling(num_units / 2), replace = FALSE)
  ah_model_1 <- pgmm(ah_formula, data = input_data[as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways", robust = TRUE)
  ah_model_2 <- pgmm(ah_formula, data = input_data[!as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways", robust = TRUE)
  lr_1 <- coef(ah_model_1)[1] / (1 - sum(coef(ah_model_1)[2:5]))
  lr_2 <- coef(ah_model_2)[1] / (1 - sum(coef(ah_model_2)[2:5]))
  ah_split_coeffs_5 <- ah_split_coeffs_5 + ((coef(ah_model_1) + coef(ah_model_2)) / 2) / num_splits_5
  ah_split_lr_5 <- ah_split_lr_5 + ((lr_1 + lr_2) / 2) / num_splits_5
}
ah_jbc5_coeffs <- 2 * ah_coeffs - ah_split_coeffs_5
ah_jbc5_lr <- 2 * ah_lr - ah_split_lr_5

# Arellano-Bond GMM Estimation
ab_formula <- lgdp ~ dem + lag(lgdp, 1:4) | lag(lgdp, 2:99) + lag(dem, 1:99)
ab_model <- pgmm(ab_formula, data = input_data, model = "twosteps", effect = "twoways")
ab_coeffs <- coef(ab_model)
ab_vcov <- vcovHC(ab_model, cluster = "group")
ab_cse <- sqrt(diag(ab_vcov))
ab_jtest <- sargan(ab_model)$statistic
ab_jdof <- sargan(ab_model)$parameter
ab_lr <- ab_coeffs[1] / (1 - sum(ab_coeffs[2:5]))
ab_jac_lr <- c(1, rep(ab_lr, 4)) / (1 - sum(ab_coeffs[2:5]))
ab_lr_cse <- sqrt(t(ab_jac_lr) %*% ab_vcov[1:5, 1:5] %*% ab_jac_lr)

# Arellano-Bond Split-Sample Bias Correction (1 partition)
num_splits_1 <- 1
ab_split_coeffs <- rep(0, length(ab_coeffs))
ab_split_lr <- 0
for (s in 1:num_splits_1) {
  split_sample <- sample(num_units, ceiling(num_units / 2), replace = FALSE)
  ab_model_1 <- pgmm(ab_formula, data = input_data[as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways")
  ab_model_2 <- pgmm(ab_formula, data = input_data[!as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways")
  lr_1 <- coef(ab_model_1)[1] / (1 - sum(coef(ab_model_1)[2:5]))
  lr_2 <- coef(ab_model_2)[1] / (1 - sum(coef(ab_model_2)[2:5]))
  ab_split_coeffs <- ab_split_coeffs + ((coef(ab_model_1) + coef(ab_model_2)) / 2) / num_splits_1
  ab_split_lr <- ab_split_lr + ((lr_1 + lr_2) / 2) / num_splits_1
}
ab_jbc_coeffs <- 2 * ab_coeffs - ab_split_coeffs
ab_jbc_lr <- 2 * ab_lr - ab_split_lr

# Arellano-Bond Split-Sample Bias Correction (5 partitions)
num_splits_5 <- 5
ab_split_coeffs_5 <- rep(0, length(ab_coeffs))
ab_split_lr_5 <- 0
for (s in 1:num_splits_5) {
  split_sample <- sample(num_units, ceiling(num_units / 2), replace = FALSE)
  ab_model_1 <- pgmm(ab_formula, data = input_data[as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways")
  ab_model_2 <- pgmm(ab_formula, data = input_data[!as.numeric(id) %in% split_sample, ], 
                     model = "twosteps", effect = "twoways")
  lr_1 <- coef(ab_model_1)[1] / (1 - sum(coef(ab_model_1)[2:5]))
  lr_2 <- coef(ab_model_2)[1] / (1 - sum(coef(ab_model_2)[2:5]))
  ab_split_coeffs_5 <- ab_split_coeffs_5 + ((coef(ab_model_1) + coef(ab_model_2)) / 2) / num_splits_5
  ab_split_lr_5 <- ab_split_lr_5 + ((lr_1 + lr_2) / 2) / num_splits_5
}
ab_jbc5_coeffs <- 2 * ab_coeffs - ab_split_coeffs_5
ab_jbc5_lr <- 2 * ab_lr - ab_split_lr_5

# Bootstrap Standard Errors
boot_reps <- 500

# Function to generate bootstrap datasets
generate_bootstrap <- function(data, params) {
  num_ids <- length(unique(data$id))
  num_years <- length(unique(data$year))
  sampled_ids <- kronecker(sample.int(num_ids, num_ids, replace = TRUE), rep(1, num_years))
  indices <- (sampled_ids - 1) * num_years + rep(seq_len(num_years), num_ids)
  boot_data <- data[indices, ]
  boot_data$id <- kronecker(seq_len(num_ids), rep(1, num_years))
  boot_data$year <- rep(1987:2009, num_ids)
  boot_data <- as.data.frame(boot_data)
  boot_data <- pdata.frame(boot_data, index = c("id", "year"))
  return(boot_data)
}

# Function to compute statistics for each bootstrap draw
bootstrap_stats <- function(data, form_fe, form_ah, form_ab) {
  # Anderson-Hsiao
  ah_fit <- pgmm(form_ah, data, model = "twosteps", effect = "twoways", robust = TRUE)
  ah_coeffs <- coef(ah_fit)
  ah_lr_val <- ah_coeffs[1] / (1 - sum(ah_coeffs[2:5]))
  
  num_ids <- length(unique(data$id))
  ah_split_1 <- rep(0, length(ah_coeffs))
  ah_lr_1 <- 0
  for (s in 1:1) {
    split_ids <- sample(num_ids, ceiling(num_ids / 2), replace = FALSE)
    ah_sub1 <- pgmm(form_ah, data[as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways", robust = TRUE)
    ah_sub2 <- pgmm(form_ah, data[!as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways", robust = TRUE)
    lr_sub1 <- coef(ah_sub1)[1] / (1 - sum(coef(ah_sub1)[2:5]))
    lr_sub2 <- coef(ah_sub2)[1] / (1 - sum(coef(ah_sub2)[2:5]))
    ah_split_1 <- ah_split_1 + ((coef(ah_sub1) + coef(ah_sub2)) / 2)
    ah_lr_1 <- ah_lr_1 + ((lr_sub1 + lr_sub2) / 2)
  }
  ah_jbc <- 2 * ah_coeffs - ah_split_1
  ah_jbc_lr <- 2 * ah_lr_val - ah_lr_1
  
  ah_split_5 <- rep(0, length(ah_coeffs))
  ah_lr_5 <- 0
  for (s in 1:5) {
    split_ids <- sample(num_ids, ceiling(num_ids / 2), replace = FALSE)
    ah_sub1 <- pgmm(form_ah, data[as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways", robust = TRUE)
    ah_sub2 <- pgmm(form_ah, data[!as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways", robust = TRUE)
    lr_sub1 <- coef(ah_sub1)[1] / (1 - sum(coef(ah_sub1)[2:5]))
    lr_sub2 <- coef(ah_sub2)[1] / (1 - sum(coef(ah_sub2)[2:5]))
    ah_split_5 <- ah_split_5 + ((coef(ah_sub1) + coef(ah_sub2)) / 2) / 5
    ah_lr_5 <- ah_lr_5 + ((lr_sub1 + lr_sub2) / 2) / 5
  }
  ah_jbc5 <- 2 * ah_coeffs - ah_split_5
  ah_jbc5_lr <- 2 * ah_lr_val - ah_lr_5
  
  # Arellano-Bond
  ab_fit <- pgmm(form_ab, data, model = "twosteps", effect = "twoways")
  ab_coeffs <- coef(ab_fit)
  ab_lr_val <- ab_coeffs[1] / (1 - sum(ab_coeffs[2:5]))
  
  ab_split_1 <- rep(0, length(ab_coeffs))
  ab_lr_1 <- 0
  for (s in 1:1) {
    split_ids <- sample(num_ids, ceiling(num_ids / 2), replace = FALSE)
    ab_sub1 <- pgmm(form_ab, data[as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways")
    ab_sub2 <- pgmm(form_ab, data[!as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways")
    lr_sub1 <- coef(ab_sub1)[1] / (1 - sum(coef(ab_sub1)[2:5]))
    lr_sub2 <- coef(ab_sub2)[1] / (1 - sum(coef(ab_sub2)[2:5]))
    ab_split_1 <- ab_split_1 + ((coef(ab_sub1) + coef(ab_sub2)) / 2)
    ab_lr_1 <- ab_lr_1 + ((lr_sub1 + lr_sub2) / 2)
  }
  ab_jbc <- 2 * ab_coeffs - ab_split_1
  ab_jbc_lr <- 2 * ab_lr_val - ab_lr_1
  
  ab_split_5 <- rep(0, length(ab_coeffs))
  ab_lr_5 <- 0
  for (s in 1:5) {
    split_ids <- sample(num_ids, ceiling(num_ids / 2), replace = FALSE)
    ab_sub1 <- pgmm(form_ab, data[as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways")
    ab_sub2 <- pgmm(form_ab, data[!as.numeric(id) %in% split_ids, ], 
                    model = "twosteps", effect = "twoways")
    lr_sub1 <- coef(ab_sub1)[1] / (1 - sum(coef(ab_sub1)[2:5]))
    lr_sub2 <- coef(ab_sub2)[1] / (1 - sum(coef(ab_sub2)[2:5]))
    ab_split_5 <- ab_split_5 + ((coef(ab_sub1) + coef(ab_sub2)) / 2) / 5
    ab_lr_5 <- ab_lr_5 + ((lr_sub1 + lr_sub2) / 2) / 5
  }
  ab_jbc5 <- 2 * ab_coeffs - ab_split_5
  ab_jbc5_lr <- 2 * ab_lr_val - ab_lr_5
  
  return(c(ah_coeffs[1:5], ah_jbc[1:5], ah_jbc5[1:5], ab_coeffs[1:5], 
           ab_jbc[1:5], ab_jbc5[1:5], ah_lr_val, ah_jbc_lr, ah_jbc5_lr, 
           ab_lr_val, ab_jbc_lr, ab_jbc5_lr))
}

# Run bootstrap
library(boot)
boot_output <- boot(data = input_data, statistic = bootstrap_stats, 
                    sim = "parametric", ran.gen = generate_bootstrap, 
                    mle = 0, form_fe = form.fe, form_ah = ah_formula, 
                    form_ab = ab_formula, parallel = "multicore", 
                    ncpus = 20, R = boot_reps)

# Compute robust bootstrap standard errors
robust_sd <- function(x) {
  (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)) / 
    (qnorm(0.75) - qnorm(0.25))
}

boot_results <- matrix(as.numeric(boot_output$t), nrow = nrow(boot_output$t))
bse_ah <- apply(boot_results[, 1:5], 2, robust_sd)
bse_ah_jbc <- apply(boot_results[, 6:10], 2, robust_sd)
bse_ah_jbc5 <- apply(boot_results[, 11:15], 2, robust_sd)
bse_ab <- apply(boot_results[, 16:20], 2, robust_sd)
bse_ab_jbc <- apply(boot_results[, 21:25], 2, robust_sd)
bse_ab_jbc5 <- apply(boot_results[, 26:30], 2, robust_sd)
bse_lr_ah <- robust_sd(boot_results[, 31])
bse_lr_ah_jbc <- robust_sd(boot_results[, 32])
bse_lr_ah_jbc5 <- robust_sd(boot_results[, 33])
bse_lr_ab <- robust_sd(boot_results[, 34])
bse_lr_ab_jbc <- robust_sd(boot_results[, 35])
bse_lr_ab_jbc5 <- robust_sd(boot_results[, 36])

# Construct results table
results_table <- matrix(NA, nrow = 18, ncol = 6,
                        dimnames = list(c("Democracy", "CSE", "BSE",
                                          "Lag 1 GDP", "CSE Lag 1", "BSE Lag 1",
                                          "Lag 2 GDP", "CSE Lag 2", "BSE Lag 2",
                                          "Lag 3 GDP", "CSE Lag 3", "BSE Lag 3",
                                          "Lag 4 GDP", "CSE Lag 4", "BSE Lag 4",
                                          "Long-Run Effect", "CSE Long-Run", "BSE Long-Run"),
                                        c("GMM1", "GMM1-BC", "GMM1-BC5", "GMM2", "GMM2-BC", "GMM2-BC5")))

# Populate GMM1 (Anderson-Hsiao) results
results_table[c(1, 4, 7, 10, 13), 1] <- ah_coeffs[1:5]
results_table[c(2, 5, 8, 11, 14), 1] <- ah_cse[1:5]
results_table[c(3, 6, 9, 12, 15), 1] <- bse_ah
results_table[16, 1] <- ah_lr
results_table[17, 1] <- ah_lr_cse
results_table[18, 1] <- bse_lr_ah

# Populate GMM1-BC (1 partition) results
results_table[c(1, 4, 7, 10, 13), 2] <- ah_jbc_coeffs[1:5]
results_table[c(3, 6, 9, 12, 15), 2] <- bse_ah_jbc
results_table[16, 2] <- ah_jbc_lr
results_table[18, 2] <- bse_lr_ah_jbc

# Populate GMM1-BC5 (5 partitions) results
results_table[c(1, 4, 7, 10, 13), 3] <- ah_jbc5_coeffs[1:5]
results_table[c(3, 6, 9, 12, 15), 3] <- bse_ah_jbc5
results_table[16, 3] <- ah_jbc5_lr
results_table[18, 3] <- bse_lr_ah_jbc5

# Populate GMM2 (Arellano-Bond) results
results_table[c(1, 4, 7, 10, 13), 4] <- ab_coeffs[1:5]
results_table[c(2, 5, 8, 11, 14), 4] <- ab_cse[1:5]
results_table[c(3, 6, 9, 12, 15), 4] <- bse_ab
results_table[16, 4] <- ab_lr
results_table[17, 4] <- ab_lr_cse
results_table[18, 4] <- bse_lr_ab

# Populate GMM2-BC (1 partition) results
results_table[c(1, 4, 7, 10, 13), 5] <- ab_jbc_coeffs[1:5]
results_table[c(3, 6, 9, 12, 15), 5] <- bse_ab_jbc
results_table[16, 5] <- ab_jbc_lr
results_table[18, 5] <- bse_lr_ab_jbc

# Populate GMM2-BC5 (5 partitions) results
results_table[c(1, 4, 7, 10, 13), 6] <- ab_jbc5_coeffs[1:5]
results_table[c(3, 6, 9, 12, 15), 6] <- bse_ab_jbc5
results_table[16, 6] <- ab_jbc5_lr
results_table[18, 6] <- bse_lr_ab_jbc5

# Scale specific rows
results_table[c(1, 2, 3, 16, 17, 18), ] <- 100 * results_table[c(1, 2, 3, 16, 17, 18), ]

# Output results table
xtable(results_table, digits = 2)
