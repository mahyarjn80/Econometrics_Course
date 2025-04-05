# Clear all objects from memory
rm(list = ls())

# Define the path to the mortgage dataset
data_path <- "../../Data/mortgage.dta"

# Load required libraries
library(haven)           # For reading Stata (.dta) files
library(xtable)          # For creating LaTeX tables
library(SortedEffects)   # For sorted effects analysis, includes mortgage data

# Import the dataset
mortgage_data <- read_dta(data_path)

# Make dataset variables accessible directly
attach(mortgage_data)

# Compute descriptive statistics
options(digits = 3)
stats_summary <- cbind(
  colMeans(mortgage_data), 
  colMeans(mortgage_data[mortgage_data$black == 1, ]), 
  colMeans(mortgage_data[mortgage_data$black == 0, ])
)
xtable(stats_summary)

# Fit OLS regression models
model_formula_simple <- deny ~ black
model_formula_full <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med + ltv_high + denpmi + selfemp + single + hischl

ols_simple <- lm(model_formula_simple, data = mortgage_data)
ols_full <- lm(model_formula_full, data = mortgage_data)

xtable(summary(ols_simple, digits = 3))
xtable(summary(ols_full, digits = 3))

# Split data into training and validation sets
set.seed(1)
total_rows <- nrow(mortgage_data)
validation_indices <- sample(1:total_rows, floor(total_rows / 3), replace = FALSE)

training_data <- mortgage_data[-validation_indices, ]
validation_data <- mortgage_data[validation_indices, ]

# Attach training data for direct variable access
attach(training_data)

# Part 1: Compare Basic Models
# Linear Probability Model
linear_model <- lm(model_formula_full)

# Models with different link functions
logistic_model <- glm(model_formula_full, family = binomial(link = "logit"))
probit_model <- glm(model_formula_full, family = binomial(link = "probit"))
cauchit_model <- glm(model_formula_full, family = binomial(link = "cauchit"))

# Generate comparison plots
# Plot 1: Logit vs Linear and Cauchit
plot(predict(logistic_model, type = "response"), predict(linear_model, type = "response"),
     xlim = c(0, 1), ylim = c(-0.1, 1.20), xlab = "Logit Predictions", ylab = "Linear & Cauchit Predictions",
     type = "p", pch = 1, col = "red")
points(predict(logistic_model, type = "response"), predict(cauchit_model, type = "response"),
       pch = 1, col = "black")
abline(0, 1)
abline(h = 1)
abline(h = 0)
legend(x = -0.05, y = 1.25, bty = "n",
       legend = c("Logit vs Linear", "Logit vs Cauchit"),
       col = c("red", "black"), pch = c(1, 1), cex = 0.7)

# Plot 2: Logit vs Probit
plot(predict(logistic_model, type = "response"), predict(probit_model, type = "response"),
     xlim = c(0, 1), ylim = c(0, 1), xlab = "Logit Predictions", ylab = "Probit Predictions",
     col = 4)
abline(0, 1)

# Evaluate out-of-sample prediction performance
prediction_errors <- matrix(0, nrow = 1, ncol = 4)
prediction_errors[1, 1] <- sqrt(mean((validation_data$deny - predict(logistic_model, validation_data, type = "response"))^2))
prediction_errors[1, 2] <- sqrt(mean((validation_data$deny - predict(probit_model, validation_data, type = "response"))^2))
prediction_errors[1, 3] <- sqrt(mean((validation_data$deny - predict(cauchit_model, validation_data, type = "response"))^2))
prediction_errors[1, 4] <- sqrt(mean((validation_data$deny - predict(linear_model, validation_data, type = "response"))^2))

colnames(prediction_errors) <- c("Logit", "Probit", "Cauchit", "Linear")
rownames(prediction_errors) <- "Mean Square Prediction Error"
xtable(prediction_errors, digits = 4, align = rep("c", 5))

# Detach training data
detach(training_data)

# Part 2: Predicted Effects of 'black'
percentiles <- seq(2, 98) / 100  # Percentiles for SPE
sig_level <- 0.1                 # Significance level
tail_quantile <- 0.05            # Tail quantile for classification
boot_reps <- 500                 # Bootstrap replications

# SPE and APE for entire population using Probit
spe_analysis <- spe(fm = model_formula_full, data = mortgage_data, method = "probit", 
                    var = "black", us = percentiles, alpha = sig_level, b = boot_reps, 
                    parallel = TRUE, ncores = 24)
xtable(summary(spe_analysis, result = "average"), digits = 3)
plot(spe_analysis, ylim = c(0, 0.25), ylab = "Change in Probability", main = "", sub = "Probit Model")

# Classification Analysis for entire population
vars_of_interest <- c("deny", "black", "p_irat", "hse_inc", "ccred", "mcred", "pubrec", 
                      "denpmi", "ltv_med", "ltv_high", "selfemp", "single", "hischl")
ca_analysis <- ca(fm = model_formula_full, data = mortgage_data, var = "black", 
                  method = "logit", u = tail_quantile, alpha = sig_level, t = vars_of_interest, 
                  b = boot_reps, parallel = TRUE, ncores = 24)
xtable(summary(ca_analysis))

# SPE and APE for Black subpopulation using Probit
spe_black_analysis <- spe(fm = model_formula_full, data = mortgage_data, method = "probit", 
                          var = "black", us = percentiles, subgroup = (mortgage_data$black == 1), 
                          alpha = sig_level, b = boot_reps, parallel = TRUE, ncores = 24)
xtable(summary(spe_black_analysis, result = "average"), digits = 3)
plot(spe_black_analysis, ylim = c(0, 0.25), ylab = "Change in Probability", main = "", sub = "Probit Model")