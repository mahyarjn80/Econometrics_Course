# 
# 
# # Part I: Functions
# 
# # Function to compute White standard errors with small sample adjustment;
# 
# hc <- function(x) 
# {;
#   vcovHC(x, type = "HC3");
# };
# 
# ####################
# 
# # Part II: Main program
# library(lmtest)
# library(sandwich);
# library(AER);
# library(systemfit);
# library(gmm);
# library(systemfit)
# # Reading the data;
# data <- as.data.frame(read.csv("../../../Data/bigfish.csv", row.names = 1,header = TRUE, stringsAsFactors = FALSE));
# attach(data);
# alpha <- .05; # Significance level
# 
# # In this example it is more convenient not to partial-out the controls because the sample size is small. If you partial out, 
# # you would need to manually correct the small sample adjustments of the  White standard errors
# 
# # ols results;
# 
# # Demand:
# formula    <- qty ~ price + day1 + day2 + day3 + day4 + rainy + cold;
# ols.fit    <- lm(formula);
# ols.coef.d <- ols.fit$coef[2];
# ols.se.d   <- coeftest(ols.fit, vcov = hc)[2,2];
# ols.lci.d  <- ols.coef.d + qnorm(alpha/2)*ols.se.d;
# ols.uci.d  <- ols.coef.d + qnorm(1-alpha/2)*ols.se.d;
# 
# # Supply:
# formula    <- qty ~ price + stormy + mixed;
# ols.fit    <- lm(formula);
# ols.coef.s <- ols.fit$coef[2];
# ols.se.s   <- coeftest(ols.fit, vcov = hc)[2,2];
# ols.lci.s  <- ols.coef.s + qnorm(alpha/2)*ols.se.s;
# ols.uci.s  <- ols.coef.s + qnorm(1-alpha/2)*ols.se.s;
# 
# 
# # First stage;
# 
# # Demand;
# formula   <- price ~ stormy + mixed + day1 + day2 + day3 + day4 + rainy + cold;
# fs.fit   <- lm(formula);
# fs.Fstat <- waldtest(fs.fit, . ~ . - stormy - mixed, vcov = hc,  test = "F")$F[2];
# print(paste('F-stat: ', fs.Fstat));
# 
# # Supply;
# fs.Fstat <- waldtest(fs.fit, . ~ .  - day1 - day2 - day3 - day4 - rainy - cold, vcov = hc,  test = "F")$F[2];
# print(paste('F-stat: ', fs.Fstat));
# 
# # tsls results
# 
# # Demand;
# formula     <- qty ~ price + day1 + day2 + day3 + day4  + rainy + cold| stormy + mixed + day1 + day2 + day3 + day4  + rainy + cold;
# tsls.fit    <- ivreg(formula);
# tsls.coef.d <- tsls.fit$coef[2];
# tsls.se.d   <- coeftest(tsls.fit, vcov = hc)[2,2]; 
# tsls.lci.d  <- tsls.coef.d + qnorm(alpha/2)*tsls.se.d;
# tsls.uci.d  <- tsls.coef.d + qnorm(1-alpha/2)*tsls.se.d;
# 
# # Supply;
# formula     <- qty ~ price + stormy + mixed| day1 + day2 + day3 + day4  + rainy + cold + stormy + mixed;
# tsls.fit    <- ivreg(formula);
# tsls.coef.s <- tsls.fit$coef[2];
# tsls.se.s   <- coeftest(tsls.fit, vcov = hc)[2,2]; 
# tsls.lci.s  <- tsls.coef.s + qnorm(alpha/2)*tsls.se.s;
# tsls.uci.s  <- tsls.coef.s + qnorm(1-alpha/2)*tsls.se.s;
# 
# # 2-stage gmm results: equation by equation;
# 
# # Demand;
# instr   <- cbind(stormy, mixed,day1,day2,day3,day4,rainy,cold);
# formula <- qty ~ price + day1 + day2 + day3 + day4  + rainy + cold;
# gmm.fit <- gmm(formula, x = instr, vcov = "iid");
# gmm.coef.d.all<-gmm.fit$coefficients
# gmm.coef.d <- summary(gmm.fit)$coefficients[2,1];
# gmm.se.d   <- summary(gmm.fit)$coefficients[2,2];
# gmm.lci.d  <- gmm.coef.d + qnorm(alpha/2)*gmm.se.d;
# gmm.uci.d  <- gmm.coef.d + qnorm(1-alpha/2)*gmm.se.d;
# print(summary(gmm.fit)$stest);
# 
# 
# 
# # Supply;
# instr   <- cbind(stormy, mixed,day1,day2,day3,day4,rainy,cold);
# formula <- qty ~ price + stormy + mixed;
# gmm.fit <- gmm(formula, x = instr, vcov = "iid");
# gmm.coef.s.all<-gmm.fit$coefficients
# gmm.coef.s <- summary(gmm.fit)$coefficients[2,1];
# gmm.se.s   <- summary(gmm.fit)$coefficients[2,2];
# gmm.lci.s  <- gmm.coef.s + qnorm(alpha/2)*gmm.se.s;
# gmm.uci.s  <- gmm.coef.s + qnorm(1-alpha/2)*gmm.se.s;
# # print(summary(gmm.fit)$stest);
# 
# 
# 
# 
# ## 3SLS: systemfit
# 
# 
# formula.d    <- qty ~ price + day1 + day2 + day3 + day4 + rainy + cold;
# formula.s    <- qty ~ price + stormy + mixed;
# formula      <- list(demand = formula.d, supply = formula.s);
# sols.fit     <- systemfit(formula);
# s3sls.fit     <- systemfit(formula, method = "3SLS", inst = ~ + day1 + day2 + day3 + day4 + rainy + cold + stormy + mixed);
# 
# # Demand;
# s3sls.coef.d <- s3sls.fit$coefficients[2];
# s3sls.se.d   <- summary(s3sls.fit)$coefficients[2,2]; 
# s3sls.lci.d  <- s3sls.coef.d + qnorm(alpha/2)*s3sls.se.d;
# s3sls.uci.d  <- s3sls.coef.d + qnorm(1-alpha/2)*s3sls.se.d;
# 
# # Supply;
# s3sls.coef.s <- s3sls.fit$coefficients[10];
# s3sls.se.s   <- summary(s3sls.fit)$coefficients[10,2]; 
# s3sls.lci.s  <- s3sls.coef.s + qnorm(alpha/2)*s3sls.se.s;
# s3sls.uci.s  <- s3sls.coef.s + qnorm(1-alpha/2)*s3sls.se.s;
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# options(digits=2);
# 
# table <- matrix(0, ncol = 4, nrow = 4, dimnames = list(c('LI-OLS', 'LI-TSLS', 'LI-3SLS', 'FI-3SLS'), c('Est.', 'Std. Error', '95% LCI','95% UCI')));
# 
# table[1,1] <- ols.coef.d;
# table[1,2] <- ols.se.d;
# table[1,3] <- ols.lci.d;
# table[1,4] <- ols.uci.d;
# 
# table[2,1] <- tsls.coef.d;
# table[2,2] <- tsls.se.d;
# table[2,3] <- tsls.lci.d;
# table[2,4] <- tsls.uci.d;
# 
# table[3,1] <- gmm.coef.d;
# table[3,2] <- gmm.se.d;
# table[3,3] <- gmm.lci.d;
# table[3,4] <- gmm.uci.d;
# 
# table[4,1] <- s3sls.coef.d;
# table[4,2] <- s3sls.se.d;
# table[4,3] <- s3sls.lci.d;
# table[4,4] <- s3sls.uci.d;
# 
# 
# print(table);
# 
# 
# table <- matrix(0, ncol = 4, nrow = 4, dimnames = list(c('LI-OLS', 'LI-TSLS', 'LI-3SLS', 'FI-3SLS'), c('Est.', 'Std. Error', '95% LCI','95% UCI')));
# 
# table[1,1] <- ols.coef.s;
# table[1,2] <- ols.se.s;
# table[1,3] <- ols.lci.s;
# table[1,4] <- ols.uci.s;
# 
# table[2,1] <- tsls.coef.s;
# table[2,2] <- tsls.se.s;
# table[2,3] <- tsls.lci.s;
# table[2,4] <- tsls.uci.s;
# 
# table[3,1] <- gmm.coef.s;
# table[3,2] <- gmm.se.s;
# table[3,3] <- gmm.lci.s;
# table[3,4] <- gmm.uci.s;
# 
# table[4,1] <- s3sls.coef.s;
# table[4,2] <- s3sls.se.s;
# table[4,3] <- s3sls.lci.s;
# table[4,4] <- s3sls.uci.s;
# 
# 
# 
# print(table);


# Part I: Functions

# Function to compute White standard errors with small sample adjustment
hc <- function(x) {
  vcovHC(x, type = "HC3")
}

# Part II: Main Program
library(lmtest)
library(sandwich)
library(AER)
library(systemfit)
library(gmm)
library(boot)  # Added for bootstrap

# Reading the data
data <- as.data.frame(read.csv("../../../Data/bigfish.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE))
attach(data)
alpha <- 0.05  # Significance level

# OLS results

# Demand
formula <- qty ~ price + day1 + day2 + day3 + day4 + rainy + cold
ols.fit <- lm(formula)
ols.coef.d <- ols.fit$coef[2]
ols.se.d <- coeftest(ols.fit, vcov = hc)[2, 2]
ols.lci.d <- ols.coef.d + qnorm(alpha / 2) * ols.se.d
ols.uci.d <- ols.coef.d + qnorm(1 - alpha / 2) * ols.se.d

# Supply
formula <- qty ~ price + stormy + mixed
ols.fit <- lm(formula)
ols.coef.s <- ols.fit$coef[2]
ols.se.s <- coeftest(ols.fit, vcov = hc)[2, 2]
ols.lci.s <- ols.coef.s + qnorm(alpha / 2) * ols.se.s
ols.uci.s <- ols.coef.s + qnorm(1 - alpha / 2) * ols.se.s

# First stage

# Demand
formula <- price ~ stormy + mixed + day1 + day2 + day3 + day4 + rainy + cold
fs.fit <- lm(formula)
fs.Fstat <- waldtest(fs.fit, . ~ . - stormy - mixed, vcov = hc, test = "F")$F[2]
print(paste('F-stat Demand: ', fs.Fstat))

# Supply
fs.Fstat <- waldtest(fs.fit, . ~ . - day1 - day2 - day3 - day4 - rainy - cold, vcov = hc, test = "F")$F[2]
print(paste('F-stat Supply: ', fs.Fstat))

# TSLS results

# Demand
formula <- qty ~ price + day1 + day2 + day3 + day4 + rainy + cold | stormy + mixed + day1 + day2 + day3 + day4 + rainy + cold
tsls.fit <- ivreg(formula)
tsls.coef.d <- tsls.fit$coef[2]
tsls.se.d <- coeftest(tsls.fit, vcov = hc)[2, 2]
tsls.lci.d <- tsls.coef.d + qnorm(alpha / 2) * tsls.se.d
tsls.uci.d <- tsls.coef.d + qnorm(1 - alpha / 2) * tsls.se.d

# Supply
formula <- qty ~ price + stormy + mixed | day1 + day2 + day3 + day4 + rainy + cold + stormy + mixed
tsls.fit <- ivreg(formula)
tsls.coef.s <- tsls.fit$coef[2]
tsls.se.s <- coeftest(tsls.fit, vcov = hc)[2, 2]
tsls.lci.s <- tsls.coef.s + qnorm(alpha / 2) * tsls.se.s
tsls.uci.s <- tsls.coef.s + qnorm(1 - alpha / 2) * tsls.se.s

# GMM results: equation by equation (LI-3SLS)

# Demand
instr <- cbind(stormy, mixed, day1, day2, day3, day4, rainy, cold)
formula <- qty ~ price + day1 + day2 + day3 + day4 + rainy + cold
gmm.fit.d <- gmm(formula, x = instr, vcov = "iid")
gmm.coef.d <- summary(gmm.fit.d)$coefficients[2, 1]
gmm.se.d <- summary(gmm.fit.d)$coefficients[2, 2]
gmm.lci.d <- gmm.coef.d + qnorm(alpha / 2) * gmm.se.d
gmm.uci.d <- gmm.coef.d + qnorm(1 - alpha / 2) * gmm.se.d
print(summary(gmm.fit.d)$stest)

# Supply
instr <- cbind(stormy, mixed, day1, day2, day3, day4, rainy, cold)
formula <- qty ~ price + stormy + mixed
gmm.fit.s <- gmm(formula, x = instr, vcov = "iid")
gmm.coef.s <- summary(gmm.fit.s)$coefficients[2, 1]
gmm.se.s <- summary(gmm.fit.s)$coefficients[2, 2]
gmm.lci.s <- gmm.coef.s + qnorm(alpha / 2) * gmm.se.s
gmm.uci.s <- gmm.coef.s + qnorm(1 - alpha / 2) * gmm.se.s

# 3SLS: systemfit (FI-3SLS)
formula.d <- qty ~ price + day1 + day2 + day3 + day4 + rainy + cold
formula.s <- qty ~ price + stormy + mixed
formula <- list(demand = formula.d, supply = formula.s)
s3sls.fit <- systemfit(formula, method = "3SLS", inst = ~ day1 + day2 + day3 + day4 + rainy + cold + stormy + mixed)

# Demand
s3sls.coef.d <- s3sls.fit$coefficients[2]
s3sls.se.d <- summary(s3sls.fit)$coefficients[2, 2]
s3sls.lci.d <- s3sls.coef.d + qnorm(alpha / 2) * s3sls.se.d
s3sls.uci.d <- s3sls.coef.d + qnorm(1 - alpha / 2) * s3sls.se.d

# Supply
s3sls.coef.s <- s3sls.fit$coefficients[10]
s3sls.se.s <- summary(s3sls.fit)$coefficients[10, 2]
s3sls.lci.s <- s3sls.coef.s + qnorm(alpha / 2) * s3sls.se.s
s3sls.uci.s <- s3sls.coef.s + qnorm(1 - alpha / 2) * s3sls.se.s

# Bootstrap for GMM (LI-3SLS)

# Bootstrap function for demand
gmm_demand_boot <- function(data, indices) {
  d <- data[indices, ]
  instr <- with(d, cbind(stormy, mixed, day1, day2, day3, day4, rainy, cold))
  formula <- qty ~ price + day1 + day2 + day3 + day4 + rainy + cold
  gmm.fit <- gmm(formula, x = instr, vcov = "iid")
  return(gmm.fit$coefficients["price"])
}

# Bootstrap function for supply
gmm_supply_boot <- function(data, indices) {
  d <- data[indices, ]
  instr <- with(d, cbind(stormy, mixed, day1, day2, day3, day4, rainy, cold))
  formula <- qty ~ price + stormy + mixed
  gmm.fit <- gmm(formula, x = instr, vcov = "iid")
  return(gmm.fit$coefficients["price"])
}

# Run bootstrap with 400 replications
set.seed(123)  # For reproducibility
boot_demand <- boot(data, gmm_demand_boot, R = 400)
boot_supply <- boot(data, gmm_supply_boot, R = 400)

# Bootstrap standard errors
boot_se_demand <- sd(boot_demand$t)
boot_se_supply <- sd(boot_supply$t)

# Bootstrap 95% CIs (percentile method)
boot_ci_demand <- quantile(boot_demand$t, c(0.025, 0.975))
boot_ci_supply <- quantile(boot_supply$t, c(0.025, 0.975))

# Output bootstrap results
cat("\nBootstrap Results for Demand Elasticity (LI-3SLS):\n")
cat("Standard Error:", boot_se_demand, "\n")
cat("95% CI:", boot_ci_demand[1], "to", boot_ci_demand[2], "\n")

cat("\nBootstrap Results for Supply Elasticity (LI-3SLS):\n")
cat("Standard Error:", boot_se_supply, "\n")
cat("95% CI:", boot_ci_supply[1], "to", boot_ci_supply[2], "\n")

# Existing result tables
options(digits = 2)

# Demand table
table_demand <- matrix(0, ncol = 4, nrow = 4, dimnames = list(c('LI-OLS', 'LI-TSLS', 'LI-3SLS', 'FI-3SLS'), c('Est.', 'Std. Error', '95% LCI', '95% UCI')))
table_demand[1, ] <- c(ols.coef.d, ols.se.d, ols.lci.d, ols.uci.d)
table_demand[2, ] <- c(tsls.coef.d, tsls.se.d, tsls.lci.d, tsls.uci.d)
table_demand[3, ] <- c(gmm.coef.d, gmm.se.d, gmm.lci.d, gmm.uci.d)
table_demand[4, ] <- c(s3sls.coef.d, s3sls.se.d, s3sls.lci.d, s3sls.uci.d)
print("Demand Elasticity Estimates:")
print(table_demand)

# Supply table
table_supply <- matrix(0, ncol = 4, nrow = 4, dimnames = list(c('LI-OLS', 'LI-TSLS', 'LI-3SLS', 'FI-3SLS'), c('Est.', 'Std. Error', '95% LCI', '95% UCI')))
table_supply[1, ] <- c(ols.coef.s, ols.se.s, ols.lci.s, ols.uci.s)
table_supply[2, ] <- c(tsls.coef.s, tsls.se.s, tsls.lci.s, tsls.uci.s)
table_supply[3, ] <- c(gmm.coef.s, gmm.se.s, gmm.lci.s, gmm.uci.s)
table_supply[4, ] <- c(s3sls.coef.s, s3sls.se.s, s3sls.lci.s, s3sls.uci.s)
print("Supply Elasticity Estimates:")
print(table_supply)