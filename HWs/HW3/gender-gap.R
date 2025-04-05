# Clear the workspace
rm(list = ls());

# Load required libraries
library(xtable);    # For generating well-formatted tables
library(boot);      # For bootstrap resampling
library(Hmisc);     # For weighted statistical functions
library(SortedEffects); # For sorted effects analysis with CPS data

# Suppress warnings to keep output clean
options(warn = -1);

# Load the CPS wage dataset
data(wage2015);

# Define indicator variables for categorical covariates

# Occupation indicators
wage2015$job_mgr <- ifelse(wage2015$occ == "manager", 1, 0);
wage2015$job_serv <- ifelse(wage2015$occ == "service", 1, 0);
wage2015$job_sales <- ifelse(wage2015$occ == "sales", 1, 0);
wage2015$job_const <- ifelse(wage2015$occ == "construction", 1, 0);
wage2015$job_prod <- ifelse(wage2015$occ == "production", 1, 0);

# Industry indicators
wage2015$sec_mining <- ifelse(wage2015$ind == "minery", 1, 0);
wage2015$sec_build <- ifelse(wage2015$ind == "construction", 1, 0);
wage2015$sec_mfg <- ifelse(wage2015$ind == "manufacture", 1, 0);
wage2015$sec_retail <- ifelse(wage2015$ind == "retail", 1, 0);
wage2015$sec_trans <- ifelse(wage2015$ind == "transport", 1, 0);
wage2015$sec_info <- ifelse(wage2015$ind == "information", 1, 0);
wage2015$sec_fin <- ifelse(wage2015$ind == "finance", 1, 0);
wage2015$sec_prof <- ifelse(wage2015$ind == "professional", 1, 0);
wage2015$sec_edu <- ifelse(wage2015$ind == "education", 1, 0);
wage2015$sec_leis <- ifelse(wage2015$ind == "leisure", 1, 0);
wage2015$sec_serv <- ifelse(wage2015$ind == "services", 1, 0);
wage2015$sec_pub <- ifelse(wage2015$ind == "public", 1, 0);

# Education indicators
wage2015$edu_lesshs <- ifelse(wage2015$educ == "lhs", 1, 0);
wage2015$edu_hschool <- ifelse(wage2015$educ == "hsg", 1, 0);
wage2015$edu_somecol <- ifelse(wage2015$educ == "sc", 1, 0);
wage2015$edu_collg <- ifelse(wage2015$educ == "cg", 1, 0);
wage2015$edu_advdeg <- ifelse(wage2015$educ == "ad", 1, 0);

# Marital status indicators
wage2015$mar_never <- ifelse(wage2015$ms == "nevermarried", 1, 0);
wage2015$mar_married <- ifelse(wage2015$ms == "married", 1, 0);
wage2015$mar_separ <- ifelse(wage2015$ms == "separated", 1, 0);
wage2015$mar_div <- ifelse(wage2015$ms == "divorced", 1, 0);
wage2015$mar_wid <- ifelse(wage2015$ms == "widowed", 1, 0);

# Region indicators
wage2015$reg_northeast <- ifelse(wage2015$region == "ne", 1, 0);
wage2015$reg_midwest <- ifelse(wage2015$region == "mw", 1, 0);
wage2015$reg_south <- ifelse(wage2015$region == "so", 1, 0);
wage2015$reg_west <- ifelse(wage2015$region == "we", 1, 0);

# Attach the dataset for easier variable access
attach(wage2015);

# Define parameters for analysis
percentiles <- seq(0.02, 0.98, by = 0.01); # Percentile range for SPE
conf_level <- 0.10; # 90% confidence level
tail_quant <- 0.05; # Tail quantile for classification
boot_reps <- 500; # Number of bootstrap replications

# List of variables for descriptive statistics
vars_list <- c("lnw", "female", "mar_married", "mar_wid", "mar_separ", "mar_div", "mar_never", 
               "edu_lesshs", "edu_hschool", "edu_somecol", "edu_collg", "edu_advdeg", 
               "reg_northeast", "reg_midwest", "reg_south", "reg_west", "exp1", 
               "job_mgr", "job_serv", "job_sales", "job_const", "job_prod", 
               "sec_mining", "sec_build", "sec_mfg", "sec_retail", "sec_trans", 
               "sec_info", "sec_fin", "sec_prof", "sec_edu", "sec_leis", "sec_serv", "sec_pub");

# Compute weighted descriptive statistics
options(digits = 2);
desc_stats <- data.frame(
  All = sapply(wage2015[, vars_list], function(x) weighted.mean(x, weight)),
  Women = sapply(wage2015[wage2015$female == 1, vars_list], function(x) weighted.mean(x, weight[wage2015$female == 1])),
  Men = sapply(wage2015[wage2015$female == 0, vars_list], function(x) weighted.mean(x, weight[wage2015$female == 0]))
);
print(xtable(desc_stats, digits = 2), type = "latex");

# Gender Wage Gap Analysis for All Women
model_formula <- lnw ~ female * (mar_married + mar_never + mar_separ + mar_div + mar_wid + 
                                   (exp1 + I(exp1^2/100) + I(exp1^3/1000) + I(exp1^4/10000)) * 
                                   (edu_lesshs + edu_hschool + edu_somecol + edu_collg + edu_advdeg) + 
                                   job_mgr + job_serv + job_sales + job_const + job_prod + 
                                   sec_mining + sec_build + sec_mfg + sec_retail + sec_trans + 
                                   sec_info + sec_fin + sec_prof + sec_edu + sec_leis + sec_serv + sec_pub + 
                                   reg_northeast + reg_midwest + reg_south + reg_west);

set.seed(123);
spe_all_fit <- spe(formula = model_formula, data = wage2015, sample_weights = weight, 
                   treatment = "female", method = "ols", subgroup = (female == 1), 
                   percentiles = percentiles, conf_level = conf_level, bc = FALSE, 
                   B = boot_reps, parallel = TRUE, n_cores = 24, boot_type = "nonpar");
print(xtable(summary(spe_all_fit, result = "average"), digits = 3), type = "latex");
plot(spe_all_fit, main = "", sub = "", xlab = "Percentile Index", ylab = "Gender Wage Gap", 
     ylim = c(-0.45, 0.1));

# Classification Analysis for Women
class_vars <- c("lnw", "mar_married", "mar_never", "mar_separ", "mar_div", "mar_wid", 
                "edu_lesshs", "edu_hschool", "edu_somecol", "edu_collg", "edu_advdeg", 
                "reg_northeast", "reg_midwest", "reg_south", "reg_west", "exp1", 
                "job_mgr", "job_serv", "job_sales", "job_const", "job_prod", 
                "sec_mining", "sec_build", "sec_mfg", "sec_retail", "sec_trans", 
                "sec_info", "sec_fin", "sec_prof", "sec_edu", "sec_leis", "sec_serv", "sec_pub");
class_fit <- ca(formula = model_formula, data = wage2015, sample_weights = weight, 
                treatment = "female", method = "ols", variables = class_vars, 
                classification = "both", quantile = tail_quant, conf_level = conf_level, 
                B = boot_reps, subgroup = (female == 1), boot_type = "nonpar", 
                parallel = TRUE, n_cores = 24, bc = FALSE);
print(xtable(summary(class_fit)), type = "latex");

# SPE for Women by Marital Status
married_group <- (female == 1 & mar_married == 1);
nevermarried_group <- (female == 1 & mar_never == 1);

set.seed(123);
spe_married_fit <- spe(formula = model_formula, data = wage2015, sample_weights = weight, 
                       treatment = "female", method = "ols", subgroup = married_group, 
                       percentiles = percentiles, conf_level = conf_level, bc = FALSE, 
                       B = boot_reps, parallel = TRUE, n_cores = 24, boot_type = "nonpar");
spe_nevermarried_fit <- spe(formula = model_formula, data = wage2015, sample_weights = weight, 
                            treatment = "female", method = "ols", subgroup = nevermarried_group, 
                            percentiles = percentiles, conf_level = conf_level, bc = FALSE, 
                            B = boot_reps, parallel = TRUE, n_cores = 24, boot_type = "nonpar");

par(mfrow = c(1, 2));
plot(spe_nevermarried_fit, main = "Never Married", sub = "", xlab = "Percentile Index", 
     ylab = "Gender Wage Gap", ylim = c(-0.45, 0.2));
plot(spe_married_fit, main = "Married", sub = "", xlab = "Percentile Index", 
     ylab = "Gender Wage Gap", ylim = c(-0.45, 0.2));

# SPE for Women by Experience Group
exp_quantiles <- wtd.quantile(exp1[female == 1], weights = weight[female == 1], probs = c(0.25, 0.75));
low_exp_group <- (female == 1 & exp1 <= exp_quantiles[1]);
mid_exp_group <- (female == 1 & exp1 > exp_quantiles[1] & exp1 < exp_quantiles[2]);
high_exp_group <- (female == 1 & exp1 >= exp_quantiles[2]);

set.seed(123);
spe_low_fit <- spe(formula = model_formula, data = wage2015, sample_weights = weight, 
                   treatment = "female", method = "ols", subgroup = low_exp_group, 
                   percentiles = percentiles, conf_level = conf_level, bc = FALSE, 
                   B = boot_reps, parallel = TRUE, n_cores = 24, boot_type = "nonpar");
spe_mid_fit <- spe(formula = model_formula, data = wage2015, sample_weights = weight, 
                   treatment = "female", method = "ols", subgroup = mid_exp_group, 
                   percentiles = percentiles, conf_level = conf_level, bc = FALSE, 
                   B = boot_reps, parallel = TRUE, n_cores = 24, boot_type = "nonpar");
spe_high_fit <- spe(formula = model_formula, data = wage2015, sample_weights = weight, 
                    treatment = "female", method = "ols", subgroup = high_exp_group, 
                    percentiles = percentiles, conf_level = conf_level, bc = FALSE, 
                    B = boot_reps, parallel = TRUE, n_cores = 24, boot_type = "nonpar");

par(mfrow = c(1, 3));
plot(spe_low_fit, main = "Low Experience", sub = "", xlab = "Percentile Index", 
     ylab = "Gender Wage Gap", ylim = c(-0.45, 0.2));
plot(spe_mid_fit, main = "Mid Experience", sub = "", xlab = "Percentile Index", 
     ylab = "Gender Wage Gap", ylim = c(-0.45, 0.2));
plot(spe_high_fit, main = "High Experience", sub = "", xlab = "Percentile Index", 
     ylab = "Gender Wage Gap", ylim = c(-0.45, 0.2));