install.packages("xtable")
install.packages("hdm")
install.packages("randomForest")
install.packages("glmnet")
install.packages("sandwich")

library(xtable)
library(randomForest)
library(hdm)
library(glmnet)
library(sandwich)

set.seed(1)

# file <- "https://raw.githubusercontent.com/CausalAIBook/MetricsMLNotebooks/main/data/GrowthData.csv"
# data <- read.csv(file)
# data <- subset(data, select = -1) # get rid of index column


# Penn<- as.data.frame(read.table("../../Data/penn_jae.dat", header=T ));
# Penn <- na.omit(Penn)
# 
# Penn$age <- as.numeric(factor(2 - Penn$agelt35 + Penn$agegt54, levels = c(1, 2, 3)))
# 
# # Create unemployment level factor (numeric levels: 1, 2, 3)
# Penn$unemp <- as.numeric(factor(Penn$lusd + 2*Penn$muld + 3*Penn$husd, levels = c(1, 2, 3)))
# 
# # Create race factor (numeric levels: 1, 2, 3, 4)
# Penn$race <- as.numeric(factor(1 + Penn$black + 2*Penn$hispanic + 3*Penn$othrace, levels = c(1, 2, 3, 4)))
# 
# # Create sex factor (numeric levels: 0, 1)
# Penn$sex <- as.numeric(factor(Penn$female, levels = c(0, 1))) - 1  # Subtract 1 to match levels 0, 1
# 
# # Create quarter factor (numeric levels: 1, 2, 3, 4, 5, 6)
# Penn$quarter <- as.numeric(factor(Penn$q1 + 2*Penn$q2 + 3*Penn$q3 + 4*Penn$q4 + 5*Penn$q5 + 6*Penn$q6, levels = c(1, 2, 3, 4, 5, 6)))
# 
# # Create log-duration variable (unchanged)
# Penn$lndur <- log(Penn$inuidur1)
# 
# 
# # Merge treatment groups 4 & 6; 
# Penn$tg[Penn$tg == 6]   <- 4;  # Assign observations in treatment group 6 to treatment group 4 following Bilias (2000);
# 
# 
# # create treatment indicator;
# Penn$treated <- 1*(Penn$tg == 4);
# 
# 
# 
# Penn_filtered <- Penn[Penn$tg %in% c(0, 4), ]
# attach(Penn);




# Create y, d, x matrices
# y: log(inuidur1)
# y <- as.matrix(Penn_filtered$lndur)
# 
# # d: factor(tg)
# d <- as.matrix(as.numeric(as.factor(Penn_filtered$tg)))
# 
# # x: controls (sex, race, factor(dep), quarter, age, durable, nondurable, unemp)
# # Assuming 'dep', 'durable', 'nondurable' are already in the dataset
# x <- as.matrix(Penn_filtered[, c("sex", "race", "dep", "quarter", "age", "durable", "nondurable")])

# y <- as.matrix(data[, 1]) # outcome: growth rate
# d <- as.matrix(data[, 3]) # treatment: initial wealth
# x <- as.matrix(data[, -c(1, 2, 3)]) # controls: country characteristics

y <- as.matrix(my_data$GDP) # outcome: growth rate
d <- as.matrix(my_data$Exprop) # treatment: initial wealth
x <- as.matrix(my_data[, c("Latitude", "Location")]) # controls: country characteristics



dml2_for_plm <- function(x, d, y, dreg, yreg, nfold = 2) {
  nobs <- nrow(x) # number of observations
  foldid <- rep.int(1:nfold, times = ceiling(nobs / nfold))[sample.int(nobs)] # define folds indices
  I <- split(1:nobs, foldid) # split observation indices into folds
  ytil <- dtil <- rep(NA, nobs)
  cat("fold: ")
  for (b in seq_along(I)) {
    dfit <- dreg(x[-I[[b]], ], d[-I[[b]]]) # take a fold out
    yfit <- yreg(x[-I[[b]], ], y[-I[[b]]]) # take a foldt out
    dhat <- predict(dfit, x[I[[b]], ], type = "response") # predict the left-out fold
    yhat <- predict(yfit, x[I[[b]], ], type = "response") # predict the left-out fold
    dtil[I[[b]]] <- (d[I[[b]]] - dhat) # record residual for the left-out fold
    ytil[I[[b]]] <- (y[I[[b]]] - yhat) # record residial for the left-out fold
    cat(b, " ")
  }
  rfit <- lm(ytil ~ dtil) # estimate the main parameter by regressing one residual on the other
  coef.est <- coef(rfit)[2] # extract coefficient
  se <- sqrt(vcovHC(rfit)[2, 2]) # record robust standard error
  cat(sprintf("\ncoef (se) = %g (%g)\n", coef.est, se)) # printing output
  return(list(coef.est = coef.est, se = se, dtil = dtil, ytil = ytil)) # save output and residuals
}





# DML with OLS
cat(sprintf("\nDML with OLS w/o feature selection \n"))
dreg <- function(x, d) {
  glmnet(x, d, lambda = 0)
} # ML method= OLS using glmnet; using lm gives bugs
yreg <- function(x, y) {
  glmnet(x, y, lambda = 0)
} # ML method = OLS
dml2_ols <- dml2_for_plm(x, d, y, dreg, yreg, nfold = 10)


# DML with Lasso:
cat(sprintf("\nDML with Lasso \n"))
dreg <- function(x, d) {
  rlasso(x, d, post = FALSE)
} # ML method= lasso from hdm
yreg <- function(x, y) {
  rlasso(x, y, post = FALSE)
} # ML method = lasso from hdm
dml2_lasso <- dml2_for_plm(x, d, y, dreg, yreg, nfold = 10)


# DML with Random Forest:
cat(sprintf("\nDML with Random Forest \n"))
dreg <- function(x, d) {
  randomForest(x, d)
} # ML method=Forest
yreg <- function(x, y) {
  randomForest(x, y)
} # ML method=Forest
dml2_rf <- dml2_for_plm(x, d, y, dreg, yreg, nfold = 10)

# DML MIX:
cat(sprintf("\nDML with Lasso for D and Random Forest for Y \n"))
dreg <- function(x, d) {
  rlasso(x, d, post = FALSE)
} # ML method=Forest
yreg <- function(x, y) {
  randomForest(x, y)
} # ML method=Forest
dml2_mix <- dml2_for_plm(x, d, y, dreg, yreg, nfold = 10)


pr_res_d <- c(mean((dml2_ols$dtil)^2), mean((dml2_lasso$dtil)^2), mean((dml2_rf$dtil)^2), mean((dml2_mix$dtil)^2))
pr_res_y <- c(mean((dml2_ols$ytil)^2), mean((dml2_lasso$ytil)^2), mean((dml2_rf$ytil)^2), mean((dml2_mix$ytil)^2))
pr_res <- rbind(sqrt(pr_res_d), sqrt(pr_res_y))
rownames(pr_res) <- c("RMSE D", "RMSE Y")
colnames(pr_res) <- c("OLS", "Lasso", "RF", "Mix")

table <- matrix(0, 4, 4)

# Point Estimate
table[1, 1] <- as.numeric(dml2_ols$coef.est)
table[2, 1] <- as.numeric(dml2_lasso$coef.est)
table[3, 1] <- as.numeric(dml2_rf$coef.est)
table[4, 1] <- as.numeric(dml2_mix$coef.est)

# SE
table[1, 2] <- as.numeric(dml2_ols$se)
table[2, 2] <- as.numeric(dml2_lasso$se)
table[3, 2] <- as.numeric(dml2_rf$se)
table[4, 2] <- as.numeric(dml2_mix$se)

# RMSE Y
table[1, 3] <- as.numeric(pr_res[2, 1])
table[2, 3] <- as.numeric(pr_res[2, 2])
table[3, 3] <- as.numeric(pr_res[2, 3])
table[4, 3] <- as.numeric(pr_res[2, 4])

# RMSE D
table[1, 4] <- as.numeric(pr_res[1, 1])
table[2, 4] <- as.numeric(pr_res[1, 2])
table[3, 4] <- as.numeric(pr_res[1, 3])
table[4, 4] <- as.numeric(pr_res[1, 4])



# print results
colnames(table) <- c("Estimate", "Standard Error", "RMSE Y", "RMSE D")
rownames(table) <- c("OLS", "Lasso", "RF", "RF/Lasso Mix")
table

print(table, digit = 3)

tab <- xtable(table, digits = 3)
print(tab, type = "latex")