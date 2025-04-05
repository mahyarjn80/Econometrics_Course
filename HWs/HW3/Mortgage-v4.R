
# Description of the data: the sample selection and variable contruction follow

# Sample selection: single-residences (excluding data on multifamily homes), 2,380 observation

# The variables in the data set include:

# deny:  = 1 if mortgage application denied
# p_irat: montly debt to income ratio
# black: = 1 if applicant is black
# hse_inc: montly housing expenses to income ratio
# loan_val: loan to assessed property value ratio (not used in analysis)
# ccred: consumer credit score 
#        = 1 if no "slow" payments or delinquencies
#        = 2 if one or two "slow" payments or delinquencies
#        = 3 if more than two "slow" payments or delinquencies
#        = 4 if insufficient credit history for determination
#        = 5 if delinquent credit history with payments 60 days overdue
#        = 6 if delinquent credit history with payments 90 days overdue
# mcred: mortgage credit score
#        = 1 if no late mortgage payments
#        = 2 if no mortgage payment history
#        = 3 if one or two late mortgage payments
#        = 4 if more than two late mortgage payments
# pubrec: = 1 if any public record of credit problems (bankruptcy, charge-offs, collection actions)
# denpmi: = 1 if applicant applied for mortgage insurance and was denied
# selfemp: = 1 if self-employed
# single: = 1 if single
# hischl: = 1 if high school graduated
# probunmp: 1989 Massachusetts unemployment rate in the applicant's industry (not used in analysis)
# condo: = 1 if unit is a condominium (not used in analysis)
# ltv_med: = 1 if medium loan to property value ratio [.80, .95]
# ltv_high: = 1 if high loan to property value ratio > .95


rm(list=ls(all=TRUE)) ## eliminating everything in  memory



filepath<-"../../Data/mortgage.dta"
###  read-in Mortgage data

library(foreign);
library(xtable);
# Load the haven library
library(haven)
library(SortedEffects); # Package for sorted effects, includes mortgage data


data <- read_dta(filepath)

# Attach the data to your workspace
attach(data)

########## Descriptive statistics

options(digits=3);
dstats <- cbind(sapply(mortgage, mean), apply(mortgage[mortgage$black==1,], 2, mean), apply(mortgage[mortgage$black==0,], 2, mean));
xtable(dstats);

# Ols regressions

fmla1 =  deny ~ black + p_irat + hse_inc + ccred + mcred  + pubrec + ltv_med + ltv_high + denpmi + selfemp + single + hischl
fit.1=lm(fmla1, data=mortgage)

fmla2 =  deny ~ black
fit.2= lm(fmla2, data=mortgage)

xtable(summary(fit.2, digits=3))
xtable(summary(fit.1, digits=3))


########## here we will split the data in two parts, data 1 and data 2

set.seed(1)
n1<- dim(mortgage)[1]
validate<- sample(1:n1, floor(n1/3), replace=F)

data1<- mortgage[-c(validate),]

### data 1 is the main part

data2<- mortgage[c(validate),]

### data 2 is the validation part, which will be used to select the link function
### you can also use to choose better predictive models


### the following command attaches data1, making its columns directly accessible by their names

attach(data1)

################ Part 1: Compare some basic models #########################################################
#basic linear probability model, linear index

fmla1 =  deny ~ black + p_irat + hse_inc + ccred + mcred  + pubrec + ltv_med + ltv_high + denpmi + selfemp + single + hischl
fit.1=lm(fmla1)

#### basic models with logistic, probit, and cauchit links

fit.lgt.1=glm(fmla1, family=binomial(link="logit"))
fit.prt.1=glm(fmla1, family=binomial(link="probit"))
fit.cat.1=glm(fmla1, family=binomial(link="cauchit"))

# pdf("L6/Results/LogitvsLinear1_mort.pdf", pointsize=15,width=6,height=6)
plot(predict(fit.lgt.1, type="response"), predict(fit.1, type="response"),
     xlim=c(0,1), ylim=c(-0.1,1.20), xlab="logit prediction", ylab="linear & cauchit prediction",
     main=" ", type="p", pch=1, lty=2 ,lwd=1, col="red")
lines(predict(fit.lgt.1, type="response"),predict(fit.cat.1, type="response"),
     type="p", pch=1, col="black")
abline(0, 1)
abline(h=1)
abline(h=0)
legend(x=-0.05, y=1.25, bty = "n",
        legend=c("Logit vs Linear","Logit vs Cauchy"),
                 col=c("red", "black"), lwd=c(1, 1), lty=c(0, 0), pch=c(1,1), cex = 0.7)




# pdf("L6/Results/LogitvsProbit1_mort.pdf", pointsize=15,width=6,height=6)
plot(predict(fit.lgt.1, type="response"),predict(fit.prt.1, type="response"),
     xlim=c(0,1), ylim=c(0,1), xlab="logit prediction", ylab="probit prediction",
     main=" ",col=4)
abline(0, 1)




# compare out-of-sample prediction scores

table.P<- matrix(0, ncol=4, nrow=1)
table.P[1,1]<- sqrt(mean((data2[,"deny"]- predict(fit.lgt.1, data2, type="response"))^2))
table.P[1,2]<- sqrt(mean((data2[,"deny"]- predict(fit.prt.1, data2, type="response"))^2))
table.P[1,3]<- sqrt(mean((data2[,"deny"]- predict(fit.cat.1, data2, type="response"))^2))
table.P[1,4]<- sqrt(mean((data2[,"deny"]- predict(fit.1, data2, type="response"))^2))

colnames(table.P)<- c("Logit", "Probit", "Cauchit", "Linear")
rownames(table.P)<- c("Mean Square Prediction Error")


xtable(table.P, digits=4, align=c(rep("c", 5)))

detach(data1);

# logistic seems to work weakly better than others in terms of predicting
# 
############## Part 2. Predicted Effects of Black ########################################################

us    <- c(2:98)/100; # percentiles of interest for SPE
alpha  <- .1; #Significance level
beta   <- .05; #Tail quantile index for classication analysis
R      <- 500; #number of bootstrap replications

# Estimates and confidence bands for APE and SPE for population using Logit;

spe.fit <- spe(fm = fmla1, data = mortgage, method = c("probit"), var = "black", us = us, alpha = alpha, b = R, parallel = TRUE, ncores = 24);
xtable(summary(spe.fit, result = "average"), digits=3);
plot(x = spe.fit, ylim = c(0, 0.25), ylab = "Change in Probability",main = " ",sub = "probit Model")

# Classification analysis for population;

# Specify variables of interest;
t <- c("deny", "black", "p_irat", "hse_inc", "ccred", "mcred", "pubrec", "denpmi", "ltv_med", "ltv_high","selfemp", "single", "hischl");
ca.fit <- ca(fm = fmla1, data = mortgage, var = "black", method = "logit", u = beta, alpha = alpha, t = t, b = R, parallel = TRUE, ncores = 24);
xtable(summary(ca.fit));


# Estimates and confidence bands for APE and SPE for black subpopulation using Logit;
spe.black.fit <- spe(fm = fmla1, data = mortgage, method = c("probit"), var = "black", us = us, subgroup = (mortgage[,"black"]==1), alpha = alpha, b = R, parallel = TRUE, ncores = 24);
xtable(summary(spe.black.fit, result = "average"), digits=3);
plot(x = spe.black.fit, ylim = c(0, 0.25), ylab = "Change in Probability",main = " ",sub = "probit Model")

# 
