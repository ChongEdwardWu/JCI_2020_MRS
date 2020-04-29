########################
### Install packages ###
########################

#### Install necessary R packages ###
## if you have not installed the necessary R packages, run the following codes ##
# if (!requireNamespace("rms", quietly = TRUE))
#   install.packages("rms")
# if (!requireNamespace("glmnet", quietly = TRUE))
#   install.packages("glmnet")


##################
### Load  data ###
##################


## Set your working directory (or use Control + Shift + H to choose it) ##
## The working directory should contain the RData file "MRS construction with lasso.RData" ##
setwd("...")
load("MRS construction with lasso.RData")
#Training <- read.csv("Training.csv")
View(Training)

############################
###-----Lasso Cox-TTR----###
############################
require(rms)
require(glmnet)

# prepare data
max.time <- 96
Training$TTR.outcome[Training$TTR.time>max.time ] <- 0
Training$TTR.time[Training$TTR.time>max.time ] <- max.time
TTR <- Surv(Training$TTR.time, Training$TTR.outcome!=0)

##identify preditors
x <- model.matrix( ~ CD33T + CD11bT + CD15T + CD68T + CD204T + CD206T +
                     CD169T + CD163T + S100T +
                     CD33N + CD11bN +  CD15N + CD68N + CD204N + CD206N +
                     CD169N + CD163N + S100N, data=Training)



# L1-penalized (Lasso) Cox regression model to select prognostic myeloid features using the 1-standard error (SE) criteria
TTR.fit1 <- glmnet(x, TTR, family = "cox", alpha=1)
set.seed(349)
cv.TTR.fit1 <- cv.glmnet(x, TTR, family = "cox", alpha=1, nfolds=10)

# plot the results
ThisPal <- c("#E78AC3","#FC8D62","#A6D854","#BEBADA",
             "#E5C494","#CCEBC5","#8DA0CB","#E5C494","#E78AC3",
             "#BEBADA","#FDB462","#D9D9D9","#A6D854","#E5C494",
             "#8DA0CB", "#BEBADA","#D9D9D9","#FDB462")
plot(TTR.fit1, xvar="lambda", label=TRUE,lwd=1,col=ThisPal)
abline(v=log(cv.TTR.fit1$lambda.1se),lty=3, col="red",lwd=1)
plot(cv.TTR.fit1)
log(cv.TTR.fit1$lambda.1se)

# print co-efficients
coef <- as.matrix(coef(TTR.fit1))
write.csv(coef,"Lasso-coef.csv")
write.csv(log(TTR.fit1$lambda),"Lasso-coef.csv")

# clean environment
graphics.off()
rm(list = ls())
