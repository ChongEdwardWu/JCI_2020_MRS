########################
### Install packages ###
########################

#### Install necessary R packages ###
## if you have not installed the necessary R packages, run the following codes ##
# if (!requireNamespace("rms", quietly = TRUE))
#   install.packages("rms")
# if (!requireNamespace("survivalROC", quietly = TRUE))
#   install.packages("survivalROC")
# if (!requireNamespace("boot", quietly = TRUE))
#   install.packages("boot")


#############################
### load data and prepare ###
#############################

require(survivalROC)
require(rms)
require(boot)

## Set your working directory (or use Control + Shift + H to choose it) ##
## The working directory should contain the RData file "ROC analysis of MRS in cohorts.RData" 
setwd("...")


#set the coloar palette
MyPal <- c(brewer.pal(12,"Set3")[c(3,6,8,10,11)],"#A6CEE3",brewer.pal(8,"Set2")[-6])


# load data from a csv table, or alternatively, from a RData file
# WholeCohort <- read.csv("CombinedWholeCohort.csv",header=TRUE, as.is=TRUE, na.strings=c("","NA"))
load("ROC analysis of MRS in cohorts.RData")
palette(MyPal) 

# preprocess data
max.time <- 96
units(WholeCohort$TTR.time) <- 'Month'
units(WholeCohort$OS.time) <- 'Month'

WholeCohort$TTR.outcome[WholeCohort$TTR.time>max.time ] <- 0
WholeCohort$TTR.time[WholeCohort$TTR.time>max.time ] <- max.time
WholeCohort$TTR <- Surv(WholeCohort$TTR.time, WholeCohort$TTR.outcome!=0)

WholeCohort$OS.outcome[WholeCohort$OS.time>max.time ] <- 0
WholeCohort$OS.time[WholeCohort$OS.time>max.time ] <- max.time
WholeCohort$OS <- Surv(WholeCohort$OS.time, WholeCohort$OS.outcome!=0)


################################
### MRS in different cohorts ###
################################

### Define cohorts
# Training set
Training <- WholeCohort[WholeCohort$Cohort==0,]
# Test set
Test <- WholeCohort[WholeCohort$Cohort==1,]
# Primary cohort
Primary <- WholeCohort[(WholeCohort$Cohort==0) | (WholeCohort$Cohort==1),]
# Internal Cohort 
Internal <- WholeCohort[(WholeCohort$Cohort==2),]
# External Cohort 
External <- WholeCohort[(WholeCohort$Cohort==3) | (WholeCohort$Cohort==4),]
# External-THZP Cohort
THZP <- WholeCohort[(WholeCohort$Cohort==3),]
# External-THZP Cohort
FUZH <- WholeCohort[(WholeCohort$Cohort==4),]

### define a function to calculate the AUC and the 95% CI
boot.TTR.AUC <- function(data,indices,cutoff){ #define a function for bootstrap
  d <- data[indices,] # allows boot to select sample 
  TTR.ROC = survivalROC(Stime=d$TTR.time, 
                        status=d$TTR.outcome,     
                        marker=d$MRS, 
                        predict.time=cutoff,
                        method="KM")
  return(TTR.ROC$AUC)
}
boot.OS.AUC <- function(data,indices,cutoff){ #define a function for bootstrap
  d <- data[indices,] # allows boot to select sample 
  OS.ROC = survivalROC(Stime=d$OS.time, 
                       status=d$OS.outcome,     
                       marker=d$MRS, 
                       predict.time=cutoff,
                       method="KM")
  return(OS.ROC$AUC)
}

##-----Training set----
#TTR
Training$TTR <- Surv(Training$TTR.time, Training$TTR.outcome!=0)
fit0 <- coxph(TTR ~ MRS,data=Training,na.action=na.omit )
TTR.12= survivalROC(Stime=Training$TTR.time, status=Training$TTR.outcome,     
                    marker=Training$MRS, predict.time=12,method="KM")   
TTR.24= survivalROC(Stime=Training$TTR.time, status=Training$TTR.outcome,     
                    marker=Training$MRS, predict.time=24,method="KM")
#OS
Training$OS <- Surv(Training$OS.time, Training$OS.outcome!=0)
fit0 <- coxph(OS ~ MRS,data=Training,na.action=na.omit )
OS.24= survivalROC(Stime=Training$OS.time, status=Training$OS.outcome,     
                   marker=Training$MRS, predict.time=24,method="KM")   
OS.36= survivalROC(Stime=Training$OS.time, status=Training$OS.outcome,     
                   marker=Training$MRS, predict.time=36,method="KM")
#Plot ROC
plot(TTR.12$FP, TTR.12$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#E41A1C",lwd=1)
lines(TTR.24$FP, TTR.24$TP,col="#377EB8",lwd=1)
lines(OS.24$FP, OS.24$TP,col="#4DAF4A",lwd=1)
lines(OS.36$FP, OS.36$TP,col="#FF7F00",lwd=1)
abline(0,1,lty=2,lwd=0.5)

#Export ROC
i <- 0.5
TTR.AUC <- c()
OS.AUC <- c()
time <- c()
repeat{if (i>96) break else{ 
  TTR.AUC.i= survivalROC(Stime=Training$TTR.time, status=Training$TTR.outcome,     
                         marker=Training$MRS, predict.time=i,method="KM")
  OS.AUC.i= survivalROC(Stime=Training$OS.time, status=Training$OS.outcome,     
                        marker=Training$MRS, predict.time=i,method="KM")
  TTR.AUC <- c(TTR.AUC, TTR.AUC.i$AUC)
  OS.AUC <- c(OS.AUC, OS.AUC.i$AUC)
  time <- c(time,i)
  i <- i+0.5
} 
}
TD.AUC <- cbind.data.frame(time,TTR.AUC,OS.AUC)
write.csv(TD.AUC,"TDROC-Training.csv")


# Calculate the AUC and the 95% CI
# 1-year TTR
TTR.12.Training.AUC <- boot(data = Training, cutoff=12, 
                            statistic = boot.TTR.AUC, R = 1000)
TTR.12$AUC
boot.ci(TTR.12.Training.AUC, conf = 0.95, type = "perc")
# 2-year TTR
TTR.24.Training.AUC <- boot(data = Training, cutoff=24, 
                            statistic = boot.TTR.AUC, R = 1000)
TTR.24$AUC
boot.ci(TTR.24.Training.AUC, conf = 0.95, type = "perc")
# 2-year OS
OS.24.Training.AUC <- boot(data = Training, cutoff=24, 
                           statistic = boot.OS.AUC, R = 1000)
OS.24$AUC
boot.ci(OS.24.Training.AUC, conf = 0.95, type = "perc")
# 3-year OS
OS.36.Training.AUC <- boot(data = Training, cutoff=36, 
                           statistic = boot.OS.AUC, R = 1000)
OS.36$AUC
boot.ci(OS.36.Training.AUC, conf = 0.95, type = "perc")


##-----Test set----
#TTR
Test$TTR <- Surv(Test$TTR.time, Test$TTR.outcome!=0)
fit0 <- coxph(TTR ~ MRS,data=Test,na.action=na.omit )
TTR.12= survivalROC(Stime=Test$TTR.time, status=Test$TTR.outcome,     
                    marker=Test$MRS, predict.time=12,method="KM")   
TTR.24= survivalROC(Stime=Test$TTR.time, status=Test$TTR.outcome,     
                    marker=Test$MRS, predict.time=24,method="KM")
#OS
Test$OS <- Surv(Test$OS.time, Test$OS.outcome!=0)
fit0 <- coxph(OS ~ MRS,data=Test,na.action=na.omit )
OS.24= survivalROC(Stime=Test$OS.time, status=Test$OS.outcome,     
                   marker=Test$MRS, predict.time=24,method="KM")   
OS.36= survivalROC(Stime=Test$OS.time, status=Test$OS.outcome,     
                   marker=Test$MRS, predict.time=36,method="KM")
#Plot ROC
plot(TTR.12$FP, TTR.12$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#E41A1C",lwd=1)
lines(TTR.24$FP, TTR.24$TP,col="#377EB8",lwd=1)
lines(OS.24$FP, OS.24$TP,col="#4DAF4A",lwd=1)
lines(OS.36$FP, OS.36$TP,col="#FF7F00",lwd=1)
abline(0,1,lty=2,lwd=0.5)

#Export ROC
i <- 0.5
TTR.AUC <- c()
OS.AUC <- c()
time <- c()
repeat{if (i>96) break else{ 
  TTR.AUC.i= survivalROC(Stime=Test$TTR.time, status=Test$TTR.outcome,     
                         marker=Test$MRS, predict.time=i,method="KM")
  OS.AUC.i= survivalROC(Stime=Test$OS.time, status=Test$OS.outcome,     
                        marker=Test$MRS, predict.time=i,method="KM")
  TTR.AUC <- c(TTR.AUC, TTR.AUC.i$AUC)
  OS.AUC <- c(OS.AUC, OS.AUC.i$AUC)
  time <- c(time,i)
  i <- i+0.5
} 
}
TD.AUC <- cbind.data.frame(time,TTR.AUC,OS.AUC)
write.csv(TD.AUC,"TDROC-Test.csv")

#Calculate the AUC and the 95% CI
# 1-year TTR
TTR.12.Test.AUC <- boot(data = Test, cutoff=12, 
                        statistic = boot.TTR.AUC, R = 1000)
TTR.12$AUC
boot.ci(TTR.12.Test.AUC, conf = 0.95, type = "perc")
# 2-year TTR
TTR.24.Test.AUC <- boot(data = Test, cutoff=24, 
                        statistic = boot.TTR.AUC, R = 1000)
TTR.24$AUC
boot.ci(TTR.24.Test.AUC, conf = 0.95, type = "perc")
# 2-year OS
OS.24.Test.AUC <- boot(data = Test, cutoff=24, 
                       statistic = boot.OS.AUC, R = 1000)
OS.24$AUC
boot.ci(OS.24.Test.AUC, conf = 0.95, type = "perc")
# 3-year OS
OS.36.Test.AUC <- boot(data = Test, cutoff=36, 
                       statistic = boot.OS.AUC, R = 1000)
OS.36$AUC
boot.ci(OS.36.Test.AUC, conf = 0.95, type = "perc")

##-----Primary Cohort----
#TTR
Primary$TTR <- Surv(Primary$TTR.time, Primary$TTR.outcome!=0)
fit0 <- coxph(TTR ~ MRS,data=Primary,na.action=na.omit )
TTR.12= survivalROC(Stime=Primary$TTR.time, status=Primary$TTR.outcome,     
                    marker=Primary$MRS, predict.time=12,method="KM")   
TTR.24= survivalROC(Stime=Primary$TTR.time, status=Primary$TTR.outcome,     
                    marker=Primary$MRS, predict.time=24,method="KM")
#OS
Primary$OS <- Surv(Primary$OS.time, Primary$OS.outcome!=0)
fit0 <- coxph(OS ~ MRS,data=Primary,na.action=na.omit )
OS.24= survivalROC(Stime=Primary$OS.time, status=Primary$OS.outcome,     
                   marker=Primary$MRS, predict.time=24,method="KM")   
OS.36= survivalROC(Stime=Primary$OS.time, status=Primary$OS.outcome,     
                   marker=Primary$MRS, predict.time=36,method="KM")
#Plot ROC
plot(TTR.12$FP, TTR.12$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#E41A1C",lwd=1)
lines(TTR.24$FP, TTR.24$TP,col="#377EB8",lwd=1)
lines(OS.24$FP, OS.24$TP,col="#4DAF4A",lwd=1)
lines(OS.36$FP, OS.36$TP,col="#FF7F00",lwd=1)
abline(0,1,lty=2,lwd=0.5)

#Export ROC
i <- 0.5
TTR.AUC <- c()
OS.AUC <- c()
time <- c()
repeat{if (i>96) break else{ 
  TTR.AUC.i= survivalROC(Stime=Primary$TTR.time, status=Primary$TTR.outcome,     
                         marker=Primary$MRS, predict.time=i,method="KM")
  OS.AUC.i= survivalROC(Stime=Primary$OS.time, status=Primary$OS.outcome,     
                        marker=Primary$MRS, predict.time=i,method="KM")
  TTR.AUC <- c(TTR.AUC, TTR.AUC.i$AUC)
  OS.AUC <- c(OS.AUC, OS.AUC.i$AUC)
  time <- c(time,i)
  i <- i+0.5
} 
}
TD.AUC <- cbind.data.frame(time,TTR.AUC,OS.AUC)
write.csv(TD.AUC,"TDROC-Primary.csv")

#Calculate the AUC and the 95% CI
# 1-year TTR
TTR.12.Primary.AUC <- boot(data = Primary, cutoff=12, 
                           statistic = boot.TTR.AUC, R = 1000)
TTR.12$AUC
boot.ci(TTR.12.Primary.AUC, conf = 0.95, type = "perc")
# 2-year TTR
TTR.24.Primary.AUC <- boot(data = Primary, cutoff=24, 
                           statistic = boot.TTR.AUC, R = 1000)
TTR.24$AUC
boot.ci(TTR.24.Primary.AUC, conf = 0.95, type = "perc")
# 2-year OS
OS.24.Primary.AUC <- boot(data = Primary, cutoff=24, 
                          statistic = boot.OS.AUC, R = 1000)
OS.24$AUC
boot.ci(OS.24.Primary.AUC, conf = 0.95, type = "perc")
# 3-year OS
OS.36.Primary.AUC <- boot(data = Primary, cutoff=36, 
                          statistic = boot.OS.AUC, R = 1000)
OS.36$AUC
boot.ci(OS.36.Primary.AUC, conf = 0.95, type = "perc")


##-----Internal validation cohort----
#TTR
Internal$TTR <- Surv(Internal$TTR.time, Internal$TTR.outcome!=0)
fit0 <- coxph(TTR ~ MRS,data=Internal,na.action=na.omit )
TTR.12= survivalROC(Stime=Internal$TTR.time, status=Internal$TTR.outcome,     
                    marker=Internal$MRS, predict.time=12,method="KM")   
TTR.24= survivalROC(Stime=Internal$TTR.time, status=Internal$TTR.outcome,     
                    marker=Internal$MRS, predict.time=24,method="KM")
#OS
Internal$OS <- Surv(Internal$OS.time, Internal$OS.outcome!=0)
fit0 <- coxph(OS ~ MRS,data=Internal,na.action=na.omit )
OS.24= survivalROC(Stime=Internal$OS.time, status=Internal$OS.outcome,     
                   marker=Internal$MRS, predict.time=24,method="KM")   
OS.36= survivalROC(Stime=Internal$OS.time, status=Internal$OS.outcome,     
                   marker=Internal$MRS, predict.time=36,method="KM")
#Plot ROC
plot(TTR.12$FP, TTR.12$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#E41A1C",lwd=1)
lines(TTR.24$FP, TTR.24$TP,col="#377EB8",lwd=1)
lines(OS.24$FP, OS.24$TP,col="#4DAF4A",lwd=1)
lines(OS.36$FP, OS.36$TP,col="#FF7F00",lwd=1)
abline(0,1,lty=2,lwd=0.5)

#Export ROC
i <- 0.5
TTR.AUC <- c()
OS.AUC <- c()
time <- c()
repeat{if (i>96) break else{ 
  TTR.AUC.i= survivalROC(Stime=Internal$TTR.time, status=Internal$TTR.outcome,     
                         marker=Internal$MRS, predict.time=i,method="KM")
  OS.AUC.i= survivalROC(Stime=Internal$OS.time, status=Internal$OS.outcome,     
                        marker=Internal$MRS, predict.time=i,method="KM")
  TTR.AUC <- c(TTR.AUC, TTR.AUC.i$AUC)
  OS.AUC <- c(OS.AUC, OS.AUC.i$AUC)
  time <- c(time,i)
  i <- i+0.5
} 
}
TD.AUC <- cbind.data.frame(time,TTR.AUC,OS.AUC)
write.csv(TD.AUC,"TDROC-Internal.csv")

#Calculate the AUC and the 95% CI
# 1-year TTR
TTR.12.Internal.AUC <- boot(data = Internal, cutoff=12, 
                            statistic = boot.TTR.AUC, R = 1000)
TTR.12$AUC
boot.ci(TTR.12.Internal.AUC, conf = 0.95, type = "perc")
# 2-year TTR
TTR.24.Internal.AUC <- boot(data = Internal, cutoff=24, 
                            statistic = boot.TTR.AUC, R = 1000)
TTR.24$AUC
boot.ci(TTR.24.Internal.AUC, conf = 0.95, type = "perc")
# 2-year OS
OS.24.Internal.AUC <- boot(data = Internal, cutoff=24, 
                           statistic = boot.OS.AUC, R = 1000)
OS.24$AUC
boot.ci(OS.24.Internal.AUC, conf = 0.95, type = "perc")
# 3-year OS
OS.36.Internal.AUC <- boot(data = Internal, cutoff=36, 
                           statistic = boot.OS.AUC, R = 1000)
OS.36$AUC
boot.ci(OS.36.Internal.AUC, conf = 0.95, type = "perc")



##-----External cohorts----
#TTR
External$TTR <- Surv(External$TTR.time, External$TTR.outcome!=0)
fit0 <- coxph(TTR ~ MRS,data=External,na.action=na.omit )
TTR.12= survivalROC(Stime=External$TTR.time, status=External$TTR.outcome,     
                    marker=External$MRS, predict.time=12,method="KM")   
TTR.24= survivalROC(Stime=External$TTR.time, status=External$TTR.outcome,     
                    marker=External$MRS, predict.time=24,method="KM")
#OS
External$OS <- Surv(External$OS.time, External$OS.outcome!=0)
fit0 <- coxph(OS ~ MRS,data=External,na.action=na.omit )
OS.24= survivalROC(Stime=External$OS.time, status=External$OS.outcome,     
                   marker=External$MRS, predict.time=24,method="KM")   
OS.36= survivalROC(Stime=External$OS.time, status=External$OS.outcome,     
                   marker=External$MRS, predict.time=36,method="KM")
#Plot ROC
plot(TTR.12$FP, TTR.12$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#E41A1C",lwd=1)
lines(TTR.24$FP, TTR.24$TP,col="#377EB8",lwd=1)
lines(OS.24$FP, OS.24$TP,col="#4DAF4A",lwd=1)
lines(OS.36$FP, OS.36$TP,col="#FF7F00",lwd=1)
abline(0,1,lty=2,lwd=0.5)

#Export ROC
i <- 0.5
TTR.AUC <- c()
OS.AUC <- c()
time <- c()
repeat{if (i>96) break else{ 
  TTR.AUC.i= survivalROC(Stime=External$TTR.time, status=External$TTR.outcome,     
                         marker=External$MRS, predict.time=i,method="KM")
  OS.AUC.i= survivalROC(Stime=External$OS.time, status=External$OS.outcome,     
                        marker=External$MRS, predict.time=i,method="KM")
  TTR.AUC <- c(TTR.AUC, TTR.AUC.i$AUC)
  OS.AUC <- c(OS.AUC, OS.AUC.i$AUC)
  time <- c(time,i)
  i <- i+0.5
} 
}
TD.AUC <- cbind.data.frame(time,TTR.AUC,OS.AUC)
write.csv(TD.AUC,"TDROC-External.csv")


#Calculate the AUC and the 95% CI
# 1-year TTR
TTR.12.External.AUC <- boot(data = External, cutoff=12, 
                            statistic = boot.TTR.AUC, R = 1000)
TTR.12$AUC
boot.ci(TTR.12.External.AUC, conf = 0.95, type = "perc")
# 2-year TTR
TTR.24.External.AUC <- boot(data = External, cutoff=24, 
                            statistic = boot.TTR.AUC, R = 1000)
TTR.24$AUC
boot.ci(TTR.24.External.AUC, conf = 0.95, type = "perc")
# 2-year OS
OS.24.External.AUC <- boot(data = External, cutoff=24, 
                           statistic = boot.OS.AUC, R = 1000)
OS.24$AUC
boot.ci(OS.24.External.AUC, conf = 0.95, type = "perc")
# 3-year OS
OS.36.External.AUC <- boot(data = External, cutoff=36, 
                           statistic = boot.OS.AUC, R = 1000)
OS.36$AUC
boot.ci(OS.36.External.AUC, conf = 0.95, type = "perc")

##-----THZP set----
#TTR
THZP$TTR <- Surv(THZP$TTR.time, THZP$TTR.outcome!=0)
fit0 <- coxph(TTR ~ MRS,data=THZP,na.action=na.omit )
TTR.12= survivalROC(Stime=THZP$TTR.time, status=THZP$TTR.outcome,     
                    marker=THZP$MRS, predict.time=12,method="KM")   
TTR.24= survivalROC(Stime=THZP$TTR.time, status=THZP$TTR.outcome,     
                    marker=THZP$MRS, predict.time=24,method="KM")
#OS
THZP$OS <- Surv(THZP$OS.time, THZP$OS.outcome!=0)
fit0 <- coxph(OS ~ MRS,data=THZP,na.action=na.omit )
OS.24= survivalROC(Stime=THZP$OS.time, status=THZP$OS.outcome,     
                   marker=THZP$MRS, predict.time=24,method="KM")   
OS.36= survivalROC(Stime=THZP$OS.time, status=THZP$OS.outcome,     
                   marker=THZP$MRS, predict.time=36,method="KM")
#Plot ROC
plot(TTR.12$FP, TTR.12$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#E41A1C",lwd=1)
lines(TTR.24$FP, TTR.24$TP,col="#377EB8",lwd=1)
lines(OS.24$FP, OS.24$TP,col="#4DAF4A",lwd=1)
lines(OS.36$FP, OS.36$TP,col="#FF7F00",lwd=1)
abline(0,1,lty=2,lwd=0.5)

#Export ROC
i <- 0.5
TTR.AUC <- c()
OS.AUC <- c()
time <- c()
repeat{if (i>96) break else{ 
  TTR.AUC.i= survivalROC(Stime=THZP$TTR.time, status=THZP$TTR.outcome,     
                         marker=THZP$MRS, predict.time=i,method="KM")
  OS.AUC.i= survivalROC(Stime=THZP$OS.time, status=THZP$OS.outcome,     
                        marker=THZP$MRS, predict.time=i,method="KM")
  TTR.AUC <- c(TTR.AUC, TTR.AUC.i$AUC)
  OS.AUC <- c(OS.AUC, OS.AUC.i$AUC)
  time <- c(time,i)
  i <- i+0.5
} 
}
TD.AUC <- cbind.data.frame(time,TTR.AUC,OS.AUC)
write.csv(TD.AUC,"TDROC-THZP.csv")


#Calculate the AUC and the 95% CI
# 1-year TTR
TTR.12.THZP.AUC <- boot(data = THZP, cutoff=12, 
                        statistic = boot.TTR.AUC, R = 1000)
TTR.12$AUC
boot.ci(TTR.12.THZP.AUC, conf = 0.95, type = "perc")
# 2-year TTR
TTR.24.THZP.AUC <- boot(data = THZP, cutoff=24, 
                        statistic = boot.TTR.AUC, R = 1000)
TTR.24$AUC
boot.ci(TTR.24.THZP.AUC, conf = 0.95, type = "perc")
# 2-year OS
OS.24.THZP.AUC <- boot(data = THZP, cutoff=24, 
                       statistic = boot.OS.AUC, R = 1000)
OS.24$AUC
boot.ci(OS.24.THZP.AUC, conf = 0.95, type = "perc")
# 3-year OS
OS.36.THZP.AUC <- boot(data = THZP, cutoff=36, 
                       statistic = boot.OS.AUC, R = 1000)
OS.36$AUC
boot.ci(OS.36.THZP.AUC, conf = 0.95, type = "perc")

##-----FUZH set----
#TTR
FUZH$TTR <- Surv(FUZH$TTR.time, FUZH$TTR.outcome!=0)
fit0 <- coxph(TTR ~ MRS,data=FUZH,na.action=na.omit )
TTR.12= survivalROC(Stime=FUZH$TTR.time, status=FUZH$TTR.outcome,     
                    marker=FUZH$MRS, predict.time=12,method="KM")   
TTR.24= survivalROC(Stime=FUZH$TTR.time, status=FUZH$TTR.outcome,     
                    marker=FUZH$MRS, predict.time=24,method="KM")
#OS
FUZH$OS <- Surv(FUZH$OS.time, FUZH$OS.outcome!=0)
fit0 <- coxph(OS ~ MRS,data=FUZH,na.action=na.omit )
OS.24= survivalROC(Stime=FUZH$OS.time, status=FUZH$OS.outcome,     
                   marker=FUZH$MRS, predict.time=24,method="KM")   
OS.36= survivalROC(Stime=FUZH$OS.time, status=FUZH$OS.outcome,     
                   marker=FUZH$MRS, predict.time=36,method="KM")
#Plot ROC
plot(TTR.12$FP, TTR.12$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#E41A1C",lwd=1)
lines(TTR.24$FP, TTR.24$TP,col="#377EB8",lwd=1)
lines(OS.24$FP, OS.24$TP,col="#4DAF4A",lwd=1)
lines(OS.36$FP, OS.36$TP,col="#FF7F00",lwd=1)
abline(0,1,lty=2,lwd=0.5)

#Export ROC
i <- 0.5
TTR.AUC <- c()
OS.AUC <- c()
time <- c()
repeat{if (i>96) break else{ 
  TTR.AUC.i= survivalROC(Stime=FUZH$TTR.time, status=FUZH$TTR.outcome,     
                         marker=FUZH$MRS, predict.time=i,method="KM")
  OS.AUC.i= survivalROC(Stime=FUZH$OS.time, status=FUZH$OS.outcome,     
                        marker=FUZH$MRS, predict.time=i,method="KM")
  TTR.AUC <- c(TTR.AUC, TTR.AUC.i$AUC)
  OS.AUC <- c(OS.AUC, OS.AUC.i$AUC)
  time <- c(time,i)
  i <- i+0.5
} 
}
TD.AUC <- cbind.data.frame(time,TTR.AUC,OS.AUC)
write.csv(TD.AUC,"TDROC-FUZH.csv")


#Calculate the AUC and the 95% CI
# 1-year TTR
TTR.12.FUZH.AUC <- boot(data = FUZH, cutoff=12, 
                        statistic = boot.TTR.AUC, R = 1000)
TTR.12$AUC
boot.ci(TTR.12.FUZH.AUC, conf = 0.95, type = "perc")
# 2-year TTR
TTR.24.FUZH.AUC <- boot(data = FUZH, cutoff=24, 
                        statistic = boot.TTR.AUC, R = 1000)
TTR.24$AUC
boot.ci(TTR.24.FUZH.AUC, conf = 0.95, type = "perc")
# 2-year OS
OS.24.FUZH.AUC <- boot(data = FUZH, cutoff=24, 
                       statistic = boot.OS.AUC, R = 1000)
OS.24$AUC
boot.ci(OS.24.FUZH.AUC, conf = 0.95, type = "perc")
# 3-year OS
OS.36.FUZH.AUC <- boot(data = FUZH, cutoff=36, 
                       statistic = boot.OS.AUC, R = 1000)
OS.36$AUC
boot.ci(OS.36.FUZH.AUC, conf = 0.95, type = "perc")




graphics.off()
rm(list = ls())