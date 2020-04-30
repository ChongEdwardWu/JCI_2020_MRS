########################
### Install packages ###
########################

#### Install necessary R packages ###
## if you have not installed the necessary R packages, run the following codes ##
# if (!requireNamespace("rms", quietly = TRUE))
#   install.packages("rms")
# if (!requireNamespace("RColorBrewer", quietly = TRUE))
#   install.packages("RColorBrewer")
# if (!requireNamespace("survivalROC", quietly = TRUE))
#   install.packages("survivalROC")
# if (!requireNamespace("nomogramEx", quietly = TRUE))
#   install.packages("nomogramEx")
# if (!requireNamespace("boot", quietly = TRUE))
#   install.packages("boot")


### load packages ###
require(rms)
require(survivalROC)
require(nomogramEx)
require(boot)


#############################
### load data and prepare ###
#############################

## Set your working directory (or use Control + Shift + H to choose it) ##
## The working directory should contain the RData file "Nomogram construction.RData", which contains the cohort data, and a R source file stdca.R" for DCA 
setwd("...")


# load data from a csv table, or alternatively, from a RData file
# Primary <- read.csv("PrimaryCohort.csv",header=TRUE, as.is=TRUE, na.strings=c("","NA"))
load("Nomogram construction.RData")

#set the coloar palette
palette(MyPal) 

# preprocess data
max.time <- 96
units(Primary$TTR.time) <- 'Month'
units(Primary$OS.time) <- 'Month'

Primary$TTR.outcome[Primary$TTR.time>max.time ] <- 0
Primary$TTR.time[Primary$TTR.time>max.time ] <- max.time
Primary$TTR <- Surv(Primary$TTR.time, Primary$TTR.outcome!=0)

Primary$OS.outcome[Primary$OS.time>max.time ] <- 0
Primary$OS.time[Primary$OS.time>max.time ] <- max.time
Primary$OS <- Surv(Primary$OS.time, Primary$OS.outcome!=0)



### determine cuf-points
Primary$age <- cut(Primary$Age, breaks=c(0,55,100))
Primary$sex <- as.factor(Primary$Sex)
Primary$dm <- cut(Primary$Diameter,breaks=c(0,4,10,Inf))
Primary$afp <- cut(Primary$AFP, breaks=c(0,25,Inf))
Primary$alt <- cut(Primary$ALT, breaks=c(0,40,Inf))
Primary$ps <- as.factor(Primary$No.path.stage)
Primary$cp <- as.factor(Primary$Child.Pugh)
Primary$tn <- cut(Primary$No.Tumors,breaks=c(0,1,Inf))
Primary$vi <- as.factor(Primary$Vasc.invad)
Primary$bclc <- as.factor(Primary$BCLC)
Primary$cir <- as.factor(Primary$Cirhosis)
Primary$cap <- as.factor(Primary$Capsule)
Primary$hbv <- as.factor(Primary$HBsAgn)
Primary$hcv <- as.factor(Primary$HCVAb)
Primary$albi <- cut(Primary$ALBI.grade, breaks=c(0,1,Inf))
Primary$t <- as.factor(Primary$T)
Primary$n <- as.factor(Primary$N)
Primary$m <- as.factor(Primary$M)
Primary$tnm <- as.factor(Primary$TNM)


### set data distribution
options(datadist=NULL)
dd <- datadist(Primary); options(datadist='dd')

#####################################
### Primary Cox regression models ###
#####################################

### construct survival data
Primary1 <- na.omit(Primary)
PrimaryTTR <- Surv(Primary1$TTR.time, Primary1$TTR.outcome!=0)
PrimaryOS <- Surv(Primary1$OS.time, Primary1$OS.outcome!=0)


### Primary TTR regression model
f0.TTR <- cph(PrimaryTTR ~ MRS + age + sex + dm + afp + alt +
ps + cp + tn + vi + cir+ cap + hbv + hcv + albi,
data=Primary1, x=TRUE,y=TRUE,surv=TRUE,time.inc=12)

### Primary OS regression model
f0.OS <- cph(PrimaryOS ~ MRS + age + sex + dm + afp + alt +
ps + cp + tn + vi + cir + cap + hbv + hcv + albi,
 data=Primary1, x=TRUE,y=TRUE,surv=TRUE,time.inc=12)

print(f0.TTR, coefs=FALSE)
print(f0.OS, coefs=FALSE)

### select predictors using AIC
set.seed(222)# so can reproduce results
v0.TTR <- validate(f0.TTR, B=300,bw=TRUE)
v0.TTR
v0.OS <- validate(f0.OS, B=300,bw=TRUE)
v0.OS



######################################################
###################  TTR Nomogram  ###################
######################################################

f.TTR <- cph(TTR ~ MRS + dm + vi + cir, 
 data=Primary, x=TRUE,y=TRUE,surv=TRUE,time.inc=12*2)


## -----Checking PH---------------------------------------------------------
z.TTR <- predict(f.TTR, type='terms')
# required x=T above to store design matrix
f.short.TTR <- cph(TTR ~ z.TTR, x=TRUE, y=TRUE, data = Primary)
# store raw x, y so can get residuals
phtest.TTR <- cox.zph(f.short.TTR, transform='identity')
phtest.TTR
plot(phtest.TTR, var='cir')

## -----Describing Predictor Effects-----------------------------------------
ggplot(Predict(f.TTR), sepdiscrete='vertical', nlevels=1,
 vnames='names') 

## -----Validating the Model-------------------------------------------------
set.seed(324)# so can reproduce results
v.TTR <- validate(f.TTR, B=300)
v.TTR
cal.TTR <- calibrate(f.TTR, B=30, u=2*12, cmethod='KM',
 cuts = c(.2,.4,.6,.8), surv=TRUE,time.inc=u)
plot(cal.TTR, subtitles=TRUE, riskdist=FALSE,xlim=c(0,1),col="#E41A1C",lwd=1)


## -----Hazard ratios and multi-level confidence bars for effects of predictors 
## in model------------------------------------------------------------------
plot(summary(f.TTR),log=TRUE, main='')


## -----Presenting the Model with nomogram-----------------------------------
RF.surv<- Survival(f.TTR)
RF.surv1 <- function(x) RF.surv(1*12,lp=x)
RF.surv2 <- function(x) RF.surv(2*12,lp=x)
quan<- Quantile(f.TTR)
med.TTR <- function(x) quan(lp=x)/12
ss<- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95)

nom.TTR <- nomogram(f.TTR, lp=FALSE, fun=list(RF.surv1, RF.surv2, med.TTR),
funlabel=c('1-year TTR','2-year TTR',
 'Median TTR Time (years)'),
fun.at=list(ss, ss, c(.5,1:7)))
plot(nom.TTR, xfrac=.55, lmgp=.4)



## -----Export nomogram parameters-----------------------------------
TTR.nomoEx <- nomogramEx(nomo=nom.TTR,np=3,digit=9)

Primary$NomoR.MRS <- Primary$MRS
Primary$NomoR.dm <-TTR.nomoEx[3][[1]][as.numeric(Primary$dm),]
Primary$NomoR.vi <-TTR.nomoEx[4][[1]][as.numeric(Primary$vi),]
Primary$NomoR.cir <-TTR.nomoEx[5][[1]][as.numeric(Primary$cir),]
Primary$NomoR.total <- Primary$NomoR.MRS + Primary$NomoR.dm +
                        Primary$NomoR.vi + Primary$NomoR.cir 
hist(Primary$NomoR.total, breaks = 24)
Primary$NomoR.group <- as.numeric(cut(Primary$NomoR.total,
                                          breaks=c(-Inf,80,160,Inf)))


#####################################################
###################  OS Nomogram  ###################
#####################################################

f.OS <- cph(OS ~ MRS + dm + vi + cir + tn + cp, 
data=Primary, x=TRUE,y=TRUE,surv=TRUE,time.inc=12*3)



## -----Checking PH---------------------------------------------------------
z.OS <- predict(f.OS, type='terms')
# required x=T above to store design matrix
f.short.OS <- cph(OS ~ z.OS, x=TRUE, y=TRUE, data = Primary)
# store raw x, y so can get residuals
phtest.OS <- cox.zph(f.short.OS, transform='identity')
phtest.OS
plot(phtest.OS, var='cir')

## -----Describing Predictor Effects-----------------------------------------
ggplot(Predict(f.OS), sepdiscrete='vertical', nlevels=1,
 vnames='names') 

## -----Validating the Model-------------------------------------------------
set.seed(324)# so can reproduce results
v.OS <- validate(f.OS, B=300)
v.OS
cal.OS <- calibrate(f.OS, B=30, u=3*12, cmethod='KM',
 cuts = c(.2,.4,.6,.8), surv=TRUE,time.inc=u)
plot(cal.OS, subtitles=TRUE, riskdist=FALSE,xlim=c(0,1),col="#E41A1C",lwd=1)


## -----Hazard ratios and multi-level confidence bars for effects of predictors 
## in model------------------------------------------------------------------
plot(summary(f.OS),log=TRUE, main='')


## -----Presenting the Model with nomogram-----------------------------------
OS.surv<- Survival(f.OS)
OS.surv2 <- function(x) OS.surv(2*12,lp=x)
OS.surv3 <- function(x) OS.surv(3*12,lp=x)
quan<- Quantile(f.OS)
med.OS <- function(x) quan(lp=x)/12
ss<- c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95)

nom.OS <- nomogram(f.OS, lp=FALSE, fun=list(OS.surv2, OS.surv3, med.OS),
                   funlabel=c('2-year OS','3-year OS',
                              'Median OS Time (years)'),
                   fun.at=list(ss, ss, c(.5,1:7)))
plot(nom.OS, xfrac=.55, lmgp=.4)


## -----Export nomogram parameters-----------------------------------
OS.nomoEx <- nomogramEx(nomo=nom.OS,np=3,digit=9)

Primary$NomoO.MRS <- Primary$MRS
Primary$NomoO.dm <-OS.nomoEx[3][[1]][as.numeric(Primary$dm),]
Primary$NomoO.vi <-OS.nomoEx[4][[1]][as.numeric(Primary$vi),]
Primary$NomoO.cir <-OS.nomoEx[5][[1]][as.numeric(Primary$cir),]
Primary$NomoO.tn <-OS.nomoEx[6][[1]][as.numeric(Primary$tn),]
Primary$NomoO.cp <-OS.nomoEx[7][[1]][as.numeric(Primary$cp),]
Primary$NomoO.total <- Primary$NomoO.MRS + Primary$NomoO.dm +
                      Primary$NomoO.vi + Primary$NomoO.cir +
                      Primary$NomoO.cp + Primary$NomoO.tn

hist(Primary$NomoO.total,28)
Primary$NomoO.group <- as.numeric(cut(Primary$NomoO.total,
                                          breaks=c(-Inf,80,160,Inf)))
length(Primary$NomoO.group[Primary$NomoO.group==3])
write.csv(Primary,"Primary-Nomo.csv")

graphics.off()
rm(list = ls())

