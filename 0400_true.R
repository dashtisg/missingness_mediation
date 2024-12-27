#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name:      0400_true 
#Purpose:	       This is the R script for simulating a large dataset (1M) and obtaining the true values 
#Author:         Ghazaleh Dashti 
#Final revision on 19/08/2024
#Code further annotated on 24/12/2024

#############################VAHCS EXAMPLE###################################### 
# Simulation study based on an example from VAHCS 
# exposure X: Depression and anxiety at least 2 waves during adolescence 
#              (as a proxy for the severity of depression and anxiety)

# intermediate confounder L: any live birth by age 24 

# mediator M: Depression and anxiety at least 2 waves in young adulthood 
#              (as a proxy for the severity of depression and anxiety)

# outcome Y: depression and anxiety in adulthood 

# baseline confounders Zs: Z1: sex registered at birth 
#                          Z2: cannabis use 
#                          Z3: adolescent (wave 2-6) antisocial behaviour 

# Auxiliary variable A: parental smoking 

#############################IMPORTANT NOTES#################################### 
# The general set-up of the codes is for two scenarios, simple and complex. 
# This code simulates data under the complex scenario. 
# In the final version of the manuscript only results from the complex scenario are presented.
# The difference between the two scenarios is in the presence of an intermediate confounder: 
# In the "simple" scenario there are no intermediate confounders  
# In the "complex" scenario any live birth by age 24 is considered as an intermediate confounder 

#############################GENERAL SET UP#####################################

rm(list=ls()) 
start_time <- Sys.time() #time the process
#setwd("xx")

#load related packages
packages <- c("parallel","foreach", "boot", "mice","devtools","data.table", "dplyr", "readstata13", "tidyverse", "zoo", "bootImpute","smcfcs","mitools","plyr","micemd","doParallel","doRNG")

for (package in packages) {
  library(package, character.only=T)
}

set.seed(1235659848)
nsim <- 1000000

## Create directories
if(!dir.exists(paste0(getwd(),"/SIM/4_true/complex"))){
  dir.create(paste0(getwd(),"/SIM/4_true/complex"))
}

#############################SIMULATE THE DATA##################################
parms<-read.csv("parms_complex.csv")
attach(parms) 

A<-rbinom(nsim,1,inv.logit(a_int))

Z1<-rbinom(nsim,1,inv.logit(z1_int))

Z2<-rbinom(nsim,1,inv.logit(z2_int+z2_a*A))

Z3<-rbinom(nsim,1,inv.logit(z3_int+z3_a*A))

X<-rbinom(nsim,1,inv.logit(x_int+x_z1*Z1+x_z2*Z2+x_z3*Z3+x_z1z3*Z1*Z3))

L<-rbinom(nsim,1,inv.logit(l_int+l_x*X+l_z1*Z1+l_z2*Z2+l_z3*Z3+l_z1z3*Z1*Z3))

M<-rbinom(nsim,1,inv.logit(m_int+m_x*X+m_l*L+m_z1*Z1+m_z2*Z2+m_z3*Z3+m_xl*X*L+m_z1z3*Z1*Z3))

Y<-rbinom(nsim,1,inv.logit(y_int+y_x*X+y_z1*Z1+y_z2*Z2+y_z3*Z3+y_l*L+y_m*M+y_xl*X*L+y_xm*X*M+y_z1z3*Z1*Z3))

data <-data.frame(A,Z1,Z2,Z3,X,L,M,Y)  

#############################RUN THE ANALYSIS ##################################
mcsim<-1000
# Model for the intermediate confounder on exposure and confounders  
fitL<-glm(L~X+Z1+Z2+Z3+Z1*Z3, family="binomial", data=data) 

# Model for the mediator on exposure and confounders  
fitM<-glm(M~X+L+Z1+Z2+Z3+X*L+Z1*Z3, family="binomial", data=data)

# Model for the outcome on expousre, mediator, and confounders 
fitY<-glm(Y~X+L+M+Z1+Z2+Z3+X*M+X*L+Z1*Z3, family="binomial", data=data)

# Replicate dataset 
dat2<-data
dat2[,c("X","L","M","Y")]<-NA_integer_ # replace X, L, M and Y to missing in replicated dataset 
dat2<-zoo::coredata(dat2)[rep(seq(nrow(dat2)),mcsim),]
n<-nrow(dat2)

# E[Y1M1]  
dat2$X<-1 
dat2$L<-rbinom(n,1,predict(fitL,newdata=dat2,type="response",na.rm = T))  
dat2$M<-rbinom(n,1,predict(fitM,newdata=dat2,type="response",na.rm = T)) 
Y11<-mean(predict(fitY,newdata=dat2,type="response",na.rm = T)) 

# E[Y0M0]  
dat2$X<-0
dat2$L<-rbinom(n,1,predict(fitL,newdata=dat2,type="response",na.rm = T))   
dat2$M<-rbinom(n,1,predict(fitM,newdata=dat2,type="response",na.rm = T)) 
Y00<-mean(predict(fitY,newdata=dat2,type="response",na.rm = T))  

# E[Y1M0]
dat2$X<-0 
dat2$M<-rbinom(n,1,predict(fitM,newdata=dat2,type="response",na.rm = T)) 
dat2$X<-1 
dat2$L<-rbinom(n,1,predict(fitL,newdata=dat2,type="response",na.rm = T)) 
Y10<-mean(predict(fitY,newdata=dat2,type="response",na.rm = T))  

# Estimate effects 
DE_RD <-(Y10 - Y00)
IE_RD <-(Y11 - Y10) 

res<-c(DE_RD,IE_RD)   

intmed_out_complex <- data.frame(DE = DE_RD, IE = IE_RD) 
comp_file <- paste0(getwd(),"/SIM/4_true/complex/true_complex.csv") 
write.csv(intmed_out_complex, file = comp_file, row.names=FALSE)  