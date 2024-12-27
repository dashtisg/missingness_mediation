#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name:    0600_applied_ana 
#Purpose:	     This is the R script for analyses of real VAHCS data 
#Author:       Ghazaleh Dashti 
#Last updated on 19/08/2024
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

#############################GENERAL SET UP#####################################
rm(list=ls()) 
#setwd("xx") 
#load related packages
packages <- c("parallel","foreach", "boot", "mice","devtools","data.table", "dplyr", "readstata13", "tidyverse", "zoo", "bootImpute","smcfcs","mitools","plyr","micemd","doParallel","doRNG","rstream")

for (package in packages) {
  library(package, character.only=T)
}
### create the folder to save the results 
if(!dir.exists(paste0(getwd(),"/SIM/9_applied"))){
  dir.create(paste0(getwd(),"/SIM/9_applied"))
}  
 

### load the data 
dat<-read.dta13("vahcs_cmd_med.dta",nonint.factors=T)
names(dat) <- toupper(names(dat)) 

prop.table(table(is.na(dat$A))) # exclude 
dat<- dat[!is.na(dat$A),]   
dat<- dat[!is.na(dat$Z1),]  
dat <- dat[,1:8]
dat$Z1 <- ifelse(dat$Z1 == levels(dat$Z1)[1], 0, 1) 

# analytic sample 1828 

# source the necessary functions 
source("func_intmed_complex.R")   

#### Fixed parameters
# bootstrap 
nboot<-200
# number of imputations 
m<-50
M<-2
imp_runs <- seq(1, m)   
# number of cores for parallel computing 
outer_cores<- 5
ncores<- 5
inner_cores<- 5
### 
nseeds<-15

initial_seed <- 7864523

set.seed(initial_seed)

# Initialize the rstream object with the initial seed
stream <- new("rstream.mrg32k3a")

# Function to get a seed from the stream
generate_seed <- function() {
  as.integer(runif(1, 1, .Machine$integer.max))
}

# Initialize a vector to store seeds
seeds <- integer(nseeds)

# Generate the seeds
for (i in 1:nseeds) {
  rstream.nextsubstream(stream)
  seeds[i] <- generate_seed()
}
#############################DATA ANALYSIS#####################################
## CC 
dat_cc<-na.omit(dat) 
seed <- seeds[1]
set.seed(seed) 
int_CC<-intmed_fun(dat = dat_cc)  
comp_file <- paste0(getwd(),"/SIM/9_applied/int_CC.csv")
write.csv(int_CC, file = comp_file, row.names=FALSE)   

## noLZYMI 
setseed<-seeds[2]
int_noLZYMIBoot <- int_noLZYMIBoot_fun(dat = dat, setseed = setseed) 
comp_file <- paste0(getwd(),"/SIM/9_applied/int_noLZYMIBoot.csv")
write.csv(int_noLZYMIBoot, file = comp_file, row.names=FALSE)   

setseed<-seeds[3]
int_noLZYBootMI <- int_noLZYBootMI_fun(dat = dat, setseed = setseed)
comp_file <- paste0(getwd(),"/SIM/9_applied/int_noLZYBootMI.csv")
write.csv(int_noLZYBootMI, file = comp_file, row.names=FALSE)  

## noYMI 
setseed<-seeds[4]
int_noYMIBoot <- int_noYMIBoot_fun(dat = dat, setseed = setseed)   
comp_file <- paste0(getwd(),"/SIM/9_applied/int_noYMIBoot.csv")
write.csv(int_noYMIBoot , file = comp_file, row.names=FALSE)   

setseed<-seeds[5]
int_noYBootMI <- int_noYBootMI_fun(dat = dat, setseed = setseed)  
comp_file <- paste0(getwd(),"/SIM/9_applied/int_noYBootMI.csv")
write.csv(int_noYBootMI, file = comp_file, row.names=FALSE) 

## noint 
setseed<-seeds[6]
int_nointMIBoot <- int_nointMIBoot_fun(dat = dat, setseed = setseed)  
comp_file <- paste0(getwd(),"/SIM/9_applied/int_nointMIBoot.csv")
write.csv(int_nointMIBoot, file = comp_file, row.names=FALSE)    

setseed<-seeds[7]
int_nointBootMI <- int_nointBootMI_fun(dat = dat, setseed = setseed) 
comp_file <- paste0(getwd(),"/SIM/9_applied/int_nointBootMI.csv")
write.csv(int_nointBootMI, file = comp_file, row.names=FALSE)   

## Mcomp 
setseed<-seeds[8]
int_McompMIBoot <- int_McompMIBoot_fun(dat = dat, setseed = setseed) 
comp_file <- paste0(getwd(),"/SIM/9_applied/int_McompMIBoot.csv")
write.csv(int_McompMIBoot, file = comp_file, row.names=FALSE)     

setseed<-seeds[9]
int_McompBootMI <- int_McompBootMI_fun(dat = dat, setseed = setseed)
comp_file <- paste0(getwd(),"/SIM/9_applied/int_McompBootMI .csv")
write.csv(int_McompBootMI, file = comp_file, row.names=FALSE)  

## Ycomp 
setseed<-seeds[10]
int_YcompMIBoot <- int_YcompMIBoot_fun(dat = dat, setseed = setseed) 
comp_file <- paste0(getwd(),"/SIM/9_applied/int_YcompMIBoot.csv")
write.csv(int_YcompMIBoot, file = comp_file, row.names=FALSE)  

setseed<-seeds[11]
int_YcompBootMI <- int_YcompBootMI_fun(dat = dat, setseed = setseed)  
comp_file <- paste0(getwd(),"/SIM/9_applied/int_YcompBootMI.csv")
write.csv(int_YcompBootMI, file = comp_file, row.names=FALSE)   

## higherint 
setseed<-seeds[12]
int_higherintMIBoot <- int_higherintMIBoot_fun(dat = dat, setseed = setseed) 
comp_file <- paste0(getwd(),"/SIM/9_applied/int_higherintMIBoot.csv")
write.csv(int_higherintMIBoot , file = comp_file, row.names=FALSE)  

setseed<-seeds[13]
int_higherintBootMI <- int_higherintBootMI_fun(dat = dat, setseed = setseed) 
comp_file <- paste0(getwd(),"/SIM/9_applied/int_higherintBootMI.csv")
write.csv(int_higherintBootMI, file = comp_file, row.names=FALSE) 

## SMC FCS 
setseed<-seeds[14] 
int_smcfcsMIBoot <-int_smcfcsMIBoot_fun(dat = dat, setseed = setseed)
comp_file <- paste0(getwd(),"/SIM/9_applied/int_smcfcsMIBoot.csv")
write.csv(int_smcfcsMIBoot, file = comp_file, row.names=FALSE)  

setseed<-seeds[15]
int_smcfcsBootMI <-int_smcfcsBootMI_fun(dat = dat, setseed = setseed)  
comp_file <- paste0(getwd(),"/SIM/9_applied/int_smcfcsBootMI.csv")
write.csv(int_smcfcsBootMI, file = comp_file, row.names=FALSE)   


### forest plot for applied results 
#setwd("xx")
results_CC<-read.csv("int_CC.csv")

results_higherintBootMI<-read.csv("int_higherintBootMI.csv")
results_higherintBootMI$SE<-sqrt(results_higherintBootMI$SE) 
results_higherintMIBoot<-read.csv("int_higherintMIBoot.csv")

results_McompBootMI<-read.csv("int_McompBootMI.csv")
results_McompBootMI$SE<-sqrt(results_McompBootMI$SE) 
results_McompMIBoot<-read.csv("int_McompMIBoot.csv")

results_nointBootMI<-read.csv("int_nointBootMI.csv")
results_nointBootMI$SE<-sqrt(results_nointBootMI$SE) 
results_nointMIBoot<-read.csv("int_nointMIBoot.csv")

results_noMYBootMI<-read.csv("int_noMYBootMI.csv")
results_noMYBootMI$SE<-sqrt(results_noMYBootMI$SE) 
results_noMYMIBoot<-read.csv("int_noMYMIBoot.csv")

results_noYBootMI<-read.csv("int_noYBootMI.csv")
results_noYBootMI$SE<-sqrt(results_noYBootMI$SE) 
results_noYMIBoot<-read.csv("int_noYMIBoot.csv")

results_smcfcsBootMI<-read.csv("int_smcfcsBootMI.csv")
results_smcfcsBootMI$SE<-sqrt(results_smcfcsBootMI$SE) 
results_smcfcsMIBoot<-read.csv("int_smcfcsMIBoot.csv")

results_YcompBootMI<-read.csv("int_YcompBootMI.csv")
results_YcompBootMI$SE<-sqrt(results_YcompBootMI$SE) 
results_YcompMIBoot<-read.csv("int_YcompMIBoot.csv")


results_de <- rbind(results_CC[1,], 
                    results_noMYBootMI[1,], 
                    results_noMYMIBoot[1,], 
                    results_noYBootMI[1,], 
                    results_noYMIBoot[1,],
                    results_nointBootMI[1,],
                    results_nointMIBoot[1,],
                    results_McompBootMI[1,], 
                    results_McompMIBoot[1,],
                    results_YcompBootMI[1,],
                    results_YcompMIBoot[1,], 
                    results_higherintBootMI[1,],
                    results_higherintMIBoot[1,],
                    results_smcfcsMIBoot[1,],
                    results_smcfcsBootMI[1,]) 

results_ie <- rbind(results_CC[2,], 
                    results_noMYBootMI[2,], 
                    results_noMYMIBoot[2,], 
                    results_noYBootMI[2,], 
                    results_noYMIBoot[2,],
                    results_nointBootMI[2,],
                    results_nointMIBoot[2,],
                    results_McompBootMI[2,], 
                    results_McompMIBoot[2,],
                    results_YcompBootMI[2,],
                    results_YcompMIBoot[2,], 
                    results_higherintBootMI[2,],
                    results_higherintMIBoot[2,],
                    results_smcfcsMIBoot[2,],
                    results_smcfcsBootMI[2,]) 


write.csv(results_de, file="results_de.csv", row.names=FALSE) 
write.csv(results_ie, file="results_ie.csv", row.names=FALSE) 


### forest plot for applied results 
library(ggplot2) 
library(ggpubr) 

results_de$method<- c("CCA",
                   "MI-noLZY (BootMI)",
                   "MI-noLZY (MIBoot)", 
                   "MI-noY (BootMI)",
                   "MI-noY (MIBoot)",
                   "MI-noint (BootMI)", 
                   "MI-noint (MIBoot)", 
                   "MI-LZint (BootMI)", 
                   "MI-LZint (MIBoot)", 
                   "MI-Yint (BootMI)",
                   "MI-Yint (MIBoot)", 
                   "MI-allint (BootMI)",
                   "MI-allint (MIBoot)", 
                   "MI-SMCFCS (BootMI)",
                   "MI-SMCFCS (MIBoot)") 

results_ie$method<- c("CCA",
                      "MI-noLZY (BootMI)",
                      "MI-noLZY (MIBoot)", 
                      "MI-noY (BootMI)",
                      "MI-noY (MIBoot)",
                      "MI-noint (BootMI)", 
                      "MI-noint (MIBoot)", 
                      "MI-LZint (BootMI)", 
                      "MI-LZint (MIBoot)", 
                      "MI-Yint (BootMI)",
                      "MI-Yint (MIBoot)", 
                      "MI-allint (BootMI)",
                      "MI-allint (MIBoot)", 
                      "MI-SMCFCS (BootMI)",
                      "MI-SMCFCS (MIBoot)") 

results_de$method <- factor(results_de$method, levels = c("MI-SMCFCS (MIBoot)",
                                                          "MI-SMCFCS (BootMI)",
                                                          "MI-allint (MIBoot)",
                                                          "MI-allint (BootMI)",
                                                          "MI-Yint (MIBoot)", 
                                                          "MI-Yint (BootMI)",
                                                          "MI-LZint (MIBoot)",
                                                          "MI-LZint (BootMI)", 
                                                          "MI-noint (MIBoot)",
                                                          "MI-noint (BootMI)", 
                                                          "MI-noY (MIBoot)",
                                                          "MI-noY (BootMI)",
                                                          "MI-noLZY (MIBoot)", 
                                                          "MI-noLZY (BootMI)",
                                                          "CCA"))

results_ie$method <- factor(results_de$method, levels = c("MI-SMCFCS (MIBoot)",
                                                          "MI-SMCFCS (BootMI)",
                                                          "MI-allint (MIBoot)",
                                                          "MI-allint (BootMI)",
                                                          "MI-Yint (MIBoot)", 
                                                          "MI-Yint (BootMI)",
                                                          "MI-LZint (MIBoot)",
                                                          "MI-LZint (BootMI)", 
                                                          "MI-noint (MIBoot)",
                                                          "MI-noint (BootMI)", 
                                                          "MI-noY (MIBoot)",
                                                          "MI-noY (BootMI)",
                                                          "MI-noLZY (MIBoot)", 
                                                          "MI-noLZY (BootMI)",
                                                          "CCA"))

results_ie$Estimate<-results_ie$Estimate*100
results_ie$`X95CIlow`<-(results_ie$`X95CIlow`)*100
results_ie$`X95CIupp`<-(results_ie$`X95CIupp`)*100
results_ie$SE<-(results_ie$SE)*100

results_de$Estimate<-results_de$Estimate*100
results_de$`X95CIlow`<-(results_de$`X95CIlow`)*100
results_de$`X95CIupp`<-(results_de$`X95CIupp`)*100
results_de$SE<-(results_de$SE)*100

library(ggplot2)
library(dplyr)
library(ggpubr)  # Ensure ggpubr is loaded for ggarrange

library(ggplot2)
library(dplyr)
library(ggpubr)  # Ensure ggpubr is loaded for ggarrange

dev.off()

# Calculate nudge values for SEs based on the maximum estimate
max_estimate <- max(results_de$Estimate)
nudge_value <- 5 # Increase the nudge value 

IE<- ggplot(data=results_ie, aes(y=method, x=Estimate, xmin=`X95CIlow`, xmax=`X95CIupp` )) +
  geom_point(aes(color = method))   + 
  geom_errorbarh(height=.1,aes(color = method)) +
  labs(title='Indirect effect', x='Risk difference (95%CI) per 100', y = 'Missing data method') +
  geom_text(aes(label = sprintf("(SE %.2f)", SE), x = max_estimate + nudge_value), hjust = 0, size = 5)+  
  scale_x_continuous(limits=c(-1, 17)) + theme(text = element_text(size =15))    

DE<- ggplot(data=results_de, aes(y=method, x=Estimate, xmin=`X95CIlow`, xmax=`X95CIupp` )) +
  geom_point(aes(color = method))   + 
  geom_errorbarh(height=.1,aes(color = method)) +
  labs(title='Direct effect', x='Risk difference (95%CI) per 100', y = '') +
  geom_text(aes(label = sprintf("(SE %.2f)", SE), x = max_estimate + nudge_value), hjust = 0, size = 5)+  
  scale_x_continuous(limits=c(-1, 17)) + theme(text = element_text(size =15))     

applied <- ggarrange(IE, DE,
                ncol = 2, nrow = 1, common.legend = T) 





