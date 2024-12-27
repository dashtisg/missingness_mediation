#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects

#############################R SCRIPT DESCRIPTION############################### 
#File name: 0200_sim_data 
#Purpose: This is the R script for simulating the data 
#Author: Ghazaleh Dashti 
#Code last updated on 09/07/2024 
#Code further annotated on 24/12/2024

#############################VAHCS EXAMPLE###################################### 
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

# NOTE 2: 
# The code simulated data for m-DAGs "T","A","B","C","D","E","F" 
# In the manuscript these m-DAGs have been renamed as follows: 
# m-DAG T is renamed to m-DAG A 
# m-DAG A is renamed to m-DAG B 
# m-DAG F is renamed to m-DAG C 
# m-DAG B is renamed to m-DAG D 
# m-DAG D is renamed to m-DAG F 
# m-DAG C is renamed to m-DAG E  

#############################GENERAL SET UP#####################################
##Set up
rm(list=ls())
#setwd("xx")  

#Load libraries
library(boot)
library(rstream)

# Set a fixed initial seed for reproducibility
initial_seed <- 9082024
set.seed(initial_seed) 

# Initialize the rstream object with the initial seed
stream <- new("rstream.mrg32k3a")

# Set the number of seeds required
n_seeds <- 2 # 2 seeds: one for the simple (not presented in the manuscript) and one for the complex scenario

# Function to get a seed from the stream
generate_seed <- function() {
  as.integer(runif(1, 1, .Machine$integer.max))
}

# Initialize a vector to store seeds
seeds <- integer(n_seeds)

# Generate the seeds
for (i in 1:n_seeds) {
  rstream.nextsubstream(stream)
  seeds[i] <- generate_seed()
}

set.seed<- seeds[2] # set seed for complex scenario 

# fixed parameters 
sim <- 2000
n   <- 2000  

#############################DATA CHECK FUNCTIONS###############################
## Missing prop functions to check distributions and proportions with missing data across simulated datasets 
missprop<-function(dat)
{
  missingprop <- c(sum(is.na(dat$Z2)),sum(is.na(dat$Z3)),sum(is.na(dat$X)),sum(is.na(dat$L)),sum(is.na(dat$M)),sum(is.na(dat$Y)),
                   sum(is.na(dat$Z2)|is.na(dat$Z3)), 
                   sum(is.na(dat$Z2)|is.na(dat$Z3)|is.na(dat$X)),
                   sum(is.na(dat$Z2)|is.na(dat$Z3)|is.na(dat$X)|is.na(dat$L)|is.na(dat$M)), 
                   sum(is.na(dat$Z2)|is.na(dat$Z3)|is.na(dat$X)|is.na(dat$L)|is.na(dat$M)|is.na(dat$Y)))/nrow(dat)
  
  return(missingprop)
} 
## Proportion function 
prop<-function(dat)
{
  proportions <- c(prop.table(table(dat$X))[2],prop.table(table(dat$L))[2],prop.table(table(dat$M))[2],prop.table(table(dat$Y))[2],
                   prop.table(table(dat$Z1))[2],prop.table(table(dat$Z2))[2],prop.table(table(dat$Z3))[2]) 
  return(proportions)
} 

#############################CREATE DIRECTORIES################################
if(!dir.exists(paste0(getwd(),"/SIM/1_simdata"))){
  dir.create(paste0(getwd(),"/SIM/1_simdata"))
}
if(!dir.exists(paste0(getwd(),"/SIM/1_simdata/complex"))){
  dir.create(paste0(getwd(),"/SIM/1_simdata/complex"))
}
if(!dir.exists(paste0(getwd(),"/SIM/1_simdata/complex/completedata"))){
  dir.create(paste0(getwd(),"/SIM/1_simdata/complex/completedata"))
}

if(!dir.exists(paste0(getwd(),"/SIM/2_datachecks"))){
  dir.create(paste0(getwd(),"/SIM/2_datachecks"))
}
if(!dir.exists(paste0(getwd(),"/SIM/2_datachecks/complex"))){
  dir.create(paste0(getwd(),"/SIM/2_datachecks/complex"))
} 

if(!dir.exists(paste0(getwd(),"/SIM/2_datachecks/complex/missingprop"))){
  dir.create(paste0(getwd(),"/SIM/2_datachecks/complex/missingprop"))
} 
if(!dir.exists(paste0(getwd(),"/SIM/2_datachecks/complex/compdataprop"))){
  dir.create(paste0(getwd(),"/SIM/2_datachecks/complex/compdataprop"))
} 

#############################SIMULATE DATA######################################
#Load data
parms<-read.csv("parms_complex.csv")
attach(parms) 
for(i in 1:sim){
  A<-rbinom(n,1,inv.logit(a_int))
  
  Z1<-rbinom(n,1,inv.logit(z1_int))
  
  Z2<-rbinom(n,1,inv.logit(z2_int+z2_a*A))
  
  Z3<-rbinom(n,1,inv.logit(z3_int+z3_a*A))
  
  X<-rbinom(n,1,inv.logit(x_int+x_z1*Z1+x_z2*Z2+x_z3*Z3+x_z1z3*Z1*Z3))
  
  L<-rbinom(n,1,inv.logit(l_int+l_x*X+l_z1*Z1+l_z2*Z2+l_z3*Z3+l_z1z3*Z1*Z3))
  
  M<-rbinom(n,1,inv.logit(m_int+m_x*X+m_l*L+m_z1*Z1+m_z2*Z2+m_z3*Z3+m_xl*X*L+m_z1z3*Z1*Z3))
  
  Y<-rbinom(n,1,inv.logit(y_int+y_x*X+y_z1*Z1+y_z2*Z2+y_z3*Z3+y_l*L+y_m*M+y_xl*X*L+y_xm*X*M+y_z1z3*Z1*Z3))
  
  dat<-data.frame(A,Z1,Z2,Z3,X,L,M,Y)  
  
  # Save the data 
  comp_file <- paste0(getwd(),"/SIM/1_simdata/complex/completedata/data_i",i,".csv")
  write.csv(dat, file = comp_file, row.names=FALSE) 
  
  # create missingness indicators (R) for different m-DAGs and impose missing data
  for(DAG in c("T","A","B","C","D","E","F")){
    
    thedir<-paste0(getwd(),"/SIM/1_simdata/complex/data_DAG",DAG)
    if(!dir.exists(thedir)){ #create directory if it doesn't exist
      dir.create(thedir) #Create folder to store data
    }
    if (DAG=="T"){ # m-DAG A in manuscript 
      W <-rbinom(n,1,inv.logit(-1.1))
      RZ2<-rbinom(n,1,inv.logit(Trz2_int+Trz2_w*W+miss_coeff*Z1)) 
      RZ3<-rbinom(n,1,inv.logit(Trz3_int+Trz3_w*W+miss_coeff*Z1))  
      RX<-rbinom(n,1,inv.logit(Trx_int+Trx_w*W+miss_coeff*Z1))  
      RL<-rbinom(n,1,inv.logit(Trl_int+Trl_w*W+miss_coeff*Z1))
      RM<-rbinom(n,1,inv.logit(Trm_int+Trm_w*W+miss_coeff*Z1))
      RY<-rbinom(n,1,inv.logit(Try_int+Try_w*W+miss_coeff*Z1))
    } 
    if (DAG=="A"){# m-DAG B in manuscript  
      W <-rbinom(n,1,inv.logit(-1.1))
      RZ2<-rbinom(n,1,inv.logit(Arz2_int+Arz2_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*X+miss_coeff*L+miss_coeff*M)) 
      RZ3<-rbinom(n,1,inv.logit(Arz3_int+Arz3_w*W+miss_coeff*Z1+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M))  
      RX<-rbinom(n,1,inv.logit(Arx_int+Arx_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M))  
      RL<-rbinom(n,1,inv.logit(Arl_int+Arl_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M))
      RM<-rbinom(n,1,inv.logit(Arm_int+Arm_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M))
      RY<-rbinom(n,1,inv.logit(Ary_int+Ary_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M)) 
    }
    if (DAG=="B"){# m-DAG D in manuscript  
      W <-rbinom(n,1,inv.logit(-1.1))
      RZ2<-rbinom(n,1,inv.logit(Brz2_int+Brz2_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*X)) 
      RZ3<-rbinom(n,1,inv.logit(Brz3_int+Brz3_w*W+miss_coeff*Z1+miss_coeff*Z3+miss_coeff*X))  
      RX<-rbinom(n,1,inv.logit(Brx_int+Brx_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))  
      RL<-rbinom(n,1,inv.logit(Brl_int+Brl_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))
      RM<-rbinom(n,1,inv.logit(Brm_int+Brm_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))
      RY<-rbinom(n,1,inv.logit(Bry_int+Bry_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))
      
    } 
    if (DAG=="C"){# m-DAG E in manuscript  
      W <-rbinom(n,1,inv.logit(-1.1))
      RZ2<-rbinom(n,1,inv.logit(Crz2_int+Crz2_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y)) 
      RZ3<-rbinom(n,1,inv.logit(Crz3_int+Crz3_w*W+miss_coeff*Z1+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y))  
      RX<-rbinom(n,1,inv.logit(Crx_int+Crx_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y))  
      RL<-rbinom(n,1,inv.logit(Crl_int+Crl_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y))
      RM<-rbinom(n,1,inv.logit(Crm_int+Crm_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y))
      RY<-rbinom(n,1,inv.logit(Cry_int+Cry_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M)) 
    } 
    if (DAG=="D"){# m-DAG F in manuscript  
      W <-rbinom(n,1,inv.logit(-1.1))
      RZ2<-rbinom(n,1,inv.logit(Drz2_int+miss_coeff*Z1+Drz2_w*W+miss_coeff*Z2+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y)) 
      RZ3<-rbinom(n,1,inv.logit(Drz3_int+Drz3_w*W+miss_coeff*Z1+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y))  
      RX<-rbinom(n,1,inv.logit(Drx_int+Drx_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M+miss_coeff*Y))  
      RL<-rbinom(n,1,inv.logit(Drl_int+Drl_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))
      RM<-rbinom(n,1,inv.logit(Drm_int+Drm_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))
      RY<-rbinom(n,1,inv.logit(Dry_int+Dry_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))
    } 
    if (DAG=="E"){# m-DAG G in manuscript  
      W <-rbinom(n,1,inv.logit(-1.1))
      RZ2<-rbinom(n,1,inv.logit(Erz2_int+miss_coeff*Z1+Erz2_w*W+miss_coeff*Z2+miss_coeff*X)) 
      RZ3<-rbinom(n,1,inv.logit(Erz3_int+Erz3_w*W+miss_coeff*Z1+miss_coeff*Z3+miss_coeff*X))  
      RX<-rbinom(n,1,inv.logit(Erx_int+Erx_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X))  
      RL<-rbinom(n,1,inv.logit(Erl_int+Erl_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M))
      RM<-rbinom(n,1,inv.logit(Erm_int+Erm_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*L+miss_coeff*M))
      RY<-rbinom(n,1,inv.logit(Ery_int+Ery_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*Y))
    }
    if (DAG=="F"){# m-DAG C in manuscript   
      W <-rbinom(n,1,inv.logit(-1.1))
      RZ2<-rbinom(n,1,inv.logit(Frz2_int+Frz2_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*X+miss_coeff*Y)) 
      RZ3<-rbinom(n,1,inv.logit(Frz3_int+Frz3_w*W+miss_coeff*Z1+miss_coeff*Z3+miss_coeff*X+miss_coeff*Y))  
      RX<-rbinom(n,1,inv.logit(Frx_int+Frx_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*Y))  
      RL<-rbinom(n,1,inv.logit(Frl_int+Frl_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*Y))
      RM<-rbinom(n,1,inv.logit(Frm_int+Frm_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*Y))
      RY<-rbinom(n,1,inv.logit(Fry_int+Fry_w*W+miss_coeff*Z1+miss_coeff*Z2+miss_coeff*Z3+miss_coeff*X+miss_coeff*Y)) 
    }
    
    #Apply missing values according to DAG
    datDAG <- dat  
    datDAG$Z2[RZ2==1]<-NA 
    datDAG$Z3[RZ3==1]<-NA 
    datDAG$X[RX==1]<-NA 
    datDAG$L[RL==1]<-NA 
    datDAG$M[RM==1]<-NA 
    datDAG$Y[RY==1]<-NA 
    
    #Save data (with missing data)
    write.csv(datDAG, file = paste0(thedir,"/data_DAG",DAG,"_i",i,".csv"), 
              row.names=FALSE) #write data named after DAG 
    
  } 
} 

#############################DATA CHECKS########################################

### Missingness proportions
savedir <- paste0(getwd(),"/SIM/2_datachecks/complex/missingprop")

RPM <- data.frame()

for(DAG in c("T","A","B","C","D","E","F")){
  rpmDAG <- data.frame()
  
  for(i in 1:sim){
    file <- paste0(getwd(),"/SIM/1_simdata/complex/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
    missing_dat <- read.csv(file)
    
    rpmDAG <- rbind(rpmDAG, missprop(dat = missing_dat))
  }
  RPM<-rbind(RPM,colMeans(rpmDAG))
} 

names(RPM)<-c("Z2","Z3","X","L","M","Y","Zs","X,Z","Z,X,M","Any")
RPM <- cbind(DAG=c("T","A","B","C","D","E","F"), RPM)
write.csv(RPM, paste0(savedir,"/missingprop.csv"), row.names=FALSE)

### proportions in complete data 
savedir <- paste0(getwd(),"/SIM/2_datachecks/complex/compdataprop")
props <- data.frame()
PROPS <- data.frame() 

for(i in 1:sim){
  file <- paste0(getwd(),"/SIM/1_simdata/complex/completedata/data_i",i,".csv")
  complete_dat <- read.csv(file)
  
  props <- rbind(props, prop(dat = complete_dat))
}
PROPS<-rbind(PROPS,colMeans(props)) 
names(PROPS)<-c("X","L","M","Y","Z1","Z2","Z3")
write.csv(PROPS, paste0(savedir,"/compdataprop.csv"), row.names=FALSE)