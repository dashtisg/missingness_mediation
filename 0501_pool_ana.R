#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name: 0501_pool_ana
#Purpose: This is the R script for calculating performance measures 
#Author: Ghazaleh Dashti 
#Last updated on 19/08/2024
#Code further annotated on 24/12/2024

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
rm(list=ls()) 
#setwd("xx")

# install and load related packages
packages <- c("parallel","foreach", "data.table", "dplyr", "readstata13", "tidyverse")

for (package in packages) {
  library(package, character.only=T)
} 

# source the necessary functions 
source("func_pool.R")   

#############################CREATE DIRECTORIES#################################
if(!dir.exists(paste0(getwd(),"/SIM/5_perfmsr"))){
  dir.create(paste0(getwd(),"/SIM/5_perfmsr"))
}

#############################PERFORMANCE MEASURES############################### 
methods <- c("complex")  
misses <- c("CC") 
DAGs <- c("T", "A", "B", "C", "D", "E", "F")
 
for (method in methods) {
  for (miss in misses) {
    for (DAG in DAGs) { 
      
      if(!dir.exists(paste0(getwd(),"/SIM/5_perfmsr/",method))){
        dir.create(paste0(getwd(),"/SIM/5_perfmsr/",method))
      }
      
      if(!dir.exists(paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG))){
        dir.create(paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG))
      }  
      
      # Load the true estimates
      if (method == "simple") {
        true_est <- read.csv(paste0(getwd(), "/SIM/4_true/simple/true_simp.csv"))
      } else if (method == "complex") {
        true_est <- read.csv(paste0(getwd(), "/SIM/4_true/complex/true_complex.csv"))
      } 
      
      causal_diagram<-paste(DAG) 
      analysis<-paste(miss)   
      
      pooled_dat <- multmerge_DERD(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG)) 
      nsim<- nrow(pooled_dat)
      true<- (true_est[1,1]) 
      DERD_sim_para<- sim_para(dat=pooled_dat) 
      DERD_sim_para<-data.frame(causal_diagram,analysis,DERD_sim_para,nsim)
      
      pooled_dat <- multmerge_IERD(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG)) 
      nsim<- nrow(pooled_dat)
      true<- (true_est[1,2]) 
      IERD_sim_para<- sim_para(dat=pooled_dat) 
      IERD_sim_para<-data.frame(causal_diagram,analysis,IERD_sim_para,nsim)
      
      perf_measures<- rbind(DERD_sim_para, IERD_sim_para)
      
      effect<- data.frame(Effect = c("DE_RD","IE_RD"))
      perf_measures<-cbind(effect, perf_measures)
      
      perf_meas <- paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG,"/perfmsr_",method,"_",miss,"_",DAG,".csv")
      write.csv(perf_measures, file=perf_meas, row.names=FALSE) 
      
    }
  }
}

### FOR BOOTMI methods the SE is actually the variance so need to take the sqrt first

methods <- c("complex")  

misses <- c("noMYBootMI",
            "noYBootMI",
            "nointBootMI",
            "McompBootMI",
            "YcompBootMI",
            "higherintBootMI",
            "smcfcsBootMI")

DAGs <- c("T", "A", "B", "C", "D", "E", "F")

for (method in methods) {
  for (miss in misses) {
    for (DAG in DAGs) { 

if(!dir.exists(paste0(getwd(),"/SIM/5_perfmsr/",method))){
        dir.create(paste0(getwd(),"/SIM/5_perfmsr/",method))
 }
      
if(!dir.exists(paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG))){
      dir.create(paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG))
    }  
  
# Load the true estimates
 if (method == "simple") {
        true_est <- read.csv(paste0(getwd(), "/SIM/4_true/simple/true_simp.csv"))
      } else if (method == "complex") {
        true_est <- read.csv(paste0(getwd(), "/SIM/4_true/complex/true_complex.csv"))
      } 
      
    causal_diagram<-paste(DAG) 
    analysis<-paste(miss)   
    
    pooled_dat <- multmerge_DERD(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG)) 
    pooled_dat[,2]<-sqrt(pooled_dat[,2])  
    nsim<- nrow(pooled_dat)
    true<- (true_est[1,1]) 
    DERD_sim_para<- sim_para(dat=pooled_dat) 
    DERD_sim_para<-data.frame(causal_diagram,analysis,DERD_sim_para,nsim)
    
    pooled_dat <- multmerge_IERD(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG)) 
    pooled_dat[,2]<-sqrt(pooled_dat[,2]) 
    nsim<- nrow(pooled_dat)
    true<- (true_est[1,2]) 
    IERD_sim_para<- sim_para(dat=pooled_dat) 
    IERD_sim_para<-data.frame(causal_diagram,analysis,IERD_sim_para,nsim)
    
    perf_measures<- rbind(DERD_sim_para, IERD_sim_para)
    
    effect<- data.frame(Effect = c("DE_RD","IE_RD"))
    perf_measures<-cbind(effect, perf_measures)
    
    perf_meas <- paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG,"/perfmsr_",method,"_",miss,"_",DAG,".csv")
    write.csv(perf_measures, file=perf_meas, row.names=FALSE) 
    
    }
  }
}

### 
methods <- c("complex")  

misses <- c("noMYMIBoot",
            "noYMIBoot",
            "nointMIBoot",
            "McompMIBoot",
            "YcompMIBoot",
            "higherintMIBoot",
            "smcfcsMIBoot")

DAGs <- c("T", "A", "B", "C", "D", "E", "F")

for (method in methods) {
  for (miss in misses) {
    for (DAG in DAGs) { 
      
      if(!dir.exists(paste0(getwd(),"/SIM/5_perfmsr/",method))){
        dir.create(paste0(getwd(),"/SIM/5_perfmsr/",method))
      }
      
      if(!dir.exists(paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG))){
        dir.create(paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG))
      }  
      
      # Load the true estimates
      if (method == "simple") {
        true_est <- read.csv(paste0(getwd(), "/SIM/4_true/simple/true_simp.csv"))
      } else if (method == "complex") {
        true_est <- read.csv(paste0(getwd(), "/SIM/4_true/complex/true_complex.csv"))
      } 
      
      causal_diagram<-paste(DAG) 
      analysis<-paste(miss)   
      
      pooled_dat <- multmerge_DERD(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG)) 
      nsim<- nrow(pooled_dat)
      true<- (true_est[1,1]) 
      DERD_sim_para<- sim_para(dat=pooled_dat) 
      DERD_sim_para<-data.frame(causal_diagram,analysis,DERD_sim_para,nsim)
      
      pooled_dat <- multmerge_IERD(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG)) 
      nsim<- nrow(pooled_dat)
      true<- (true_est[1,2]) 
      IERD_sim_para<- sim_para(dat=pooled_dat) 
      IERD_sim_para<-data.frame(causal_diagram,analysis,IERD_sim_para,nsim)
      
      perf_measures<- rbind(DERD_sim_para, IERD_sim_para)
      
      effect<- data.frame(Effect = c("DE_RD","IE_RD"))
      perf_measures<-cbind(effect, perf_measures)
      
      perf_meas <- paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG,"/perfmsr_",method,"_",miss,"_",DAG,".csv")
      write.csv(perf_measures, file=perf_meas, row.names=FALSE) 
      
    }
  }
}
