#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name: 0502_pool_prfmsr
#Purpose: This is the R script for pooling performing measures 
#Author: Ghazaleh Dashti 
#Last updated on 19/08/2024
#Code further annotated on 24/12/2024 

#############################GENERAL SET UP#####################################
rm(list=ls()) 
#setwd("xx") 

# install and load related packages
packages <- c("foreach", "boot", "devtools")

for (package in packages) {
  library(package, character.only=T)
} 

# source the necessary functions 
source("func_pool.R")  

methods <- c("complex")   
DAGs <- c("T", "A", "B", "C", "D", "E", "F")
 
#############################CREATE DIRECTORIES AND POOL PERFORMANCE MEASURES####
if(!dir.exists(paste0(getwd(),"/SIM/6_poolperfmsr"))){
  dir.create(paste0(getwd(),"/SIM/6_poolperfmsr"))
} 

for (method in methods) {
    for (DAG in DAGs) {  
      
if(!dir.exists(paste0(getwd(),"/SIM/6_poolperfmsr/",method))){
   dir.create(paste0(getwd(),"/SIM/6_poolperfmsr/",method))
      }        

      pool_file <- multmerge(paste0(getwd(),"/SIM/5_perfmsr/",method,"/data_DAG",DAG))  
pooled_dat <- paste0(getwd(),"/SIM/6_poolperfmsr/",method,"/pool_perfmsr_DAG",DAG,".csv")
write.csv(pool_file, file=pooled_dat, row.names=FALSE)

    } 
} 
