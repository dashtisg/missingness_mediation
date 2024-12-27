#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#file name: 0503_final_res 
#Purpose: This is the R script for combining pooled performance measures across all causal diagrams 
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


#Load arguments
#args <- commandArgs(trailingOnly=TRUE)  
#method <- args[1] #

#############################CREATE DIRECTORIES AND POOL ALL RESULTS############
if(!dir.exists(paste0(getwd(),"/SIM/7_finalres"))){
  dir.create(paste0(getwd(),"/SIM/7_finalres"))
}

pool_file <- multmerge(paste0(getwd(),"/SIM/6_poolperfmsr/simple")) 
pooled_dat <- paste0(getwd(),"/SIM/7_finalres/finalres_simple.csv")
write.csv(pool_file, file=pooled_dat, row.names=FALSE)

pool_file <- multmerge(paste0(getwd(),"/SIM/6_poolperfmsr/complex")) 
pooled_dat <- paste0(getwd(),"/SIM/7_finalres/finalres_complex.csv")
write.csv(pool_file, file=pooled_dat, row.names=FALSE)

