#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name: 0300_sim_ana
#Purpose: This is the R script for the analysis of data using different missing data methods 
#Author: Ghazaleh Dashti 
#Last updated on 13/08/2024
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
setwd("/group/cebu1/BACKUP/Ghazaleh/simstudy_med_revision")

#load related packages
packages <- c("parallel","foreach", "boot", "mice","devtools","data.table", "dplyr", "readstata13", "tidyverse", "zoo", "bootImpute","smcfcs","mitools","plyr","micemd","doParallel","doRNG")

for (package in packages) {
  library(package, character.only=T)
}

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
outer_cores<- 2
ncores<- 5
inner_cores<- 5

#Load arguments
args <- commandArgs(trailingOnly=TRUE)  

method<-args[1] #simple, complex

miss<-args[2]

if(miss=="CC"|
   miss=="noMYMIBoot"|miss=="noMYBootMI"| 
   miss=="noYMIBoot"|miss=="noYBootMI"|
   miss=="nointMIBoot"|miss=="nointBootMI"|
   miss=="McompMIBoot"|miss=="McompBootMI"|
   miss=="YcompMIBoot"|miss=="YcompBootMI"|
   miss=="higherintMIBoot"|miss=="higherintBootMI"|
   miss=="smcfcsMIBoot"|miss=="smcfcsBootMI"){ 
  DAG <- args[3] 
} 

batch <- as.numeric(args[4])   # 1:50
# generate batch_i_values 
total_batches <- 50
iterations_per_batch <- 40

# Function to generate i values for each batch
generate_i_values <- function(batch_number, iterations_per_batch) {
  start_i <- (batch_number - 1) * iterations_per_batch + 1
  end_i <- batch_number * iterations_per_batch
  return(start_i:end_i)
}

# Create the batch_i_values list
batch_i_values <- lapply(1:total_batches, generate_i_values, iterations_per_batch = iterations_per_batch)

#############################CREARE DIRECTORIES#################################
if(!dir.exists(paste0(getwd(),"/SIM/3_simres"))){
  dir.create(paste0(getwd(),"/SIM/3_simres"))
}

if(!dir.exists(paste0(getwd(),"/SIM/3_simres/",method))){
  dir.create(paste0(getwd(),"/SIM/3_simres/",method))
}

if(!dir.exists(paste0(getwd(),"/SIM/3_simres/",method,"/",miss))){
  dir.create(paste0(getwd(),"/SIM/3_simres/",method,"/",miss))
} 

if(miss=="CC"|
   miss=="noMYMIBoot"|miss=="noMYBootMI"| 
   miss=="noYMIBoot"|miss=="noYBootMI"|
   miss=="nointMIBoot"|miss=="nointBootMI"|
   miss=="McompMIBoot"|miss=="McompBootMI"|
   miss=="YcompMIBoot"|miss=="YcompBootMI"|
   miss=="higherintMIBoot"|miss=="higherintBootMI"|
   miss=="smcfcsMIBoot"|miss=="smcfcsBootMI"){ 
  if(!dir.exists(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG))){
    dir.create(paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG))
  }  
} 


#Print at beginning 
if(miss=="CC"|
   miss=="noMYMIBoot"|miss=="noMYBootMI"| 
   miss=="noYMIBoot"|miss=="noYBootMI"|
   miss=="nointMIBoot"|miss=="nointBootMI"|
   miss=="McompMIBoot"|miss=="McompBootMI"|
   miss=="YcompMIBoot"|miss=="YcompBootMI"|
   miss=="higherintMIBoot"|miss=="higherintBootMI"|
   miss=="smcfcsMIBoot"|miss=="smcfcsBootMI"){ 
  print(paste0("This is running: ", method, "using ", miss , "under DAG ", DAG , ", under batch ", batch))
} 

#############################DATA ANALYSIS#####################################
if (!is.null(batch_i_values[[batch]])) {
  if (miss=="CC") {
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      # create complete data 
      dat_cc<-na.omit(dat) 
      # run analysis
      RNGkind("L'Ecuyer-CMRG")
      set.seed(setseed)
      int_CC <- intmed_fun(dat = dat_cc)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_CC, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="noMYMIBoot") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      seed<-setseed
      # run analysis
      int_noMYMIBoot <- int_noMYMIBoot_fun(dat = dat, setseed = setseed) 
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_noMYMIBoot, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="noMYBootMI") {
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- as.numeric(seed_dat$seed)
      print(setseed) 
      
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_noMYBootMI <- int_noMYBootMI_fun(dat = dat, setseed = setseed)  
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_noMYBootMI, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="noYMIBoot") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_noYMIBoot <- int_noYMIBoot_fun(dat = dat, setseed = setseed)  
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_noYMIBoot, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="noYBootMI") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_noYBootMI <- int_noYBootMI_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_noYBootMI, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="nointMIBoot") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_nointMIBoot <- int_nointMIBoot_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_nointMIBoot, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="nointBootMI") {
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_nointBootMI <- int_nointBootMI_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_nointBootMI, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  } 
  else if (miss=="McompMIBoot") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_McompMIBoot <- int_McompMIBoot_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_McompMIBoot, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="McompBootMI") {
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_McompBootMI <- int_McompBootMI_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_McompBootMI, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  } 
  else if (miss=="YcompMIBoot") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_YcompMIBoot <- int_YcompMIBoot_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_YcompMIBoot, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="YcompBootMI") {
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_YcompBootMI <- int_YcompBootMI_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_YcompBootMI, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="higherintMIBoot") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_higherintMIBoot <- int_higherintMIBoot_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_higherintMIBoot, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="higherintBootMI") {
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_higherintBootMI <- int_higherintBootMI_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_higherintBootMI, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  } 
  else if (miss=="smcfcsMIBoot") {
    
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # load the data
      thefile <- paste0(getwd(),"/SIM/1_simdata/",method,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,".csv")
      dat <- read.csv(thefile)
      
      # run analysis
      int_smcfcsMIBoot <- int_smcfcsMIBoot_fun(dat = dat, setseed = setseed)
      
      # save the results
      comp_file <- paste0(getwd(),"/SIM/3_simres/",method,"/",miss,"/data_DAG",DAG,"/data_DAG",DAG,"_i",i,"_results.csv")
      write.csv(int_smcfcsMIBoot, file = comp_file, row.names = FALSE)
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }
  else if (miss=="smcfcsBootMI") {
    # Function to run the analysis and save the results
    run_int_analysis <- function(i) {
      # get the seed 
      thefile <- paste0(getwd(),"/seeds.csv")
      seeds <- read.csv(thefile)
      meth<- paste0(method)
      missing<-paste0(miss)
      DAG_name<-paste0(DAG)
      sim<-paste0(i)
      seed_dat<-seeds %>% filter(method==meth,miss==missing,DAG==DAG_name,simrun==sim)  
      setseed <- seed_dat$seed
      print(setseed) 
      # Define a flag to track whether an error occurred
      error_occurred <- FALSE
      
      # load the data
      thefile <- paste0(getwd(), "/SIM/1_simdata/",method,"/data_DAG", DAG, "/data_DAG", DAG, "_i", i, ".csv")
      dat <- read.csv(thefile)
      
      # Use tryCatch to run the analysis and handle errors
      int_smcfcsBootMI <- tryCatch({
        # run analysis
        int_smcfcsBootMI_fun(dat = dat, setseed = setseed)
      }, error = function(err) {
        # Error handling code
        cat("Error occurred while processing data_i", i, ":", conditionMessage(err), "\n")
        error_occurred <- TRUE
        # You can choose to return or save an error object here if needed
        return(NULL)
      })
      
      # Check if an error occurred and only save results if no error
      if (!error_occurred) {
        # save the results
        comp_file <- paste0(getwd(), "/SIM/3_simres/", method, "/", miss, "/data_DAG", DAG, "/data_DAG", DAG, "_i", i, "_results.csv")
        write.csv(int_smcfcsBootMI, file = comp_file, row.names = FALSE)
      }
    }
    # Get the appropriate set of i values for the batch
    i_values <- batch_i_values[[batch]]
    # Use lapply to loop over i values and perform the analysis for each i
    lapply(i_values, run_int_analysis)
  }  
}
time <- Sys.time() - start_time ; print(time) 

