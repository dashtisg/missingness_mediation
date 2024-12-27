#############################PROJECT############################################ 
# Handling multivariable missing data in causal mediation analysis using interventional effects

#############################R SCRIPT DESCRIPTION###############################
#File name: 000_seeds 
#Purpose: This is the R script for creating seeds for running the analyses using 
# Different missing data methods on the HPC 
#Author: Ghazaleh Dashti 
#Code last updated on 13/08/2024  
#Code further annotated on 24/12/2024

#############################IMPORTANT NOTES#################################### 
# NOTE 1: 
# The seeds are created for two scenarios (methods), referred to as "simple" and "complex" 
# In the final version of the manuscript only results from the complex scenario are presented.  
# The difference between the two scenarios is in the presence of an intermediate confounder: 
# In the "simple" scenario there are no intermediate confounders  
# In the "complex" scenario any live birth by age 24 is considered as an intermediate confounder 


# NOTE 2: 
# The seeds are created for m-DAGs"T","A","B","C","D","E","F" 
# In the manuscript these m-DAGs have been renamed as follows: 
# m-DAG T is renamed to m-DAG A 
# m-DAG A is renamed to m-DAG B 
# m-DAG F is renamed to m-DAG C 
# m-DAG B is renamed to m-DAG D 
# m-DAG D is renamed to m-DAG F 
# m-DAG C is renamed to m-DAG E 

#############################GENERATE SEEDS##################################### 
# load the rstream package 
library(rstream)
rm(list=ls()) 

# set working directory 
#setwd("xx")  

miss<-c("CC",
        "noMYBootMI", 
        "noYBootMI", 
        "nointBootMI", 
        "McompBootMI",
        "YcompBootMI",
        "higherintBootMI", 
        "smcfcsBootMI",
        "noMYMIBoot",
        "noYMIBoot", 
        "nointMIBoot",
        "McompMIBoot",
        "YcompMIBoot", 
        "higherintMIBoot",
        "smcfcsMIBoot")

method<-c("simple","complex")

DAG <- c("T","A","B","C","D","E","F")

simrun<- (c(1:2000))

nseeds<-length(miss)*length(method)*length(DAG)*length(simrun)

initial_seed <- 12345

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

dat<-cbind(method=rep(method, each=nseeds/length(method)),
           miss=rep(rep(miss, each=nseeds/(length(method)*length(miss))),length(method)),
           DAG=rep(rep(DAG, each=nseeds/(length(method)*length(miss)*length(DAG))),length(miss)), 
           simrun=rep(simrun, (length(method)*length(miss)*length(DAG))),seed=seeds)
dat<-as.data.frame(dat)
write.csv(dat,"seeds.csv",row.names=F, quote=FALSE)
