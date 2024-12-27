#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name: pool_func  
#Purpose:	This is R script for functions for pooling results 
#Author:Ghazaleh Dashti 
#Final revision on 14/06/2023    
#Code further annotated on 24/12/2024   

#-------------------------------------------------------------------------------
#POOL DATASET FUNCTIONS######################################################### 
multmerge_DERD <- function(mypath) {
  # Get list of files
  filenames <- list.files(path = mypath, full.names = TRUE)
  
  # Function to process each file
  process_file <- function(x) {
    if (file.info(x)$size == 0) {
      return(NULL)  # Skip empty files
    }
    
    tryCatch({
      # Read the file
      data <- read.csv(file = x, header = TRUE)
      
      # Check if the data frame has the required number of rows and columns
      if (nrow(data) > 0 && ncol(data) >= 5) {
        return(data[1, 1:5])
      } else {
        warning(paste("File", x, "does not have enough data or columns. Skipping."))
        return(NULL)
      }
    }, error = function(e) {
      # Handle read errors
      warning(paste("Error reading file", x, ":", e$message))
      return(NULL)
    })
  }
  
  # Process all files and filter out NULL results
  datalist <- Filter(Negate(is.null), lapply(filenames, process_file))
  
  if (length(datalist) == 0) {
    stop("No valid data to merge.")
  }
  
  # Merge remaining valid data
  Reduce(function(x, y) merge(x, y, all = TRUE), datalist)
}



multmerge_IERD<- function(mypath){
  # Get list of files
  filenames <- list.files(path = mypath, full.names = TRUE)
  
  # Function to process each file
  process_file <- function(x) {
    if (file.info(x)$size == 0) {
      return(NULL)  # Skip empty files
    }
    
    tryCatch({
      # Read the file
      data <- read.csv(file = x, header = TRUE)
      
      # Check if the data frame has the required number of rows and columns
      if (nrow(data) > 0 && ncol(data) >= 5) {
        return(data[2, 1:5])
      } else {
        warning(paste("File", x, "does not have enough data or columns. Skipping."))
        return(NULL)
      }
    }, error = function(e) {
      # Handle read errors
      warning(paste("Error reading file", x, ":", e$message))
      return(NULL)
    })
  }
  
  # Process all files and filter out NULL results
  datalist <- Filter(Negate(is.null), lapply(filenames, process_file))
  
  if (length(datalist) == 0) {
    stop("No valid data to merge.")
  }
  
  # Merge remaining valid data
  Reduce(function(x, y) merge(x, y, all = TRUE), datalist)
}

multmerge <- function(mypath){
  filenames <- list.files(path=mypath, full.names=TRUE)
  datalist <- lapply(filenames, function(x){read.csv(file=x,header=T,)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
} 
#-------------------------------------------------------------------------------
#PERFORMANCE MEASURES ##########################################################
sim_para<-function(dat){
mean_est  <- mean(pooled_dat$Estimate)  
# Absolute bias 
bias_abs  <- mean_est - true
# Absolute bias MC SE 
bias_abs_se  <- sqrt((sum((pooled_dat$Estimate - mean_est)^2))/(nsim*(nsim-1)))
# Relative bias 
bias_rel <- 100*(bias_abs/true) 
# Relative bias MC SE 
bias_rel_se <- 100*(1/true)*bias_abs_se
# Empirical SE 
emp_se <- sqrt((sum((pooled_dat$Estimate- mean_est)^2))/(nsim-1)) 
# Empirical SE MC SE 
emp_se_se  <- emp_se/sqrt(2*(nsim-1))
# Average Model SE 
mod_se <- sqrt(mean((pooled_dat$SE)^2)) 
# Average Model SE MC SE 
mod_se_se  <- sqrt((var((((pooled_dat$SE)^2)) - (mod_se)))/(4*nsim*(mod_se)^2)) 
# Relative % error in Model SE 
relmodse_err <- 100*((mod_se/emp_se) - 1) 
# Relative % error in Model SE MC SE 
relmodse_err_se  <- 100*(mod_se/emp_se)*sqrt(((var((((pooled_dat$SE)^2)) - (mod_se)))/(4*nsim*(mod_se)^4)) + (1/(2*(nsim-1)))) 
# Coverage 
int <- cbind(pooled_dat$X95CIlow, pooled_dat$X95CIupp) 
cov<-mean(apply(int, 1, findInterval, x = true) == 1)
# Coverage MC SE 
cov_se  <- sqrt((cov*(1-cov))/nsim)
# Bias Eliminated Coverage 
cov_bc<-mean(apply(int, 1, findInterval, x = mean_est) == 1)
# Bias Eliminated Coverage MC SE 
cov_bc_se  <- sqrt((cov_bc*(1-cov_bc))/nsim) 

fails<- (sum(pooled_dat$fail)/nsim)*100

sim_parameters <-data.frame(mean_est, 
                            bias_abs, bias_abs_se, 
                            bias_rel, bias_rel_se, 
                            cov, cov_se, 
                            cov_bc, cov_bc_se, 
                            emp_se, emp_se_se,
                            mod_se, mod_se_se, 
                            relmodse_err,relmodse_err_se,fails)
sim_parameters<-head(sim_parameters, n = 1)   
sim_parameters <-as.data.frame(sim_parameters) 
} 

