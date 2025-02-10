#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name: intmed_func_complex
#Purpose: This is R script for functions for intventional mediation analysis:
#           noMYMIBoot (noLZYMIBoot in the manuscript) 
#           noMYBootMI (noLZYBootMI in the manuscript) 
#           noYMIBoot 
#           noYBootMI
#           nointMIBoot 
#           nointBootMI
#           McompMIBoot (LZintMIBoot in the manuscript)
#           McompBootMI (LZintBootMI in the manuscript)
#           YcompMIBoot (YintMIBoot in the manuscript)
#           YcompBootMI (YintBootMI in the manuscript)
#           higherintMIBoot (allintMIBoot in the manuscript)
#           higherintBootMI (allintBootMI in the manuscript)
#           smcfcsMIBoot  
#           smcfcsBootMI  
#Author:    Ghazaleh Dashti 
#Final revision on 08/07/2024 
#Code further annotated on 24/12/2024  
# All models are correctly specified (i.e, same as data generating models)

#############################VAHCS EXAMPLE###################################### 
# Simulation study based on an example from VAHCS 
# exposure X: Depression and anxiety at least 2 waves during adolescence 
#              (as a proxy for severity of depression and anxiety)

# intermediate confounder L: any live birth by age 24 

# mediator M: Depression and anxiety at least 2 wave in young adulthood 
#              (as a proxy for severity of depression and anxiety)

# outcome Y: depression and anxiety in adulthood 

# baseline confounders Zs: Z1: sex registered at birth 
#                          Z2: cannabis use 
#                          Z3: adolescent (wave 2-6) antisocial behaviour 

# Auxiliary variable A: parental smoking 

#-------------------------------------------------------------------------------
###########int-based mediation analysis with bootstrap for CI####
intmed_fun<-function(dat)
{
  intmedest_fun<-function(dat, ind=1:nrow(dat))
  { 
    mcsim<-50
    # take bootstrap sample 
    data <- dat[ind,] 
    
    #Set flag to capture bootstrap samples to reject
    flag<-FALSE 
    
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
    
    if(!flag)return(res) else
      return(rep(NA,length(res)))  
  }
  
  #Run mediation with the bootstrap
  bstrap_intmed<-boot(data=dat, statistic=intmedest_fun,stype="i", R=nboot,parallel="multicore",ncpus=inner_cores)
  
  #Recover bootstrap results
  res<-data.frame()
  for(i in 1:length(bstrap_intmed$t0))
  { 
    est<-(bstrap_intmed$t0[i]) 
    se<-sd(bstrap_intmed$t[,i],na.rm=T)
    bt<-boot.ci(bstrap_intmed,index=i,type=c("perc"))  
    cilow<-(bt$percent[4]) 
    ciupp<-(bt$percent[5])
    pval<-2*(1-pnorm(q=abs((est)/se)))
    res<-rbind(res,c(est,se,cilow,ciupp,pval))
  }
  names(res)<-c("Estimate","SE","95CIlow","95CIupp","pvalue")
  row.names(res)<-c("DE_RD", "IE_RD")
  intmed_out<-as.data.frame(res)
}  

#-------------------------------------------------------------------------------
# wrapper function to analyse an imputed dataset (without bootstrap - this is used with the boot mi functions)####
intmedbootimpute_fun <-function(dat) 
{ 
  data <- dat 
  mcsim<-50
  
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
} 

#-------------------------------------------------------------------------------
# Function to pool MI boot results##############################################  
pool_MIBoot_fun <- function(results_mi)
{
  
  DERD_df <- results_mi[seq(1, nrow(results_mi), 8), ] 
  IERD_df <- results_mi[seq(2, nrow(results_mi), 8), ]  
  
   #### 
  DERD_est <- mean(DERD_df$Estimate)
  var_inter <- (1/(m-1))*sum((DERD_df$Estimate-DERD_est)^2)
  var_intra <- mean((DERD_df$SE)^2)
  DERD_var <- var_intra + (1 + 1/m)*var_inter
  r<-(1+1/m)*var_inter/var_intra
  vvv<-(m-1)*(1+1/r)^2 
  tal<-qt(0.025,df=vvv,lower.tail=F) 
  DERD_lCI <- DERD_est - tal*sqrt(DERD_var)
  DERD_uCI <- DERD_est + tal*sqrt(DERD_var) 
  DERD_pval<-2*(1-pt(q=abs(DERD_est/sqrt(DERD_var)),df=vvv))  
  DERD_se<-sqrt(DERD_var)   
  
  #### 
  IERD_est <- mean(IERD_df$Estimate)
  var_inter <- (1/(m-1))*sum((IERD_df$Estimate-IERD_est)^2)
  var_intra <- mean((IERD_df$SE)^2)
  IERD_var <- var_intra + (1 + 1/m)*var_inter
  r<-(1+1/m)*var_inter/var_intra
  vvv<-(m-1)*(1+1/r)^2 
  tal<-qt(0.025,df=vvv,lower.tail=F) 
  IERD_lCI <- IERD_est - tal*sqrt(IERD_var)
  IERD_uCI <- IERD_est + tal*sqrt(IERD_var) 
  IERD_pval<-2*(1-pt(q=abs(IERD_est/sqrt(IERD_var)),df=vvv))  
  IERD_se<-sqrt(IERD_var)    
  
  DERD<- cbind(DERD_est,DERD_se,DERD_lCI,DERD_uCI,DERD_pval)
  IERD<- cbind(IERD_est,IERD_se,IERD_lCI,IERD_uCI,IERD_pval)  
  
  res<-data.frame(rbind(DERD,IERD)) 
  names(res)<-c("Estimate","SE","95CIlow","95CIupp","pvalue")
  row.names(res)<-c("DERD","IERD")
  res<-as.data.frame(res) 
}

#-------------------------------------------------------------------------------
#Function to get p values for bootMI############################################ 
pvals_BootMI_fun <- function(ests) 
{ 
  df <- as.data.frame(ests)
  
  DERD_pval <- 2 * stats::pt(abs(df[1,1]/df[1,2]^0.5), df = df[1,5], 
                             lower.tail = FALSE)  
  
  IERD_pval <- 2 * stats::pt(abs(df[2,1]/df[2,2]^0.5), df = df[2,5], 
                             lower.tail = FALSE)  
  
  df$pval <- c(DERD_pval, IERD_pval)
  df<-df[-c(5)]
  names(df)<-c("Estimate","SE","95CIlow","95CIupp","pvalue")
  row.names(df)<-c("DE_RD", "IE_RD") 
  res_bootMI <-as.data.frame(df)  
} 

#-------------------------------------------------------------------------------
# noMYMIBOOT####################################################################   
int_noMYMIBoot_fun <- function(dat,setseed) {
  
  # Start the outer parallel cluster
  outer_cl <- makeCluster(outer_cores, type = "FORK")
  registerDoParallel(outer_cl) 
  registerDoRNG(seed = setseed)  
  
  dat<-dat
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[,"Y"] <-0 
  predMat[,"M"] <-0 
  predMat[,"L"] <-0 
  predMat[c("X","L","M","Y"),"A"] <-0    
  
  ###MI 
  MI<-mice::mice(dat,m=m,method=c(" "," ","logreg","logreg","logreg","logreg","logreg","logreg"),
                 predictorMatrix=predMat,maxit=5,printFlag=F)
  
  impdatasets <- vector("list", m)    
  for (j in 1:m) {
    impdatasets[[j]] <- mice::complete(MI, j)
  }
  
  # Save the MI and impdatasets in a list
  result_list <- list(MI = MI, impdatasets = impdatasets, outer_cl = outer_cl)
  
  # Define the function for parallel processing 
  parallel_processing <- function(j) {
    impdataset <- as.data.table(impdatasets[[j]])
    out <- intmed_fun(dat = impdataset)
    out
  }
  
  # Perform the parallel processing
  results_mi <- foreach(j = 1:m, .packages = c("mice", "boot", "parallel"), .combine = rbind) %dopar% {
    # Call the parallel processing function and pass impdatasets
    parallel_processing(j)
  } 
  
  # Stop the outer parallel cluster
  stopCluster(outer_cl) 
  
  res_MIboot <- pool_MIBoot_fun(results_mi) 
}

#-------------------------------------------------------------------------------
# noMYBootMI####################################################################   
# wrapper function to call mice generating M imputations 
impM_noMY <- function(dat,M) {
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[,"Y"] <-0 
  predMat[,"M"] <-0 
  predMat[,"L"] <-0     
  predMat[c("X","L","M","Y"),"A"] <-0      
  
  miceImps<-mice::mice(dat,m=M,method=c(" "," ","logreg","logreg","logreg","logreg","logreg","logreg"),
                       predictorMatrix=predMat,maxit=5,printFlag=F) 
  implist <- vector("list",M)
  for (j in 1:M) {
    implist[[j]] <- mice::complete(miceImps,j)
  }
  implist
}  

int_noMYBootMI_fun <- function(dat,setseed) {  
  impute <- bootImpute(dat, impM_noMY, nBoot=nboot, nImp=M, M=M,nCores=ncores, seed=setseed) 
  ests <- bootImputeAnalyse(impute, intmedbootimpute_fun,nCores=ncores)  
  
  df <- as.data.frame(ests)
  res_bootMI<-pvals_BootMI_fun(df)   
}


#-------------------------------------------------------------------------------
# noYMIBOOT#####################################################################   
int_noYMIBoot_fun <- function(dat,setseed) {
  
  # Start the outer parallel cluster
  outer_cl <- makeCluster(outer_cores, type = "FORK")
  registerDoParallel(outer_cl) 
  registerDoRNG(seed = setseed)  
  
  dat<-dat
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[,"Y"] <-0    
  predMat[c("X","L","M","Y"),"A"] <-0    
  
  ###MI 
  MI<-mice::mice(dat,m=m,method=c(" "," ","logreg","logreg","logreg","logreg","logreg","logreg"),
                 predictorMatrix=predMat,maxit=5,printFlag=F)
  
  impdatasets <- vector("list", m)    
  for (j in 1:m) {
    impdatasets[[j]] <- mice::complete(MI, j)
  }
  
  # Save the MI and impdatasets in a list
  result_list <- list(MI = MI, impdatasets = impdatasets, outer_cl = outer_cl)
  
  # Define the function for parallel processing 
  parallel_processing <- function(j) {
    impdataset <- as.data.table(impdatasets[[j]])
    out <- intmed_fun(dat = impdataset)
    out
  }
  
  # Perform the parallel processing
  results_mi <- foreach(j = 1:m, .packages = c("mice", "boot", "parallel"), .combine = rbind) %dopar% {
    # Call the parallel processing function and pass impdatasets
    parallel_processing(j)
  } 
  
  # Stop the outer parallel cluster
  stopCluster(outer_cl) 
  
  res_MIboot <- pool_MIBoot_fun(results_mi) 
}

#-------------------------------------------------------------------------------
# noYBootMI#####################################################################   
# wrapper function to call mice generating M imputations 
impM_noY <- function(dat,M) {
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[,"Y"] <-0    
  predMat[c("X","L","M","Y"),"A"] <-0      
  
  miceImps<-mice::mice(dat,m=M,method=c(" "," ","logreg","logreg","logreg","logreg","logreg","logreg"),
                       predictorMatrix=predMat,maxit=5,printFlag=F) 
  implist <- vector("list",M)
  for (j in 1:M) {
    implist[[j]] <- mice::complete(miceImps,j)
  }
  implist
}  

int_noYBootMI_fun <- function(dat,setseed) {  
  impute <- bootImpute(dat, impM_noY, nBoot=nboot, nImp=M, M=M,nCores=ncores, seed=setseed) 
  ests <- bootImputeAnalyse(impute, intmedbootimpute_fun,nCores=ncores)  
  
  df <- as.data.frame(ests)
  res_bootMI<-pvals_BootMI_fun(df)   
}

#-------------------------------------------------------------------------------
# nointMIBOOT####################################################################    
int_nointMIBoot_fun <- function(dat,setseed) {
  
  # Start the outer parallel cluster
  outer_cl <- makeCluster(outer_cores, type = "FORK")
  registerDoParallel(outer_cl) 
  registerDoRNG(seed = setseed)  
  
  dat<-dat
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[c("X","L","M","Y"),"A"] <-0    
  
  ###MI 
  MI<-mice::mice(dat,m=m,method=c(" "," ","logreg","logreg","logreg","logreg","logreg","logreg"),
                 predictorMatrix=predMat,maxit=5,printFlag=F)
  
  impdatasets <- vector("list", m)    
  for (j in 1:m) {
    impdatasets[[j]] <- mice::complete(MI, j)
  }
  
  # Save the MI and impdatasets in a list
  result_list <- list(MI = MI, impdatasets = impdatasets, outer_cl = outer_cl)
  
  # Define the function for parallel processing 
  parallel_processing <- function(j) {
    impdataset <- as.data.table(impdatasets[[j]])
    out <- intmed_fun(dat = impdataset)
    out
  }
  
  # Perform the parallel processing
  results_mi <- foreach(j = 1:m, .packages = c("mice", "boot", "parallel"), .combine = rbind) %dopar% {
    # Call the parallel processing function and pass impdatasets
    parallel_processing(j)
  } 
  
  # Stop the outer parallel cluster
  stopCluster(outer_cl) 
  
  res_MIboot <- pool_MIBoot_fun(results_mi) 
}

#-------------------------------------------------------------------------------
# nointBootMI####################################################################   
# wrapper function to call mice generating M imputations 
impM_noint <- function(dat,M) {
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[c("X","L","M","Y"),"A"] <-0   
  
  miceImps<-mice::mice(dat,m=M,method=c(" "," ","logreg","logreg","logreg","logreg","logreg","logreg"),
                       predictorMatrix=predMat,maxit=5,printFlag=F) 
  implist <- vector("list",M)
  for (j in 1:M) {
    implist[[j]] <- mice::complete(miceImps,j)
  }
  implist
}  

int_nointBootMI_fun <- function(dat,setseed) {  
  impute <- bootImpute(dat, impM_noint, nBoot=nboot, nImp=M, M=M,nCores=ncores, seed=setseed) 
  ests <- bootImputeAnalyse(impute, intmedbootimpute_fun,nCores=ncores)  
  
  df <- as.data.frame(ests)
  res_bootMI<-pvals_BootMI_fun(df)   
}

#-------------------------------------------------------------------------------
# McompMIBoot####################################################################     
int_McompMIBoot_fun <- function(dat,setseed) {
  
  # Start the outer parallel cluster
  outer_cl <- makeCluster(outer_cores, type = "FORK")
  registerDoParallel(outer_cl) 
  registerDoRNG(seed = setseed)  
  
  dat<-dat
  
  # add relevant interactions 
  dat$LZ1 <- dat$L*dat$Z1 
  dat$LZ3 <- dat$L*dat$Z3 
  dat$MZ1 <- dat$M*dat$Z1 
  dat$MZ3 <- dat$M*dat$Z3 
  dat$XL <- dat$X*dat$L  
  dat$XM <- dat$X*dat$M 
  dat$Z1Z3 <- dat$Z1*dat$Z3
  
  # modify the matrix 
  predMat <- mice::make.predictorMatrix(dat)
  predMat[c("X"),c("XL","XM")] <-0 
  predMat[c("L"),c("LZ1","LZ3","XL")] <-0  
  predMat[c("M"),c("MZ1","MZ3","XM")] <-0   
  predMat["Z1",c("MZ1","Z1Z3")]<-0 
  predMat["Z3",c("MZ3","Z1Z3")]<-0 
  predMat[c("X","L","M","Y"),"A"] <-0       
  
  ### specify the imputation methods 
  impMethod <- mice::make.method(dat)
  impMethod[c("Z2","Z3","X","L","M","Y")]<-"logreg"
  impMethod["LZ1"] <- "~I(L*Z1)" 
  impMethod["LZ3"] <- "~I(L*Z3)" 
  impMethod["MZ1"] <- "~I(M*Z1)" 
  impMethod["MZ3"] <- "~I(M*Z3)" 
  impMethod["XL"] <- "~I(X*L)" 
  impMethod["XM"] <- "~I(X*M)" 
  impMethod["Z1Z3"] <- "~I(Z1*Z3)"   
  
  ###MI 
  MI<-mice::mice(dat,m=m,method=impMethod,
                 predictorMatrix=predMat,maxit=5,printFlag=F)
  
  impdatasets <- vector("list", m)    
  for (j in 1:m) {
    impdatasets[[j]] <- mice::complete(MI, j)
  }
  
  # Save the MI and impdatasets in a list
  result_list <- list(MI = MI, impdatasets = impdatasets, outer_cl = outer_cl)
  
  # Define the function for parallel processing 
  parallel_processing <- function(j) {
    impdataset <- as.data.table(impdatasets[[j]])
    out <- intmed_fun(dat = impdataset)
    out
  }
  
  # Perform the parallel processing
  results_mi <- foreach(j = 1:m, .packages = c("mice", "boot", "parallel"), .combine = rbind) %dopar% {
    # Call the parallel processing function and pass impdatasets
    parallel_processing(j)
  } 
  
  # Stop the outer parallel cluster
  stopCluster(outer_cl) 
  
  res_MIboot <- pool_MIBoot_fun(results_mi) 
}

#-------------------------------------------------------------------------------
# McompBootMI####################################################################       
# wrapper function to call mice generating M imputations 
impM_mcomp <- function(dat,M) {
  
  # add relevant interactions 
  dat$LZ1 <- dat$L*dat$Z1 
  dat$LZ3 <- dat$L*dat$Z3 
  dat$MZ1 <- dat$M*dat$Z1 
  dat$MZ3 <- dat$M*dat$Z3 
  dat$XL <- dat$X*dat$L  
  dat$XM <- dat$X*dat$M 
  dat$Z1Z3 <- dat$Z1*dat$Z3
  
  # modify the matrix 
  predMat <- mice::make.predictorMatrix(dat)
  predMat[c("X"),c("XL","XM")] <-0 
  predMat[c("L"),c("LZ1","LZ3","XL")] <-0  
  predMat[c("M"),c("MZ1","MZ3","XM")] <-0   
  predMat["Z1",c("MZ1","Z1Z3")]<-0 
  predMat["Z3",c("MZ3","Z1Z3")]<-0 
  predMat[c("X","L","M","Y"),"A"] <-0       
  
  ### specify the imputation methods 
  impMethod <- mice::make.method(dat)
  impMethod[c("Z2","Z3","X","L","M","Y")]<-"logreg"
  impMethod["LZ1"] <- "~I(L*Z1)" 
  impMethod["LZ3"] <- "~I(L*Z3)" 
  impMethod["MZ1"] <- "~I(M*Z1)" 
  impMethod["MZ3"] <- "~I(M*Z3)" 
  impMethod["XL"] <- "~I(X*L)" 
  impMethod["XM"] <- "~I(X*M)" 
  
  impMethod["Z1Z3"] <- "~I(Z1*Z3)"  
  
  #MI
  miceImps<- mice::mice(dat,m=M,method=impMethod,predictorMatrix=predMat,maxit=5,printFlag=F) 
  implist <- vector("list",M) 
  
  for (j in 1:M) {
    implist[[j]] <- mice::complete(miceImps,j)
  }
  implist
}  

int_McompBootMI_fun <- function(dat,setseed) {  
  impute <- bootImpute(dat, impM_mcomp, nBoot=nboot, nImp=M, M=M,nCores=ncores, seed=setseed) 
  ests <- bootImputeAnalyse(impute, intmedbootimpute_fun,nCores=ncores)  
  
  df <- as.data.frame(ests)
  res_bootMI<-pvals_BootMI_fun(df)   
}

#-------------------------------------------------------------------------------
# ycompMIBoot####################################################################    
int_YcompMIBoot_fun <- function(dat,setseed) {
  
  # Start the outer parallel cluster
  outer_cl <- makeCluster(outer_cores, type = "FORK")
  registerDoParallel(outer_cl) 
  registerDoRNG(seed = setseed)  
  
  dat<-dat
  
  # add relevant interactions 
  dat$YX <- dat$Y*dat$X 
  dat$YM <- dat$Y*dat$M 
  dat$YL <- dat$Y*dat$L  
  dat$XL   <- dat$X*dat$L
  dat$XM   <- dat$X*dat$M
  dat$LM   <- dat$L*dat$M 
  
  dat$YZ1 <- dat$Y*dat$Z1 
  dat$YZ3 <- dat$Y*dat$Z3 

  dat$Z1Z3 <- dat$Z1*dat$Z3
  
  # modify the matrix 
  predMat <- mice::make.predictorMatrix(dat)
  predMat[c("Y"),c("YX","YL","YM","YZ1","YZ3")] <-0 
  predMat[c("X"),c("YX","XL","XM")] <-0
  predMat[c("L"),c("YL","XL")] <-0   
  predMat[c("M"),c("YM","XM")] <-0   
  predMat["Z1",c("YZ1","Z1Z3")]<-0 
  predMat["Z3",c("YZ3","Z1Z3")]<-0 
  predMat[c("X","L","M","Y"),"A"] <-0       
  
  ### specify the imputation methods 
  impMethod <- mice::make.method(dat)
  impMethod[c("Z2","Z3","X","L","M","Y")]<-"logreg"
  
  impMethod["YX"] <- "~I(Y*X)" 
  impMethod["YM"] <- "~I(Y*M)" 
  impMethod["YL"] <- "~I(Y*L)"
  
  impMethod["XL"] <- "~I(X*L)" 
  impMethod["XM"] <- "~I(X*M)" 
   
  impMethod["YZ1"] <- "~I(Y*Z1)" 
   
  impMethod["YZ3"] <- "~I(Y*Z3)" 

  impMethod["Z1Z3"] <- "~I(Z1*Z3)"   
  
  ###MI 
  MI<-mice::mice(dat,m=m,method=impMethod,
                 predictorMatrix=predMat,maxit=5,printFlag=F) 
  
  impdatasets <- vector("list", m)    
  for (j in 1:m) {
    impdatasets[[j]] <- mice::complete(MI, j)
  }
  
  # Save the MI and impdatasets in a list
  result_list <- list(MI = MI, impdatasets = impdatasets, outer_cl = outer_cl)
  
  # Define the function for parallel processing 
  parallel_processing <- function(j) {
    impdataset <- as.data.table(impdatasets[[j]])
    out <- intmed_fun(dat = impdataset)
    out
  }
  
  # Perform the parallel processing
  results_mi <- foreach(j = 1:m, .packages = c("mice", "boot", "parallel"), .combine = rbind) %dopar% {
    # Call the parallel processing function and pass impdatasets
    parallel_processing(j)
  } 
  
  # Stop the outer parallel cluster
  stopCluster(outer_cl) 
  
  res_MIboot <- pool_MIBoot_fun(results_mi) 
}

#-------------------------------------------------------------------------------
# ycompBootMI####################################################################       
# wrapper function to call mice generating M imputations 
impM_ycomp <- function(dat,M) {
  
  # add relevant interactions 
  dat$YX <- dat$Y*dat$X 
  dat$YM <- dat$Y*dat$M 
  dat$YL <- dat$Y*dat$L  
  dat$XL   <- dat$X*dat$L
  dat$XM   <- dat$X*dat$M
  dat$LM   <- dat$L*dat$M 
  
  dat$YZ1 <- dat$Y*dat$Z1 
  dat$YZ3 <- dat$Y*dat$Z3 
  
  dat$Z1Z3 <- dat$Z1*dat$Z3
  
  # modify the matrix 
  predMat <- mice::make.predictorMatrix(dat)
  predMat[c("Y"),c("YX","YL","YM","YZ1","YZ3")] <-0 
  predMat[c("X"),c("YX","XL","XM")] <-0
  predMat[c("L"),c("YL","XL")] <-0   
  predMat[c("M"),c("YM","XM")] <-0   
  predMat["Z1",c("YZ1","Z1Z3")]<-0 
  predMat["Z3",c("YZ3","Z1Z3")]<-0 
  predMat[c("X","L","M","Y"),"A"] <-0   
  
  ### specify the imputation methods 
  impMethod <- mice::make.method(dat)
  impMethod[c("Z2","Z3","X","L","M","Y")]<-"logreg"
  
  impMethod["YX"] <- "~I(Y*X)" 
  impMethod["YM"] <- "~I(Y*M)" 
  impMethod["YL"] <- "~I(Y*L)"
  impMethod["XL"] <- "~I(X*L)" 
  impMethod["XM"] <- "~I(X*M)" 

  impMethod["YZ1"] <- "~I(Y*Z1)" 
  impMethod["YZ3"] <- "~I(Y*Z3)" 
  
  impMethod["Z1Z3"] <- "~I(Z1*Z3)" 
  
  #MI
  miceImps<- mice::mice(dat,m=M,method=impMethod,predictorMatrix=predMat,maxit=5,printFlag=F) 
  implist <- vector("list",M) 
  
  for (j in 1:M) {
    implist[[j]] <- mice::complete(miceImps,j)
  }
  implist
}  

int_YcompBootMI_fun <- function(dat,setseed) {  
  impute <- bootImpute(dat, impM_ycomp, nBoot=nboot, nImp=M, M=M,nCores=ncores, seed=setseed) 
  ests <- bootImputeAnalyse(impute, intmedbootimpute_fun,nCores=ncores)  
  
  df <- as.data.frame(ests)
  res_bootMI<-pvals_BootMI_fun(df)   
}

#-------------------------------------------------------------------------------
# higherintMIBoot####################################################################    
int_higherintMIBoot_fun <- function(dat,setseed) {
  
  # Start the outer parallel cluster
  outer_cl <- makeCluster(outer_cores, type = "FORK")
  registerDoParallel(outer_cl) 
  registerDoRNG(seed = setseed)  
  
  dat<-dat
  
  # add relevant interactions 
  dat$LZ1 <- dat$L*dat$Z1 
  dat$LZ3 <- dat$L*dat$Z3 
  dat$MZ1 <- dat$M*dat$Z1 
  dat$MZ3 <- dat$M*dat$Z3 
  dat$YZ1 <- dat$Y*dat$Z1 
  dat$YZ3 <- dat$Y*dat$Z3 
  
  dat$YX <- dat$Y*dat$X 
  dat$YM <- dat$Y*dat$M 
  dat$YL <- dat$Y*dat$L  
  
  dat$XL <- dat$X*dat$L  
  dat$XM <- dat$X*dat$M 

  dat$Z1Z3 <- dat$Z1*dat$Z3
  
  # modify the matrix 
  predMat <- mice::make.predictorMatrix(dat)
  predMat[c("Y"),c("YX","YM","YL","YZ1","YZ3")] <-0 
  predMat[c("X"),c("YX","XL","XM")] <-0   
  predMat[c("L"),c("YL","LZ1","LZ3","XL")] <-0  
  predMat[c("M"),c("YM","MZ1","MZ3","XM")] <-0   
  predMat["Z1",c("YZ1","LZ1","MZ1","Z1Z3")]<-0 
  predMat["Z3",c("YZ3","LZ3","MZ3","Z1Z3")]<-0 
  predMat[c("X","L","M","Y"),"A"] <-0       
  
  ### specify the imputation methods 
  impMethod <- mice::make.method(dat)
  impMethod[c("Z2","Z3","X","L","M","Y")]<-"logreg"
  
  impMethod["LZ1"] <- "~I(L*Z1)" 
  impMethod["LZ3"] <- "~I(L*Z3)" 
  impMethod["MZ1"] <- "~I(M*Z1)" 
  impMethod["MZ3"] <- "~I(M*Z3)" 
  impMethod["YZ1"] <- "~I(Y*Z1)" 
  impMethod["YZ3"] <- "~I(Y*Z3)" 
  
  impMethod["YX"] <- "~I(Y*X)" 
  impMethod["YM"] <- "~I(Y*M)" 
  impMethod["YL"] <- "~I(Y*L)"
  
  impMethod["XL"] <- "~I(X*L)"
  impMethod["XM"] <- "~I(X*M)"

  impMethod["Z1Z3"] <- "~I(Z1*Z3)"   
  
  ###MI 
  MI<-mice::mice(dat,m=m,method=impMethod,
                 predictorMatrix=predMat,maxit=5,printFlag=F) 
  
  impdatasets <- vector("list", m)    
  for (j in 1:m) {
    impdatasets[[j]] <- mice::complete(MI, j)
  }
  
  # Save the MI and impdatasets in a list
  result_list <- list(MI = MI, impdatasets = impdatasets, outer_cl = outer_cl)
  
  # Define the function for parallel processing 
  parallel_processing <- function(j) {
    impdataset <- as.data.table(impdatasets[[j]])
    out <- intmed_fun(dat = impdataset)
    out
  }
  
  # Perform the parallel processing
  results_mi <- foreach(j = 1:m, .packages = c("mice", "boot", "parallel"), .combine = rbind) %dopar% {
    # Call the parallel processing function and pass impdatasets
    parallel_processing(j)
  } 
  
  # Stop the outer parallel cluster
  stopCluster(outer_cl) 
  
  res_MIboot <- pool_MIBoot_fun(results_mi) 
}

#-------------------------------------------------------------------------------
# higherintBootMI###############################################################     
# wrapper function to call mice generating M imputations 
impM_higherint <- function(dat,M) {
  
  # add relevant interactions 
  dat$LZ1 <- dat$L*dat$Z1 
  dat$LZ3 <- dat$L*dat$Z3 
  dat$MZ1 <- dat$M*dat$Z1 
  dat$MZ3 <- dat$M*dat$Z3 
  dat$YZ1 <- dat$Y*dat$Z1 
  dat$YZ3 <- dat$Y*dat$Z3 
  dat$YX <- dat$Y*dat$X 
  dat$YM <- dat$Y*dat$M 
  dat$YL <- dat$Y*dat$L  
  dat$XL <- dat$X*dat$L  
  dat$XM <- dat$X*dat$M 
  dat$Z1Z3 <- dat$Z1*dat$Z3
  
  # modify the matrix 
  predMat <- mice::make.predictorMatrix(dat)
  predMat[c("Y"),c("YX","YM","YL","YZ1","YZ3")] <-0 
  predMat[c("X"),c("YX","XL","XM")] <-0   
  predMat[c("L"),c("YL","LZ1","LZ3","XL")] <-0  
  predMat[c("M"),c("YM","MZ1","MZ3","XM")] <-0   
  predMat["Z1",c("YZ1","LZ1","MZ1","Z1Z3")]<-0 
  predMat["Z3",c("YZ3","LZ3","MZ3","Z1Z3")]<-0 
  predMat[c("X","L","M","Y"),"A"] <-0    
  
  ### specify the imputation methods 
  impMethod <- mice::make.method(dat)
  impMethod[c("Z2","Z3","X","L","M","Y")]<-"logreg"
  
  impMethod["LZ1"] <- "~I(L*Z1)" 
  impMethod["LZ3"] <- "~I(L*Z3)" 
  
  impMethod["MZ1"] <- "~I(M*Z1)" 
  impMethod["MZ3"] <- "~I(M*Z3)" 
  
  impMethod["YZ1"] <- "~I(Y*Z1)" 
  impMethod["YZ3"] <- "~I(Y*Z3)" 
  
  impMethod["YX"] <- "~I(Y*X)" 
  impMethod["YM"] <- "~I(Y*M)" 
  impMethod["YL"] <- "~I(Y*L)"
  
  impMethod["XL"] <- "~I(X*L)"
  impMethod["XM"] <- "~I(X*M)"
  
  impMethod["Z1Z3"] <- "~I(Z1*Z3)"   
  #MI
  miceImps<- mice::mice(dat,m=M,method=impMethod,predictorMatrix=predMat,maxit=5,printFlag=F) 
  implist <- vector("list",M) 
  
  for (j in 1:M) {
    implist[[j]] <- mice::complete(miceImps,j)
  }
  implist
}  

int_higherintBootMI_fun <- function(dat,setseed) {  
  impute <- bootImpute(dat, impM_higherint, nBoot=nboot, nImp=M, M=M,nCores=ncores, seed=setseed) 
  ests <- bootImputeAnalyse(impute, intmedbootimpute_fun,nCores=ncores)  
  
  df <- as.data.frame(ests)
  res_bootMI<-pvals_BootMI_fun(df)   
}

#-------------------------------------------------------------------------------
# smcfcsMIBoot####################################################################  
int_smcfcsMIBoot_fun <- function(dat,setseed) {
  
  # Start the outer parallel cluster
  outer_cl <- makeCluster(outer_cores, type = "FORK")
  registerDoParallel(outer_cl) 
  registerDoRNG(seed = setseed)  
  
  dat<-data.frame(dat) #making sure the data frame is in correct format for smcfcs 
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[c("X","L","M","Y"),"A"] <-0     
  
  ###MI 
  MI<- smcfcs(dat, smtype="logistic",
              smformula="Y~X+L+M+Z1+Z2+Z3+X*M+X*L+Z1*Z3",
              method=c("","","brlogreg","brlogreg","brlogreg","brlogreg","brlogreg",""), predictorMatrix=predMat, m=m)
  
  impdatasets <- vector("list", m)    
  for (j in 1:m) {
    impdatasets[[j]] <- MI$impDatasets[j]
  }
  
  # Save the MI and impdatasets in a list
  result_list <- list(MI = MI, impdatasets = impdatasets, outer_cl = outer_cl)
  
  # Define the function for parallel processing 
  parallel_processing <- function(j) {
    impdataset <- as.data.table(impdatasets[[j]])
    out <- intmed_fun(dat = impdataset)
    out
  }
  
  # Perform the parallel processing
  results_mi <- foreach(j = 1:m, .packages = c("mice", "boot", "parallel"), .combine = rbind) %dopar% {
    # Call the parallel processing function and pass impdatasets
    parallel_processing(j)
  } 
  
  # Stop the outer parallel cluster
  stopCluster(outer_cl) 
  
  res_MIboot <- pool_MIBoot_fun(results_mi) 
}

#-------------------------------------------------------------------------------
# smcfcsBootMI####################################################################     
int_smcfcsBootMI_fun<- function(dat,setseed) { 
  
  dat<-data.frame(dat) #making sure the data frame is in correct format for smcfcs 
  
  predMat <- mice::make.predictorMatrix(dat) 
  predMat[c("X","L","M","Y"),"A"] <-0      
  
  impute <- bootSmcfcs(dat, nBoot=nboot, nImp=M,predictorMatrix=predMat, 
                       smtype="brlogistic", smformula="Y~X+L+M+Z1+Z2+Z3+X*M+X*L+Z1*Z3",
                       method=c("","","brlogreg","brlogreg","brlogreg","brlogreg","brlogreg",""),nCores=ncores, seed=setseed) 
  
  
  
  ests <- bootImputeAnalyse(impute, intmedbootimpute_fun)  
  df <- as.data.frame(ests)
  res_bootMI<-pvals_BootMI_fun(df)   
}   

#-------------------------------------------------------------------------------
