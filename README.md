# missingness_mediation
Simulation code for the paper "Handling multivariable missing data in causal mediation analysis using interventional effects"

# The simulation study is based on a motivating example from a previously published investigation using data from the Victorian Adolescent Health Cohort Study (VAHCS) - DOI: 10.1192/bjp.2022.3 
    
# Important notes: 
# NOTE 1: 
    # The general set-up of the codes is for two scenarios, simple and complex. However, only the codes for the complex scenario are uploaded to the repository and only results from the complex scenario are presented in the manuscript. 
    # The difference between the two scenarios is in the presence of an intermediate confounder: 
    # In the "simple" scenario there are no intermediate confounders.  
    # In the "complex" scenario any live birth by age 24 is considered as an intermediate confounder. 

# NOTE 2: 
    # In the simulation codes, different missingness mechanisms are simulated under seven m-DAGs titled "T","A","B","C","D","E","F" 
    # In the manuscript these m-DAGs have been renamed as follows: 
    # m-DAG T is renamed to m-DAG A 
    # m-DAG A is renamed to m-DAG B 
    # m-DAG F is renamed to m-DAG C 
    # m-DAG B is renamed to m-DAG D 
    # m-DAG D is renamed to m-DAG F 
    # m-DAG C is renamed to m-DAG E  

# NOTE 3: 
The parameter values used for simulating the data were determined by fitting similar models to the available VAHCS data. These parameter values are available from the file "parms" uploaded to the repository. 


    
