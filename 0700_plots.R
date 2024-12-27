#############################PROJECT############################################
# Handling multivariable missing data in causal mediation analysis using interventional effects 

#############################R SCRIPT DESCRIPTION############################### 
#File name: 0700_plot
#Purpose: This is R script for creating the plots 
#Author: Ghazaleh Dashti 
#Last updated: 19/08/2024
#Code further annotated on 24/12/2024 

#############################IMPORTANT NOTES#################################### 
# The general set up of the codes is for two scenarios, simple and complex. 
# This code simulates data under complex scenario. 
# In the final version of the manuscript only results from the complex scenario are presented.
# The difference between the two scenarios is in presence of an intermediate confounder: 
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

# THE m-DAGS ARE RENAMED IN THIS SCRIPT  

#############################GENERAL SET UP#####################################
#load related packages
packages <- c("magrittr","plyr", "ggpubr", "ggthemes","ggsci","ggplot2", "readr", "gridExtra", "cowplot", "dplyr", "patchwork")

for (package in packages) {
  library(package, character.only=T)
}

rm(list=ls())
#setwd("xx") 

res_complex<- read_csv("SIM/7_finalres/finalres_complex.csv") 

res_complex$causal_diagram[res_complex$causal_diagram=="TRUE"]<- "T" # rename true to T 
res_complex$causal_diagram[res_complex$causal_diagram=="FALSE"]<- "F" # rename true to F  

#############################RENAME m-DAGS BASED ON WHAT IS PRESENTED IN THE MANUSCRIPT#### 
res_complex$causal_diagram[res_complex$causal_diagram=="E"]<- "G" 
res_complex$causal_diagram[res_complex$causal_diagram=="C"]<- "E" 
res_complex$causal_diagram[res_complex$causal_diagram=="F"]<- "C" 
res_complex$causal_diagram[res_complex$causal_diagram=="D"]<- "F" 
res_complex$causal_diagram[res_complex$causal_diagram=="B"]<- "D" 
res_complex$causal_diagram[res_complex$causal_diagram=="A"]<- "B" 
res_complex$causal_diagram[res_complex$causal_diagram=="T"]<- "A" 
res_complex$causal_diagram <- factor(res_complex$causal_diagram, 
                                 levels = c("A","B","C","D","E","F","G"))

#############################RENAME MISSING METHODS BASED ON WHAT IS PRESENTED IN THE MANUSCRIPT#### 
res_complex$analysis[res_complex$analysis=="smcfcsBootMI"]<- "MI-SMCFCS (BootMI)"  
res_complex$analysis[res_complex$analysis=="smcfcsMIBoot"]<- "MI-SMCFCS (MIBoot)"  
res_complex$analysis[res_complex$analysis=="higherintBootMI"]<- "MI-allint (BootMI)"  
res_complex$analysis[res_complex$analysis=="higherintMIBoot"]<- "MI-allint (MIBoot)"  
res_complex$analysis[res_complex$analysis=="YcompBootMI"]<- "MI-Yint (BootMI)"  
res_complex$analysis[res_complex$analysis=="YcompMIBoot"]<- "MI-Yint (MIBoot)"  
res_complex$analysis[res_complex$analysis=="McompBootMI"]<- "MI-LZint (BootMI)"  
res_complex$analysis[res_complex$analysis=="McompMIBoot"]<- "MI-LZint (MIBoot)"  
res_complex$analysis[res_complex$analysis=="nointBootMI"]<- "MI-noint (BootMI)"  
res_complex$analysis[res_complex$analysis=="nointMIBoot"]<- "MI-noint (MIBoot)"  
res_complex$analysis[res_complex$analysis=="noYBootMI"]<- "MI-noY (BootMI)"  
res_complex$analysis[res_complex$analysis=="noYMIBoot"]<- "MI-noY (MIBoot)"  
res_complex$analysis[res_complex$analysis=="noMYBootMI"]<- "MI-noLZY (BootMI)"  
res_complex$analysis[res_complex$analysis=="noMYMIBoot"]<- "MI-noLZY (MIBoot)"  
res_complex$analysis[res_complex$analysis=="CC"]<- "CCA" 

#############################PREPARE FOR PLOTTING###############################
### order the methods 
methods_order_int <- c("MI-SMCFCS (BootMI)",
                       "MI-SMCFCS (MIBoot)",
                       "MI-allint (BootMI)",
                       "MI-allint (MIBoot)",
                       "MI-Yint (BootMI)", 
                       "MI-Yint (MIBoot)", 
                       "MI-LZint (BootMI)",
                       "MI-LZint (MIBoot)",
                       "MI-noint (BootMI)",
                       "MI-noint (MIBoot)", 
                       "MI-noY (BootMI)",
                       "MI-noY (MIBoot)", 
                       "MI-noLZY (BootMI)",
                       "MI-noLZY (MIBoot)", 
                       "CCA")  



res_complex$analysis <- factor(res_complex$analysis, levels = methods_order_int)

res_complex$cov<- (res_complex$cov)*100
res_complex$cov_bc<- (res_complex$cov_bc)*100
res_complex$emp_se<- (res_complex$emp_se)*100
res_complex$emp_se_se<- (res_complex$emp_se_se)*100

### save the reordered results 
res_complex<-ddply(res_complex, c("causal_diagram", "analysis")) 

res_complex_derd <- res_complex[which(res_complex$Effect=='DE_RD'),] 
path <- paste0(getwd(),"/SIM/8_plots/res_complex_derd.csv")
write.csv(res_complex_derd, file=path, row.names=FALSE)   

res_complex_ierd <- res_complex[which(res_complex$Effect=='IE_RD'),] 
path <- paste0(getwd(),"/SIM/8_plots/res_complex_ierd.csv")
write.csv(res_complex_ierd, file=path, row.names=FALSE)    

miss_method <- c("MI-SMCFCS (BootMI)",
                 "MI-SMCFCS (MIBoot)",
                 "MI-allint (BootMI)",
                 "MI-allint (MIBoot)",
                 "MI-Yint (BootMI)", 
                 "MI-Yint (MIBoot)", 
                 "MI-LZint (BootMI)",
                 "MI-LZint (MIBoot)",
                 "MI-noint (BootMI)",
                 "MI-noint (MIBoot)", 
                 "MI-noY (BootMI)",
                 "MI-noY (MIBoot)", 
                 "MI-noLZY (BootMI)",
                 "MI-noLZY (MIBoot)", 
                 "CCA")    


res_complex_ierd$analysis <- factor(res_complex_ierd$analysis, levels = paste0(miss_method)) 
res_complex_derd$analysis <- factor(res_complex_derd$analysis, levels = paste0(miss_method))

#############################RELATIVE BIAS######################################
# create dataset to only include 1 of the two MI approaches 
res_complex_ierd_relbias<- res_complex_ierd[(res_complex_ierd$analysis=="CCA" | 
                                       res_complex_ierd$analysis=="MI-noLZY (BootMI)" |
                                       res_complex_ierd$analysis=="MI-noY (BootMI)" |
                                       res_complex_ierd$analysis=="MI-noint (BootMI)" | 
                                       res_complex_ierd$analysis=="MI-LZint (BootMI)" |
                                       res_complex_ierd$analysis=="MI-Yint (BootMI)" |
                                       res_complex_ierd$analysis=="MI-allint (BootMI)" |
                                       res_complex_ierd$analysis=="MI-SMCFCS (BootMI)"),]

res_complex_derd_relbias<- res_complex_derd[(res_complex_derd$analysis=="CCA" | 
                                       res_complex_derd$analysis=="MI-noLZY (BootMI)" |
                                       res_complex_derd$analysis=="MI-noY (BootMI)" |
                                       res_complex_derd$analysis=="MI-noint (BootMI)" | 
                                       res_complex_derd$analysis=="MI-LZint (BootMI)" |
                                       res_complex_derd$analysis=="MI-Yint (BootMI)" |
                                       res_complex_derd$analysis=="MI-allint (BootMI)" |
                                       res_complex_derd$analysis=="MI-SMCFCS (BootMI)"),]


### separate the DAGS 
res_complex_ierd_relbiasA <- res_complex_ierd_relbias[res_complex_ierd_relbias$causal_diagram == "A", ] 
res_complex_ierd_relbiasB <- res_complex_ierd_relbias[res_complex_ierd_relbias$causal_diagram == "B", ] 
res_complex_ierd_relbiasC <- res_complex_ierd_relbias[res_complex_ierd_relbias$causal_diagram == "C", ] 
res_complex_ierd_relbiasD <- res_complex_ierd_relbias[res_complex_ierd_relbias$causal_diagram == "D", ] 
res_complex_ierd_relbiasE <- res_complex_ierd_relbias[res_complex_ierd_relbias$causal_diagram == "E", ] 
res_complex_ierd_relbiasF <- res_complex_ierd_relbias[res_complex_ierd_relbias$causal_diagram == "F", ] 
res_complex_ierd_relbiasG <- res_complex_ierd_relbias[res_complex_ierd_relbias$causal_diagram == "G", ] 


res_complex_derd_relbiasA <- res_complex_derd_relbias[res_complex_derd_relbias$causal_diagram == "A", ] 
res_complex_derd_relbiasB <- res_complex_derd_relbias[res_complex_derd_relbias$causal_diagram == "B", ] 
res_complex_derd_relbiasC <- res_complex_derd_relbias[res_complex_derd_relbias$causal_diagram == "C", ] 
res_complex_derd_relbiasD <- res_complex_derd_relbias[res_complex_derd_relbias$causal_diagram == "D", ] 
res_complex_derd_relbiasE <- res_complex_derd_relbias[res_complex_derd_relbias$causal_diagram == "E", ] 
res_complex_derd_relbiasF <- res_complex_derd_relbias[res_complex_derd_relbias$causal_diagram == "F", ] 
res_complex_derd_relbiasG <- res_complex_derd_relbias[res_complex_derd_relbias$causal_diagram == "G", ] 


dev.off()
int_iebias_A <- ggplot(res_complex_ierd_relbiasA, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_iebias_B <- ggplot(res_complex_ierd_relbiasB, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_iebias_C <- ggplot(res_complex_ierd_relbiasC, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_iebias_D <- ggplot(res_complex_ierd_relbiasD, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_iebias_E <- ggplot(res_complex_ierd_relbiasE, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_iebias_F <- ggplot(res_complex_ierd_relbiasF, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))


int_iebias_G <- ggplot(res_complex_ierd_relbiasG, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

# Add margin adjustments to each plot individually
int_iebias_A <- int_iebias_A + theme(plot.margin = margin(0, 0, 0, 0))
int_iebias_B <- int_iebias_B + theme(plot.margin = margin(0, 0, 0, 0))
int_iebias_C <- int_iebias_C + theme(plot.margin = margin(0, 0, 0, 0))
int_iebias_D <- int_iebias_D + theme(plot.margin = margin(0, 0, 0, 0))
int_iebias_E <- int_iebias_E + theme(plot.margin = margin(0, 0, 0, 0))
int_iebias_F <- int_iebias_F + theme(plot.margin = margin(0, 0, 0, 0))
int_iebias_G <- int_iebias_G + theme(plot.margin = margin(0, 0, 0, 0))

# Use ggarrange to arrange the plots with reduced space
all_bias_ie <- ggarrange(
  int_iebias_A, int_iebias_B, int_iebias_C,
  int_iebias_D, int_iebias_E, int_iebias_F,
  int_iebias_G, ncol = 2, nrow = 4,
  common.legend = TRUE, legend = "bottom",
  align = "hv",  # Align plots horizontally and vertically
  widths = c(1, 1, 1),  # Relative widths of columns
  heights = c(1, 1, 1),  # Relative heights of rows
  padding = 0  # Reduce space between plots
)

# Print the arranged plots
all_bias_ie


dev.off()
int_debias_A <- ggplot(res_complex_derd_relbiasA, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_debias_B <- ggplot(res_complex_derd_relbiasB, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_debias_C <- ggplot(res_complex_derd_relbiasC, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_debias_D <- ggplot(res_complex_derd_relbiasD, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_debias_E <- ggplot(res_complex_derd_relbiasE, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_debias_F <- ggplot(res_complex_derd_relbiasF, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))

int_debias_G <- ggplot(res_complex_derd_relbiasG, aes(bias_rel,analysis)) + 
  geom_errorbar(
    aes(xmin = bias_rel-bias_rel_se, xmax = bias_rel+bias_rel_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(-85, 40)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +  
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10))


# Add margin adjustments to each plot individually
int_debias_A <- int_debias_A + theme(plot.margin = margin(0, 0, 0, 0))
int_debias_B <- int_debias_B + theme(plot.margin = margin(0, 0, 0, 0))
int_debias_C <- int_debias_C + theme(plot.margin = margin(0, 0, 0, 0))
int_debias_D <- int_debias_D + theme(plot.margin = margin(0, 0, 0, 0))
int_debias_E <- int_debias_E + theme(plot.margin = margin(0, 0, 0, 0))
int_debias_F <- int_debias_F + theme(plot.margin = margin(0, 0, 0, 0))
int_debias_G <- int_debias_G + theme(plot.margin = margin(0, 0, 0, 0))

# Use ggarrange to arrange the plots with reduced space
all_bias_de <- ggarrange(
  int_debias_A, int_debias_B, int_debias_C,
  int_debias_D, int_debias_E, int_debias_F,
  int_debias_G, ncol = 2, nrow = 4,
  common.legend = TRUE, legend = "bottom",
  align = "hv",  # Align plots horizontally and vertically
  widths = c(1, 1, 1),  # Relative widths of columns
  heights = c(1, 1, 1),  # Relative heights of rows
  padding = 0  # Reduce space between plots
)

# Print the arranged plots
all_bias_de

#############################EMP SE#############################################
### separate the DAGS 
res_complex_ierdA <- res_complex_ierd[res_complex_ierd$causal_diagram == "A", ] 
res_complex_ierdB <- res_complex_ierd[res_complex_ierd$causal_diagram == "B", ] 
res_complex_ierdC <- res_complex_ierd[res_complex_ierd$causal_diagram == "C", ] 
res_complex_ierdD <- res_complex_ierd[res_complex_ierd$causal_diagram == "D", ] 
res_complex_ierdE <- res_complex_ierd[res_complex_ierd$causal_diagram == "E", ] 
res_complex_ierdF <- res_complex_ierd[res_complex_ierd$causal_diagram == "F", ] 
res_complex_ierdG <- res_complex_ierd[res_complex_ierd$causal_diagram == "G", ] 

res_complex_derdA <- res_complex_derd[res_complex_derd$causal_diagram == "A", ] 
res_complex_derdB <- res_complex_derd[res_complex_derd$causal_diagram == "B", ] 
res_complex_derdC <- res_complex_derd[res_complex_derd$causal_diagram == "C", ] 
res_complex_derdD <- res_complex_derd[res_complex_derd$causal_diagram == "D", ] 
res_complex_derdE <- res_complex_derd[res_complex_derd$causal_diagram == "E", ] 
res_complex_derdF <- res_complex_derd[res_complex_derd$causal_diagram == "F", ] 
res_complex_derdG <- res_complex_derd[res_complex_derd$causal_diagram == "G", ] 


int_ieempse_A <- ggplot(res_complex_ierdA, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieempse_B <- ggplot(res_complex_ierdB, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieempse_C <- ggplot(res_complex_ierdC, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieempse_D <- ggplot(res_complex_ierdD, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieempse_E <- ggplot(res_complex_ierdE, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieempse_F <- ggplot(res_complex_ierdF, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 


int_ieempse_G <- ggplot(res_complex_ierdG, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

all_empse_ie <- ggarrange(int_ieempse_A , int_ieempse_B,int_ieempse_C,int_ieempse_D,int_ieempse_E,int_ieempse_F,int_ieempse_G,   
                          ncol = 4, nrow = 2, common.legend = T,legend="bottom") 

all_empse_ie 


int_deempse_A <- ggplot(res_complex_derdA, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deempse_B <- ggplot(res_complex_derdB, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deempse_C <- ggplot(res_complex_derdC, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deempse_D <- ggplot(res_complex_derdD, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deempse_E <- ggplot(res_complex_derdE, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deempse_F <- ggplot(res_complex_derdF, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deempse_G <- ggplot(res_complex_derdG, aes(emp_se,analysis)) + 
  geom_errorbar(
    aes(xmin = emp_se-emp_se_se, xmax = emp_se+emp_se_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(0.04,4.0)) + 
  scale_y_discrete(name="") + 
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 


all_empse_de <- ggarrange(int_deempse_A , int_deempse_B,int_deempse_C,int_deempse_D,int_deempse_E,int_deempse_F,int_deempse_G,   
                          ncol = 4, nrow = 2, common.legend = T,legend="bottom")  

all_empse_de

#############################SE ERROR###########################################
int_ieseerror_A <- ggplot(res_complex_ierdA, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieseerror_B <- ggplot(res_complex_ierdB, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieseerror_C <- ggplot(res_complex_ierdC, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieseerror_D <- ggplot(res_complex_ierdD, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieseerror_E <- ggplot(res_complex_ierdE, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieseerror_F <- ggplot(res_complex_ierdF, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_ieseerror_G <- ggplot(res_complex_ierdG, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 


# Add margin adjustments to each plot individually
int_ieseerror_A <- int_ieseerror_A + theme(plot.margin = margin(0, 0, 0, 0))
int_ieseerror_B <- int_ieseerror_B + theme(plot.margin = margin(0, 0, 0, 0))
int_ieseerror_C <- int_ieseerror_C + theme(plot.margin = margin(0, 0, 0, 0))
int_ieseerror_D <- int_ieseerror_D + theme(plot.margin = margin(0, 0, 0, 0))
int_ieseerror_E <- int_ieseerror_E + theme(plot.margin = margin(0, 0, 0, 0))
int_ieseerror_F <- int_ieseerror_F + theme(plot.margin = margin(0, 0, 0, 0))
int_ieseerror_G <- int_ieseerror_G + theme(plot.margin = margin(0, 0, 0, 0))

# Use ggarrange to arrange the plots with reduced space
all_ieseerror_ie <- ggarrange(
  int_ieseerror_A, int_ieseerror_B, int_ieseerror_C,
  int_ieseerror_D, int_ieseerror_E, int_ieseerror_F,
  int_ieseerror_G, ncol = 2, nrow = 4,
  common.legend = TRUE, legend = "bottom",
  align = "hv",  # Align plots horizontally and vertically
  widths = c(1, 1, 1),  # Relative widths of columns
  heights = c(1, 1, 1),  # Relative heights of rows
  padding = 0  # Reduce space between plots
)

# Print the arranged plots
all_ieseerror_ie

int_deseerror_A <- ggplot(res_complex_derdA, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deseerror_B <- ggplot(res_complex_derdB, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deseerror_C <- ggplot(res_complex_derdC, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deseerror_D <- ggplot(res_complex_derdD, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deseerror_E <- ggplot(res_complex_derdE, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_deseerror_F <- ggplot(res_complex_derdF, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 


int_deseerror_G <- ggplot(res_complex_derdG, aes(relmodse_err,analysis)) + 
  geom_errorbar(
    aes(xmin = relmodse_err-relmodse_err_se, xmax = relmodse_err+relmodse_err_se, color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(-30,60)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

# Add margin adjustments to each plot individually
int_deseerror_A <- int_deseerror_A + theme(plot.margin = margin(0, 0, 0, 0))
int_deseerror_B <- int_deseerror_B + theme(plot.margin = margin(0, 0, 0, 0))
int_deseerror_C <- int_deseerror_C + theme(plot.margin = margin(0, 0, 0, 0))
int_deseerror_D <- int_deseerror_D + theme(plot.margin = margin(0, 0, 0, 0))
int_deseerror_E <- int_deseerror_E + theme(plot.margin = margin(0, 0, 0, 0))
int_deseerror_F <- int_deseerror_F + theme(plot.margin = margin(0, 0, 0, 0))
int_deseerror_G <- int_deseerror_G + theme(plot.margin = margin(0, 0, 0, 0))

# Use ggarrange to arrange the plots with reduced space
all_deseerror_de <- ggarrange(
  int_deseerror_A, int_deseerror_B, int_deseerror_C,
  int_deseerror_D, int_deseerror_E, int_deseerror_F,
  int_deseerror_G, ncol = 2, nrow = 4,
  common.legend = TRUE, legend = "bottom",
  align = "hv",  # Align plots horizontally and vertically
  widths = c(1, 1, 1),  # Relative widths of columns
  heights = c(1, 1, 1),  # Relative heights of rows
  padding = 0  # Reduce space between plots
)

# Print the arranged plots
all_deseerror_de

#############################COVERAGE PROBABILITY###############################
int_iecov_A <- ggplot(res_complex_ierdA, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_iecov_B <- ggplot(res_complex_ierdB, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_iecov_C <- ggplot(res_complex_ierdC, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_iecov_D <- ggplot(res_complex_ierdD, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_iecov_E <- ggplot(res_complex_ierdE, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_iecov_F <- ggplot(res_complex_ierdF, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_iecov_G <- ggplot(res_complex_ierdG, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 


# Add margin adjustments to each plot individually
int_iecov_A <- int_iecov_A + theme(plot.margin = margin(0, 0, 0, 0))
int_iecov_B <- int_iecov_B + theme(plot.margin = margin(0, 0, 0, 0))
int_iecov_C <- int_iecov_C + theme(plot.margin = margin(0, 0, 0, 0))
int_iecov_D <- int_iecov_D + theme(plot.margin = margin(0, 0, 0, 0))
int_iecov_E <- int_iecov_E + theme(plot.margin = margin(0, 0, 0, 0))
int_iecov_F <- int_iecov_F + theme(plot.margin = margin(0, 0, 0, 0))
int_iecov_G <- int_iecov_G + theme(plot.margin = margin(0, 0, 0, 0))

# Use ggarrange to arrange the plots with reduced space
all_iecov_ie <- ggarrange(
  int_iecov_A, int_iecov_B, int_iecov_C,
  int_iecov_D, int_iecov_E, int_iecov_F,
  int_iecov_G, ncol = 2, nrow = 4,
  common.legend = TRUE, legend = "bottom",
  align = "hv",  # Align plots horizontally and vertically
  widths = c(1, 1, 1),  # Relative widths of columns
  heights = c(1, 1, 1),  # Relative heights of rows
  padding = 0  # Reduce space between plots
)

# Print the arranged plots
all_iecov_ie


int_decov_A <- ggplot(res_complex_derdA, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG A", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_decov_B <- ggplot(res_complex_derdB, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG B", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_decov_C <- ggplot(res_complex_derdC, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG C", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_decov_D <- ggplot(res_complex_derdD, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG D", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_decov_E <- ggplot(res_complex_derdE, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG E", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 

int_decov_F <- ggplot(res_complex_derdF, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG F", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 


int_decov_G <- ggplot(res_complex_derdG, aes(cov,analysis)) + 
  geom_errorbar(
    aes(xmin = cov-(cov_se*100), xmax = cov+(cov_se*100), color = analysis),
    position = position_dodge(0.3), width = 0.5) + geom_point(size = 1, aes(color = analysis)) + 
  scale_x_continuous(name="m-DAG G", limits=c(0,100)) + 
  scale_y_discrete(name="") + geom_vline(xintercept = 95, linetype = "dashed", color = "red", size = 0.3) +   
  theme(text = element_text(size = 13), plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),
        legend.title=element_blank(),axis.title.x = element_text(size = 10)) 



# Add margin adjustments to each plot individually
int_decov_A <- int_decov_A + theme(plot.margin = margin(0, 0, 0, 0))
int_decov_B <- int_decov_B + theme(plot.margin = margin(0, 0, 0, 0))
int_decov_C <- int_decov_C + theme(plot.margin = margin(0, 0, 0, 0))
int_decov_D <- int_decov_D + theme(plot.margin = margin(0, 0, 0, 0))
int_decov_E <- int_decov_E + theme(plot.margin = margin(0, 0, 0, 0))
int_decov_F <- int_decov_F + theme(plot.margin = margin(0, 0, 0, 0))
int_decov_G <- int_decov_G + theme(plot.margin = margin(0, 0, 0, 0))

# Use ggarrange to arrange the plots with reduced space
all_decov_ie <- ggarrange(
  int_decov_A, int_decov_B, int_decov_C,
  int_decov_D, int_decov_E, int_decov_F,
  int_decov_G, ncol = 2, nrow = 4,
  common.legend = TRUE, legend = "bottom",
  align = "hv",  # Align plots horizontally and vertically
  widths = c(1, 1, 1),  # Relative widths of columns
  heights = c(1, 1, 1),  # Relative heights of rows
  padding = 0  # Reduce space between plots
)

# Print the arranged plots
all_decov_ie
