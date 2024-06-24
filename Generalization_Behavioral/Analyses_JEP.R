#DATA ANLYSIS JEP:General paper

#clear console, environment and plots
cat("\014")
rm(list = ls())

#Load libraries
library(tidyverse)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(effectsize)

#Set a seed to replicate same results every time
set.seed(9999)

#Folder paths
Main_path <- "/Users/pieter/Documents/Generalization_study/JEP/"
Plot_path <- paste(Main_path, "Plots/", sep="") 

#Load in data for all experiments and sessions
Experiments = 3
Sessions = 2

for (exp in 1:Experiments){
  Exp_path <- paste(Main_path, "Experiment", exp,"/", sep="") 
  for (ses in 1:Sessions){
    Sess_path <- paste(Exp_path, "Session", ses,"/", sep="") 
    file_paths <- list.files(path = Sess_path, full.names = TRUE, pattern = "\\.csv$")
    
    session_data <- file_paths %>% map_df(read.csv)
    
    subjects <- c()
    for (s in file_paths){
      if (exp ==1){
        subjects <- append(subjects, substr(s, nchar(Sess_path)+19, nchar(Sess_path)+24))
      }else{
        subjects <- append(subjects, substr(s, nchar(Sess_path)+19, nchar(Sess_path)+23))
      }
    }
    
    trials <- nrow(read.csv(file_paths[1]))
    session_data$Subject <- rep(subjects, each = trials)
    
    Dat_name <- paste("Data_Exp_", exp, "_Sess_", ses, sep = "")
    assign(Dat_name, session_data)
    
    Sub_name <- paste("Subs_Exp_", exp, "_Sess_", ses, sep = "")
    assign(Sub_name, subjects)
    
  }
}
#############################################################################
#################################  Experiment 1 #############################
#############################################################################

#First experiment was special since it had 3 sessions (third session in MRI)
#we deal with this here
Exp1_sess3_path <- paste(Main_path, "Experiment1/Session3/", sep="") 

for (run in 0:2){
  file_paths <- list.files(path = Exp1_sess3_path, full.names = TRUE, pattern = paste("\\Run_", run, ".csv", sep = ""))
  Run_data <- file_paths %>% map_df(read.csv)
  if (run ==0){
    Run_data$round <- Run_data$round +9
  }else{
    Run_data$round <- Run_data$round +11
  }
  
  subjects <- c()
  for (s in file_paths){
    subjects <- append(subjects, substr(s, nchar(Exp1_sess3_path)+19, nchar(Exp1_sess3_path)+24))
  }
  
  trials <- nrow(read.csv(file_paths[1]))
  Run_data$Subject <- rep(subjects, each = trials)
  
  Dat_name <- paste("Data_Exp_1_Run_", run, sep = "")
  assign(Dat_name, Run_data)
}

#different columns in run0 (out of scanner) compared to other runs
sess3_dat <- rbind(Data_Exp_1_Run_0[,c(1:19,26:28,30)], Data_Exp_1_Run_1[,c(1:15,24:27,29:31,33)], Data_Exp_1_Run_2[,c(1:15,24:27,29:31,33)])
Data_Exp_1_Sess_3 <-data.frame(matrix(nrow = nrow(sess3_dat), ncol = ncol(Data_Exp_1_Sess_1)))
colnames(Data_Exp_1_Sess_3) <-colnames(Data_Exp_1_Sess_1)

#also different columns for first two days compared to third day
for (n in colnames(Data_Exp_1_Sess_3)){
  c1 <- which(n == colnames(Data_Exp_1_Sess_3)) 
  if (n %in% colnames(sess3_dat)){
    c2 <- which(n == colnames(sess3_dat)) 
    Data_Exp_1_Sess_3[, c1] <- sess3_dat[,c2]
  }else{
    Data_Exp_1_Sess_3[, c1] <- NaN
  }
}

#Check subjects in third session
Subs_Exp_1_Sess_3 <- unique(Data_Exp_1_Sess_3$Subject)

#Check all subjects and which ones did all sessions
Exp_1_allSubs <- unique(c(Subs_Exp_1_Sess_1,Subs_Exp_1_Sess_2,Subs_Exp_1_Sess_3))
Exp_1_selected_subs <- intersect(Subs_Exp_1_Sess_1,Subs_Exp_1_Sess_2)
Exp_1_selected_subs <- intersect(Exp_1_selected_subs,Subs_Exp_1_Sess_3)

#No subjects removed
removed_subs_Exp1 <- Exp_1_allSubs[! Exp_1_allSubs %in% Exp_1_selected_subs]

#Combine all data for subjects that did all sessions
Exp1_allData <- rbind(Data_Exp_1_Sess_1[Data_Exp_1_Sess_1$Subject %in% Exp_1_selected_subs,], Data_Exp_1_Sess_2[Data_Exp_1_Sess_2$Subject %in% Exp_1_selected_subs,], Data_Exp_1_Sess_3[Data_Exp_1_Sess_3$Subject %in% Exp_1_selected_subs,])

#Compute errrorscore
Exp1_allData$Errorscore <- abs(Exp1_allData$coordID -Exp1_allData$respID)*100

#Compute errorscore based on random behavior 
for (iterations in seq(1:1000)){
  randombaseline <- runif(n = nrow(Exp1_allData), min = 0, max = 1)
  random_error1 <- abs(Exp1_allData$coordID - randombaseline)*100
  if (iterations ==1){
    random_error <- random_error1
  }else{
    random_error<-random_error+random_error1
  }
}

#Mean over 1000 repetitions
Exp1_allData$randombaseline<-random_error/1000

#Baseline participant errorscore base on random baseline
Exp1_allData$BaselinedScore <- Exp1_allData$Errorscore / Exp1_allData$randombaseline

#Subjects that were worse than random in the last test block of the old animals were also removed
round8<- aggregate(Exp1_allData$BaselinedScore[Exp1_allData$round ==8], list(Exp1_allData$Subject[Exp1_allData$round ==8]), FUN = mean)
good_enough <- round8$Group.1[round8$x < 1]
#No rejections
Exp1_allData <- Exp1_allData[Exp1_allData$Subject %in% good_enough, ]

#Check performance across blocks/rounds
Blocked_initialaggregate <- aggregate(list(Exp1_allData$BaselinedScore), list(Exp1_allData$Subject, Exp1_allData$round, Exp1_allData$blocktype), FUN = mean)
colnames(Blocked_initialaggregate)<- c("Subject", "Round", "Blocktype", "Errorscore")

Blocked_errorscores <- summarySE(data = Blocked_initialaggregate, measurevar = c("Errorscore"), groupvars = c("Round", "Blocktype"))

#make plot
Block_plot <- ggplot(Blocked_errorscores, aes(x = Round, y = Errorscore, color = factor(Blocktype), shape = factor(Blocktype))) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", size = 1.5)+
  geom_errorbar(aes(ymin = Errorscore - ci, ymax = Errorscore + ci), width = 0.2) +
  scale_shape_manual(values = c("0" = 4, "1" = 16, "2" = 4, "3"= 16), labels = c("Initial learning", "Initial test", "Generalization learning", "Generalization test")) +
  scale_color_manual(values = c("0" = "blue", "1" = "blue", "2"= "red", "3" = "red"), labels = c("Initial learning", "Initial test", "Generalization learning", "Generalization test")) +
  geom_vline(xintercept = 4.33, linetype="dotted", color = "grey", size=1.5)+
  geom_vline(xintercept = 8.33, linetype="dotted", color = "grey", size=1.5)+
  #geom_point(size = 4, shape = 18, color = "black", aes(y = Blocked_baseline$Baseline, x = Blocked_baseline$Round ))+
  scale_x_continuous(breaks = 1:max(Blocked_errorscores$Round)) +  
  ylim(c(0,1))+
  theme_classic() +
  labs(x = "Round number",
       y = "Baselined errorscore") +
  guides(color = guide_legend("Round"), shape = guide_legend("Round"))

#Check plot and save it
Block_plot
ggsave("Block_plot_Exp1.jpg", Block_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Get all data for Anova
Block_data <- aggregate(list(Exp1_allData$BaselinedScore), list(Exp1_allData$Subject,Exp1_allData$round, Exp1_allData$blocktype), FUN = mean)
colnames(Block_data)<- c("Subject", "Blocknumber", "Blocktype", "Errorscore")
Block_data$Test <- Block_data$Blocktype %%2
Block_data$Animals <- floor(Block_data$Blocktype /2)

#Perform Anova
Model_errorscore <- aov(Errorscore ~ Blocknumber*factor(Test)*factor(Animals)+Error(factor(Subject)/(Blocknumber*factor(Test)*factor(Animals))), data = Block_data)
summary(Model_errorscore)

# Error: factor(Subject)
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 48  19.76  0.4117               
# 
# Error: factor(Subject):Blocknumber
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Blocknumber  1  2.851  2.8514   55.73 1.43e-09 ***
#   Residuals   48  2.456  0.0512                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Test)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(Test)  1 11.513   11.51   114.6 2.58e-14 ***
#   Residuals    48  4.821    0.10                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Animals)
# Df Sum Sq Mean Sq F value Pr(>F)  
# factor(Animals)  1 0.2444 0.24443   7.129 0.0103 *
#   Residuals       48 1.6458 0.03429                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(Test)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Blocknumber:factor(Test)  1 0.2972 0.29717    12.3 0.000994 ***
#   Residuals                48 1.1598 0.02416                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(Animals)
# Df Sum Sq  Mean Sq F value Pr(>F)
# Blocknumber:factor(Animals)  1 0.0005 0.000549   0.036   0.85
# Residuals                   48 0.7323 0.015255               
# 
# Error: factor(Subject):factor(Test):factor(Animals)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(Test):factor(Animals)  1  1.486   1.486   47.98 9.52e-09 ***
#   Residuals                    48  1.486   0.031                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(Test):factor(Animals)
# Df Sum Sq Mean Sq F value  Pr(>F)    
# Blocknumber:factor(Test):factor(Animals)  1  1.771  1.7712    63.6 2.4e-10 ***
#   Residuals                                48  1.337  0.0278                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: Within
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 686  9.274 0.01352                

eta_squared(Model_errorscore)
# # Effect Size for ANOVA (Type I)
# 
# Group                                                    |                                Parameter | Eta2 (partial) |       95% CI
# -----------------------------------------------------------------------------------------------------------------------------------
#   factor(Subject):Blocknumber                              |                              Blocknumber |           0.54 | [0.38, 1.00]
# factor(Subject):factor(Test)                             |                             factor(Test) |           0.70 | [0.59, 1.00]
# factor(Subject):factor(Animals)                          |                          factor(Animals) |           0.13 | [0.02, 1.00]
# factor(Subject):Blocknumber:factor(Test)                 |                 Blocknumber:factor(Test) |           0.20 | [0.06, 1.00]
# factor(Subject):Blocknumber:factor(Animals)              |              Blocknumber:factor(Animals) |       7.49e-04 | [0.00, 1.00]
# factor(Subject):factor(Test):factor(Animals)             |             factor(Test):factor(Animals) |           0.50 | [0.33, 1.00]
# factor(Subject):Blocknumber:factor(Test):factor(Animals) | Blocknumber:factor(Test):factor(Animals) |           0.57 | [0.41, 1.00]

# To get more insight in blocktype differences
Test_animals_interactie <- aggregate(list(Exp1_allData$BaselinedScore), list(Exp1_allData$blocktype, Exp1_allData$Subject), FUN = mean)
colnames(Test_animals_interactie) <- c("Blocktype", "Subject","Errorscore")
pairwise.t.test(Test_animals_interactie$Errorscore, Test_animals_interactie$Blocktype, paired = TRUE)

# Pairwise comparisons using paired t tests 
# 
# data:  Test_animals_interactie$Errorscore and Test_animals_interactie$Blocktype 
# 
# 0       1       2      
# 1 8.5e-11 -       -      
#   2 0.0126  8.5e-11 -      
#   3 5.3e-07 0.0071  2.1e-09
# 
# P value adjustment method: holm 

# For experiment 1 there was context information that we need to read in and add to data of generalization test blocks
cont_info <- read.csv(paste(Main_path, "Experiment1/ContextInfo.csv", sep=""))

#Now we add for each context, which type it was
Generalization_test_data <- Exp1_allData[Exp1_allData$blocktype==3,]
for (x in 1:nrow(cont_info)){
  Generalization_test_data$ContextType[((Generalization_test_data$Subject == cont_info$Subject[x]) + (Generalization_test_data$locationID==cont_info$Training_1[x]))==2]<-"Train"
  Generalization_test_data$ContextType[((Generalization_test_data$Subject == cont_info$Subject[x]) + (Generalization_test_data$locationID==cont_info$Training_2[x]))==2]<-"Train"
  Generalization_test_data$ContextType[((Generalization_test_data$Subject == cont_info$Subject[x]) + (Generalization_test_data$locationID==cont_info$Full_1[x]))==2]<-"Full"
  Generalization_test_data$ContextType[((Generalization_test_data$Subject == cont_info$Subject[x]) + (Generalization_test_data$locationID==cont_info$Full_2[x]))==2]<-"Full"
  Generalization_test_data$ContextType[((Generalization_test_data$Subject == cont_info$Subject[x]) + (Generalization_test_data$locationID==cont_info$Anti_1[x]))==2]<-"Anti"
  Generalization_test_data$ContextType[((Generalization_test_data$Subject == cont_info$Subject[x]) + (Generalization_test_data$locationID==cont_info$Anti_2[x]))==2]<-"Anti"
}

#aggregate data of generalization test blocks
Gen_data <- aggregate(list(Generalization_test_data$BaselinedScore), list(Generalization_test_data$Subject,Generalization_test_data$round, Generalization_test_data$ContextType), FUN = mean)
colnames(Gen_data)<- c("Subject", "Blocknumber", "ContextType", "Errorscore")

#check whether there were subjects that were below chance level on the trained contexts in last generalization test block
check <- Gen_data[Gen_data$ContextType == "Train",]
check <- check[check$Blocknumber ==22,]
check$Subject[check$Baselined_error< 0]

#Perform anova for contexttypes
Model_generalization <- aov(Errorscore ~ Blocknumber*factor(ContextType)+Error(factor(Subject)/(Blocknumber*factor(ContextType))), data = Gen_data)
summary(Model_generalization)

# Error: factor(Subject)
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 48  49.47   1.031               
# 
# Error: factor(Subject):Blocknumber
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Blocknumber  1  3.038  3.0383   26.28 5.24e-06 ***
#   Residuals   48  5.550  0.1156                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(ContextType)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(ContextType)  2  31.86  15.929   39.74 2.66e-13 ***
#   Residuals           96  38.48   0.401                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(ContextType)
# Df Sum Sq Mean Sq F value Pr(>F)
# Blocknumber:factor(ContextType)  2  0.003 0.00130   0.022  0.978
# Residuals                       96  5.613 0.05847               
# 
# Error: Within
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 588  22.26 0.03786   

eta_squared(Model_generalization)
# # Effect Size for ANOVA (Type I)
# 
# Group                                           |                       Parameter | Eta2 (partial) |       95% CI
# -----------------------------------------------------------------------------------------------------------------
#   factor(Subject):Blocknumber                     |                     Blocknumber |           0.35 | [0.18, 1.00]
# factor(Subject):factor(ContextType)             |             factor(ContextType) |           0.45 | [0.33, 1.00]
# factor(Subject):Blocknumber:factor(ContextType) | Blocknumber:factor(ContextType) |       4.62e-04 | [0.00, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].

#Pairwise t-tests across contexttypes
Gen_type_dat <- aggregate(Gen_data$Errorscore, list(Gen_data$Subject, Gen_data$ContextType), FUN = mean)
colnames(Gen_type_dat)<-c("Subject", "ContextType", "Errorscore")
pairwise.t.test(Gen_type_dat$Errorscore, Gen_type_dat$ContextType, paired= TRUE)
# Pairwise comparisons using paired t tests 
# 
# data:  Gen_type_dat$Errorscore and Gen_type_dat$ContextType 
# 
# Anti    Full   
# Full  0.00023 -      
#   Train 8.6e-11 1.6e-06
# 
# P value adjustment method: holm 

#Here, again comparison by comparison to check t-values
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)

#Also check cohen's d
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)

#Check mean and sd for each contexttype
mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"])
#0.1872871
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"])
#0.1624696

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"])
# 0.4115802
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"])
# 0.3430434

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"])
#0.652722
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"])
#0.4016165

#check for each context whether they were better than baseline (below 1)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Train"] -     1
# t = -35.016, df = 48, p-value < 2.2e-16
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.7737847
# sample estimates:
#   mean of x 
# -0.8127129 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Full"] -     1
# t = -12.007, df = 48, p-value = 2.285e-16
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.5062255
# sample estimates:
#   mean of x 
# -0.5884198 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Anti"] -     1
# t = -6.0529, df = 48, p-value = 1.041e-07
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.2510493
# sample estimates:
#   mean of x 
# -0.347278 

#Make plot of generalization performance
Generalisation_toplot<- summarySE(data = Gen_data, measurevar = "Errorscore", groupvars = c("ContextType", "Blocknumber"))

Generalisation_plot <- ggplot(Generalisation_toplot, aes(x = Blocknumber, y = Errorscore, color = factor(ContextType))) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", size = 1.5)+
  geom_point(size = 4, shape = 16, position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = Errorscore- ci, ymax = Errorscore + ci), width = 0.2, position = position_dodge(width = .5)) +
  scale_color_manual(labels = c("Anti", "Full", "Trained"), values = c("Anti" = "magenta", "Full" = "cyan", "Train"= "gold")) +
  scale_x_continuous("Round number", breaks = c(11,13, 16,17,19,22), labels = c(11,13, 16,17,19,22))+
  ylim(c(0,1))+
  theme_classic() +
  labs(x = "Round",
       y = "Baselined errorscore") +
  guides(color = guide_legend("Type"), shape = guide_legend("Type"))

Generalisation_plot
ggsave("Generalisation_plot_Exp1.jpg", Generalisation_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Check in each block what they learned
CDat_try <- aggregate(x = Generalization_test_data$BaselinedScore, by = list(Generalization_test_data$Subject, Generalization_test_data$round, Generalization_test_data$ContextType), FUN = mean)
colnames(CDat_try) <- c("Subject", "Round", "Type", "BaselinedScore")
CDat_try$Rand <- CDat_try$BaselinedScore < 1
CDat_try$Good <- CDat_try$BaselinedScore < .5

CDat_try$RandNone <- FALSE
CDat_try$GoodNone <- FALSE
CDat_try$RandTrained <- FALSE
CDat_try$GoodTrained <- FALSE
CDat_try$RandFull <- FALSE
CDat_try$GoodFull <- FALSE
CDat_try$RandAnti <- FALSE
CDat_try$GoodAnti <- FALSE

CDat_try$RandCombined <- FALSE
CDat_try$GoodCombined <- FALSE
CDat_try$RandTFull <- FALSE
CDat_try$RandTAnti <- FALSE
CDat_try$GoodTFull <- FALSE
CDat_try$GoodTAnti <- FALSE

CDat_try$RandAll <- FALSE
CDat_try$GoodAll <- FALSE

for (s in Exp_1_selected_subs){
  d <- CDat_try[CDat_try$Subject == s,]
  for (r in unique(CDat_try$Round)){
    d2 <- d[d$Round == r,]
    
    if (sum(d2$Rand)==3){
      CDat_try$RandAll[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Rand)==0){
      CDat_try$RandNone[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Rand)==1){
      if (d2$Rand[d2$Type=="Anti"]){
        CDat_try$RandAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]){
        CDat_try$RandFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]){
        CDat_try$RandTrained[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }else if (sum(d2$Rand)==2){
      if (d2$Rand[d2$Type=="Anti"]==FALSE){
        CDat_try$RandTFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]==FALSE){
        CDat_try$RandTAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]==FALSE){
        CDat_try$RandCombined[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }
    
    if (sum(d2$Good)==3){
      CDat_try$GoodAll[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Good)==0){
      CDat_try$GoodNone[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Good)==1){
      if (d2$Good[d2$Type=="Anti"]){
        CDat_try$GoodAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Good[d2$Type=="Full"]){
        CDat_try$GoodFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Good[d2$Type=="Train"]){
        CDat_try$GoodTrained[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }else if (sum(d2$Good)==2){
      if (d2$Good[d2$Type=="Anti"]==FALSE){
        CDat_try$GoodTFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Good[d2$Type=="Full"]==FALSE){
        CDat_try$GoodTAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Good[d2$Type=="Train"]==FALSE){
        CDat_try$GoodCombined[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }
  }
}

#aggregate
aggregate_Rand <- aggregate(list(CDat_try$RandNone, CDat_try$RandTrained, CDat_try$RandFull, CDat_try$RandAnti, CDat_try$RandCombined, CDat_try$RandTFull, CDat_try$RandTAnti, CDat_try$RandAll), list(CDat_try$Round), FUN = sum )
colnames(aggregate_Rand) <- c("Block", "None", "OnlyTrain", "OnlyFull", "OnlyAnti", "AntiFull", "TrainedFull", "TrainedAnti", "All")
aggregate_Rand[,2:9] <-aggregate_Rand[,2:9]/3

#see who learned full and full without anti
aggregate_Rand$full <- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull + aggregate_Rand$AntiFull + aggregate_Rand$All
aggregate_Rand$fullnoanti <- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull

#compute proportions
sum(aggregate_Rand$fullnoanti)/(49*6)
sum(aggregate_Rand$fullnoanti)/sum(aggregate_Rand$full)

#see who learned anti and anti without full
aggregate_Rand$anti <- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti + + aggregate_Rand$AntiFull + aggregate_Rand$All
aggregate_Rand$antinofull <- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti

#compute proportions
sum(aggregate_Rand$antinofull)/(49*6)
sum(aggregate_Rand$antinofull)/sum(aggregate_Rand$anti)

#Do proportional test to see whether it is more likely to have learned full without anti than inverse
prop.test(x = c(sum(aggregate_Rand$fullnoanti), sum(aggregate_Rand$antinofull)), n = c(sum(aggregate_Rand$full), sum(aggregate_Rand$anti)))
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(sum(aggregate_Rand$fullnoanti), sum(aggregate_Rand$antinofull)) out of c(sum(aggregate_Rand$full), sum(aggregate_Rand$anti))
# X-squared = 23.581, df = 1, p-value = 1.198e-06
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   0.1227796 0.2796633
# sample estimates:
#   prop 1    prop 2 
# 0.3281250 0.1269036 

prop.test(x = c(0, sum(aggregate_Rand$antinofull)), n = c(49*6, 49*6))

#prepare data for plotting
Rand_long <- gather(aggregate_Rand, Type, Proportion, None:All, factor_key = TRUE)
plot_learnseq_random <-ggplot(data = Rand_long, aes(x = Block, y = (Proportion/49)*100, fill = Type, group = Type))+
  geom_bar(stat="identity", position = "dodge")+
  scale_fill_manual("Learned Contexts", values=c("None" = "black", "OnlyTrain"="gold", "OnlyFull"= "cyan", "OnlyAnti"="magenta", "AntiFull" = "purple", "TrainedFull"= "orange", "TrainedAnti"= "green", "All"= "red"),labels=c("None" = "None", "OnlyTrain"="Only Trained", "OnlyFull"= "Only Full", "OnlyAnti"="Only Anti", "AntiFull" = "Anti & Full", "TrainedFull"= "Trained & Full", "TrainedAnti"= "Trained & Anti", "All"= "All"))+
  scale_y_continuous("Proportion subjects below baseline")+
  scale_x_continuous("Round number", breaks = c(0:5), labels = c(0:5))+
  theme_classic()+
  theme(text = element_text(size=10))

ggsave("Learning_Baseline_Exp1.jpg", plot_learnseq_random, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#more specific plot (used in paper)
Rand_bis <- gather(aggregate_Rand, Type, Proportion, full:antinofull, factor_key = TRUE)
plot_learnseq_random_bis <-ggplot(data = Rand_bis, aes(x = Block, y = (Proportion/49)*100, fill = Type, group = Type))+
  geom_bar(stat="identity", position = "dodge")+
  scale_fill_manual("Learned", values=c("full" = "cyan", "fullnoanti"="darkblue", "anti"= "magenta", "antinofull"="purple"),labels=c("full" = "Full", "fullnoanti"="Full but not Anti", "anti"= "Anti", "antinofull"="Anti but not Full"))+
  scale_y_continuous("Proportion participants")+
  scale_x_continuous("Round number", breaks = c(11,13, 16,17,19,22), labels = c(11,13, 16,17,19,22))+
  theme_classic()+
  theme(text = element_text(size=10))

ggsave("Learning_Baseline_Exp1_bis.jpg", plot_learnseq_random_bis, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#we double check robustness of this finding by doing it again with baseline/2 as threshold (definitely not random anymore)
aggregate_Good <- aggregate(list(CDat_try$GoodNone, CDat_try$GoodTrained, CDat_try$GoodFull, CDat_try$GoodAnti, CDat_try$GoodCombined, CDat_try$GoodTFull, CDat_try$GoodTAnti, CDat_try$GoodAll), list(CDat_try$Round), FUN = sum )
colnames(aggregate_Good) <- c("Block", "None", "OnlyTrain", "OnlyFull", "OnlyAnti", "AntiFull", "TrainedFull", "TrainedAnti", "All")
aggregate_Good[,2:9] <-aggregate_Good[,2:9]/3

aggregate_Good$full <- aggregate_Good$OnlyFull + aggregate_Good$TrainedFull + aggregate_Good$AntiFull + aggregate_Good$All
aggregate_Good$fullnoanti <- aggregate_Good$OnlyFull + aggregate_Good$TrainedFull

sum(aggregate_Good$fullnoanti)/(49*6)
sum(aggregate_Good$fullnoanti)/sum(aggregate_Good$full)

aggregate_Good$anti <- aggregate_Good$OnlyAnti + aggregate_Good$TrainedAnti + + aggregate_Good$AntiFull + aggregate_Good$All
aggregate_Good$antinofull <- aggregate_Good$OnlyAnti + aggregate_Good$TrainedAnti

sum(aggregate_Good$antinofull)/(49*6)
sum(aggregate_Good$antinofull)/sum(aggregate_Good$anti)

#Do proportional test to see whether it is more likely to have learned full without anti than inverse
prop.test(x = c(sum(aggregate_Good$fullnoanti), sum(aggregate_Good$antinofull)), n = c(sum(aggregate_Good$full), sum(aggregate_Good$anti)))
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(sum(aggregate_Good$fullnoanti), sum(aggregate_Good$antinofull)) out of c(sum(aggregate_Good$full), sum(aggregate_Good$anti))
# X-squared = 23.363, df = 1, p-value = 1.342e-06
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   0.1790174 0.4021418
# sample estimates:
#   prop 1    prop 2 
# 0.4797688 0.1891892 

Good_long <- gather(aggregate_Good, Type, Proportion, None:All, factor_key = TRUE)
plot_learnseq_good <-ggplot(data = Good_long, aes(x = Block, y = (Proportion/49)*100, fill = Type, group = Type))+
  geom_bar(stat="identity", position = "dodge")+
  scale_fill_manual("Learned Contexts", values=c("None" = "black", "OnlyTrain"="gold", "OnlyFull"= "cyan", "OnlyAnti"="magenta", "AntiFull" = "purple", "TrainedFull"= "orange", "TrainedAnti"= "green", "All"= "red"),labels=c("None" = "None", "OnlyTrain"="Only Trained", "OnlyFull"= "Only Full", "OnlyAnti"="Only Anti", "AntiFull" = "Anti & Full", "TrainedFull"= "Trained & Full", "TrainedAnti"= "Trained & Anti", "All"= "All"))+
  scale_y_continuous("Proportion subjects below baseline/2")+
  scale_x_continuous("Round number", breaks = c(0:5), labels = c(0:5))+
  theme_minimal()+
  theme(text = element_text(size=10))

ggsave("Learning_Good_Exp1.jpg", plot_learnseq_good, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Combine all relevant plots in one
Exp1_plot<-ggarrange(Block_plot, Generalisation_plot, plot_learnseq_random_bis, nrow = 3)
ggsave("All_Exp1.jpg", Exp1_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 15, units = "cm")

#############################################################################
#################################  Experiment 2 #############################
#############################################################################
Exp_2_allSubs <- unique(c(Subs_Exp_2_Sess_1,Subs_Exp_2_Sess_2))
Exp_2_selected_subs <- intersect(Subs_Exp_2_Sess_1,Subs_Exp_2_Sess_2)

removed_subs_Exp2 <- Exp_2_allSubs[! Exp_2_allSubs %in% Exp_2_selected_subs]
#4 participants did not complete all sessions

Exp2_allData <- rbind(Data_Exp_2_Sess_1[Data_Exp_2_Sess_1$Subject %in% Exp_2_selected_subs,], Data_Exp_2_Sess_2[Data_Exp_2_Sess_2$Subject %in% Exp_2_selected_subs,])
Exp2_allData$Errorscore <- abs(Exp2_allData$coordID -Exp2_allData$respID)*100

for (iterations in seq(1:1000)){
  randombaseline <- runif(n = nrow(Exp2_allData), min = 0, max = 1)
  random_error1 <- abs(Exp2_allData$coordID - randombaseline)*100
  if (iterations ==1){
    random_error <- random_error1
  }else{
    random_error<-random_error+random_error1
  }
}

Exp2_allData$randombaseline<-random_error/1000

Exp2_allData$BaselinedScore <- Exp2_allData$Errorscore / Exp2_allData$randombaseline

#Subjects that were worse than random in the last test block of the old animals were also removed
round4<- aggregate(Exp2_allData$BaselinedScore[Exp2_allData$round ==4], list(Exp2_allData$Subject[Exp2_allData$round ==4]), FUN = mean)
good_enough <- round4$Group.1[round4$x < 1]
#4 subjects were not better than random in final test block with old animals, we exclude them below

Exp2_allData <- Exp2_allData[Exp2_allData$Subject %in% good_enough, ]
Exp_2_selected_subs <- Exp_2_selected_subs[Exp_2_selected_subs%in%good_enough]

#Check performance across blocks/rounds
Blocked_initialaggregate <- aggregate(list(Exp2_allData$BaselinedScore), list(Exp2_allData$Subject, Exp2_allData$round, Exp2_allData$blocktype), FUN = mean)
colnames(Blocked_initialaggregate)<- c("Subject", "Round", "Blocktype", "Errorscore")

Blocked_errorscores <- summarySE(data = Blocked_initialaggregate, measurevar = c("Errorscore"), groupvars = c("Round", "Blocktype"))

#make plot
Block_plot <- ggplot(Blocked_errorscores, aes(x = Round, y = Errorscore, color = factor(Blocktype), shape = factor(Blocktype))) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype="dashed", color = "black", size=1.5)+
  geom_errorbar(aes(ymin = Errorscore - ci, ymax = Errorscore + ci), width = 0.2) +
  scale_shape_manual(values = c("0" = 4, "1" = 16, "2" = 4, "3"= 16), labels = c("Initial learning", "Initial test", "Generalization learning", "Generalization test")) +
  scale_color_manual(values = c("0" = "blue", "1" = "blue", "2"= "red", "3" = "red"), labels = c("Initial learning", "Initial test", "Generalization learning", "Generalization test")) +
  geom_vline(xintercept = 4.33, linetype="dotted", color = "grey", size=1.5)+
  #geom_point(size = 4, shape = 18, color = "black", aes(y = Blocked_baseline$Baseline, x = Blocked_baseline$Round ))+
  scale_x_continuous(breaks = 1:max(Blocked_errorscores$Round)) +  
  ylim(c(0,1))+
  theme_classic() +
  labs(x = "Round number",
       y = "Baselined errorscore") +
  guides(color = guide_legend("Round"), shape = guide_legend("Round"))

#show and save plot
Block_plot
ggsave("Block_plot_Exp2.jpg", Block_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Get all data for Anova
Block_data <- aggregate(list(Exp2_allData$BaselinedScore), list(Exp2_allData$Subject,Exp2_allData$round, Exp2_allData$blocktype), FUN = mean)
colnames(Block_data)<- c("Subject", "Blocknumber", "Blocktype", "Errorscore")
Block_data$Test <- Block_data$Blocktype %%2
Block_data$Animals <- floor(Block_data$Blocktype /2)

#Perform Anova
Model_errorscore <- aov(Errorscore ~ Blocknumber*factor(Test)*factor(Animals)+Error(factor(Subject)/(Blocknumber*factor(Test)*factor(Animals))), data = Block_data)
summary(Model_errorscore)

# Error: factor(Subject)
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 36  9.468   0.263               
# 
# Error: factor(Subject):Blocknumber
# Df Sum Sq Mean Sq F value Pr(>F)  
# Blocknumber  1 0.2022 0.20222   6.238 0.0172 *
#   Residuals   36 1.1670 0.03242                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Test)
# Df Sum Sq Mean Sq F value Pr(>F)    
# factor(Test)  1 19.736  19.736     223 <2e-16 ***
#   Residuals    36  3.185   0.088                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Animals)
# Df Sum Sq  Mean Sq F value Pr(>F)
# factor(Animals)  1 0.0028 0.002796   0.105  0.748
# Residuals       36 0.9613 0.026703               
# 
# Error: factor(Subject):Blocknumber:factor(Test)
# Df Sum Sq Mean Sq F value Pr(>F)
# Blocknumber:factor(Test)  1 0.0314 0.03136   1.189  0.283
# Residuals                36 0.9497 0.02638               
# 
# Error: factor(Subject):Blocknumber:factor(Animals)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# Blocknumber:factor(Animals)  1 0.1787 0.17872   11.72 0.00156 **
#   Residuals                   36 0.5490 0.01525                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Test):factor(Animals)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(Test):factor(Animals)  1 0.6350  0.6350   30.81 2.77e-06 ***
#   Residuals                    36 0.7419  0.0206                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(Test):factor(Animals)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Blocknumber:factor(Test):factor(Animals)  1 0.5555  0.5555   23.68 2.26e-05 ***
#   Residuals                                36 0.8446  0.0235                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: Within
# Df Sum Sq  Mean Sq F value Pr(>F)
# Residuals 185  1.467 0.007931  

eta_squared(Model_errorscore)
# # Effect Size for ANOVA (Type I)
# 
# Group                                                    |                                Parameter | Eta2 (partial) |       95% CI
# -----------------------------------------------------------------------------------------------------------------------------------
#   factor(Subject):Blocknumber                              |                              Blocknumber |           0.15 | [0.02, 1.00]
# factor(Subject):factor(Test)                             |                             factor(Test) |           0.86 | [0.79, 1.00]
# factor(Subject):factor(Animals)                          |                          factor(Animals) |       2.90e-03 | [0.00, 1.00]
# factor(Subject):Blocknumber:factor(Test)                 |                 Blocknumber:factor(Test) |           0.03 | [0.00, 1.00]
# factor(Subject):Blocknumber:factor(Animals)              |              Blocknumber:factor(Animals) |           0.25 | [0.07, 1.00]
# factor(Subject):factor(Test):factor(Animals)             |             factor(Test):factor(Animals) |           0.46 | [0.26, 1.00]
# factor(Subject):Blocknumber:factor(Test):factor(Animals) | Blocknumber:factor(Test):factor(Animals) |           0.40 | [0.19, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].

# To get more insight in blocktype differences
Test_animals_interactie <- aggregate(list(Exp2_allData$BaselinedScore), list(Exp2_allData$blocktype, Exp2_allData$Subject), FUN = mean)
colnames(Test_animals_interactie) <- c("Blocktype", "Subject","Errorscore")
pairwise.t.test(Test_animals_interactie$Errorscore, Test_animals_interactie$Blocktype, paired = TRUE)

# Pairwise comparisons using paired t tests 
# 
# data:  Test_animals_interactie$Errorscore and Test_animals_interactie$Blocktype 
# 
# 0       1       2      
# 1 5.3e-16 -       -      
#   2 2.5e-06 < 2e-16 -      
#   3 5.1e-11 0.45    4.0e-14
# 
# P value adjustment method: holm 

#Assign contexttype information
Generalization_test_data <- Exp2_allData[Exp2_allData$blocktype==3,]
Generalization_test_data$ContextType[Generalization_test_data$locationID==0]<-"Train"
Generalization_test_data$ContextType[Generalization_test_data$locationID==1]<-"Full"
Generalization_test_data$ContextType[Generalization_test_data$locationID==2]<-"Full"
Generalization_test_data$ContextType[Generalization_test_data$locationID==3]<-"Train"
Generalization_test_data$ContextType[Generalization_test_data$locationID==4]<-"Compositional"
Generalization_test_data$ContextType[Generalization_test_data$locationID==5]<-"Anti"
Generalization_test_data$ContextType[Generalization_test_data$locationID==6]<-"Compositional"
Generalization_test_data$ContextType[Generalization_test_data$locationID==7]<-"Anti"

#aggregate data of generalization test blocks
Gen_data <- aggregate(list(Generalization_test_data$BaselinedScore), list(Generalization_test_data$Subject,Generalization_test_data$round, Generalization_test_data$ContextType), FUN = mean)
colnames(Gen_data)<- c("Subject", "Blocknumber", "ContextType", "Errorscore")

#check whether there were subjects that were below chance level on the trained contexts in last generalization test block
check <- Gen_data[Gen_data$ContextType == "Train",]
check <- check[check$Blocknumber ==13,]
check$Subject[check$Baselined_error< 0]

#Perform anova for contexttypes
Model_generalization <- aov(Errorscore ~ Blocknumber*factor(ContextType)+Error(factor(Subject)/(Blocknumber*factor(ContextType))), data = Gen_data)
summary(Model_generalization)

# Error: factor(Subject)
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 36  34.96  0.9712               
# 
# Error: factor(Subject):Blocknumber
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Blocknumber  1  1.375  1.3750   13.21 0.000864 ***
#   Residuals   36  3.748  0.1041                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(ContextType)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(ContextType)   3  5.867  1.9558   16.28 8.45e-09 ***
#   Residuals           108 12.974  0.1201                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(ContextType)
# Df Sum Sq Mean Sq F value Pr(>F)
# Blocknumber:factor(ContextType)   3  0.067 0.02223   0.333  0.801
# Residuals                       108  7.204 0.06670               
# 
# Error: Within
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 148  7.438 0.05025   

eta_squared(Model_generalization)
# # Effect Size for ANOVA (Type I)
# 
# Group                                           |                       Parameter | Eta2 (partial) |       95% CI
# -----------------------------------------------------------------------------------------------------------------
#   factor(Subject):Blocknumber                     |                     Blocknumber |           0.27 | [0.08, 1.00]
# factor(Subject):factor(ContextType)             |             factor(ContextType) |           0.31 | [0.19, 1.00]
# factor(Subject):Blocknumber:factor(ContextType) | Blocknumber:factor(ContextType) |       9.17e-03 | [0.00, 1.00]

#Pairwise t-tests across contexttypes
Gen_type_dat <- aggregate(Gen_data$Errorscore, list(Gen_data$Subject, Gen_data$ContextType), FUN = mean)
colnames(Gen_type_dat)<-c("Subject", "ContextType", "Errorscore")
pairwise.t.test(Gen_type_dat$Errorscore, Gen_type_dat$ContextType, paired= TRUE)
# Pairwise comparisons using paired t tests 
# 
# data:  Gen_type_dat$Errorscore and Gen_type_dat$ContextType 
# 
# Anti    Compositional Full   
# Compositional 0.75402 -             -      
#   Full          0.05096 0.12295       -      
#   Train         1.9e-05 1.9e-05       0.00045
# 
# P value adjustment method: holm 

#Here, again to check t-values
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)

cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)

#Check mean and sd for each contexttype
mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"])
#0.5019641
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"])
#0.3505684

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"])
#0.6825204
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"])
#0.3407935

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"])
#0.789924
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"])
#0.3314545

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"])
#0.7764952
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"])
#0.3081388

#check for each context whether they were better than baseline (below 1)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Train"] -     1
# t = -8.6415, df = 36, p-value = 1.319e-10
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.4007341
# sample estimates:
#   mean of x 
# -0.4980359 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Full"] -     1
# t = -5.6666, df = 36, p-value = 9.695e-07
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.2228908
# sample estimates:
#   mean of x 
# -0.3174796 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Anti"] -     1
# t = -3.8553, df = 36, p-value = 0.0002295
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.1180793
# sample estimates:
#   mean of x 
# -0.210076 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"]-1, alternative = "less")

# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Compositional"] -     1
# t = -4.4121, df = 36, p-value = 4.456e-05
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.1379795
# sample estimates:
#   mean of x 
# -0.2235048

#Make plot of generalization performance
Generalisation_toplot<- summarySE(data = Gen_data, measurevar = "Errorscore", groupvars = c("ContextType", "Blocknumber"))

Generalisation_plot <- ggplot(Generalisation_toplot, aes(x = Blocknumber, y = Errorscore, color = factor(ContextType))) +
  geom_point(size = 4, shape = 16, position = position_dodge(width = .5)) +
  geom_hline(yintercept = 1, linetype="dashed", color = "black", size=1.5)+
  geom_errorbar(aes(ymin = Errorscore- ci, ymax = Errorscore + ci), width = 0.2, position = position_dodge(width = .5)) +
  scale_color_manual(labels = c("Anti", "Compositional", "Full", "Trained"),values = c("Anti" = "magenta", "Full" = "cyan", "Train"= "gold", "Compositional"="green")) +
  scale_x_continuous(breaks = seq(7,13,3))+
  ylim(c(0,1))+
  theme_classic() +
  labs(x = "Round number",
       y = "Baselined errorscore") +
  guides(color = guide_legend("Type"), shape = guide_legend("Type"))

Generalisation_plot
ggsave("Generalisation_plot_Exp2.jpg", Generalisation_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Now check for compositions: we divide them in their components
Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(4,7)) +(Generalization_test_data$locationID %in% c(0,1,6)))==2] <-"A"
Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(4,7)) +(Generalization_test_data$locationID %in% c(2,3)))==2] <-"C"
Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(4,7)) +(Generalization_test_data$locationID %in% c(4,5)))==2] <-"-A"
Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(4,7)) +(Generalization_test_data$locationID %in% c(7)))==2] <-"-C"

Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(5,6)) +(Generalization_test_data$locationID %in% c(0,1,4)))==2] <-"B"
Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(5,6)) +(Generalization_test_data$locationID %in% c(2,3)))==2] <-"D"
Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(5,6)) +(Generalization_test_data$locationID %in% c(5,6)))==2] <-"-B"
Generalization_test_data$Divide[((Generalization_test_data$StimID %in% c(5,6)) +(Generalization_test_data$locationID %in% c(7)))==2] <-"-D"

#aggregate data
Gen_specific_data <- aggregate(list(Generalization_test_data$BaselinedScore), list(Generalization_test_data$Subject,Generalization_test_data$round, Generalization_test_data$ContextType, Generalization_test_data$Divide), FUN = mean)
colnames(Gen_specific_data)<- c("Subject", "Blocknumber", "ContextType", "Rule", "Errorscore")

#we only look at blocks that required generalization so we exclude trained contexts
Gen_specific_data<-Gen_specific_data[Gen_specific_data$ContextType!="Trained",]

#and now aggregate per rule
Gen_rule_data<- aggregate(list(Gen_specific_data$Errorscore), list(Gen_specific_data$Subject,Gen_specific_data$Rule), FUN = mean)
colnames(Gen_rule_data)<- c("Subject", "Rule", "Errorscore")

pairwise.t.test(Gen_rule_data$Errorscore, Gen_rule_data$Rule, paired = TRUE)
#-D and -A were slightly worse than the others

# Pairwise comparisons using paired t tests 
# 
# data:  Gen_rule_data$Errorscore and Gen_rule_data$Rule 
# 
#         -A      -B      -C      -D      A       B       C      
# -B 1.00000 -       -       -       -       -       -      
#   -C 1.00000 1.00000 -       -       -       -       -      
#   -D 1.00000 0.31258 0.39105 -       -       -       -      
#   A  0.43083 1.00000 1.00000 0.01226 -       -       -      
#   B  0.01249 0.03639 0.06155 0.00121 0.43083 -       -      
#   C  0.06672 0.25246 0.16598 0.00183 1.00000 1.00000 -      
#   D  0.05357 0.25246 0.22451 0.00079 0.39105 1.00000 1.00000
# 
# P value adjustment method: holm 

Generalisation_specific_toplot<- summarySE(data = Gen_rule_data, measurevar = "Errorscore", groupvars = c("Rule"))
Generalisation_specific_plot <- ggplot(Gen_rule_data, aes(x = Rule, y = Errorscore) )+
  geom_violin(fill = "grey") +
  geom_boxplot(width = .2) +
  ylim(c(0,1.5))+
  theme_classic() +
  labs(x = "Rule",
       y = "Baselined Errorscore")

Generalisation_specific_plot
ggsave("Generalisation_rule_plot_Exp2.jpg", Generalisation_specific_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Now prepare the additional anova
Gen_specific_data$inversion<-FALSE
Gen_specific_data$inversion[Gen_specific_data$Rule %in% c("-A", "-B","-C", "-D")]<-TRUE
Gen_specific_data$comp<-FALSE
Gen_specific_data$comp[Gen_specific_data$ContextType=="Compositional"]<-TRUE

aov_comp_data<- aggregate(list(Gen_specific_data$Errorscore), list(Gen_specific_data$Subject,Gen_specific_data$inversion, Gen_specific_data$comp), FUN = mean)
colnames(aov_comp_data)<- c("Subject", "Inversion", "Compositionality","Errorscore")

Model_comp <- aov(Errorscore ~ factor(Inversion)*factor(Compositionality)+Error(factor(Subject)/(factor(Inversion)*factor(Compositionality))), data = aov_comp_data)
summary(Model_comp)

# Error: factor(Subject)
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 36  11.36  0.3156               
# 
# Error: factor(Subject):factor(Inversion)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(Inversion)  1 0.6008  0.6008   13.85 0.000674 ***
#   Residuals         36 1.5615  0.0434                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Compositionality)
# Df Sum Sq Mean Sq F value Pr(>F)  
# factor(Compositionality)  1 0.2699 0.26992   5.053 0.0308 *
#   Residuals                36 1.9231 0.05342                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Inversion):factor(Compositionality)
# Df Sum Sq Mean Sq F value Pr(>F)  
# factor(Inversion):factor(Compositionality)  1 0.1826 0.18264   5.621 0.0232 *
#   Residuals                                  36 1.1697 0.03249                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

eta_squared(Model_comp)
# # Effect Size for ANOVA (Type I)
# 
# Group                                                      |                                  Parameter | Eta2 (partial) |       95% CI
# ---------------------------------------------------------------------------------------------------------------------------------------
#   factor(Subject):factor(Inversion)                          |                          factor(Inversion) |           0.28 | [0.09, 1.00]
# factor(Subject):factor(Compositionality)                   |                   factor(Compositionality) |           0.12 | [0.01, 1.00]
# factor(Subject):factor(Inversion):factor(Compositionality) | factor(Inversion):factor(Compositionality) |           0.14 | [0.01, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].

#prepare and make plot
Generalisation_comp_toplot<- summarySE(data = aov_comp_data, measurevar = "Errorscore", groupvars = c("Inversion", "Compositionality"))

Generalisation_comp_plot <- ggplot(Generalisation_comp_toplot, aes(x = Inversion, y = Errorscore, color = factor(Compositionality))) +
  geom_point(size = 4, shape = 16, position = position_dodge(width = .5)) +
  geom_hline(yintercept = 1, linetype="dashed", color = "black", size=1.5)+
  geom_errorbar(aes(ymin = Errorscore- ci, ymax = Errorscore + ci), width = 0.2, position = position_dodge(width = .5)) +
  scale_color_manual(labels = c("No", "Yes"), values = c("FALSE" = "darkgreen", "TRUE" = "lightgreen")) +
  scale_x_discrete(labels = c("No", "Yes"))+
  ylim(c(0,1))+
  theme_classic() +
  labs(x = "Inversion",
       y = "Baselined errorscore") +
  guides(color = guide_legend("Compositionality"))

Generalisation_comp_plot
ggsave("comp_inverse_plot_Exp2.jpg", Generalisation_comp_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Check in each block what they learned again
CDat_try <- aggregate(x = Generalization_test_data$BaselinedScore, by = list(Generalization_test_data$Subject, Generalization_test_data$round, Generalization_test_data$ContextType), FUN = mean)
colnames(CDat_try) <- c("Subject", "Round", "Type", "BaselinedScore")
CDat_try$Rand <- CDat_try$BaselinedScore < 1

CDat_try$RandNone <- FALSE

CDat_try$RandTrained <- FALSE
CDat_try$RandFull <- FALSE
CDat_try$RandAnti <- FALSE
CDat_try$RandComp <- FALSE

CDat_try$RandCombined <- FALSE
CDat_try$RandTFull <- FALSE
CDat_try$RandTAnti <- FALSE
CDat_try$RandTComp <- FALSE
CDat_try$RandAComp <- FALSE
CDat_try$RandFComp <- FALSE

CDat_try$RandTFComp <- FALSE
CDat_try$RandTAComp <- FALSE
CDat_try$RandAFComp <- FALSE
CDat_try$RandTAF <- FALSE

CDat_try$RandAll <- FALSE

for (s in Exp_2_selected_subs){
  d <- CDat_try[CDat_try$Subject == s,]
  for (r in unique(CDat_try$Round)){
    d2 <- d[d$Round == r,]
    
    if (sum(d2$Rand)==4){
      CDat_try$RandAll[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Rand)==0){
      CDat_try$RandNone[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Rand)==1){
      if (d2$Rand[d2$Type=="Anti"]){
        CDat_try$RandAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]){
        CDat_try$RandFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]){
        CDat_try$RandTrained[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Compositional"]){
        CDat_try$RandComp[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }else if (sum(d2$Rand)==2){
      if (d2$Rand[d2$Type=="Anti"]==FALSE & d2$Rand[d2$Type=="Compositional"]==FALSE){
        CDat_try$RandTFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]==FALSE & d2$Rand[d2$Type=="Compositional"]==FALSE){
        CDat_try$RandTAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]==FALSE & d2$Rand[d2$Type=="Compositional"]==FALSE){
        CDat_try$RandCombined[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]==FALSE & d2$Rand[d2$Type=="Anti"]==FALSE){
        CDat_try$RandTComp[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]==FALSE & d2$Rand[d2$Type=="Anti"]==FALSE){
          CDat_try$RandFComp[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      } else if (d2$Rand[d2$Type=="Train"]==FALSE & d2$Rand[d2$Type=="Full"]==FALSE){
        CDat_try$RandAComp[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }else if (sum(d2$Rand)==3){
      if (d2$Rand[d2$Type=="Anti"]==FALSE){
        CDat_try$RandTFComp[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]==FALSE){
        CDat_try$RandTAComp[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]==FALSE){
        CDat_try$RandCombined[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Compositional"]==FALSE){
        CDat_try$RandTAF[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }
  }
}

#aggregate
aggregate_Rand <- aggregate(list(CDat_try$RandNone, CDat_try$RandTrained, CDat_try$RandFull, CDat_try$RandAnti, CDat_try$RandComp, CDat_try$RandCombined, CDat_try$RandTFull, CDat_try$RandTAnti, CDat_try$RandTComp,CDat_try$RandFComp, CDat_try$RandAComp,CDat_try$RandTAF,CDat_try$RandTFComp,CDat_try$RandTAComp,CDat_try$RandAFComp,CDat_try$RandAll), list(CDat_try$Round), FUN = sum )
colnames(aggregate_Rand) <- c("Block", "None", "OnlyTrain", "OnlyFull", "OnlyAnti", "OnlyComp", "AntiFull", "TrainedFull", "TrainedAnti", "TrainedComp", "FullComp", "AntiComp","TrainedAntiFull", "TrainedFullComp", "TrainedAntiComp","AntiFullComp","All")
aggregate_Rand[,2:17] <-aggregate_Rand[,2:17]/4

#check different counts for full
aggregate_Rand$full <- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull+ aggregate_Rand$AntiFull+ aggregate_Rand$FullComp+ aggregate_Rand$TrainedAntiFull+ aggregate_Rand$TrainedFullComp + aggregate_Rand$AntiFullComp + aggregate_Rand$All
aggregate_Rand$fullnoanti <- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull+ aggregate_Rand$FullComp+ aggregate_Rand$TrainedFullComp
aggregate_Rand$fullnocomp <- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull+ aggregate_Rand$AntiFull+ aggregate_Rand$TrainedAntiFull
aggregate_Rand$onlyfull<- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull

#and their proportions
sum(aggregate_Rand$onlyfull)/(37*3)
sum(aggregate_Rand$fullnoanti)/sum(aggregate_Rand$full)
sum(aggregate_Rand$fullnocomp)/sum(aggregate_Rand$full)

#check different counts for anti
aggregate_Rand$anti <- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti+ aggregate_Rand$AntiFull+ aggregate_Rand$AntiComp+ aggregate_Rand$TrainedAntiFull+ aggregate_Rand$TrainedAntiComp + aggregate_Rand$AntiFullComp + aggregate_Rand$All
aggregate_Rand$antinofull <- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti+ aggregate_Rand$AntiComp+ aggregate_Rand$TrainedAntiComp
aggregate_Rand$antinocomp <- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti+ aggregate_Rand$AntiFull+ aggregate_Rand$TrainedAntiFull
aggregate_Rand$onlyanti<- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti

#and their proportions
sum(aggregate_Rand$onlyanti)/(37*3)
sum(aggregate_Rand$antinofull)/sum(aggregate_Rand$anti)
sum(aggregate_Rand$antinocomp)/sum(aggregate_Rand$anti)

#check different counts for compositional
aggregate_Rand$comp <- aggregate_Rand$OnlyComp + aggregate_Rand$TrainedComp+ aggregate_Rand$FullComp+ aggregate_Rand$AntiComp+ aggregate_Rand$TrainedFullComp+ aggregate_Rand$TrainedAntiComp + aggregate_Rand$AntiFullComp + aggregate_Rand$All
aggregate_Rand$compnofull <- aggregate_Rand$OnlyComp + aggregate_Rand$TrainedComp+ aggregate_Rand$AntiComp+ aggregate_Rand$TrainedAntiComp
aggregate_Rand$compnoanti <- aggregate_Rand$OnlyComp + aggregate_Rand$TrainedComp+ aggregate_Rand$FullComp+ aggregate_Rand$TrainedFullComp
aggregate_Rand$onlycomp<- aggregate_Rand$OnlyComp + aggregate_Rand$TrainedComp

#and their proportions
sum(aggregate_Rand$onlycomp)/(37*3)
sum(aggregate_Rand$compnofull)/sum(aggregate_Rand$comp)
sum(aggregate_Rand$compnoanti)/sum(aggregate_Rand$comp)

#Do proportional tests
prop.test(x = c(sum(aggregate_Rand$fullnoanti), sum(aggregate_Rand$antinofull)), n = c(sum(aggregate_Rand$full), sum(aggregate_Rand$anti)))
prop.test(x = c(sum(aggregate_Rand$fullnocomp), sum(aggregate_Rand$compnofull)), n = c(sum(aggregate_Rand$full), sum(aggregate_Rand$comp)))
prop.test(x = c(sum(aggregate_Rand$antinocomp), sum(aggregate_Rand$compnoanti)), n = c(sum(aggregate_Rand$anti), sum(aggregate_Rand$comp)))

#make plots
Rand_long <- gather(aggregate_Rand, Type, Proportion, None:All, factor_key = TRUE)
plot_learnseq_random <-ggplot(data = Rand_long, aes(x = Block, y = (Proportion/37)*100, fill = Type, group = Type))+
  geom_bar(stat="identity", position = "dodge")+
  #scale_fill_manual("Learned Contexts", values=c("None" = "black", "OnlyTrain"="gold", "OnlyFull"= "cyan", "OnlyAnti"="magenta", "AntiFull" = "purple", "TrainedFull"= "orange", "TrainedAnti"= "green", "All"= "red"),labels=c("None" = "None", "OnlyTrain"="Only Trained", "OnlyFull"= "Only Full", "OnlyAnti"="Only Anti", "AntiFull" = "Anti & Full", "TrainedFull"= "Trained & Full", "TrainedAnti"= "Trained & Anti", "All"= "All"))+
  scale_y_continuous("Proportion subjects below baseline")+
  scale_x_continuous("Round", breaks = c(0:5), labels = c(0:5))+
  theme_classic()+
  theme(text = element_text(size=10))

plot_learnseq_random 
ggsave("Learning_Baseline_Exp2.jpg", plot_learnseq_random, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#now we will do this test only within the compositional contexts
comp_test<- Gen_specific_data[Gen_specific_data$comp,]
for (bl in c(7, 10, 13)){
  for (sub in Exp_2_selected_subs){
    comp_test$negscore[comp_test$Blocknumber==bl & comp_test$Subject==sub] <- mean(comp_test$Errorscore[comp_test$Blocknumber==bl & comp_test$Subject==sub][3:4])
    comp_test$posscore[comp_test$Blocknumber==bl & comp_test$Subject==sub] <- mean(comp_test$Errorscore[comp_test$Blocknumber==bl & comp_test$Subject==sub][1:2])
  }
}
comp_test$passed_full<- comp_test$posscore<1
comp_test$passed_anti<- comp_test$negscore<1

comp_test$None <- ((comp_test$passed_full + comp_test$passed_anti)==0)*1
comp_test$All <- ((comp_test$passed_full + comp_test$passed_anti)==2)*1
comp_test$Anti <- ((comp_test$passed_anti + comp_test$All)==1)*1
comp_test$Full <- ((comp_test$passed_full + comp_test$All)==1)*1

aggregate_comp <- aggregate(list(comp_test$None, comp_test$Full, comp_test$Anti, comp_test$All), list(comp_test$Blocknumber), FUN = sum )
colnames(aggregate_comp) <- c("Block", "None", "Full", "Anti", "All")
aggregate_comp[,2:5] <-aggregate_comp[,2:5]/4

aggregate_comp$fullLearned <- aggregate_comp$Full + aggregate_comp$All 
aggregate_comp$antiLearned <- aggregate_comp$Anti + aggregate_comp$All 

#Do proportional test to see whether it is more likely to have learned full without anti than inverse
prop.test(x = c(sum(aggregate_comp$Full), sum(aggregate_comp$Anti)), n = c(sum(aggregate_comp$fullLearned), sum(aggregate_comp$antiLearned)))

#make plot on this
comp_long <- gather(aggregate_comp, Type, Proportion, None:All, factor_key = TRUE)
plot_learnseq_comp <-ggplot(data = comp_long, aes(x = Block, y = (Proportion/37)*100, fill = Type, group = Type))+
  geom_bar(stat="identity", position = "dodge")+
  scale_fill_manual("Learned", values = c("None" = "black", "Full"="cyan","Anti"="magenta", "All" = "red"), labels=c("None" = "None", "Full"="Only full","Anti"="Only Anti", "All" = "Both"))+
  scale_y_continuous("Proportion subjects below random baseline")+
  scale_x_continuous("Round number", breaks = c(7, 10, 13), labels = c(7, 10, 13))+
  theme_classic()+
  theme(text = element_text(size=10))+
  guides(fill = guide_legend("Learned rules"))

plot_learnseq_comp
ggsave("Learning_Baseline_Comp_Exp2.jpg", plot_learnseq_comp, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#combine all in one plot
Exp2_plot<-ggarrange(Block_plot, Generalisation_plot, Generalisation_comp_plot, plot_learnseq_comp, nrow = 2, ncol = 2)
ggsave("All_Exp2.jpg", Exp2_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 20, height = 12, units = "cm")

Exp2_plot_bis<-ggarrange(Block_plot, Generalisation_plot, Generalisation_comp_plot, nrow = 3)
ggsave("All_Exp2_bis.jpg", Exp2_plot_bis, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 15, units = "cm")
#############################################################################
#################################  Experiment 3 #############################
#############################################################################
#check subjects
Exp_3_allSubs <- unique(c(Subs_Exp_3_Sess_1,Subs_Exp_3_Sess_2))
Exp_3_selected_subs <- intersect(Subs_Exp_3_Sess_1,Subs_Exp_3_Sess_2)

removed_subs_Exp3 <- Exp_3_allSubs[! Exp_3_allSubs %in% Exp_3_selected_subs]
#5 got removed because they did not perform both sessions

Exp3_allData <- rbind(Data_Exp_3_Sess_1[Data_Exp_3_Sess_1$Subject %in% Exp_3_selected_subs,], Data_Exp_3_Sess_2[Data_Exp_3_Sess_2$Subject %in% Exp_3_selected_subs,])
#compute errorscore
Exp3_allData$Errorscore <- abs(Exp3_allData$coordID -Exp3_allData$respID)*100

for (iterations in seq(1:1000)){
  randombaseline <- runif(n = nrow(Exp3_allData), min = 0, max = 1)
  random_error1 <- abs(Exp3_allData$coordID - randombaseline)*100
  if (iterations ==1){
    random_error <- random_error1
  }else{
    random_error<-random_error+random_error1
  }
}

Exp3_allData$randombaseline<-random_error/1000

Exp3_allData$BaselinedScore <- Exp3_allData$Errorscore / Exp3_allData$randombaseline

#Subjects that were worse than random in the last test block of the old animals were also removed
round4<- aggregate(Exp3_allData$BaselinedScore[Exp3_allData$round ==4], list(Exp3_allData$Subject[Exp3_allData$round ==4]), FUN = mean)
good_enough <- round4$Group.1[round4$x < 1]
#1 subjects was not better than random in final test block with old animals

Exp3_allData <- Exp3_allData[Exp3_allData$Subject %in% good_enough, ]
Exp_3_selected_subs <- Exp_3_selected_subs[Exp_3_selected_subs%in%good_enough]

#Check performance across blocks/rounds
Blocked_initialaggregate <- aggregate(list(Exp3_allData$BaselinedScore), list(Exp3_allData$Subject, Exp3_allData$round, Exp3_allData$blocktype), FUN = mean)
colnames(Blocked_initialaggregate)<- c("Subject", "Round", "Blocktype", "Errorscore")

Blocked_errorscores <- summarySE(data = Blocked_initialaggregate, measurevar = c("Errorscore"), groupvars = c("Round", "Blocktype"))

#make plot
Block_plot <- ggplot(Blocked_errorscores, aes(x = Round, y = Errorscore, color = factor(Blocktype), shape = factor(Blocktype))) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1, linetype="dashed", color = "black", size=1.5)+
  geom_errorbar(aes(ymin = Errorscore - ci, ymax = Errorscore + ci), width = 0.2) +
  scale_shape_manual(values = c("0" = 4, "1" = 16, "2" = 4, "3"= 16), labels = c("Initial learning", "Initial test", "Generalization learning", "Generalization test")) +
  scale_color_manual(values = c("0" = "blue", "1" = "blue", "2"= "red", "3" = "red"), labels = c("Initial learning", "Initial test", "Generalization learning", "Generalization test")) +
  geom_vline(xintercept = 4.33, linetype="dotted", color = "grey", size=1.5)+
  #geom_point(size = 4, shape = 18, color = "black", aes(y = Blocked_baseline$Baseline, x = Blocked_baseline$Round ))+
  scale_x_continuous(breaks = 1:max(Blocked_errorscores$Round)) +  
  ylim(c(0,1))+
  theme_classic() +
  labs(x = "Round number",
       y = "Baselined errorscore") +
  guides(color = guide_legend("Round"), shape = guide_legend("Round"))

Block_plot
ggsave("Block_plot_Exp3.jpg", Block_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Get all data for Anova
Block_data <- aggregate(list(Exp3_allData$BaselinedScore), list(Exp3_allData$Subject,Exp3_allData$round, Exp3_allData$blocktype), FUN = mean)
colnames(Block_data)<- c("Subject", "Blocknumber", "Blocktype", "Errorscore")
Block_data$Test <- Block_data$Blocktype %%2
Block_data$Animals <- floor(Block_data$Blocktype /2)

#Perform Anova
Model_errorscore <- aov(Errorscore ~ Blocknumber*factor(Test)*factor(Animals)+Error(factor(Subject)/(Blocknumber*factor(Test)*factor(Animals))), data = Block_data)
summary(Model_errorscore)
# Error: factor(Subject)
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 38  14.35  0.3775               
# 
# Error: factor(Subject):Blocknumber
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Blocknumber  1  1.260  1.2600   40.77 1.69e-07 ***
#   Residuals   38  1.174  0.0309                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Test)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(Test)  1 10.366   10.37   115.6 4.38e-13 ***
#   Residuals    38  3.407    0.09                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(Animals)
# Df Sum Sq  Mean Sq F value Pr(>F)
# factor(Animals)  1 0.0001 0.000087    0.01  0.919
# Residuals       38 0.3175 0.008355               
# 
# Error: factor(Subject):Blocknumber:factor(Test)
# Df Sum Sq Mean Sq F value Pr(>F)
# Blocknumber:factor(Test)  1 0.0565 0.05648   2.106  0.155
# Residuals                38 1.0194 0.02683               
# 
# Error: factor(Subject):Blocknumber:factor(Animals)
# Df Sum Sq  Mean Sq F value Pr(>F)
# Blocknumber:factor(Animals)  1 0.0018 0.001838     0.2  0.658
# Residuals                   38 0.3500 0.009211               
# 
# Error: factor(Subject):factor(Test):factor(Animals)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(Test):factor(Animals)  1 0.5996  0.5996   25.45 1.15e-05 ***
#   Residuals                    38 0.8952  0.0236                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(Test):factor(Animals)
# Df Sum Sq Mean Sq F value  Pr(>F)    
# Blocknumber:factor(Test):factor(Animals)  1  1.085  1.0849   38.89 2.7e-07 ***
#   Residuals                                38  1.060  0.0279                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: Within
# Df Sum Sq  Mean Sq F value Pr(>F)
# Residuals 195  1.058 0.005427     

eta_squared(Model_errorscore)
# # Effect Size for ANOVA (Type I)
# 
# Group                                                    |                                Parameter | Eta2 (partial) |       95% CI
# -----------------------------------------------------------------------------------------------------------------------------------
#   factor(Subject):Blocknumber                              |                              Blocknumber |           0.52 | [0.33, 1.00]
# factor(Subject):factor(Test)                             |                             factor(Test) |           0.75 | [0.63, 1.00]
# factor(Subject):factor(Animals)                          |                          factor(Animals) |       2.74e-04 | [0.00, 1.00]
# factor(Subject):Blocknumber:factor(Test)                 |                 Blocknumber:factor(Test) |           0.05 | [0.00, 1.00]
# factor(Subject):Blocknumber:factor(Animals)              |              Blocknumber:factor(Animals) |       5.22e-03 | [0.00, 1.00]
# factor(Subject):factor(Test):factor(Animals)             |             factor(Test):factor(Animals) |           0.40 | [0.20, 1.00]
# factor(Subject):Blocknumber:factor(Test):factor(Animals) | Blocknumber:factor(Test):factor(Animals) |           0.51 | [0.32, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].

# To get more insight in blocktype differences
Test_animals_interactie <- aggregate(list(Exp3_allData$BaselinedScore), list(Exp3_allData$blocktype, Exp3_allData$Subject), FUN = mean)
colnames(Test_animals_interactie) <- c("Blocktype", "Subject","Errorscore")
pairwise.t.test(Test_animals_interactie$Errorscore, Test_animals_interactie$Blocktype, paired = TRUE)

# Pairwise comparisons using paired t tests 
# 
# data:  Test_animals_interactie$Errorscore and Test_animals_interactie$Blocktype 
# 
# 0       1       2      
# 1 1.2e-12 -       -      
#   2 2.1e-11 8.0e-16 -      
#   3 1.5e-06 0.017   1.1e-09
# 
# P value adjustment method: holm 

#Assign contexttype information
Generalization_test_data <- Exp3_allData[Exp3_allData$blocktype==3,]
Generalization_test_data$ContextType[Generalization_test_data$locationID==0]<-"Train"
Generalization_test_data$ContextType[Generalization_test_data$locationID==1]<-"Full"
Generalization_test_data$ContextType[Generalization_test_data$locationID==2]<-"Shrinkage"
Generalization_test_data$ContextType[Generalization_test_data$locationID==3]<-"Train"
Generalization_test_data$ContextType[Generalization_test_data$locationID==4]<-"Compositional"
Generalization_test_data$ContextType[Generalization_test_data$locationID==5]<-"Anti"
Generalization_test_data$ContextType[Generalization_test_data$locationID==6]<-"Full"
Generalization_test_data$ContextType[Generalization_test_data$locationID==7]<-"Expansion"

#aggregate data of generalization test blocks
Gen_data <- aggregate(list(Generalization_test_data$BaselinedScore), list(Generalization_test_data$Subject,Generalization_test_data$round, Generalization_test_data$ContextType), FUN = mean)
colnames(Gen_data)<- c("Subject", "Blocknumber", "ContextType", "Errorscore")

#check whether there were subjects that were below chance level on the trained contexts in last generalization test block
check <- Gen_data[Gen_data$ContextType == "Train",]
check <- check[check$Blocknumber ==13,]
check$Subject[check$Baselined_error< 0]

#Perform anova for contexttypes
Model_generalization <- aov(Errorscore ~ Blocknumber*factor(ContextType)+Error(factor(Subject)/(Blocknumber*factor(ContextType))), data = Gen_data)
summary(Model_generalization)

# Error: factor(Subject)
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 38  57.06   1.502               
# 
# Error: factor(Subject):Blocknumber
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Blocknumber  1  3.387   3.387   18.71 0.000106 ***
#   Residuals   38  6.877   0.181                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):factor(ContextType)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# factor(ContextType)   5   6.66  1.3326   6.854 6.65e-06 ***
#   Residuals           190  36.94  0.1944                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: factor(Subject):Blocknumber:factor(ContextType)
# Df Sum Sq Mean Sq F value Pr(>F)  
# Blocknumber:factor(ContextType)   5  0.911 0.18227   2.047 0.0739 .
# Residuals                       190 16.917 0.08904                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: Within
# Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 234  18.36 0.07846   

eta_squared(Model_generalization)
# # Effect Size for ANOVA (Type I)
# 
# Group                                           |                       Parameter | Eta2 (partial) |       95% CI
# -----------------------------------------------------------------------------------------------------------------
#   factor(Subject):Blocknumber                     |                     Blocknumber |           0.33 | [0.14, 1.00]
# factor(Subject):factor(ContextType)             |             factor(ContextType) |           0.15 | [0.07, 1.00]
# factor(Subject):Blocknumber:factor(ContextType) | Blocknumber:factor(ContextType) |           0.05 | [0.00, 1.00]
#
# - One-sided CIs: upper bound fixed at [1.00].

#Pairwise t-tests across contexttypes
Gen_type_dat <- aggregate(Gen_data$Errorscore, list(Gen_data$Subject, Gen_data$ContextType), FUN = mean)
colnames(Gen_type_dat)<-c("Subject", "ContextType", "Errorscore")
pairwise.t.test(Gen_type_dat$Errorscore, Gen_type_dat$ContextType, paired= TRUE)
# Pairwise comparisons using paired t tests 
# 
# data:  Gen_type_dat$Errorscore and Gen_type_dat$ContextType 
# 
# Anti   Compositional Expansion Full   Shrinkage
# Compositional 0.0422 -             -         -      -        
#   Expansion     1.0000 1.0000        -         -      -        
#   Full          1.0000 0.0073        1.0000    -      -        
#   Shrinkage     1.0000 0.0409        1.0000    1.0000 -        
#   Train         0.0073 3.6e-07       0.0409    0.0214 0.0026   
# 
# P value adjustment method: holm 

#Here, again to check t-values
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"], paired = TRUE)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"], paired = TRUE)

cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"], paired = TRUE)
cohens_d(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"], Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"], paired = TRUE)

#Check mean and sd for each contexttype
mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"])
#0.3320237
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"])
#0.2746905

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"])
#0.4821115
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"])
#0.3273685

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"])
#0.5213985
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"])
#0.34758

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"])
#0.6523572
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"])
#0.2912755

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"])
#0.5742289
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"])
#0.5690515

mean(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"])
#0.5103984
sd(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"])
#0.3353748

#check for each context whether they were better than baseline (below 1)
t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Train"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Train"] -     1
# t = -15.186, df = 38, p-value < 2.2e-16
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.5938185
# sample estimates:
#   mean of x 
# -0.6679763 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Full"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Full"] -     1
# t = -9.8794, df = 38, p-value = 2.389e-12
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.4295093
# sample estimates:
#   mean of x 
# -0.5178885 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Anti"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Anti"] -     1
# t = -8.5991, df = 38, p-value = 9.499e-11
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.3847657
# sample estimates:
#   mean of x 
# -0.4786015 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Compositional"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Compositional"] -     1
# t = -7.4535, df = 38, p-value = 3.019e-09
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.2690075
# sample estimates:
#   mean of x 
# -0.3476428

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Expansion"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Expansion"] -     1
# t = -4.6726, df = 38, p-value = 1.836e-05
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.272145
# sample estimates:
#   mean of x 
# -0.4257711 

t.test(Gen_type_dat$Errorscore[Gen_type_dat$ContextType=="Shrinkage"]-1, alternative = "less")
# One Sample t-test
# 
# data:  Gen_type_dat$Errorscore[Gen_type_dat$ContextType == "Shrinkage"] -     1
# t = -9.1168, df = 38, p-value = 2.09e-11
# alternative hypothesis: true mean is less than 0
# 95 percent confidence interval:
#   -Inf -0.3990608
# sample estimates:
#   mean of x 
# -0.4896016

#Make plot of generalization performance
Generalisation_toplot<- summarySE(data = Gen_data, measurevar = "Errorscore", groupvars = c("ContextType", "Blocknumber"))

Generalisation_plot <- ggplot(Generalisation_toplot, aes(x = Blocknumber, y = Errorscore, color = factor(ContextType))) +
  geom_hline(yintercept = 1, linetype="dashed", color = "black", size=1.5)+
  geom_point(size = 4, shape = 16, position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = Errorscore- ci, ymax = Errorscore + ci), width = 0.2, position = position_dodge(width = .5)) +
  scale_color_manual(labels = c("Anti", "Compositional", "Expansion", "Full", "Shrinkage", "Trained"), values = c("Anti" = "magenta", "Full" = "cyan", "Train"= "gold", "Compositional"="green", "Shrinkage"="burlywood", "Expansion" = "darkorange")) +
  scale_x_continuous(breaks = seq(7,13,3))+
  ylim(c(0,1))+
  theme_classic() +
  labs(x = "Round number",
       y = "Baselined errorscore") +
  guides(color = guide_legend("Type"), shape = guide_legend("Type"))

Generalisation_plot
ggsave("Generalisation_plot_Exp3.jpg", Generalisation_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#Now we create a novel baseline, based on the trained context
for (s in Exp_3_selected_subs){
  d <- Generalization_test_data[Generalization_test_data$Subject == s,]
  for (stim in unique(d$StimID)){
    d$Base2[d$StimID==stim] <- d$coordID[(d$StimID==stim & d$locationID ==3)][1]
  }
  Generalization_test_data$Base2[Generalization_test_data$Subject == s] <-d$Base2
}

#Compute the confusion score
Generalization_test_data$Based2_error <- abs(Generalization_test_data$Base2 - Generalization_test_data$respID)*100

#aggregate data
transform_test <- aggregate(list(Generalization_test_data$Errorscore, Generalization_test_data$Based2_error ), list(Generalization_test_data$Subject, Generalization_test_data$ContextType), FUN = mean)
colnames(transform_test) <- c("Subject", "ContextType", "Errorscore", "Score2")

#Perform t-tests
t.test(transform_test$Errorscore[transform_test$ContextType=="Shrinkage"], transform_test$Score2[transform_test$ContextType=="Shrinkage"], paired = TRUE)
# Paired t-test
# 
# data:  transform_test$Errorscore[transform_test$ContextType == "Shrinkage"] and transform_test$Score2[transform_test$ContextType == "Shrinkage"]
# t = -3.8558, df = 38, p-value = 0.000432
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -11.916663  -3.711542
# sample estimates:
#   mean of the differences 
# -7.814103 

t.test(transform_test$Errorscore[transform_test$ContextType=="Expansion"], transform_test$Score2[transform_test$ContextType=="Expansion"], paired = TRUE)
# Paired t-test
# 
# data:  transform_test$Errorscore[transform_test$ContextType == "Expansion"] and transform_test$Score2[transform_test$ContextType == "Expansion"]
# t = -3.1856, df = 38, p-value = 0.002885
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -30.357858  -6.766073
# sample estimates:
#   mean of the differences 
# -18.56197 

mean(transform_test$Errorscore[transform_test$ContextType=="Shrinkage"])
#13.26389
sd(transform_test$Errorscore[transform_test$ContextType=="Shrinkage"])
#8.709028

mean(transform_test$Errorscore[transform_test$ContextType=="Expansion"])
#23.56517
sd(transform_test$Errorscore[transform_test$ContextType=="Expansion"])
#23.34897

mean(transform_test$Score2[transform_test$ContextType=="Shrinkage"])
#21.07799
sd(transform_test$Score2[transform_test$ContextType=="Shrinkage"])
#6.235115

mean(transform_test$Score2[transform_test$ContextType=="Expansion"])
#42.12714
sd(transform_test$Score2[transform_test$ContextType=="Expansion"])
#14.44108

#now for plotting, the confusion and error scores per block
transform_test_bis <- aggregate(list(Generalization_test_data$Errorscore, Generalization_test_data$Based2_error ), list(Generalization_test_data$Subject, Generalization_test_data$round, Generalization_test_data$ContextType), FUN = mean)
colnames(transform_test_bis) <- c("Subject", "Round", "ContextType", "Errorscore", "Score2")
transform_test_bis<- transform_test_bis[transform_test_bis$ContextType != "Train",]

transform_toplot_confusion<- summarySE(data = transform_test_bis[((transform_test_bis$ContextType=="Shrinkage") + (transform_test_bis$ContextType=="Expansion"))>0,], measurevar = "Score2", groupvars = c("ContextType", "Round"))
transform_toplot_real<- summarySE(data = transform_test_bis[((transform_test_bis$ContextType=="Shrinkage") + (transform_test_bis$ContextType=="Expansion"))>0,], measurevar = "Errorscore", groupvars = c("ContextType", "Round"))

colnames(transform_toplot_confusion) <- c("ContextType", "Round", "N", "Score", "SD", "SE", "ci")
colnames(transform_toplot_real) <- c("ContextType", "Round", "N", "Score", "SD", "SE", "ci")

transform_toplot_confusion$ScoreType <- "Confusion"
transform_toplot_real$ScoreType <- "Error"

transform_toplot<- rbind(transform_toplot_confusion, transform_toplot_real)

transform_plot <- ggplot(transform_toplot, aes(x = Round, y = Score, color = factor(ContextType), fill = factor(ScoreType))) +
  geom_bar(width = .75, position = position_dodge(1.5), stat = "identity", size = 2) +
  geom_errorbar(aes(ymin = Score- ci, ymax = Score + ci), width = 0.75, position = position_dodge(1.5)) +
  scale_color_manual(values = c("Shrinkage"="burlywood", "Expansion" = "darkorange")) +
  scale_fill_manual(values = c("Confusion"= "grey", "Error" = "black")) +
  scale_x_continuous(breaks = seq(7,13,3))+
  #ylim(c(0,1))+
  theme_classic() +
  labs(x = "Round number",
       y = "Score") +
  guides(color = guide_legend("Type"), fill = guide_legend("Score"))

transform_plot
ggsave("Confusion_plot_Exp3.jpg", transform_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 7, units = "cm")

#make combined plot
Exp3_plot<-ggarrange(Block_plot, Generalisation_plot, transform_plot, nrow = 3)
ggsave("All_Exp3.jpg", Exp3_plot, device = "jpeg", path = Plot_path, dpi = 300, width = 15, height = 15, units = "cm")

#Check in each block what they learned
CDat_try <- aggregate(x = Generalization_test_data$BaselinedScore, by = list(Generalization_test_data$Subject, Generalization_test_data$round, Generalization_test_data$ContextType), FUN = mean)
colnames(CDat_try) <- c("Subject", "Round", "Type", "BaselinedScore")
CDat_try$Rand <- CDat_try$BaselinedScore < 1

CDat_try$RandNone <- FALSE
CDat_try$RandTrained <- FALSE
CDat_try$RandFull <- FALSE
CDat_try$RandAnti <- FALSE

CDat_try$RandCombined <- FALSE
CDat_try$RandTFull <- FALSE
CDat_try$RandTAnti <- FALSE

CDat_try$RandAll <- FALSE

CDat_try <- CDat_try[CDat_try$Type %in% c("Full", "Anti", "Train"),]

for (s in Exp_3_selected_subs){
  d <- CDat_try[CDat_try$Subject == s,]
  for (r in unique(CDat_try$Round)){
    d2 <- d[d$Round == r,]
    
    if (sum(d2$Rand)==3){
      CDat_try$RandAll[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Rand)==0){
      CDat_try$RandNone[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
    }else if (sum(d2$Rand)==1){
      if (d2$Rand[d2$Type=="Anti"]){
        CDat_try$RandAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]){
        CDat_try$RandFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]){
        CDat_try$RandTrained[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }else if (sum(d2$Rand)==2){
      if (d2$Rand[d2$Type=="Anti"]==FALSE){
        CDat_try$RandTFull[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Full"]==FALSE){
        CDat_try$RandTAnti[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }else if (d2$Rand[d2$Type=="Train"]==FALSE){
        CDat_try$RandCombined[(CDat_try$Subject == s & CDat_try$Round == r)] <- TRUE
      }
    }
  }
}

#aggregate
aggregate_Rand <- aggregate(list(CDat_try$RandNone, CDat_try$RandTrained, CDat_try$RandFull, CDat_try$RandAnti, CDat_try$RandCombined, CDat_try$RandTFull, CDat_try$RandTAnti, CDat_try$RandAll), list(CDat_try$Round), FUN = sum )
colnames(aggregate_Rand) <- c("Block", "None", "OnlyTrain", "OnlyFull", "OnlyAnti", "AntiFull", "TrainedFull", "TrainedAnti", "All")
aggregate_Rand[,2:9] <-aggregate_Rand[,2:9]/3

#see who learned full and full without anti
aggregate_Rand$full <- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull + aggregate_Rand$AntiFull + aggregate_Rand$All
aggregate_Rand$fullnoanti <- aggregate_Rand$OnlyFull + aggregate_Rand$TrainedFull

#compute proportions
sum(aggregate_Rand$fullnoanti)/(39*3)
sum(aggregate_Rand$fullnoanti)/sum(aggregate_Rand$full)

#see who learned anti and anti without full
aggregate_Rand$anti <- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti + + aggregate_Rand$AntiFull + aggregate_Rand$All
aggregate_Rand$antinofull <- aggregate_Rand$OnlyAnti + aggregate_Rand$TrainedAnti

#compute proportions
sum(aggregate_Rand$antinofull)/(39*3)
sum(aggregate_Rand$antinofull)/sum(aggregate_Rand$anti)

#Do proportional test to see whether it is more likely to have learned full without anti than inverse
prop.test(x = c(sum(aggregate_Rand$fullnoanti), sum(aggregate_Rand$antinofull)), n = c(sum(aggregate_Rand$full), sum(aggregate_Rand$anti)))
# 
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(sum(aggregate_Rand$fullnoanti), sum(aggregate_Rand$antinofull)) out of c(sum(aggregate_Rand$full), sum(aggregate_Rand$anti))
# X-squared = 0.78799, df = 1, p-value = 0.3747
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.05025736  0.15475140
# sample estimates:
#   prop 1    prop 2 
# 0.1553398 0.1030928 
