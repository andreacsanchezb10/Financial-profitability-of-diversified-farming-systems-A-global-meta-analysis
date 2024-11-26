

### FINANCIAL PROFITABILITY OF DIVERSIFIED FARMING SYSTEMS: A GLOBAL META-ANALYSIS

# Andrea C. Sánchez*; Hannah N. Kamau2; Francesca Grazioli1; Sarah K. Jones1
# 1Bioversity International, Parc Scientifique d'Agropolis II, 34397 Montpellier, France, (andrea.sanchez@cgiar.org); (f.Grazioli@cgiar.org); (s.jones@cgiar.org).
# 2Department of Ecology and Natural Resources Management, Center for Development Research (ZEF), Genscherallee 3-D 53113 Bonn, Germany, (hkamau@uni-bonn.de)
# *Corresponding author: (andrea.sanchez@cgiar.org)

##Library###
library(devtools)
library("magrittr") # for piping (%>%)
library(dplyr) # for data manipulation
library(stringr) # string manipulationstringr
library(tidyr)
library(reshape2)
library(Matrix)
library(data.table)
library(metafor) #meta-analysis##
library(Formula)
library(funModeling)
library(frequency)
library(ggnewscale)
library(magick)
library(rgdal)
library(paletteer)
library(forestplot)
library(rlang)
library(doParallel)
library(parallel)
library(iterators)
library(foreach)
library(ggpubr)
library("viridis")
library("wesanderson")
library(ggsci)
library(tidyverse)
library(grid)
library(MuMIn)
library(magick)
library(rgdal)
library(scales)
library(PerformanceAnalytics)
library(vcd)
library("DescTools")

##Set working directory##
setwd("C:/User/Meta-analysis")

###############################################################################################################################
### Appendix B: Dataset 2 Sheet: Data ##
meta <- read.csv("Dataset_2_Data.csv", header = TRUE,  sep = ",")%>%

#Number of articles included in the original dataset
length(sort(unique(meta$ID)))#119
#Number of effect sizes included in the original dataset
length(meta$ES_ID) #3192
sort(unique(meta$Financial_outcome))

#Number of articles and effect sizes by Financial_measure
meta %>%  group_by(Financial_outcome) %>% summarise(n_distinct(ID)) #Number of articles per financial metric
meta %>%  group_by(Financial_outcome) %>% summarise(n_distinct(ES_ID)) #Number of effect sizes per financial metric

## Filter the the Financial_output equals to "B/C ratio", "gross income", "partial cost", "total cost"
## To calculate LOG RESPONSE RATIO effect size -> zero values have to be replaced by 0.0001
meta_logRR<- meta%>%filter(Financial_outcome != "Net income"& Financial_outcome != "Gross margin")%>%
  filter(!is.na(Financial_outcome))%>%
  mutate(Financial_SD_raw_C = if_else(Financial_SD_raw_C==0,0.001,Financial_SD_raw_C))%>%
  mutate(Financial_SD_raw_T = if_else(Financial_SD_raw_T==0,0.001,Financial_SD_raw_T))%>%
  mutate(Financial_value_raw_C= if_else(Financial_value_raw_C == 0, 0.0001, Financial_value_raw_C))%>%
  mutate(Financial_value_raw_T= if_else(Financial_value_raw_T == 0, 0.0001, Financial_value_raw_T))

## Filter the the Financial_outcome equals to "Net income", "Gross margin"
## To calculate STANDARDIZED MEAN DIFFERENCE effect size (because these metrics have negative values)
meta_SMD <- meta%>%filter(Financial_outcome == "Net income" | Financial_outcome == "Gross margin")


################################ EFFECT SIZE CALCULATION ####################################################################33
##########   LOG RESPONSE RATIO
###Calculate the percentage of effect sizes presented sqrt((N_samples_C*mean_C)/SD_C) values >3 (Hedges et al. 1999)
(length(which(meta$prueba_logRR_C > 3))/nrow(meta))*100 #56.45363
(length(which(meta$prueba_logRR_T > 3))/nrow(meta))*100 #72.86967

###Calculate the effect size using log response ratio Reference: Hedge et al. 1999 (vtype=LS)#####
### see = https://rstudio-pubs-static.s3.amazonaws.com/28456_ea0b1faf0f4645cc8af81d81aaf0c1af.html####
##LRR= mean response ratio, LRR_var= variance###
###Calculate the effect size LOG RESPONSE RATIO for the Financial_output == to "benefit-cost ratio", "gross income", "partial cost", "total cost"
effectsize_logRR <- escalc(measure = "ROM", m1i= Financial_value_raw_T, m2i= Financial_value_raw_C,
                           sd1i= Financial_SD_raw_T,sd2i= Financial_SD_raw_C,
                           n1i= Financial_N_T, n2i= Financial_N_C, 
                           data= meta_logRR, var.names=c("Financial_mean","Financial_var"),vtype="LS",digits=4)%>%
  mutate(Financial_se = sqrt(Financial_var),
         Financial_precision = (1/Financial_se))%>%
  mutate(Financial_percent = (100*(exp(Financial_mean)-1)))%>%
  filter(!is.na(Financial_mean))%>%
  mutate(effect_size= "log_RR")

hist(effectsize_logRR$Financial_mean) ##Frequency of Effect sizes##
addmargins(table(effectsize_logRR$Financial_outcome,effectsize_logRR$Validity_overall))

##########   STANDARDIZED MEAN DIFFERENCE
effectsize_SMD <- escalc(measure = "SMD", m1i= Financial_value_raw_T, m2i= Financial_value_raw_C, sd1i= Financial_SD_raw_T,sd2i= Financial_SD_raw_C,n1i= Financial_N_T, 
                         n2i= Financial_N_C, data= meta_SMD, var.names=c("Financial_mean","Financial_var"),vtype="LS",digits=4)%>%
  filter(!is.na(Financial_mean))%>%
  mutate(Financial_se = sqrt(Financial_var),
         Financial_precision = (1/Financial_se))%>%
  mutate(effect_size="SMD")

hist(effectsize_SMD$Financial_mean) ##Frequency of Effect sizes##
addmargins(table(effectsize_SMD$Financial_outcome,effectsize_SMD$Validity_overall))

##Combine log-RR and SMD databases
effectsize<- rbind(effectsize_SMD, effectsize_logRR)%>%
  mutate(Financial_mean_percentage = if_else(Financial_outcome!="Gross margin"&Financial_outcome!="Net income",(100*(exp(Financial_mean)-1)),Financial_mean),
         Financial_var_percentage = if_else(Financial_outcome!="Gross margin"&Financial_outcome!="Net income",(100*(exp(Financial_var)-1)),Financial_var))%>%
  mutate(ES_ID = rownames(.),
         ES_ID_numeric= as.numeric(ES_ID))
View(effectsize)

####################################################################################################################
#---------------------------##########  META-ANALYSES ################----------------------------------------------------#

##### Equation: Cheung (2014) Formula to calculate the estimate sampling variance (formula 14)####
#b= LRR_var
estimated.sampling.variance.func <- function (b) {  
  result<- ((length(b)-1) * sum(1/b))/ (((sum(1/b))^2)-(sum(1/(b^2))))
  return(result)
}

###############################################################################################################################3
############################## Benefit-cost ratio ######################################################################
########---- INTERCEPT ONLY MODEL -------###########
#Estimate the overall mean effect size by fitting an intercept-only model.
#See:  Lopez-Lopez et al. 2018; Assink and Wibbelink (2016)
#tdist=TRUE = the argument specifying that test statistics and confidence intervals must be based on the t-distribution (Assink and Wibbelink 2016)
#For the methodology read Li, Y., Shi, L., & Roth, D. (1994) to see the problem of not use tdis .
#For the methodology read Knapp, G. & Hartung, J. (2003), they proposed to use tdist.
#LEVEL 2: ES_ID (the variable containing the unique identifiers of all effect sizes in the data set)
#LEVEL 3: ID (article ID)
#RESULTS: k= number of ES comprised; 
#REML(REstricted Maximum Likelihood estimation method): method is superior to other methods (see, Hox, 2010; Viechtbauer, 2005), but has restrictions (see Cheung, 2014; Van den Noortgate et al., 2013).
#RESULTS: sigma^2.1 (estim): variance between ES within studies (level 2) if it is small indicates that the ES are similar within studies (Cheung 2014); 
#sigma^2.2(estim): variance between studies (level 3): if it is large, indicates the population effect sizes vary across Level 3, so study characteristics can be included to explain the heterogeneity at level 3 (Cheung 2014)
#RESULTS: Test for heterogeneity: p-val<0.001 significant variation between all effect sizes in the data set.
#RESULTS: estimate = the overall effect size; se = standard error; tval = t value; pval = p value; ci.lb = lower bound of the confidence in- terval; and ci.ub = upper bound of the confidence interval.
BCR.overall <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist= TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="B/C ratio"))
summary(BCR.overall, digits=3)

#####Heterogeneity of within-study variance (level 2)###
#Build a two-level model without within-study variance.
#If the test results provide support for rejecting the null hypothesis, we can conclude that the fit of the original three-level model is statistically better than the fit of the two-level model, and consequently, that there is significant variability between effect sizes within studies.
#sigma2=c(0,NA) = the argument is taken by the rma.mv function when the user wants to fix a specific variance component to a user-defined value. The first parameter (0) states that the within-study variance is fixed to zero (i.e., no within-study variance will be modeled), and the second parameter (NA) states that the between-study variance is estimated.
#The variance at the first level (sampling variance) was not included in the model, because it is assumed to be known.
BCR.modelnovar2 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          sigma2=c(0,NA), tdist=TRUE, data=effectsize, method="REML",
                          verbose=TRUE, subset=(Financial_outcome=="B/C ratio"))
summary(BCR.modelnovar2)

# Perform a likelihood-ratio-test to determine the significance of the within-study variance (level 2).
#RESULTS: Full= represents the three-level model stored in the object abun.overall; 
#Reduced= represents the two-level model stored in the object abun.modelnovar2
#df = degrees of freedom= The reduced model has one degree less than the full model, since within-study variance is not present in the reduced model;
#LRT = likelihood-ratio test. In this column, the value of the test statistic is presented;
#pval = the two-sided p value of the test statistic
#QE= resembles the test for heterogeneity in all effect sizes in the data set, and the value of the test statistic is given in this column
#CONCLUSION: LRT<pval = We found significant variability between effect sizes within studies
anova(BCR.overall,BCR.modelnovar2) #benefit-cost ratio

###Heterogeneity of between-study variance (level 3)
# Build a two-level model without between-study variance;
#In the last model, between-study variance is not modeled. 
#If the null hypothesis should be rejected based on the test results, we can conclude that the fit of the original three-level model is statistically better than the fit of the two-level model, and consequently, that there is significant variability between studies.
#sigma2=c(NA,0): Since we want to fix the between-study variance to zero and freely estimate the within-study variance
BCR.modelnovar3 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID),
                          sigma2=c(NA,0), tdist=TRUE, data=effectsize, method="REML", 
                          subset=(Financial_outcome=="B/C ratio"))
summary(BCR.modelnovar3)

# Perform a likelihood-ratio-test to determine the significance of the between-study variance.
#CONCLUSION: LRT<pval = We found significant variability between studies
anova(BCR.overall,BCR.modelnovar3) 

###The distribution of the variance over the three levels of the meta-analytic model
#Recall that three different sources of variance are modeled in our meta-analytic model: sampling variance at the first level; within-study variance at the second level; and between-study variance at the third level.
#To determine how much variance can be attributed to differences between effect sizes within studies (level 2) and to differences between studies (level 3), formulas given by Cheung (2014 - formula 14 on page 2015) can be used to determine how the total variance is distributed over the three levels of the meta-analytic model;
#Print the results in percentages on screen.
BCR.effectsize<- effectsize%>%filter(Financial_outcome=="B/C ratio")
BCR.estimated.sampling.variance<-  estimated.sampling.variance.func(BCR.effectsize$Financial_var)

###Each of the three variance components (I2_1, I2_2, I2_3) is divided by the total amount of variance, so that a proportional estimate of each variance component is stored in an object.
###overall$sigma2[1]: refers to the amount of within-study variance in the object BCR.overall 
###overall$sigma2[2]: refers to the amount of between-study variance in the object BCR.overall
###The proportional estimates of the three variance components are multiplied by 100 (%), so that a percentage 
#estimate of each variance component is stored in an object

#Sampling variance (Amount of variance at level 1) 
((BCR.estimated.sampling.variance)/(BCR.overall$sigma2[1]+BCR.overall$sigma2[2]+BCR.estimated.sampling.variance))*100

#Within-study variance (Amount of variance at level 2) 
((BCR.overall$sigma2[1]) / (BCR.overall$sigma2[1] + BCR.overall$sigma2[2] + BCR.estimated.sampling.variance))*100

#Between-study variance (Amount of variance at level 3)
((BCR.overall$sigma2[2]) / (BCR.overall$sigma2[1] + BCR.overall$sigma2[2] + BCR.estimated.sampling.variance))*100

########---- META-REGRESSION MODELS -------###########
###Because the variance within and between studies was substantial, we fit a moderator analysis (meta-regression models) to determine the potential moderating effect of:

#-------------------------##########  META-REGRESSION (MODERATOR: FAO_group_T_recla) ################----------------------------------------------------#
BCR.FAO_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ (FAO_group_T_recla)-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.FAO_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: System_T) ################----------------------------------------------------#
BCR.system_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~System_T -1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                       tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.system_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_duration_T) ################----------------------------------------------------#
BCR.Crop_duration_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_duration_T-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                               tdist=TRUE, data=effectsize, method="REML",
                               subset=(Financial_outcome=="B/C ratio"))
summary(BCR.Crop_duration_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Input_type_CT) ################----------------------------------------------------#
BCR.input_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Input_type_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                       tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.input_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_Development_status) ################----------------------------------------------------#
BCR.UN_development <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_Development_status-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                             tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.UN_development, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Fertiliser_chem_CT) ################----------------------------------------------------#
BCR.fertiliser_chem_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Fertiliser_chem_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                 tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.fertiliser_chem_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Pesticide_CT) ################----------------------------------------------------#
BCR.pesticide_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Pesticide_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.pesticide_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Irrigation_CT) ################----------------------------------------------------#
BCR.irrigation_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Irrigation_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.irrigation_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Soil_management_CT) ################----------------------------------------------------#
BCR.soil_management_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Soil_management_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                              tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.soil_management_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_woodiness_T_recla) ################----------------------------------------------------#
BCR.Crop_woodiness_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_woodiness_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.Crop_woodiness_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_subregion) ################----------------------------------------------------#
BCR.UN_subregion <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_subregion-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.UN_subregion, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: labour_difference_TC) ################----------------------------------------------------#
BCR.Labour_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ labour_difference_TC, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.Labour_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Year_assessment_range) ################----------------------------------------------------#
BCR.Year_assessment <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Year_assessment_range-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.Year_assessment, digits=3)

###########################################################################################################3
############################## Gross income ######################################################################
########---- INTERCEPT ONLY MODEL -------###########
GI.overall <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                     tdist= TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Gross income"))
summary(GI.overall, digits=3)

#####Heterogeneity of within-study variance (level 2)###
GI.modelnovar2 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                         sigma2=c(0,NA), tdist=TRUE, data=effectsize, method="REML", verbose=TRUE, 
                         subset=(Financial_outcome=="Gross income"))
summary(GI.modelnovar2)

# Perform a likelihood-ratio-test to determine the
anova(GI.overall,GI.modelnovar2) 

###Heterogeneity of between-study variance (level 3)
GI.modelnovar3 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID),
                         sigma2=c(NA,0), tdist=TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Gross income"))
summary(GI.modelnovar3)

# Perform a likelihood-ratio-test to determine the significance of the between-study variance.
anova(GI.overall,GI.modelnovar3) 

###The distribution of the variance over the three levels of the meta-analytic model
GI.effectsize<- effectsize%>%filter(Financial_outcome=="Gross income")
GI.estimated.sampling.variance<-  estimated.sampling.variance.func(GI.effectsize$Financial_var)

#Sampling variance (Amount of variance at level 1) 
((GI.estimated.sampling.variance)/(GI.overall$sigma2[1]+GI.overall$sigma2[2]+GI.estimated.sampling.variance))*100

#Within-study variance (Amount of variance at level 2)
((GI.overall$sigma2[1]) / (GI.overall$sigma2[1] + GI.overall$sigma2[2] + GI.estimated.sampling.variance))*100

#Between-study variance (Amount of variance at level 3) 
((GI.overall$sigma2[2]) / (GI.overall$sigma2[1] + GI.overall$sigma2[2] + GI.estimated.sampling.variance))*100

########---- META-REGRESSION MODELS -------###########

#-------------------------##########  META-REGRESSION (MODERATOR: FAO_group_T_recla) ################----------------------------------------------------#
GI.FAO_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ FAO_group_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                     tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.FAO_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: System_T) ################----------------------------------------------------#
GI.system_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Diversified_system-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.system_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_duration_T) ################----------------------------------------------------#
GI.Crop_duration_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_duration_T-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                              tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.Crop_duration_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Input_type_CT) ################----------------------------------------------------#
GI.input_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Input_type_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.input_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_Development_status) ################----------------------------------------------------#
GI.UN_development <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_Development_status-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.UN_development, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Fertiliser_chem_CT) ################----------------------------------------------------#
GI.fertiliser_chem_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Fertiliser_chem_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.fertiliser_chem_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Pesticide_CT) ################----------------------------------------------------#
GI.pesticide_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Pesticide_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.pesticide_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Irrigation_CT) ################----------------------------------------------------#
GI.irrigation_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Irrigation_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.irrigation_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Soil_management_CT) ################----------------------------------------------------#
GI.soil_management_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Soil_management_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                             tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.soil_management_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_woodiness_T_recla) ################----------------------------------------------------#
GI.Crop_woodiness_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_woodiness_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                               tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.Crop_woodiness_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_subregion) ################----------------------------------------------------#
GI.UN_subregion <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_subregion-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.UN_subregion, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: labour_difference_TC) ################----------------------------------------------------#
GI.Labour_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ labour_difference_TC, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.Labour_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Year_assessment) ################----------------------------------------------------#
GI.Year_assessment <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Year_assessment_range-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                              tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.Year_assessment, digits=3)

###########################################################################################################3
############################## Total cost ######################################################################
########---- INTERCEPT ONLY MODEL -------###########
TC.overall <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                     tdist= TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Total cost"))
summary(TC.overall, digits=3)

#####Heterogeneity of within-study variance (level 2)###
TC.modelnovar2 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                         sigma2=c(0,NA), tdist=TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Total cost"))
summary(TC.modelnovar2)

#Perform a likelihood-ratio-test to determine the
anova(TC.overall,TC.modelnovar2) 

###Heterogeneity of between-study variance (level 3)
TC.modelnovar3 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID),
                         sigma2=c(NA,0), tdist=TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Total cost"))
summary(GI.modelnovar3)

# Perform a likelihood-ratio-test to determine the significance of the between-study variance.
anova(TC.overall,TC.modelnovar3) #ABUNDANCE

###The distribution of the variance over the three levels of the meta-analytic model
TC.effectsize<- effectsize%>%filter(Financial_outcome=="Total cost")
TC.estimated.sampling.variance<-  estimated.sampling.variance.func(TC.effectsize$Financial_var)

#Sampling variance (Amount of variance at level 1) 
((TC.estimated.sampling.variance)/(TC.overall$sigma2[1]+TC.overall$sigma2[2]+TC.estimated.sampling.variance))*100

#Within-study variance (Amount of variance at level 2) 
((TC.overall$sigma2[1]) / (TC.overall$sigma2[1] + TC.overall$sigma2[2] + TC.estimated.sampling.variance))*100

#Between-study variance (Amount of variance at level 3) 
((TC.overall$sigma2[2]) / (TC.overall$sigma2[1] + TC.overall$sigma2[2] + TC.estimated.sampling.variance))*100

########---- META-REGRESSION MODELS -------###########

#-------------------------##########  META-REGRESSION (MODERATOR: FAO_group_T_recla) ################----------------------------------------------------#
TC.FAO_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ FAO_group_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                     tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.FAO_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: System_T) ################----------------------------------------------------#
TC.system_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Diversified_system-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.system_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_duration_T) ################----------------------------------------------------#
TC.Crop_duration_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_duration_T-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                              tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.Crop_duration_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Input_type_CT) ################----------------------------------------------------#
TC.input_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Input_type_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.input_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_Development_status) ################----------------------------------------------------#
TC.UN_development <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_Development_status-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.UN_development, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Fertiliser_chem_CT) ################----------------------------------------------------#
TC.fertiliser_chem_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Fertiliser_chem_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.fertiliser_chem_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Pesticide_CT) ################----------------------------------------------------#
TC.pesticide_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Pesticide_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.pesticide_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Irrigation_CT) ################----------------------------------------------------#
TC.irrigation_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Irrigation_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.irrigation_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Soil_management_CT) ################----------------------------------------------------#
addmargins(table(TC.effectsize$Soil_management_CT))
TC.soil_management_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Soil_management_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                             tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.soil_management_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_woodiness_T_recla) ################----------------------------------------------------#
TC.Crop_woodiness_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_woodiness_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.Crop_woodiness_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_subregion) ################----------------------------------------------------#
TC.UN_subregion <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_subregion-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.UN_subregion, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: labour_difference_TC) ################----------------------------------------------------#
TC.Labour_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ labour_difference_TC, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                       tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.Labour_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Year_assessment_range) ################----------------------------------------------------#
addmargins(table(TC.effectsize$Time_state_T))
TC.Year_assessment <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Year_assessment_range-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.Year_assessment, digits=3)

###########################################################################################################3
############################## Net income ######################################################################

########---- INTERCEPT ONLY MODEL -------###########
NI.overall <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                     tdist= TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Net income"))
summary(NI.overall, digits=3)

#####Heterogeneity of within-study variance (level 2)###
NI.modelnovar2 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                         sigma2=c(0,NA), tdist=TRUE, data=effectsize, method="REML", verbose=TRUE, subset=(Financial_outcome=="Net income"))
summary(NI.modelnovar2)

# Perform a likelihood-ratio-test to determine the
anova(NI.overall,NI.modelnovar2) #B/C ratio

###Heterogeneity of between-study variance (level 3)
NI.modelnovar3 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID),
                         sigma2=c(NA,0), tdist=TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Net income"))
summary(NI.modelnovar3)

# Perform a likelihood-ratio-test to determine the significance of the between-study variance.
anova(NI.overall,NI.modelnovar3) #ABUNDANCE

###The distribution of the variance over the three levels of the meta-analytic model
NI.effectsize<- effectsize%>%filter(Financial_outcome=="Net income")
NI.estimated.sampling.variance<-  estimated.sampling.variance.func(NI.effectsize$Financial_var)

#Sampling variance (Amount of variance at level 1) #11.45035
((NI.estimated.sampling.variance)/(NI.overall$sigma2[1]+NI.overall$sigma2[2]+NI.estimated.sampling.variance))*100

#Within-study variance (Amount of variance at level 2) #19.42363
((NI.overall$sigma2[1]) / (NI.overall$sigma2[1] + NI.overall$sigma2[2] + NI.estimated.sampling.variance))*100

#Between-study variance (Amount of variance at level 3) 69.12602
((NI.overall$sigma2[2]) / (NI.overall$sigma2[1] + NI.overall$sigma2[2] + NI.estimated.sampling.variance))*100

########---- META-REGRESSION MODELS -------###########

#-------------------------##########  META-REGRESSION (MODERATOR: FAO_group_T_recla) ################----------------------------------------------------#
NI.FAO_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ FAO_group_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                     tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.FAO_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: System_T) ################----------------------------------------------------#
NI.system_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ System_T-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.system_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_duration_T) ################----------------------------------------------------#
NI.Crop_duration_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_duration_T-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                              tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.Crop_duration_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Input_type_CT) ################----------------------------------------------------#
NI.input_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Input_type_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.input_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_Development_status) ################----------------------------------------------------#
NI.UN_development <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_Development_status-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.UN_development, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Fertiliser_chem_CT) ################----------------------------------------------------#
NI.fertiliser_chem_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Fertiliser_chem_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.fertiliser_chem_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Pesticide_CT) ################----------------------------------------------------#
NI.pesticide_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Pesticide_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.pesticide_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Irrigation_CT) ################----------------------------------------------------#
NI.irrigation_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Irrigation_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.irrigation_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Soil_management_CT) ################----------------------------------------------------#
NI.soil_management_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Soil_management_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                             tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.soil_management_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_woodiness_T) ################----------------------------------------------------#
NI.Crop_woodiness_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_woodiness_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                               tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.Crop_woodiness_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_subregion) ################----------------------------------------------------#
NI.UN_subregion <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_subregion-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.UN_subregion, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: labour_difference_TC) ################----------------------------------------------------#
NI.Labour_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ labour_difference_TC, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.Labour_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Year_assessment_range) ################----------------------------------------------------#
NI.Year_assessment <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Year_assessment_range-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.Year_assessment, digits=3)

##########################################################################################################3
############################## Gross margin ######################################################################
########---- INTERCEPT ONLY MODEL -------###########
GM.overall <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                     tdist= TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Gross margin"))
summary(GM.overall, digits=3)

#####Heterogeneity of within-study variance (level 2)###
GM.modelnovar2 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                         sigma2=c(0,NA), tdist=TRUE, data=effectsize, method="REML", verbose=TRUE, subset=(Financial_outcome=="Gross margin"))
summary(GM.modelnovar2)

# Perform a likelihood-ratio-test to determine the
anova(GM.overall,GM.modelnovar2)

###Heterogeneity of between-study variance (level 3)
GM.modelnovar3 <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID),
                         sigma2=c(NA,0), tdist=TRUE, data=effectsize, method="REML", subset=(Financial_outcome=="Gross margin"))
summary(GM.modelnovar3)

# Perform a likelihood-ratio-test to determine the significance of the between-study variance.
anova(GM.overall,GM.modelnovar3) #ABUNDANCE

###The distribution of the variance over the three levels of the meta-analytic model
GM.effectsize<- effectsize%>%filter(Financial_outcome=="Gross margin")
GM.estimated.sampling.variance<-  estimated.sampling.variance.func(GM.effectsize$Financial_var)

#Sampling variance (Amount of variance at level 1) # 6.204258
((GM.estimated.sampling.variance)/(GM.overall$sigma2[1]+GM.overall$sigma2[2]+GM.estimated.sampling.variance))*100

#Within-study variance (Amount of variance at level 2) #39.48765
((GM.overall$sigma2[1]) / (GM.overall$sigma2[1] + GM.overall$sigma2[2] + GM.estimated.sampling.variance))*100

#Between-study variance (Amount of variance at level 3) #49.85954
((GM.overall$sigma2[2]) / (GM.overall$sigma2[1] + GM.overall$sigma2[2] + GM.estimated.sampling.variance))*100

########---- META-REGRESSION MODELS -------###########

#-------------------------##########  META-REGRESSION (MODERATOR: FAO_group_T_recla) ################----------------------------------------------------#
GM.FAO_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ FAO_group_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.FAO_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: "System_T",) ################----------------------------------------------------#
GM.system_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ System_T-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.system_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_duration_CT) ################----------------------------------------------------#
GM.Crop_duration_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_duration_T-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                              tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.Crop_duration_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Input_type_CT) ################----------------------------------------------------#
GM.input_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Input_type_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.input_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_Development_status) ################----------------------------------------------------#
GM.UN_development <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_Development_status-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                            tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.UN_development, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Fertiliser_chem_CT) ################----------------------------------------------------#
GM.fertiliser_chem_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Fertiliser_chem_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.fertiliser_chem_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Pesticide_CT) ################----------------------------------------------------#
GM.pesticide_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Pesticide_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.pesticide_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Irrigation_CT) ################----------------------------------------------------#
GM.irrigation_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Irrigation_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                           tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.irrigation_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Soil_management_CT) ################----------------------------------------------------#
GM.soil_management_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Soil_management_CT-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                             tdist=TRUE, data=effectsize_SMD, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.soil_management_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Crop_woodiness_T_recla) ################----------------------------------------------------#
GM.Crop_woodiness_T <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Crop_woodiness_T_recla-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                               tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.Crop_woodiness_T, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: UN_subregion) ################----------------------------------------------------#
addmargins(table(GM.effectsize$UN_subregion))
GM.UN_subregion <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ UN_subregion-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.UN_subregion, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: labour_difference_TC) ################----------------------------------------------------#
GM.Labour_CT <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ labour_difference_TC, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                       tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.Labour_CT, digits=3)

#-------------------------##########  META-REGRESSION (MODERATOR: Year_assessment_range) ################----------------------------------------------------#
addmargins(table(GM.effectsize$Time_state_T))
GM.Year_assessment <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Year_assessment_range-1, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                          tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.Year_assessment, digits=3)

##  MULTIVARIATE META-REGRESSION MODEL
## CHECK FOR MULTICOLLINEARITY BETWEEN EXPLANATORY VARIABLES 
#http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
#Correlation between categorical variables (i.e., Functional_group_recla and min.distance.categorical)
old_column<-c("V1",  "V2",  "V3",  "V4",  "V5",  "V6",  "V7",  "V8",  "V9",  "V10", "V11","V12")
new_column<- c("FAO_group_T_recla_2","Crop_duration_T" ,"Crop_woodiness_T_recla", "Diversified_system",
               "Input_type_CT" , "Fertiliser_chem_CT","Pesticide_CT" ,"Soil_management_CT" ,
               "Irrigation_CT", "UN_Development_status" ,"UN_subregion","Year_assessment_range")
  
cramer_data<- effectsize%>%
  select(Financial_outcome, FAO_group_T_recla_2, Crop_duration_T,Crop_woodiness_T_recla,
         Diversified_system, Input_type_CT, Fertiliser_chem_CT,
         Pesticide_CT, Soil_management_CT, Irrigation_CT, UN_Development_status,UN_subregion,Year_assessment_range)

cor_cramer <- function(x) {
  require(vcd)
  nc <- ncol(x)
  v <- expand.grid(1:nc, 1:nc)
  matrix(mapply(function(i1, i2) assocstats(table(x[,i1],
                                                  x[,i2]))$cramer, v[,1], v[,2]), nc, nc)
  }

BCR_cramer<- cramer_data%>%
  filter(Financial_outcome=="B/C ratio")%>%select(2:13)

BCR_cramer<-  as.data.frame(cor_cramer(BCR_cramer))%>%
  mutate_all(funs(round(., digits =1)))%>%
  setnames(., old_column, new_column)%>%
  mutate(variables = new_column)%>%
  select(13,1:12)
write.csv(BCR_cramer, "C:/Users/Meta-analysis/Table_A6_a.csv")

GI_cramer<- cramer_data%>%
  filter(Financial_outcome=="Gross income")%>%select(2:13)
GI_cramer<-  as.data.frame(cor_cramer(GI_cramer))%>%
  mutate_all(funs(round(., digits =1)))%>%
  setnames(., old_column, new_column)%>%
  mutate(variables = new_column)%>%
  select(13,1:12)
write.csv(GI_cramer, "C:/Users/Meta-analysis/Table_A6_b.csv")

TC_cramer<- cramer_data%>%
  filter(Financial_outcome=="Total cost")%>%select(2:13)
TC_cramer<-  as.data.frame(cor_cramer(TC_cramer))%>%
  mutate_all(funs(round(., digits =1)))%>%
  setnames(., old_column, new_column)%>%
  mutate(variables = new_column)%>%
  select(13,1:12)
write.csv(TC_cramer, "C:/Users/Meta-analysis/Table_A6_c.csv")

GM_cramer<- cramer_data%>%
  filter(Financial_outcome=="Gross margin")%>%select(2:13)
GM_cramer<-  as.data.frame(cor_cramer(GM_cramer))%>%
  mutate_all(funs(round(., digits =1)))%>%
  setnames(., old_column, new_column)%>%
  mutate(variables = new_column)%>%
  select(13,1:12)
write.csv(GM_cramer, "C:/Users/Meta-analysis/Table_A6_d.csv")

NI_cramer<- cramer_data%>%
  filter(Financial_outcome=="Net income")%>%select(2:13)
NI_cramer<-  as.data.frame(cor_cramer(NI_cramer))%>%
  mutate_all(funs(round(., digits =1)))%>%
  setnames(., old_column, new_column)%>%
  mutate(variables = new_column)%>%
  select(13,1:12)
write.csv(NI_cramer, "C:/Users/Meta-analysis/Table_A6_e.csv")

#Correlation between continuous (i.e., labour_difference_TC) and categorical variables (i.e., the other variables)
kruskal_data<- effectsize%>%
  filter(!is.na(labour_difference_TC))%>%
  select(Financial_outcome, FAO_group_T_recla_2, Crop_duration_T,Crop_woodiness_T_recla,
         Diversified_system, Input_type_CT, Fertiliser_chem_CT,
         Pesticide_CT, Soil_management_CT, Irrigation_CT, UN_Development_status,UN_subregion, Year_assessment_range, labour_difference_TC)%>%
  mutate_at(c(2:13),as.factor)

BCR.kruskal_data<- kruskal_data%>%filter(Financial_outcome=="B/C ratio")
GI.kruskal_data<- kruskal_data%>%filter(Financial_outcome=="Gross income")
TC.kruskal_data<- kruskal_data%>%filter(Financial_outcome=="Total cost")
GM.kruskal_data<- kruskal_data%>%filter(Financial_outcome=="Gross margin")
NI.kruskal_data<- kruskal_data%>%filter(Financial_outcome=="Net income")
View(NI.kruskal_data)

kruskal.test(data=GI.kruskal_data,Crop_duration_T~labour_difference_TC)

cor_kruscal <- function (b) {  
  result<- as.data.frame(c(kruskal.test(FAO_group_T_recla_2~labour_difference_TC, data =b )))%>%
    rbind(as.data.frame(c(kruskal.test(Crop_duration_T~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Crop_woodiness_T_recla~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Diversified_system~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Input_type_CT~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Fertiliser_chem_CT~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Pesticide_CT~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Soil_management_CT~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Irrigation_CT~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(UN_Development_status~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(UN_subregion~labour_difference_TC, data =b ))))%>%
    rbind(as.data.frame(c(kruskal.test(Year_assessment_range~labour_difference_TC, data =b ))))
  return(result)
}

kruskal_test<- cor_kruscal(BCR.kruskal_data)%>%
  rbind(cor_kruscal(GI.kruskal_data))%>%
  rbind(cor_kruscal(TC.kruskal_data))%>%
  rbind(cor_kruscal(GM.kruskal_data))%>%
  rbind(cor_kruscal(NI.kruskal_data))%>%
  mutate(Financial_outcome=c(rep("B/C ratio",12),rep("Gross income",12), rep("Total cost",12),
                             rep("Gross margin",12), rep("Net income",12)),
         statistic = round(statistic,2),
         label= paste("H(",parameter, ") = ", statistic, ", p = ",p.value,sep = ""),
         observation = 1:n())
View(kruskal_test)
write.csv(kruskal_test, "C:/Users/Meta-analysis/Table_A6_kruskal_test.csv")

## Multivariate meta-regression: Including no-collinear moderators found to be significant during the univariate meta-regression analysis
#---- B/C ratio
## Model 1
BCR.multivariate_1<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)
                          ,tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize,
                          method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.multivariate_1, digits=5)

summary(BCR.FAO_T, digits=3)

## Model 2
BCR.multivariate_2<- rma.mv(Financial_mean, Financial_var, mods = ~factor(UN_subregion)
                          ,tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize,
                          method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.multivariate_2, digits=5)

## Multi-model inference: Including all no-collinear moderators and check model goodness of fit
#https://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti_and_mumin
eval(metafor:::.MuMIn)

#---- Gross income
## Model 1
GI.multivariate_1<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_duration_T)+factor(System_T)+
                             factor(Fertiliser_chem_CT)+ factor(Pesticide_CT)+
                             factor(Soil_management_CT)+factor(Irrigation_CT)+
                           factor(UN_subregion)+factor(Year_assessment_range),
                         tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",
                         subset=(Financial_outcome=="Gross income"))
summary(GI.multivariate_1, digits=5)

GI.multivariate_1_AICc <- dredge(GI.multivariate_1, trace=2, rank = "AICc")
subset(GI.multivariate_1_AICc, delta <= 2, recalc.weights=FALSE)
write.csv(GI.multivariate_1_AICc, "C:/Users/Andrea/Documents/Bioversity/Cost_benefits_analysis/Meta-analysis/Results_2022.06.21/GI.multivariate_1_AICc.csv")


#Multivariate model omnibus test ARREGLAR

GI.multivariate_1_results<- data.frame(df_1 = c(anova(GI.multivariate_1, btt=2:11)$m,
                                                anova(GI.multivariate_1, btt=12:13)$m,
                                                anova(GI.multivariate_1, btt=14:18)$m,
                                                anova(GI.multivariate_1, btt=19:21)$m,
                                                anova(GI.multivariate_1, btt=22:25)$m,
                                                anova(GI.multivariate_1, btt=26:33)$m,
                                                anova(GI.multivariate_1, btt=34:36)$m,
                                                anova(GI.multivariate_1, btt=37:42)$m,
                                                anova(GI.multivariate_1, btt=43:47)$m),
                                       df_2 = c(anova(GI.multivariate_1, btt=2:11)$k,anova(GI.multivariate_1, btt=12:13)$k,
                                                anova(GI.multivariate_1, btt=14:18)$k,anova(GI.multivariate_1, btt=19:21)$k,
                                                anova(GI.multivariate_1, btt=22:25)$k,anova(GI.multivariate_1, btt=26:33)$k,
                                                anova(GI.multivariate_1, btt=34:36)$k,
                                                anova(GI.multivariate_1, btt=37:42)$k,
                                                anova(GI.multivariate_1, btt=43:47)$k),
                                       QM = c(anova(GI.multivariate_1, btt=2:11)$QM,anova(GI.multivariate_1, btt=12:13)$QM,
                                              anova(GI.multivariate_1, btt=14:18)$QM,anova(GI.multivariate_1, btt=19:21)$QM,
                                              anova(GI.multivariate_1, btt=22:25)$QM,anova(GI.multivariate_1, btt=26:33)$QM,
                                              anova(GI.multivariate_1, btt=34:36)$QM,
                                              anova(GI.multivariate_1, btt=37:42)$QM,
                                              anova(GI.multivariate_1, btt=43:47)$QM),
                                       QMp = c(anova(GI.multivariate_1, btt=2:11)$QMp,anova(GI.multivariate_1, btt=12:13)$QMp,
                                               anova(GI.multivariate_1, btt=14:18)$QMp,anova(GI.multivariate_1, btt=19:21)$QMp,
                                               anova(GI.multivariate_1, btt=22:25)$QMp,anova(GI.multivariate_1, btt=26:33)$QMp,
                                               anova(GI.multivariate_1, btt=34:36)$QMp,
                                               anova(GI.multivariate_1, btt=37:42)$QMp,
                                               anova(GI.multivariate_1, btt=43:47)$QMp),
                                       moderator = c("FAO_crop","Crop_duration","System_T","Fertiliser_chem_CT",
                                                     "Pesticide_CT","Soil_management_CT","Irrigation_CT","UN_subregion","Year_assessment_range" ))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(GI.multivariate_1_results, "C:/Users/Meta-analysis/Table_A10_GI_m_1.csv")

## Model 2
GI.multivariate_2<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_duration_T)+factor(System_T)+
                             factor(Fertiliser_chem_CT)+factor(Pesticide_CT)+
                             factor(Soil_management_CT)+
                             factor(Irrigation_CT)+factor(UN_Development_status),
                         tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.multivariate_2, digits=5)

#Multivariate model omnibus test
GI.multivariate_2_results<- data.frame(df_1 = c(anova(GI.multivariate_2, btt=2:11)$m,
                                                anova(GI.multivariate_2, btt=12:13)$m,
                                                anova(GI.multivariate_2, btt=14:18)$m,
                                                anova(GI.multivariate_2, btt=19:21)$m,
                                                anova(GI.multivariate_2, btt=22:25)$m,
                                                anova(GI.multivariate_2, btt=26:33)$m,
                                                anova(GI.multivariate_2, btt=34:36)$m,
                                                anova(GI.multivariate_2, btt=37:37)$m),
                                       df_2 = c(anova(GI.multivariate_2, btt=2:11)$k,anova(GI.multivariate_2, btt=12:13)$k,
                                                anova(GI.multivariate_2, btt=14:18)$k,anova(GI.multivariate_2, btt=19:21)$k,
                                                anova(GI.multivariate_2, btt=22:25)$k,anova(GI.multivariate_2, btt=26:33)$k,
                                                anova(GI.multivariate_2, btt=34:36)$k,
                                                anova(GI.multivariate_2, btt=37:37)$k),
                                       QM = c(anova(GI.multivariate_2, btt=2:11)$QM,anova(GI.multivariate_2, btt=12:13)$QM,
                                                anova(GI.multivariate_2, btt=14:18)$QM,anova(GI.multivariate_2, btt=19:21)$QM,
                                                anova(GI.multivariate_2, btt=22:25)$QM,anova(GI.multivariate_2, btt=26:33)$QM,
                                                anova(GI.multivariate_2, btt=34:36)$QM,
                                                anova(GI.multivariate_2, btt=37:37)$QM),
                                       QMp = c(anova(GI.multivariate_2, btt=2:11)$QMp,anova(GI.multivariate_2, btt=12:13)$QMp,
                                                anova(GI.multivariate_2, btt=14:18)$QMp,anova(GI.multivariate_2, btt=19:21)$QMp,
                                                anova(GI.multivariate_2, btt=22:25)$QMp,anova(GI.multivariate_2, btt=26:33)$QMp,
                                                anova(GI.multivariate_2, btt=34:36)$QMp,
                                                anova(GI.multivariate_2, btt=37:37)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_duration_T","System_T","Fertiliser_chem_CT",
                                                     "Pesticide_CT","Soil_management_CT","Irrigation_CT","UN_Development_status"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))
  
write.csv(GI.multivariate_2_results, "C:/Users/Meta-analysis/Table_A10_GI_m_2.csv")

## Model 3
GI.multivariate_3<- rma.mv(Financial_mean, Financial_var, mods = ~factor(Crop_woodiness_T_recla)+
                             factor(System_T)+factor(Fertiliser_chem_CT)+
                             factor(Pesticide_CT)+
                             factor(Soil_management_CT)+factor(Irrigation_CT)+
                             factor(UN_subregion)+factor(Year_assessment_range),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize,
                           method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.multivariate_3, digits=5)

#Multivariate model omnibus test
GI.multivariate_3_results<- data.frame(df_1 = c(anova(GI.multivariate_3, btt=2:3)$m,anova(GI.multivariate_3, btt=4:8)$m,
                                                anova(GI.multivariate_3, btt=9:11)$m,anova(GI.multivariate_3, btt=12:15)$m,
                                                anova(GI.multivariate_3, btt=16:23)$m,anova(GI.multivariate_3, btt=24:26)$m,
                                                anova(GI.multivariate_3, btt=27:32)$m,anova(GI.multivariate_3, btt=33:37)$m),
                                       df_2 = c(anova(GI.multivariate_3, btt=2:3)$k,anova(GI.multivariate_3, btt=4:8)$k,
                                                anova(GI.multivariate_3, btt=9:11)$k,anova(GI.multivariate_3, btt=12:15)$k,
                                                anova(GI.multivariate_3, btt=16:23)$k,anova(GI.multivariate_3, btt=24:26)$k,
                                                anova(GI.multivariate_3, btt=27:32)$k,anova(GI.multivariate_3, btt=33:37)$k),
                                       QM = c(anova(GI.multivariate_3, btt=2:3)$QM,anova(GI.multivariate_3, btt=4:8)$QM,
                                              anova(GI.multivariate_3, btt=9:11)$QM,anova(GI.multivariate_3, btt=12:15)$QM,
                                              anova(GI.multivariate_3, btt=16:23)$QM,anova(GI.multivariate_3, btt=24:26)$QM,
                                              anova(GI.multivariate_3, btt=27:32)$QM,anova(GI.multivariate_3, btt=33:37)$QM),
                                       QMp = c(anova(GI.multivariate_3, btt=2:3)$QMp,anova(GI.multivariate_3, btt=4:8)$QMp,
                                               anova(GI.multivariate_3, btt=9:11)$QMp,anova(GI.multivariate_3, btt=12:15)$QMp,
                                               anova(GI.multivariate_3, btt=16:23)$QMp,anova(GI.multivariate_3, btt=24:26)$QMp,
                                               anova(GI.multivariate_3, btt=27:32)$QMp,anova(GI.multivariate_3, btt=33:37)$QMp),
                                       moderator = c("Crop_woodiness_T_recla","System_T","Fertiliser_chem_CT","Pesticide_CT",
                                                     "Soil_management_CT","Irrigation_CT",
                                                     "UN_subregion","Year_assessment_range"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(GI.multivariate_3_results, "C:/Users/Meta-analysis/Table_A10_GI_m_3.csv")

## Model 4
GI.multivariate_4<- rma.mv(Financial_mean, Financial_var, mods = ~factor(Crop_woodiness_T_recla)+
                             factor(System_T)+factor(Fertiliser_chem_CT)+
                             factor(Pesticide_CT)+factor(Soil_management_CT)+
                             factor(Irrigation_CT)+factor(UN_Development_status),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize,
                           method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.multivariate_4, digits=5)

#Multivariate model omnibus test
GI.multivariate_4_results<- data.frame(df_1 = c(anova(GI.multivariate_4, btt=2:3)$m,anova(GI.multivariate_4, btt=4:8)$m,
                                                anova(GI.multivariate_4, btt=9:11)$m,anova(GI.multivariate_4, btt=12:15)$m,
                                                anova(GI.multivariate_4, btt=16:23)$m,anova(GI.multivariate_4, btt=24:26)$m,
                                                anova(GI.multivariate_4, btt=27:27)$m),
                                       df_2 = c(anova(GI.multivariate_4, btt=2:3)$k,anova(GI.multivariate_4, btt=4:8)$k,
                                                anova(GI.multivariate_4, btt=9:11)$k,anova(GI.multivariate_4, btt=12:15)$k,
                                                anova(GI.multivariate_4, btt=16:23)$k,anova(GI.multivariate_4, btt=24:26)$k,
                                                anova(GI.multivariate_4, btt=27:27)$k),
                                       QM = c(anova(GI.multivariate_4, btt=2:3)$QM,anova(GI.multivariate_4, btt=4:8)$QM,
                                              anova(GI.multivariate_4, btt=9:11)$QM,anova(GI.multivariate_4, btt=12:15)$QM,
                                              anova(GI.multivariate_4, btt=16:23)$QM,anova(GI.multivariate_4, btt=24:26)$QM,
                                              anova(GI.multivariate_4, btt=27:27)$QM),
                                       QMp = c(anova(GI.multivariate_4, btt=2:3)$QMp,anova(GI.multivariate_4, btt=4:8)$QMp,
                                               anova(GI.multivariate_4, btt=9:11)$QMp,anova(GI.multivariate_4, btt=12:15)$QMp,
                                               anova(GI.multivariate_4, btt=16:23)$QMp,anova(GI.multivariate_4, btt=24:26)$QMp,
                                               anova(GI.multivariate_4, btt=27:27)$QMp),
                                       moderator = c("Crop_woodiness_T_recla","System_T","Fertiliser_chem_CT","Pesticide_CT",
                                                     "Soil_management_CT","Irrigation_CT","UN_Development_status"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(GI.multivariate_4_results, "C:/Users/Meta-analysis/Table_A10_GI_m_4.csv")

#---- Total cost
## Model 1
TC.multivariate_1<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_duration_T)+factor(Fertiliser_chem_CT)+
                             factor(Pesticide_CT)+factor(Irrigation_CT),
                         tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",
                         subset=(Financial_outcome=="Total cost"))
summary(TC.multivariate_1, digits=5)

#Multivariate model omnibus test
TC.multivariate_1_results<- data.frame(df_1 = c(anova(TC.multivariate_1, btt=2:10)$m,anova(TC.multivariate_1, btt=11:12)$m,
                                                anova(TC.multivariate_1, btt=13:16)$m,anova(TC.multivariate_1, btt=17:20)$m,
                                                anova(TC.multivariate_1, btt=21:23)$m),
                                       df_2 = c(anova(TC.multivariate_1, btt=2:10)$k,anova(TC.multivariate_1, btt=11:12)$k,
                                                anova(TC.multivariate_1, btt=13:16)$k,anova(TC.multivariate_1, btt=17:20)$k,
                                                anova(TC.multivariate_1, btt=21:23)$k),
                                       QM = c(anova(TC.multivariate_1, btt=2:10)$QM,anova(TC.multivariate_1, btt=11:12)$QM,
                                              anova(TC.multivariate_1, btt=13:16)$QM,anova(TC.multivariate_1, btt=17:20)$QM,
                                              anova(TC.multivariate_1, btt=21:23)$QM),
                                       QMp = c(anova(TC.multivariate_1, btt=2:10)$QMp,anova(TC.multivariate_1, btt=11:12)$QMp,
                                               anova(TC.multivariate_1, btt=13:16)$QMp,anova(TC.multivariate_1, btt=17:20)$QMp,
                                               anova(TC.multivariate_1, btt=21:23)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_duration","Fertiliser_chem_CT","Pesticide_CT","Irrigation_CT"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(TC.multivariate_1_results, "C:/Users/Meta-analysis/Table_A10_TC_m_1.csv")

## Model 2
TC.multivariate_2<- rma.mv(Financial_mean, Financial_var, mods = ~factor(Crop_woodiness_T)+
                             factor(System_T)+factor(Fertiliser_chem_CT)+
                             factor(Pesticide_CT)+factor(Irrigation_CT),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",
                           subset=(Financial_outcome=="Total cost"))
summary(TC.multivariate_2, digits=5)

#Multivariate model omnibus test
TC.multivariate_2_results<- data.frame(df_1 = c(anova(TC.multivariate_2, btt=2:7)$m,anova(TC.multivariate_2, btt=8:12)$m,
                                                anova(TC.multivariate_2, btt=13:16)$m,anova(TC.multivariate_2, btt=17:20)$m,
                                                anova(TC.multivariate_2, btt=21:23)$m),
                                       df_2 = c(anova(TC.multivariate_2, btt=2:7)$k,anova(TC.multivariate_2, btt=8:12)$k,
                                                anova(TC.multivariate_2, btt=13:16)$k,anova(TC.multivariate_2, btt=17:20)$k,
                                                anova(TC.multivariate_2, btt=21:23)$k),
                                       QM = c(anova(TC.multivariate_2, btt=2:7)$QM,anova(TC.multivariate_2, btt=8:12)$QM,
                                              anova(TC.multivariate_2, btt=13:16)$QM,anova(TC.multivariate_2, btt=17:20)$QM,
                                              anova(TC.multivariate_2, btt=21:23)$QM),
                                       QMp = c(anova(TC.multivariate_2, btt=2:7)$QMp,anova(TC.multivariate_2, btt=8:12)$QMp,
                                               anova(TC.multivariate_2, btt=13:16)$QMp,anova(TC.multivariate_2, btt=17:20)$QMp,
                                               anova(TC.multivariate_2, btt=21:23)$QMp),
                                       moderator = c("Crop_woodiness_T","System_T","Fertiliser_chem_CT","Pesticide_CT","Irrigation_CT"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(TC.multivariate_2_results, "C:/Users/Meta-analysis/Table_A10_TC_m_2.csv")

#---- Gross margin
## Model 1
GM.multivariate_1<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_duration_T)+factor(System_T)+
                             factor(Fertiliser_chem_CT)+factor(UN_subregion),
                         tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",
                         subset=(Financial_outcome=="Gross margin"))
summary(GM.multivariate_1, digits=5)

#Multivariate model omnibus test
GM.multivariate_1_results<- data.frame(df_1 = c(anova(GM.multivariate_1, btt=2:8)$m,anova(GM.multivariate_1, btt=9:9)$m,
                                                anova(GM.multivariate_1, btt=10:12)$m,anova(GM.multivariate_1, btt=13:15)$m,
                                                anova(GM.multivariate_1, btt=16:20)$m),
                                       df_2 = c(anova(GM.multivariate_1, btt=2:8)$k,anova(GM.multivariate_1, btt=9:9)$k,
                                                anova(GM.multivariate_1, btt=10:12)$k,anova(GM.multivariate_1, btt=13:15)$k,
                                                anova(GM.multivariate_1, btt=16:20)$k),
                                       QM = c(anova(GM.multivariate_1, btt=2:8)$QM,anova(GM.multivariate_1, btt=9:9)$QM,
                                              anova(GM.multivariate_1, btt=10:12)$QM,anova(GM.multivariate_1, btt=13:15)$QM,
                                              anova(GM.multivariate_1, btt=16:20)$QM),
                                       QMp = c(anova(GM.multivariate_1, btt=2:8)$QMp,anova(GM.multivariate_1, btt=9:9)$QMp,
                                               anova(GM.multivariate_1, btt=10:12)$QMp,anova(GM.multivariate_1, btt=13:15)$QMp,
                                               anova(GM.multivariate_1, btt=16:20)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_duration_T","System_T","Fertiliser_chem_CT","UN_subregion"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))
write.csv(GM.multivariate_1_results, "C:/Users/Meta-analysis/Table_A10_GM_m_1.csv")

## Model 2
GM.multivariate_2<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_duration_T)+factor(Fertiliser_chem_CT)+
                             factor(UN_Development_status),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.multivariate_2, digits=5)

#Multivariate model omnibus test
GM.multivariate_2_results<- data.frame(df_1 = c(anova(GM.multivariate_2, btt=2:8)$m,anova(GM.multivariate_2, btt=9:9)$m,
                                                anova(GM.multivariate_2, btt=10:12)$m,anova(GM.multivariate_2, btt=13:13)$m),
                                       df_2 = c(anova(GM.multivariate_2, btt=2:8)$k,anova(GM.multivariate_2, btt=9:9)$k,
                                                anova(GM.multivariate_2, btt=10:12)$k,anova(GM.multivariate_2, btt=7:7)$k),
                                       QM = c(anova(GM.multivariate_2, btt=2:2)$QM,anova(GM.multivariate_2, btt=13:13)$QM,
                                              anova(GM.multivariate_2, btt=4:6)$QM,anova(GM.multivariate_2, btt=7:7)$QM),
                                       QMp = c(anova(GM.multivariate_2, btt=2:8)$QMp,anova(GM.multivariate_2, btt=9:9)$QMp,
                                               anova(GM.multivariate_2, btt=10:12)$QMp,anova(GM.multivariate_2, btt=13:13)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_duration_T","Fertiliser_chem_CT","UN_Development_status"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(GM.multivariate_2_results, "C:/Users/Meta-analysis/Table_A10_GM_m_2.csv")

## Model 3
GM.multivariate_3<- rma.mv(Financial_mean, Financial_var, mods = ~factor(Crop_woodiness_T)+
                             factor(System_T)+factor(Fertiliser_chem_CT)+
                             factor(UN_subregion),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.multivariate_3, digits=5)

#Multivariate model omnibus test
GM.multivariate_3_results<- data.frame(df_1 = c(anova(GM.multivariate_3, btt=2:3)$m,anova(GM.multivariate_3, btt=9:6)$m,
                                                anova(GM.multivariate_3, btt=10:9)$m,anova(GM.multivariate_3, btt=13:15)$m),
                                       df_2 = c(anova(GM.multivariate_3, btt=2:8)$k,anova(GM.multivariate_3, btt=9:9)$k,
                                                anova(GM.multivariate_3, btt=10:12)$k,anova(GM.multivariate_3, btt=7:7)$k),
                                       QM = c(anova(GM.multivariate_3, btt=2:2)$QM,anova(GM.multivariate_3, btt=13:13)$QM,
                                              anova(GM.multivariate_3, btt=4:6)$QM,anova(GM.multivariate_3, btt=7:7)$QM),
                                       QMp = c(anova(GM.multivariate_3, btt=2:8)$QMp,anova(GM.multivariate_3, btt=9:9)$QMp,
                                               anova(GM.multivariate_3, btt=10:12)$QMp,anova(GM.multivariate_3, btt=13:13)$QMp),
                                       moderator = c("Crop_woodiness_T","System_T","Fertiliser_chem_CT","UN_subregion"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))
write.csv(GM.multivariate_3_results, "C:/Users/Meta-analysis/Table_A10_GM_m_3.csv")

## Model 4
GM.multivariate_4<- rma.mv(Financial_mean, Financial_var, mods = ~factor(Crop_woodiness_T)+
                             factor(Fertiliser_chem_CT)+factor(UN_Development_status),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.multivariate_4, digits=5)

#Multivariate model omnibus test
GM.multivariate_4_results<- data.frame(df_1 = c(anova(GM.multivariate_4, btt=2:3)$m,anova(GM.multivariate_4, btt=4:6)$m,
                                                anova(GM.multivariate_4, btt=7:7)$m),
                                       df_2 = c(anova(GM.multivariate_4, btt=2:3)$k,anova(GM.multivariate_4, btt=4:6)$k,
                                                anova(GM.multivariate_4, btt=7:7)$k),
                                       QM = c(anova(GM.multivariate_4, btt=2:3)$QM,anova(GM.multivariate_4, btt=4:6)$QM,
                                              anova(GM.multivariate_4, btt=7:7)$QM),
                                       QMp = c(anova(GM.multivariate_4, btt=2:3)$QMp,anova(GM.multivariate_4, btt=4:6)$QMp,
                                               anova(GM.multivariate_4, btt=7:7)$QMp),
                                       moderator = c("Crop_woodiness_T","Fertiliser_chem_CT","UN_Development_status"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(GM.multivariate_4_results, "C:/Users/Andrea/Documents/Bioversity/Cost_benefits_analysis/Meta-analysis/Results_2022.06.21/Table_A10_GM_m_4.csv")

#---- Net income
## Model 1
NI.multivariate_1<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_duration_T)+factor(System_T)+
                             factor(Fertiliser_chem_CT)+factor(Pesticide_CT)+
                             factor(UN_Development_status),
                         tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.multivariate_1, digits=5)

#Multivariate model omnibus test
NI.multivariate_1_results<- data.frame(df_1 = c(anova(NI.multivariate_1, btt=2:13)$m,anova(NI.multivariate_1, btt=14:15)$m,
                                                anova(NI.multivariate_1, btt=16:20)$m,anova(NI.multivariate_1, btt=21:24)$m,
                                                anova(NI.multivariate_1, btt=25:28)$m,anova(NI.multivariate_1, btt=29:29)$m),
                                       df_2 = c(anova(NI.multivariate_1, btt=2:13)$k,anova(NI.multivariate_1, btt=14:15)$k,
                                                anova(NI.multivariate_1, btt=16:20)$k,anova(NI.multivariate_1, btt=21:24)$k,
                                                anova(NI.multivariate_1, btt=25:28)$k,anova(NI.multivariate_1, btt=29:29)$k),
                                       QM = c(anova(NI.multivariate_1, btt=2:13)$QM,anova(NI.multivariate_1, btt=14:15)$QM,
                                              anova(NI.multivariate_1, btt=16:20)$QM,anova(NI.multivariate_1, btt=21:24)$QM,
                                              anova(NI.multivariate_1, btt=25:28)$QM,anova(NI.multivariate_1, btt=29:29)$QM),
                                       QMp = c(anova(NI.multivariate_1, btt=2:13)$QMp,anova(NI.multivariate_1, btt=14:15)$QMp,
                                               anova(NI.multivariate_1, btt=16:20)$QMp,anova(NI.multivariate_1, btt=21:24)$QMp,
                                               anova(NI.multivariate_1, btt=25:28)$QMp,anova(NI.multivariate_1, btt=29:29)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_duration_T","System_T","Fertiliser_chem_CT",
                                                     "Pesticide_CT","UN_Development_status"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(NI.multivariate_1_results, "C:/Users/Meta-analysis/Table_A10_NI_m_1.csv")

## Model 2
NI.multivariate_2<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_duration_T)+factor(System_T)+
                             factor(Fertiliser_chem_CT)+factor(Pesticide_CT)+
                             factor(Year_assessment_range),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",
                           subset=(Financial_outcome=="Net income"))
summary(NI.multivariate_2, digits=5)

#Multivariate model omnibus test
NI.multivariate_2_results<- data.frame(df_1 = c(anova(NI.multivariate_2, btt=2:13)$m,anova(NI.multivariate_2, btt=14:15)$m,
                                                anova(NI.multivariate_2, btt=16:20)$m,anova(NI.multivariate_2, btt=21:24)$m,
                                                anova(NI.multivariate_2, btt=25:28)$m,anova(NI.multivariate_2, btt=29:33)$m),
                                       df_2 = c(anova(NI.multivariate_2, btt=2:11)$k,anova(NI.multivariate_2, btt=14:15)$k,
                                                anova(NI.multivariate_2, btt=16:20)$k,anova(NI.multivariate_2, btt=21:24)$k,
                                                anova(NI.multivariate_2, btt=25:28)$k,anova(NI.multivariate_2, btt=29:33)$k),
                                       QM = c(anova(NI.multivariate_2, btt=2:13)$QM,anova(NI.multivariate_2, btt=14:15)$QM,
                                              anova(NI.multivariate_2, btt=16:20)$QM,anova(NI.multivariate_2, btt=21:24)$QM,
                                              anova(NI.multivariate_2, btt=25:28)$QM,anova(NI.multivariate_2, btt=29:33)$QM),
                                       QMp = c(anova(NI.multivariate_2, btt=2:13)$QMp,anova(NI.multivariate_2, btt=14:15)$QMp,
                                               anova(NI.multivariate_2, btt=16:20)$QMp,anova(NI.multivariate_2, btt=21:24)$QMp,
                                               anova(NI.multivariate_2, btt=25:28)$QMp,anova(NI.multivariate_2, btt=29:33)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_duration_T","System_T","Fertiliser_chem_CT",
                                                     "Pesticide_CT","Year_assessment_range"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

View(NI.multivariate_2_results)
write.csv(NI.multivariate_2_results, "C:/Users/Meta-analysis/Table_A10_NI_m_2.csv")

## Model 3
NI.multivariate_3<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_woodiness_T_recla)+factor(System_T)+
                             factor(Fertiliser_chem_CT)+factor(Pesticide_CT)+
                             factor(UN_Development_status),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.multivariate_3, digits=5)

#Multivariate model omnibus test
NI.multivariate_3_results<- data.frame(df_1 = c(anova(NI.multivariate_3, btt=2:13)$m,anova(NI.multivariate_3, btt=14:15)$m,
                                                anova(NI.multivariate_3, btt=16:20)$m,anova(NI.multivariate_3, btt=21:24)$m,
                                                anova(NI.multivariate_3, btt=25:28)$m,anova(NI.multivariate_3, btt=29:29)$m),
                                       df_2 = c(anova(NI.multivariate_3, btt=2:13)$k,anova(NI.multivariate_3, btt=14:15)$k,
                                                anova(NI.multivariate_3, btt=16:20)$k,anova(NI.multivariate_3, btt=21:24)$k,
                                                anova(NI.multivariate_3, btt=25:28)$k,anova(NI.multivariate_3, btt=29:29)$k),
                                       QM = c(anova(NI.multivariate_3, btt=2:13)$QM,anova(NI.multivariate_3, btt=14:15)$QM,
                                              anova(NI.multivariate_3, btt=16:20)$QM,anova(NI.multivariate_3, btt=21:24)$QM,
                                              anova(NI.multivariate_3, btt=25:28)$QM,anova(NI.multivariate_3, btt=29:29)$QM),
                                       QMp = c(anova(NI.multivariate_3, btt=2:13)$QMp,anova(NI.multivariate_3, btt=14:15)$QMp,
                                               anova(NI.multivariate_3, btt=16:20)$QMp,anova(NI.multivariate_3, btt=21:24)$QMp,
                                               anova(NI.multivariate_3, btt=25:28)$QMp,anova(NI.multivariate_3, btt=29:29)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_woodiness_T_recla","System_T","Fertiliser_chem_CT",
                                                     "Pesticide_CT","UN_Development_status"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(NI.multivariate_3_results, "C:/Users/Meta-analysis/Table_A10_NI_m_3.csv")

## Model 4
NI.multivariate_4<- rma.mv(Financial_mean, Financial_var, mods = ~factor(FAO_group_T_recla)+
                             factor(Crop_woodiness_T_recla)+factor(System_T)+
                             factor(Fertiliser_chem_CT)+factor(Pesticide_CT)+
                             factor(Year_assessment_range),
                           tdist=TRUE,random = list(~ 1 | ES_ID, ~ 1 | ID),data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.multivariate_4, digits=5)


#Multivariate model omnibus test
NI.multivariate_4_results<- data.frame(df_1 = c(anova(NI.multivariate_4, btt=2:13)$m,anova(NI.multivariate_4, btt=14:15)$m,
                                                anova(NI.multivariate_4, btt=16:20)$m,anova(NI.multivariate_4, btt=21:24)$m,
                                                anova(NI.multivariate_4, btt=25:28)$m,anova(NI.multivariate_4, btt=29:33)$m),
                                       df_2 = c(anova(NI.multivariate_4, btt=2:13)$k,anova(NI.multivariate_4, btt=14:15)$k,
                                                anova(NI.multivariate_4, btt=16:20)$k,anova(NI.multivariate_4, btt=21:24)$k,
                                                anova(NI.multivariate_4, btt=25:28)$k,anova(NI.multivariate_4, btt=29:33)$k),
                                       QM = c(anova(NI.multivariate_4, btt=2:13)$QM,anova(NI.multivariate_4, btt=14:15)$QM,
                                              anova(NI.multivariate_4, btt=16:20)$QM,anova(NI.multivariate_4, btt=21:24)$QM,
                                              anova(NI.multivariate_4, btt=25:28)$QM,anova(NI.multivariate_4, btt=29:33)$QM),
                                       QMp = c(anova(NI.multivariate_4, btt=2:13)$QMp,anova(NI.multivariate_4, btt=14:15)$QMp,
                                               anova(NI.multivariate_4, btt=16:20)$QMp,anova(NI.multivariate_4, btt=21:24)$QMp,
                                               anova(NI.multivariate_4, btt=25:28)$QMp,anova(NI.multivariate_4, btt=29:33)$QMp),
                                       moderator = c("FAO_group_T_recla","Crop_woodiness_T_recla","System_T","Fertiliser_chem_CT",
                                                     "Pesticide_CT","Year_assessment_range"))%>%
  mutate(QM =round(QM, digits = 3),
         QMp =round(QMp, digits = 4),
         QMp_1 = if_else(QMp < 0.001, 0, QMp),
         QMp_1 = as.character(QMp_1),
         QMp_1 = if_else(QMp_1=="0","< 0.001",QMp_1))%>%
  mutate(significance = if_else(QMp <=0.001,"***",
                                if_else(QMp>0.001&QMp<0.01,"**",
                                        if_else(QMp>0.01&QMp<=0.05,"*",
                                                if_else(QMp>0.05&QMp<=0.1,"","")))))%>%
  mutate(label=paste("F (",df_1,", ",df_2,") = ", QM,", p", QMp_1, significance, sep=""))

write.csv(NI.multivariate_4_results, "C:/Users/Meta-analysis/Table_A10_NI_m_4.csv")

#__________________________________________________________________________________________________________________________________________________________
################################  PUBLICATION BIAS ###################################################################3----------------------------------------------------#
##--------- Egger regression test
#B/C ratio
BCR.bias.egger <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Financial_precision, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                             tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="B/C ratio"))
summary(BCR.bias.egger, digits=3)

#Gross income
GI.bias.egger <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Financial_precision, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                         tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.bias.egger, digits=3)

#Total cost
TC.bias.egger <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Financial_precision, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.bias.egger, digits=3)

sort(unique(prueba$Financial_precision))

#Gross margin
GM.bias.egger <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Financial_precision, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.bias.egger, digits=3)

#Net income
NI.bias.egger <- rma.mv(y=Financial_mean, V=Financial_var, mods = ~ Financial_precision, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=effectsize, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.bias.egger, digits=3)

#__________________________________________________________________________________________________________________________________________________________
################################  SENSITIVITY ANALYSIS ###################################################################3----------------------------------------------------#
getMKLthreads()
##COOK´S DISTANCE, IDENTIFICATION OF EFFECT SIZE OUTLIERS
# B/C ratio, Gross income, Total cost
system.time(BCR.overall_cook <- cooks.distance(BCR.overall, reestimate=FALSE))
system.time(GI.overall_cook <- cooks.distance(GI.overall, reestimate=FALSE))
system.time(TC.overall_cook <- cooks.distance(TC.overall, reestimate=FALSE))
system.time(GM.overall_cook <- cooks.distance(GM.overall, reestimate=FALSE))
system.time(NI.overall_cook <- cooks.distance(NI.overall, reestimate=FALSE))
BCR.overall_cook
overall_cook_df<- as.data.frame(BCR.overall_cook)%>%
  rename("GM.overall_cook"="BCR.overall_cook")%>%
  rbind(as.data.frame(GM.overall_cook))%>%
  rename("NI.overall_cook"="GM.overall_cook")%>%
  rbind(as.data.frame(NI.overall_cook))%>%
  rename("GI.overall_cook"="NI.overall_cook")%>%
  rbind(as.data.frame(GI.overall_cook))%>%
  rename("TC.overall_cook"="GI.overall_cook")%>%
  rbind(as.data.frame(TC.overall_cook))%>%
  rename("overall_cook"="TC.overall_cook")%>%
  mutate(ES_ID=rownames(.))%>%
    mutate(Financial_outcome=c(rep("B/C ratio",90),rep("Gross margin",800),
                               rep("Net income",1265),rep("Gross income",534),
                               rep("Total cost",503)),
           ES_ID = as.numeric(ES_ID))%>%
  group_by(Financial_outcome)%>%
  mutate(n =length(Financial_outcome),
         cut_off_y = (4/(n-1)))%>%
  #mutate(cut_off_y = if_else(Financial_outcome =="Gross margin",0.01,cut_off_y))%>%
  mutate(overall_cook_outlier = if_else(overall_cook >cut_off_y ,overall_cook,200))%>%
  mutate(ES_ID_outlier = if_else(overall_cook > cut_off_y,ES_ID,5000))%>%
  mutate(Financial_group = sapply(as.character(Financial_outcome), switch, "B/C ratio" = 1, "Gross margin" = 2, "Net income" = 3, "Gross income" = 4, "Total cost"=5,
                           USE.NAMES = F))
View(overall_cook_df)

#Recomended cutoff for Cooks distance 4/(n-p) Bollen, K. A., & Jackman, R. W. Regression diagnostics: An expository treatment of
sensitivity_cook_outlier<- overall_cook_df%>%
  mutate(ES_ID=as.character(ES_ID))%>%
  left_join(effectsize, by="ES_ID")%>%
  filter(overall_cook_outlier==200)

View(sensitivity_cook_outlier)

## Cook´s distance sensitivity analysis
#B/C ratio
BCR.cook <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                                 tdist= TRUE, data=sensitivity_cook_outlier, method="REML", subset=(Financial_outcome.x=="B/C ratio"))
summary(BCR.cook, digits=3)

#Gross income
GI.cook <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=sensitivity_cook_outlier, method="REML",subset=(Financial_outcome.x=="Gross income"))
summary(GI.cook, digits=3)

#Total cost
TC.cook <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=sensitivity_cook_outlier, method="REML",subset=(Financial_outcome.x=="Total cost"))
summary(TC.cook, digits=3)

#Gross margin
GM.cook <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=sensitivity_cook_outlier, method="REML",subset=(Financial_outcome.x=="Gross margin"))
summary(GM.cook, digits=3)

#Net income
NI.cook <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                        tdist=TRUE, data=sensitivity_cook_outlier, method="REML",subset=(Financial_outcome.x=="Net income"))
summary(NI.cook, digits=3)

### IDENTIFICATION OF EFFECT SIZES WITH HIGH AND MEDIUM RISK OF BIAS
risk_identification<- effectsize%>%group_by(Financial_outcome,Validity_overall)%>%
  mutate(N_articles = n_distinct(ID))%>%ungroup()%>%
  group_by(Financial_outcome,Validity_overall)%>%
  mutate(N_articles = mean(N_articles))%>%
  group_by(Financial_outcome,Validity_overall,N_articles)%>%tally()%>%
  rename("N_effect_size" = "n")%>%group_by()%>%
  mutate(Validity_overall=if_else(Validity_overall==5.5,6,Validity_overall))%>%
  add_row(Financial_outcome="Gross margin", Validity_overall = 0, N_articles=0,N_effect_size=0)%>%
  add_row(Financial_outcome="Gross margin", Validity_overall = 1, N_articles=0,N_effect_size=0)%>%
  add_row(Financial_outcome="Gross margin", Validity_overall = 2, N_articles=0,N_effect_size=0)%>%
  complete(., Financial_outcome, Validity_overall)%>%
  mutate(N_articles= if_else(is.na(N_articles),0,N_articles),
         N_effect_size=as.numeric(N_effect_size),
         N_effect_size = if_else(is.na(N_effect_size),0,N_effect_size))%>%
  mutate(y_axis = paste(Validity_overall, " (",N_effect_size,", ",N_articles,")",sep=""))%>%
  left_join(effectsize, by=c("Financial_outcome", "Validity_overall"))%>%
  mutate(Financial_group = sapply(as.character(Financial_outcome), switch, "B/C ratio" = 1, "Gross margin" = 2, "Net income" = 3, "Gross income" = 4, "Total cost"=5,
                                  USE.NAMES = F))
View(risk_identification)




## SENSITIVITY ANALYSIS EXCLUDING EFFECT SIZES WITH HIGH RISK OF BIAS
sensitivity_risk_bias <- effectsize%>%
  filter(Validity_overall != 3 & Validity_overall !=4)
  
sort(unique(sensitivity_risk_bias$Validity_overall))
sort(unique(sensitivity_risk_bias$Financial_outcome))

#B/C ratio
BCR.risk <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                   tdist= TRUE, data=sensitivity_risk_bias, method="REML", subset=(Financial_outcome=="B/C ratio"))
summary(BCR.risk, digits=3)

#Gross income
GI.risk <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist=TRUE, data=sensitivity_risk_bias, method="REML",subset=(Financial_outcome=="Gross income"))
summary(GI.risk, digits=3)

#Total cost
TC.risk <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist=TRUE, data=sensitivity_risk_bias, method="REML",subset=(Financial_outcome=="Total cost"))
summary(TC.risk, digits=3)

#Gross margin
GM.risk <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist=TRUE, data=sensitivity_risk_bias, method="REML",subset=(Financial_outcome=="Gross margin"))
summary(GM.risk, digits=3)

#Net income
NI.risk <- rma.mv(y=Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist=TRUE, data=sensitivity_risk_bias, method="REML",subset=(Financial_outcome=="Net income"))
summary(NI.risk, digits=3)

###### EXCLUDING THE ARTICLES WITH SIMPLIFIED SYSTEMS DIFFERENT THAN MONOCULTURE
sensitivity_monoculture<- effectsize%>%
  filter(System_raw_C== "Monoculture")

sort(unique(effectsize$System_raw_C))
BCR.mono <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                      tdist= TRUE, data=sensitivity_monoculture, method="REML", subset=(Financial_outcome=="B/C ratio"))
summary(BCR.mono, digits=3)
summary(BCR.overall, digits=3)

GI.mono <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                   tdist= TRUE, data=sensitivity_monoculture, method="REML", subset=(Financial_outcome=="Gross income"))
summary(GI.overall, digits=3)

TC.mono <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist= TRUE, data=sensitivity_monoculture, method="REML", subset=(Financial_outcome=="Total cost"))
summary(TC.mono, digits=3)
summary(TC.overall, digits=3)


NI.mono <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist= TRUE, data=sensitivity_monoculture, method="REML", subset=(Financial_outcome=="Net income"))
summary(NI.mono, digits=3)

GM.mono <- rma.mv(y= Financial_mean, V=Financial_var, random = list(~ 1 | ES_ID, ~ 1 | ID), 
                  tdist= TRUE, data=sensitivity_monoculture, method="REML", subset=(Financial_outcome=="Gross margin"))
summary(GM.mono, digits=3)