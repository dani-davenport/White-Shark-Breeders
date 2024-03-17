---
title: "Age_At_Length"
output: html_document
author: "Danielle Davenport"
---
# setup
  ## brew install jags
  ## install.packages("rjags")
  ## library(rjags)

########## 
# The von Bertalanffy growth function is defined as:
#Length ~ Linf - ((Linf - L0) * exp(-k * Age)) 
VBGF2 <-  
  #Length ~ Linf*(1-exp(-K*(Age-t0)))
  #solve for age 
VBGF3 <-   
  ### E[L|t] = L∞ − (L∞ − L0) e−Kt
  ### the model proposed by von Bertalanffy, relatively (as compared to the typical VBGM) rarely used in the literature. However, Cailliet et al. (2006) recommended its use with chondrychthians.
  
## solve for age
#t(L) = (t0 - (1/k)ln(1-(L/Linf))) #implemented in tropfish package 
  
# Estimating age at length
  # info:
  # Best fitting parmaters for Males
Linf<- 798.94# cm TL
k<- 0.047 # 
L0<- 140 #cm 
T0 <- -3.8 #years 
Males<- list(Linf, k, L0, T0)
  # Best fitting paramters for Females 
Linf<- 719.02# cm TL
k<- 0.056 # 
L0<- 140 #cm 
T0 <- -3.8 #years 
Females<- list(Linf, k, L0, T0)
  # Best fitting paramter for Combined sexes
Linf<- 746.66# cm TL
k<- 0.053# 
L0<- 140 #cm 
T0 <- -3.8 #years 
Combined<- list(Linf, k, L0, T0)
# Estimate age directly using values from O'Connor
library("TropFishR") 
# calculating ages
Sample_Info = read.csv("Sample_TL_Raw.csv", header = T)
L <- Sample_Info$TL
Sex<- Sample_Info$Sex
DOC_Year<- Sample_Info$DOC_Year

# combined
t_combined <- TropFishR::VBGF(L = L, param = list(Linf= Combined[[1]], K= Combined[[2]], t0 = Combined[[4]], L0 = Combined[[3]]))
plot(t_combined, L)

# female 
idx_f<- which(Sex == "F")
t_female <- TropFishR::VBGF(L = L[idx_f], param = list(Linf= Females[[1]], K= Females[[2]], t0 =Females[[4]], L0 = Females[[3]]))
plot(t_female, L[idx_f])
# male
idx_m<- which(Sex == "M")
t_male <- TropFishR::VBGF(L = L[idx_m], param = list(Linf= Males[[1]], K= Males[[2]], t0 = Males[[4]], L0 = Males[[3]]))
plot(t_male, L[idx_m])
# other
idx_none<- which(is.na(Sex))
t_none <- TropFishR::VBGF(L = L[idx_none], param = list(Linf= Combined[[1]], K= Combined[[2]], t0 = Combined[[4]], L0 = Combined[[3]]))
plot(t_none, L[idx_none])
# write out
sample_specific<- as.data.frame(rbind(cbind(idx_f, t_female), cbind(idx_m, t_male), cbind(idx_none, t_none)))
colnames(sample_specific)<- c("index", "sex_model_age")
out<- cbind(Sample_Info, sample_specific[order(sample_specific$index),], t_combined)

out_2<- cbind(out, sex_model_cohort = mapply(function(x,y) trunc(x-y), out$DOC_Year, out$sex_model_age))
out_2<- cbind(out_2, combined_cohort = mapply(function(x,y) trunc(x-y), out$DOC_Year, out$t_combined))
write.csv(out_2, file = "Cohort_Info.csv", quote = F,row.names = F, sep = ",")
# To Do
#
#
## Bayseain and also NLS modelling with O'Connor Data - when it arrives to check it and also to get CI around age
# setwd("")
# Read data
setwd("~/Google Drive/UQ_PhD/Lab_Book/ECWS_NSW_Fisheries/Data/Analysis/Age_Samples")
oconnor<- read.csv("../../Sample_TL_Age.csv", header = T)
# Plot relationship between length and age
dataset <- data.frame(oconnor$Age, oconnor$Length)
colnames(dataset) <- c("Age", "Total Length")

#draw graph
ggplot(dataset, aes(x = Age, y = TL)) +
  geom_point()
# Estimation of starting parameters, informative priors based on previous studies 
# Priors based on X?
########## 

## O'Connor 
## Best fitting paramters combined sexes
Linf<- 746.66 # cm TL
k<- 0.053 # is the so-called Brody growth rate coefficient, (units are yr−1) that is time^−1
L0<- 140 #cm as in Francis 1996 and Malcom et al., 2001 - length at birth ## in three parameter model this is like the T0 parameter which is the theoretical age at length 0
t0 <- -3.8 #years 

# make non linear model from known paramters
# NLS
sapply(as.numeric(oconnor$TL), VBGF2)
n2<- nls(VBGF2, start = list(Linf = Linf, k = k , t0 = t0 ))
n3<- nls(VBGF3, start = list(Linf = Linf, k = k , L0 = L0 ))

predict(n2, oconnor$TL,t=3)


## TO DO
# Bayesian Modeling of Age at Length 
# Hillary 2018
## Cohort is equal to year sampled minus age

# rjags package
c <-  y - age 



## prepping sara file 

sara<- read.csv("../../../Data/NSWDPI_Sent_Sara.csv", header = T)
colnames(sara)<- c("PlateID","Row","Column","Organism","Species","TARGET_ID", "Comments")
strata<- read.table("../Full_Strata.tsv", header = T)
cohort<- read.csv("../Age_Samples/Cohort_Info.csv")[,1:10]
colnames(cohort)[1]<- c("TARGET_ID")
sara2<- dplyr::inner_join(sara, cohort, by = "TARGET_ID")
write.csv(sara2, file = "../../../Data/NSWDPI_Sara_Full_Info.csv",quote = F, row.names = F)

