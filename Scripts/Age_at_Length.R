# Age samples
# read in sample metadata, using TL, estimate the age of samples using selected model, assign a cohort 


meta_data <- readr::read_csv("Data/Raw/Copy of UPDTAED LAT LONGS UQ genetics PhD181018 Paul Butcher NSW DPI.csv") %>%
  #select(MBB_Code, FL, TL) %>%
  mutate(date_column = as.Date(`Date tagged`, format="%d/%m/%y"),
         modified_year = year(date_column))

# 1. relationship between TL and FL 
lm1 <- lm(meta_data$TL ~ meta_data$FL)

Tl_intercept <- lm1$coefficients[1]
Tl_slope <- lm1$coefficients[2]

TL <- sapply(meta_data$FL, function(x) {Tl_intercept + x * Tl_slope }) # eq (1) manuscript 

# model age 
# The von Bertalanffy growth function is defined as:
#Length ~ Linf - ((Linf - L0) * exp(-k * Age)) 
#VBGF2 <-  
  #Length ~ Linf*(1-exp(-K*(Age-t0)))
  #solve for age 
# VBGF3 <-   
  ### E[L|t] = L∞ − (L∞ − L0) e−Kt  # Linf - (Linf - L0) * exp(-K * t_0)
  ### the model proposed by von Bertalanffy, relatively (as compared to the typical VBGM) rarely used in the literature. However, Cailliet et al. (2006) recommended its use with chondrychthians.
  
vbgf3 <- function(L, L_infinity, K_value, t_0) {
  result <- t_0 - (1 / K_value) * log(1 - (L / L_infinity))
  return(result)
}


# vbgf3(t_value, L_infinity_value, K_value, t_0_value)

# Estimate age directly using values from O'Connor 2011
#solve for age

# Estimating age at length
# info:
# Best fitting parmaters for Males
Linf <- 798.94# cm TL
k <- 0.047 # 
L0 <- 140 #cm 
T0 <- -3.8 #years 
Males <- list(L_infinity_value = Linf, K_value = k, L0 = L0, t_0_value = T0)

# Best fitting paramters for Females 
Linf<- 719.02# cm TL
k<- 0.056 # 
L0<- 140 #cm 
T0 <- -3.8 #years 
Females <- list(L_infinity_value = Linf, K_value = k, L0 = L0, t_0_value = T0)

Linf<- 746.66# cm TL
k<- 0.053# 
L0<- 140 #cm 
T0 <- -3.8 #years 
Combined <- list(L_infinity_value = Linf, K_value = k, L0 = L0, t_0_value = T0)

# female 
idx_f<- which(meta_data$Sex == "F")
idx_m<- which(meta_data$Sex == "M")

male_t = sapply(TL[idx_m], function(x) {vbgf3(x, Males[["L_infinity_value"]], Males[["K_value"]], Males[["t_0_value"]])})
female_t = sapply(TL[idx_f], function(x) {vbgf3(x, Females[["L_infinity_value"]], Females[["K_value"]], Females[["t_0_value"]])})

# Cohort is equal to year sampled minus age
male_age = trunc(meta_data$modified_year[idx_m] -  male_t)
female_age = trunc(meta_data$modified_year[idx_f] -  female_t) # including some samples, ie. MBB1348 depends on how you round the age or the year
#male_age = round(meta_data$modified_year[idx_m] -  male_t)
#female_age = round(meta_data$modified_year[idx_f] -  female_t) # including some samples, ie. MBB1348 depends on how you round the age or the year
output = tibble::tibble(STRATA = c(male_age, female_age), MBB_CODE = c(meta_data$MBB_Code[idx_m], meta_data$MBB_Code[idx_f]))
readr::write_tsv(output, "Data/Processed/New_Cohort_Assignment.tsv")

                  
# check with tropfishR
#remotes::install_github("tokami/TropFishR")
t_male_trop = TropFishR::VBGF(L = TL[idx_m], param = list(Linf= Males[[1]], K= Males[[2]], t0 = Males[[4]], L0 = Males[[3]])) # thats age?
f_female_top = TropFishR::VBGF(L = TL[idx_f], param = list(Linf= Females[[1]], K= Females[[2]], t0 = Females[[4]], L0 = Females[[3]])) # thats age?
male_age_trop = round(meta_data$modified_year[idx_m] -  t_male_trop)
female_age_trop = round(meta_data$modified_year[idx_f] -  f_female_top)
output_trop = tibble::tibble(STRATA = c(male_age_trop, female_age_trop), MBB_CODE = c(meta_data$MBB_Code[idx_m], meta_data$MBB_Code[idx_f]))

# the tropfish results are slightly different????



