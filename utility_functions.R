#some functions for later

Ne_to_r <- function(Ne, S){
  r <- (-69*S^2 + sqrt(10000 * (S^4) * (Ne^2) + 4761 * (S^4) - 248400*(S^3)*Ne) + 1800*S*(Ne^2))/( 1800*(S^2)*(Ne^2))  
  return(r)
}

Ne_to_r_star <- function(Ne, S){
  r <- Ne_to_r(Ne, S)
  r <- r - 1/S
  return(r)
}


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

# This file will pull the Weighted Mean of r^2-drift
# To get this value, you need to run Ne estimator with the check-box marked for Output Burrows coefficient, then pull the line that says "Weighted Mean of r^2-drift" using this script

bur_stripper <- function(x){
  neinfo<- tibble(lines = readLines(x))
  p = grep('\\POPULATION+', neinfo$lines, value = TRUE)
  r = grep('\\# Weighted Mean of r\\^2\\-drift+', neinfo$lines, value = TRUE)
  rnew = r[-grep('\\Revised', r)]
  return(tibble::tibble("Pop_name" = str_match(p, "POPULATION\\s*(.*?)\\s*\\(")[,2],
                        "Pop_samplesize" = str_match(p, "\\= \\s*(.*?)\\s*\\)")[,2], 
                        "Rstar" = str_match(rnew, "\\= \\s*(.*?)\\s*\\(")[,2]))
}


# file = "newBur.txt"
# output = bur_stripper(file)
# readr::write_tsv(output, file = "rstar_extract.tsv")