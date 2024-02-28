library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)


# This file will pull the Weighted Mean of r^2-drift
# To get this value, you need to run Ne estimator with the check-box marked for Output Burrows coefficient, then pull the line that says "Weighted Mean of r^2-drift" using this script

file = "newBur.txt"
 
bur_stripper<- function(x){
  neinfo<- tibble(lines = readLines(x))
  p = grep('\\POPULATION+', neinfo$lines, value = TRUE)
  r = grep('\\# Weighted Mean of r\\^2\\-drift+', neinfo$lines, value = TRUE)
  rnew = r[-grep('\\Revised', r)]
  return(tibble::tibble("Pop_name" = str_match(p, "POPULATION\\s*(.*?)\\s*\\(")[,2],
          "Pop_samplesize" = str_match(p, "\\= \\s*(.*?)\\s*\\)")[,2], 
          "Rstar" = str_match(rnew, "\\= \\s*(.*?)\\s*\\(")[,2]))
}

output = bur_stripper(file)
readr::write_tsv(output, file = "rstar_extract.tsv")
