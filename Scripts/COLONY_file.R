all <- genomic_converter("//NASv2/Transfers2/new.gen")

for(i in 5:9){
  year <- i + 2005
  all$tidy.data[all$tidy.data$STRATA == i, ] %>% write_colony(filename = paste(paste0("//NASv2/Transfers2/COLONY_input_", year,"_v2.txt")))
}



devtools::install_github("thierrygosselin/radiator")

library(radiator)
library(tidyr)

all <- genomic_converter("//NASv2/Transfers2/new.gen")

for(i in 5:9){
  year <- i + 2005
  all$tidy.data[all$tidy.data$STRATA == i, ] %>% write_colony(filename = paste(paste0("//NASv2/Transfers2/COLONY_input_", year,"_v2.txt")))
}












write_COLONY <- function(agenepop, current_path = ".") {
  
  cohort_pops <- adegenet::seppop(agenepop)
  marker.num = length(agenepop@all.names)
  
  n_pops <- length(cohort_pops)
  # write a colony file - (see website for specs) https://www.zsl.org/science/software/colony
  for (i in seq_len(n_pops)){
    
    data <- cohort_pops[[i]]
    obj <- data
    year <- str_split(names(cohort_pops)[i], "_")[[1]][1]
    filename = file.path(current_path, paste0("COLONY_input_", year, ".txt"))
    # line 1 - dataset name
    line_1 = paste0(dQuote(paste0("Colony Input Pop: ", pop)),"                 ! Dataset name\n")
    readr::write_file(line_1, filename, append = FALSE)
    # line 2 - output file name 
    line_2 = paste0(dQuote(paste0("Output_", pop)),"                 ! Output file name\n")
    readr::write_file(line_2, filename, append = T)
    # line 3 - offspring number 
    line_3 <- paste(nrow(data), "                                  ! Number of offspring in the sample\n", sep = "")
    readr::write_file(line_3, filename, append = T)
    # line 4 - number of markers
    line_4 <- paste(marker.num, "                                 ! Number of loci\n", sep = "")
    readr::write_file(line_4, filename, append = T)
    # line 5 - random seed
    line_5 <- "1234                                 ! Seed for random number generator\n"
    readr::write_file(line_5, filename, append = T)
    # line 6 - allele frequency update?
    line_6 <- "0                                    ! 0/1=Not updating/updating allele frequency\n"
    readr::write_file(x = line_6, file = filename, append = TRUE)
    # line 7 - dioecious/monecious
    line_7 <- "2                                    ! 2/1=Dioecious/Monoecious species\n"
    readr::write_file(x = line_7, file = filename, append = TRUE)
    # line 8 - inbreeding
    line_8 <- "0                                    ! 0/1=No inbreeding/inbreeding\n"
    readr::write_file(x = line_8, file = filename, append = TRUE)
    # line 9 - ploidy
    line_9 <- "0                                    ! 0/1=Diploid species/HaploDiploid species\n"
    readr::write_file(x = line_9, file = filename, append = TRUE)
    # line 10 - mating system 
    line_10 <- paste("0","  ", "0", "                                 ! 0/1=Polygamy/Monogamy for males & females\n", sep = "")
    readr::write_file(x = line_10, file = filename, append = TRUE)
    # line 11 - clone inference 
    line_11 <- paste("0                                    ! 0/1=Clone inference =No/Yes\n", sep = "")
    readr::write_file(x = line_11, file = filename, append = TRUE)
    # line 12 - sibship size scaling 
    line_12 <- "1                                    ! 0/1=Full sibship size scaling =No/Yes\n"
    readr::write_file(x = line_12, file = filename, append = TRUE)
    # line 13 - sibship prior (0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size)
    line_13 <- "0 0 0                                ! 0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size\n"
    readr::write_file(x = line_13, file = filename, append = TRUE)
    # line 14 - pop allele freq indiciator **note - if you set this to 1 you must add other lines inbetween line 14 and line 15 
    line_14 = "0                                    ! 0/1=Unknown/Known population allele frequency\n"
    readr::write_file(x = line_14, file = filename, append = TRUE)
    # line 15 - number of runs
    line_15 <- "1                                    ! Number of runs\n"
    readr::write_file(x = line_15, file = filename, append = TRUE)
    # line 16 - length of run
    line_16 <- "2                                ! Length of run\n"
    readr::write_file(x = line_16, file = filename, append = TRUE)
    # line 17 - 0/1=Monitor method by Iterate
    line_17 <- "0                                    ! 0/1=Monitor method by Iterate\n"
    readr::write_file(x = line_17, file = filename, append = TRUE)
    # line 18 - monitor interval in Iterate
    line_18 <- "10000                                ! Monitor interval in Iterate\n"
    readr::write_file(x = line_18, file = filename, append = TRUE)
    # line 19 - non-Windows version
    line_19 <- "0                                    ! non-Windows version\n"
    readr::write_file(x = line_19, file = filename, append = TRUE)
    # line 20 - analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)
    line_20 <- "1                                    ! Analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)\n"
    readr::write_file(x = line_20, file = filename, append = TRUE)
    # line 21 - 1/2/3=low/medium/high Precision for Full likelihood
    line_21 <- "3                                    ! 1/2/3=low/medium/high Precision for Full likelihood\n"
    readr::write_file(x = line_21, file = filename, append = TRUE)
    # line 22 - marker IDs/Names
    #snp.id = stringi::stri_join(seq(from = 1, to = marker.num, by = 1), collapse = " ")
    line_22 <- paste0("\nmk@","                                    ! Marker IDs\n")
    readr::write_file(x = line_22, file = filename, append = TRUE)
    # line 23 - marker type (codominant/dominant)
    #line_23 = stringi::stri_join(rep(0, marker.num), collapse = " ")
    line_23 = "0@                                    ! Marker Type\n"
    readr::write_file(x = line_23, file = filename, append = TRUE)
    # line 24 - alellelic dropout rate 
    #line_24 <- stringi::stri_join(rep(0, marker.num), collapse = " ")
    line_24 = "0@                                    ! Dropout Rate\n"
    readr::write_file(x = line_24, file = filename, append = TRUE)
    # line 25 - error rate (two new lines)
    #error <- stringi::stri_join(rep(0.2, marker.num), collapse = " ")
    error = paste0("0.02@                                    ! False allele rate\n", "\n")
    #error <- paste0(error, "     ! False allele rate\n\n")
    readr::write_file(x = error, file = filename, append = TRUE)
    # line 26 - offspring id and genotype (two new lines after last inds)
    # change homozygotes (0) to 11
    # change heterozygote (1) to 12
    # change homozygote (2) to 22 
    # change NA to 00
    data = as.matrix(obj)
    data[data == 0] <- "1,1"
    data[data == 1] <- "1,2"
    data[data == 2] <- "2,2"
    data[is.na(data)]<- "0,0"
    
    df = data.frame(data)
    df = purrr:::reduce(seq_along(df), 
                        .init = df, 
                        ~ .x %>% separate(names(df)[.y], 
                                          sep = ',', 
                                          into = paste0(names(df)[.y], '_col_' , seq(1 + max(str_count(df[[.y]], ',')))),
                                          fill = 'right'
                        )
    ) 
    
    utils::write.table(x = df, file = filename, sep = "\t", append = TRUE, col.names = FALSE, row.names = TRUE, quote = FALSE)
    
    # line 27 - Prob. of dad/mum included in the candidates
    line_27 <- "\n\n0  0                                 ! Prob. of dad/mum included in the candidates\n"
    readr::write_file(x = line_27, file = filename, append = TRUE)
    # line 28 - Numbers of candidate males & females
    line_28 <- "0  0                                ! Numbers of candidate males & females\n"
    readr::write_file(x = line_28, file = filename, append = TRUE)
    # line 29 - Number of offspring with known father
    line_29 <- "0  0                                 ! Number of offspring with known father\n"
    readr::write_file(x = line_29, file = filename, append = TRUE)
    # line 30 - Number of offspring with known mother
    line_30 <- "0  0                                 ! Number of offspring with known mother\n"
    readr::write_file(x = line_30, file = filename, append = TRUE)
    # line 31 - Number of known paternal sibships
    line_31 <- "0                                    ! Number of known paternal sibships\n"
    readr::write_file(x = line_31, file = filename, append = TRUE)
    # line 32 - Number of known maternal sibships
    line_32 <- "0                                    ! Number of known paternal sibships\n"
    readr::write_file(x = line_32, file = filename, append = TRUE)
    # line 33 - Number of offspring with known excluded fathers
    line_33 <- "0                                    ! Number of offspring with known excluded fathers\n"
    readr::write_file(x = line_33, file = filename, append = TRUE)
    # line 34 - Number of offspring with known excluded mothers
    line_34 <- "0                                    ! Number of offspring with known excluded fathers\n"
    readr::write_file(x = line_34, file = filename, append = TRUE)
    # line 35 - Number of offspring with known excluded paternal sibships
    line_35 <- "0                                    ! Number of offspring with known excluded fathers\n"
    readr::write_file(x = line_35, file = filename, append = TRUE)
    # line 36 - Number of offspring with known excluded maternal sibships
    line_36 <- "0                                    ! Number of offspring with known excluded fathers\n"
    readr::write_file(x = line_36, file = filename, append = TRUE)
    
    message(paste('Writing', filename, "complete."))
  }
  

}


