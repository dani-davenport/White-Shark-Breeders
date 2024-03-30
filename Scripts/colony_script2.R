radiator_colony2 <- function (data, filename, allele.freq = NULL, inbreeding = 0, pop.select=NULL,
          mating.sys.males = 0, mating.sys.females = 0, clone = 0, 
          run.length = 2, analysis = 1, allelic.dropout = 0, error.rate = 0.02, 
          print.all.colony.opt = FALSE, ...) 
{
  data %<>% dplyr::select(MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
    dplyr::mutate(A1 = stringi::stri_sub(GT, 1, 3), A2 = stringi::stri_sub(GT, 
                                                                           4, 6), GT = NULL) %>% radiator::rad_long(x = ., cols = c("POP_ID", 
                                                                                                                                    "INDIVIDUALS", "MARKERS"), names_to = "ALLELE_GROUP", 
                                                                                                                    values_to = "ALLELES") %>% dplyr::mutate(ALLELES = as.numeric(ALLELES)) %>% 
    dplyr::arrange(MARKERS, POP_ID, INDIVIDUALS)
  if (!is.null(allele.freq)) {
    message("Computing allele frequency")
    if (allele.freq != "overall") {
      input.alleles <- dplyr::filter(.data = data, POP_ID %in% 
                                       allele.freq) %>% dplyr::filter(ALLELES != 0)
    }
    else {
      input.alleles <- dplyr::filter(.data = data, ALLELES !=  0)
    }
    allele.per.locus <- dplyr::distinct(input.alleles, MARKERS, 
                                        ALLELES) %>% dplyr::count(x = ., MARKERS) %>% dplyr::arrange(MARKERS) %>% 
      dplyr::select(n) %>% purrr::flatten_chr(.)
    freq <- input.alleles %>% dplyr::group_by(MARKERS, ALLELES) %>% 
      dplyr::tally(.) %>% dplyr::ungroup() %>% dplyr::group_by(MARKERS) %>% 
      dplyr::mutate(FREQ = round(n/sum(n), 2)) %>% dplyr::select(MARKERS, 
                                                                 ALLELES, FREQ) %>% dplyr::arrange(MARKERS) %>% dplyr::mutate(GROUP = seq(1, 
                                                                                                                                          n(), by = 1)) %>% dplyr::mutate(dplyr::across(tidyselect::everything(), 
                                                                                                                                                                                        .fns = as.character)) %>% radiator::rad_long(x = ., 
                                                                                                                                                                                                                                     cols = c("GROUP", "MARKERS"), names_to = "ALLELES_FREQ", 
                                                                                                                                                                                                                                     values_to = "VALUE") %>% dplyr::mutate(ALLELES_FREQ = factor(ALLELES_FREQ, 
                                                                                                                                                                                                                                                                                                  levels = c("ALLELES", "FREQ"), ordered = TRUE)) %>% 
      radiator::rad_wide(x = ., formula = "MARKERS + ALLELES_FREQ ~ GROUP", 
                         values_from = "VALUE") %>% tidyr::unite(data = ., 
                                                                 col = INFO, -c(MARKERS, ALLELES_FREQ), sep = " ") %>% 
      dplyr::mutate(INFO = stringi::stri_replace_all_regex(str = INFO, 
                                                           pattern = "NA", replacement = "", vectorize_all = FALSE), 
                    INFO = stringi::stri_trim_right(str = INFO, pattern = "\\P{Wspace}")) %>% 
      dplyr::select(-MARKERS, -ALLELES_FREQ)
    #input.alleles <- NULL
  }else{
    input.alleles <- dplyr::filter(.data = data, ALLELES !=  0)
  }
  
  if (!is.null(pop.select)) {
    data %<>% dplyr::filter(POP_ID %in% pop.select)
  }
  
  markers.name <- input.alleles %>% dplyr::distinct(MARKERS) %>% dplyr::arrange(MARKERS) %$% MARKERS
  marker.num <- length(markers.name)
  data <- tidyr::unite(data = data, col = MARKERS.ALLELE_GROUP, 
                       MARKERS, ALLELE_GROUP, sep = ".") %>% radiator::rad_wide(x = ., 
                                                                                formula = "POP_ID + INDIVIDUALS ~ MARKERS.ALLELE_GROUP", 
                                                                                values_from = "ALLELES") %>% dplyr::arrange(POP_ID, INDIVIDUALS) %>% 
    dplyr::select(-POP_ID) %>% dplyr::mutate(dplyr::across(tidyselect::everything(), 
                                                           .fns = as.character))
  dataset.opt <- "`My first COLONY run`                ! Dataset name\n"
  readr::write_file(x = dataset.opt, file = filename, append = FALSE)
  colony.output.filename <- stringi::stri_replace_all_fixed(filename, 
                                                            pattern = ".dat", replacement = "")
  colony.output.filename <- paste(colony.output.filename, "         ! Output file name\n")
  readr::write_file(x = colony.output.filename, file = filename, 
                    append = TRUE)
  off.num.opt <- paste(nrow(data), "                                  ! Number of offspring in the sample\n", 
                       sep = "")
  readr::write_file(x = off.num.opt, file = filename, append = TRUE)
  marker.num.opt <- paste(marker.num, "                                 ! Number of loci\n", 
                          sep = "")
  readr::write_file(x = marker.num.opt, file = filename, append = TRUE)
  seed.opt <- "1234                                 ! Seed for random number generator\n"
  readr::write_file(x = seed.opt, file = filename, append = TRUE)
  update.allele.freq.opt <- "0                                    ! 0/1=Not updating/updating allele frequency\n"
  readr::write_file(x = update.allele.freq.opt, file = filename, 
                    append = TRUE)
  dioecious.opt <- "2                                    ! 2/1=Dioecious/Monoecious species\n"
  readr::write_file(x = dioecious.opt, file = filename, append = TRUE)
  inbreeding.opt <- paste(inbreeding, "                                    ! 0/1=No inbreeding/inbreeding\n", 
                          sep = "")
  readr::write_file(x = inbreeding.opt, file = filename, append = TRUE)
  ploidy.opt <- "0                                    ! 0/1=Diploid species/HaploDiploid species\n"
  readr::write_file(x = ploidy.opt, file = filename, append = TRUE)
  mating.opt <- paste(mating.sys.males, "  ", mating.sys.females, 
                      "                                 ! 0/1=Polygamy/Monogamy for males & females\n", 
                      sep = "")
  readr::write_file(x = mating.opt, file = filename, append = TRUE)
  clone.opt <- paste(clone, "                                    ! 0/1=Clone inference =No/Yes\n", 
                     sep = "")
  readr::write_file(x = clone.opt, file = filename, append = TRUE)
  sib.size.scal.opt <- "1                                    ! 0/1=Full sibship size scaling =No/Yes\n"
  readr::write_file(x = sib.size.scal.opt, file = filename, 
                    append = TRUE)
  sib.prior.opt <- "0 0 0                                ! 0, 1, 2, 3 = No, weak, medium, strong sibship size prior; mean paternal & maternal sibship size\n"
  readr::write_file(x = sib.prior.opt, file = filename, append = TRUE)
  if (is.null(allele.freq)) {
    allele.freq.ind.opt <- "0                                    ! 0/1=Unknown/Known population allele frequency\n"
  }
  else {
    allele.freq.ind.opt <- "1                                    ! 0/1=Unknown/Known population allele frequency\n"
  }
  readr::write_file(x = allele.freq.ind.opt, file = filename, 
                    append = TRUE)
  if (!is.null(allele.freq)) {
    readr::write_file(x = paste0(paste0(allele.per.locus, 
                                        collapse = " "), "  !Number of alleles per locus\n"), 
                      file = filename, append = TRUE)
    utils::write.table(x = freq, file = filename, append = TRUE, 
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  num.run.opt <- "1                                    ! Number of runs\n"
  readr::write_file(x = num.run.opt, file = filename, append = TRUE)
  run.length.opt <- paste(run.length, "                                    ! Length of run\n", 
                          sep = "")
  readr::write_file(x = run.length.opt, file = filename, append = TRUE)
  monitor.met.opt <- "0                                    ! 0/1=Monitor method by Iterate\n"
  readr::write_file(x = monitor.met.opt, file = filename, append = TRUE)
  monitor.int.opt <- "10000                                ! Monitor interval in Iterate\n"
  readr::write_file(x = monitor.int.opt, file = filename, append = TRUE)
  windows.gui.opt <- "0                                    ! non-Windows version\n"
  readr::write_file(x = windows.gui.opt, file = filename, append = TRUE)
  analysis.opt <- paste(analysis, "                                    ! Analysis 0 (Pairwise-Likelihood Score), 1 (Full Likelihood), 2 (combined Pairwise-Likelihood Score and Full Likelihood)\n", 
                        sep = "")
  readr::write_file(x = analysis.opt, file = filename, append = TRUE)
  precision.opt <- "3                                    ! 1/2/3=low/medium/high Precision for Full likelihood\n"
  readr::write_file(x = precision.opt, file = filename, append = TRUE)
  #markers.name.opt <- stringi::stri_join(markers.name, collapse = " ")
  #markers.name.opt <- stringi::stri_join(markers.name.opt,  " ! Marker IDs\n")
  markers.name.opt <- "MK@ ! Marker IDs\n" #stringi::stri_join(markers.name.opt,  " ! Marker IDs\n")
  readr::write_file(x = markers.name.opt, file = filename,   append = TRUE)
  marker.type.opt <- stringi::stri_join(rep(0, marker.num), 
                                        collapse = " ")
  marker.type.opt <- stringi::stri_join(marker.type.opt, "  ! Marker types, 0/1=Codominant/Dominant\n")
  readr::write_file(x = marker.type.opt, file = filename, append = TRUE)
  dropout <- stringi::stri_join(rep(allelic.dropout, marker.num), 
                                collapse = " ")
  dropout <- stringi::stri_join(dropout, "     ! Allelic dropout rate at each locus\n")
  readr::write_file(x = dropout, file = filename, append = TRUE)
  error <- stringi::stri_join(rep(error.rate, marker.num), 
                              collapse = " ")
  error <- stringi::stri_join(error, "     ! False allele rate\n\n")
  readr::write_file(x = error, file = filename, append = TRUE)
  utils::write.table(x = data, file = filename, sep = " ", 
                     append = TRUE, col.names = FALSE, row.names = FALSE, 
                     quote = FALSE)
  prob.opt <- "\n\n0  0                                 ! Prob. of dad/mum included in the candidates\n"
  readr::write_file(x = prob.opt, file = filename, append = TRUE)
  candidate.opt <- "0  0                                ! Numbers of candidate males & females\n"
  readr::write_file(x = candidate.opt, file = filename, append = TRUE)
  if (print.all.colony.opt) {
    message("Printing all COLONY options...")
    candidate.male.id.geno.opt <- "!Candidate male ID and genotypes\n"
    readr::write_file(x = candidate.male.id.geno.opt, file = filename, 
                      append = TRUE)
    candidate.female.id.geno.opt <- "!Candidate female ID and genotypes\n"
    readr::write_file(x = candidate.female.id.geno.opt, file = filename, 
                      append = TRUE)
  }
  known.paternity.opt <- "0  0                                 ! Number of offspring with known father\n"
  readr::write_file(x = known.paternity.opt, file = filename, 
                    append = TRUE)
  if (print.all.colony.opt) {
    known.father.dyad.opt <- "! Offspring ID and known father ID (Known offspring-father dyad)\n"
    readr::write_file(x = known.father.dyad.opt, file = filename, 
                      append = TRUE)
  }
  known.maternity.opt <- "0  0                                 ! Number of offspring with known mother\n"
  readr::write_file(x = known.maternity.opt, file = filename, 
                    append = TRUE)
  if (print.all.colony.opt) {
    known.mother.dyad.opt <- "! Offspring ID and known mother ID (Known offspring-mother dyad)\n"
    readr::write_file(x = known.mother.dyad.opt, file = filename, 
                      append = TRUE)
  }
  known.paternal.sibships.opt <- "0                                    ! Number of known paternal sibships\n"
  readr::write_file(x = known.paternal.sibships.opt, file = filename, 
                    append = TRUE)
  if (print.all.colony.opt) {
    known.paternal.sibships.size.opt <- "! Paternal sibship size and members\n"
    readr::write_file(x = known.paternal.sibships.size.opt, 
                      file = filename, append = TRUE)
  }
  known.maternal.sibships.opt <- "0                                    ! Number of known maternal sibships\n"
  readr::write_file(x = known.maternal.sibships.opt, file = filename, 
                    append = TRUE)
  if (print.all.colony.opt) {
    known.maternal.sibships.size.opt <- "! Maternal sibship size and members\n"
    readr::write_file(x = known.maternal.sibships.size.opt, 
                      file = filename, append = TRUE)
  }
  offspring.known.excl.paternity.opt <- "0                                    ! Number of offspring with known excluded fathers\n"
  readr::write_file(x = offspring.known.excl.paternity.opt, 
                    file = filename, append = TRUE)
  if (print.all.colony.opt) {
    excl.paternity.opt <- "! Offspring ID, number of excluded fathers, and excluded father IDs\n"
    readr::write_file(x = excl.paternity.opt, file = filename, 
                      append = TRUE)
  }
  offspring.known.excl.maternity.opt <- "0                                    ! Number of offspring with known excluded mothers\n"
  readr::write_file(x = offspring.known.excl.maternity.opt, 
                    file = filename, append = TRUE)
  if (print.all.colony.opt) {
    excl.maternity.opt <- "! Offspring ID, number of excluded mothers, and excluded father IDs\n"
    readr::write_file(x = excl.maternity.opt, file = filename, 
                      append = TRUE)
  }
  offspring.known.excl.paternal.sibships.opt <- "0                                    ! Number of offspring with known excluded paternal sibships\n"
  readr::write_file(x = offspring.known.excl.paternal.sibships.opt, 
                    file = filename, append = TRUE)
  if (print.all.colony.opt) {
    excluded.paternal.siblings.opt <- "! Excluded paternal siblings\n"
    readr::write_file(x = excluded.paternal.siblings.opt, 
                      file = filename, append = TRUE)
  }
  offspring.known.excl.maternal.sibships.opt <- "0                                    ! Number of offspring with known excluded maternal sibships\n"
  readr::write_file(x = offspring.known.excl.maternal.sibships.opt, 
                    file = filename, append = TRUE)
  if (print.all.colony.opt) {
    excluded.maternal.siblings.opt <- "! Excluded maternal siblings\n"
    readr::write_file(x = excluded.maternal.siblings.opt, 
                      file = filename, append = TRUE)
  }
}