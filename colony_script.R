write_colony2 <- function (data, strata = NULL, sample.markers = NULL, pop.select = NULL, 
          allele.freq = NULL, inbreeding = 0, mating.sys.males = 0, 
          mating.sys.females = 0, clone = 0, run.length = 2, analysis = 1, 
          allelic.dropout = 0, error.rate = 0.02, print.all.colony.opt = FALSE, 
          random.seed = NULL, verbose = FALSE, parallel.core = parallel::detectCores() - 
            1, filename = NULL, ...) 
{
  radiator_function_header(f.name = "write_colony", verbose = verbose)
  timing <- radiator_tic()
  on.exit(radiator_toc(timing), add = TRUE)
  on.exit(radiator_function_header(f.name = "write_colony", 
                                   start = FALSE, verbose = verbose), add = TRUE)
  if (missing(data)) 
    rlang::abort("Input file missing")
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_colony_", file.date)
  }
  message("Importing data...")
  data.type <- radiator::detect_genomic_format(data)
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") 
      data %<>% radiator::read_rad(data = .)
    data <- gds2tidy(gds = data, parallel.core = parallel.core)
    data.type <- "tbl_df"
  }
  else {
    if (is.vector(data)) 
      data %<>% radiator::tidy_wide(data = ., import.metadata = TRUE)
  }
  if (rlang::has_name(data, "STRATA") && !rlang::has_name(data, 
                                                          "POP_ID")) {
    data %<>% dplyr::rename(POP_ID = STRATA)
  }
  if (!is.null(strata)) 
    data %<>% join_strata(data = ., strata = strata, pop.id = TRUE)
  if (!rlang::has_name(data, "GT")) {
    data %<>% calibrate_alleles(data = ., gt = TRUE) %$% 
      input
  }
  if (!is.null(pop.select)) {
    data %<>% dplyr::filter(POP_ID %in% pop.select)
    data %<>% filter_monomorphic(data = .)
  }
  if (!is.null(sample.markers)) {
    message("Randomly subsampling ", sample.markers, " markers...")
    markers.list <- dplyr::distinct(data, MARKERS) %>% dplyr::sample_n(tbl = ., 
                                                                       size = sample.markers, replace = FALSE)
    data <- suppressWarnings(dplyr::semi_join(data, markers.list, 
                                              by = "MARKERS"))
    markers.list <- NULL
  }
  message("Generating COLONY file...\n")
  radiator:::radiator_colony(data = data, filename = filename, allele.freq = allele.freq , inbreeding = inbreeding, 
                             mating.sys.males = mating.sys.males, mating.sys.females = mating.sys.females, clone = clone, 
                             run.length = run.length, analysis = analysis, allelic.dropout = allelic.dropout, error.rate = error.rate, 
                             print.all.colony.opt = print.all.colony.opt )
  message("COLONY file(s) written in the working directory")
  return("COLONY file(s) written in the working directory")
}