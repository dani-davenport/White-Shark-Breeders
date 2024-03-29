# adapted from CACA_DAVENPORT_DART_Filtering___Initial_Analyses_gl_filt_2_2_v5.R
# by Andrew Jones

######################################################
# parameters used in the filtering

#call rate thresholds
IND_CALLRATE_THRESHOLD <- 0.85
LOCI_CALLRATE_THRESHOLD <- 0.75
REPRODUCIBILITY_THRESHOLD <- 0.98

# MAC threshold BELOW which (< MAC_THRESHOLD) loci are removed and above which (>= MAC_THRESHOLD) loci are kept
MAC_THRESHOLD <- 3

#read depth
READ_DEPTH_LOWER <- 5
READ_DEPTH_UPPER <- 25

#individuals to remove
DUPLICATE_LIST <- c("MBB_1341_Dup","MBB_1372","MBB_1417","MBB_1455","MBB_1483_Dup","MBB_1544")
XVALDAPC_BAD_SAMPLES <- c("MBB_1412","MBB_1446")

#loci to remove based on PCAS
df_pcadapt_loci_to_remove_final <- read.table(file=file.path('./Data/Raw', "3_pops_WAS_as_SAS.pcadapt.results.outlier_loci.qvalues.CACA_2021.Davenport_Dart_Orig_Loci.gl_dart_0_0_v2.csv"), header=FALSE, sep=",")

#options for hwe tests
str_each_or_all_pops <- 'all'
float_alpha_hwe <- 0.01

#absolute, rejects >=50 and <= -0.5
FIS_THRESHOLD <- 0.5
#various file names used
# IMPORT DaRTSeq data
str_working_dir = './Data/Raw'
str_File_Name_METADATA = 'Davenport_Dart__Report_DSha18_Dups_Relabeled__dartR_METADATA_v1.csv'
str_File_Name_DATA <- 'Report_DSha18-3402_SNP_2_ReLabeled.csv'
str_File_Name_FILTERED <- 'filtered_genotypes.RData'
str_File_Name_FILTERED_GL <- 'filtered_genotypes.csv'
gl_filename <- file.path('./Data/Raw/', paste0(str_File_Name_DATA, ".gl.Rdata"))

#############################################
MAC_calcs <- function(gl){
  
  mtx_gl <- as.matrix(gl)
  
  # FULL DATASET
  df_allele_type_counts <- data.frame(colSums(mtx_gl == 0,na.rm=TRUE),colSums(mtx_gl==1,na.rm=TRUE),colSums(mtx_gl==2,na.rm=TRUE),colSums(is.na(mtx_gl)))
  colnames(df_allele_type_counts) <- c("sum_0_hom_ref","sum_1_het","sum_2_hom_alt","sum_NA_missing")
  
  df_MAC=NULL
  df_MAC <- data.frame(df_allele_type_counts
                       ,sum_ref_allele_count=(colSums(mtx_gl == 0,na.rm=TRUE)*2)+colSums(mtx_gl==1,na.rm=TRUE)
                       ,sum_alt_allele_count=colSums(mtx_gl == 1,na.rm=TRUE)+(colSums(mtx_gl==2,na.rm=TRUE)*2)
  )
  
  df_MAC <- data.frame(df_MAC
                       ,MAC__minor_allele_count=ifelse(df_MAC$sum_ref_allele_count > df_MAC$sum_alt_allele_count, df_MAC$sum_alt_allele_count, df_MAC$sum_ref_allele_count)
                       ,sum_alleles=df_MAC$sum_ref_allele_count + df_MAC$sum_alt_allele_count
  )
  
  df_MAC_MAF <- data.frame(df_MAC, MAF__minor_allele_freq=df_MAC$MAC__minor_allele_count / df_MAC$sum_alleles
  )
  
  return(df_MAC_MAF)
}
#############################################

# # An example of the function used to input data to GENLIGHT object is as follows:
gl <- gl.read.dart(filename = file.path(str_working_dir, str_File_Name_DATA), ind.metafile = file.path(str_working_dir, str_File_Name_METADATA))
# # Check data file
# gl <- gl.compliance.check(gl)
# # -----------------------------------------------------------------------
# # SAVE GENLIGHT and re-import
# # It is sensible to save your genlight object in binary form using
saveRDS(gl, file = gl_filename)

# and then read it in again with
gl <- readRDS(gl_filename)

# -----------------------------------------------------------------------
# Filtering the dataset
# -----------------------------------------------------------------------

# -------------
# DAVENPORT DATASET_Nb - Step 1
# REMOVE POPS - only want NSW 

gl_2 <- gl.drop.pop(
  gl,
  pop.list = c("SA","SAFRICA","WA"),
  as.pop = NULL,
  recalc = TRUE,
  mono.rm = TRUE,
  verbose = 5
)

# FILTER BY INDIV CALLRATE
# You can consider call rate across individuals
gl_3 <- gl.filter.callrate(gl_2, method='ind', threshold = IND_CALLRATE_THRESHOLD, mono.rm = FALSE, recalc = TRUE, plot.out = FALSE)

#filter monomorphic
gl_4 <- gl.filter.monomorphs(gl_3)

# -------------
# DAVENPORT DATASET_Nb - Step 3
# REMOVE OUTLIER SAMPLES IDENTIFIED BY XVALDAPC

#DRP INDIVS
# Samples identified as possibly divergant with filter gl_dart_0_0_v2
gl_5 <- gl.drop.ind(gl_4, XVALDAPC_BAD_SAMPLES, recalc = TRUE, mono.rm = TRUE, verbose = 5)

# -------------
# DAVENPORT DATASET_Nb - Step 4
# FILTER BY REPRODUCIBILITY
gl_6 <- gl.filter.reproducibility(gl_5, threshold = REPRODUCIBILITY_THRESHOLD, plot.out = FALSE)

# -------------
# DAVENPORT DATASET_Nb - Step 5
# Remove MONOMORPHIC LOCI

#filter monomorphic
gl_7 <- gl.filter.monomorphs(gl_6)

# -------------
# DAVENPORT DATASET_Nb - Step 6
# DCB FILTER BY MAC

#Minor allele frequency (MAF) is the frequency at which the second most common allele
#Minor allele count (MAF) is the count of the second most common allele

# 0 - Homozygous reference allele
# 1 - Heterozygous
# 2 - Homozygous alternate allele
# NA - Missing

df_MAC_MAF <- MAC_calcs(gl_7)

# BELOW MAC threshold - Loci removed when filtered
df_MAC_below_threshold <- df_MAC_MAF[df_MAC_MAF$MAC__minor_allele_count < MAC_THRESHOLD,]

# Remove loci by name
# Create env of locus names to remove
env_MAC_loci_to_remove <- c(rownames(df_MAC_below_threshold))

gl_8 <- gl.drop.loc(gl_7, env_MAC_loci_to_remove, verbose = 5)

# -------------
# DAVENPORT DATASET_Nb - Step 7
# FILTER BY LOCUS READ DEPTH
gl_9 <- gl.filter.rdepth(gl_8, lower = 5, upper = 25, verbose = 5, plot.out = FALSE)


# -------------
# DAVENPORT DATASET_Nb - Step 8
# FILTER BY LOCUS CALLRATE

gl_10 <- gl.filter.callrate(gl_9, method='loc', threshold = LOCI_CALLRATE_THRESHOLD, plot.out = FALSE)


# -------------
# DAVENPORT DATASET_Nb - Step 9
# FILTER BY SECONDARY SNPS (i.e. reduce to 1 SNP per tag)
# SNP datasets generated by DArT include fragments with more than one SNP (that is, with secondaries)
# and record them separately with the same CloneID (=AlleleID).
# These multiple SNP loci within a fragment are likely to be linked, and so you may wish to remove secondaries.
#gl.report.secondaries(gl_10)

gl_11 <- gl.filter.secondaries(gl_10)

# -------------
# DAVENPORT DATASET_Nb - Step 11
# FILTER BY INDIV CALLRATE

# Consider call rate across individuals
#gl.report.callrate(gl_temp, method='ind')

gl_12 <- gl.filter.callrate(gl_11, method='ind', threshold = IND_CALLRATE_THRESHOLD, mono.rm = FALSE, recalc = TRUE, plot.out = FALSE)

# -------------
# DAVENPORT DATASET_Nb - Step 10
# Remove KNOWN DUPLICATES

#DRP INDIVS
# What was actually removed by DP when comparing raw Dart NSW samples to DP genpop file samples
#gl_temp = gl.drop.ind(gl, c("MBB_1341_Dup","MBB_1371","MBB_1417","MBB_1455","MBB_1483_Dup"), recalc = TRUE, mono.rm = TRUE, verbose = 5)
# Duplicates (> 96.5% the same) found using radiator and DCB pairwise comparision of genotypes
gl_13 <- gl.drop.ind(gl_12, DUPLICATE_LIST, recalc = TRUE, mono.rm = TRUE, verbose = 5)

# -------------
# DAVENPORT DATASET_Nb - Step 12a
# FILTER BY ADAPTATION LOCI USING GENOTYPE MATRIX WITH PC AND Q-VALUE

# Remove loci from GENLIGHT
env_pcadapt_loci_to_remove <- c(df_pcadapt_loci_to_remove_final[,1])

#DROP LOCI
gl_14 <- gl.drop.loc(gl_13, env_pcadapt_loci_to_remove, verbose = 5)

# -------------
# DAVENPORT DATASET_Nb - Step 12b
# FILTER BY HWE

# Calculates the probabilities of agreement with H-W equilibrium based on observed frequencies
# of reference homozygotes, heterozygotes and alternate homozygotes.
# Uses the exact calculations contained in function utils.prob.hwe() as developed by Wigginton et al. (2005).
#str_each_or_all_pops= c('NSW','SAS', 'SAF')
#str_each_or_all_pops='each'
# 
# bool_multi_comp = FALSE
# 
# # This is ignored if bool_multi_comp = FALSE
# str_multi_comp_method = "bonferroni"
# 
#   df_hwe_diagnotics = gl.diagnostics.hwe(
#     gl,
#     alpha_val = float_alpha_hwe,
#     bins = 20,
#     #colors_hist = two_colors,
#     #colors_barplot = two_colors_contrast,
#     plot_theme = theme_dartR(),
#     save2tmp = FALSE,
#     n.cores = 28
#   )
# 
#   str_each_or_all_pops
#   
#   df_hwe_report = gl.report.hwe(
#     gl,
#     subset = str_each_or_all_pops,
#     method_sig = "Exact",
#     multi_comp = bool_multi_comp,
#     multi_comp_method = str_multi_comp_method,
#     alpha_val = float_alpha_hwe,
#     pvalue_type = "midp",
#     cc_val = 0.5,
#     sig_only = TRUE,
#     min_sample_size = 3,
#     plot.out = TRUE,
#     #plot_colors = two_colors_contrast,
#     max_plots = 4,
#     save2tmp = FALSE)
# }

# Remove loci that are out of HWE
#gl_temp = gl.filter.hwe(gl, alpha = 0.01, basis = "any", bon = FALSE)

gl_15 <- gl.filter.hwe(
  gl_14,
  subset = str_each_or_all_pops,
  method_sig = "Exact",
  multi_comp = FALSE,
  alpha_val = float_alpha_hwe,
  pvalue_type = "midp",
  cc_val = 0.5,
  min_sample_size = 2
)


# -------------
# DAVENPORT DATASET_Nb - Step 13
# FILTER BY FIS
# -----------------
# DCB Method:

df_basic_stats <- gl.basic.stats(gl_15, digits = 4)

#head(df_basic_stats$perloc)

df_perloc <- data.frame(df_basic_stats$perloc)
df_fis_not_within_threshold <- df_perloc[abs(df_perloc$Fis) > FIS_THRESHOLD, ]

env_loci_to_remove <- row.names(df_fis_not_within_threshold)

#find and replace certain chars
env_loci_to_remove <-  gsub("X","",env_loci_to_remove)
env_loci_to_remove <-  sub("[.]","-",env_loci_to_remove)
env_loci_to_remove <-  sub("[.]","-",env_loci_to_remove)
env_loci_to_remove <-  sub("[.]","/",env_loci_to_remove)

#DROP LOCI
gl_final <- gl.drop.loc(gl_15, env_loci_to_remove, verbose = 5)

gl_final
# -----------------------------------------------------------------------
# SAVE FILTERED GENLIGHT 

# It is sensible to save your genlight object in binary form using
saveRDS(gl_final, file = file.path('./Data/Processed', str_File_Name_FILTERED))


