# Filtering 
# This document replicates the filtering steps in in Davenport et al., 2021, but using updated methods to refelect availabe packages and versions of the R-packages (now reflected in the updated manuscripts methods). 

# outlier function 
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

# mac function
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

# # -----------------------------------------------------------------------

# STEP 1. 
# Import data 
str_working_dir = './Data/Raw'
str_File_Name_METADATA = 'Davenport_Dart__Report_DSha18_Dups_Relabeled__dartR_METADATA_v1.csv'
str_File_Name_DATA <- 'Report_DSha18-3402_SNP_2_ReLabeled.csv'
gl_filename <- file.path('./Data/Raw/', paste0(str_File_Name_DATA, ".gl.Rdata"))
# # An example of the function used to input data to GENLIGHT object is as follows:
gl <- dartR::gl.read.dart(filename = file.path(str_working_dir, str_File_Name_DATA), ind.metafile = file.path(str_working_dir, str_File_Name_METADATA))
# # Check data file
# gl <- gl.compliance.check(gl)
# # -----------------------------------------------------------------------
# # SAVE GENLIGHT and re-import
# # It is sensible to save your genlight object in binary form using
saveRDS(gl, file = gl_filename)

# and then read it in again with
gl <- readRDS(gl_filename)
# SA is South Australia and presumed to be part of West Aus pop, these are capture locations 
table(gl@pop) 
#NSW      SA SAFRICA      WA 
#245       9      21       3 
# Dataset A 
# This is the full dataset, with all samples and individuals - we need to use this dataset to find possible outlier loci and assign samples to populations/groups
# Filter data the following steps:
#call rate thresholds
IND_MISSING_THRESHOLD <- 1-0.95
LOCI_MISSING_THRESHOLD <- 1-0.70
LOCI_CALLRATE_THRESHOLD <- 0.75
REPRODUCIBILITY_THRESHOLD <- 0.95
MAC_THRESHOLD <- 3 # MAC threshold BELOW which (< MAC_THRESHOLD) loci are removed and above which (>= MAC_THRESHOLD) loci are kept
#read depth
READ_DEPTH_LOWER <- 5
READ_DEPTH_UPPER <- 25
# known Duplicate samples
DUPLICATE_LIST <- c("MBB_1341_Dup","MBB_1574_Dup","MBB_1483_Dup") # known
# identified duplicates RDS name
DUPLICATE_ALL_RDS = "./Data/Raw/duplicate_samples.RDS"
#options for hwe tests
str_each_or_all_pops <- 'all'
n.pop.threshold_val <- 4
float_alpha_hwe <- 0.01 # HWE <0.01 
# pcadapt 
PCADAPT_LOCI_RDS = "./Data/Processed/pcadapt_loci.RDS"
LEA_FILENAME = "./Data/Processed/pcadapt.loci.lfmm"
## -- 
# apply filters 
adegenet::nLoc(gl) # 9841
nInd(gl) # 278
# #--
# Step 1 Reproducibility 
gl_a1 = dartR::gl.filter.reproducibility(gl, threshold = REPRODUCIBILITY_THRESHOLD)
#a1_indx = which(gl@other$loc.metrics$RepAvg >= REPRODUCIBILITY_THRESHOLD)
#gl_1 = gl[,a1_indx] #dartR is not intuitive you have to be careful here, you have to remove loci, samples by name otherwise all the extra data that it stores gets messed up 
adegenet::nLoc(gl_a1) # 9525
gl_1 <- gl.recalc.metrics(gl_a1,  mono.rm = TRUE) 
# -- 
# Step 2 Call Rate 
#a2_indx = which(gl_1@other$loc.metrics$CallRate >= LOCI_CALLRATE_THRESHOLD)
#gl_2 = gl_1[,a2_indx] 
gl_2 = dartR::gl.filter.callrate(gl_1, LOCI_CALLRATE_THRESHOLD)
adegenet::nLoc(gl_2) #  8114
# -- 
# a_3 Read Depth 
gl_3 = dartR::gl.filter.rdepth(gl_2, lower = READ_DEPTH_LOWER, upper = READ_DEPTH_UPPER)
#a3_indx = which(gl_2@other$loc.metrics$rdepth <= READ_DEPTH_UPPER & gl_2@other$loc.metrics$rdepth >= READ_DEPTH_LOWER)
#gl_3 = gl_2[,a3_indx] 
adegenet::nLoc(gl_3) # 7473
# -- 
# Filter MAC 
df_MAC_MAF <- MAC_calcs(gl_3)
# BELOW MAC threshold - Loci removed when filtered
df_MAC_below_threshold <- df_MAC_MAF[df_MAC_MAF$MAC__minor_allele_count < MAC_THRESHOLD,]
#gl_4 = gl_3[,-as.numeric(rownames(df_MAC_below_threshold))]
df_MAC_below_threshold_names = rownames(df_MAC_below_threshold)
keep_MAC_above_threshold_names = locNames(gl_3)[-which(locNames(gl_3) %in% df_MAC_below_threshold_names)]
gl_4 = dartR::gl.keep.loc(gl_3, keep_MAC_above_threshold_names)
gl_4 = dartR::gl.recalc.metrics(gl_4)
nLoc(gl_4) # 7008
# # -- 
# Step 4 Missingness - Individuals 
indmissing_tbl = tibble::tibble(factor = gl_4@pop, sample_name = indNames(gl_4), 
                                value = rowMeans(is.na(tab(gl_4, NA.method = "asis"))))
indmissing_tbl$outlier <- ave(indmissing_tbl$value, indmissing_tbl$factor, FUN = is_outlier) 
ggplot2::ggplot(indmissing_tbl, aes(x = factor, y = value)) +
  geom_boxplot() +
  ggrepel::geom_text_repel(data = indmissing_tbl[indmissing_tbl$outlier == TRUE,], aes(label = sample_name), size = 3, color = "red") +
  theme_minimal()
inds_to_drop_missingness = which(rowMeans(is.na(tab(gl_4, NA.method = "asis"))) >= IND_MISSING_THRESHOLD)
inds_to_drop_missingness
inds_to_keep_missingness = indNames(gl_4)[-which(indNames(gl_4) %in% names(inds_to_drop_missingness))]
# these will be removed: MBB_1338 MBB_1348 MBB_1336 MBB_1455 MBB_1431 
gl_5 = dartR::gl.keep.ind(gl_4, inds_to_keep_missingness, recalc = T, mono.rm = T)
nInd(gl_5) # 273
adegenet::nLoc(gl_5) # 7008
# -- 
# Step 5 Identify (unknown) duplicates and contaminated samples remove them 
# Filter duplicated samples (these were intentionally included, known) an dwe know which ones they are becuase they are labelled as so with "_DUP"
# Also check for accidental or mixed/contaminated samples using heterozygosity ect
# cal het per sample
ind_het = gl.report.heterozygosity(gl_5,method = 'ind')
# I prefer a boxplot, lets zoom in just by making a new one 
df_het = tibble::tibble(value = ind_het$Ho, sample_name = ind_het$ind.name, factor = gl_5@pop[which(indNames(gl_5) %in% ind_het$ind.name)])
# Function to add a column for outliers
df_het$outlier <- ave(df_het$value, df_het$factor, FUN = is_outlier)
p <- ggplot2::ggplot(df_het, aes(x = factor, y = value)) +
  geom_boxplot() +
  ggrepel::geom_text_repel(data = df_het[df_het$outlier == TRUE,], aes(label = sample_name), size = 3, color = "red") +
  theme_minimal()
p
# here we can see some samples have really low heterozygosity relative to other samples
# At this point you need to distinguish between an artifact of poor polymorphism discovery or a biological reason (highly inbred individual, etc.).
missing_tbl = tibble::tibble(missing = rowMeans(is.na(tab(gl_5, NA.method = "asis")))[which(df_het$value <= 0.18)],
                             df_het[which(df_het$value <= 0.18),])
missing_tbl
# so ill leave all of these 
# 6b Remove known duplicates and other identified (possible) duplicates. 
gl_dist_ind <- as.matrix(stats::dist(tab(gl_5, NA.method= "zero"), method = "euclidean"))
df <- reshape2::melt(gl_dist_ind,  id.vars = "Row")
names(df) <- c("row", "col", "Value")
ggplot(df, aes(x=factor(row), y=Value)) +
  geom_boxplot() +
  theme_minimal()
indx = df[df$Value < 30 & df$Value >0,]
indx = indx[!grepl("_Dup", indx$row), ]
indx = indx[!grepl("_Dup", indx$col), ]
# this is clumsy but gets the job done, keep one from each pair, no duplicates
indx$pair <- apply(indx[, c("row", "col")], 1, function(x) paste(sort(x), collapse = "_"))
df_unique <- indx[!duplicated(indx$pair), ] # "MBB_1544"     "MBB_1372"     "MBB_1516"     "MBB_1554" 
DUPLICATE_LIST_ALL = c(as.character(df_unique$row), DUPLICATE_LIST) # just keep one of row or column from df_unique
saveRDS(DUPLICATE_LIST_ALL, file = DUPLICATE_ALL_RDS)
length(DUPLICATE_LIST_ALL) # 7  
no_duplicates_keep = indNames(gl_5)[-which(adegenet::indNames(gl_5) %in% DUPLICATE_LIST_ALL)]
gl_6 = dartR::gl.keep.ind(gl_5, no_duplicates_keep,recalc = T, mono.rm = T)
nLoc(gl_6) # 7008
nInd(gl_6) # 266
table(gl_6@pop)
#NSW      SA SAFRICA      WA 
#235       8      20       3 
# # - - 
# Step 7 HWE
gl_7 = dartR::gl.filter.hwe(
  gl_6,
  subset = str_each_or_all_pops,
  n.pop.threshold = n.pop.threshold_val,
  method_sig = "Exact",
  multi_comp = FALSE,
  alpha_val = float_alpha_hwe,
  pvalue_type = "midp",
  cc_val = 0.5,
  min_sample_size = 2)
nLoc(gl_7) # 7008

# # -- Check population structure, assignment 
# DAPC 
gl_dapc = gl_7
gl_dapc = gl_dapc[-which(gl_dapc@pop == "WA"),]
#gl_dapc@pop[which(gl_dapc@pop == "WA")]<- "SA"
pramx = adegenet::xvalDapc(tab(gl_dapc, NA.method= "mean") , pop(gl_dapc),
                           n.rep = 100,
                           parallel = "multicore", 
                           ncpus = 4L,
                           center = T, scale = F) 
saveRDS(pramx, file = "./Data/Processed/dapcobj.RDS")
# plot
png(file="./Data/Processed/DAPC.plot.png")
scatter.dapc(pramx$DAPC, cex = 2, legend = TRUE,
             clabel = FALSE, posi.leg = "bottomleft", scree.pca = FALSE, cleg = 0.75, xax = 1, yax = 2, 
             inset.solid = 1)
dev.off()
# inds 1446 and 1412 removed, close/group with SA
DAPC_OUTTIES = c("MBB_1412", "MBB_1446") 
# # -- 
# # # PCAdapt
# for this lets not included very related samples. 
# remove samples with relatedness in sibling category (Kinship is relative, not absolute) beta method see Weir BS,Goudet J. 2017. Genetics.
m8 <- as.data.frame(as.matrix(gl_7))
pb8 = hierfstat::beta.dosage(m8) # you can use any kinship or distance metric, must be careful not to remove real related samples (i.e. MZ twins = 0.3) 
diag(pb8)<- 0
indx8 = data.frame(row = rownames(pb8)[which(pb8 >= 0.1, arr.ind = TRUE)[, 1]], 
                   col = colnames(pb8)[which(pb8 >= 0.1, arr.ind = TRUE)[, 2]]) # this cut off capture all known dups, above this not all known dups were included
# this is clumsy but gets the job done, keep one from each pair, no duplicates
indx8$pair <- apply(indx8[, c("row", "col")], 1, function(x) paste(sort(x), collapse = "_"))
df_unique <- indx8[!duplicated(indx8$pair), ]
gl_pcadapt_in = gl_8[-which(adegenet::indNames(gl_8) %in% df_unique$row),]
dos = tab(gl_pcadapt_in, NA.method= "asis")
dos[is.na(dos)] <- 9
LEA::write.lfmm(dos,LEA_FILENAME)
pcadaptin <- pcadapt::read.pcadapt(LEA_FILENAME, type = "lfmm")
pcadaptout <- pcadapt::pcadapt(input = pcadaptin, K = 20, min.maf = 0.05)
scree = gridExtra::grid.arrange(plot(pcadaptout, option = "screeplot", K = 20) +  ggplot2::theme_classic())
as_tibble(pcadaptout$scores)[,1:2] %>% 
  ggplot(aes(x = .[[1]], y = .[[2]], color = gl_pcadapt_in@pop))  +  
  ggplot2::geom_point(size = 4) + ggplot2::theme_classic() 
# looks like only pc1 seperates SA and others, maybe two
k = 1 #  an integer specifying the number of principal components to retain
pcadaptout<- pcadapt::pcadapt(input = pcadaptin, K = k, min.maf = 0.05)
PVALS = pcadaptout$pvalues # loci with MAF below threshold are not computed, they return NA
PCMAP = as.data.frame(PVALS)
QVALS = qvalue::qvalue(as.numeric(PCMAP$PVALS), pi0 = 1) # pi0 = 1, using pi0 will apply the BH procedure and is conservative
PCMAP$QVALS = QVALS$qvalues
# # find outlier 
# # Benjamini-Hochberg Procedure
padj <- p.adjust(PVALS, method="BH")
# #print (head(padj))
alpha <- 0.05
outliers <- which(padj < alpha) # this is an index, only 521 loci
pcadapt_remove_loci = gl_pcadapt_in@loc.names[outliers] # keep these for later 
saveRDS(pcadapt_remove_loci, file = PCADAPT_LOCI_RDS)

# # -- 
# STEP 2 Get Final dataset for Nb analysis 
# Filter again, this time we only keep the target samples (NSW samples). 
# Use code in filter_script_datasetb.R