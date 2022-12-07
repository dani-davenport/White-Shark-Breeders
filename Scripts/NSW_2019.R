#############################
#setup
library(dplyr)
library(ggplot2)
#############################
# make a STRATA file with only NSW samples
strata<- read.csv("/Users/danielle/Google Drive/UQ_PhD/Lab_Book/ECWS_NSW_Fisheries/Data/Analysis/ByPutPop_Strata.tsv", header = T, sep = "\t")
# make a file using all individuals in NSW ONLY, grouped by YOB and written into a single genepop with "multiple" populations/cohorts represented in the file 
NSW<- strata[which(strata$STRATA == "NSW"),]
write.table(NSW, file = "Strata_NSW_ONLY.tsv", quote = F, sep = "\t", col.names = T, row.names = F)
#############################
# calculate stats for filtering SNP data
##{r Calculting Filtering Stat D}
x<- read.csv("/Users/danielle/Google Drive/UQ_PhD/Lab_Book/ECWS_NSW_Fisheries/Data/Report-DSha18-3402/Report_DSha18-3402_SNP_singlerow_2.csv")
dart_check<- x[7:nrow(x),]## number is the row number which corresponds to where the AlleleSequence data begins
x.names<- as.character(unlist(x[6,]))
colnames(dart_check)<- x.names
## calculate your coverage statistic d
## A maximum read depth equal to  d+4*√d, whereby d is the average read depth. This step reduces the number of false heterozygotes due to sequencing errors or due to the presence of paralogs (see Li, 2014). 
#' Having a higher than "normal" SNP number is usually the results of
#' assembly artifacts or bad assembly parameters.
# filter on AvgDepthCoverage
coverage<- dart_check[,1:21] %>%
  dplyr::rename(REF = AvgCountRef, ALT = AvgCountSnp) %>%
  dplyr::mutate(READS = as.numeric(REF) + as.numeric(ALT))
# find average
av<- mean(coverage$READS)
## [1] 11.58773
# find max
max(as.numeric(as.character(coverage$READS)))
## [1] 117.0637
# find min
min(as.numeric(as.character(coverage$READS)))
## [1] 5.04819
# calculate d, you need this value later for filtering 
d<- av+4*sqrt((av))
d
##[1] 25.20403
##

# filter with RADIATOR
## check read-depth (see section A above)
radiator::filter_rad(interactive.filter = T, data = "/Users/danielle/Google Drive/UQ_PhD/Lab_Book/ECWS_NSW_Fisheries/Data/Report-DSha18-3402/Report_DSha18_Dups_Relabeled.csv", strata = "Strata_NSW_ONLY.tsv") # this strata file contains only inds sampled from NSW that did not group out in the DAPC (inds 1446 and 1412 removed)
# didn't do HWE in radiator. I do below using the HW package
# individuals removed who were intentionally duplciated, or looked to have a higher distance score then duplicated indivudals (see list "Pairwise_Dists" in main folder)
# some samples on plates looked very simialr (close pairwise distance) and were confirmed to be plated in ajoining wells, may have been an error in lab work so samples were removed. 
########### 
# HW filter
# we cannot have loci out of HWE for NEEstimator or COLONY 
radiator::filter_hwe(data = tidy_genomic_data("filter_rad_20190624@1024/15_filtered/radiator_data_20190624@1038.rad", strata = "Strata_NSW_ONLY.tsv"))
# remove outlier loci id'd in ALL dataset 
outlier<- read.table("~/Google Drive/UQ_PhD/Lab_Book/ECWS_NSW_Fisheries/Data/Analysis/2_Filtering/Filter_1/Filter_1.1/ALL_Pops/pcadapt/PCAdapt_Add_TO_BLACKLIST_LOCI2019.tsv", header = T)
whitelist<- read.table("24_filter_hwe_20190624@1133/whitelist.markers.hwe.0.05.mid.p.value.1.hw.pop.threshold.tsv", header = T)
writeout<- whitelist %>%
   # also remove outlier loci if remaning id'd by PCadapt in "All" Filtering 
  filter(!MARKERS %in% outlier$MARKERS) 
nrow(writeout)
colnames(writeout)<- c("MARKERS","CHROM","LOCUS","POS")
#[1] 3736
write.table(writeout, file = "Whitelist_HWD_0.1.nooutlier.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
# write genepop
#radiator::write_genepop(filename = "Genepop_Files/NSW_NO_HWD_OUTLIER", radiator::tidy_genomic_data("24_filter_hwe_20190624@1133/tidy.filtered.hwe.0.05.mid.p.value.1.hw.pop.threshold.rad", strata = "Strata_Sex_Model_COHORT_REGRESSION.tsv", whitelist.markers = "Whitelist_HWD_0.1.nooutlier.tsv"))
## check for FIS
library(hierfstat)
library(tibble)
library(dplyr)
td2.dart<- radiator::tidy_genomic_data(radiator::tidy_genomic_data("filter_rad_20190624@1024/15_filtered/radiator_data_20190624@1038.rad", strata = "Strata_Sex_Model_COHORT_REGRESSION.tsv", whitelist.markers = "Whitelist_HWD_0.1.nooutlier.tsv"))
radiator::write_genepop(td2.dart, filename = "Genepop_Files/NSW_NO_HWD_OUTLIER_")

## PCA
tidy_gi<- radiator::write_genind(td2.dart) # blacklisted ids have been removed from the STRATA file 
gi.scaled<- adegenet::scaleGen(tidy_gi, NA.method="zero")
pca<- ade4::dudi.pca(df = gi.scaled, scannf = FALSE, nf = 100, center = T, scale = T) # number of axis kept

factoextra::fviz_eig(pca)
factoextra::fviz_pca_ind(pca, axes = c(1,2))

#Tidy genomic data:
 #Number of markers: 3693
 #Number of chromosome/contig/scaffold: 1
 #Number of individuals: 241
spall<- radiator::write_hierfstat(radiator::tidy_genepop("Genepop_Files/NSW_NO_HWD_OUTLIER__genepop.gen"))
bootfis<- hierfstat::basic.stats(spall)
fisvals<- as.data.frame(bootfis$perloc$Fis)
rownames(fisvals)<- rownames(bootfis$perloc)
row_order <- order(fisvals$`bootfis$perloc$Fis`, decreasing = FALSE)
x<- tibble::rownames_to_column(fisvals[row_order, , drop = FALSE])
colnames(x)<- c("MARKERS", "fisvals")
#lociouthwe<- read.table("NSW/14_filter_hwe_20190112@0003/blacklist.markers.hwd.0.001.mid.p.value.1.hw.pop.threshold.tsv", header = T)
#to_plot_fis<- dplyr::left_join(x, lociouthwe, "MARKERS") %>%
#  mutate(LOCUS = factor(ifelse(is.na(LOCUS), 1, 2)))
##to_plot_fis$loci = factor(to_plot_fis$loci, levels=to_plot_fis[order(to_plot_fis$loci), "loci"])
## plot
to_plot_fis<- x
library(ggplot2)
dplyr::filter(to_plot_fis, !is.na(fisvals)) %>%
  ggplot(aes(reorder(MARKERS, fisvals), fisvals)) + #color = LOCUS
  geom_point(aes(alpha = 1/10)) + 
  #scale_size(range=c(0,0.2)) +
  #scale_x_discrete(expand = c(.1, .1)) + # padding x axis which is chopping off points
  geom_hline(yintercept=0)+ 
  #labs(size = "Missing Data") + scale_alpha(guide = 'none') +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), panel.background = element_blank(), axis.line.y = element_line(color="black"))
# some greater then and less than 0.5
idx<- which(x$fisvals <=-0.5 | x$fisvals >=0.5)
blkFIS<- x$MARKERS[idx]
b1<- do.call("rbind", strsplit(blkFIS , "__" ), quote = F)
markers<- do.call("rbind", strsplit(blkFIS , " " ), quote = F)
b2<- cbind(markers[,2], b1)
colnames(b2)<- c("MARKERS","CHROM","LOCUS","POS")
b2<- as_tibble(b2)
nrow(b2)
#25
which(td2.dart$MARKERS %in% b2$MARKERS)
write.table(b2, file = "Blacklist_FIS", quote = F, row.names = F, col.names = T, sep = "\t")
b2<- read.table("Blacklist_FIS", header = T)
filtered_dart2<- td2.dart %>%
  dplyr::filter(!MARKERS %in% b2$MARKERS)
radiator::data_info(filtered_dart2)
# loci 3668
radiator::write_genepop(filtered_dart2, filename = "Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier") # DATASET for Ne Estimator
##########################
# how many in each group 
strata = as_tibble(read.table("Strata_Sex_Model_COHORT_REGRESSION.tsv", header = T))
write.table(strata %>%
  count(STRATA), row.names = F, col.names = T, sep = "\t", file = "counts_per_cohort.tsv", quote = F)

#Interquartile Range for LDNe, COLONY ect

##summarise SNP DATA for publication
#p<- poppr::poppr(td2.dart)
gi.dart<- radiator::write_genind(radiator::tidy_genepop("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop.gen"))
#div<- adegenet::summary(gi.dart)
#div = tibble::as_tibble(cbind(Hobs = div$Hobs, Hexp = div$Hexp))
#ggplot(div, aes(x= Hobs, y = Hexp)) + 
#  geom_point() + 
#  theme_bw()# title = "Expected heterozygosity as a function of observed heterozygosity per locus"
htd2<- radiator::write_hierfstat(radiator::tidy_genepop("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop.gen"))
bs<- hierfstat::basic.stats(htd2)
diveRsity::divBasic("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop.gen", outfile = "diveBasic_Stats")
# check 
radiator::allele_count(radiator::tidy_genepop("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop.gen"))
# write out stats
write.table(bs$overall, file = "hierfstat_summary_stats_snp", sep = "\t", quote = F, row.names = T)

## WRITE COLONY FILES 
# write one with all allele freqs (calculated from all samples)
# then write each seperatley (per cohort)
# edit mannually 
# write cohorts seperatley 
pop_info<- tibble::tibble(pop = filtered_dart2$POP_ID)
pop_info_series<- sort(levels(as.factor(pop_info$pop)),decreasing = T)

for (i in 1:length(pop_info_series)){
  #z = septd2[[i]]
  fname = pop_info_series[i]
  print(fname)
  #file.create(paste0("SepCohorts/", fname))
  #radiator::write_colony(radiator::tidy_genind(z), filename = paste0("COLONY/SepCohorts/", fname))
  # can try to subset with allelle freqs
  radiator::write_colony(data = filtered_dart2, filename = paste0("COLONY/SepCohorts/Allele_Freqs/", fname), print.all.colony.opt = T, allele.freq = pop_info_series, pop.select = c(fname))
}
#this removes loci and doesnt calculate allele freqs
# get for all samples as one "pop" and mannually curate later
radiator::write_colony(filter_monomorphic(data =radiator::tidy_genepop("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop-ALL.gen")), filename = paste0("COLONY/SepCohorts/Allele_Freqs/", "ALL"), allele.freq = "overall")
af<- hierfstat::pop.freq(radiator::write_hierfstat(filter_monomorphic(data =radiator::tidy_genepop("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop-ALL.gen"))))
## format for coancestory 
# format is: two rows per locus, one coloumn for each allele, where the name of the alelle (an integer) is in row 1, the frequecy in row 2
file.date <- format(Sys.time(), "%Y%m%d@%H%M")
filepath <- stringi::stri_join("allele_freqs_coancestory_", file.date)
# get write format 
df<- do.call(rbind, lapply(af, data.frame))
df<- rownames_to_column(df)
colnames(df)<- c("LOCUS", "Allele", "NA", "Freq")
same_locus<- stringr::str_split(df$LOCUS, "[.]")
df$LOCUS_NAME<- sapply(same_locus, "[[", 2)

# this function is applied row-wise, it will write out each loci with 2 alleles with name in one line (eg row 1), followed by frequency data in following line (eg row2 )
# make data wide
df2<- df[,c("LOCUS_NAME", "Allele", "Freq")]
# not working
# df3<- df2 %>%
#   group_by(LOCUS_NAME) %>%
#   mutate(id = row_number()) %>%
#   gather(key, value, -(LOCUS_NAME:Freq), id) %>%
#   unite(id_key, key, value) %>%
#   gather(Allele_VAR, Allele_Value, -c(LOCUS_NAME, Freq)) %>%
#   gather(Freq_VAR, Freq_Value, -c(LOCUS_NAME, Allele_VAR, Allele_Value)) %>%
#   spread(key = id_key, value = c(Allele_Var, Freq_VAR))

x<- as_tibble(df2) %>%
    mutate(Allele = as.numeric(as.character(Allele)))
    for(i in seq(1, length(x$LOCUS_NAME), by = 2)){
      z = i + 1
      print(i)
      locusName <- x$LOCUS_NAME[i]
      locusNameVector<- c(paste0(locusName, "1"), paste0(locusName, "2"))
      freqs<- c(x$Freq[[i]], x$Freq[[i+1]])
      allele_id<- c(x$Allele[[i]], x$Allele[[z]])
      cat(paste(allele_id), file = paste(file.date, filepath, sep= "_"), append = T, fill = T)
      cat(paste(freqs), file= paste(file.date, filepath, sep= "_"), append = T, fill = T)}

# remove weird end of line in vim ":%s/\r//g" or sed $'s/$/\r/'
## write out emprical error and missing rate
out<- filtered_dart2 %>% 
  group_by(MARKERS) %>% 
  summarize(callRate=1-mean(CALL_RATE), avgError=1-mean(REP_AVG)) %>% ## call rate accounting for missingness, repavg accounting for genotyping error
  rowwise() %>%
  mutate(dropOut = c(rep(0)))
file.date <- format(Sys.time(), "%Y%m%d@%H%M")
filepath <- stringi::stri_join("allele_error_", file.date)
write.table(out[,2:4], quote = F, row.names = F, file = paste0(filepath, ".tsv"), col.names = F)










# FOR LDNb...get all samples within Q1 and Q3 and check the effects of including sample from another YOB in the set.....

age<- as.tibble(read.csv("../../../../Age_Samples/Cohort_Info_No_L0_Regression.csv"))
# split into cohorts
x<- age %>%
  group_by(sex_model_cohort) %>%
  mutate(Q2 = quantile(Predicted_TL)[2]) %>%
  mutate(Q3 = quantile(Predicted_TL)[4]) %>%
  filter(Predicted_TL >= Q2 & Predicted_TL <= Q3) 


ages<- x %>%
  group_by(sex_model_cohort) %>%
  tally()

# plot 
bp <- ggplot(x, aes(y=Predicted_TL, x=as.factor(sex_model_cohort))) + 
  geom_boxplot() + 
  xlab(label = "\nCohort") + 
  ylab(label = "Estimated Total Length\n") + 
  theme_bw() + 
  theme(legend.position="none", text = element_text(size=12), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

bp

# make a strata file
strataz<- x %>%
  select(MBB_Code, sex_model_cohort) %>%
  mutate(INDIVIDUALS = MBB_Code) %>%
  select(MBB_Code, INDIVIDUALS, sex_model_cohort) %>%
  ungroup()
colnames(strataz)<- c("TARGET_ID",	"INDIVIDUALS",	"STRATA")
strataz$INDIVIDUALS<- c(stringr::str_replace_all(strataz$INDIVIDUALS, " ", "-"))
strataz$TARGET_ID<- c(stringr::str_replace_all(strataz$TARGET_ID, " ", "-"))
dart_ids<- radiator::extract_dart_target_id("../../../../Report-DSha18-3402/Report_DSha18-3402_SNP_2_ReLabeled.csv")
out<- inner_join(strataz, dart_ids)
write.table(out, file = "Q2_Q3_STRATA.tsv", quote = F, sep = "\t", row.names = F)

ids_keep<- tibble::as_tibble(read.table("Q2_Q3_STRATA.tsv", header = T))
filtered_dart <- filtered_dart2 %>%
  left_join(ids_keep, ., by = "INDIVIDUALS") 

radiator::data_info(filtered_dart) #2566 loci, 126 inds
# # make new strata
# x = c()
# for(i in ids_keep$STRATA){
#   x = c(x, rep(i, radiator::data_info(filtered_dart)$n.locus))
# }
# #add in the year strata 
radiator::write_genepop(filtered_dart, filename = "Q2_Q3")

###################
# Sensitivity analsysis to show we believe there is little effect of linked loci, where we would expect that linked loci included in a sample would downwardly bias the estimated of Nb (?) 
# subsample the loci N=1000, 1000 times, only for the dataset with the most individuals, and reesimate Nb for that
## following on from above
setwd("~/Google Drive/UQ_PhD/Lab_Book/ECWS_NSW_Fisheries/Data/Analysis/2_Filtering/Filter_1/Filter_1.1/")
dart<- tidy_genepop("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop.gen")
dart_markers<- tidy_genepop("Genepop_Files/NSW_Sex_Regression_No_FIS_HWD_Outlier_genepop.gen")$MARKERS
age<- tibble::as_tibble(read.csv("Strata_Sex_Model_COHORT_REGRESSION.tsv"))

input<- tibble::enframe(dart_markers)
rand_gi_list<- list()
for(i in c(1:1000)){
  rand_gi_list[[i]]<-  input %>%
    dplyr::sample_n(., size = 500, replace = F)
}
counter = 1
for( l in rand_gi_list){
  out= dart %>%
    dplyr::filter(MARKERS %in% l$value)
  
  write_genepop(out, filename = paste0("SubSample_Genepop_", counter))
  counter = counter + 1
}
## then run 
~/Desktop/NeEstimator/Ne2-1M c:batch_file_1000.txt










##-------- MSATS
# Prep files for COLONY
# We will use average missing rate 0, error rate 0.02

## Check MATURITY 
## Pratt (1996) studied male white sharks in the WNA and suggested that the smallest mature male in his sample 
##(352 cm FL) should be considered size at maturity. 
## Francis (1996) suggested that female size at maturity occurred over a broad size range (female 450–500 cm TL, 417–464cm FL). 
## HERE **males were considered mature when >3.6 m, and females were considered mature when >4.5 m (Francis, 1996; Pratt, 1996).

age<- as.tibble(read.csv("../../../Age_Samples/Cohort_Info_No_L0_Regression.csv"))
# which are m and bigger then 3.5?
# which are f and bigger then 4.5?
male<- age %>%
  dplyr::filter(Sex == "M") %>%
  dplyr::select(MBB_Code, Sex, Predicted_TL) %>%
  dplyr::filter(Predicted_TL >= 360)
# none
female<- age %>%
  dplyr::filter(Sex == "F") %>%
  dplyr::select(MBB_Code, Sex, Predicted_TL) %>%
  dplyr::filter(Predicted_TL >= 450)
# none 

# check summary ages sent in msat 
age<- as.tibble(read.csv("../../../../Age_Samples/Cohort_Info_No_L0_Regression.csv"))
msat.id<- read.csv("../../../../../Sara_Andreotti/NSWDPI_Sara_Full_Info.csv")$TARGET_ID

ages<- age %>%
  filter(sub(" ", "_", MBB_Code) %in% msat.id) %>%
  group_by(sex_model_cohort) %>%
  tally()