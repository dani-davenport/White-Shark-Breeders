library(ggplot2)
library(descriptr)
library(dplyr)
#---- Subsampled loci 
# fix up output files bash
## sed -i '' '/^[-]/d' output_sublocixLD.txt
## sed -i '' '/genepop/d' output_sublocixLD.txt
# Plot
setwd("~/Google Drive/UQ_PhD/Lab_Book/ECWS_NSW_Fisheries/Data/Data_Filtering_2019/SubSample_Loci")
dat500<- readxl::read_xlsx("Loci_500/Loci_500_OutputSubLociLD.xlsx")
dat1000<- readxl::read_xlsx("Loci_1000/Loci_1000_OutputSubLociLD.xlsx")
dat2000<- readxl::read_xlsx("Loci_2000/Loci_2000_OutputSubLociLD.xlsx")
# get the cohorts we want
ckeep<- dat500$Pop[1:5]
dat5 <- dat500 %>% 
  filter(Pop %in% ckeep) %>% # change the cohort name/factor
  mutate(Pop_Fact = as.factor(Pop)) %>%
  mutate(LCI_Jackknife = as.double(LCI_Jackknife)) %>%
  mutate(UCI_Jackknife = as.double(UCI_Jackknife))
dat5<- dat5 %>% add_column(Loci_Number = rep(500, nrow(dat5)))

dat1 <- dat1000 %>% 
  filter(Pop %in% ckeep) %>% # change the cohort name/factor
  mutate(Pop_Fact = as.factor(Pop)) %>%
  mutate(LCI_Jackknife = as.double(LCI_Jackknife)) %>%
  mutate(UCI_Jackknife = as.double(UCI_Jackknife)) 
dat1<- dat1 %>% add_column(Loci_Number = rep(1000, nrow(dat1)))

dat2 <- dat2000 %>% 
  filter(Pop %in% ckeep) %>% # change the cohort name/factor
  mutate(Pop_Fact = as.factor(Pop)) %>%
  mutate(LCI_Jackknife = as.double(LCI_Jackknife)) %>%
  mutate(UCI_Jackknife = as.double(UCI_Jackknife))
dat2<- dat2 %>% add_column(Loci_Number = rep(2000, nrow(dat2)))

dat<- bind_rows(dat5, dat1, dat2)

levels(dat$Pop_Fact)
pops = c("2013", #1:MBB-1160
         "2014", #2:MBB-1161
         "2012", #3:MBB-1163
         "2010", #4:MBB-1162
         "2011") #5:MBB-1158
levels(dat$Pop_Fact) <- pops
dat$Pop_Fact<- factor(dat$Pop_Fact, levels =c("2010", "2011", "2012", "2013", "2014"))
# Plot denisity of Nb(point) estimates 
p<- dat %>%
  #filter(Pop != "2:MBB-1161") %>%
  ggplot(aes(x=Ne, group = Pop_Fact, color = Pop_Fact)) + 
  stat_density(geom="line", size = 1.5) +
  facet_wrap(~Loci_Number, nrow = 1, strip.position = "bottom") +
  scale_color_manual(values=wes_palette(n=5, name="Darjeeling1")) +  
  coord_cartesian(xlim = c(-500, 1000)) + 
  #geom_vline(xintercept= 0 , color = "black", linetype = 5) +
  xlab("Estimated Nb") + ylab("Density") + 
  theme_bw() + 
  theme(text = element_text(size=35), #legend.position="none"
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) + 
  labs(color = "YOB Cohort")
  
  
# Calculate the H mean of Nb and Jacknife CIs
datx<- dat
datx$UCI_Parametric[dat$UCI_Parametric == "Infinite"]<- NA
hmean<- datx %>%
  filter(Pop != "2:MBB-1161") %>%
  group_by(.dots=c("Pop_Fact", "Loci_Number")) %>%
  summarize(HMeanNe = descriptr::ds_hmean(Ne, na.rm = F), ULCI = descriptr::ds_hmean(LCI_Jackknife), UCI = descriptr::ds_hmean(UCI_Jackknife))
write.table(hmean, "Results_Summary_SubSampled.tsv", sep = "\t", quote = F, row.names = F)
