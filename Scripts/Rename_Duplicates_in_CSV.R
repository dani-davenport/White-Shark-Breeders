library(stringr)
#paul<- read.csv("../../Data/Paul's master sheet- UQ genetics PhD070218.csv")
order<- read.csv(file = "../Master_sheet_for_MFL_database.csv", header = T)

#doesnt matter
# label has spaces
newcode<- str_replace_all(order$, "_", "-")
#newcode<- str_replace_all(newcode, " ", "")
#Shark.ID.code<- newcode

# label has - but also random spaces
#newcodepaul<- str_replace_all(paul$Shark.ID.code, " ", "")
#newcodepaul2<- str_replace_all(newcodepaul, "_", "-")
#paul$Shark.ID.code<- newcodepaul2
#neworder<- as.data.frame(cbind(MFL.Reg.Num = as.character(order$MFL.Registration.number), Shark.ID.code))
#pauldata<- left_join(paul, neworder[1:247,], by = "Shark.ID.code" )
# check not missing any
#'%!in%' <- function(x,y)!('%in%'(x,y))
#which(order$Shark.ID.code %!in% Shark.ID.code)
#which(paul$Shark.ID.code %!in% Shark.ID.code)


#282 samples to dart
dartsort<- read.csv("../../DArT_Order/DArT-extract-White_shark-t-2018-03-15.csv", header = T, as.is = T) # this file contains the samples added to the dart plate
colnames(dartsort)[6]<- "MFL.Registration.number"
# join the two datasets
order1<- left_join(dartsort, order, by = "MFL.Registration.number")
# there duplicates in the data file that need to be renamed!
which(duplicated(dartsort$MFL.Registration.number))
# [1]  94 188 282
newnames<- lapply(order1$MFL.Registration.number[c(94, 188, 282)], function(x) paste(x, "Dup", sep = "_"))
order1$MFL.Registration.number[c(94)]<- newnames[[1]]
order1$MFL.Registration.number[c(188)]<- newnames[[2]]
order1$MFL.Registration.number[c(282)]<- newnames[[3]]
order1$MFL.Registration.number<- str_replace_all(order1$MFL.Registration.number, " ", "_")
#$Location<- tidyr::separate(order1, Location, sep = ",", into = c("Location_Beach", "Location_Town") ) 
#order1$Location<- str_replace_all(order1$Location, " ", "_")

# change in datafile as well 
dartfile<- read.csv("../Report-DSha18-3402/Report_DSha18-3402_SNP_singlerow_2.csv", header = F, as.is = T)
coldets<- unlist(dartfile[7,], use.names=FALSE)
coldets2<- str_replace_all(coldets, " ", "_")
which(duplicated(unlist(dartfile[7,], use.names=FALSE)))
#[1]  85 177 287
dartfile[7,c(85, 177, 287)]
# MBB 1341 MBB 1483 MBB 1574

newnames<- lapply(coldets2[c(85,177,287)], function(x) paste(x, "Dup", sep = "_"))
coldets2[85]<- newnames[[1]]
coldets2[177]<- newnames[[2]]
coldets2[287]<- newnames[[3]]
# save new dart file
dartfile[7,]<- coldets2
write.csv(dartfile, quote = F, file = "../Report-DSha18-3402/Report_DSha18_Dups_Relabeled.csv")
# three samples missing
length(coldets2[22:length(coldets2)])
# 278
length(order1$MFL.Registration.number)
# 282
which(order1$MFL.Registration.number %!in% coldets2[22:length(coldets2)])
# [1] 131 177 195 198
order1$MFL.Registration.number[c(131, 177, 195, 198)]
#[1] "MBB_1430" "MBB_1481" "MBB_1498" "MBB_1501"
#### MANNUALLY REMOVED FROM FILE
# make a strata file
makestrata<- as.data.frame(cbind(order1$MFL.Registration.number, order1$MFL.Registration.number, as.character(order1$Location), as.character(order1$State.of.Australia),as.character(order1$Sex), as.character(order1$Total.Length..mm.)))
colnames(makestrata)<- c("TARGET_ID", "INDIVIDUALS", "STRATA_LOC", "STATE", "SEX", "TL")
makestrata$SEX[makestrata$SEX == "Female"]<- "F"
makestrata$SEX[makestrata$SEX == "Male"]<- "M"
makestrata$SEX[makestrata$SEX == "Unknown"]<- "NA"

# size data total length is not corect for some in mm some in cm 

# the duplicate are missing the strata data... need to add
# MBB 1341 MBB 1483 MBB 1574
dets1314<- makestrata[makestrata$TARGET_ID == "MBB_1341",][3:6]
# add
makestrata[makestrata$TARGET_ID == "MBB_1341_Dup",][3:6]<- dets1314

dets1483<- makestrata[makestrata$TARGET_ID == "MBB_1483",][3:6]
# add
makestrata[makestrata$TARGET_ID == "MBB_1483_Dup",][3:6]<- dets1483

dets1574<- makestrata[makestrata$TARGET_ID == "MBB_1574",][3:6]
# add
makestrata[makestrata$TARGET_ID == "MBB_1574_Dup",][3:6]<- dets1574



write.table(makestrata, sep = "\t", quote = F, row.names = F, file = "Full_Strata.tsv")
write.table(makestrata[,c(1,2,4)], sep = "\t", quote = F, row.names = F, file = "ByState_Strata.tsv")
