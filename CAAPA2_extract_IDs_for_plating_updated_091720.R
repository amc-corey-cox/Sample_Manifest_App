library(gmodels)
library(psych)
library(doBy)
library(dplyr)
library(lubridate)
library(stringi)
library(stringr)
library(readxl)

data1<-read.csv("CAAPA2NewApproachesF-AsthmaAgeGender_DATA_2020-05-21_1558.csv",stringsAsFactors = F)
data2 <- read.csv("CAAPA2NewApproachesF-Asthmaagegender_DATA_2020-05-21_1559.csv",stringsAsFactors = F)
data2$recruit_asthma <- gsub(1,2,data2$recruit_asthma)
data2$recruit_asthma <- gsub(0,1,data2$recruit_asthma)
data2$arhq_asthma_confirm <- gsub(1,2,data2$arhq_asthma_confirm)
data2$arhq_asthma_confirm <- gsub(0,1,data2$arhq_asthma_confirm)

data <- rbind(data1, data2)
data$sid<-gsub(",","",data$sid)
data$sid<-as.numeric(data$sid)
data_sub<-data[which(data$sid != "" & data$sid>=10201001 & data$sid <=10207200),] # remove rows with empty sid  
colnames(data_sub)[2] <- "Site"
table(data_sub$Site)

table(data_sub$calculated_age)


#Add age category (make sure for sites Baltimore , washingtonDC, Barbados, Brazil Age category is all Adults and Denver is all Child and Chicago and Nigeria are 001-050 sid are Adults and 051-100 sids are Peds)/
#data_sub$Age<-ifelse(data_sub$arhq_date != "" | data_sub$crhq_date == "","Adult","Child")

data_sub$Age<-data_sub$calculated_age
data_sub$Age<-ifelse(is.na(data_sub$Age) == T , ifelse(data_sub$arhq_date != "", as.period(interval(data_sub$arhq_birthdate,data_sub$arhq_date),unit = "year")$year, as.period(interval(data_sub$crhq_birthdate,data_sub$crhq_date),unit = "year")$year),data_sub$Age)

data_sub$Age_category<-ifelse(data_sub$Age>=18,"Adult","Child")

#Asthma status
data_sub$Asthma_Self <- ifelse(data_sub$Age_category=="Adult",data_sub$arhq_asthma, data_sub$crhq_asthma)
data_sub$Asthma_MD <- ifelse(data_sub$Age_category=="Adult",data_sub$arhq_asthma_confirm, data_sub$crhq_asthma_confirm)
data_sub$Asthma <- ifelse(data_sub$Asthma_Self==1,1,ifelse(data_sub$Asthma_Self==2,ifelse(data_sub$Asthma_MD==2,2,0),"NA"))


#remove duplicate sids 
data_sub<-data_sub[which(!duplicated(data_sub$sid)),] ################## NOTE: Check duplicated sid and keep sid that has more data across all columns #### 

#Gender
data_sub$Gender <- ifelse(is.na(data_sub$recruit_gender == T),ifelse(data_sub$Age_category=="Adult", data_sub$arhq_gender, data_sub$crhq_gender),data_sub$recruit_gender)


### Overall numbers for the field "Enrolled" in the final recruitment table #########
CrossTable(data_sub[data_sub$Site=="03_baltimore",]$Asthma, data_sub[data_sub$Site=="03_baltimore",]$Age_category)
CrossTable(data_sub[data_sub$Site=="01_denver",]$Asthma, data_sub[data_sub$Site=="01_denver",]$Age_category)
CrossTable(data_sub[data_sub$Site=="02_chicago",]$Asthma, data_sub[data_sub$Site=="02_chicago",]$Age_category)
CrossTable(data_sub[data_sub$Site=="07_nigeria",]$Asthma, data_sub[data_sub$Site=="07_nigeria",]$Age_category)
CrossTable(data_sub[data_sub$Site=="04_washingtondc",]$Asthma, data_sub[data_sub$Site=="04_washingtondc",]$Age_category)
CrossTable(data_sub[data_sub$Site=="05_brazil",]$Asthma, data_sub[data_sub$Site=="05_brazil",]$Age_category)

colnames(data_sub)
data_sub.condensed<-data_sub[,c(2,3,19,22,23)]

write.table(data_sub.condensed,file="CAAPA2_phenotypes_site_age_gender_asthma_052120.txt",sep="\t",row.names=F,quote=F)

phenotype <- read.delim("CAAPA2_phenotypes_site_age_gender_asthma_052120.txt", stringsAsFactors = F)

######## Download the master files from google drive #######

########## Read in DC brazil & baltimore file 
pbmc_dna1 <- read_excel("CAAPA2_DC_Brazil_Baltimore_Master_file_072120_edited.xlsx",sheet=1,skip=1)
colnames(pbmc_dna1)<- gsub("...[0-9]+$","",colnames(pbmc_dna1))

pbmc_rna1 <- read_excel("CAAPA2_DC_Brazil_Baltimore_Master_file_072120_edited.xlsx",sheet=2,skip=1)
colnames(pbmc_rna1)<- gsub("...[0-9]+$","",colnames(pbmc_rna1))

nasal_dna1 <- read_excel("CAAPA2_DC_Brazil_Baltimore_Master_file_072120_edited.xlsx",sheet=3,skip=1)
colnames(nasal_dna1)<- gsub("...[0-9]+$","",colnames(nasal_dna1))

nasal_rna1 <- read_excel("CAAPA2_DC_Brazil_Baltimore_Master_file_072120_edited.xlsx",sheet=4,skip=1)
colnames(nasal_rna1)<- gsub("...[0-9]+$","",colnames(nasal_rna1))

######## Read in Denver Chicago & Nigeria file 
pbmc_dna2 <- read_excel("CAAPA2_Denver_Chicago_Nigeria_Master_file_072120_edited.xlsx",sheet=3,skip=1)
colnames(pbmc_dna2)<- gsub("...[0-9]+$","",colnames(pbmc_dna2))

pbmc_rna2 <- read_excel("CAAPA2_Denver_Chicago_Nigeria_Master_file_072120_edited.xlsx",sheet=4,skip=1)
colnames(pbmc_rna2)<- gsub("...[0-9]+$","",colnames(pbmc_rna2))

nasal_dna2 <- read_excel("CAAPA2_Denver_Chicago_Nigeria_Master_file_072120_edited.xlsx",sheet=5,skip=1)
colnames(nasal_dna2)<- gsub("...[0-9]+$","",colnames(nasal_dna2))

nasal_rna2 <- read_excel("CAAPA2_Denver_Chicago_Nigeria_Master_file_072120_edited.xlsx",sheet=6,skip=1)
colnames(nasal_rna2)<- gsub("...[0-9]+$","",colnames(nasal_rna2))

######### Look for columns between each replicate data entry and add Rep2, Rep3 and Rerecruit #####
######### in the column names to distinguish data entries #######

pbmc_dna <- rbind(pbmc_dna1, pbmc_dna2)
colnames(pbmc_dna) <- make.unique(colnames(pbmc_dna))
colnames(pbmc_dna)<- ifelse(str_detect(colnames(pbmc_dna), "\\.1$")==TRUE, paste("Rep2",colnames(pbmc_dna),sep="_"),ifelse(str_detect(colnames(pbmc_dna), "\\.2$")==TRUE,paste("Rep3",colnames(pbmc_dna),sep="_"),ifelse(str_detect(colnames(pbmc_dna), "\\.3$")==TRUE,paste("ReRecruit",colnames(pbmc_dna),sep="_"),ifelse(str_detect(colnames(pbmc_dna), "\\.4$")==TRUE,paste("ReRecruit2",colnames(pbmc_dna),sep="_"),paste("Rep1",colnames(pbmc_dna),sep="_")))))
colnames(pbmc_dna) <- gsub("\\.[0-9]$","",colnames(pbmc_dna))
colnames(pbmc_dna) <- gsub("/","_",colnames(pbmc_dna))
colnames(pbmc_dna) <- gsub("-","_",colnames(pbmc_dna))

pbmc_rna <- rbind(pbmc_rna1, pbmc_rna2)
colnames(pbmc_rna) <- make.unique(colnames(pbmc_rna))
colnames(pbmc_rna)<- ifelse(str_detect(colnames(pbmc_rna), "\\.1$")==TRUE, paste("Rep2",colnames(pbmc_rna),sep="_"),ifelse(str_detect(colnames(pbmc_rna), "\\.2$")==TRUE,paste("Rep3",colnames(pbmc_rna),sep="_"),ifelse(str_detect(colnames(pbmc_rna), "\\.3$")==TRUE,paste("ReRecruit",colnames(pbmc_rna),sep="_"),ifelse(str_detect(colnames(pbmc_rna), "\\.4$")==TRUE,paste("ReRecruit2",colnames(pbmc_rna),sep="_"),paste("Rep1",colnames(pbmc_rna),sep="_")))))
colnames(pbmc_rna) <- gsub("\\.[0-9]$","",colnames(pbmc_rna))
colnames(pbmc_rna) <- gsub("/","_",colnames(pbmc_rna))
colnames(pbmc_rna) <- gsub("-","_",colnames(pbmc_rna))

nasal_dna <- rbind(nasal_dna1, nasal_dna2)
colnames(nasal_dna) <- make.unique(colnames(nasal_dna))
colnames(nasal_dna)<- ifelse(str_detect(colnames(nasal_dna), "\\.1$")==TRUE, paste("Rep2",colnames(nasal_dna),sep="_"),ifelse(str_detect(colnames(nasal_dna), "\\.2$")==TRUE,paste("Rep3",colnames(nasal_dna),sep="_"),ifelse(str_detect(colnames(nasal_dna), "\\.3$")==TRUE,paste("Rep4",colnames(nasal_dna),sep="_"),ifelse(str_detect(colnames(nasal_dna), "\\.4$")==TRUE,paste("ReRecruit",colnames(nasal_dna),sep="_"),ifelse(str_detect(colnames(nasal_dna), "\\.5$")==TRUE,paste("ReRecruit2",colnames(nasal_dna),sep="_"),paste("Rep1",colnames(nasal_dna),sep="_"))))))
colnames(nasal_dna) <- gsub("\\.[0-9]$","",colnames(nasal_dna))
colnames(nasal_dna) <- gsub("/","_",colnames(nasal_dna))
colnames(nasal_dna) <- gsub("-","_",colnames(nasal_dna))

nasal_rna <- rbind(nasal_rna1, nasal_rna2)
colnames(nasal_rna) <- make.unique(colnames(nasal_rna))
colnames(nasal_rna)<- ifelse(str_detect(colnames(nasal_rna), "\\.1$")==TRUE, paste("Rep2",colnames(nasal_rna),sep="_"),ifelse(str_detect(colnames(nasal_rna), "\\.2$")==TRUE,paste("Rep3",colnames(nasal_rna),sep="_"),ifelse(str_detect(colnames(nasal_rna), "\\.3$")==TRUE,paste("Rep4",colnames(nasal_rna),sep="_"),ifelse(str_detect(colnames(nasal_rna), "\\.4$")==TRUE,paste("ReRecruit",colnames(nasal_rna),sep="_"),ifelse(str_detect(colnames(nasal_rna), "\\.5$")==TRUE,paste("ReRecruit2",colnames(nasal_rna),sep="_"),paste("Rep1",colnames(nasal_rna),sep="_"))))))
colnames(nasal_rna) <- gsub("\\.[0-9]$","",colnames(nasal_rna))
colnames(nasal_rna) <- gsub("/","_",colnames(nasal_rna))
colnames(nasal_rna) <- gsub("-","_",colnames(nasal_rna))

######## Nigeria had some cell count problems due to overdilution that were fixed manually by Tonya and labelled as pass #######
nigeria_overdiluted_samples <- c("10207001","10207002","10207003","10207004","10207005","10207006","10207007","10207008","10207009","10207010","10207011","10207012","10207013","10207014","10207015","10207016","10207017","10207018","10207019","10207020","10207022","10207051","10207052","10207053","10207054","10207055","10207056","10207057","10207058","10207059","10207060","10207061","10207063","10207065","10207066")

########### PBMC DNA #############
######### Repeat the same for pbmc RNA
######### Nasal DNA & RNA have the same conditions for thresholds Except that PBMC assessment considers cell counts
####### and nasal assessment considers Slide status which is already assigned by the lab as "pass/fail
##### Rep1, Rep2, Rep3 get the slide status from "Slide" column and Rerecruits have a separate assessment since
##### a slide is re-collected along with a new nasal sample

###### Total amount in nanograms (ng) = (ng_ul * volume (ul)) 
###### For example, 15ng_ul in 50ul = 750ng

###### Add missing volumes for automated extractions #######
pbmc_dna$Rep1_Volume <- ifelse(is.na(pbmc_dna$Rep1_Volume)==TRUE,50,pbmc_dna$Rep1_Volume)
pbmc_dna$Rep2_Volume <- ifelse(is.na(pbmc_dna$Rep2_Volume)==TRUE,50,pbmc_dna$Rep2_Volume)
pbmc_dna$Rep3_Volume <- ifelse(is.na(pbmc_dna$Rep3_Volume)==TRUE,50,pbmc_dna$Rep3_Volume)
pbmc_dna$ReRecruit_Volume <- ifelse(is.na(pbmc_dna$ReRecruit_Volume)==TRUE,50,pbmc_dna$ReRecruit_Volume)
pbmc_dna$ReRecruit2_Volume <- ifelse(is.na(pbmc_dna$ReRecruit2_Volume)==TRUE,50,pbmc_dna$ReRecruit2_Volume)

####### Lab is going to add volume column to the master spreadsheet
####### Once we have that info, calculate Total by multiplying conc (ng_ul) and volume (ul) to get total in ng

pbmc_dna$Rep1_TOTAL <- pbmc_dna$Rep1_ng_ul*pbmc_dna$Rep1_Volume
pbmc_dna$Rep2_TOTAL<- pbmc_dna$Rep2_ng_ul*pbmc_dna$Rep2_Volume 
pbmc_dna$Rep3_TOTAL<- pbmc_dna$Rep3_ng_ul*pbmc_dna$Rep3_Volume
pbmc_dna$ReRecruit_TOTAL<- pbmc_dna$ReRecruit_ng_ul*pbmc_dna$ReRecruit_Volume
pbmc_dna$ReRecruit2_TOTAL<- pbmc_dna$ReRecruit2_ng_ul*pbmc_dna$ReRecruit2_Volume

pbmc_dna$Rep1_TOTAL_Qubit <- pbmc_dna$Rep1_Qubit_ng_ul*pbmc_dna$Rep1_Volume
pbmc_dna$Rep2_Qubit_ng_ul <- as.numeric(pbmc_dna$Rep2_Qubit_ng_ul)
pbmc_dna$Rep2_TOTAL_Qubit<- pbmc_dna$Rep2_Qubit_ng_ul*pbmc_dna$Rep2_Volume 
pbmc_dna$Rep3_TOTAL_Qubit<- pbmc_dna$Rep3_Qubit_ng_ul*pbmc_dna$Rep3_Volume
pbmc_dna$ReRecruit_TOTAL_Qubit<- pbmc_dna$ReRecruit_Qubit_ng_ul*pbmc_dna$ReRecruit_Volume
pbmc_dna$ReRecruit2_TOTAL_Qubit<- pbmc_dna$ReRecruit2_Qubit_ng_ul*pbmc_dna$ReRecruit2_Volume

######## PBMC cell count pass fail status #########

#### Since rep1, rep2 and rep3 are from the same blood sample, cell counts used are also the same
#### For Re-recruit, since the blood sample is new, use cell counts from Re-recruit columns

pbmc_dna$PBMC_status <- ifelse((is.na(pbmc_dna$Rep1_PercentLive_PBMC_1to200)=="FALSE"), ifelse(pbmc_dna$Rep1_PercentLive_PBMC_1to200>=85,"PASS","FAIL"),ifelse((is.na(pbmc_dna$Rep1_TotalCellCount_default_1to1)=="FALSE"),ifelse(pbmc_dna$Rep1_TotalCellCount_default_1to1>=58700000,"PASS","FAIL"),"Missing"))
pbmc_dna$PBMC_rerecruit_status <- ifelse(is.na(pbmc_dna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(pbmc_dna$ReRecruit_PercentLive_PBMC_1to200)=="FALSE"), ifelse(pbmc_dna$ReRecruit_PercentLive_PBMC_1to200>=85,"PASS","FAIL"),ifelse((is.na(pbmc_dna$ReRecruit_TotalCellCount_default_1to1)=="FALSE"), ifelse(pbmc_dna$ReRecruit_TotalCellCount_default_1to1>=58700000,"PASS","FAIL"),"Missing")),NA)
pbmc_dna$PBMC_rerecruit2_status <- ifelse(is.na(pbmc_dna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(pbmc_dna$ReRecruit2_PercentLive_PBMC_1to200)=="FALSE"), ifelse(pbmc_dna$ReRecruit2_PercentLive_PBMC_1to200>=85,"PASS","FAIL"),ifelse((is.na(pbmc_dna$ReRecruit2_TotalCellCount_default_1to1)=="FALSE"),ifelse(pbmc_dna$ReRecruit2_TotalCellCount_default_1to1>=58700000,"PASS","FAIL"),"Missing")),NA)

pbmc_dna$PBMC_status <- ifelse(pbmc_dna$Rep1_subject %in% nigeria_overdiluted_samples,"PASS",pbmc_dna$PBMC_status)

###### Change these thresholds if needed, consult with Monica
###### ADD Qubit & Tape station QC below in the conditions
# Agilent Tape station: Rep1_Agilent_DIN>=6
# Qubit Concentration: Rep1_Qubit_ng_ul
# Calculate total for Qubit concentration (ng_ul * Vol(ul) to get total amount in ng)
# Qubit total should also be greater than 750ng

### Add in Qubit thresholds also to call the subject pass or fail ####

pbmc_dna$Rep1_status_Nanodrop_amt <- ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep1_TOTAL)==FALSE,ifelse(pbmc_dna$Rep1_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))
pbmc_dna$Rep1_status_Nanodrop_qc <- ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep1_260_280)==FALSE,ifelse((pbmc_dna$Rep1_260_280>=1.4 & pbmc_dna$Rep1_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))
pbmc_dna$Rep1_status_Qubit_amt <- ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep1_TOTAL_Qubit)==FALSE, ifelse(pbmc_dna$Rep1_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))
pbmc_dna$Rep1_status_DIN <- ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep1_Agilent_DIN)==FALSE,ifelse(pbmc_dna$Rep1_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))

pbmc_dna$Rep2_status_Nanodrop_amt <- ifelse(is.na(pbmc_dna$Rep2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep2_TOTAL)==FALSE,ifelse(pbmc_dna$Rep2_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$Rep2_status_Nanodrop_qc <- ifelse(is.na(pbmc_dna$Rep2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep2_260_280)==FALSE,ifelse((pbmc_dna$Rep2_260_280>=1.4 & pbmc_dna$Rep2_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$Rep2_status_Qubit_amt <- ifelse(is.na(pbmc_dna$Rep2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep2_TOTAL_Qubit)==FALSE, ifelse(pbmc_dna$Rep2_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$Rep2_status_DIN <- ifelse(is.na(pbmc_dna$Rep2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep2_Agilent_DIN)==FALSE,ifelse(pbmc_dna$Rep2_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_dna$Rep3_status_Nanodrop_amt <- ifelse(is.na(pbmc_dna$Rep3_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep3_TOTAL)==FALSE,ifelse(pbmc_dna$Rep3_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$Rep3_status_Nanodrop_qc <- ifelse(is.na(pbmc_dna$Rep3_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep3_260_280)==FALSE,ifelse((pbmc_dna$Rep3_260_280>=1.4 & pbmc_dna$Rep3_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$Rep3_status_Qubit_amt <- ifelse(is.na(pbmc_dna$Rep3_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep3_TOTAL_Qubit)==FALSE, ifelse(pbmc_dna$Rep3_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$Rep3_status_DIN <- ifelse(is.na(pbmc_dna$Rep3_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_status =="PASS",ifelse(is.na(pbmc_dna$Rep3_Agilent_DIN)==FALSE,ifelse(pbmc_dna$Rep3_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse(pbmc_dna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_dna$ReRecruit_status_Nanodrop_amt <- ifelse(is.na(pbmc_dna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit_TOTAL)==FALSE,ifelse(pbmc_dna$ReRecruit_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_dna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$ReRecruit_status_Nanodrop_qc <- ifelse(is.na(pbmc_dna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit_260_280)==FALSE,ifelse((pbmc_dna$ReRecruit_260_280>=1.4 & pbmc_dna$ReRecruit_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_dna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$ReRecruit_status_Qubit_amt <- ifelse(is.na(pbmc_dna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit_TOTAL_Qubit)==FALSE, ifelse(pbmc_dna$ReRecruit_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_dna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$ReRecruit_status_DIN <- ifelse(is.na(pbmc_dna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit_Agilent_DIN)==FALSE,ifelse(pbmc_dna$ReRecruit_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse(pbmc_dna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_dna$ReRecruit2_status_Nanodrop_amt <- ifelse(is.na(pbmc_dna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit2_TOTAL)==FALSE,ifelse(pbmc_dna$ReRecruit2_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_dna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$ReRecruit2_status_Nanodrop_qc <- ifelse(is.na(pbmc_dna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit2_260_280)==FALSE,ifelse((pbmc_dna$ReRecruit2_260_280>=1.4 & pbmc_dna$ReRecruit2_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_dna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$ReRecruit2_status_Qubit_amt <- ifelse(is.na(pbmc_dna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit2_TOTAL_Qubit)==FALSE, ifelse(pbmc_dna$ReRecruit2_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_dna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_dna$ReRecruit2_status_DIN <- ifelse(is.na(pbmc_dna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_dna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_dna$ReRecruit2_Agilent_DIN)==FALSE,ifelse(pbmc_dna$ReRecruit2_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse(pbmc_dna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_dna$Rep1_Info <- paste(pbmc_dna$Rep1_status_Nanodrop_amt,pbmc_dna$Rep1_status_Nanodrop_qc,pbmc_dna$Rep1_status_Qubit_amt,pbmc_dna$Rep1_status_DIN,sep=";")
pbmc_dna$Rep1_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_dna$Rep1_Info)
pbmc_dna$Rep1_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_dna$Rep1_Info)
pbmc_dna$Rep1_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_dna$Rep1_Info)
pbmc_dna$Rep1_Info <- gsub(("PASS;"),"",pbmc_dna$Rep1_Info)
pbmc_dna$Rep1_Info <- gsub((";PASS"),"",pbmc_dna$Rep1_Info)
pbmc_dna$Rep1_Info <- gsub("NA.*",NA,pbmc_dna$Rep1_Info)

pbmc_dna$Rep2_Info <- paste(pbmc_dna$Rep2_status_Nanodrop_amt,pbmc_dna$Rep2_status_Nanodrop_qc,pbmc_dna$Rep2_status_Qubit_amt,pbmc_dna$Rep2_status_DIN,sep=";")
pbmc_dna$Rep2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_dna$Rep2_Info)
pbmc_dna$Rep2_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_dna$Rep2_Info)
pbmc_dna$Rep2_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_dna$Rep2_Info)
pbmc_dna$Rep2_Info <- gsub(("PASS;"),"",pbmc_dna$Rep2_Info)
pbmc_dna$Rep2_Info <- gsub((";PASS"),"",pbmc_dna$Rep2_Info)
pbmc_dna$Rep2_Info <- gsub("NA.*",NA,pbmc_dna$Rep2_Info)

pbmc_dna$Rep3_Info <- paste(pbmc_dna$Rep3_status_Nanodrop_amt,pbmc_dna$Rep3_status_Nanodrop_qc,pbmc_dna$Rep3_status_Qubit_amt,pbmc_dna$Rep3_status_DIN,sep=";")
pbmc_dna$Rep3_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_dna$Rep3_Info)
pbmc_dna$Rep3_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_dna$Rep3_Info)
pbmc_dna$Rep3_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_dna$Rep3_Info)
pbmc_dna$Rep3_Info <- gsub(("PASS;"),"",pbmc_dna$Rep3_Info)
pbmc_dna$Rep3_Info <- gsub((";PASS"),"",pbmc_dna$Rep3_Info)
pbmc_dna$Rep3_Info <- gsub("NA.*",NA,pbmc_dna$Rep3_Info)

pbmc_dna$ReRecruit_Info <- paste(pbmc_dna$ReRecruit_status_Nanodrop_amt,pbmc_dna$ReRecruit_status_Nanodrop_qc,pbmc_dna$ReRecruit_status_Qubit_amt,pbmc_dna$ReRecruit_status_DIN,sep=";")
pbmc_dna$ReRecruit_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_dna$ReRecruit_Info)
pbmc_dna$ReRecruit_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_dna$ReRecruit_Info)
pbmc_dna$ReRecruit_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_dna$ReRecruit_Info)
pbmc_dna$ReRecruit_Info <- gsub(("PASS;"),"",pbmc_dna$ReRecruit_Info)
pbmc_dna$ReRecruit_Info <- gsub((";PASS"),"",pbmc_dna$ReRecruit_Info)
pbmc_dna$ReRecruit_Info <- gsub("NA.*",NA,pbmc_dna$ReRecruit_Info)

pbmc_dna$ReRecruit2_Info <- paste(pbmc_dna$ReRecruit2_status_Nanodrop_amt,pbmc_dna$ReRecruit2_status_Nanodrop_qc,pbmc_dna$ReRecruit2_status_Qubit_amt,pbmc_dna$ReRecruit2_status_DIN,sep=";")
pbmc_dna$ReRecruit2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_dna$ReRecruit2_Info)
pbmc_dna$ReRecruit2_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_dna$ReRecruit2_Info)
pbmc_dna$ReRecruit2_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_dna$ReRecruit2_Info)
pbmc_dna$ReRecruit2_Info <- gsub(("PASS;"),"",pbmc_dna$ReRecruit2_Info)
pbmc_dna$ReRecruit2_Info <- gsub((";PASS"),"",pbmc_dna$ReRecruit2_Info)
pbmc_dna$ReRecruit2_Info <- gsub("NA.*",NA,pbmc_dna$ReRecruit2_Info)

pbmc_dna$Rep1_status <- ifelse((!(pbmc_dna$Rep1_Info=="PASS")),"Missing_Fail","PASS")
pbmc_dna$Rep2_status <- ifelse((!(pbmc_dna$Rep2_Info=="PASS")),"Missing_Fail","PASS")
pbmc_dna$Rep3_status <- ifelse((!(pbmc_dna$Rep3_Info=="PASS")),"Missing_Fail","PASS")
pbmc_dna$ReRecruit_status <- ifelse((!(pbmc_dna$ReRecruit_Info=="PASS")),"Missing_Fail","PASS")
pbmc_dna$ReRecruit2_status <- ifelse((!(pbmc_dna$ReRecruit2_Info=="PASS")),"Missing_Fail","PASS")

pbmc_dna$sample_to_use <- ifelse(pbmc_dna$Rep1_status=="Missing_Fail",ifelse((is.na(pbmc_dna$Rep2_status)==FALSE & pbmc_dna$Rep2_status=="PASS"), "Rep2",ifelse((is.na(pbmc_dna$Rep3_status)==FALSE & pbmc_dna$Rep3_status=="PASS"),"Rep3",ifelse((is.na(pbmc_dna$ReRecruit_status)==FALSE & pbmc_dna$ReRecruit_status=="PASS"),"ReRecruit",ifelse((is.na(pbmc_dna$ReRecruit2_status)==FALSE & pbmc_dna$ReRecruit2_status=="PASS"),"ReRecruit2",NA)))),ifelse(pbmc_dna$Rep1_status=="PASS","Rep1",NA))

pbmc_dna$Overall_status <- ifelse(pbmc_dna$Rep1_status=="Missing_Fail",ifelse(((is.na(pbmc_dna$Rep2_status)==FALSE & pbmc_dna$Rep2_status=="PASS")| (is.na(pbmc_dna$Rep3_status)==FALSE & pbmc_dna$Rep3_status=="PASS")|(is.na(pbmc_dna$ReRecruit_status)==FALSE & pbmc_dna$ReRecruit_status=="PASS")|(is.na(pbmc_dna$ReRecruit2_status)==FALSE & pbmc_dna$ReRecruit2_status=="PASS")),"PASS","Missing_Fail"),ifelse(pbmc_dna$Rep1_status=="PASS","PASS",NA))

############## PBMC RNA ###########
######### Look for columns between each replicate data entry and add Rep1, Rep2, Rep3 and Rerecruit #####
######### in the column names to distinguish data entries #######

####### CHECK with Monica if we should do the total amount or concentration and volume 
###### Total amount in nanograms (ng) = (ng_ul * volume (ul)) 
###### For example, 15ng_ul in 50ul = 750ng

###### Add missing volumes for automated extractions #######
pbmc_rna$Rep1_Volume <- ifelse(is.na(pbmc_rna$Rep1_Volume)==TRUE,20,pbmc_rna$Rep1_Volume)
pbmc_rna$Rep2_Volume <- ifelse(is.na(pbmc_rna$Rep2_Volume)==TRUE,20,pbmc_rna$Rep2_Volume)
pbmc_rna$Rep3_Volume <- ifelse(is.na(pbmc_rna$Rep3_Volume)==TRUE,20,pbmc_rna$Rep3_Volume)
pbmc_rna$ReRecruit_Volume <- ifelse(is.na(pbmc_rna$ReRecruit_Volume)==TRUE,20,pbmc_rna$ReRecruit_Volume)
pbmc_rna$ReRecruit2_Volume <- ifelse(is.na(pbmc_rna$ReRecruit2_Volume)==TRUE,20,pbmc_rna$ReRecruit2_Volume)

####### Lab is going to add volume column to the master spreadsheet
####### Once we have that info, calculate Total by multiplying conc (ng_ul) and volume (ul) to get total in ng

pbmc_rna$Rep1_TOTAL <- pbmc_rna$Rep1_ng_ul*pbmc_rna$Rep1_Volume
pbmc_rna$Rep2_TOTAL<- pbmc_rna$Rep2_ng_ul*pbmc_rna$Rep2_Volume 
pbmc_rna$Rep3_TOTAL<- pbmc_rna$Rep3_ng_ul*pbmc_rna$Rep3_Volume
pbmc_rna$ReRecruit_TOTAL<- pbmc_rna$ReRecruit_ng_ul*pbmc_rna$ReRecruit_Volume
pbmc_rna$ReRecruit2_TOTAL<- pbmc_rna$ReRecruit2_ng_ul*pbmc_rna$ReRecruit2_Volume

pbmc_rna$Rep1_TOTAL_Qubit <- pbmc_rna$Rep1_Qubit_ng_ul*pbmc_rna$Rep1_Volume
pbmc_rna$Rep2_TOTAL_Qubit<- pbmc_rna$Rep2_Qubit_ng_ul*pbmc_rna$Rep2_Volume 
pbmc_rna$Rep3_TOTAL_Qubit<- pbmc_rna$Rep3_Qubit_ng_ul*pbmc_rna$Rep3_Volume
pbmc_rna$ReRecruit_TOTAL_Qubit<- pbmc_rna$ReRecruit_Qubit_ng_ul*pbmc_rna$ReRecruit_Volume
pbmc_rna$ReRecruit2_TOTAL_Qubit<- pbmc_rna$ReRecruit2_Qubit_ng_ul*pbmc_rna$ReRecruit2_Volume

######## PBMC cell count pass fail status #########

#### Since rep1, rep2 and rep3 are from the same blood sample, cell counts used are also the same
#### For Re-recruit, since the blood sample is new, use cell counts from Re-recruit columns

pbmc_rna$PBMC_status <- ifelse((is.na(pbmc_rna$Rep1_PercentLive_PBMC_1to200)=="FALSE"), ifelse(pbmc_rna$Rep1_PercentLive_PBMC_1to200>=85,"PASS","FAIL"),ifelse((is.na(pbmc_rna$Rep1_TotalCellCount_default_1to1)=="FALSE"),ifelse(pbmc_rna$Rep1_TotalCellCount_default_1to1>=58700000,"PASS","FAIL"),"Missing"))
pbmc_rna$PBMC_rerecruit_status <- ifelse(is.na(pbmc_rna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(pbmc_rna$ReRecruit_PercentLive_PBMC_1to200)=="FALSE"), ifelse(pbmc_rna$ReRecruit_PercentLive_PBMC_1to200>=85,"PASS","FAIL"),ifelse((is.na(pbmc_rna$ReRecruit_TotalCellCount_default_1to1)=="FALSE"), ifelse(pbmc_rna$ReRecruit_TotalCellCount_default_1to1>=58700000,"PASS","FAIL"),"Missing")),NA)
pbmc_rna$PBMC_rerecruit2_status <- ifelse(is.na(pbmc_rna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(pbmc_rna$ReRecruit2_PercentLive_PBMC_1to200)=="FALSE"), ifelse(pbmc_rna$ReRecruit2_PercentLive_PBMC_1to200>=85,"PASS","FAIL"),ifelse((is.na(pbmc_rna$ReRecruit2_TotalCellCount_default_1to1)=="FALSE"),ifelse(pbmc_rna$ReRecruit2_TotalCellCount_default_1to1>=58700000,"PASS","FAIL"),"Missing")),NA)

pbmc_rna$PBMC_status <- ifelse(pbmc_rna$Rep1_subject %in% nigeria_overdiluted_samples,"PASS",pbmc_rna$PBMC_status)

###### Change these thresholds if needed, consult with Monica
###### ADD Qubit & Tape station QC below in the conditions
# Agilent Tape station: Rep1_Agilent_DIN>=6
# Qubit Concentration: Rep1_Qubit_ng_ul
# Calculate total for Qubit concentration (ng_ul * Vol(ul) to get total amount in ng)
# Qubit total should also be greater than 750ng

### Add in Qubit thresholds also to call the subject pass or fail ####

pbmc_rna$Rep1_status_Nanodrop_amt <- ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep1_TOTAL)==FALSE,ifelse(pbmc_rna$Rep1_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))
pbmc_rna$Rep1_status_Nanodrop_qc <- ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep1_260_280)==FALSE,ifelse((pbmc_rna$Rep1_260_280>=1.7 & pbmc_rna$Rep1_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))
pbmc_rna$Rep1_status_Qubit_amt <- ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep1_TOTAL_Qubit)==FALSE, ifelse(pbmc_rna$Rep1_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))
pbmc_rna$Rep1_status_RIN <- ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep1_Agilent_RINe)==FALSE,ifelse(pbmc_rna$Rep1_Agilent_RINe>=6,"PASS","Fail_RIN"),"Missing_RIN"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts"))

pbmc_rna$Rep2_status_Nanodrop_amt <- ifelse(is.na(pbmc_rna$Rep2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep2_TOTAL)==FALSE,ifelse(pbmc_rna$Rep2_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$Rep2_status_Nanodrop_qc <- ifelse(is.na(pbmc_rna$Rep2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep2_260_280)==FALSE,ifelse((pbmc_rna$Rep2_260_280>=1.7 & pbmc_rna$Rep2_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$Rep2_status_Qubit_amt <- ifelse(is.na(pbmc_rna$Rep2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep2_TOTAL_Qubit)==FALSE, ifelse(pbmc_rna$Rep2_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$Rep2_status_RIN <- ifelse(is.na(pbmc_rna$Rep2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep2_Agilent_RINe)==FALSE,ifelse(pbmc_rna$Rep2_Agilent_RINe>=6,"PASS","Fail_RIN"),"Missing_RIN"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_rna$Rep3_status_Nanodrop_amt <- ifelse(is.na(pbmc_rna$Rep3_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep3_TOTAL)==FALSE,ifelse(pbmc_rna$Rep3_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$Rep3_status_Nanodrop_qc <- ifelse(is.na(pbmc_rna$Rep3_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep3_260_280)==FALSE,ifelse((pbmc_rna$Rep3_260_280>=1.7 & pbmc_rna$Rep3_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$Rep3_status_Qubit_amt <- ifelse(is.na(pbmc_rna$Rep3_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep3_TOTAL_Qubit)==FALSE, ifelse(pbmc_rna$Rep3_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$Rep3_status_RIN <- ifelse(is.na(pbmc_rna$Rep3_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_status =="PASS",ifelse(is.na(pbmc_rna$Rep3_Agilent_RINe)==FALSE,ifelse(pbmc_rna$Rep3_Agilent_RINe>=6,"PASS","Fail_RIN"),"Missing_RIN"),ifelse(pbmc_rna$PBMC_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_rna$ReRecruit_status_Nanodrop_amt <- ifelse(is.na(pbmc_rna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit_TOTAL)==FALSE,ifelse(pbmc_rna$ReRecruit_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_rna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$ReRecruit_status_Nanodrop_qc <- ifelse(is.na(pbmc_rna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit_260_280)==FALSE,ifelse((pbmc_rna$ReRecruit_260_280>=1.7 & pbmc_rna$ReRecruit_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_rna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$ReRecruit_status_Qubit_amt <- ifelse(is.na(pbmc_rna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit_TOTAL_Qubit)==FALSE, ifelse(pbmc_rna$ReRecruit_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_rna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$ReRecruit_status_RIN <- ifelse(is.na(pbmc_rna$ReRecruit_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit_Agilent_RINe)==FALSE,ifelse(pbmc_rna$ReRecruit_Agilent_RINe>=6,"PASS","Fail_RIN"),"Missing_RIN"),ifelse(pbmc_rna$PBMC_rerecruit_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_rna$ReRecruit2_status_Nanodrop_amt <- ifelse(is.na(pbmc_rna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit2_TOTAL)==FALSE,ifelse(pbmc_rna$ReRecruit2_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse(pbmc_rna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$ReRecruit2_status_Nanodrop_qc <- ifelse(is.na(pbmc_rna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit2_260_280)==FALSE,ifelse((pbmc_rna$ReRecruit2_260_280>=1.7 & pbmc_rna$ReRecruit2_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse(pbmc_rna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$ReRecruit2_status_Qubit_amt <- ifelse(is.na(pbmc_rna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit2_TOTAL_Qubit)==FALSE, ifelse(pbmc_rna$ReRecruit2_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse(pbmc_rna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)
pbmc_rna$ReRecruit2_status_RIN <- ifelse(is.na(pbmc_rna$ReRecruit2_ng_ul)==FALSE,ifelse(pbmc_rna$PBMC_rerecruit2_status =="PASS",ifelse(is.na(pbmc_rna$ReRecruit2_Agilent_RINe)==FALSE,ifelse(pbmc_rna$ReRecruit2_Agilent_RINe>=6,"PASS","Fail_RIN"),"Missing_RIN"),ifelse(pbmc_rna$PBMC_rerecruit2_status =="FAIL","Fail_Cellcounts","Missing_Cellcounts")),NA)

pbmc_rna$Rep1_Info <- paste(pbmc_rna$Rep1_status_Nanodrop_amt,pbmc_rna$Rep1_status_Nanodrop_qc,pbmc_rna$Rep1_status_Qubit_amt,pbmc_rna$Rep1_status_RIN,sep=";")
pbmc_rna$Rep1_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_rna$Rep1_Info)
pbmc_rna$Rep1_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_rna$Rep1_Info)
pbmc_rna$Rep1_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_rna$Rep1_Info)
pbmc_rna$Rep1_Info <- gsub(("PASS;"),"",pbmc_rna$Rep1_Info)
pbmc_rna$Rep1_Info <- gsub((";PASS"),"",pbmc_rna$Rep1_Info)
pbmc_rna$Rep1_Info <- gsub("NA.*",NA,pbmc_rna$Rep1_Info)

pbmc_rna$Rep2_Info <- paste(pbmc_rna$Rep2_status_Nanodrop_amt,pbmc_rna$Rep2_status_Nanodrop_qc,pbmc_rna$Rep2_status_Qubit_amt,pbmc_rna$Rep2_status_RIN,sep=";")
pbmc_rna$Rep2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_rna$Rep2_Info)
pbmc_rna$Rep2_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_rna$Rep2_Info)
pbmc_rna$Rep2_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_rna$Rep2_Info)
pbmc_rna$Rep2_Info <- gsub(("PASS;"),"",pbmc_rna$Rep2_Info)
pbmc_rna$Rep2_Info <- gsub((";PASS"),"",pbmc_rna$Rep2_Info)
pbmc_rna$Rep2_Info <- gsub("NA.*",NA,pbmc_rna$Rep2_Info)

pbmc_rna$Rep3_Info <- paste(pbmc_rna$Rep3_status_Nanodrop_amt,pbmc_rna$Rep3_status_Nanodrop_qc,pbmc_rna$Rep3_status_Qubit_amt,pbmc_rna$Rep3_status_RIN,sep=";")
pbmc_rna$Rep3_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_rna$Rep3_Info)
pbmc_rna$Rep3_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_rna$Rep3_Info)
pbmc_rna$Rep3_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_rna$Rep3_Info)
pbmc_rna$Rep3_Info <- gsub(("PASS;"),"",pbmc_rna$Rep3_Info)
pbmc_rna$Rep3_Info <- gsub((";PASS"),"",pbmc_rna$Rep3_Info)
pbmc_rna$Rep3_Info <- gsub("NA.*",NA,pbmc_rna$Rep3_Info)

pbmc_rna$ReRecruit_Info <- paste(pbmc_rna$ReRecruit_status_Nanodrop_amt,pbmc_rna$ReRecruit_status_Nanodrop_qc,pbmc_rna$ReRecruit_status_Qubit_amt,pbmc_rna$ReRecruit_status_RIN,sep=";")
pbmc_rna$ReRecruit_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_rna$ReRecruit_Info)
pbmc_rna$ReRecruit_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_rna$ReRecruit_Info)
pbmc_rna$ReRecruit_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_rna$ReRecruit_Info)
pbmc_rna$ReRecruit_Info <- gsub(("PASS;"),"",pbmc_rna$ReRecruit_Info)
pbmc_rna$ReRecruit_Info <- gsub((";PASS"),"",pbmc_rna$ReRecruit_Info)
pbmc_rna$ReRecruit_Info <- gsub("NA.*",NA,pbmc_rna$ReRecruit_Info)

pbmc_rna$ReRecruit2_Info <- paste(pbmc_rna$ReRecruit2_status_Nanodrop_amt,pbmc_rna$ReRecruit2_status_Nanodrop_qc,pbmc_rna$ReRecruit2_status_Qubit_amt,pbmc_rna$ReRecruit2_status_RIN,sep=";")
pbmc_rna$ReRecruit2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",pbmc_rna$ReRecruit2_Info)
pbmc_rna$ReRecruit2_Info <- gsub("Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts;Fail_Cellcounts","Fail_Cellcounts",pbmc_rna$ReRecruit2_Info)
pbmc_rna$ReRecruit2_Info <- gsub("Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts;Missing_Cellcounts","Missing_Cellcounts",pbmc_rna$ReRecruit2_Info)
pbmc_rna$ReRecruit2_Info <- gsub(("PASS;"),"",pbmc_rna$ReRecruit2_Info)
pbmc_rna$ReRecruit2_Info <- gsub((";PASS"),"",pbmc_rna$ReRecruit2_Info)
pbmc_rna$ReRecruit2_Info <- gsub("NA.*",NA,pbmc_rna$ReRecruit2_Info)

pbmc_rna$Rep1_status <- ifelse((!(pbmc_rna$Rep1_Info=="PASS")),"Missing_Fail","PASS")
pbmc_rna$Rep2_status <- ifelse((!(pbmc_rna$Rep2_Info=="PASS")),"Missing_Fail","PASS")
pbmc_rna$Rep3_status <- ifelse((!(pbmc_rna$Rep3_Info=="PASS")),"Missing_Fail","PASS")
pbmc_rna$ReRecruit_status <- ifelse((!(pbmc_rna$ReRecruit_Info=="PASS")),"Missing_Fail","PASS")
pbmc_rna$ReRecruit2_status <- ifelse((!(pbmc_rna$ReRecruit2_Info=="PASS")),"Missing_Fail","PASS")

pbmc_rna$sample_to_use <- ifelse(pbmc_rna$Rep1_status=="Missing_Fail",ifelse((is.na(pbmc_rna$Rep2_status)==FALSE & pbmc_rna$Rep2_status=="PASS"), "Rep2",ifelse((is.na(pbmc_rna$Rep3_status)==FALSE & pbmc_rna$Rep3_status=="PASS"),"Rep3",ifelse((is.na(pbmc_rna$ReRecruit_status)==FALSE & pbmc_rna$ReRecruit_status=="PASS"),"ReRecruit",ifelse((is.na(pbmc_rna$ReRecruit2_status)==FALSE & pbmc_rna$ReRecruit2_status=="PASS"),"ReRecruit2",NA)))),ifelse(pbmc_rna$Rep1_status=="PASS","Rep1",NA))

pbmc_rna$Overall_status <- ifelse(pbmc_rna$Rep1_status=="Missing_Fail",ifelse(((is.na(pbmc_rna$Rep2_status)==FALSE & pbmc_rna$Rep2_status=="PASS")| (is.na(pbmc_rna$Rep3_status)==FALSE & pbmc_rna$Rep3_status=="PASS")|(is.na(pbmc_rna$ReRecruit_status)==FALSE & pbmc_rna$ReRecruit_status=="PASS")|(is.na(pbmc_rna$ReRecruit2_status)==FALSE & pbmc_rna$ReRecruit2_status=="PASS")),"PASS","Missing_Fail"),ifelse(pbmc_rna$Rep1_status=="PASS","PASS",NA))

########### Nasal DNA #############
######### Look for columns between each replicate data entry and add Rep1, Rep2, Rep3 and Rerecruit #####
######### in the column names to distinguish data entries #######

####### CHECK with Monica if we should do the total amount or concentration and volume 
###### Total amount in nanograms (ng) = (ng_ul * volume (ul)) 
###### For example, 15ng_ul in 50ul = 750ng

###### Add missing volumes for automated extractions #######
nasal_dna$Rep1_Volume <- ifelse(is.na(nasal_dna$Rep1_Volume)==TRUE,50,nasal_dna$Rep1_Volume)
nasal_dna$Rep2_Volume <- ifelse(is.na(nasal_dna$Rep2_Volume)==TRUE,50,nasal_dna$Rep2_Volume)
nasal_dna$Rep3_Volume <- ifelse(is.na(nasal_dna$Rep3_Volume)==TRUE,50,nasal_dna$Rep3_Volume)
nasal_dna$Rep4_Volume <- ifelse(is.na(nasal_dna$Rep4_Volume)==TRUE,50,nasal_dna$Rep4_Volume)
nasal_dna$ReRecruit_Volume <- ifelse(is.na(nasal_dna$ReRecruit_Volume)==TRUE,50,nasal_dna$ReRecruit_Volume)
nasal_dna$ReRecruit2_Volume <- ifelse(is.na(nasal_dna$ReRecruit2_Volume)==TRUE,50,nasal_dna$ReRecruit2_Volume)

####### Lab is going to add volume column to the master spreadsheet
####### Once we have that info, calculate Total by multiplying conc (ng_ul) and volume (ul) to get total in ng

nasal_dna$Rep1_TOTAL <- nasal_dna$Rep1_ng_ul*nasal_dna$Rep1_Volume
nasal_dna$Rep2_TOTAL<- nasal_dna$Rep2_ng_ul*nasal_dna$Rep2_Volume 
nasal_dna$Rep3_TOTAL<- nasal_dna$Rep3_ng_ul*nasal_dna$Rep3_Volume
nasal_dna$Rep4_TOTAL<- nasal_dna$Rep4_ng_ul*nasal_dna$Rep4_Volume
nasal_dna$ReRecruit_TOTAL<- nasal_dna$ReRecruit_ng_ul*nasal_dna$ReRecruit_Volume
nasal_dna$ReRecruit2_TOTAL<- nasal_dna$ReRecruit2_ng_ul*nasal_dna$ReRecruit2_Volume

nasal_dna$Rep1_TOTAL_Qubit <- nasal_dna$Rep1_Qubit_ng_ul*nasal_dna$Rep1_Volume
nasal_dna$Rep2_TOTAL_Qubit<- nasal_dna$Rep2_Qubit_ng_ul*nasal_dna$Rep2_Volume
nasal_dna$Rep3_TOTAL_Qubit<- nasal_dna$Rep3_Qubit_ng_ul*nasal_dna$Rep3_Volume
nasal_dna$Rep4_TOTAL_Qubit<- nasal_dna$Rep4_Qubit_ng_ul*nasal_dna$Rep4_Volume
nasal_dna$ReRecruit_TOTAL_Qubit<- nasal_dna$ReRecruit_Qubit_ng_ul*nasal_dna$ReRecruit_Volume
nasal_dna$ReRecruit2_TOTAL_Qubit<- nasal_dna$ReRecruit2_Qubit_ng_ul*nasal_dna$ReRecruit2_Volume


###### Change these thresholds if needed, consult with Monica
###### ADD Qubit & Tape station QC below in the conditions
# Agilent Tape station: Rep1_Agilent_DIN>=6
# Qubit Concentration: Rep1_Qubit_ng_ul
# Calculate total for Qubit concentration (ng_ul * Vol(ul) to get total amount in ng)
# Qubit total should also be greater than 750ng

### Add in Qubit thresholds also to call the subject pass or fail ####
nasal_dna$Rep1_status_Nanodrop_amt <- ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep1_TOTAL)==FALSE,ifelse(nasal_dna$Rep1_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))
nasal_dna$Rep1_status_Nanodrop_qc <- ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep1_260_280)==FALSE,ifelse((nasal_dna$Rep1_260_280>=1.4 & nasal_dna$Rep1_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))
nasal_dna$Rep1_status_Qubit_amt <- ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep1_TOTAL_Qubit)==FALSE, ifelse(nasal_dna$Rep1_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))
nasal_dna$Rep1_status_DIN <- ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep1_Agilent_DIN)==FALSE,ifelse(nasal_dna$Rep1_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse((is.na(nasal_dna$Rep1_Slide)==FALSE & nasal_dna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))

nasal_dna$Rep2_status_Nanodrop_amt <- ifelse(is.na(nasal_dna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep2_TOTAL)==FALSE,ifelse(nasal_dna$Rep2_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep2_status_Nanodrop_qc <- ifelse(is.na(nasal_dna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep2_260_280)==FALSE,ifelse((nasal_dna$Rep2_260_280>=1.4 & nasal_dna$Rep2_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep2_status_Qubit_amt <- ifelse(is.na(nasal_dna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep2_TOTAL_Qubit)==FALSE, ifelse(nasal_dna$Rep2_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep2_status_DIN <- ifelse(is.na(nasal_dna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep2_Agilent_DIN)==FALSE,ifelse(nasal_dna$Rep2_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse((is.na(nasal_dna$Rep2_Slide)==FALSE & nasal_dna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_dna$Rep3_status_Nanodrop_amt <- ifelse(is.na(nasal_dna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep3_TOTAL)==FALSE,ifelse(nasal_dna$Rep3_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep3_status_Nanodrop_qc <- ifelse(is.na(nasal_dna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep3_260_280)==FALSE,ifelse((nasal_dna$Rep3_260_280>=1.4 & nasal_dna$Rep3_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep3_status_Qubit_amt <- ifelse(is.na(nasal_dna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep3_TOTAL_Qubit)==FALSE, ifelse(nasal_dna$Rep3_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep3_status_DIN <- ifelse(is.na(nasal_dna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep3_Agilent_DIN)==FALSE,ifelse(nasal_dna$Rep3_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse((is.na(nasal_dna$Rep3_Slide)==FALSE & nasal_dna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_dna$Rep4_status_Nanodrop_amt <- ifelse(is.na(nasal_dna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep4_TOTAL)==FALSE,ifelse(nasal_dna$Rep4_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep4_status_Nanodrop_qc <- ifelse(is.na(nasal_dna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep4_260_280)==FALSE,ifelse((nasal_dna$Rep4_260_280>=1.4 & nasal_dna$Rep4_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep4_status_Qubit_amt <- ifelse(is.na(nasal_dna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep4_TOTAL_Qubit)==FALSE, ifelse(nasal_dna$Rep4_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$Rep4_status_DIN <- ifelse(is.na(nasal_dna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_dna$Rep4_Agilent_DIN)==FALSE,ifelse(nasal_dna$Rep4_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse((is.na(nasal_dna$Rep4_Slide)==FALSE & nasal_dna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_dna$ReRecruit_status_Nanodrop_amt <- ifelse(is.na(nasal_dna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit_TOTAL)==FALSE,ifelse(nasal_dna$ReRecruit_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$ReRecruit_status_Nanodrop_qc <- ifelse(is.na(nasal_dna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit_260_280)==FALSE,ifelse((nasal_dna$ReRecruit_260_280>=1.4 & nasal_dna$ReRecruit_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$ReRecruit_status_Qubit_amt <- ifelse(is.na(nasal_dna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit_TOTAL_Qubit)==FALSE, ifelse(nasal_dna$ReRecruit_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$ReRecruit_status_DIN <- ifelse(is.na(nasal_dna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit_Agilent_DIN)==FALSE,ifelse(nasal_dna$ReRecruit_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse((is.na(nasal_dna$ReRecruit_Slide)==FALSE & nasal_dna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_dna$ReRecruit2_status_Nanodrop_amt <- ifelse(is.na(nasal_dna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit2_TOTAL)==FALSE,ifelse(nasal_dna$ReRecruit2_TOTAL>=750,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$ReRecruit2_status_Nanodrop_qc <- ifelse(is.na(nasal_dna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit2_260_280)==FALSE,ifelse((nasal_dna$ReRecruit2_260_280>=1.4 & nasal_dna$ReRecruit2_260_280<=2.15),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$ReRecruit2_status_Qubit_amt <- ifelse(is.na(nasal_dna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit2_TOTAL_Qubit)==FALSE, ifelse(nasal_dna$ReRecruit2_TOTAL_Qubit>=750,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_dna$ReRecruit2_status_DIN <- ifelse(is.na(nasal_dna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_dna$ReRecruit2_Agilent_DIN)==FALSE,ifelse(nasal_dna$ReRecruit2_Agilent_DIN>=6,"PASS","Fail_DIN"),"Missing_DIN"),ifelse((is.na(nasal_dna$ReRecruit2_Slide)==FALSE & nasal_dna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_dna$Rep1_Info <- paste(nasal_dna$Rep1_status_Nanodrop_amt,nasal_dna$Rep1_status_Nanodrop_qc,nasal_dna$Rep1_status_Qubit_amt,nasal_dna$Rep1_status_DIN,sep=";")
nasal_dna$Rep1_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_dna$Rep1_Info)
nasal_dna$Rep1_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_dna$Rep1_Info)
nasal_dna$Rep1_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_dna$Rep1_Info)
nasal_dna$Rep1_Info <- gsub(("PASS;"),"",nasal_dna$Rep1_Info)
nasal_dna$Rep1_Info <- gsub((";PASS"),"",nasal_dna$Rep1_Info)
nasal_dna$Rep1_Info <- gsub("NA.*",NA,nasal_dna$Rep1_Info)

nasal_dna$Rep2_Info <- paste(nasal_dna$Rep2_status_Nanodrop_amt,nasal_dna$Rep2_status_Nanodrop_qc,nasal_dna$Rep2_status_Qubit_amt,nasal_dna$Rep2_status_DIN,sep=";")
nasal_dna$Rep2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_dna$Rep2_Info)
nasal_dna$Rep2_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_dna$Rep2_Info)
nasal_dna$Rep2_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_dna$Rep2_Info)
nasal_dna$Rep2_Info <- gsub(("PASS;"),"",nasal_dna$Rep2_Info)
nasal_dna$Rep2_Info <- gsub((";PASS"),"",nasal_dna$Rep2_Info)
nasal_dna$Rep2_Info <- gsub("NA.*",NA,nasal_dna$Rep2_Info)

nasal_dna$Rep3_Info <- paste(nasal_dna$Rep3_status_Nanodrop_amt,nasal_dna$Rep3_status_Nanodrop_qc,nasal_dna$Rep3_status_Qubit_amt,nasal_dna$Rep3_status_DIN,sep=";")
nasal_dna$Rep3_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_dna$Rep3_Info)
nasal_dna$Rep3_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_dna$Rep3_Info)
nasal_dna$Rep3_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_dna$Rep3_Info)
nasal_dna$Rep3_Info <- gsub(("PASS;"),"",nasal_dna$Rep3_Info)
nasal_dna$Rep3_Info <- gsub((";PASS"),"",nasal_dna$Rep3_Info)
nasal_dna$Rep3_Info <- gsub("NA.*",NA,nasal_dna$Rep3_Info)

nasal_dna$Rep4_Info <- paste(nasal_dna$Rep4_status_Nanodrop_amt,nasal_dna$Rep4_status_Nanodrop_qc,nasal_dna$Rep4_status_Qubit_amt,nasal_dna$Rep4_status_DIN,sep=";")
nasal_dna$Rep4_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_dna$Rep4_Info)
nasal_dna$Rep4_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_dna$Rep4_Info)
nasal_dna$Rep4_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_dna$Rep4_Info)
nasal_dna$Rep4_Info <- gsub(("PASS;"),"",nasal_dna$Rep4_Info)
nasal_dna$Rep4_Info <- gsub((";PASS"),"",nasal_dna$Rep4_Info)
nasal_dna$Rep4_Info <- gsub("NA.*",NA,nasal_dna$Rep4_Info)

nasal_dna$ReRecruit_Info <- paste(nasal_dna$ReRecruit_status_Nanodrop_amt,nasal_dna$ReRecruit_status_Nanodrop_qc,nasal_dna$ReRecruit_status_Qubit_amt,nasal_dna$ReRecruit_status_DIN,sep=";")
nasal_dna$ReRecruit_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_dna$ReRecruit_Info)
nasal_dna$ReRecruit_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_dna$ReRecruit_Info)
nasal_dna$ReRecruit_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_dna$ReRecruit_Info)
nasal_dna$ReRecruit_Info <- gsub(("PASS;"),"",nasal_dna$ReRecruit_Info)
nasal_dna$ReRecruit_Info <- gsub((";PASS"),"",nasal_dna$ReRecruit_Info)
nasal_dna$ReRecruit_Info <- gsub("NA.*",NA,nasal_dna$ReRecruit_Info)

nasal_dna$ReRecruit2_Info <- paste(nasal_dna$ReRecruit2_status_Nanodrop_amt,nasal_dna$ReRecruit2_status_Nanodrop_qc,nasal_dna$ReRecruit2_status_Qubit_amt,nasal_dna$ReRecruit2_status_DIN,sep=";")
nasal_dna$ReRecruit2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_dna$ReRecruit2_Info)
nasal_dna$ReRecruit2_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_dna$ReRecruit2_Info)
nasal_dna$ReRecruit2_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_dna$ReRecruit2_Info)
nasal_dna$ReRecruit2_Info <- gsub(("PASS;"),"",nasal_dna$ReRecruit2_Info)
nasal_dna$ReRecruit2_Info <- gsub((";PASS"),"",nasal_dna$ReRecruit2_Info)
nasal_dna$ReRecruit2_Info <- gsub("NA.*",NA,nasal_dna$ReRecruit2_Info)

nasal_dna$Rep1_status <- ifelse((!(nasal_dna$Rep1_Info=="PASS")),"Missing_Fail","PASS")
nasal_dna$Rep2_status <- ifelse((!(nasal_dna$Rep2_Info=="PASS")),"Missing_Fail","PASS")
nasal_dna$Rep3_status <- ifelse((!(nasal_dna$Rep3_Info=="PASS")),"Missing_Fail","PASS")
nasal_dna$Rep4_status <- ifelse((!(nasal_dna$Rep4_Info=="PASS")),"Missing_Fail","PASS")
nasal_dna$ReRecruit_status <- ifelse((!(nasal_dna$ReRecruit_Info=="PASS")),"Missing_Fail","PASS")
nasal_dna$ReRecruit2_status <- ifelse((!(nasal_dna$ReRecruit2_Info=="PASS")),"Missing_Fail","PASS")

nasal_dna$sample_to_use <- ifelse(nasal_dna$Rep1_status=="Missing_Fail",ifelse((is.na(nasal_dna$Rep2_status)==FALSE & nasal_dna$Rep2_status=="PASS"), "Rep2",ifelse((is.na(nasal_dna$Rep3_status)==FALSE & nasal_dna$Rep3_status=="PASS"),"Rep3",ifelse((is.na(nasal_dna$Rep4_status)==FALSE & nasal_dna$Rep4_status=="PASS"),"Rep4", ifelse((is.na(nasal_dna$ReRecruit_status)==FALSE & nasal_dna$ReRecruit_status=="PASS"),"ReRecruit",ifelse((is.na(nasal_dna$ReRecruit2_status)==FALSE & nasal_dna$ReRecruit2_status=="PASS"),"ReRecruit2",NA))))),ifelse(nasal_dna$Rep1_status=="PASS","Rep1",NA))

nasal_dna$Overall_status <- ifelse(nasal_dna$Rep1_status=="Missing_Fail",ifelse(((is.na(nasal_dna$Rep2_status)==FALSE & nasal_dna$Rep2_status=="PASS")| (is.na(nasal_dna$Rep3_status)==FALSE & nasal_dna$Rep3_status=="PASS")|(is.na(nasal_dna$Rep4_status)==FALSE & nasal_dna$Rep4_status=="PASS")|(is.na(nasal_dna$ReRecruit_status)==FALSE & nasal_dna$ReRecruit_status=="PASS")|(is.na(nasal_dna$ReRecruit2_status)==FALSE & nasal_dna$ReRecruit2_status=="PASS")),"PASS","Missing_Fail"),ifelse(nasal_dna$Rep1_status=="PASS","PASS",NA))

#nasal_dna$Overall_status <- ifelse((nasal_dna$Overall_status=="Missing_Fail") & 
########### Nasal RNA #############
######### Look for columns between each replicate data entry and add Rep1, Rep2, Rep3 and Rerecruit #####
######### in the column names to distinguish data entries #######


####### CHECK with Monica if we should do the total amount or concentration and volume 
###### Total amount in nanograms (ng) = (ng_ul * volume (ul)) 
###### For example, 15ng_ul in 50ul = 750ng

###### Add missing volumes for automated extractions #######
nasal_rna$Rep1_Volume <- ifelse(is.na(nasal_rna$Rep1_Volume)==TRUE,20,nasal_rna$Rep1_Volume)
nasal_rna$Rep2_Volume <- ifelse(is.na(nasal_rna$Rep2_Volume)==TRUE,20,nasal_rna$Rep2_Volume)
nasal_rna$Rep3_Volume <- ifelse(is.na(nasal_rna$Rep3_Volume)==TRUE,20,nasal_rna$Rep3_Volume)
nasal_rna$Rep4_Volume <- ifelse(is.na(nasal_rna$Rep4_Volume)==TRUE,20,nasal_rna$Rep4_Volume)
nasal_rna$ReRecruit_Volume <- ifelse(is.na(nasal_rna$ReRecruit_Volume)==TRUE,20,nasal_rna$ReRecruit_Volume)
nasal_rna$ReRecruit2_Volume <- ifelse(is.na(nasal_rna$ReRecruit2_Volume)==TRUE,20,nasal_rna$ReRecruit2_Volume)

####### Lab is going to add Volume column to the master spreadsheet
####### Once we have that info, calculate Total by multiplying conc (ng_ul) and Volume (ul) to get total in ng

nasal_rna$Rep1_TOTAL <- nasal_rna$Rep1_ng_ul*nasal_rna$Rep1_Volume
nasal_rna$Rep2_TOTAL<- nasal_rna$Rep2_ng_ul*nasal_rna$Rep2_Volume 
nasal_rna$Rep3_TOTAL<- nasal_rna$Rep3_ng_ul*nasal_rna$Rep3_Volume
nasal_rna$Rep4_TOTAL<- nasal_rna$Rep4_ng_ul*nasal_rna$Rep4_Volume
nasal_rna$ReRecruit_TOTAL<- nasal_rna$ReRecruit_ng_ul*nasal_rna$ReRecruit_Volume
nasal_rna$ReRecruit2_TOTAL<- nasal_rna$ReRecruit2_ng_ul*nasal_rna$ReRecruit2_Volume

nasal_rna$Rep1_TOTAL_Qubit <- nasal_rna$Rep1_Qubit_ng_ul*nasal_rna$Rep1_Volume
nasal_rna$Rep2_TOTAL_Qubit<- nasal_rna$Rep2_Qubit_ng_ul*nasal_rna$Rep2_Volume 
nasal_rna$Rep3_TOTAL_Qubit<- nasal_rna$Rep3_Qubit_ng_ul*nasal_rna$Rep3_Volume
nasal_rna$Rep4_TOTAL_Qubit<- nasal_rna$Rep4_Qubit_ng_ul*nasal_rna$Rep4_Volume
nasal_rna$ReRecruit_TOTAL_Qubit<- nasal_rna$ReRecruit_Qubit_ng_ul*nasal_rna$ReRecruit_Volume
nasal_rna$ReRecruit2_TOTAL_Qubit<- nasal_rna$ReRecruit2_Qubit_ng_ul*nasal_rna$ReRecruit2_Volume


###### Change these thresholds if needed, consult with Monica
###### ADD Qubit & Tape station QC below in the conditions
# Agilent Tape station: Rep1_Agilent.RIN>=6
# Qubit Concentration: Rep1_Qubit_ng_ul
# Calculate total for Qubit concentration (ng_ul * Vol(ul) to get total amount in ng)
# Qubit total should also be greater than 600ng

### Add in Qubit thresholds also to call the subject pass or fail ####
nasal_rna$Rep1_status_Nanodrop_amt <- ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep1_TOTAL)==FALSE,ifelse(nasal_rna$Rep1_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))
nasal_rna$Rep1_status_Nanodrop_qc <- ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep1_260_280)==FALSE,ifelse((nasal_rna$Rep1_260_280>=1.7 & nasal_rna$Rep1_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))
nasal_rna$Rep1_status_Qubit_amt <- ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep1_TOTAL_Qubit)==FALSE, ifelse(nasal_rna$Rep1_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))
nasal_rna$Rep1_status_RINe <- ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep1_Agilent_RINe)==FALSE,ifelse(nasal_rna$Rep1_Agilent_RINe>=6,"PASS","Fail_RINe"),"Missing_RINe"),ifelse((is.na(nasal_rna$Rep1_Slide)==FALSE & nasal_rna$Rep1_Slide =="FAIL"),"Fail_Slide","Missing_Slide"))

nasal_rna$Rep2_status_Nanodrop_amt <- ifelse(is.na(nasal_rna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep2_TOTAL)==FALSE,ifelse(nasal_rna$Rep2_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep2_status_Nanodrop_qc <- ifelse(is.na(nasal_rna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep2_260_280)==FALSE,ifelse((nasal_rna$Rep2_260_280>=1.7 & nasal_rna$Rep2_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep2_status_Qubit_amt <- ifelse(is.na(nasal_rna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep2_TOTAL_Qubit)==FALSE, ifelse(nasal_rna$Rep2_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep2_status_RINe <- ifelse(is.na(nasal_rna$Rep2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep2_Agilent_RINe)==FALSE,ifelse(nasal_rna$Rep2_Agilent_RINe>=6,"PASS","Fail_RINe"),"Missing_RINe"),ifelse((is.na(nasal_rna$Rep2_Slide)==FALSE & nasal_rna$Rep2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_rna$Rep3_status_Nanodrop_amt <- ifelse(is.na(nasal_rna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep3_TOTAL)==FALSE,ifelse(nasal_rna$Rep3_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep3_status_Nanodrop_qc <- ifelse(is.na(nasal_rna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep3_260_280)==FALSE,ifelse((nasal_rna$Rep3_260_280>=1.7 & nasal_rna$Rep3_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep3_status_Qubit_amt <- ifelse(is.na(nasal_rna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep3_TOTAL_Qubit)==FALSE, ifelse(nasal_rna$Rep3_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep3_status_RINe <- ifelse(is.na(nasal_rna$Rep3_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep3_Agilent_RINe)==FALSE,ifelse(nasal_rna$Rep3_Agilent_RINe>=6,"PASS","Fail_RINe"),"Missing_RINe"),ifelse((is.na(nasal_rna$Rep3_Slide)==FALSE & nasal_rna$Rep3_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_rna$Rep4_status_Nanodrop_amt <- ifelse(is.na(nasal_rna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep4_TOTAL)==FALSE,ifelse(nasal_rna$Rep4_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep4_status_Nanodrop_qc <- ifelse(is.na(nasal_rna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep4_260_280)==FALSE,ifelse((nasal_rna$Rep4_260_280>=1.7 & nasal_rna$Rep4_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep4_status_Qubit_amt <- ifelse(is.na(nasal_rna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep4_TOTAL_Qubit)==FALSE, ifelse(nasal_rna$Rep4_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$Rep4_status_RINe <- ifelse(is.na(nasal_rna$Rep4_ng_ul)==FALSE,ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="PASS"),ifelse(is.na(nasal_rna$Rep4_Agilent_RINe)==FALSE,ifelse(nasal_rna$Rep4_Agilent_RINe>=6,"PASS","Fail_RINe"),"Missing_RINe"),ifelse((is.na(nasal_rna$Rep4_Slide)==FALSE & nasal_rna$Rep4_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_rna$ReRecruit_status_Nanodrop_amt <- ifelse(is.na(nasal_rna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit_TOTAL)==FALSE,ifelse(nasal_rna$ReRecruit_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$ReRecruit_status_Nanodrop_qc <- ifelse(is.na(nasal_rna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit_260_280)==FALSE,ifelse((nasal_rna$ReRecruit_260_280>=1.7 & nasal_rna$ReRecruit_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$ReRecruit_status_Qubit_amt <- ifelse(is.na(nasal_rna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit_TOTAL_Qubit)==FALSE, ifelse(nasal_rna$ReRecruit_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$ReRecruit_status_RINe <- ifelse(is.na(nasal_rna$ReRecruit_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit_Agilent_RINe)==FALSE,ifelse(nasal_rna$ReRecruit_Agilent_RINe>=6,"PASS","Fail_RINe"),"Missing_RINe"),ifelse((is.na(nasal_rna$ReRecruit_Slide)==FALSE & nasal_rna$ReRecruit_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_rna$ReRecruit2_status_Nanodrop_amt <- ifelse(is.na(nasal_rna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit2_TOTAL)==FALSE,ifelse(nasal_rna$ReRecruit2_TOTAL>=600,"PASS","Fail_Nanodrop_Qty"),"Missing_Nanodrop_Qty"),ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$ReRecruit2_status_Nanodrop_qc <- ifelse(is.na(nasal_rna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit2_260_280)==FALSE,ifelse((nasal_rna$ReRecruit2_260_280>=1.7 & nasal_rna$ReRecruit2_260_280<=2.2),"PASS","Fail_260_280"),"Missing_260_280"),ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$ReRecruit2_status_Qubit_amt <- ifelse(is.na(nasal_rna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit2_TOTAL_Qubit)==FALSE, ifelse(nasal_rna$ReRecruit2_TOTAL_Qubit>=600,"PASS","Fail_Qubit_Qty"),"Missing_Qubit"),ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)
nasal_rna$ReRecruit2_status_RINe <- ifelse(is.na(nasal_rna$ReRecruit2_ng_ul)==FALSE,ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="PASS"),ifelse(is.na(nasal_rna$ReRecruit2_Agilent_RINe)==FALSE,ifelse(nasal_rna$ReRecruit2_Agilent_RINe>=6,"PASS","Fail_RINe"),"Missing_RINe"),ifelse((is.na(nasal_rna$ReRecruit2_Slide)==FALSE & nasal_rna$ReRecruit2_Slide =="FAIL"),"Fail_Slide","Missing_Slide")),NA)

nasal_rna$Rep1_Info <- paste(nasal_rna$Rep1_status_Nanodrop_amt,nasal_rna$Rep1_status_Nanodrop_qc,nasal_rna$Rep1_status_Qubit_amt,nasal_rna$Rep1_status_RINe,sep=";")
nasal_rna$Rep1_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_rna$Rep1_Info)
nasal_rna$Rep1_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_rna$Rep1_Info)
nasal_rna$Rep1_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_rna$Rep1_Info)
nasal_rna$Rep1_Info <- gsub(("PASS;"),"",nasal_rna$Rep1_Info)
nasal_rna$Rep1_Info <- gsub((";PASS"),"",nasal_rna$Rep1_Info)
nasal_rna$Rep1_Info <- gsub("NA.*",NA,nasal_rna$Rep1_Info)

nasal_rna$Rep2_Info <- paste(nasal_rna$Rep2_status_Nanodrop_amt,nasal_rna$Rep2_status_Nanodrop_qc,nasal_rna$Rep2_status_Qubit_amt,nasal_rna$Rep2_status_RINe,sep=";")
nasal_rna$Rep2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_rna$Rep2_Info)
nasal_rna$Rep2_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_rna$Rep2_Info)
nasal_rna$Rep2_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_rna$Rep2_Info)
nasal_rna$Rep2_Info <- gsub(("PASS;"),"",nasal_rna$Rep2_Info)
nasal_rna$Rep2_Info <- gsub((";PASS"),"",nasal_rna$Rep2_Info)
nasal_rna$Rep2_Info <- gsub("NA.*",NA,nasal_rna$Rep2_Info)

nasal_rna$Rep3_Info <- paste(nasal_rna$Rep3_status_Nanodrop_amt,nasal_rna$Rep3_status_Nanodrop_qc,nasal_rna$Rep3_status_Qubit_amt,nasal_rna$Rep3_status_RINe,sep=";")
nasal_rna$Rep3_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_rna$Rep3_Info)
nasal_rna$Rep3_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_rna$Rep3_Info)
nasal_rna$Rep3_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_rna$Rep3_Info)
nasal_rna$Rep3_Info <- gsub(("PASS;"),"",nasal_rna$Rep3_Info)
nasal_rna$Rep3_Info <- gsub((";PASS"),"",nasal_rna$Rep3_Info)
nasal_rna$Rep3_Info <- gsub("NA.*",NA,nasal_rna$Rep3_Info)

nasal_rna$Rep4_Info <- paste(nasal_rna$Rep4_status_Nanodrop_amt,nasal_rna$Rep4_status_Nanodrop_qc,nasal_rna$Rep4_status_Qubit_amt,nasal_rna$Rep4_status_RINe,sep=";")
nasal_rna$Rep4_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_rna$Rep4_Info)
nasal_rna$Rep4_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_rna$Rep4_Info)
nasal_rna$Rep4_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_rna$Rep4_Info)
nasal_rna$Rep4_Info <- gsub(("PASS;"),"",nasal_rna$Rep4_Info)
nasal_rna$Rep4_Info <- gsub((";PASS"),"",nasal_rna$Rep4_Info)
nasal_rna$Rep4_Info <- gsub("NA.*",NA,nasal_rna$Rep4_Info)

nasal_rna$ReRecruit_Info <- paste(nasal_rna$ReRecruit_status_Nanodrop_amt,nasal_rna$ReRecruit_status_Nanodrop_qc,nasal_rna$ReRecruit_status_Qubit_amt,nasal_rna$ReRecruit_status_RINe,sep=";")
nasal_rna$ReRecruit_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_rna$ReRecruit_Info)
nasal_rna$ReRecruit_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_rna$ReRecruit_Info)
nasal_rna$ReRecruit_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_rna$ReRecruit_Info)
nasal_rna$ReRecruit_Info <- gsub(("PASS;"),"",nasal_rna$ReRecruit_Info)
nasal_rna$ReRecruit_Info <- gsub((";PASS"),"",nasal_rna$ReRecruit_Info)
nasal_rna$ReRecruit_Info <- gsub("NA.*",NA,nasal_rna$ReRecruit_Info)

nasal_rna$ReRecruit2_Info <- paste(nasal_rna$ReRecruit2_status_Nanodrop_amt,nasal_rna$ReRecruit2_status_Nanodrop_qc,nasal_rna$ReRecruit2_status_Qubit_amt,nasal_rna$ReRecruit2_status_RINe,sep=";")
nasal_rna$ReRecruit2_Info <- gsub("PASS;PASS;PASS;PASS","PASS",nasal_rna$ReRecruit2_Info)
nasal_rna$ReRecruit2_Info <- gsub("Fail_Slide;Fail_Slide;Fail_Slide;Fail_Slide","Fail_Slide",nasal_rna$ReRecruit2_Info)
nasal_rna$ReRecruit2_Info <- gsub("Missing_Slide;Missing_Slide;Missing_Slide;Missing_Slide","Missing_Slide",nasal_rna$ReRecruit2_Info)
nasal_rna$ReRecruit2_Info <- gsub(("PASS;"),"",nasal_rna$ReRecruit2_Info)
nasal_rna$ReRecruit2_Info <- gsub((";PASS"),"",nasal_rna$ReRecruit2_Info)
nasal_rna$ReRecruit2_Info <- gsub("NA.*",NA,nasal_rna$ReRecruit2_Info)

nasal_rna$Rep1_status <- ifelse((!(nasal_rna$Rep1_Info=="PASS")),"Missing_Fail","PASS")
nasal_rna$Rep2_status <- ifelse((!(nasal_rna$Rep2_Info=="PASS")),"Missing_Fail","PASS")
nasal_rna$Rep3_status <- ifelse((!(nasal_rna$Rep3_Info=="PASS")),"Missing_Fail","PASS")
nasal_rna$Rep4_status <- ifelse((!(nasal_rna$Rep4_Info=="PASS")),"Missing_Fail","PASS")
nasal_rna$ReRecruit_status <- ifelse((!(nasal_rna$ReRecruit_Info=="PASS")),"Missing_Fail","PASS")
nasal_rna$ReRecruit2_status <- ifelse((!(nasal_rna$ReRecruit2_Info=="PASS")),"Missing_Fail","PASS")

nasal_rna$sample_to_use <- ifelse(nasal_rna$Rep1_status=="Missing_Fail",ifelse((is.na(nasal_rna$Rep2_status)==FALSE & nasal_rna$Rep2_status=="PASS"), "Rep2",ifelse((is.na(nasal_rna$Rep3_status)==FALSE & nasal_rna$Rep3_status=="PASS"),"Rep3",ifelse((is.na(nasal_rna$Rep4_status)==FALSE & nasal_rna$Rep4_status=="PASS"),"Rep4", ifelse((is.na(nasal_rna$ReRecruit_status)==FALSE & nasal_rna$ReRecruit_status=="PASS"),"ReRecruit",ifelse((is.na(nasal_rna$ReRecruit2_status)==FALSE & nasal_rna$ReRecruit2_status=="PASS"),"ReRecruit2",NA))))),ifelse(nasal_rna$Rep1_status=="PASS","Rep1",NA))

nasal_rna$Overall_status <- ifelse(nasal_rna$Rep1_status=="Missing_Fail",ifelse(((is.na(nasal_rna$Rep2_status)==FALSE & nasal_rna$Rep2_status=="PASS")| (is.na(nasal_rna$Rep3_status)==FALSE & nasal_rna$Rep3_status=="PASS")|(is.na(nasal_rna$Rep4_status)==FALSE & nasal_rna$Rep4_status=="PASS")|(is.na(nasal_rna$ReRecruit_status)==FALSE & nasal_rna$ReRecruit_status=="PASS")|(is.na(nasal_rna$ReRecruit2_status)==FALSE & nasal_rna$ReRecruit2_status=="PASS")),"PASS","Missing_Fail"),ifelse(nasal_rna$Rep1_status=="PASS","PASS",NA))


############### List failed & missing samples for the lab ############

nasal_dna_failed <- nasal_dna[which(nasal_dna$Overall_status=="Missing_Fail"),]
nasal_dna_failed <- nasal_dna_failed[c("Rep1_subject","Rep1_site","Rep1_type","Rep1_Info", "Rep2_Info", "Rep3_Info", "Rep4_Info","ReRecruit_Info", "ReRecruit2_Info")]

nasal_rna_failed <- nasal_rna[which(nasal_rna$Overall_status=="Missing_Fail"),]
nasal_rna_failed <- nasal_rna_failed[c("Rep1_subject","Rep1_site","Rep1_type","Rep1_Info", "Rep2_Info", "Rep3_Info", "Rep4_Info","ReRecruit_Info", "ReRecruit2_Info")]

pbmc_dna_failed <- pbmc_dna[which(pbmc_dna$Overall_status=="Missing_Fail"),]
pbmc_dna_failed <- pbmc_dna_failed[c("Rep1_subject","Rep1_site","Rep1_type","Rep1_Info", "Rep2_Info", "Rep3_Info", "ReRecruit_Info", "ReRecruit2_Info")]

pbmc_rna_failed <- pbmc_rna[which(pbmc_rna$Overall_status=="Missing_Fail"),]
pbmc_rna_failed <- pbmc_rna_failed[c("Rep1_subject","Rep1_site","Rep1_type","Rep1_Info", "Rep2_Info", "Rep3_Info", "ReRecruit_Info", "ReRecruit2_Info")]

write.table(nasal_dna_failed, file="Nasal_DNA_Missing_failed_072320.txt",sep="\t",row.names=F,quote=F)
write.table(nasal_rna_failed, file="Nasal_RNA_Missing_failed_072320.txt",sep="\t",row.names=F,quote=F)
write.table(pbmc_dna_failed, file="PBMC_DNA_Missing_failed_072320.txt",sep="\t",row.names=F,quote=F)
write.table(pbmc_rna_failed, file="PBMC_RNA_Missing_failed_072320.txt",sep="\t",row.names=F,quote=F)

############# All data for plots ##########
nasal_dna$Container_ID <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_container_id,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_container_id,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_container_id,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_container_id,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_container_id,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_container_id,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_container_id,NA)))))))
nasal_dna$Conc_Nanodrop <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_ng_ul,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_ng_ul,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_ng_ul,NA)))))))
nasal_dna$Volume <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_Volume,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_Volume,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_Volume,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_Volume,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_Volume,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_Volume,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_Volume,NA)))))))
nasal_dna$Conc_Qubit <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_Qubit_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_Qubit_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_Qubit_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_Qubit_ng_ul,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_Qubit_ng_ul,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_Qubit_ng_ul,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_Qubit_ng_ul,NA)))))))
nasal_dna$TOTAL_Nanodrop <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_TOTAL,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_TOTAL,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_TOTAL,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_TOTAL,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_TOTAL,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_TOTAL,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_TOTAL,NA)))))))
nasal_dna$TOTAL_Qubit <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_TOTAL_Qubit,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_TOTAL_Qubit,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_TOTAL_Qubit,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_TOTAL_Qubit,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_TOTAL_Qubit,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_TOTAL_Qubit,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_TOTAL_Qubit,NA)))))))
nasal_dna$Nanodrop_260_280 <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_260_280,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_260_280,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_260_280,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_260_280,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_260_280,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_260_280,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_260_280,NA)))))))
nasal_dna$Agilent_DIN <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_Agilent_DIN,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_Agilent_DIN,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_Agilent_DIN,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_Agilent_DIN,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_Agilent_DIN,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_Agilent_DIN,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_Agilent_DIN,NA)))))))
nasal_dna$Slide_status <- ifelse(is.na(nasal_dna$sample_to_use)==TRUE,nasal_dna$Rep1_Slide,ifelse(nasal_dna$sample_to_use=="Rep1",nasal_dna$Rep1_Slide,ifelse(nasal_dna$sample_to_use=="Rep2",nasal_dna$Rep2_Slide,ifelse(nasal_dna$sample_to_use=="Rep3",nasal_dna$Rep3_Slide,ifelse(nasal_dna$sample_to_use=="Rep4",nasal_dna$Rep4_Slide,ifelse(nasal_dna$sample_to_use=="ReRecruit",nasal_dna$ReRecruit_Slide,ifelse(nasal_dna$sample_to_use=="ReRecruit2",nasal_dna$ReRecruit2_Slide,NA)))))))

nasal_rna$Container_ID <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_container_id,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_container_id,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_container_id,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_container_id,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_container_id,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_container_id,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_container_id,NA)))))))
nasal_rna$Conc_Nanodrop <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_ng_ul,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_ng_ul,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_ng_ul,NA)))))))
nasal_rna$Volume <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_Volume,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_Volume,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_Volume,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_Volume,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_Volume,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_Volume,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_Volume,NA)))))))
nasal_rna$Conc_Qubit <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_Qubit_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_Qubit_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_Qubit_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_Qubit_ng_ul,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_Qubit_ng_ul,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_Qubit_ng_ul,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_Qubit_ng_ul,NA)))))))
nasal_rna$TOTAL_Nanodrop <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_TOTAL,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_TOTAL,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_TOTAL,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_TOTAL,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_TOTAL,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_TOTAL,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_TOTAL,NA)))))))
nasal_rna$TOTAL_Qubit <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_TOTAL_Qubit,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_TOTAL_Qubit,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_TOTAL_Qubit,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_TOTAL_Qubit,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_TOTAL_Qubit,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_TOTAL_Qubit,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_TOTAL_Qubit,NA)))))))
nasal_rna$Nanodrop_260_280 <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_260_280,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_260_280,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_260_280,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_260_280,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_260_280,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_260_280,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_260_280,NA)))))))
nasal_rna$Agilent_RINe <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_Agilent_RINe,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_Agilent_RINe,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_Agilent_RINe,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_Agilent_RINe,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_Agilent_RINe,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_Agilent_RINe,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_Agilent_RINe,NA)))))))
nasal_rna$Slide_status <- ifelse(is.na(nasal_rna$sample_to_use)==TRUE,nasal_rna$Rep1_Slide,ifelse(nasal_rna$sample_to_use=="Rep1",nasal_rna$Rep1_Slide,ifelse(nasal_rna$sample_to_use=="Rep2",nasal_rna$Rep2_Slide,ifelse(nasal_rna$sample_to_use=="Rep3",nasal_rna$Rep3_Slide,ifelse(nasal_rna$sample_to_use=="Rep4",nasal_rna$Rep4_Slide,ifelse(nasal_rna$sample_to_use=="ReRecruit",nasal_rna$ReRecruit_Slide,ifelse(nasal_rna$sample_to_use=="ReRecruit2",nasal_rna$ReRecruit2_Slide,NA)))))))

pbmc_dna$Container_ID <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_container_id,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_container_id,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_container_id,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_container_id,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_container_id,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_container_id,NA))))))
pbmc_dna$Conc_Nanodrop <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_ng_ul,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_ng_ul,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_ng_ul,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_ng_ul,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_ng_ul,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_ng_ul,NA))))))
pbmc_dna$Volume <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_Volume,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_Volume,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_Volume,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_Volume,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_Volume,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_Volume,NA))))))
pbmc_dna$Conc_Qubit <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_Qubit_ng_ul,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_Qubit_ng_ul,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_Qubit_ng_ul,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_Qubit_ng_ul,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_Qubit_ng_ul,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_Qubit_ng_ul,NA))))))
pbmc_dna$TOTAL_Nanodrop <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_TOTAL,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_TOTAL,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_TOTAL,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_TOTAL,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_TOTAL,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_TOTAL,NA))))))
pbmc_dna$TOTAL_Qubit <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_TOTAL_Qubit,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_TOTAL_Qubit,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_TOTAL_Qubit,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_TOTAL_Qubit,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_TOTAL_Qubit,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_TOTAL_Qubit,NA))))))
pbmc_dna$Nanodrop_260_280 <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_260_280,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_260_280,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_260_280,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_260_280,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_260_280,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_260_280,NA))))))
pbmc_dna$Agilent_DIN <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$Rep1_Agilent_DIN,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$Rep1_Agilent_DIN,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$Rep2_Agilent_DIN,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$Rep3_Agilent_DIN,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$ReRecruit_Agilent_DIN,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$ReRecruit2_Agilent_DIN,NA))))))
pbmc_dna$PBMC_cellcounts <- ifelse(is.na(pbmc_dna$sample_to_use)==TRUE,pbmc_dna$PBMC_status,ifelse(pbmc_dna$sample_to_use=="Rep1",pbmc_dna$PBMC_status,ifelse(pbmc_dna$sample_to_use=="Rep2",pbmc_dna$PBMC_status,ifelse(pbmc_dna$sample_to_use=="Rep3",pbmc_dna$PBMC_status,ifelse(pbmc_dna$sample_to_use=="ReRecruit",pbmc_dna$PBMC_rerecruit_status,ifelse(pbmc_dna$sample_to_use=="ReRecruit2",pbmc_dna$PBMC_rerecruit2_status,NA))))))

pbmc_rna$Container_ID <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_container_id,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_container_id,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_container_id,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_container_id,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_container_id,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_container_id,NA))))))
pbmc_rna$Conc_Nanodrop <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_ng_ul,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_ng_ul,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_ng_ul,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_ng_ul,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_ng_ul,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_ng_ul,NA))))))
pbmc_rna$Volume <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_Volume,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_Volume,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_Volume,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_Volume,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_Volume,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_Volume,NA))))))
pbmc_rna$Conc_Qubit <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_Qubit_ng_ul,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_Qubit_ng_ul,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_Qubit_ng_ul,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_Qubit_ng_ul,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_Qubit_ng_ul,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_Qubit_ng_ul,NA))))))
pbmc_rna$TOTAL_Nanodrop <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_TOTAL,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_TOTAL,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_TOTAL,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_TOTAL,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_TOTAL,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_TOTAL,NA))))))
pbmc_rna$TOTAL_Qubit <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_TOTAL_Qubit,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_TOTAL_Qubit,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_TOTAL_Qubit,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_TOTAL_Qubit,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_TOTAL_Qubit,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_TOTAL_Qubit,NA))))))
pbmc_rna$Nanodrop_260_280 <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_260_280,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_260_280,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_260_280,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_260_280,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_260_280,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_260_280,NA))))))
pbmc_rna$Agilent_RINe <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$Rep1_Agilent_RINe,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$Rep1_Agilent_RINe,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$Rep2_Agilent_RINe,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$Rep3_Agilent_RINe,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$ReRecruit_Agilent_RINe,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$ReRecruit2_Agilent_RINe,NA))))))
pbmc_rna$PBMC_cellcounts <- ifelse(is.na(pbmc_rna$sample_to_use)==TRUE,pbmc_rna$PBMC_status,ifelse(pbmc_rna$sample_to_use=="Rep1",pbmc_rna$PBMC_status,ifelse(pbmc_rna$sample_to_use=="Rep2",pbmc_rna$PBMC_status,ifelse(pbmc_rna$sample_to_use=="Rep3",pbmc_rna$PBMC_status,ifelse(pbmc_rna$sample_to_use=="ReRecruit",pbmc_rna$PBMC_rerecruit_status,ifelse(pbmc_rna$sample_to_use=="ReRecruit2",pbmc_rna$PBMC_rerecruit2_status,NA))))))

pbmc_dna_pheno <- merge(pbmc_dna, phenotype, by.x="Rep1_subject",by.y="sid")
pbmc_rna_pheno <- merge(pbmc_rna, phenotype, by.x="Rep1_subject",by.y="sid")
nasal_dna_pheno <- merge(nasal_dna, phenotype, by.x="Rep1_subject",by.y="sid")
nasal_rna_pheno <- merge(nasal_rna, phenotype, by.x="Rep1_subject",by.y="sid")

pbmc_dna_for_plots <- pbmc_dna_pheno[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop","Conc_Qubit","TOTAL_Qubit","Nanodrop_260_280","Agilent_DIN","PBMC_cellcounts")]
pbmc_rna_for_plots <- pbmc_rna_pheno[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop","Conc_Qubit","TOTAL_Qubit","Nanodrop_260_280","Agilent_RINe","PBMC_cellcounts")]

nasal_dna_for_plots <- nasal_dna_pheno[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop","Conc_Qubit","TOTAL_Qubit","Nanodrop_260_280","Agilent_DIN","Slide_status")]
nasal_rna_for_plots <- nasal_rna_pheno[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop","Conc_Qubit","TOTAL_Qubit","Nanodrop_260_280","Agilent_RINe","Slide_status")]

pbmc_dna_for_plots$TOTAL_Nanodrop_ug <- (pbmc_dna_for_plots$TOTAL_Nanodrop)/1000
pbmc_dna_for_plots$TOTAL_Qubit_ug <- (pbmc_dna_for_plots$TOTAL_Qubit)/1000
pbmc_dna_for_plots$TOTAL_Nanodrop_ug <- ifelse(is.na(pbmc_dna_for_plots$TOTAL_Nanodrop_ug)==TRUE,"-5",pbmc_dna_for_plots$TOTAL_Nanodrop_ug)
pbmc_dna_for_plots$TOTAL_Qubit_ug <- ifelse(is.na(pbmc_dna_for_plots$TOTAL_Qubit_ug)==TRUE,"-5",pbmc_dna_for_plots$TOTAL_Qubit_ug)
pbmc_dna_for_plots$Nanodrop_260_280 <- ifelse(is.na(pbmc_dna_for_plots$Nanodrop_260_280)==TRUE,"0",pbmc_dna_for_plots$Nanodrop_260_280)
pbmc_dna_for_plots$Agilent_DIN <- ifelse(is.na(pbmc_dna_for_plots$Agilent_DIN)==TRUE,"-1",pbmc_dna_for_plots$Agilent_DIN)
pbmc_dna_final <- pbmc_dna_for_plots[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop_ug","Conc_Qubit","TOTAL_Qubit_ug","Nanodrop_260_280","Agilent_DIN","PBMC_cellcounts")]
colnames(pbmc_dna_final)[1] <- "subject"
colnames(pbmc_dna_final)[2] <- "site"
colnames(pbmc_dna_final)[3] <- "type"
write.table(pbmc_dna_final,file="PBMC_DNA_Final_for_plots_091720.csv",sep=",",row.names=F,quote=F)

pbmc_rna_for_plots$TOTAL_Nanodrop_ug <- (pbmc_rna_for_plots$TOTAL_Nanodrop)/1000
pbmc_rna_for_plots$TOTAL_Qubit_ug <- (pbmc_rna_for_plots$TOTAL_Qubit)/1000
pbmc_rna_for_plots$TOTAL_Nanodrop_ug <- ifelse(is.na(pbmc_rna_for_plots$TOTAL_Nanodrop_ug)==TRUE,"-5",pbmc_rna_for_plots$TOTAL_Nanodrop_ug)
pbmc_rna_for_plots$TOTAL_Qubit_ug <- ifelse(is.na(pbmc_rna_for_plots$TOTAL_Qubit_ug)==TRUE,"-5",pbmc_rna_for_plots$TOTAL_Qubit_ug)
pbmc_rna_for_plots$Nanodrop_260_280 <- ifelse(is.na(pbmc_rna_for_plots$Nanodrop_260_280)==TRUE,"0",pbmc_rna_for_plots$Nanodrop_260_280)
pbmc_rna_for_plots$Agilent_RINe <- ifelse(is.na(pbmc_rna_for_plots$Agilent_RINe)==TRUE,"-1",pbmc_rna_for_plots$Agilent_RINe)
pbmc_rna_final <- pbmc_rna_for_plots[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop_ug","Conc_Qubit","TOTAL_Qubit_ug","Nanodrop_260_280","Agilent_RINe","PBMC_cellcounts")]
colnames(pbmc_rna_final)[1] <- "subject"
colnames(pbmc_rna_final)[2] <- "site"
colnames(pbmc_rna_final)[3] <- "type"
write.table(pbmc_rna_final,file="PBMC_RNA_Final_for_plots_091720.csv",sep=",",row.names=F,quote=F)

nasal_dna_for_plots$TOTAL_Nanodrop_ug <- (nasal_dna_for_plots$TOTAL_Nanodrop)/1000
nasal_dna_for_plots$TOTAL_Qubit_ug <- (nasal_dna_for_plots$TOTAL_Qubit)/1000
nasal_dna_for_plots$TOTAL_Nanodrop_ug <- ifelse(is.na(nasal_dna_for_plots$TOTAL_Nanodrop_ug)==TRUE,"-5",nasal_dna_for_plots$TOTAL_Nanodrop_ug)
nasal_dna_for_plots$TOTAL_Qubit_ug <- ifelse(is.na(nasal_dna_for_plots$TOTAL_Qubit_ug)==TRUE,"-5",nasal_dna_for_plots$TOTAL_Qubit_ug)
nasal_dna_for_plots$Nanodrop_260_280 <- ifelse(is.na(nasal_dna_for_plots$Nanodrop_260_280)==TRUE,"0",nasal_dna_for_plots$Nanodrop_260_280)
nasal_dna_for_plots$Agilent_DIN <- ifelse(is.na(nasal_dna_for_plots$Agilent_DIN)==TRUE,"-1",nasal_dna_for_plots$Agilent_DIN)
nasal_dna_final <- nasal_dna_for_plots[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop_ug","Conc_Qubit","TOTAL_Qubit_ug","Nanodrop_260_280","Agilent_DIN","Slide_status")]
colnames(nasal_dna_final)[1] <- "subject"
colnames(nasal_dna_final)[2] <- "site"
colnames(nasal_dna_final)[3] <- "type"
write.table(nasal_dna_final,file="Nasal_DNA_Final_for_plots_091720.csv",sep=",",row.names=F,quote=F)

nasal_rna_for_plots$TOTAL_Nanodrop_ug <- (nasal_rna_for_plots$TOTAL_Nanodrop)/1000
nasal_rna_for_plots$TOTAL_Qubit_ug <- (nasal_rna_for_plots$TOTAL_Qubit)/1000
nasal_rna_for_plots$TOTAL_Nanodrop_ug <- ifelse(is.na(nasal_rna_for_plots$TOTAL_Nanodrop_ug)==TRUE,"-5",nasal_rna_for_plots$TOTAL_Nanodrop_ug)
nasal_rna_for_plots$TOTAL_Qubit_ug <- ifelse(is.na(nasal_rna_for_plots$TOTAL_Qubit_ug)==TRUE,"-5",nasal_rna_for_plots$TOTAL_Qubit_ug)
nasal_rna_for_plots$Nanodrop_260_280 <- ifelse(is.na(nasal_rna_for_plots$Nanodrop_260_280)==TRUE,"0",nasal_rna_for_plots$Nanodrop_260_280)
nasal_rna_for_plots$Agilent_RINe <- ifelse(is.na(nasal_rna_for_plots$Agilent_RINe)==TRUE,"-1",nasal_rna_for_plots$Agilent_RINe)
nasal_rna_final <- nasal_rna_for_plots[c("Rep1_subject","Rep1_site","Rep1_type","Age_category","Asthma","Gender","Container_ID","Volume","Conc_Nanodrop","TOTAL_Nanodrop_ug","Conc_Qubit","TOTAL_Qubit_ug","Nanodrop_260_280","Agilent_RINe","Slide_status")]
colnames(nasal_rna_final)[1] <- "subject"
colnames(nasal_rna_final)[2] <- "site"
colnames(nasal_rna_final)[3] <- "type"
write.table(nasal_rna_final,file="Nasal_RNA_Final_for_plots_091720.csv",sep=",",row.names=F,quote=F)

