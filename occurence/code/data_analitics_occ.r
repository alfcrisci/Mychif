#############################################################################################################################
# Setup working directory

#setwd("/home/alf/Scrivania/lav_michyf/occurence")
setwd("..")

##############################################################################################################################
# Load libraries & functions . Check if all packages was installed.

library(XLConnect)
library(doBy)

source("load_lib.r")
source("aux_mycosources.r")

##############################################################################################################################
# Read data

DB_DATA_TRUE=readRDS("data/DB_DATA_TRUE.rds")

DB_DATA_TRUE_occ=DB_DATA_TRUE[which(DB_DATA_TRUE$Co_occurrence==1),]
DB_DATA_TRUE_occ$refunique=gsub("_(.+)","",DB_DATA_TRUE_occ$Ref)
DB_DATA_TRUE_occ$meanTot=as.numeric(DB_DATA_TRUE_occ$meanTot)
DB_DATA_TRUE_occ$Concentration=as.numeric(DB_DATA_TRUE_occ$Concentration)

DB_DATA_TRUE_occ_v=DB_DATA_TRUE_occ[unique(c(which(DB_DATA_TRUE_occ$meanTot>-1),which(DB_DATA_TRUE_occ$Concentration>-1))),]
DB_DATA_TRUE_occ_v$data=ifelse(is.na(DB_DATA_TRUE_occ_v$meanTot),DB_DATA_TRUE_occ_v$Concentration,DB_DATA_TRUE_occ_v$meanTot)

file.remove("PlantTox_occurence_data.xls")

##########################################################################################################################################

XLConnect::writeWorksheetToFile("PlantTox_occurence_data.xls",DB_DATA_TRUE_occ_v,"data_clean")

file.remove("PlantTox_occurence_stat.xls")
full_plants=as.data.frame.matrix(xtabs(~Ref+paramType,data=DB_DATA_TRUE_occ_v))
full_plants=data.frame(ref=row.names(full_plants),full_plants)
XLConnect::writeWorksheetToFile("PlantTox_occurence_stat.xls",full_plants,"full_plants")

plants=unique(DB_DATA_TRUE_occ_v$sampMatbased)
label_plants=paste0("data_",plants)
refs=unique(DB_DATA_TRUE_occ_v$Ref)

for ( i in 1:length(plants)) {
temp=as.data.frame.matrix(xtabs(~Ref+paramType,data=subset(DB_DATA_TRUE_occ_v,sampMatbased==plants[i])))
temp=data.frame(ref=row.names(temp),temp)
XLConnect::writeWorksheetToFile("PlantTox_occurence_stat.xls",temp,label_plants[i])
}













res=list()
z=1
for ( j in unique(DB_DATA_TRUE_occ$sampMatbased)) {
  xtabs(~Ref+paramType,suset(DB_DATA_TRUE_occ,para[which(paramType)])
}  
  
XLConnect::writeWorksheetToFile("PlantTox_occurence_stat.xls",j)
