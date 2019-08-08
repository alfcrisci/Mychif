################################################################################
library("XLConnect")
library("doBy")
################################################################################

setwd("")

con=read.csv("Consumption_data.csv",sep=";")
mat=con
attach(con)

################################################################################################################################################
# definitions

mycotoxs=c("AcDONs","15AcDON","3AcDON","15AcDON+3AcDON","aZEL","aZEL4G ","DON","DON3G","FB1","FB1+FB2","FB2","FB3","FBs","NIV","OTA","?ZEL","?ZEL4G","2+HT2","ZEN","ZEN4G","ZEN4S","ZEN16G")
sampMat=c("wheat","oat","barley","triticale","spelt","sorghum","rye","cereals","rice","buckwheat","maize")

sumfun <- function(x, ...){
  c(m=mean(x, na.rm=TRUE, ...),len=length(x))
}
################################################################################################################################################

res_cons=data.frame(country,COUNTRY_cons_average_adolescent=Adolescent_mean/Adolescent_weight,
                            COUNTRY_cons_p95_adolescent=Adolescent_95/Adolescent_weight,
                            COUNTRY_cons_average_adult=Adult_mean/Adult_weight,
                            COUNTRY_cons_p95_adult=Adult_95/Adult_weight,
                            COUNTRY_cons_average_elderly=Elderly_mean/Elderly_weight,
                            COUNTRY_cons_p95_elderly=Elderly_95/Elderly_weight)


################################################################################

ALL_cons_average_adolescent=mean(Adolescent_mean,na.rm=T)/mean(Adolescent_weight,na.rm=T)
ALL_cons_p95_adolescent=mean(Adolescent_95,na.rm=T)/mean(Adolescent_weight,na.rm=T)
ALL_cons_average_adult=mean(Adult_mean,na.rm=T)/mean(Adult_weight,na.rm=T)
ALL_cons_p95_adult=mean(Adult_95,na.rm=T)/mean(Adult_weight,na.rm=T)
ALL_cons_average_elderly=mean(Elderly_mean,na.rm=T)/mean(Elderly_weight,na.rm=T)
ALL_cons_p95_elderly=mean(Elderly_95,na.rm=T)/mean(Elderly_weight,na.rm=T)

mat_all=data.frame(cbind(ALL_cons_average_adolescent,
          ALL_cons_p95_adolescent,
          ALL_cons_average_adult,
          ALL_cons_p95_adult,
          ALL_cons_average_elderly,
          ALL_cons_p95_elderly))
################################################################################
#  read occurrence&co_final_RP 


wb <- loadWorkbook("all_cereals_FINALE_RP_29_07.xls")
db_data=XLConnect::readWorksheet(wb,1)
var2num=c("resLOD","resLOQ","meanTot","meanPos","median","min","max","MAJ_LOD","IQRmin","IQRmax","Concentration","Co_occurrence") 
for ( i in var2num){db_data[,i]=as.numeric(gsub(",",".",db_data[,i]))}
saveRDS(db_data,"all_cereals_FINALE_RP_29_07.rds")

##############################################################################################################################################à



db_food=db_data[which(db_data$sampMatType == "food"),]
db_food_occ=db_food[which(db_food$Co_occurrence == 0),]
db_food_occ_clean=db_food_occ[which(nchar(db_food_occ$sampCountry) == 2),]

db_food_coocc=db_food[which(db_food$Co_occurrence == 1),]
db_food_coocc_clean=db_food_coocc[which(nchar(db_food_coocc$sampCountry) == 2),]


########################################################################################################################################################################################
###  WRITE ALL DATA MEDIUM LOW AND UPPER BOUND (MB, LB, UB) occurence & co-occurence

db_food_occ_clean$LB=db_food_occ_clean$MB=db_food_occ_clean$UB=NA
db_food_coocc_clean$LB=db_food_coocc_clean$MB=db_food_coocc_clean$UB=NA


id=which(db_food_occ_clean$meanTot>0 | 
         db_food_occ_clean$meanPos>0 | 
         db_food_occ_clean$Concentration>0)

db_food_occ_clean$LB[id]=db_food_occ_clean$MB[id]=db_food_occ_clean$UB[id]=apply(db_food_occ_clean[id,c("meanTot","meanPos","Concentration")],1,FUN=function(x) max(x,na.rm=T))


id_A=which(db_food_occ_clean$meanTot==-1 | 
           db_food_occ_clean$meanPos==-1 | 
           db_food_occ_clean$Concentration==-1)

db_food_occ_clean$LB[id_A]=0
db_food_occ_clean$MB[id_A]=0.5*db_food_occ_clean$resLOD[id_A]
db_food_occ_clean$UB[id_A]=db_food_occ_clean$resLOD[id_A]


id_B=which(db_food_occ_clean$meanTot==-2 | 
             db_food_occ_clean$meanPos==-2 | 
             db_food_occ_clean$Concentration==-2)


db_food_occ_clean$LB[id_B]=db_food_occ_clean$resLOD[id_B]
db_food_occ_clean$MB[id_B]=0.5*db_food_occ_clean$resLOQ[id_B]
db_food_occ_clean$UB[id_B]=db_food_occ_clean$resLOQ[id_B]

##########################################################################################################################################################
# co-occurences

db_food_coocc_clean$LB=db_food_coocc_clean$MB=db_food_coocc_clean$UB=NA
db_food_coocc_clean$LB=db_food_coocc_clean$MB=db_food_coocc_clean$UB=NA


id=which(db_food_coocc_clean$meanTot>0 | 
           db_food_coocc_clean$meanPos>0 | 
           db_food_coocc_clean$Concentration>0)

db_food_coocc_clean$LB[id]=db_food_coocc_clean$MB[id]=db_food_coocc_clean$UB[id]=apply(db_food_coocc_clean[id,c("meanTot","meanPos","Concentration")],1,FUN=function(x) max(x,na.rm=T))


id_A=which(db_food_coocc_clean$meanTot==-1 | 
             db_food_coocc_clean$meanPos==-1 | 
             db_food_coocc_clean$Concentration==-1)

db_food_coocc_clean$LB[id_A]=0
db_food_coocc_clean$MB[id_A]=0.5*db_food_coocc_clean$resLOD[id_A]
db_food_coocc_clean$UB[id_A]=db_food_coocc_clean$resLOD[id_A]


id_B=which(db_food_coocc_clean$meanTot==-2 | 
             db_food_coocc_clean$meanPos==-2 | 
             db_food_coocc_clean$Concentration==-2)


db_food_coocc_clean$LB[id_B]=db_food_coocc_clean$resLOD[id_B]
db_food_coocc_clean$MB[id_B]=0.5*db_food_coocc_clean$resLOQ[id_B]
db_food_coocc_clean$UB[id_B]=db_food_coocc_clean$resLOQ[id_B]

##########################################################################################################################################################

BALKGREE=c("HR","RS","BA","MK","SI","BG","GR","CY")
BALTIC=c("EE","LV","LT")
EEUROPA=c("SK","HU","RO","PL")
ENG=c("IE","GB")
IBERIA=c("ES","PT")
ITALIA=c("IT","MT")
MIDDLEEU=c("FR","CH","AT","BE","NL","DE","CZ","DK")
SCAND=c("NO","SE","FI")
ALL=unique(c(BALKGREE,BALTIC,EEUROPA,ENG,IBERIA,ITALIA,MIDDLEEU,SCAND))


countries_occ=unlist(strsplit(gsub("NaN;","",paste(unique(db_food_occ_clean$sampCountry),collapse=";")),";"))
countries_coocc=unlist(strsplit(gsub("NaN;","",paste(unique(db_food_coocc_clean$sampCountry),collapse=";")),";"))

########################################################################################################################################
# Occurence workflow


db_food_occ_clean=db_food_occ_clean[which(db_food_occ_clean$sampMatbased %in% sampMat ==T),]
table_occ_clean=summaryBy(UB+MB+LB~ paramType+sampMatbased+sampCountry,data=db_food_occ_clean,FUN=sumfun)
table_occ_clean$defparamType=table_occ_clean$paramType
table_occ_clean$defparamType[grep("DON",table_occ_clean$defparamType)]="DONeq"

file.remove("Table_occurence_full.xls")
XLConnect::writeWorksheetToFile("Table_occurence_full.xls",table_occ_clean,"full")

selected_multiple=c("DONeq")
tempM=table_occ_clean[which(table_occ_clean$defparamType %in% selected_multiple ==T),]
tempMs=summaryBy(UB.m+MB.m+LB.m+LB.len~sampCountry+defparamType+sampMatbased,data=tempM,FUN=sum)
names(tempMs)=gsub(".sum","",names(tempMs))
tempMs=na.omit(tempMs)

selected_single=c("NIV","T2+HT2")
table_occ_clean_def=table_occ_clean[which(table_occ_clean$defparamType %in% selected_single ==T),]
table_occ_clean_def=na.omit(table_occ_clean_def)[,names(tempMs)]
table_occ_clean_def=rbind(table_occ_clean_def,tempMs)
XLConnect::writeWorksheetToFile("Table_occurence_full.xls",table_occ_clean_def,"def")



############################################################################################################

findcountry_occ=function(x,mat=table_occ_clean_def) grep(x,mat$sampCountry)

macroaree_ls_occ=list(
  id_balkgree=as.numeric(unlist(sapply(BALKGREE,FUN=findcountry_occ))),
  id_baltic=as.numeric(unlist(sapply(BALTIC,FUN=findcountry_occ))),
  id_eeuropa=as.numeric(unlist(sapply(EEUROPA,FUN=findcountry_occ))),
  id_eng=as.numeric(unlist(sapply(ENG,FUN=findcountry_occ))),
  id_italia=as.numeric(unlist(sapply(ITALIA,FUN=findcountry_occ))),
  id_iberia=as.numeric(unlist(sapply(IBERIA,FUN=findcountry_occ))),
  id_middleeu=as.numeric(unlist(sapply(MIDDLEEU,FUN=findcountry_occ))),
  id_scand=as.numeric(unlist(sapply(SCAND,FUN=findcountry_occ))),
  id_all=as.numeric(unlist(sapply(ALL,FUN=findcountry_occ)))
)

macroaree_ls_occ=macroaree_ls_occ[as.numeric(which(unlist(lapply(macroaree_ls_occ,length))>0))]

res_macro_exp_ls=list()


for ( i in 1:length(macroaree_ls_occ)) 
  
  {
  tempexp=table_occ_clean_def[macroaree_ls_occ[[i]],]
  res_macro_exp_ls[[i]]=cbind(name=names(macroaree_ls_occ)[i],summaryBy(UB.m+MB.m+LB.m~defparamType,data=tempexp,FUN = mean))
  
  }

#####################################################################################################


res_macro_exp_df=do.call("rbind",res_macro_exp_ls)

res_macro_exp=cbind(res_macro_exp_df,mat_all,mat_all,mat_all)
names(res_macro_exp)[6:11]=paste0("UB_EXP_",names(res_macro_exp)[6:11])
names(res_macro_exp)[12:17]=paste0("MB_EXP_",names(res_macro_exp)[12:17])
names(res_macro_exp)[18:23]=paste0("LB_EXP_",names(res_macro_exp)[18:23])


res_macro_exp[,6:11]=res_macro_exp[,6:11]*res_macro_exp[,3]/1000
res_macro_exp[,12:17]=res_macro_exp[,12:17]*res_macro_exp[,4]/1000
res_macro_exp[,18:23]=res_macro_exp[,18:23]*res_macro_exp[,5]/1000

res_macro_exp$EF=res_macro_exp$PDI=NA

id_DONeq=which(res_macro_exp$defparamType=="DONeq")
id_NIV=which(res_macro_exp$defparamType=="NIV")
id_T2_HT2=which(res_macro_exp$defparamType=="T2+HT2")

#############################################################################################################
###PDI per DON, NIV, T2+HT2

PDI_T2_HT2 = 3.33
PDI_DON = 110
PDI_NIV = 350

### EQUIVALENT FACTOR , sarebbero i Relative potency factor

EF_T_HT2=1
EF_DON=0.03
EF_NIV=0.0095

res_macro_exp$EF[id_DONeq]=EF_DON
res_macro_exp$PDI[id_DONeq]=PDI_DON

res_macro_exp$EF[id_NIV]=EF_NIV
res_macro_exp$PDI[id_NIV]=PDI_NIV
res_macro_exp$EF[id_T2_HT2]=EF_T_HT2
res_macro_exp$PDI[id_T2_HT2]=PDI_T2_HT2

res_macro_exp[,gsub("EXP","EF",names(res_macro_exp)[6:23])]=NA
res_macro_exp[, grep("_EF",names(res_macro_exp))]=res_macro_exp[,6:23]*res_macro_exp$EF # 26:43


final=res_macro_exp[,c(1:2,26:43)]
final_res_sum=summaryBy(.~name,data=final,FUN=sum)
final_moe=3.33/final_res_sum[,2:length(final_res_sum)]
final_moe=cbind(final_res_sum$name,final_moe)


names(final_res_sum)=gsub("EF","SUM",names(final_moe))
names(final_moe)=gsub("EF","MOE",names(final_moe))
names(final_moe)=gsub(".sum","",names(final_moe))


XLConnect::writeWorksheetToFile("Table_occurence_full.xls",final,"exposure per macroaree")
XLConnect::writeWorksheetToFile("Table_occurence_full.xls",final_res_sum,"SUM per macroaree")
XLConnect::writeWorksheetToFile("Table_occurence_full.xls",final_moe,"MOE per macroaree")

###########################################################################################################################################ààà
# Co - occurence work flows

db_food_coocc_clean=db_food_coocc_clean[which(db_food_coocc_clean$sampMatbased %in% sampMat ==T),]
table_coocc_clean=summaryBy(UB+MB+LB~ paramType+sampMatbased+sampCountry,data=db_food_coocc_clean,FUN=sumfun)
table_coocc_clean$defparamType=table_coocc_clean$paramType
table_coocc_clean$defparamType[grep("DON",table_coocc_clean$defparamType)]="DONeq"

file.remove("Table_cooccurence_full.xls")
XLConnect::writeWorksheetToFile("Table_cooccurence_full.xls",table_coocc_clean,"full")

selected_multiple=c("DONeq")
tempM=table_coocc_clean[which(table_coocc_clean$defparamType %in% selected_multiple ==T),]
tempMs=summaryBy(UB.m+MB.m+LB.m+LB.len~sampCountry+defparamType+sampMatbased,data=tempM,FUN=sum)
names(tempMs)=gsub(".sum","",names(tempMs))
tempMs=na.omit(tempMs)

selected_single=c("NIV","T2+HT2")
table_coocc_clean_def=table_coocc_clean[which(table_coocc_clean$defparamType %in% selected_single ==T),]
table_coocc_clean_def=na.omit(table_coocc_clean_def)[,names(tempMs)]
table_coocc_clean_def=rbind(table_coocc_clean_def,tempMs)
XLConnect::writeWorksheetToFile("Table_cooccurence_full.xls",table_coocc_clean_def,"def")


########################################################################################################
findcountry_coocc=function(x,mat=table_coocc_clean_def) grep(x,mat$sampCountry)

macroaree_ls_coocc=list(
  id_balkgree=as.numeric(unlist(sapply(BALKGREE,FUN=findcountry_coocc))),
  id_baltic=as.numeric(unlist(sapply(BALTIC,FUN=findcountry_coocc))),
  id_eeuropa=as.numeric(unlist(sapply(EEUROPA,FUN=findcountry_occ))),
  id_eng=as.numeric(unlist(sapply(ENG,FUN=findcountry_coocc))),
  id_iberia=as.numeric(unlist(sapply(IBERIA,FUN=findcountry_coocc))),
  id_italia=as.numeric(unlist(sapply(ITALIA,FUN=findcountry_coocc))),
  id_middleeu=as.numeric(unlist(sapply(MIDDLEEU,FUN=findcountry_coocc))),
  id_scand=as.numeric(unlist(sapply(SCAND,FUN=findcountry_coocc))),
  id_all=as.numeric(unlist(sapply(ALL,FUN=findcountry_coocc)))
)

macroaree_ls_coocc=macroaree_ls_coocc[as.numeric(which(unlist(lapply(macroaree_ls_coocc,length))>0))]

res_macro_exp_ls=list()


for ( i in 1:length(macroaree_ls_coocc)) 
  
{
  tempexp=table_coocc_clean_def[macroaree_ls_coocc[[i]],]
  res_macro_exp_ls[[i]]=cbind(name=names(macroaree_ls_coocc)[i],summaryBy(UB.m+MB.m+LB.m~defparamType,data=tempexp,FUN = mean))
  
}

#####################################################################################################

res_macro_exp_ls[[2]]=res_macro_exp_ls[[2]][1:3,]

res_macro_exp_df=do.call("rbind",res_macro_exp_ls)

res_macro_exp=cbind(res_macro_exp_df,mat_all,mat_all,mat_all)
names(res_macro_exp)[6:11]=paste0("UB_EXP_",names(res_macro_exp)[6:11])
names(res_macro_exp)[12:17]=paste0("MB_EXP_",names(res_macro_exp)[12:17])
names(res_macro_exp)[18:23]=paste0("LB_EXP_",names(res_macro_exp)[18:23])


res_macro_exp[,6:11]=res_macro_exp[,6:11]*res_macro_exp[,3]/1000 # UB
res_macro_exp[,12:17]=res_macro_exp[,12:17]*res_macro_exp[,4]/1000 # MB
res_macro_exp[,18:23]=res_macro_exp[,18:23]*res_macro_exp[,5]/1000 # LB

res_macro_exp$EF=res_macro_exp$PDI=NA

id_DONeq=which(res_macro_exp$defparamType=="DONeq")
id_NIV=which(res_macro_exp$defparamType=="NIV")
id_T2_HT2=which(res_macro_exp$defparamType=="T2+HT2")

#############################################################################################################
###PDI per DON, NIV, T2+HT2

PDI_T2_HT2 = 3.33
PDI_DON = 110
PDI_NIV = 350

### EQUIVALENT FACTOR , sarebbero i Relative potency factor

EF_T_HT2=1
EF_DON=0.03
EF_NIV=0.0095

res_macro_exp$EF[id_DONeq]=EF_DON
res_macro_exp$PDI[id_DONeq]=PDI_DON

res_macro_exp$EF[id_NIV]=EF_NIV
res_macro_exp$PDI[id_NIV]=PDI_NIV
res_macro_exp$EF[id_T2_HT2]=EF_T_HT2
res_macro_exp$PDI[id_T2_HT2]=PDI_T2_HT2

res_macro_exp[,gsub("EXP","EF",names(res_macro_exp)[6:23])]=NA
res_macro_exp[, grep("_EF",names(res_macro_exp))]=res_macro_exp[,6:23]*res_macro_exp$EF # 26:43


final=res_macro_exp[,c(1:2,26:43)]
final_res_sum=summaryBy(.~name,data=final,FUN=sum)
final_moe=3.33/final_res_sum[,2:length(final_res_sum)]
final_moe=cbind(final_res_sum$name,final_moe)
names(final_res_sum)=gsub("EF","SUM",names(final_moe))
names(final_moe)=gsub("EF","MOE",names(final_moe))
names(final_moe)=gsub(".sum","",names(final_moe))

XLConnect::writeWorksheetToFile("Table_cooccurence_full.xls",final,"exposure per macroaree")
XLConnect::writeWorksheetToFile("Table_cooccurence_full.xls",final_res_sum,"SUM per macroaree")
XLConnect::writeWorksheetToFile("Table_cooccurence_full.xls",final_moe,"MOE per macroaree")
###########################################################################################################################