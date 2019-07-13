library(readxl)
library(fitdistrplus)
library(doBy)
library(XLConnect)
library(ggplot2)

setwd("/home/alf/Scrivania/lav_michyf/lav_PDI/child_DON")

micotoxins=excel_sheets("DON_BIOMARKER_CHILD.xlsx")

pdi_func_DON=function(x,escr=70,vol_urine=2,weight=70) {return(x*(vol_urine/weight)*(100/escr))}
pdi_func=function(x,escr=70,vol_urine=2,weight=70) {return(x*(vol_urine/weight)*(100/escr))}
sumfit=function(x) {return(data.frame(par=x$estimate[1],errpar=x$estimate[2],aic=x$aic,names=x$distname))}

####################################################################################
# from Piero Toscano
# C = i valori di concentrazione che trovi nel file xls
# V = quantità urine per adulti 2Litri (EFSA, NDA) ma abbiamo ref con range [1, 2.4] Litri
# W = adulti 70 kg, bambini fino a 9 anni = 30kg. 
# 
# Nel file xls per tutti i dati di Brera 2015 abbiamo pure il peso
# 

####################################################################################
# DON analisys
# E (DON) = 72.3% Turner, 2010; 68% Warth 2013; 70% EFSA, 2017

EDON=70

DON_data=as.data.frame(read_excel("DON_BIOMARKER_CHILD.xlsx",1))
names(DON_data)[5:6]=c("weight","biomark")
table(DON_data$Sex,DON_data$Country)


Full_mat_DON=DON_data[which(DON_data$Country %in% c("IT","NO","UK")==T),]
Full_mat_DON$weight=as.numeric(Full_mat_DON$weight)
DON_data_nation=summaryBy(biomark+weight~Sex+Country,data=Full_mat_DON, FUN=c(mean,sd), na.rm=TRUE)
DON_data_nation_full=summaryBy(biomark+weight~Sex,data=Full_mat_DON, FUN=c(mean,sd), na.rm=TRUE)

aa=cbind(DON_data_nation_full[1],c("full","full"),DON_data_nation_full[,c(2:5)])
bb=c("full","full",apply(DON_data_nation_full[,2:5],2,function(x) mean(x,na.rm=T)))
names(bb)=names(DON_data_nation)
names(aa)=names(DON_data_nation)
mat_meta=rbind(bb,aa,DON_data_nation)




M_data_DON=Full_mat_DON[Full_mat_DON$Sex=="Male",]
F_data_DON=Full_mat_DON[Full_mat_DON$Sex=="Female",]
mat_meta$biomark.sd[1]=sd(Full_mat_DON$biomark,na.rm=T)
mat_meta$biomark.sd[2]=sd(F_data_DON$biomark,na.rm=T)
mat_meta$biomark.sd[3]=sd(M_data_DON$biomark,na.rm=T)
mat_meta$weight.sd[1]=sd(Full_mat_DON$weight,na.rm=T)
mat_meta$weight.sd[2]=sd(F_data_DON$weight,na.rm=T)
mat_meta$weight.sd[3]=sd(M_data_DON$weight,na.rm=T)
mat_meta$biomark.mean=as.numeric(mat_meta$biomark.mean)
mat_meta$weight.mean=as.numeric(mat_meta$weight.mean)
mat_meta$weight.sd=as.numeric(mat_meta$weight.sd)
mat_meta$biomark.sd=as.numeric(mat_meta$biomark.sd)
M_data_DON_nation=split(M_data_DON,M_data_DON$Country) # IT 1 NO 2 UK 3
F_data_DON_nation=split(F_data_DON,F_data_DON$Country)

##############################################################################
# collect DON data

data_DON_ls=list()
data_DON_ls$data_DON_tot=as.numeric(na.omit(Full_mat_DON$biomark))
data_DON_ls$data_M_DON_tot=as.numeric(na.omit(M_data_DON$biomark))
data_DON_ls$data_F_DON_tot=as.numeric(na.omit(F_data_DON$biomark))
data_DON_ls$data_M_DON_IT=as.numeric(na.omit(M_data_DON_nation[[1]]$biomark))
data_DON_ls$data_F_DON_IT=as.numeric(na.omit(F_data_DON_nation[[1]]$biomark))
data_DON_ls$data_M_DON_NO=as.numeric(na.omit(M_data_DON_nation[[2]]$biomark))
data_DON_ls$data_F_DON_NO=as.numeric(na.omit(F_data_DON_nation[[2]]$biomark))
data_DON_ls$data_M_DON_UK=as.numeric(na.omit(M_data_DON_nation[[3]]$biomark))
data_DON_ls$data_F_DON_UK=as.numeric(na.omit(F_data_DON_nation[[3]]$biomark))
saveRDS(data_DON_ls,"data_DON_ls.rds")

set.seed(2)

##############################################################################################################################
names_DON=c("DON_tot",
          "F_DON_tot",
          "M_DON_tot",
          "F_DON_IT",
          "F_DON_NO",
          "F_DON_UK",
          "M_DON_IT",
          "M_DON_NO",
          "M_DON_UK")


res_DON=list(
DON_tot=fitdist(data_DON_ls$data_DON_tot,"norm"),
F_DON_tot=fitdist(data_DON_ls$data_F_DON_tot,"norm"),
M_DON_tot=fitdist(data_DON_ls$data_M_DON_tot,"norm"),
F_DON_IT=fitdist(data_DON_ls$data_F_DON_IT,"norm"),
F_DON_NO=fitdist(data_DON_ls$data_F_DON_NO,"norm"),
F_DON_UK=fitdist(data_DON_ls$data_F_DON_UK,"norm"),
M_DON_IT=fitdist(data_DON_ls$data_M_DON_IT,"norm"),
M_DON_NO=fitdist(data_DON_ls$data_M_DON_NO,"norm"),
M_DON_UK=fitdist(data_DON_ls$data_M_DON_UK,"norm"),
DON_tote=fitdist(data_DON_ls$data_DON_tot,"exp"),
F_DON_tote=fitdist(data_DON_ls$data_F_DON_tot,"exp"),
M_DON_tote=fitdist(data_DON_ls$data_M_DON_tot,"exp"),
F_DON_ITe=fitdist(data_DON_ls$data_F_DON_IT,"exp"),
F_DON_NOe=fitdist(data_DON_ls$data_F_DON_NO,"exp"),
F_DON_UKe=fitdist(data_DON_ls$data_F_DON_UK,"exp"),
M_DON_ITe=fitdist(data_DON_ls$data_M_DON_IT,"exp"),
M_DON_NOe=fitdist(data_DON_ls$data_M_DON_NO,"exp"),
M_DON_UKe=fitdist(data_DON_ls$data_M_DON_UK,"exp")
)

df_res_DON=data.frame(names=names(res_DON),do.call("rbind",lapply(res_DON,FUN=sumfit)))
row.names(df_res_DON)=NULL

for (i in 1:9) {

png(paste0(names_DON[i],".png"))
f=i+9
denscomp(res_DON[c(i,f)],legendtext = c("Normal", "Exponential"),
         main = paste("Fitting",names_DON[i],"biomarker data"), xlab = "microg/L")

dev.off()

}



exp_par_DON=df_res_DON$par[10:18] 
res_pdi_DON_05=list()
for ( i in 1:length(exp_par_DON)) {
                                  temp_PDI=sapply(rexp(10000,exp_par_DON[i]),FUN=function(x){ pdi_func(x,vol_urine = 0.5,weight=mat_meta$weight.mean[i])})
                                  res_pdi_DON_05[[i]]=as.numeric(c(t.test(temp_PDI)$estimate,t.test(temp_PDI)$conf.int))
}

df_pdi_DON_05=data.frame(name=names_DON,do.call("rbind",res_pdi_DON_05))
names(df_pdi_DON_05) [2:4]=c("mean","conf.int.inf","conf.int.sup")




res_pdi_DON_17=list()
for ( i in 1:length(exp_par_DON)) {
  temp_PDI=sapply(rexp(10000,exp_par_DON[i]),FUN=function(x){ pdi_func(x,vol_urine = 1.7,weight=mat_meta$weight.mean[i])})
  res_pdi_DON_17[[i]]=as.numeric(c(t.test(temp_PDI)$estimate,t.test(temp_PDI)$conf.int))
}

df_pdi_DON_17=data.frame(name=names_DON,do.call("rbind",res_pdi_DON_17))
names(df_pdi_DON_17) [2:4]=c("mean","conf.int.inf","conf.int.sup")


df_pdi_DON_05$vol_urine=0.5
df_pdi_DON_17$vol_urine=1.7

temp=rbind(df_pdi_DON_05,df_pdi_DON_17)
temp_fit=rbind(df_res_DON)

file.remove("PDI_stat_bimbi.xls")
XLConnect::writeWorksheetToFile("PDI_stat_bimbi.xls",temp,"PDI")
XLConnect::writeWorksheetToFile("PDI_stat_bimbi.xls",temp_fit,"fit stats")
XLConnect::writeWorksheetToFile("PDI_stat_bimbi.xls",mat_meta,"metadati")

##############################################################################################################################



