#####################################################################################################################################
setwd("")

setwd("/home/alf/Scrivania/lav_michyf/repositories/Mychif/occurence")

source("load_lib.r")
source("aux_mycosources.r")

#####################################################################################################################################
# read data

data_occ=readRDS("data/data_occ.rds")

#####################################################################################################################################
# data_occ=read_excel("data/all_cereals_FINALE_RP.xls",1)
# data_occ[data_occ=="NaN"]=NA
# data_occ$meanTot=as.numeric(data_occ$meanTot)
# data_occ$Concentration=as.numeric(data_occ$Concentration)
# data_occ=as.data.frame(data_occ)
# saveRDS(data_occ,"data/data_occ.rds")

#####################################################################################################################################
#' @Roberta
#' 
#' - le ref sono 219, inizialmente avevo buttato una tra Co-occurrence of type A e Co-occurrence of type-A, 
#'   ma in realtà sono due pubblicazioni differenti. Trovi la lista sempre in questo file. Chiedi ad Ana di risistemare tutte le sue ref in modo univoco (titolo primo autore ed anno).
#' 
#' - le crop per ora sono 16, eventualmente decidiamo in post-analisi quali accorpare o escludere
#' 
#' - le tossine e parenti sono 121, eventualmente decidiamo in post-analisi quali accorpare o escludere
#' 
#' @Alfonso
#' 
#' Prima analisi
#' 
#' per ogni crop ed ogni tossina estrarre e graficare (ECD e tutto quello che ti sembra coerente, 
#' con i vari check sulla tipologia di distribuzione) i valori di concentrazione
#' Nei valori di concentrazione devi contemplare sia il valore di concentrazione 
#' se il datasamplelevel=1 che il valore di meantot se il datasamplelevel=0
#' Se per il datasamplelevel=0 non esiste il valore di meantot 
#' puoi includere i record che hanno il valore di meanpos
#' 
#' Se il valore di concentrazione o di meantot o di meanpos è =-1, tieni traccia sia di questo che del LOD. 
#' Significa che i valori sono inferiori a LOD. In qualche modo dovremo graficarli visto che i LOD sono sempre differenti 
#' (dipendono dalla tossina e dallo strumento di analisi)
#' 
#' Nei dati sono presenti anche i -2 (<LOQ) che per adesso puoi considerare come NaN.

#####################################################################################################################################

#####################################################################################################################################
# [1] "DB_origin"           "sampCountryorigin"   "sampContinentorigin" "sampContinent"       "paramType"           "sampCountry"        
# [7] "sampRegion"          "sampInfo_latitude"   "sampInfo_longitude"  "sampInfo_altitude"   "sampYear"            "sampYearIncreas"    
# [13] "sampSize"            "sampMethod"          "sampPoint"           "sampPointInfo"       "sampMatType"         "sampMatbased"       
# [19] "sampMatInfo"         "sampMatCode"         "growingSystem"       "anMethRefId"         "resUnit"             "resLOD"             
# [25] "resLODinfo"          "resLOQ"              "resLOQinfo"          "exprResPerc"         "ToTresValUncertSD"   "meanTot"            
# [31] "POSresValUncertSD"   "meanPos"             "median"              "medianinfo"          "min"                 "max"                
# [37] "MAJ_LOD"            "IQRmin"              "IQRmax"              "Concentration"       "DataSampLevel"       "Co_occurrence"      
# [43] "Ref"
###########################################################################
##########################################################################
# working on GEO

data_auth_geo=c("DB_origin","sampContinent","sampCountryorigin","sampContinentorigin","sampCountry","sampRegion")
DB_AUTH_GEO=data_occ[data_auth_geo]

#   to uppercase

DB_AUTH_GEO=as.data.frame(apply(DB_AUTH_GEO,2,toupper))
DB_AUTH_GEO$sampCountryorigin=gsub("UK","GB",DB_AUTH_GEO$sampCountryorigin)
DB_AUTH_GEO$sampCountry=gsub("UK","GB",DB_AUTH_GEO$sampCountry)
saveRDS(DB_AUTH_GEO,"data/DB_AUTH_GEO.rds")

continents=c("AFRICA","EUROPE","SOUTHAMERICA","ASIA")
countries_orig=unique(as.vector(do.call("rbind",sapply(levels(as.factor(c(DB_AUTH_GEO$sampCountryorigin))),FUN=function(x) strsplit(x, ";")))))
countries=unique(as.vector(do.call("rbind",sapply(levels(as.factor(c(DB_AUTH_GEO$sampCountry))),FUN=function(x) strsplit(x, ";")))))
countries=unique(c(countries_orig,countries))


res=list()
res_nation=list()
for ( i in 1:length(countries)) {
  res[[i]]=length(grep(countries[i],DB_AUTH_GEO$sampCountryorigin,fixed=T))
  res_nation[[i]]=countrycode(countries[i], 'iso2c','country.name')

}

res_country_orig=data.frame(code=countries,country=unlist(res_nation),Ndata_orig=unlist(res),stringsAsFactors = F)

res=list()
res_nation=list()
for ( i in 1:length(countries)) {
  res[[i]]=length(grep(countries[i],DB_AUTH_GEO$sampCountry,fixed=T))
  res_nation[[i]]=countrycode(countries[i], 'iso2c','country.name')
  
}

res_country=data.frame(code=countries,country=unlist(res_nation),Ndata=unlist(res),stringsAsFactors = F)

res_country_orig$Ndata=res_country$Ndata


res=list()
for ( i in 1:length(continents)) {
  res[[i]]=length(grep(continents[i],DB_AUTH_GEO$sampContinent,fixed=T))
  
}

res_cont=data.frame(code=continents,Ndata=unlist(res),stringsAsFactors = F)


res=list()
for ( i in 1:length(continents)) {
  res[[i]]=length(grep(continents[i],DB_AUTH_GEO$sampContinentorigin,fixed=T))

}

res_cont_orig=data.frame(code=continents,Ndata_orig=unlist(res),Ndata=res_cont$Ndata,stringsAsFactors = F)


# results

saveRDS(res_cont_orig,"data/res_continent_orig.rds")
saveRDS(res_country_orig,"data/res_country_orig.rds")
file.remove("data/mycotoxins_geo.xls")
XLConnect::writeWorksheetToFile("data/mycotoxins_geo.xls",res_cont_orig,"continent")
XLConnect::writeWorksheetToFile("data/mycotoxins_geo.xls",res_country_orig,"nations")


#################################################################################################################################
# info on data values by selected data fields

info_data_values=c("sampSize","resUnit","resLOD","resLODinfo","exprResPerc","ToTresValUncertSD","medianinfo",
                   "POSresValUncertSD","MAJ_LOD","Co_occurrence","Ref")

DB_INFODATA=as.data.frame(data_occ[info_data_values])
saveRDS(DB_INFODATA,"data/DB_INFODATA.rds")


# DB data_stratum

data_stratum=c("paramType","sampMatbased","sampMatType","sampMatInfo","sampMatCode","growingSystem","anMethRefId","sampMethod")
DB_STRATUM=as.data.frame(data_occ[data_stratum])
saveRDS(DB_STRATUM,"data/DB_STRATUM.rds")

# DB data values

data_values=c("sampMatbased","paramType","meanTot","Concentration","meanPos","median","min","max","sampSize","Co_occurrence","resUnit","resLOD","Ref")
DB_DATA_FULL=as.data.frame(data_occ[data_values])


# working on missing data management


id_overlod=which(DB_DATA_FULL$meanTot >=-1)       # 2556
length(id_overlod)
id_overlod2=which(DB_DATA_FULL$Concentration>=-1) # 3136
length(id_overlod2)
id_overlod_full=unique(c(id_overlod,id_overlod2)) # 5598
length(id_overlod_full)
DB_DATA_TRUE=DB_DATA_FULL[id_overlod_full,]       # 
DB_DATA_TRUE$Concentration=as.numeric(DB_DATA_TRUE$Concentration)
DB_DATA_TRUE$meanPos=as.numeric(DB_DATA_TRUE$meanPos)
DB_DATA_TRUE$median=as.numeric(DB_DATA_TRUE$median)
DB_DATA_TRUE$min=as.numeric(DB_DATA_TRUE$min)                
DB_DATA_TRUE$max=as.numeric(DB_DATA_TRUE$max) 
DB_DATA_TRUE$sampSize=as.numeric(DB_DATA_TRUE$sampSize)

saveRDS(DB_DATA_TRUE,"data/DB_DATA_TRUE.rds")
#################################################################################################################################
DB_DATA_TRUE=readRDS("data/DB_DATA_TRUE.rds")


by_plant=split(DB_DATA_TRUE,DB_DATA_TRUE$sampMatbased)

plants=unique(DB_DATA_TRUE$sampMatbased) # 
mycotoxins=unique(DB_DATA_TRUE$paramType) # mycotoxins

table_plants_mycotoxin=as.data.frame.array(table(DB_DATA_TRUE$sampMatbased, 
                                                 DB_DATA_TRUE$paramType))

mycorank=sort(colSums(table_plants_mycotoxin), decreasing = T)

plantrank=sort(rowSums(table_plants_mycotoxin),decreasing = T)

aux_info_mycotox=data_occ[id_overlod_full,c("sampYearIncreas","Co_occurrence","Ref")]

aux_info_mycotox$agepaper=2019-as.numeric(unlist(lapply(regmatches(aux_info_mycotox$Ref,gregexpr("[[:digit:]]+", aux_info_mycotox$Ref)),function(x) x[1])))

years_paper=aux_info_mycotox$agepaper

png("images/Years_paper.png")
hist(years_paper,main="Distribution age of mycotoxin occurence papers \n( Ref year= 2019)")
dev.off()

DB_DATA_TRUE$agepaper=aux_info_mycotox$agepaper
DB_DATA_TRUE$havebounds=ifelse((!is.na(DB_DATA_TRUE$min) & !is.na(DB_DATA_TRUE$max)),1,0)
DB_DATA_TRUE$havemeanpos=ifelse((!is.na(DB_DATA_TRUE$meanPos)),1,0)
DB_DATA_TRUE$havedata=ifelse((!is.na(DB_DATA_TRUE$meanTot) || !is.na(DB_DATA_TRUE$Concentration)),1,0)
DB_DATA_TRUE$norm_agepaper=range01(DB_DATA_TRUE$agepaper)
DB_DATA_TRUE$norm_sampSize=range01(DB_DATA_TRUE$sampSize)



saveRDS(DB_DATA_TRUE,"data/DB_DATA_TRUE.rds")
saveRDS(temp_plan_raw,"data/template.rds")
saveRDS(mycorank,"data/mycorank.rds")
saveRDS(plantrank,"data/plantrank.rds")

DB_DATA_TRUE_co_occur=DB_DATA_TRUE[which(DB_DATA_TRUE$Co_occurrence==1),]
BY_plant_occur=split(DB_DATA_TRUE_co_occur,DB_DATA_TRUE_co_occur$sampMatbased)
plants_occur=names(BY_plant_occur)


############################################################################
# save data for CO-OCCURENCE analisys

saveRDS(DB_DATA_TRUE_co_occur,"data/DB_DATA_TRUE_co_occur.rds")
saveRDS(BY_plant_occur,"data/BY_plant_occur.rds")

##########################################################################
# working on OCCURENCES

# reload data

data_occ=readRDS("data/data_occ.rds")
DB_DATA_TRUE=readRDS("data/DB_DATA_TRUE.rds")
temp_plan_raw=readRDS("data/template.rds")
mycorank=readRDS("data/mycorank.rds")
plantrank=readRDS("data/plantrank.rds")

plants=unique(DB_DATA_TRUE$sampMatbased) # 13 
mycotoxins=unique(DB_DATA_TRUE$paramType) # mycotoxins

######################################################################################## 
# create list for OCCURENCE data analisys

res_tot=list()
res_data=list()
res_names=list()
res_pooled=list()

z=1
for (i in plants)  {
   temp_plant=DB_DATA_TRUE[DB_DATA_TRUE$sampMatbased==i,]
   for (j in mycotoxins) { temp_plant_myco=temp_plant[temp_plant$paramType==j,]
                         
                           id_data=unique(c(which(temp_plant_myco$meanTot>0),which(temp_plant_myco$Concentration>0)))
                           
                           datavalid=as.numeric(c(temp_plant_myco$meanTot[which(temp_plant_myco$meanTot>0)],temp_plant_myco$Concentration[which(temp_plant_myco$Concentration>0)]))
                           
                         if (length(id_data)>5)  # criteria:  minimum  5 data valid

                         {  res_tot[[z]]=temp_plant_myco;
                            res_data[[z]]=temp_plant_myco[id_data,];
                            res_names[[z]]=paste(i,j,sep="_")
                            res_pooled[[z]]=datavalid
                            z=z+1
                         }
                         }
}

# update data savings

saveRDS(res_tot,"data/res_tot_min5.rds")
saveRDS(res_data,"data/res_data_min5.rds")
saveRDS(res_names,"data/res_names_min5.rds")
saveRDS(res_pooled,"data/res_pooled_min5.rds")


#########################################################################################################


df_stats_valid=data.frame(plant_mycotoxins=unlist(res_names),do.call("rbind",lapply(res_pooled,summary_large)),do.call("rbind",lapply(res_tot,function(x) {paste(as.character(unique(x$Ref_unique)),collapse = ", ")})))

file.remove("data/mycotoxins_stats_valid.xls")

XLConnect::writeWorksheetToFile("data/mycotoxins_stats_valid.xls",df_stats_valid,"tab_stats")


#########################################################################################################

plant_myco_db=read.csv(textConnection(gsub("_",",",unlist(res_names))),header=F)
names(plant_myco_db)=c("plants","mycotoxins")

plant_myco_db$ndata_valid=unlist(lapply(res_pooled,length))
plant_myco_db$nrecords=unlist(lapply(res_tot,nrow))

summ_dimensions=data.frame(rbind(
  summ_nrecords=summary_large(plant_myco_db$nrecords),
  summ_nvalid=summary_large(plant_myco_db$ndata_valid),
  summ_pvalid=summary_large(plant_myco_db$ndata_valid/plant_myco_db$nrecords),
  summ_cv=summary_large(unlist(lapply(res_pooled,function(x){ sd(as.numeric(x), na.rm=TRUE)/mean(as.numeric(x), na.rm=TRUE)}))),
  summ_pnormality=summary_large(unlist(lapply(res_pooled,function(x) ifelse(shapiro.test(as.numeric(x))$p.value>0.05,1,0)))),
  summ_psampSize=summary_large(unlist(lapply(res_tot,function(x) mean(x$sampSize,na.rm=T)))),
  summ_p_bibintensity=summary_large(unlist(lapply(res_tot,function(x) length(unique((x$Ref_unique)))))),
  summ_p_havebounds=summary_large(unlist(lapply(res_tot,function(x) mean(x$havebounds))))))
summ_dimensions$par=row.names(summ_dimensions)
summ_dimensions=summ_dimensions[c("par",names(summ_dimensions)[2:16])]
saveRDS(summ_dimensions,"data/summ_dimensions.rds")
file.remove("data/DATA_summary_by_dimensions.xls")
XLConnect::writeWorksheetToFile("data/DATA_summary_by_dimensions.xls",summ_dimensions,"table par dimensions")


#########################################################################################################
# WORKING on OCCURENCE BIBLIOMETRICS
# reload data

res_tot=readRDS("data/res_tot_min5.rds")
res_data=readRDS("data/res_data_min5.rds")
res_names=readRDS("data/res_names_min5.rds")
res_pooled=readRDS("data/res_pooled_min5.rds")

plant_myco_db=read.csv(textConnection(gsub("_",",",unlist(res_names))),header=F)
names(plant_myco_db)=c("plants","mycotoxins")

plant_myco_db$ndata_valid=unlist(lapply(res_pooled,length))
plant_myco_db$nrecords=unlist(lapply(res_tot,nrow))


plant_myco_db$score_numerosity=ifelse(plant_myco_db$ndata_valid>25,1,plant_myco_db$ndata_valid/25)
plant_myco_db$score_validity=plant_myco_db$ndata_valid/plant_myco_db$nrecords
plant_myco_db$p_cv=range01(scale(unlist(lapply(res_pooled,function(x){ sd(as.numeric(x), na.rm=TRUE)/mean(as.numeric(x), na.rm=TRUE)})),center=F))
plant_myco_db$p_normality=unlist(lapply(res_pooled,function(x) ifelse(shapiro.test(as.numeric(x))$p.value>0.05,1,0)))
plant_myco_db$p_sampsize=range01(unlist(lapply(res_tot,function(x) mean(x$norm_sampSize,na.rm=T))))
plant_myco_db$p_agepaper=range01(scale(unlist(lapply(res_tot,function(x) mean(x$norm_agepaper,na.rm=T))),center=F))
plant_myco_db$p_bibintensity=range01(scale(unlist(lapply(res_tot,function(x) length(unique((x$Ref))))),center = F))
real_value_bibintensity=unlist(lapply(res_tot,function(x) length(unique((x$Ref)))))
min_max_bibintensity=c(min(real_value_bibintensity),max(real_value_bibintensity))
file.remove("data/values_bibintensity_5_updated.xls")
XLConnect::writeWorksheetToFile("data/values_bibintensity_5_updated.xls",cbind(plant_myco_db[,1:2],real_value_bibintensity),"bib_intens_more5")

plant_myco_db$p_havebounds=unlist(lapply(res_tot,function(x) mean(x$havebounds)))

plant_myco_db$scoreGEN=plant_myco_db$score_numerosity+
                       plant_myco_db$score_validity+
                       plant_myco_db$p_cv+
                       plant_myco_db$p_havebounds+
                       plant_myco_db$p_sampsize+
                       plant_myco_db$p_agepaper+
                       plant_myco_db$p_bibintensity



saveRDS(plant_myco_db,"data/plant_myco_db.rds")
file.remove("data/mycotoxins_score_5_updated.xls")
XLConnect::writeWorksheetToFile("data/mycotoxins_score_5_updated.xls",plant_myco_db,"data")



list_reliable_myco5=data.frame(name=unlist(res_names),ndata=unlist(lapply(res_data, nrow)))
file.remove("data/mycotoxins_info.xls")

XLConnect::writeWorksheetToFile("data/mycotoxins_info.xls",list_reliable_myco5,"combo_n5")



###################################################################################################################
# Working on CO-OCCURENCE

DB_DATA_TRUE_co_occur=readRDS("data/DB_DATA_TRUE_co_occur.rds")

names(DB_DATA_TRUE_co_occur)


BY_plant_occur=readRDS("data/BY_plant_occur.rds")
plants_occur=names(BY_plant_occur) # 13

table_plants_co_occur=reshape::melt(unlist(lapply(BY_plant_occur,nrow)))



res_occ_plants=list()

for ( i in 1:length(plants_occur)) {
  
  res_occ_plants[[i]]=split(BY_plant_occur[[i]],BY_plant_occur[[i]]$Ref)
  
}
names(res_occ_plants)=plants_occur
saveRDS(res_occ_plants,"data/res_occ_plants_fin.rds")

########################################################################################################
# WORKING ON DATA LIST FOR ALL CO-OCCURENCE DATA

res_tot_occ=list()
res_data_occ=list()
res_names_occ=list()
res_pooled_occ=list()

z=1

for ( i in 1:length(plants_occur)) {
  
     temp_plant=res_occ_plants[[i]]
     
     for (j in seq_along(temp_plant)) { 
       
    temp_plant_myco=temp_plant[[j]]
    id_data=unique(c(which(temp_plant_myco$meanTot>0),which(temp_plant_myco$Concentration>0)))
    datavalid=as.numeric(c(temp_plant_myco$meanTot[which(temp_plant_myco$meanTot>0)],temp_plant_myco$Concentration[which(temp_plant_myco$Concentration>0)]))
    datavalid = datavalid[!is.na(datavalid)]
    
    if (length(datavalid)>1) 
    {  res_tot_occ[[z]]=temp_plant_myco;
       res_data_occ[[z]]=temp_plant_myco[id_data,];
       res_names_occ[[z]]=paste(names(res_occ_plants)[i],paste(unique(res_occ_plants[[i]][[j]]$paramType),collapse = ":"),sep=";")
       res_pooled_occ[[z]]=datavalid
       z=z+1
    }
    }
  }

saveRDS(res_tot,"data/res_tot_occ.rds")
saveRDS(res_data,"data/res_data_occ.rds")
saveRDS(res_names,"data/res_names_occ.rds")
saveRDS(res_pooled,"data/res_pooled_occ.rds")

########################################################################################################
# WORKING on co-OCCURENCE BIBLIOMETRICS


plant_myco_db_occ=read.csv(textConnection(unlist(res_names_occ)),sep=";",header=F)
names(plant_myco_db_occ)=c("plants","mycotoxins")

plant_myco_db_occ$ndata_valid=unlist(lapply(res_pooled_occ,length))
plant_myco_db_occ$nrecords=unlist(lapply(res_tot_occ,nrow))


plant_myco_db_occ$score_numerosity=ifelse(plant_myco_db_occ$ndata_valid>25,1,plant_myco_db_occ$ndata_valid/25)
plant_myco_db_occ$score_validity=plant_myco_db_occ$ndata_valid/plant_myco_db_occ$nrecords
plant_myco_db_occ$p_sampsize=range01(unlist(lapply(res_tot_occ,function(x) mean(x$norm_sampSize,na.rm=T))))
plant_myco_db_occ$p_agepaper=range01(scale(unlist(lapply(res_tot_occ,function(x) mean(x$norm_agepaper,na.rm=T))),center=F))
plant_myco_db_occ$p_agepaper=ifelse(is.na(plant_myco_db_occ$p_agepaper),0,plant_myco_db_occ$p_agepaper)
plant_myco_db_occ$p_havebounds=unlist(lapply(res_tot_occ,function(x) mean(x$havebounds)))

plant_myco_db_occ$scoreGEN=plant_myco_db_occ$score_numerosity+
                           plant_myco_db_occ$score_validity+
                           plant_myco_db_occ$p_havebounds+
                            plant_myco_db_occ$p_sampsize+
                           plant_myco_db_occ$p_agepaper

file.remove("data/mycotoxins_score_co_occurence.xls")

XLConnect::writeWorksheetToFile("data/mycotoxins_score_co_occurence.xls",plant_myco_db_occ,"data_co_occ")

saveRDS(plant_myco_db_occ,"data/plant_myco_db_occ_tot.rds")



#################################################################################################
# WORKING on last requests

DB_occ_paola_maize=DB_DATA_TRUE_co_occur[c(which(DB_DATA_TRUE_co_occur$sampMatbased=="maize")),]
DB_occ_paola_wheat=DB_DATA_TRUE_co_occur[c(which(DB_DATA_TRUE_co_occur$sampMatbased=="wheat")),]
DB_occ_paola_oat=DB_DATA_TRUE_co_occur[c(which(DB_DATA_TRUE_co_occur$sampMatbased=="oat")),]
DB_occ_paola_barley=DB_DATA_TRUE_co_occur[c(which(DB_DATA_TRUE_co_occur$sampMatbased=="barley")),]

# merging DON & FB

DB_occ_paola_maize$paramType[grep("DON",DB_occ_paola_maize$paramType)]="DONs"
DB_occ_paola_maize$paramType[grep("FB",DB_occ_paola_maize$paramType)]="FBs"
DB_occ_paola_wheat$paramType[grep("DON",DB_occ_paola_wheat$paramType)]="DONs"
DB_occ_paola_oat$paramType[grep("DON",DB_occ_paola_oat$paramType)]="DONs"
DB_occ_paola_barley$paramType[grep("DON",DB_occ_paola_barley$paramType)]="DONs"

###########################################################################################################################

DB_occ_paola_maize=DB_occ_paola_maize[c(which(DB_occ_paola_maize$paramType=="DONs"),
                                        which(DB_occ_paola_maize$paramType=="FBs"),
                                        which(DB_occ_paola_maize$paramType=="AF")),]

DB_occ_paola_wheat=DB_occ_paola_wheat[c(which(DB_occ_paola_wheat$paramType=="DONs"),
                                        which(DB_occ_paola_wheat$paramType=="ZEN"),
                                        which(DB_occ_paola_wheat$paramType=="NIV"),
                                        which(DB_occ_paola_wheat$paramType=="T2+HT2")),]

DB_occ_paola_barley=DB_occ_paola_barley[c(which(DB_occ_paola_barley$paramType=="DONs"),
                                          which(DB_occ_paola_barley$paramType=="ZEN"),
                                          which(DB_occ_paola_wheat$paramType=="NIV"),
                                        which(DB_occ_paola_wheat$paramType=="T2+HT2")),]

DB_occ_paola_oat=DB_occ_paola_oat[c(which(DB_occ_paola_oat$paramType=="DONs"),
                                    which(DB_occ_paola_barley$paramType=="ZEN"),
                                    which(DB_occ_paola_oat$paramType=="T2+HT2"),
                                    which(DB_occ_paola_oat$paramType=="NIV")),]



BY_plant_occur=list(maize=DB_occ_paola_maize,
                    wheat=DB_occ_paola_wheat,
                    oat=DB_occ_paola_oat,
                    barley=DB_occ_paola_barley)


plants_occur=names(BY_plant_occur) # 4

table_plants_co_occur=reshape::melt(unlist(lapply(BY_plant_occur,nrow)))

file.remove("data/table_co_occur_fin.xls")
XLConnect::writeWorksheetToFile("data/table_co_occur_fin.xls",table_plants_co_occur,"table_plants_occur")


########################################################################################################
# WORKING ON DATA LIST FOR CO-OCCURENCE selected data

res_tot_occ=list()
res_data_occ=list()
res_names_occ=list()
res_pooled_occ=list()
res_stats_occ=list()
res_prop_occ=list()

z=1

for ( i in 1:length(plants_occur)) {
  
  temp_plant=res_occ_plants[[i]]
  
  for (j in seq_along(temp_plant)) { 
    
    temp_plant_myco=temp_plant[[j]]
    id_data=unique(c(which(temp_plant_myco$meanTot>0),which(temp_plant_myco$Concentration>0)))
    datavalid=as.numeric(c(temp_plant_myco$meanTot[which(temp_plant_myco$meanTot>0)],temp_plant_myco$Concentration[which(temp_plant_myco$Concentration>0)]))
    paramvalid=c(temp_plant_myco$paramType[which(temp_plant_myco$meanTot>0)],temp_plant_myco$paramType[which(temp_plant_myco$Concentration>0)])
    
     paramvalid= paramvalid[which(!is.na(datavalid))]
     datavalid = datavalid[!is.na(datavalid)]
     
    if (length(datavalid)>0) 
    {  res_tot_occ[[z]]=temp_plant_myco;
       res_data_occ[[z]]=temp_plant_myco[id_data,];
       res_names_occ[[z]]=paste(names(res_occ_plants)[i],paste(unique(res_occ_plants[[i]][[j]]$paramType),collapse = ":"),sep=";")
       res_pooled_occ[[z]]=datavalid
       names(res_pooled_occ[[z]])=paramvalid
       res_stats_occ[[z]]=c(tapply(datavalid, paste0("Mean_",paramvalid), mean),
                            tapply(datavalid, paste0("Sd_",paramvalid), sd),
                            tapply(datavalid, paste0("Max_",paramvalid), max))
       res_prop_occ[[z]]=tapply(datavalid, paste0("Count_",paramvalid),length)
       
    z=z+1
    }
  }
}


plant_myco_db_occ=read.csv(textConnection(unlist(res_names_occ)),sep=";",header=F)

names(plant_myco_db_occ)=c("plants","mycotoxins")

plant_myco_db_occ$ndata_valid=unlist(lapply(res_pooled_occ,length))
plant_myco_db_occ$nrecords=unlist(lapply(res_tot_occ,nrow))



plant_myco_db_occ$Ref=unlist(lapply(res_tot_occ,function(x)x$Ref[1]))

plant_myco_db_occ$mycotoxins=as.character(plant_myco_db_occ$mycotoxins)
plant_myco_db_occ$plants=as.character(plant_myco_db_occ$plants)
id_valid=grep(":",plant_myco_db_occ$mycotoxins)

##########################################################################################
# here co-occurence selection id_valid variables

plant_myco_db_occ=plant_myco_db_occ[id_valid,]
res=res_stats_occ[id_valid]
res_prop=res_prop_occ[id_valid]

#####################################################################################################################################################

saveRDS(res_tot_occ[id_valid],"data/res_tot_occ_sel.rds")
saveRDS(res_data_occ[id_valid],"data/res_data_occ_sel.rds")
saveRDS(res_names_occ[id_valid],"data/res_names_occ_sel.rds")
saveRDS(res_pooled_occ[id_valid],"data/res_pooled_occ_sel.rds")

#####################################################################################################################################################
# export database dump

res_pooled_occ=res_pooled_occ[id_valid]

names(res_pooled_occ)=sapply(1:nrow(plant_myco_db_occ),function(x) (paste(plant_myco_db_occ[x,1],plant_myco_db_occ[x,2],collapse=":")))

file.remove("dump/file_data_co_occurence.txt")

sink(file="dump/file_data_co_occurence.txt")
res_pooled_occ
sink(type = "message")
sink()

res_pool_data=list()
for (i in 1:nrow(plant_myco_db_occ)){
  res_pool_data[[i]]=paste(res_pooled_occ[[i]],collapse=";")
  names()
}

resls=list()
for (i in 1:nrow(plant_myco_db_occ)){
  resls[[i]]=paste0(paste(names(res[[i]]),collapse=";"),";",paste(names(res_prop[[i]]),collapse=";"),";",paste(res[[i]],collapse=";"),";",paste(res_prop[[i]],collapse=";"),collapse = "")
  
}



plant_myco_db_occ$data=do.call("rbind",resls)


#########################################################################################################################################


file.remove("data/mycotoxins_co_occurence_ONEtable_fin.xls")

XLConnect::writeWorksheetToFile("data/mycotoxins_co_occurence_ONEtable_fin.xls",plant_myco_db_occ,"table data co-occurence")


 
file.remove("data/mycotoxins_co_occurence_MULTItable_fin.xls")

plant_myco_db_occ_ls=split(plant_myco_db_occ,as.factor(paste(plant_myco_db_occ$plants,plant_myco_db_occ$mycotoxins)))



for ( i in seq_along(plant_myco_db_occ_ls)) {
                                    XLConnect::writeWorksheetToFile("data/mycotoxins_co_occurence_MULTItable_fin.xls",plant_myco_db_occ_ls[[i]],gsub(":"," ",substr(names(plant_myco_db_occ_ls)[i],1,15)))

}

# AFTER CREATE MANUALLY ANOTHER EXCEL FILE ARRANGED  mycotoxins_co_occurence_fin_arranged.xls

saveRDS(plant_myco_db_occ,"data/plant_myco_db_occ_fin.rds")





###########################################################################################################################################################


#  Reference
#  https://dabblingwithdata.wordpress.com/2018/01/02/my-favourite-r-package-for-summarising-data/
#  https://rawgraphs.io/learning/how-to-make-a-beeswarm-plot/
#  Rlist packages https://renkun-ken.github.io/rlist/
