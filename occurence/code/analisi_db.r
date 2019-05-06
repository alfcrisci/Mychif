#####################################################################################################################################

setwd("/home/alf/Scrivania/lav_michyf/occurence")

source("load_lib.r")
source("aux_mycosources.r")


#####################################################################################################################################
# data_occ=read_excel("data/all_cereals_FINALE.xls",1)
# data_occ[is.na(data_occ)]=NA
# data_occ=as.data.frame(data_occ)
# saveRDS(data_occ,"data_occ.rds")

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
# read data

data_occ=readRDS("data_occ.rds")

#####################################################################################################################################
# [1] "DB_origin"           "sampCountryorigin"   "sampContinentorigin" "sampContinent"       "paramType"           "sampCountry"        
# [7] "sampRegion"          "sampInfo_latitude"   "sampInfo_longitude"  "sampInfo_altitude"   "sampYear"            "sampYearIncreas"    
# [13] "sampSize"            "sampMethod"          "sampPoint"           "sampPointInfo"       "sampMatType"         "sampMatbased"       
# [19] "sampMatInfo"         "sampMatCode"         "growingSystem"       "anMethRefId"         "resUnit"             "resLOD"             
# [25] "resLODinfo"          "resLOQ"              "resLOQinfo"          "exprResPerc"         "ToTresValUncertSD"   "meanTot"            
# [31] "POSresValUncertSD"   "meanPos"             "median"              "medianinfo"          "min"                 "max"                
# [37] "MAJ_LOD"            "IQRmin"              "IQRmax"              "Concentration"       "DataSampLevel"       "Co_occurrence"      
# [43] "Ref"
#####################################################################################################################################
# GEODATA

data_auth_geo=c("DB_origin","sampContinent","sampCountryorigin","sampContinentorigin","sampCountry","sampRegion")

DB_AUTH_GEO=data_occ[data_auth_geo]
#   to uppercase
DB_AUTH_GEO=as.data.frame(apply(DB_AUTH_GEO,2,toupper))
DB_AUTH_GEO$sampCountryorigin=gsub("UK","GB",DB_AUTH_GEO$sampCountryorigin)
DB_AUTH_GEO$sampCountry=gsub("UK","GB",DB_AUTH_GEO$sampCountry)

#view(dfSummary(DB_AUTH_GEO))

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


########################################################################à
# Risultati

saveRDS(res_cont_orig,"data/res_continent_orig.rds")
saveRDS(res_country_orig,"data/res_country_orig.rds")


########################################################################
# data_stratum

data_stratum=c("paramType","sampMatbased","sampMatType","sampMatInfo","sampMatCode","growingSystem","anMethRefId","sampMethod")
DB_STRATUM=as.data.frame(data_occ[data_stratum])
saveRDS(DB_STRATUM,"data/DB_STRATUM.rds")

#########################################
# data values

data_values=c("sampMatbased","paramType","meanTot","Concentration","meanPos","median","min","max","sampSize","Co_occurrence","resUnit","resLOD","Ref")
DB_DATA_FULL=as.data.frame(data_occ[data_values])


#########################################
# missing data management

id_overlod=which(DB_DATA_FULL$meanTot >=-1)       # 1480
id_overlod2=which(DB_DATA_FULL$Concentration>=-1) # 1488
id_overlod_full=unique(c(id_overlod,id_overlod2)) # 5802
DB_DATA_TRUE=DB_DATA_FULL[id_overlod_full,]       # 

#########################################
# view(DB_DATA_TRUE)
by_plant=split(DB_DATA_TRUE,DB_DATA_TRUE$sampMatbased)

plants=unique(DB_DATA_TRUE$sampMatbased) # 14 
mycotoxins=unique(DB_DATA_TRUE$paramType) # mycotoxins

table_plants_mycotoxin=as.data.frame.array(table(DB_DATA_TRUE$sampMatbased, 
                                                 DB_DATA_TRUE$paramType))
mycorank=sort(colSums(table_plants_mycotoxin), decreasing = T)
plantrank=sort(rowSums(table_plants_mycotoxin),decreasing = T)

aux_info_mycotox=data_occ[id_overlod_full,c("sampYearIncreas","Co_occurrence","Ref")]
aux_info_mycotox$agepaper=2019-as.numeric(gsub("[A-z]","",gsub("\\W+","",aux_info_mycotox$Ref)))
aux_info_mycotox$agepaper[764:782]=8
DB_DATA_TRUE$agepaper=aux_info_mycotox$agepaper

DB_DATA_TRUE$havebounds=ifelse((!is.na(DB_DATA_TRUE$min) & !is.na(DB_DATA_TRUE$max)),1,0)
DB_DATA_TRUE$havemeanpos=ifelse((!is.na(DB_DATA_TRUE$meanPos)),1,0)
DB_DATA_TRUE$havedata=ifelse((!is.na(DB_DATA_TRUE$meanTot) || !is.na(DB_DATA_TRUE$Concentration)),1,0)
DB_DATA_TRUE$norm_agepaper=scale(DB_DATA_TRUE$agepaper,center=F)
DB_DATA_TRUE$norm_sampSize=scale(DB_DATA_TRUE$sampSize,center=F)

# aa=DB_DATA_TRUE$agepaper[which(DB_DATA_TRUE$agepaper<0)]
# aa=2019-as.numeric(substr(aa,1,5))*-1
# DB_DATA_TRUE$agepaper[which(DB_DATA_TRUE$agepaper<0)]=aa

########################################################################################
# sampMatbased paramType
saveRDS(DB_DATA_TRUE,"data/DB_DATA_TRUE.rds")
saveRDS(temp_plan_raw,"data/template.rds")
saveRDS(mycorank,"data/mycorank.rds")
saveRDS(plantrank,"data/plantrank.rds")
DB_DATA_TRUE_co_occur=DB_DATA_TRUE[which(DB_DATA_TRUE$Co_occurrence==1),]
BY_plant_occur=split(DB_DATA_TRUE_co_occur,DB_DATA_TRUE_co_occur$sampMatbased)
plants_occur=names(BY_plant_occur)
saveRDS(DB_DATA_TRUE_co_occur,"data/DB_DATA_TRUE_co_occur.rds")
saveRDS(BY_plant_occur,"data/BY_plant_occur.rds")

######################################################################################## 
######################################################################################## 
######################################################################################## 
data_occ=readRDS("data_occ.rds")
DB_DATA_TRUE=readRDS("data/DB_DATA_TRUE.rds")


# [1] "sampMatbased"  "paramType"     "meanTot"       "Concentration" "meanPos"       "median"        "min"          
# [8] "max"           "sampSize"      "Co_occurrence" "resUnit"       "resLOD"        "Ref"           "agepaper"     
# [15] "havebounds"    "havemeanpos"   "havedata"      "norm_agepaper" "norm_sampSize"
 


temp_plan_raw=readRDS("data/template.rds")
mycorank=readRDS("data/mycorank.rds")
plantrank=readRDS("data/plantrank.rds")
res_continent_orig=readRDS("data/res_continent_orig.rds")
res_country_orig=readRDS("data/res_country_orig.rds")

plants=unique(DB_DATA_TRUE$sampMatbased) # 14 
mycotoxins=unique(DB_DATA_TRUE$paramType) # mycotoxins



######################################################################################## 


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
                           
                         if (length(id_data)>5) 
                         {  res_tot[[z]]=temp_plant_myco;
                            res_data[[z]]=temp_plant_myco[id_data,];
                            res_names[[z]]=paste(i,j,sep="_")
                            res_pooled[[z]]=datavalid
                            z=z+1
                         }
                         }
}

#########################################################################################################
# minimum 5 data valid

saveRDS(res_tot,"data/res_tot_min5.rds")
saveRDS(res_data,"data/res_data_min5.rds")
saveRDS(res_names,"data/res_names_min5.rds")
saveRDS(res_pooled,"data/res_pooled_min5.rds")


#########################################################################################################
res_tot=readRDS("data/res_tot_min5.rds")
res_data=readRDS("data/res_data_min5.rds")
res_names=readRDS("data/res_names_min5.rds")
res_pooled=readRDS("data/res_pooled_min5.rds")

#########################################################################################################

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
plant_myco_db$p_havebounds=unlist(lapply(res_tot,function(x) mean(x$havebounds)))

plant_myco_db$scoreGEN=plant_myco_db$score_numerosity+
                       plant_myco_db$score_validity+
                       plant_myco_db$p_cv+
                       plant_myco_db$p_havebounds+
                       plant_myco_db$p_sampsize+
                       plant_myco_db$p_agepaper+
                       plant_myco_db$p_bibintensity

file.remove("data/mycotoxins_score_5_updated.xls")
XLConnect::writeWorksheetToFile("data/mycotoxins_score_5_updated.xls",plant_myco_db,"data")
saveRDS(plant_myco_db,"data/plant_myco_db.rds")

#########################################################################################################

list_reliable_myco4=data.frame(name=unlist(res_names),ndata=unlist(lapply(res_data, nrow)))
list_reliable_myco5=data.frame(name=unlist(res_names),ndata=unlist(lapply(res_data, nrow)))
list_reliable_myco6=data.frame(name=unlist(res_names),ndata=unlist(lapply(res_data, nrow)))

XLConnect::writeWorksheetToFile("data/mycotoxins_COMBO.xls",list_reliable_myco4,"combon4")
XLConnect::writeWorksheetToFile("data/mycotoxins_COMBO.xls",list_reliable_myco5,"combon5")
XLConnect::writeWorksheetToFile("data/mycotoxins_COMBO.xls",list_reliable_myco6,"combon6")




#########################################
# info on data values

info_data_values=c("sampSize","resUnit","resLOD","resLODinfo","exprResPerc","ToTresValUncertSD","medianinfo","POSresValUncertSD","MAJ_LOD","Co_occurrence","Ref")

DB_INFODATA=as.data.frame(data_occ[info_data_values])



###########################################################################

DB_DATA_TRUE_co_occur=readRDS("data/DB_DATA_TRUE_co_occur.rds")
BY_plant_occur=readRDS("data/BY_plant_occur.rds")

plants_occur=names(BY_plant_occur)# 13
table_plants_occur=reshape::melt(unlist(lapply(BY_plant_occur,nrow)))

res_occ_plants=list()



for ( i in 1:length(plants_occur)) {
         
  res_occ_plants[[i]]=split(BY_plant_occur[[i]],BY_plant_occur[[i]]$Ref)

}
names(res_occ_plants)=plants_occur
saveRDS(res_occ_plants,"data/res_occ_plants.rds")
########################################################################################################

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
    
    if (length(id_data)>5) 
    {  res_tot_occ[[z]]=temp_plant_myco;
       res_data_occ[[z]]=temp_plant_myco[id_data,];
       res_names_occ[[z]]=paste(names(res_occ_plants)[i],paste(unique(res_occ_plants[[i]][[j]]$paramType),collapse = ":"),sep=";")
       res_pooled_occ[[z]]=datavalid
       z=z+1
    }
    }
  }


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

file.remove("data/mycotoxins_score_5_occurence.xls")
XLConnect::writeWorksheetToFile("data/mycotoxins_score_5_occurence.xls",plant_myco_db_occ,"dataocc")
saveRDS(plant_myco_db_occ,"data/plant_myco_db_occ.rds")



###########################################################################################################################################################


#  Reference
#  https://dabblingwithdata.wordpress.com/2018/01/02/my-favourite-r-package-for-summarising-data/
#  https://rawgraphs.io/learning/how-to-make-a-beeswarm-plot/
#  Rlist packages https://renkun-ken.github.io/rlist/
