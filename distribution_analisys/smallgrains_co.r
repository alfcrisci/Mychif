########################################################################
library(readxl)
library(XLConnect)
library(ggplot2)
library(dplyr)
library(scales)
library(cartography)


########################################################################

setwd("/home/alf/Scrivania/lav_michyf/distrib")
source("aux_distrib.r")

#########################################################################################################

SMALLGRAINS=readRDS("SMALLGRAINS.rds")

query_vec=unlist(strsplit("almond;barley;buckwheat;caraway;cereals;durum wheat;millet;oat;peanuts;rye;soft wheat;soy;spelt;triticale;walnuts;wheat",";"))


#############################################################################################################

SMALLGRAINS=SMALLGRAINS[which(!is.na(SMALLGRAINS$sampMatType)),]
SMALLGRAINS=SMALLGRAINS[-which(SMALLGRAINS$sampMatType=="food;feed"),]
SMALLGRAINS=SMALLGRAINS[-which(SMALLGRAINS$sampMatType=="food, feed;seed"),]

#############################################################################################################


unique(SMALLGRAINS$sampMatType)
res_co_data_mat=list()

zz=1

for ( query in query_vec) {
  
SMALLGRAINS_co=SMALLGRAINS[which(SMALLGRAINS$sampMatbased==query),]
names(SMALLGRAINS_co)[42]="Co_occurrence"
SMALLGRAINS_co=as.data.frame(SMALLGRAINS_co[which(SMALLGRAINS_co$Co_occurrence==1),])
if (nrow(SMALLGRAINS_co)==0) {next}

SMALLGRAINS_co$paramType=gsub("/-","_",tolower(SMALLGRAINS_co$paramType))
SMALLGRAINS_co$paramType=tolower(SMALLGRAINS_co$paramType)
SMALLGRAINS_co$sampMatType=tolower(SMALLGRAINS_co$sampMatType)
SMALLGRAINS_co$sampMatbased=tolower(SMALLGRAINS_co$sampMatbased)
SMALLGRAINS_co$sampMatCode=tolower(SMALLGRAINS_co$sampMatCode)
U_ref_co=unique(SMALLGRAINS_co$Ref)
U_type_co=unique(SMALLGRAINS_co$sampMatType)

#############################################################################################
res_co_data=list()
res_co_tab=list()

z=1
for ( i in seq_along(U_ref_co)) {
  for ( j in seq_along(U_type_co)){
      
    tab=info_extract_co_occur(SMALLGRAINS_co,U_ref_co[i],U_type_co[j])
    U_sampsize=unique(tab$mat$sampSize)
    
    for ( ind in seq_along(U_sampsize)) {
    res_co_tab[[z]]=tab$ls_mat
    matsize=tab$mat[which(tab$mat$sampSize==U_sampsize[ind]),]  
    res_co_data[[z]]=data.frame(sampMatbased=query,
                                sampSize=U_sampsize[ind],
                                paramType=paste(unique(matsize$paramType),collapse = "+"),
                                sampCountry=paste(unique(matsize$sampCountry),collapse=";"),
                                sampCountryorigin=paste(unique(matsize$sampCountry),collapse=";"),
                                Records=nrow(matsize),
                                Ref=U_ref_co[i])
                                
     z=z+1                           
    
    } # ind Usampsize
   } #i
}
data_tab=data.frame(do.call("rbind",res_co_data))
data_tab=data_tab[which(data_tab$paramType != 0),]
res_co_data_mat[[zz]]=data_tab
zz=zz+1
}


data_tab_mat=data.frame(do.call("rbind",res_co_data_mat))

file.remove(paste0("Tabella_co_smallgrain.xls"))
XLConnect::writeWorksheetToFile(paste0("Tabella_co_smallgrain.xls"),data_tab_mat,query)



#########################################################################################################
# References

# https://stats.stackexchange.com/questions/153725/plotting-a-ecdf-in-r-and-overlay-cdf
# https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist
# http://rstudio-pubs-static.s3.amazonaws.com/5554_25ed8319163a4df6bd644e68c6fd4b21.html
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# https://www.r-bloggers.com/analysing-longitudinal-data-multilevel-growth-models-i/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
