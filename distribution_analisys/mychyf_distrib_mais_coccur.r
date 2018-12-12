########################################################################

library(readxl)
library(XLConnect)
library(ggplot2)
library(dplyr)
library(scales)
library(fitdistrplus)

setwd("/home/alf/Scrivania/lav_michyf/distrib")
source("aux_distrib.r")

#########################################################################################################

MAIZE=readRDS("MAIZE.rds")
query="maize"
MAIZE_co=MAIZE[which(MAIZE$sampMatbased==query),]
MAIZE_co=as.data.frame(MAIZE[which(MAIZE$`Co-occurrence`==1),])

####

MAIZE_co$paramType=tolower(MAIZE_co$paramType)
MAIZE_co$sampMatType=tolower(MAIZE_co$sampMatType)
MAIZE_co$sampMatbased=tolower(MAIZE_co$sampMatbased)
MAIZE_co$sampMatCode=tolower(MAIZE_co$sampMatCode)

U_ref_co=unique(MAIZE_co$Ref)
U_type_co=unique(MAIZE_co$sampMatType)

#############################################################################################

res_co_data=list()
res_co_tab=list()

z=1
for ( i in seq_along(U_ref_co)) {
  for ( j in seq_along(U_type_co)){
      
    tab=info_extract_co_occur(MAIZE_co,U_ref_co[i],U_type_co[j])
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
    
    }

  }
  }
  
#########################################################################################################


# tabella Maize+Feed
# Paramtype	Sampsize	sampCountry	sampCountryOrigin	Records number
# DON+HT2+ADON	8	ES	ES

#########################################################################################################


data_tab=data.frame(do.call("rbind",res_co_data))
data_tab=data_tab[which(data_tab$paramType != 0),]
file.remove(paste0("Tabella_co_",query,".xls"))
XLConnect::writeWorksheetToFile(paste0("Tabella_co_",query,".xls"),data_tab,query)



#########################################################################################################
# References

# https://stats.stackexchange.com/questions/153725/plotting-a-ecdf-in-r-and-overlay-cdf
# https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist
# http://rstudio-pubs-static.s3.amazonaws.com/5554_25ed8319163a4df6bd644e68c6fd4b21.html
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# https://www.r-bloggers.com/analysing-longitudinal-data-multilevel-growth-models-i/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
