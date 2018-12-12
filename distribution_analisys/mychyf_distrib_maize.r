
library(XLConnect)
library(ggplot2)
library(dplyr)
library(scales)
library(fitdistrplus)
library(mc2d)

#########################################################################################################

setwd("/home/alf/Scrivania/lav_michyf/distrib")
source("aux_distrib.r")


#########################################################################################################
MAIZE=readRDS("matrix_data/MAIZE.rds")

#MAIZE$paramType=tolower(MAIZE$paramType)
#MAIZE$sampMatType=tolower(MAIZE$sampMatType)
#MAIZE$sampMatbased=tolower(MAIZE$sampMatbased)


U_param=unique(MAIZE$paramType)
U_type=unique(MAIZE$sampMatType)
U_plant=unique(MAIZE$sampMatbased)[1]

#########################################################################################################

res_maize_tab=list()
res_maize_num=list()

res_maize_data=list()
res_maize_names=list()


#########################################################################################################
z=1
for ( i in seq_along(U_param)) {
      for ( j in seq_along(U_type[1:2])){
    
        tab=info_extract_db(MAIZE,U_param[i],U_type[j],"maize")
        res_maize_data[[z]]=tab
        res_maize_tab[[z]]=tab$ls_tab
        res_maize_names[[z]]=paste0(U_param[i],";",U_type[j],";","maize")
        z=z+1;
      }
}   

maize_tab=data.frame(names=do.call("rbind",res_maize_names),do.call("rbind",res_maize_tab))
file.remove("Tabella_maize.xls")
XLConnect::writeWorksheetToFile("Tabella_maize.xls",maize_tab,"maize")


z=1
for ( i in seq_along(U_param)) {
  for ( j in seq_along(U_type[1:2])){
        tab=info_extract_db(MAIZE,U_param[i],U_type[j],"maize")
    
         if (tab$N_conc > 0) {
                       ggplot(tab$mat_op, aes(Concentration_op))+
                              stat_ecdf(geom = "point")+
                              ylab("Density")+
                              ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = "maize")
                       ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_simple_maize.jpg"),dpi=300)
    
                       ##############################################################################################################################
                       
                       smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),labelx = "Concentration_op")
                       ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_full_maize.jpg"),dpi=300)
                       ##############################################################################################################################
                       
                       if(length(na.omit(tab$mat_op$Concentration_op)) > 4) 
                          {jpeg(filename = paste0(U_param[i],"_on_",U_type[j],"_culley_graph_mtot.jpg"),width = 800, height = 800)
                           try(print(descdist(as.numeric(na.omit(tab$mat_op$Concentration_op)),boot=500)))
                           dev.off()
                       }
                    }
        
          if (tab$N_mtot > 0) {
                       ggplot(tab$mat_op, aes(meanTot_op))+
                              stat_ecdf(geom = "point")+
                              ylab("Density")+
                              ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = "maize")
    
                      ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_simple_maize.jpg"),dpi=300)
    
                      smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),labelx = "meanTot_op")
                      ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_mtot_full_maize.jpg"),dpi=300)
                    
                      if(length(na.omit(tab$mat_op$meanTot_op)) > 4) 
                            {jpeg(filename = paste0(U_param[i],"_on_",U_type[j],"_culley_graph_mtot.jpg"),width = 800, height = 800)
                                  try(print(descdist(as.numeric(na.omit(tab$mat_op$meanTot_op)),boot=500)))
                            dev.off()
                            }
                    }
    z=z+1;
    }
}

dir.create("myco_mais_plot")
system("mv *.jpg myco_mais_plot")
#########################################################################################################
# References

# https://stats.stackexchange.com/questions/153725/plotting-a-ecdf-in-r-and-overlay-cdf
# https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist
# http://rstudio-pubs-static.s3.amazonaws.com/5554_25ed8319163a4df6bd644e68c6fd4b21.html
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# https://www.r-bloggers.com/analysing-longitudinal-data-multilevel-growth-models-i/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
