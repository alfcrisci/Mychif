########################################################################

library(XLConnect)
library(ggplot2)
library(dplyr)
library(scales)
library(fitdistrplus)

########################################################################

setwd("/home/alf/Scrivania/lav_michyf/distrib")

source("aux_distrib.r")

##################################################################################

SMALLGRAINS=readRDS("matrix_data/SMALLGRAINS.rds")

##################################################################################

query_vec=unlist(strsplit("almond;barley;buckwheat;caraway;cereals;durum wheat;millet;oat;peanuts;rye;soft wheat;soy;spelt;triticale;walnuts;wheat",";"))


##################################################################################

dir.create(paste0("SMALLGRAINS_plots"))

for ( query in query_vec) {
dir.create(paste0("SMALLGRAINS_plots/",query))
dir.create(paste0("SMALLGRAINS_plots/",query,"/plot_simple"))
dir.create(paste0("SMALLGRAINS_plots/",query,"/plot_full"))

SMALLGRAINSq=SMALLGRAINS[which(SMALLGRAINS$sampMatbased==query),]
SMALLGRAINSq$paramType=gsub("/-","_",SMALLGRAINSq$paramType)

#SMALLGRAINSq$sampMatType=tolower(SMALLGRAINSq$sampMatType)
#SMALLGRAINSq$sampMatbased=tolower(SMALLGRAINSq$sampMatbased)


U_param=unique(SMALLGRAINSq$paramType)
U_type=c("food","feed")
U_plant=unique(SMALLGRAINSq$sampMatbased)[1]

#########################################################################################################

res_SMALLGRAINS_tab=list()
res_SMALLGRAINS_data=list()
res_SMALLGRAINS_names=list()


#########################################################################################################

z=1
for ( i in seq_along(U_param)) {
      for ( j in seq_along(U_type)){
    
        tab=info_extract_db(SMALLGRAINSq,U_param[i],U_type[j],query)
        
        res_SMALLGRAINS_data[[z]]=tab
        res_SMALLGRAINS_tab[[z]]=tab$ls_tab
        res_SMALLGRAINS_names[[z]]=paste0(U_param[i],";",U_type[j],";",paste(query,"on",U_type[j]))
        z=z+1
        } 
 }
      SMALLGRAINSq_tab=data.frame(names=do.call("rbind",res_SMALLGRAINS_names),do.call("rbind",res_SMALLGRAINS_tab))
      file.remove(paste0("Tabella_SMALLGRAINS_",query,".xls"))
      XLConnect::writeWorksheetToFile(paste0("Tabella_SMALLGRAINS_",query,".xls"),SMALLGRAINSq_tab,query)
        

z=1
      for ( i in seq_along(U_param)) {
        for ( j in seq_along(U_type)){
          
          tab=info_extract_db(SMALLGRAINSq,U_param[i],U_type[j],query)
          
           if (tab$N_conc > 0) {
                       ggplot(tab$mat_op, aes(Concentration_op))+
                              stat_ecdf(geom = "point")+
                              ylab("Density")+
                              ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = query)
                              ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_simple_SMALLGRAINS_",query,".jpg"),dpi=300)
    
                       smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),subtitleg=query,labelx = "Concentration_op")
                       ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_conc_full_SMALLGRAINS_",query,".jpg"),dpi=300)
                       if(length(na.omit(tab$mat_op$Concentration_op)) > 4) 
                       {jpeg(filename = paste0(U_param[i],"_on_",U_type[j],"_culley_graph_mtot_",query,".jpg"),width = 800, height = 800)
                         try(print(descdist(as.numeric(na.omit(tab$mat_op$Concentration_op)),boot=500)))
                         dev.off()
                       }
                    }
  if (tab$N_mtot > 0) {
                       ggplot(tab$mat_op, aes(meanTot_op))+
                              stat_ecdf(geom = "point")+
                              ylab("Density")+
                              ggtitle(paste(U_param[i],"on",U_type[j]),subtitle = query)
    
                       ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_mtot_simple_SMALLGRAINS_",query,".jpg"),dpi=300)
    
                       smooth_ecd_gray(as.data.frame(tab$mat_op),titleg = paste(U_param[i],"on",U_type[j]),subtitleg=query,labelx = "meanTot_op")
                       ggsave(filename = paste0(U_param[i],"_on_",U_type[j],"_mtot_full_SMALLGRAINS_",query,".jpg"),dpi=300)
                       if(length(na.omit(tab$mat_op$meanTot_op)) > 4) 
                       {jpeg(filename = paste0(U_param[i],"_on_",U_type[j],"_culley_graph_mtot_",query,".jpg"),width = 800, height = 800)
                         try(print(descdist(as.numeric(na.omit(tab$mat_op$meanTot_op)),boot=500)))
                         dev.off()
                       }
                    }
    z=z+1;
    }
}
}
################################################################################################################################

dir.create("SMALLGRAINS_plots/culley_graphs")
system(paste0("mv *culley*.jpg ","SMALLGRAINS_plots/culley_graphs"))

for ( query in query_vec) {
  system(paste0("mv *simple*",query,"*.jpg SMALLGRAINS_plots/",query,"/plot_simple"))
}

for ( query in query_vec) {
  system(paste0("mv *full*",query,"*.jpg SMALLGRAINS_plots/",query,"/plot_full"))
}


#########################################################################################################
# References

# https://stats.stackexchange.com/questions/153725/plotting-a-ecdf-in-r-and-overlay-cdf
# https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist
# http://rstudio-pubs-static.s3.amazonaws.com/5554_25ed8319163a4df6bd644e68c6fd4b21.html
# https://cran.r-project.org/web/packages/codyn/vignettes/Temporal_Diversity_Indices.html
# https://www.r-bloggers.com/analysing-longitudinal-data-multilevel-growth-models-i/
# https://rstudio-pubs-static.s3.amazonaws.com/228019_f0c39e05758a4a51b435b19dbd321c23.html
