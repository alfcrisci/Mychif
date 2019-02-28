#############################################################################################################################
# Setup working directory

setwd("/home/alf/Scrivania/lav_michyf/occurence")
##############################################################################################################################
# Load libraries & functions . Check if packages was installed.

library(XLConnect)
source("load_lib.r")
source("aux_mycosources.r")

##############################################################################################################################
# Read data

res_tot=readRDS("data/res_tot_min5.rds")
res_data=readRDS("data/res_data_min5.rds")
res_names=readRDS("data/res_names_min5.rds")
res_pooled=readRDS("data/res_pooled_min5.rds")
plant_myco_db=readRDS("hdata/plant_myco_db.rds")
plant_myco_db_occ=readRDS("data/plant_myco_db_occ.rds")

#########################################################################################################
# filtering by general final bibliometrics score
tresh=2.5
index=which(plant_myco_db$scoreGEN>tresh)
plant_myco_db_sel=plant_myco_db[index,]
names(plant_myco_db_sel)

# create color palette

col.pal2=colorRampPalette(c("white", "yellow", "red"))(216)
for ( i in c(5:7,9:13) ){
tdata=dcast(plant_myco_db_sel, mycotoxins  ~ plants , 
            value.var = as.character(names(plant_myco_db_sel)[i]), 
            fun.aggregate = sum)
row.names(tdata)=tdata$mycotoxins
data=tdata[,2:length(tdata)]
pheatmap(data, color = col.pal2,cluster_rows = F, 
         cluster_cols = F,
         angle_col=0,
         border_color="gray",
         cellwidth = 60, cellheight = 30,
         main=paste("Heatmap  of ",firstup(names(plant_myco_db_sel)[i])),
         filename = paste0("Heatmap_",as.character(names(plant_myco_db_sel)[i]),".png"))
}


######################################################################################
# plotting fit
plant_myco_db_sel$check_names=unlist(res_names[index])
plant_myco_db_sel[,c(1:2,14)] # check only

res_dists=list()
res_dists_names=list()
res_summary=list()

z=1

for ( i in index) {
x=as.numeric(res_pooled[[i]])
outfilepdf=gsub(" ","_",paste0(as.character(res_names[i]),"_dists.pdf"))
outfilepng=gsub(" ","_",paste0(as.character(res_names[i]),"_dists.png"))
#######################################################
dists=list()
j=1
idlist=c(NULL)
fw <- try(fitdist(x,"weibull"))
if (!grepl("Error",fw)) {dists[[j]]=fw;idlist=c("weibull");j=j+1}

fg <- try(fitdist(x,"gamma"))
if (!grepl("Error",fg)) {dists[[j]]=fg;idlist=c(idlist,"gamma"); j=j+1 } 

fln <- try(fitdist(x,"lnorm"))
if (!grepl("Error",fln)) {dists[[j]]=fln;idlist=c(idlist,"lognormal");j=j+1 }

flnor <- try(fitdist(x,"norm"))
if (!grepl("Error",flnor)) {dists[[j]]=flnor;idlist=c(idlist,"normal") }

ds=denscomp(dists,legendtext=eval(idlist),plotstyle = "ggplot")
qq=qqcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot")
cd=cdfcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot") 
pp=ppcomp(dists,legendtext=eval(idlist),plotstyle = "ggplot")
out=ggarrange(ds, qq, cd,pp + rremove("x.text"), ncol = 2, nrow = 2)
options(warn=0)
ggsave(outfilepdf,plot=out,device="pdf")
ggsave(outfilepng,plot=out,device="png")
#######################################################
res_dists_names[[z]]=idlist
res_dists[[z]]=dists
res_summary[[z]]=summary_large(x)
z=z+1
}

###############################################################################################
b=melt(unlist(lapply(res_dists,length)))
b$names=rownames(b)
res=list()
for ( i in 1:nrow(b)) {res[[i]]=rep(b$names[i],b$value[i])}
names_dist=unlist(res)

names(res_dists)=NULL
list_results=do.call(c, res_dists)
list_results_type=do.call(c,res_dists_names)

resdf=data.frame(c(names_dist[1],names_dist),c(list_results_type[1],list_results_type),do.call("rbind",lapply(list_results,function(x) names(x$estimate))),
           do.call("rbind",lapply(list_results,function(x) as.numeric(x$estimate))),
           do.call("rbind",lapply(list_results,function(x) c(x$aic,x$bic)))
          )

names(resdf)=c("PlantTtox","Type_distrib","Namepar_1","Namepar_2","Valpar_1","Valpar_2","AIC","BIC")
###############################################################################################

table_describe=data.frame(PlantTox=unlist(res_names[index]),do.call("rbind",res_summary))


file.remove("PlantTox_select_score.xls")
XLConnect::writeWorksheetToFile("PlantTox_select_score.xls",plant_myco_db_sel,"data_bib_scores")
XLConnect::writeWorksheetToFile("PlantTox_select_score.xls",resdf,"Fit stasts")
XLConnect::writeWorksheetToFile("PlantTox_select_score.xls",table_describe,"data summaries")
file.remove("PlantTox_full_score.xls")
XLConnect::writeWorksheetToFile("PlantTox_full_score.xls",plant_myco_db,"data_bib_scores fiveD")


###################################################################################
# Discrete Data Analysis with R: Visualization and Modeling Techniques 
# for Categorical and Count Data 
#' @Book{FriendlyMeyer:2016:DDAR,
#'   title    = {Discrete Data Analysis with R: Visualization and Modeling Techniques for Categorical and Count Data},
#'   year     = {2016},
#'   author   = {Friendly, Michael and Meyer, David},
#'   publisher    = {Chapman \& Hall/CRC},
#'   address  = {Boca Raton, FL},
#'   isbn     = {978-1-4987-2583-5},
#' }
#' 
#' 
#' Friendly, M. & Meyer, D. (2016). Discrete Data Analysis with R: Visualization and Modeling Techniques for Categorical and Count Data. Boca Raton, FL: Chapman & Hall/CRC.

# http://ddar.datavis.ca/