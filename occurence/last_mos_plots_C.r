#############################################################################################################################
# Setup working directory

setwd("")
##############################################################################################################################
# Load libraries & functions . Check if packages were previously installed.

source("load_lib.r")
source("aux_mycosources.r")

##############################################################################################################################
# load data

res_tot=readRDS("data/res_tot_min5.rds")
res_data=readRDS("data/res_data_min5.rds")
res_names=readRDS("data/res_names_min5.rds")
res_pooled=readRDS("data/res_pooled_min5.rds")
plant_myco_db=readRDS("data/plant_myco_db.rds")
plant_myco_db_occ=readRDS("data/plant_myco_db_occ_tot.rds")

#########################################################################################################
# filtering by general final bibliometrics score

tresh=2.5

index=which(plant_myco_db$scoreGEN>tresh)

plant_myco_db_sel=plant_myco_db[index,]

names(plant_myco_db_sel)

#########################################################################################################
# create color palette

col.pal2=colorRampPalette(c("white", "yellow", "red"))(216)

for ( i in c(5:7,9:13) )  {

tdata=dcast(plant_myco_db_sel, mycotoxins  ~ plants , 
            value.var = as.character(names(plant_myco_db_sel)[i]), 
            fun.aggregate = sum)

row.names(tdata)=tdata$mycotoxins

data=tdata[,2:length(tdata)]

######################################################################################
# Plotting heatmaps


pheatmap(data, color = col.pal2,cluster_rows = F, 
         cluster_cols = F,
         angle_col=0,
         border_color="gray",
         cellwidth = 60, cellheight = 30,
         main=paste("Heatmap  of ",firstup(names(plant_myco_db_sel)[i])),
         fontsize_row=10,
         filename = paste0("Heatmap_",as.character(names(plant_myco_db_sel)[i]),".png")
         )
}

mycotoxs =unique(as.character(plant_myco_db_sel$mycotoxins))

for ( j in 1:length(mycotoxs)) {
  
  data=subset(plant_myco_db_sel,mycotoxins==mycotoxs[j])[,c(1,5:7,9:12)]
  row.names(data)=data$plants
  data=data[,2:length(data)]
  pheatmap(data, color = col.pal2,cluster_rows = F, 
           cluster_cols = F,
           angle_col=0,
           border_color="gray",
           cellwidth = 60, cellheight = 30,
           main=paste("Heatmap  of ",as.character(mycotoxs[j])),
           fontsize=8,
           height=3,
           filename = paste0("Heatmap_",as.character(mycotoxs[j]),".png")
  )
}



######################################################################################
# Plotting distrib fits

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

#ggsave(outfilepdf,plot=out, width = 8, height = 7,device="pdf")
#ggsave(outfilepng,plot=out, width = 8, height = 7,device="png")

#######################################################
res_dists_names[[z]]=idlist
res_dists[[z]]=dists
res_summary[[z]]=summary_large(x)
z=z+1
}

###########################################################################################################################################
# WRITE DATA FIT STATITICS AND SCORES

b=melt(unlist(lapply(res_dists,length)))
b$names=rownames(b)


res=list()

for ( i in 1:nrow(b)) {res[[i]]=rep(b$names[i],b$value[i])}

names_dist=unlist(res)
names(res_dists)=NULL

list_results=do.call(c, res_dists)
list_results_type=do.call(c,res_dists_names)

resdf=data.frame(plant_myco_db_sel[as.numeric(c(names_dist)),c("check_names")],
                 c(list_results_type),
                 do.call("rbind",lapply(list_results,function(x) names(x$estimate))),
                 do.call("rbind",lapply(list_results,function(x) as.numeric(x$estimate))),
                 do.call("rbind",lapply(list_results,function(x) c(x$aic,x$bic)))
          )

names(resdf)=c("Planttox","Type_distrib","Namepar_1","Namepar_2","Valpar_1","Valpar_2","AIC","BIC")


###############################################################################################

table_describe=data.frame(PlantTox=unlist(res_names[index]),do.call("rbind",res_summary))

file.remove("data/PlantTox_select_score.xls")
XLConnect::writeWorksheetToFile("data/PlantTox_select_score.xls",plant_myco_db_sel,"data_bib_scores")
XLConnect::writeWorksheetToFile("data/PlantTox_select_score.xls",resdf,"Fit stasts")
XLConnect::writeWorksheetToFile("data/PlantTox_select_score.xls",table_describe,"data summaries")
file.remove("data/PlantTox_full_score.xls")
XLConnect::writeWorksheetToFile("data/PlantTox_full_score.xls",plant_myco_db,"data_bib_scores fiveD")

plant_sel=c("barley","cereals","maize","oat","rice","rye","wheat")

png("images/boxplot_occurence_score_GEN.png")
bwplot(scoreGEN~plants, data =plant_myco_db[which(plant_myco_db$plants %in%  plant_sel ==T),])
dev.off()

png("images/boxplot_co_occurenze_score_GEN.png")
bwplot(scoreGEN~plants, data =plant_myco_db_occ[which(plant_myco_db_occ$plants %in%  plant_sel ==T),])
dev.off()

#############################################################################################################################################
# @references
#  http://ddar.datavis.ca/
#
#' Friendly, M. & Meyer, D. (2016). Discrete Data Analysis with R: Visualization and Modeling Techniques for Categorical and Count Data. Boca Raton, FL: Chapman & Hall/CRC.
#

