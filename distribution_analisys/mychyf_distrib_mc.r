library(fitdistrplus)
library(readxl)
library(mc2d)
library(VGAM)
library(car)

#######################################################################################

setwd("/home/alf/Scrivania/lav_michyf/distrib")
source("aux_distrib.r")



MAIZE=readRDS("MAIZE.rds")
MAIZE_co=as.data.frame(MAIZE[which(MAIZE$`Co-occurrence`==1),])


U_param_co=unique(MAIZE_co$paramType)
U_type_co=unique(MAIZE_co$sampMatType)
U_plant_co=unique(MAIZE_co$sampMatbased)[1]

U_param=unique(MAIZE$paramType)
U_type=unique(MAIZE$sampMatType)
U_plant=unique(MAIZE$sampMatbased)[1]




##################################################################################################################
# don_on_feed_mtot_full_maize   

res_ls=list()

res_ls[[1]]=info_extract_db(MAIZE,U_param[17],U_type[2],"maize")
res_ls[[1]]$values=as.numeric(na.omit(res_ls[[1]]$mat_op$meanTot_op))
res_ls[[1]]$max_val=max(as.numeric(na.omit(res_ls[[1]]$values)),na.rm=T)
res_ls[[1]]$min_val=min(as.numeric(na.omit(res_ls[[1]]$values)),na.rm=T)
res_ls[[1]]$mc=mcstoc(rempiricalC, type = "VU", nsv=1000,nsu=100,values = res_ls[[1]]$values, min=res_ls[[1]]$min_val, max=res_ls[[1]]$max_val)

png("don_on_feed_mtot_full_maize_MC.png")
plot(res_ls[[1]]$mc,main="\n\n\n\n\nMC Simulation CDF\nDON  on Feed -  Maize\n Nvar=1000\nNunc=100\nN=6",xlab= "MicroG/kg",ylab="CDF")
dev.off()

# don_on raw material_conc_full_maize  U_param[16]

res_ls[[2]]=info_extract_db(MAIZE,U_param[17],U_type[3],"maize")
res_ls[[2]]$values=as.numeric(na.omit(res_ls[[2]]$mat_op$Concentration_op))
res_ls[[2]]$max_val=max(as.numeric(na.omit(res_ls[[2]]$values)),na.rm=T)
res_ls[[2]]$min_val=min(as.numeric(na.omit(res_ls[[2]]$values)),na.rm=T)
res_ls[[2]]$mc=mcstoc(rempiricalC, type = "VU", nsv=1000,nsu=100,values = res_ls[[2]]$values, min=res_ls[[2]]$min_val, max=res_ls[[2]]$max_val)

png("don_raw material_conc_full_maize_MC.png")
plot(res_ls[[2]]$mc,main="\n\n\n\nMC Simulation CDF DON  on raw material -  Maize\n Nvar=1000\nNunc=100\nN=93",xlab= "MicroG/kg",ylab="CDF")
dev.off()

# fb1+fb2_on_feed_mtot_full_maize

res_ls[[3]]=info_extract_db(MAIZE,U_param[30],U_type[2],"maize")
res_ls[[3]]$values=as.numeric(na.omit(res_ls[[3]]$mat_op$meanTot_op))
res_ls[[3]]$max_val=max(as.numeric(na.omit(res_ls[[3]]$values)),na.rm=T)
res_ls[[3]]$min_val=min(as.numeric(na.omit(res_ls[[3]]$values)),na.rm=T)
res_ls[[3]]$mc=mcstoc(rempiricalC, type = "VU", nsv=1000,nsu=100,values = res_ls[[3]]$values, min=res_ls[[3]]$min_val, max=res_ls[[3]]$max_val)

png("fb1+fb2_on_feed_mtot_full_maize_MC.png")
plot(res_ls[[3]]$mc,main="\n\n\n\nMC Simulation CDF FB1+FB2  on feed -  Maize\n Nvar=1000\nNunc=100\nN=5",xlab= "MicroG/kg",ylab="CDF")
dev.off()

##################################################################################################################
# QQplots
##################################################################################################################
# don_on_feed_mtot_full_maize   


png("don_on_feed_mtot_full_maize_qqplot.png")
qqnorm(res_ls[[1]]$values,main="QQplot DON  on Feed -  Maize \nN=6")
dev.off()

png("don_on_feed_mtot_full_maize_qqplot_MC.png")
qqnorm(res_ls[[1]]$mc,main="QQplot MC DON  on Feed -  Maize\n Nvar=1000 Nunc=100 N=6")
dev.off()

# don_on raw material_conc_full_maize  U_param[16]


png("don_raw material_conc_full_maize_qplot.png")
qqnorm(res_ls[[2]]$values,main="DON  on raw material -  Maize\nN=93")
dev.off()


png("don_raw material_conc_full_maize_qqplot_MC.png")
qqnorm(res_ls[[2]]$mc,main="QQplot MC DON  on raw material -  Maize\n Nvar=1000 Nunc=100 N=93")
dev.off()

# fb1+fb2_on_feed_mtot_full_maize


png("fb1+fb2_on_feed_mtot_full_maize_qplot.png")
qqnorm(res_ls[[3]]$values,main="QQplot FB1+FB2  on feed - Maize\n N=5")
dev.off()

png("fb1+fb2_on_feed_mtot_full_maize_qplot_MC.png")
qqnorm(res_ls[[3]]$mc,main="QQplot MC FB1+FB2  on feed -  Maize\n Nvar=1000 Nunc=100 N=5")
dev.off()

##################################################################################################################





png("don_on_feed_mtot_full_maize_qqplot_poisson.png")
qqPlot(res_ls[[1]]$values,distribution="pois",lambda=mean(res_ls[[1]]$values),main="QQplot Poisson DON  on Feed -  Maize \nN=6",ylab= "MicroG/kg")
dev.off()

png("don_raw material_conc_full_maize_qqplot_poisson.png")
qqPlot(res_ls[[2]]$values,distribution="pois",lambda=mean(res_ls[[2]]$values),main="QQplot Poisson DON  on raw material -  Maize \nN=93",ylab= "MicroG/kg")
dev.off()

png("fb1+fb2_on_feed_mtot_full_qqplot_poisson.png")
qqPlot(res_ls[[3]]$values,distribution="pois",lambda=mean(res_ls[[3]]$values),main="QQplot Poisson FB1+FB2  on feed -  Maize\nN=5",ylab= "MicroG/kg")
dev.off()

#########################################################################################################
load(url("http://www.openintro.org/stat/data/bdims.RData"))
fdims = subset(bdims, bdims$sex == 0)

qqnorm(fdims$wgt, col=adjustcolor("orange", 0.4), pch=19)
qqline(fdims$wgt)

qqnormsim = function(dat, dim=c(2,2)) {
  par(mfrow=dim)
  qqnorm(dat, col=adjustcolor("orange", 0.4), 
         pch=19, cex=0.7, main="Normal QQ Plot (Data)")
  qqline(dat)
  for (i in 1:(prod(dim) - 1)) {
    simnorm = rnorm(n=length(dat), mean=mean(dat), sd=sd(dat))
    qqnorm(simnorm, col=adjustcolor("orange", 0.4), 
           pch=19, cex=0.7,
           main="Normal QQ Plot (Sim)")
    qqline(simnorm)
  }
  par(mfrow=c(1, 1))
}
qqnormsim(fdims$wgt)
#########################################################################################################








#################################################################################################################################################################
# References

# https://rpubs.com/ajlarso2/working_draft_RA_campylobacter
# http://rstudio-pubs-static.s3.amazonaws.com/179995_71c735332ca34f19b4ef78bd8f693744.html
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4423557/
# The Two-Dimensional Monte Carlo: A New Methodologic Paradigm for Dose Reconstruction for Epidemiological Studies
# https://cran.r-project.org/web/packages/mc2d/vignettes/docmcEnglish.pdf

# rempiricalC

#https://ehp.niehs.nih.gov/ehp283/

# Two-dimensional Monte Carlo simulations consist of multiple simulations addressing variability 
# (i.e., the variability dimension) nested within a single larger simulation that addresses uncertainty (i.e., the uncertainty dimension). In this context, the term uncertainty is akin to measurement error; 
# it refers to heterogeneity in QMRA inputs that could theoretically be reduced by collecting more measurements. In contrast, our use of the term variability refers to heterogeneity in QMRA inputs that is not reducible 
# because it is due to natural processes. Heterogeneity in risk factors (e.g., age, exposure time) was assigned to the variability dimension, while heterogeneity in uncertain model parameters (e.g., dose–response parameters)
# was assigned to the uncertainty dimension (see Table 1 for the dimension associated with each individual input variable).
