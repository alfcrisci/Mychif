#############################################################
# Modelling dose - responses

library(drc)
library(mixtox)
library(bmdModeling)
library(synergyfinder)
library(multcomp)
library(bmd)

# library(BIGL)

#############################################################
# Read and graphics

library(readxl)
library(ggplot2)

#############################################################
# Setup working directory

setwd("")

#############################################################
# Setup working directory

names_end=excel_sheets("data/TOX_VIVO_CHICKEN_crisci.xls")

endlen=length(names_end)

cases_ls=list()

for ( i in 1:endlen) {
  
  cases_ls[[i]]=as.data.frame(read_excel("data/TOX_VIVO_CHICKEN_crisci.xls",names_end[i]))

  }

saveRDS(cases_ls,"cases_ls.rds")


###############################################################################################
#  

cases_ls_liv=list()
cases_ls_mar=list()

for(  i in 1:endlen) {
  cases_ls_liv[[i]]=cases_ls[[i]][which(cases_ls[[i]]$Specie=="White Leghorn chickens"),]
  cases_ls_mar[[i]]=cases_ls[[i]][which(cases_ls[[i]]$Specie=="Marek chickens"),]
}

saveRDS(cases_ls_liv,"cases_ls_liv.rds")
saveRDS(cases_ls_mar,"cases_ls_mar.rds")

############################################################
# 'LL.3' and 'LL2.3' provide the three-parameter 
#  log-logistic function where the lower limit is equal to 0. '

############################################################
old.par=par()
par(old.par)
set.seed(2)
range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

############################################################
# BW  gain
# Marek

datamyc=cases_ls_mar[[1]]
dataexpand=data.frame(Bodyweight_gain=rep(datamyc$Bodyweight_gain,100)+runif(length(datamyc$Bodyweight_gain)*100,-1.5,1.5)*rep(datamyc$Bodyweight_gain.sd,100),
                      Mycotoxin=rep(datamyc$Mycotoxin,100),
                      Dose_T=rep(datamyc$Dose_T,100))



############################################################
# BW  gain
# Marek

datamyc=cases_ls_mar[[1]]
dataexpand=data.frame(Bodyweight_gain=rep(datamyc$Bodyweight_gain,100)+runif(length(datamyc$Bodyweight_gain)*100,-1.5,1.5)*rep(datamyc$Bodyweight_gain.sd,100),
                      Mycotoxin=rep(datamyc$Mycotoxin,100),
                      Dose_T=rep(datamyc$Dose_T,100))

dataexpand$response=range01(dataexpand$Bodyweight_gain)

multi.m3=drm(Bodyweight_gain~Dose_T,Mycotoxin,data=dataexpand,fct = AR.3())
#multi.s=drm(response~Dose_T,Mycotoxin,data=dataexpand,fct = LL.2())

ED(multi.m3,c(10,50), interval = "delta")

file.remove("summary_BW_gain_Marek.txt")

cat(capture.output(ED(multi.m3,c(10,50), interval = "delta")),
    file = "summary_BW_gain_Marek.txt",
    sep="\n",append = T)



cat(capture.output(summary(multi.m3)),
    file = "summary_BW_gain_Marek.txt",
    sep="\n",append = T)

png("images/BW_gain_Marek.png")

plot(multi.m3, 
     col = TRUE, 
     xlab = "Log Dose(mg/kg feed)",
     ylab = "BW_gain", 
     main="Bodyweight_gain AF & OTA Marek chickens\n Bootstrap=Yes Model=L.3")

dev.off()





############################################################

datamyc=cases_ls_liv[[1]]

multi.m3 <- drm(Terminal_bodyweight~Dose_T,Mycotoxin,data=datamyc,fct = MM.3())
print(summary(multi.m3))

file.remove("summary_TW_gain_Liv.txt")

cat(capture.output(ED(multi.m3,c(10,50), interval = "delta")),
    file = "summary_TW_gain_Liv.txt",
    sep="\n",append = T)

cat(capture.output(summary(multi.m3)),
    file = "summary_TW_gain_Liv.txt",
    sep="\n")

png("images/TW_gain_Liv.png")

plot(multi.m3, 
     col = TRUE, 
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Terminal BW", 
     main="Terminal BW AF & OTA  White Leghorn chickens\n Bootstrap=Not Model=MM.3")
dev.off()



########################################################################################à
# Feed conversion #Marek

datamyc=cases_ls_mar[[2]]
dataexpand=data.frame(Feed_conversion=rep(datamyc$Feed_conversion,100)+runif(length(datamyc$Feed_conversion)*100,-1.5,1.5)*rep(datamyc$Feed_conversion.sd,100),
                      Mycotoxin=rep(datamyc$Mycotoxin,100),
                      Dose_T=rep(datamyc$Dose_T,100))

multi.m3 <- drm(Feed_conversion~Dose_T,Mycotoxin,data=dataexpand,fct = AR.3())

print(summary(multi.m3))

file.remove("summary_Feed_conversion_mar.txt")

cat(capture.output(ED(multi.m3,c(10,50), interval = "delta")),
    file = "summary_Feed_conversion_mar.txt",
    sep="\n",append = T)

cat(capture.output(summary(multi.m3)),
    file = "summary_Feed_conversion_mar.txt",
    sep="\n")

png("images/Feed_conversion_mar.png")

plot(multi.m3, 
     col = TRUE, 
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Feed_conversion", 
     main="Feed_conversion AF & OTA  Marek chickens\n Bootstrap=Yes Model=AR.3",
     legendPos = c(0.4, 2.8))

dev.off()

#########################################################################################
# Feed intake #Marek

datamyc=cases_ls_mar[[3]]
multi.m3 <- drm(Feed_intake~Dose_T,Mycotoxin,data=datamyc,fct = AR.3())
print(summary(multi.m3))


png("images/Feed_intake_mar_bootno.png")

plot(multi.m3, 
     col = TRUE, 
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Feed_intake",
     main="Feed_intake AF & OTA Marek chickens\n Bootstrap=No Model=AR.3")

dev.off()

dataexpand=data.frame(Feed_intake=rep(datamyc$Feed_intake,100)+runif(length(datamyc$Feed_intake)*100,-1.5,1.5)*rep(datamyc$Feed_intake.sd,100),
                      Mycotoxin=rep(datamyc$Mycotoxin,100),
                      Dose_T=rep(datamyc$Dose_T,100))


multi.m3 <- drm(Feed_intake~Dose_T,Mycotoxin,data=dataexpand,fct = AR.3())
print(summary(multi.m3))

file.remove("summary_Feed_intake_mar.txt")

cat(capture.output(ED(multi.m3,c(10,50), interval = "delta")),
    file = "summary_Feed_intake_mar.txt",
    sep="\n",append = T)

cat(capture.output(summary(multi.m3)),
    file = "summary_Feed_intake_mar.txt",
    sep="\n",append = T)

png("images/Feed_intake_mar.png")

plot(multi.m3, 
     col = TRUE, 
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Feed_intake",
     main=" Feed_intake AF & OTA  Marek chickens\n Bootstrap=Yes Model=AR.3")
dev.off()

##########################################################################################

datamyc=cases_ls_liv[[3]]
multi.m3 <- drm(Feed_intake~Dose_T,Mycotoxin,data=datamyc,fct = L.3())
print(summary(multi.m3))

file.remove("summary_Feed_intake_liv.txt")

cat(capture.output(ED(multi.m3,c(10,50), interval = "delta")),
    file = "summary_Feed_intake_liv.txt",
    sep="\n",append = T)

cat(capture.output(summary(multi.m3)),
    file = "summary_Feed_intake_mar.txt",
    sep="\n",append = T)


png("images/Feed_intake_liv_bootno.png")

plot(multi.m3,
     col = TRUE,
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Feed_intake",
     xlim=c(0,10),
     main="Feed_intake AF & OTA   White Leghorn chickens\n Bootstrap=No Model=L.3",
     legendPos = c(9, 5400))

dev.off()



##########################################################################################
# Liver_relative_weight #Marek


datamyc=cases_ls_mar[[4]]

dataexpand=data.frame(Liver_relative_weight=rep(datamyc$Liver_relative_weight,100)+runif(length(datamyc$Liver_relative_weight)*100,-1.2,1.2)*rep(datamyc$Liver_relative_weight.sd,100),
                      Mycotoxin=rep(datamyc$Mycotoxin,100),
                      Dose_T=rep(datamyc$Dose_T,100))

multi.m3 <- drm(Liver_relative_weight~Dose_T,Mycotoxin,data=dataexpand,fct = AR.3())
print(summary(multi.m3))
file.remove("summary_Liver_RW_mar.txt")

cat(capture.output(ED(multi.m3,c(10,50), interval = "delta")),
    file = "summary_Liver_RW_mar.txt",
    sep="\n",append = T)

cat(capture.output(summary(multi.m3)),
    file = "summary_Liver_RW_mar.txt",
    sep="\n",append = T)


png("images/Liver_RW_mar.png")

plot(multi.m3, col = TRUE,
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Liver RW", 
     main="Liver RW AF & OTA  Marek chickens\n Bootstrap=Yes Model=AR.3",
     legendPos = c(0.5, 4.8))

dev.off()




##########################################################################################
# Kidney_relative_weight#Marek

datamyc=cases_ls_mar[[5]]

dataexpand=data.frame(Kidney_relative_weight=rep(datamyc$Kidney_relative_weight,100)+runif(length(datamyc$Kidney_relative_weight)*100,-1.5,1.5)*rep(datamyc$Kidney_relative_weight.sd,100),
                      Mycotoxin=rep(datamyc$Mycotoxin,100),
                      Dose_T=rep(datamyc$Dose_T,100))

multi.m3 <- drm(Kidney_relative_weight~Dose_T,Mycotoxin,data=dataexpand,fct = AR.3())
print(summary(multi.m3))

file.remove("summary_Kidney_RW_mar.txt")

cat(capture.output(ED(multi.m3,c(10,50), interval = "delta")),
    file = "summary_Kidney_RW_mar.txt",
    sep="\n",append = T)

cat(capture.output(summary(multi.m3)),
    file = "summary_Kidney_RW_mar.txt",
    sep="\n",append = T)


png("images/Kidney_RW_mar.png")

plot(multi.m3, col = TRUE,
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Kidney RW", 
     main="Kidney RW AF & OTA Marek chickens\n Bootstrap=Yes Model=AR.3",
     legendPos = c(0.5, 0.9))

dev.off()

##########################################################################################
# Spleen_relative weight

datamyc=cases_ls_mar[[6]]
multi.m4 <- drm(Spleen_relative_weight~Dose_T,Mycotoxin,data=datamyc[,],fct=L.3())
print(summary(multi.m4))

file.remove("summary_Spleen_RW_mar.txt")

cat(capture.output(ED(multi.m4,c(10,50), interval = "delta")),
    file = "summary_Spleen_RW_mar.txt",
    sep="\n",append = T)

cat(capture.output(summary(multi.m4)),
    file = "summary_Spleen_RW_mar.txt",
    sep="\n",append = T)


png("images/Spleen_RW_mar.png")

plot(multi.m4, col = TRUE,
     xlab = "Log Dose(mg/kg feed)",
     ylab = "Spleen RW", 
     main="Spleen RW AF & OTA Marek chickens\n Bootstrap=No Model=L.3",
     legendPos = c(0.4, 0.329))

dev.off()
