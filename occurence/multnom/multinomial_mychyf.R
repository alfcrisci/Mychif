setwd("/home/alf/Scrivania/lav_michyf/repositories/Mychif/occurence/multnom")

library(readxl)
library(ggplot2)
library(multcomp)
library(nnet)
library(XLConnect)


plants=excel_sheets("co_occurrence_multinomial.xlsx")

##################################################################################

file.remove("res_multinomial.xls")
file.remove("modeling.txt")

sink("modeling.txt")
i=1

plant=read_xlsx("co_occurrence_multinomial.xlsx",i)
plant_data=plant[which(!is.na(plant$COUNT)),]
formula_model=paste(paste(names(plant_data)[1:3],collapse = "+"),"~ COUNT")
test<-plant_data %>% group_by(DON,NIV,ZEN) %>% summarize(COUNT=sum(COUNT))
plant_data=test
prova=as.matrix(plant_data[,1:3])
plant_model=multinom(formula = prova~1,weights=COUNT, data = plant_data)
cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))
print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
XLConnect::writeWorksheetToFile("res_multinomial.xls",plant_data,paste(plants[i],"data"))

i=2

plant=read_xlsx("co_occurrence_multinomial.xlsx",i)
plant_data=plant[which(!is.na(plant$COUNT)),]
formula_model=paste(paste(names(plant_data)[1:3],collapse = "+"),"~ COUNT")
test<-plant_data %>% group_by(DON,FB,AF) %>% summarize(COUNT=sum(COUNT))
plant_data=test
prova=as.matrix(plant_data[,1:3])
plant_model=multinom(formula = prova~1,weights=COUNT, data = plant_data)
cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))
print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
XLConnect::writeWorksheetToFile("res_multinomial.xls",plant_data,paste(plants[i],"data"))

i=3
plant=read_xlsx("co_occurrence_multinomial.xlsx",i)
plant_data=plant[which(!is.na(plant$COUNT)),]
formula_model=paste(paste(names(plant_data)[1:3],collapse = "+"),"~ COUNT")
test<-plant_data %>% group_by(DON,T2_HT2,NIV) %>% summarize(COUNT=sum(COUNT))
plant_data=test
prova=as.matrix(plant_data[,1:3])
plant_model=multinom(formula = prova~1,weights=COUNT, data = plant_data)
cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))
print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
XLConnect::writeWorksheetToFile("res_multinomial.xls",plant_data,paste(plants[i],"data"))

i=4
plant=read_xlsx("co_occurrence_multinomial.xlsx",i)
plant_data=plant[which(!is.na(plant$COUNT)),]
formula_model=paste(paste(names(plant_data)[1:3],collapse = "+"),"~ COUNT")
test<-plant_data %>% group_by(DON,NIV,ZEN) %>% summarize(COUNT=sum(COUNT))
plant_data=test
prova=as.matrix(plant_data[,1:3])
plant_model=multinom(formula = prova~1,weights=COUNT, data = plant_data)
cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))
print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
XLConnect::writeWorksheetToFile("res_multinomial.xls",plant_data,paste(plants[i],"data"))

sink()
