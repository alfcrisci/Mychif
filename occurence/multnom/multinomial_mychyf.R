setwd("/home/alf/Scrivania/lav_michyf/mychif_tasks/multnom")

library(readxl)
library(ggplot2)
library(multcomp)
library(nnet)
library(XLConnect)


plants=excel_sheets("co_occurrence_multinomial.xlsx")

##################################################################################

file.remove("res_multinomial.xls")
file.remove("modeling.txt")


#sink("modeling.txt")

for ( i in 1:4) {
plant=read_xlsx("co_occurrence_multinomial.xlsx",i)
plant_data=plant[which(!is.na(plant$COUNT)),]

formula_model=paste(paste(names(plant_data)[1:3],collapse = "+"),"~ COUNT")

plant_model=multinom(formula = formula_model, data = plant_data)
cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))

print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
}
#sink()

