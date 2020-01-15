setwd("C:/Users/Robypalumb/OneDrive - Università Cattolica del Sacro Cuore/TESI/MULTINOMIAL THESIS")

library(readxl)
library(ggplot2)
library(multcomp)
library(nnet)
library(XLConnect)
library(dplyr)

plants=excel_sheets("co_occurrence_multinomial_4_12_2019.xlsx")

################################     MAIZE          ######################################


i<-2

plant=read_xlsx("co_occurrence_multinomial_4_12_2019.xlsx",i)

plant_data=plant[which(!is.na(plant$COUNT)),]

plant_data=plant_data %>% group_by(DON,FB,AF) %>% summarize(COUNT=sum(COUNT))

y<-as.matrix(plant_data[,1:3])

plant_model=multinom(formula = y~1, weights= COUNT, data = plant_data)


cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))

print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
#XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
#}
#sink()

### Simulating samples taken and the presence of one of the micotoxin ###
### Estimating the multinomial probability for maize from Literature ###
maize<-colSums(out)/dim(out)[1]

### Number of simulations ###
nsim<-1000

### Setting the seed ###
set.seed(12345)


### Generating presence abscense of each of the micotoxins in each sample ###
sampleCooccu<-(t(rmultinom(nsim,3,maize))>0)*1

### Summary table describing the amount of samples with one two or three micotoxins present ###
table(rowSums(sampleCooccu))

data.frame(sampleCooccu) %>% group_by(DON,FB,AF) %>% summarize(COUNT=n()/nsim)

sampleCooccu[1:10,]


############################   BARLEY   ##################################################

i<-1

plant=read_xlsx("co_occurrence_multinomial_4_12_2019.xlsx",i)

plant_data=plant[which(!is.na(plant$COUNT)),]

plant_data=plant_data %>% group_by(DON,NIV,ZEN) %>% summarize(COUNT=sum(COUNT))

y<-as.matrix(plant_data[,1:3])

plant_model=multinom(formula = y~1, weights= COUNT, data = plant_data)


cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))

print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
#XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
#}
#sink()

### Simulating samples taken and the presence of one of the micotoxin ###
### Estimating the multinomial probability for maize from Literature ###
barley<-colSums(out)/dim(out)[1]

### Number of simulations ###
nsim<-1000

### Setting the seed ###
set.seed(12345)


### Generating presence abscense of each of the micotoxins in each sample ###
sampleCooccu<-(t(rmultinom(nsim,3,barley))>0)*1

### Summary table describing the amount of samples with one two or three micotoxins present ###
table(rowSums(sampleCooccu))

data.frame(sampleCooccu) %>% group_by(DON,NIV,ZEN) %>% summarize(COUNT=n()/nsim)

sampleCooccu[1:10,]

############################   OAT   ##################################################

i<-3

plant=read_xlsx("co_occurrence_multinomial_4_12_2019.xlsx",i)

plant_data=plant[which(!is.na(plant$COUNT)),]

plant_data=plant_data %>% group_by(DON,T2_HT2,NIV) %>% summarize(COUNT=sum(COUNT))

y<-as.matrix(plant_data[,1:3])

plant_model=multinom(formula = y~1, weights= COUNT, data = plant_data)


cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))

print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
#XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
#}
#sink()

### Simulating samples taken and the presence of one of the micotoxin ###
### Estimating the multinomial probability for maize from Literature ###
oat<-colSums(out)/dim(out)[1]

### Number of simulations ###
nsim<-1000

### Setting the seed ###
set.seed(12345)


### Generating presence abscense of each of the micotoxins in each sample ###
sampleCooccu<-(t(rmultinom(nsim,3,oat))>0)*1

### Summary table describing the amount of samples with one two or three micotoxins present ###
table(rowSums(sampleCooccu))

data.frame(sampleCooccu) %>% group_by(DON,T2_HT2,NIV) %>% summarize(COUNT=n()/nsim)

sampleCooccu[1:10,]

############################   OAT   ##################################################

i<-3

plant=read_xlsx("co_occurrence_multinomial_4_12_2019.xlsx",i)

plant_data=plant[which(!is.na(plant$COUNT)),]

plant_data=plant_data %>% group_by(DON,T2_HT2,NIV) %>% summarize(COUNT=sum(COUNT))

y<-as.matrix(plant_data[,1:3])

plant_model=multinom(formula = y~1, weights= COUNT, data = plant_data)


cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))

print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
#XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
#}
#sink()

### Simulating samples taken and the presence of one of the micotoxin ###
### Estimating the multinomial probability for maize from Literature ###
oat<-colSums(out)/dim(out)[1]

### Number of simulations ###
nsim<-1000

### Setting the seed ###
set.seed(12345)


### Generating presence abscense of each of the micotoxins in each sample ###
sampleCooccu<-(t(rmultinom(nsim,3,oat))>0)*1

### Summary table describing the amount of samples with one two or three micotoxins present ###
table(rowSums(sampleCooccu))

data.frame(sampleCooccu) %>% group_by(DON,T2_HT2,NIV) %>% summarize(COUNT=n()/nsim)

sampleCooccu[1:10,]

############################   WHEAT   ##################################################

i<-4

plant=read_xlsx("co_occurrence_multinomial_4_12_2019.xlsx",i)

plant_data=plant[which(!is.na(plant$COUNT)),]

plant_data=plant_data %>% group_by(DON,NIV,ZEN) %>% summarize(COUNT=sum(COUNT))

y<-as.matrix(plant_data[,1:3])

plant_model=multinom(formula = y~1, weights= COUNT, data = plant_data)


cat(paste("Modeling ", plants[i],"\n"))
cat(paste(formula_model,"\n"))

print(summary(plant_model))
out=as.data.frame(fitted.values(plant_model))
names(out)=names(plant_data)[1:3]
#XLConnect::writeWorksheetToFile("res_multinomial.xls",out,plants[i])
#}
#sink()

### Simulating samples taken and the presence of one of the micotoxin ###
### Estimating the multinomial probability for maize from Literature ###
wheat<-colSums(out)/dim(out)[1]

### Number of simulations ###
nsim<-1000

### Setting the seed ###
set.seed(12345)


### Generating presence abscense of each of the micotoxins in each sample ###
sampleCooccu<-(t(rmultinom(nsim,3,wheat))>0)*1

### Summary table describing the amount of samples with one two or three micotoxins present ###
table(rowSums(sampleCooccu))

data.frame(sampleCooccu) %>% group_by(DON,NIV,ZEN) %>% summarize(COUNT=n()/nsim)

sampleCooccu[1:10,]


