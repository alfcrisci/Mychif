######################################################################################
# Load codes and libraries
######################################################################################
# 
# library(foreach)
# library(doSNOW)
# cl <- makeCluster(4, type="SOCK") # for 4 cores machine
# registerDoSNOW (cl)
######################################################################################

library(reshape)
library(ggplot2)
library(compiler)
library(XLConnect)

######################################################################################
# High-Throughput Toxicokinetics

library(httk)  
library(fitdistrplus)

######################################################################################
# devtools::install_github("nanhung/pksensi")
# devtools::install_github("nanhung/pksensi")
# devtools::install_github("nanhung/pksensi")

library(sensitivity) # Sensitivity Analysis 
library(pksensi) # Global Sensitivity Analysis in Pharmacokinetic Modeling 
library(PK) # Basic Non-Compartmental Pharmacokinetics
library(PKPDmisc) # Pharmacokinetic and Pharmacodynamic Data Management Functions

######################################################################################

library(PKPDsim) # Simulate dose regimens for PKPD models described by ODE or equatino
library(PKPDsimShiny) # 
library(shiny)
library(shinydashboard)
################################################################################################################################################################
# Set working directory
################################################################################################################################################################

setwd("/home/alf/Scrivania/lav_michyf/repositories/Mychif/pbkm_modeling")

source("R_code/PBPK_aux.R")
source("R_code/PBTK_model_core.R")


############################################################################################################################
# Loading physicochemical data

chem <- read.csv("data/chem_data_prob.csv", stringsAsFactors = FALSE)      #physicochemical parameters constant are in 1/minutes

# http://www.t3db.ca/toxins

############################################################################################################################
# Loading TK data

TK <- read.csv("data/tk_data/TK_prob.csv",stringsAsFactors=FALSE)                  # toxicokinetic parameters

############################################################################################################################
# Loading physiology data 

fBW <- read.csv("data/physio/fBW_animals_prob.csv",stringsAsFactors = FALSE)      # organ fractions, incl var-distribution parameters
fCO <- read.csv("data/physio/fCO_animals_prob.csv",stringsAsFactors = FALSE)      # blood flow fractions, incl var-distribution parameters
rates <- read.csv("data/physio/rates_animals_prob.csv",stringsAsFactors = FALSE)  # physiological rates (not chemical-dependent)
fPC <- read.csv("data/physio/PC_tissue_prob.csv",stringsAsFactors = FALSE)        # tissue composition, fixed during var-analysis

############################################################################################################################
# Defining cases for species and toxins referenced customiing population parameter

#################################################################################################################
# Saint-Cyr et al 100 μg/kg BW

idsim="Swine_DONa"
species <- "swine"                      # cat/cattle/chicken/sheep/swine
chemical <- "DON"                       # for now: 

BW_animal=27                            # kg
BW_SD=0                                 # kg standard deviation 
F_dose=0.75                             # %
Cl_don=0.52                             # 0.52 ± 0.25 L/h/kg
Vd_don=2.01                             # ± 0.67 L/kg
kel_don=0.003660697                     # 1/min  # 0.2567164 1/h
kabs_don=0.062                          # 1/min   # 3.72 1/h
E_dose <-0.1*F_dose*BW_animal           # Single dose in 2 volte


params <- c(F = F_dose, KA = kabs_don*60, KE = kel_don*60, V = Vd_don)  
t <- seq(0, 24, 0.01)
C <-FFPK(params = params, time = t,dose=100)
res=PKPDmisc::nca(t, C, dose=100, last_times = c(3, 4, 5), digits = 2)

png("case_a.png")
plot(t, C, type = "l",xlab = "Hours", ylab = "DON concentration microg/l", main="PK 1cmpt Pig DON \nSaint-Cyr et al 2015")
text(15,15,labels=paste0("Dose=100 μg/kg BW\n",
                         "F= ",F_dose*100,"%\n",
                         "Cmax= ",res$Cmax,"\n",
                         "Tmax= ",res$Tmax," h\n",
                         "half_life= ",res$half_life," h\n"),col="red")

dev.off()

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_don
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=0.00001
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=Cl_don
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD

n_sim  <- 1                         # number of iterations
E_start <- 0                        # start of exposure phase (h)
E_end <- 3                          # end of exposure phase (h)
E_int <- 80                         # interval between doses (h); 
t_start <- 0                  # start of simulation (h)
t_end <- 24                   # end of simulation (h)

source("R_code/simrun_bolus_oral.r") # res_sims_body res_sims_blood Two R list corresponding nsim simulation


t=as.numeric(res_sims_blood[[1]][,1])
C=as.numeric(res_sims_blood[[1]][,2])*1000
res=PKPDmisc::nca(t, C, dose=E_dose, last_times = c(3, 4, 5), digits = 2)

plot(t,C, type = "l",xlab = "Hours", ylab = "DON concentration microg/l", main="TD Pig multicmpt blood ven DON \nSaint-Cyr et al 2015")
text(15,2,labels=paste0("Dose=100 μg/kg BW\n",
                         "F= ",F_dose*100,"%\n",
                         "Cmax= ",res$Cmax,"\n",
                         "Tmax= ",res$Tmax," h\n",
                         "half_life= ",res$half_life," h\n"),col="red")

t=as.numeric(res_sims_body[[1]][,1])
C=as.numeric(res_sims_body[[1]][,2])*1000
res=PKPDmisc::nca(t, C, dose=E_dose*1000, last_times = c(3, 4, 5), digits = 2)


plot(t,C, type = "l",xlab = "Hours", ylab = "DON concentration microg/l", main="TD Pig multicmpt wholebody  DON \nSaint-Cyr et al 2015")
text(15,2,labels=paste0("Dose=100 μg/kg BW\n",
                        "F= ",F_dose*100,"%\n",
                        "Cmax= ",res$Cmax,"\n",
                        "Tmax= ",res$Tmax," h\n",
                        "half_life= ",res$half_life," h\n"),col="red")

#################################################################################################################
# Paulick et al 75 μg/kg BW

idsim="Swine_DONb"
species <- "swine"                      # cat/cattle/chicken/sheep/swine
chemical <- "DON"                         # for now: 

BW_animal=41                            # kg
BW_SD=6.7                               # kg standard deviation
F_dose=0.986
Cl_don=12.53                            #± 3.89 ml/kg min #0.01253 l/kg min
Vd_don=1.68                             # ± 0.31 l/kg
Tmax=4.85                               # h
Cmax =28.76                             # ng/mL
kel_don=0.002166667                     # 1/min 0.13 1/h
kabs_don=0.01016667                     # 1/min 0.61 1/h 
E_dose <-0.075*41*0.986                 # Single dose mg

params <- c(F = F_dose, KA = kabs_don*60, KE = kel_don*60, V = Vd_don)  
t <- seq(0, 24, 0.01)
C <-FFPK(params = params, time = t,dose=75)

plot(t, C, type = "l", xlab = "Hours", ylab = "DON concentration")
res=PKPDmisc::nca(t, C, dose=75, last_times = c(3, 4, 5), digits = 2)

png("case_b.png")
plot(t, C, type = "l",xlab = "Hours", ylab = "DON concentration microg/l", main="PK 1cmpt Pig DON \nPaulick et al 2015")
text(15,18,labels=paste0("Dose=100 μg/kg BW\n",
                         "F= ",F_dose*100,"%\n",
                         "Cmax= ",res$Cmax,"\n",
                         "Tmax= ",res$Tmax," h\n",
                         "half_life= ",res$half_life," h\n"),col="red")
dev.off()

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_don
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=0.0001
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=Cl_don/1000
TK$Vmax_tot[which(TK$species==species & TK$chemical==chemical)]=1.68   
TK$Km_tot[which(TK$species==species & TK$chemical==chemical)]=5.793103
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD

n_sim  <- 1                         # number of iterations
E_start <- 0                        # start of exposure phase (h)
E_end <- 3                          # end of exposure phase (h)
E_int <- 80                         # interval between doses (h); 
t_start <- 0                  # start of simulation (h)
t_end <- 24                   # end of simulation (h)

source("R_code/simrun_bolus_oral.r") # res_sims_body res_sims_blood Two R list corresponding nsim simulation

# case_b_blood_TD.png

t=as.numeric(res_sims_blood[[1]][,1])
C=as.numeric(res_sims_blood[[1]][,2])*1000
res=PKPDmisc::nca(t, C, dose=E_dose, last_times = c(3, 4, 5), digits = 2)
plot(t,C, type = "l",xlab = "Hours", ylab = "DON concentration microg/l", main="TD Pig multicmpt blood ven DON \nPaulick et al 2015")
text(15,10,labels=paste0("Dose=75 μg/kg BW\n",
                        "F= ",F_dose*100,"%\n",
                        "Cmax= ",res$Cmax,"\n",
                        "Tmax= ",res$Tmax," h\n",
                        "half_life= ",res$half_life," h\n"),col="red")
# case_b_body_TD.png

t=as.numeric(res_sims_body[[1]][,1])
C=as.numeric(res_sims_body[[1]][,2])*1000
res=PKPDmisc::nca(t, C, dose=E_dose*1000, last_times = c(3, 4, 5), digits = 2)
plot(t,C, type = "l",xlab = "Hours", ylab = "DON concentration microg/l", main="TD Pig multicmpt wholebody  DON \nPaulick et al 2015")
text(15,8,labels=paste0("Dose=75 μg/kg BW\n",
                        "F= ",F_dose*100,"%\n",
                        "Cmax= ",res$Cmax,"\n",
                        "Tmax= ",res$Tmax," h\n",
                        "half_life= ",res$half_life," h\n"),col="red")


# file.remove(paste("./Output/SimRes_",idsim,".xls",sep=""))
# writeWorksheetToFile(paste("./Output/SimRes_",idsim,".xls",sep=""), data_blood, sheet=paste0("blood_",idsim))



########################################################################################################
# Sensibility DON pigs

# params <- c(F = 0.75, KA = 3.72, KE = 0.29, Vd = 1.5074)  
# params <- c(F = 0.893, KA = 0.99, KE = 0.13, Vd = 0.98)  

q <- "qunif"
params <- c("F","KA","KE","V")

q.arg <- list(list(min = 0.6, max = 1.0),
              list(min = 0.5, max = 4),
              list(min = 0.02, max = 0.3),
              list(min = 1, max = 3))

x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, rep = 20)
time <- seq(0.1, 24, 0.1)
y <- solve_fun(x, model = FFPK, time = time, vars = "output")
tell2(x,y) # Link decoupling simulation result
check(x)
heat_check(x)
heat_check(x, index = "CI")
test=x
dimnames(test$y)[[4]]=" Fast99 sensitivity analisys - Pigs"
plot(test)

#################################################################################################################
# Osselaere et al
# DON 0.75 mg/kg BW ( 5mg/kg feed max DON), T2 0.02 mg/kg BW and ZEN 0.3 mg/kg BW ( max 2mg/kg feed)

idsim="chicken_DON"
species <- "chicken"                      # cat/cattle/chicken/sheep/swine
chemical <- "DON"                         # for now: chlorpyrifos/melamine/meloxicam/monensin/oxytetracycline/PFOS
BW_animal=0.96                             # kg 
BW_SD=0.1                                 # kg standard deviation 
F_dose=0.193
Cl_don=0.65                               # ± 0.217 L/min kg
Vd_don=35.72                              # ± 15.563 L/kg
kel_don=0.02                      # 1/min 
Tmax=35                                   # min
kabs_don=0.038                           # 1/min 0.6 1/h rate_kabs(35,0.01819709)
E_dose <-0.75*0.96                        # Single dose mg

params <- c(F = F_dose, KA = kabs_don*60, KE = kel_don*60, V = Vd_don)  
t <- seq(0, 24, 0.01)
C <-FFPK(params = params, time = t,dose=720)
plot(t, C, type = "l", xlab = "Hours", ylab = "DON concentration")
PKPDmisc::nca(t, C, dose=100, last_times = c(3, 4, 5), digits = 2)

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_don
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=0
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=Cl_don
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD

#################################################################################################################
# Buranatragool et al 2015 1.2 mg/kg of bw

idsim="chicken_ZEN"
species <- "chicken"                      # cat/cattle/chicken/sheep/swine
chemical <- "ZEN"                          # for now: chlorpyrifos/melamine/meloxicam/monensin/oxytetracycline/PFOS

F_dose=0.26
BW_animal=1.56                             # kg
BW_SD=2.32                                 # kg standard deviation
Vd_zea=6.40                                # ± 0.89 l/kg
Cl_zea=0.34                                # ± 0.03  l/h/kg
kel_zea=0.0002333333                      # 1/min 0.53
kabs_zea= 0.01133333                            # 1/min  0.68/60
Tmax_zea=10                                # min
Cmax=15.9                                  # ng/ml
                                           # 1/min 0.6 1/h rate_kabs(10,0.008833333 )
E_dose <-1.2*1.56*F_dose                    # Single dose mg

params <- c(F = F_dose, KA = kabs_zea*60, KE = kel_zea*60, V = Vd_zea/F_dose)  
t <- seq(0, 24, 0.01)
C <-FFPK(params = params, time = t,dose=1200)

plot(t, C, type = "l", xlab = "Hours", ylab = "ZEA concentration")
res=PKPDmisc::nca(t, C, dose=1200, last_times = c(3, 4, 5), digits = 2)


png("case_d.png")

plot(t, C, type = "l",xlab = "Hours", ylab = "ZEA concentration microg/l", main="PK 1cmpt Chicken Zea\nBuranatragool et al 2015")
text(15,4,labels=paste0("Dose=1200 μg/kg BW\n",
                         "F= ",F_dose*100,"%\n",
                         "Cmax= ",res$Cmax,"\n",
                         "Tmax= ",res$Tmax," h\n",
                         "half_life= ",res$half_life," h\n"),col="red")

dev.off()


BW_animal=1.56                             # kg
BW_SD=0                                    # kg standard deviation
E_dose <-1.2*1.56                          # Single dose mg

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_zea
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=(Cl_zea*F_dose)/60
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=0
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD

n_sim  <- 1                         # number of iterations
E_start <- 0                        # start of exposure phase (h)
E_int <- 2                          # interval between doses (h); 
E_end <- 1                          # end of exposure phase (h)
t_start <- 0                        # start of simulation (h)
t_end <- 24                         # end of simulation (h)

source("R_code/simrun_bolus_oral.r") # res_sims_body res_sims_blood Two R list corresponding nsim simulation

# case_d_blood_TD.png

t=as.numeric(res_sims_blood[[1]][,1])
C=as.numeric(res_sims_blood[[1]][,2])*(1000)
res=PKPDmisc::nca(t, C, dose=1200, last_times = c(3, 4, 5), digits = 2)
plot(t,C, type = "l",xlab = "Hours", ylab = "ZEN concentration microg/l", main="TD Chicken  multicmpt blood ven ZEN \nBuranatragool et al 2015")
text(18,80,labels=paste0("Dose=1200 μg/kg BW\n",
                         "F= ",F_dose*100,"%\n",
                         "Cmax= ",res$Cmax,"\n",
                         "Tmax= ",res$Tmax," h\n",
                         "half_life= ",res$half_life," h\n"),col="red")

# case_d_body_TD.png

t=as.numeric(res_sims_body[[1]][,1])
C=as.numeric(res_sims_body[[1]][,2])*1000
res=PKPDmisc::nca(t, C, dose=1200, last_times = c(3, 4, 5), digits = 2)
plot(t,C, type = "l",xlab = "Hours", ylab = "ZEN concentration microg/l", main="TD Chicken multicmpt  wholebody ven ZEN \nBuranatragool et al 2015")
text(15,200,labels=paste0("Dose=1200 μg/kg BW\n",
                        "F= ",F_dose*100,"%\n",
                        "Cmax= ",res$Cmax,"\n",
                        "Tmax= ",res$Tmax," h\n",
                        "half_life= ",res$half_life," h\n"),col="red")


# file.remove(paste("./Output/SimRes_",idsim,".xls",sep=""))
# writeWorksheetToFile(paste("./Output/SimRes_",idsim,".xls",sep=""), data_blood, sheet=paste0("blood_",idsim))

#################################################################################################################
# Sun et al 2.0 mg/kg b.w., every 12 h for 2 days N 20 

idsim="chicken_T2"
species <- "chicken"                      # cat/cattle/chicken/sheep/swine
chemical <- "T2"                          # 
BW_animal=1.3                             # kg
BW_SD=0.1                                 # kg standard deviation 
F_dose=0.177 
Cl_t2=0.12                                # ± 0.217 L/min/kg
Vd_t2=3.17                                # ± 0.43 L/kg
Tmax=13.2                                 # ± 4.80 min
Cmax=53.10                                # ± 10.42 ng/ml
kel_t2=0.03785489                         # 1/min 2.271293 1/h
kabs_t2=0.132                             # 1/min 7.92 1/h # rate_kabs(13.2,0.03785489)
E_dose <-1.3*2                            # Single dose mg

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_t2
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=0
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=Cl_t2
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD

n_sim  <- 1                         # number of iterations
E_start <- 0                        # start of exposure phase (h)
E_end <- 1                          # end of exposure phase (h)
E_int <- 80                         # interval between doses (h); 
t_start <- 0                        # start of simulation (h)
t_end <- 24                         # end of simulation (h)

source("R_code/simrun_bolus_oral.r") # res_sims_body res_sims_blood Two R list corresponding nsim simulation


#################################################################################################################
# Kongkapan et al 2016 0.8 mg/kg bw 

idsim="chicken_NIV"
species <- "chicken"                      # cat/cattle/chicken/sheep/swine
chemical <- "NIV"                         # 
BW_animal=1.3                             # kg  g
BW_SD=0.2                                 # kg standard deviation 
F_dose=0.0398 
Cl_NIV=113.57                             # ± 18.95 ml/h/kg 1.892833 ml/min/kg
Vd_NIV=853.93                             # ± 136.29 ml/kg 
Tmax=2.4                                  # ± 0.89 h 
Cmax=62.56                                # ± 30.86 ng/ml
kel_NIV=0.002216613                       # 1/min 0.1329968 1/h
kabs_NIV=0.01586667                       # 1/min rate_kabs(2.4,0.1329968) 0.952 1/h
E_dose <-0.8*2                            # Single dose mg

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_NIV
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=0
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=Cl_NIV/60
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD


###############################################################################################
# Freyman DOI: 10.1021/acs.jafc.6b02913 N6 0.2 mg/kg bw 

idsim="chicken_ENNB1"
species <- "chicken"                        # cat/cattle/chicken/sheep/swine
chemical <- "ENNB1"                         # 
BW_animal=0.96                               # kg  g
BW_SD=0.142                                   # kg standard deviation 
F_dose=0.05                                # mg/kg bw
Cl_ENNB1=6.64                             # ± 0,81 ml/h/kg 1.892833 ml/min/kg
Vd_ENNB1=14.36                             # ± 6.84 l/kg 
Tmax=0.63                                   # ± 0.89 h 
Cmax=1.37                                  # ± 30.86 ng/ml
kel_ENNB1=0.008666667                      # 1/min 0.52 ± 0.07 1/h
kabs_ENNB1=0.05971667                      # 1/min rate_kabs(0.63,0.52) 3.583 1/h
E_dose <-0.96*0.2                            # Single dose mg

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_ENNB1
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=0
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=Cl_ENNB1/60
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD

###############################################################################################
# Freyman DOI: 10.1021/acs.jafc.6b02913 N6 0.2 mg/kg bw

idsim="chicken_ENNB" 
species <- "chicken"                         # cat/cattle/chicken/sheep/swine
chemical <- "ENNB"                           # 
BW_animal=0.96                               # kg  g
BW_SD=0.142                                  # kg standard deviation 
F_dose=0.11                                  # mg/kg bw
Cl_ENNB=7.18                                 # ± 0,81 ml/h/kg 1.892833 ml/min/kg
Vd_ENNB=29.39                                # ± 6.84 l/kg 
Tmax=0.29                                    # ± 0.89 h 
Cmax=0.96                                    # ± 30.86 ng/ml
kel_ENNB=0.004666667                         # 1/min 0.25 ± 0.07 1/h
kabs_ENNB= 0.227                             # 1/min rate_kabs(17.4,0.004666667) 13.62 1/h
E_dose <-0.96*0.2                            # Single dose mg

TK$kabs[which(TK$species==species & TK$chemical==chemical)]=kabs_ENNB
TK$Cl_renal[which(TK$species==species & TK$chemical==chemical)]=0
TK$Cl_hepatic[which(TK$species==species & TK$chemical==chemical)]=Cl_ENNB
TK$fbact[which(TK$species==species & TK$chemical==chemical)]=0
fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD

###############################################################################################
# Toxicodinamics

############################################################################################################################
# Defined exposure scenario

E_dose=E_dose                       # exposure dose : mg/kgfeed at continuous; 
n_sim  <- 2                   # number of iterations

regime <- "bolus"                   # exposure regime (bolus/continuous)
route <- "oral"                     # exposure route (oral/iv)
E_start <- 0                        # start of exposure phase (h)
E_end <- 3                          # end of exposure phase (h)
E_int <- 80                         # interval between doses (h); only relevant when regime="bolus" or iv

# Simulation parameters

A_type <- "VA"                # type of probabilistic analysis ("SA" or "VA")
chem_fix <- TRUE              # fixing the chemical and TK parameters or not (TRUE/FALSE)
n_boot <- 1000                # number of boostrap iterations (for SA)
n_out  <- 2                   # number of compartments to output (blood, total body)
t_start <- 0                  # start of simulation (h)
t_end <- 12                   # end of simulation (h)
t_A <- c(seq(0.05,0.2,by=0.05),seq(0.25,t_end,by=0.25))     # all time points for analysis (h), only relevant when chem_fix=TRUE


#################################################################################################################################################
# Set up

source("R_code/simrun.r")

#################################################################################################################################################
# Write

if (chem_fix) {sim_data=na.omit(t(rbind(SimRes,t_A)))
               data_blood=data.frame(t_A=sim_data[,n_sim+1],sim_data[,1:n_sim])
               row.names(data_blood)=NULL
               PKPDmisc::nca(as.numeric(data_blood$t_A), as.numeric(data_blood$X2)*1000, dose=E_dose, last_times = c(3, 4, 5), digits = 2)
               # file.remove(paste("./Output/SimRes_",idsim,".xls",sep=""))
               # writeWorksheetToFile(paste("./Output/SimRes_",idsim,".xls",sep=""), data_body, sheet=paste0("body_",idsim))
}
# saveRDS(res_sims,paste("./Output/SimRes_",idsim,".rds",sep=""))

##############################################################################################################################
# Plotting of results
# Create data frame with results

plotting_blood <- as.data.frame(matrix(NA,nrow=length(t_A),ncol=5))


for (i in 1:length(t_A)) {
      if (!all(is.na(SimRes[,i]))) {
        plotting_blood[i,1] <- min(SimRes[,i],na.rm = TRUE)
        plotting_blood[i,2] <- quantile(SimRes[,i],0.025, na.rm = TRUE)
        plotting_blood[i,3] <- quantile(SimRes[,i],0.5, na.rm = TRUE)
        plotting_blood[i,4] <- quantile(SimRes[,i],0.975, na.rm = TRUE)
        plotting_blood[i,5] <- max(SimRes[,i],na.rm=TRUE)
      }
    }

    coltemp <- colnames(SimRes[grep('blood',colnames(SimRes))])
    plotting_blood$t <- substr(coltemp,start=8,stop=nchar(coltemp))
    plotting_blood <- subset(plotting_blood,t%in%t_A)
    plotting_blood$t <- NULL
    
# i-100 need to be ajusted if the variables in SimRes are different

plotting_wb <- as.data.frame(matrix(NA,nrow=length(t_A),ncol=5))
  for (i in (length(t_A)+1):(2*length(t_A))) {
      if (!all(is.na(SimRes[,i]))) {
        plotting_wb[i-100,1] <- min(SimRes[,i],na.rm = TRUE)
        plotting_wb[i-100,2] <- quantile(SimRes[,i],0.025, na.rm = TRUE)
        plotting_wb[i-100,3] <- quantile(SimRes[,i],0.5, na.rm = TRUE)
        plotting_wb[i-100,4] <- quantile(SimRes[,i],0.975, na.rm = TRUE)
        plotting_wb[i-100,5] <- max(SimRes[,i],na.rm=TRUE)
      }
    }
coltemp <- colnames(SimRes[grep('whole',colnames(SimRes))])
plotting_wb$t <- substr(coltemp,start=12,stop=nchar(coltemp))
plotting_wb <- subset(plotting_wb,t%in%t_A)
plotting_wb$t <- NULL
 
   
################################################################################################################################
# Plots
    
# Plot blood concentration    
# set figure margins

par(mgp=c(1.75,0.75,0),mar=c(5.5,3.5,3,1),xpd=NA) 

    #create empty plot frame

    plot(x=NULL,y=NULL, xlim=c(0,max(t_A)), ylim=c(0,1.25*max(plotting_blood,na.rm=TRUE)),xlab="Time (hr)",
    ylab="Concentration in blood (mg/L)", frame.plot=FALSE, xaxt="n", yaxt="n", cex.lab=1.25)

    #create polygon for 95% CI

    polygon(c(0,t_A,rev(t_A),0), c(0,plotting_blood[,2],rev(c(0,plotting_blood[,4]))), col="gray", border=NA) 

    #create line representing median response

    lines(c(0,plotting_blood[,3])~c(0,t_A),type="l",xlab="Time (hr)",ylab="Concentration in blood (mg/L)")
    #add median points

    points(t_A,plotting_blood[,3],pch=16,cex=0.75,col="black") 

    #add min and max points
    #points(1:(nrow(plotting_blood)),plotting_blood[,1],pch=1,cex=0.75,col="black") 
    #points(1:(nrow(plotting_blood)),plotting_blood[,5],pch=1,cex=0.75,col="black")

    #add axes
    axis(1,at=c(0:max(t_A)),pos=0,cex.axis=1)
    axis(2,at=seq(0,ceiling(1300*max(plotting_blood, na.rm=TRUE))/1000,by=0.0001),
         pos=par("usr")[1]+(par("usr")[2]-par("usr")[1])/35,cex.axis=1)

    #add legend

    #legend(x="topright", legend=c("P50","min-max","95% interval"),pch=c(16,1,NA),lty=c(NA,NA,1),seg.len=0.3,lwd=c(1,1,2),col=c("black","black","darkgrey"),bty="n",text.col="black",cex=1, y.intersp=0.8,x.intersp=c(1,1,0.9),merge=F) 

    # Adds a legend box to the plot
    legend(x="topright", legend=c("P50","min-max","95% interval"),pch=c(16,1,NA),lty=c(NA,NA,1),seg.len=0.3,lwd=c(1,1,5),col=c("black","black","grey"),bty="n",text.col="black",cex=1, y.intersp=0.8,x.intersp=c(1,1,0.9),merge=F) # Adds a legend box to the plot
    dev.off()    
    
###############################################################
# Plot whole body concentration    
# set figure margins

par(mgp=c(1.75,0.75,0),mar=c(5.5,3.5,3,1),xpd=NA)
# create empty plot frame
    plot(x=NULL,y=NULL, xlim=c(0,max(t_A)), ylim=c(0,1.25*max(plotting_wb,na.rm=TRUE)),xlab="Time (hr)", ylab="Concentration in whole body (mg/L)", frame.plot=FALSE, xaxt="n", yaxt="n", cex.lab=1.25)
    #create polygon for 95% CI
    polygon(c(0,t_A,rev(t_A),0), c(0,plotting_wb[,2],rev(c(0,plotting_wb[,4]))), col="gray", border=NA) 
    #create line representing median response
    lines(c(0,plotting_wb[,3])~c(0,t_A),type="l")
    #add median points
    points(t_A,plotting_wb[,3],pch=16,cex=0.75,col="black") 
    #add min and max points
    #points(1:(nrow(plotting_wb)),plotting_wb[,1],pch=1,cex=0.75,col="black") 
    #points(1:(nrow(plotting_wb)),plotting_wb[,5],pch=1,cex=0.75,col="black")
    #add axes
    axis(1,at=c(0:max(t_A)),pos=0,cex.axis=1)
    axis(2,at=seq(0,ceiling(1300*max(plotting_wb,na.rm=TRUE))/1000,by=0.01),pos=par("usr")[1]+(par("usr")[2]-par("usr")[1])/35,cex.axis=1)
    # add legend
    # legend(x="topright", legend=c("P50","min-max","95% interval"),pch=c(16,1,NA),lty=c(NA,NA,1),seg.len=0.3,lwd=c(1,1,2),col=c("black","black","darkgrey"),bty="n",text.col="black",cex=1, y.intersp=0.8,x.intersp=c(1,1,0.9),merge=F) # Adds a legend box to the plot
    legend(x="topright", legend=c("P50","min-max","95% interval"),pch=c(16,1,NA),lty=c(NA,NA,1),seg.len=0.3,lwd=c(1,1,5),col=c("black","black","grey"),bty="n",text.col="black",cex=1, y.intersp=0.8,x.intersp=c(1,1,0.9),merge=F) # Adds a legend box to the plot
dev.off()
    
###########################################################################################################################    
# Sobol analysis results and plots, including lowry plots 
    
###############################################################
# Sobol analysis plot blood


if (chem_fix & A_type=="SA") {
  
  
par(mfrow=c(2,2), las=2, cex=0.7)
FOI = TI = TI.borninf= TI.bornsup= matrix(NA, nrow = NP_var, ncol = length(SimRes)) 
rownames(FOI)= rownames(TI)= rownames(TI.borninf) = rownames(TI.bornsup) = Names_var
t_SA <- c(SimRes$AUC_h)
    


for(i in 1:length(SimRes)){
      print(i)
      if (SimRes$AUC_h[i] %in% t_SA) {
        tell(x = sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
        FOI[,i]       = sa$S[,1]    #First order indices
        TI[,i]        = sa$T[,1]    #Total indices
        TI.borninf[,i] = sa$T[,4]   #Lower CL total indices 
        TI.bornsup[,i] = sa$T[,5]   #Upper CL total indices
        
        plot(sa, main=colnames(SimRes)[i])  
      }
      
}

dev.off()
    
############################################################################################
# Sobol analysis whole body

par(mfrow=c(2,2), las=2, cex=0.7)

FOI = TI = TI.borninf= TI.bornsup= matrix(NA, nrow = NP_var, ncol = length(t_A))

rownames(FOI)= rownames(TI)= rownames(TI.borninf) = rownames(TI.bornsup) = Names_var

t_SA <- c(0.2,1.5,4,12)
    
for(i in 1:length(t_A)) {
      print(i)
      if (t_A[i] %in% t_SA) {
        tell(x = sa, y = SimRes[,length(t_A)+i], nboot = n_boot, conf = 0.95)
        FOI[,i]       = sa$S[,1]    #First order indices
        TI[,i]        = sa$T[,1]    #Total indices
        TI.borninf[,i] = sa$T[,4]   #Lower CL total indices 
        TI.bornsup[,i] = sa$T[,5]   #Upper CL total indices
        
        plot(sa, main=colnames(SimRes)[i])  
      }
      
    }

dev.off()
 
#############################################################################################   
# Lowry plots

i <- 20  # repeat for each point in time (column of t_A)


# sum of first order variances should not exceed 1

    for (j in 1:nrow(FOI)) {
      if (FOI[j,i]<0) {
        FOI[j,i] <- 0
      }
    }
    
    if (sum(FOI[,i])>1) {
      FOI[,i] <- FOI[,i]/sum(FOI[,i])
    } 

# Create data frame for lowry, with first order and interaction effects

lowry <- data.frame(Parameter=Names_var,Interaction=TI[,i]-FOI[,i],Main.Effect=FOI[,i])

# Make sure that interaction indices are not below 0

lowry$Interaction <- ifelse(lowry$Interaction<0,0,lowry$Interaction)
    
lowry_plot(lowry)

}

#########################################################################################################à   


