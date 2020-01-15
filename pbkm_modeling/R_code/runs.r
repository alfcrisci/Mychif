library(fitdistrplus)

######################################################################################

library(httk) # High-Throughput Toxicokinetics 
library(sensitivity) # Sensitivity Analysis 

######################################################################################

# devtools::install_github("nanhung/pksensi")

library(pksensi) # Global Sensitivity Analysis in Pharmacokinetic Modeling Hsieh 2018

library(PK) # Basic Non-Compartmental Pharmacokinetics https://cran.r-project.org/web/packages/PK/index.html

#devtools::install_github("dpastoor/PKPDmisc")

library(PKPDmisc) # Pharmacokinetic and Pharmacodynamic Data Management Functions

######################################################################################
# http://pkpdsim.ronkeizer.com/ https://github.com/InsightRX/PKPDsim

# devtools::install_github("InsightRX/PKPDsim")

library(PKPDsim) # Simulate dose regimens for PKPD models described by ODE or equation 

#devtools::install_github("ronkeizer/PKPDsimShiny")

library(PKPDsimShiny) #  

#######################################################################################################àà

source('https://alfcrisci.github.io/Mychif/plot_senstivity_aux.r')


#######################################################################################################àà
# Distribution analisys

######################################################################################################################################
# check data in httk sources

mycotoxins_data=read.csv("https://alfcrisci.github.io/Mychif/mycotoxins.csv")

mycoCAS <- mycotoxins_data$CAS
chem.lists$Tox21[which(chem.lists$Tox21$CAS %in% mycoCAS ==T),] # Zearalenone
chem.lists$ExpoCast[which(chem.lists$ExpoCast$CAS %in% mycoCAS ==T),] # Zearalenone ZEA 807, Aflatoxin B1  6469
chem.physical_and_invitro.data[which(chem.physical_and_invitro.data$CAS %in% mycoCAS ==T),] # OTA, Aflatoxin B1


#############################################################################################################################################
# Non compartimental analisys
#############################################################################################################################################
# set seeds for reproducibility

set.seed(34534)

#############################################################################################################################################

# Considering 1 comp Flip-flop kinetics (FFPK) PK Model for Sensitivity Analysis by pksensi R packages

# Inputs are generally considered for bolus regimes:

# F bioavailability :the fraction or percentage of the administrated dose that can reach the general circulation
# k_a is the first-order absorption rate constant (1/time)
# k_e is the first-order elimination rate constant (/time), 
# Vd is the distribution volume. see https://en.wikipedia.org/wiki/Volume_of_distribution
# Volume of distribution is called a “primary pharmacokinetic parameter”, which means that this parameter depends on the physiologic properties of the body and the physiochemical properties of the drug. Volume of distribution is not derived from other PK parameters, instead it is used to estimate the “secondary” PK parameters
# ndr when regime is bolus Vd/F equal 2.01 Vd=2.01*0.75

#################################################################################
#################################################################################
# toxicokinetics
#################################################################################

# po Saint Cyr et al 2015 doi:10.3390/toxins7124873

D=100 # microL/kgbw
params <- c(F = 0.75, KA = 3.72, KE = 0.29, V = 1.5074)  
t <- seq(0, 24, 0.01)
C <-pksensi::FFPK(params = params, time = t,dose=D)
plot(t, C, type = "l", xlab = "Hours", ylab = "DON concentration")
PKPDmisc::nca(t, C, dose=100, last_times = c(3, 4, 5), digits = 2)


####################################################################################
# po NDON control Paulick et al 2015 doi:10.3390/toxins7114622

D=75 # microL/kgbw
params <- c(F = 0.98, KA = 0.61, KE = 0.12, V = 1.68)  
t <- seq(0, 24, 0.01)
C <-FFPK(params = params, time = t,dose=D)
plot(t, C, type = "l", xlab = "Hours", ylab = "DON concentration")
PKPDmisc::nca(t, C, dose=75, last_times = c(3, 4, 5), digits = 2)

# sulfite DON 

D=71 # microL/kgbw
params <- c(F = 0.893, KA = 0.99, KE = 0.13, V = 1.98)  
t <- seq(0, 24, 0.01)
C <-FFPK(params = params, time = t,dose=D)
plot(t, C, type = "l", xlab = "Hours", ylab = "DON concentration")
PKPDmisc::nca(t, C, dose=71, last_times = c(3, 4, 5), digits = 2)

########################################################################################################
# Sensibility DON pigs

# params <- c(F = 0.893, KA = 0.99, KE = 0.13, V = 1.98)  
# params <- c(F = 0.75, KA = 3.72, KE = 0.29, V = 1.5074)  

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
plot(x)

##################################################################################################
# sensitivity analysis

q.arg <- list(list(min = 0.6, max = 1.0),
              list(min = 0.5, max = 4),
              list(min = 0.02, max = 0.3),
              list(min = 1, max = 3))


sobolSim <- function(r,n){
  M <- length(r)
  X <- data.frame(matrix(runif(M * n), nrow = n))
  for (m in (1:M)){
    rm <- r[[m]]
    X[,m] <- X[,m]*(rm$max-rm$min) + rm$min
  }
  return(X)
}

X1 <- sobolSim(q.arg,n=200)
X2 <- sobolSim(q.arg,n=200)

dose=100
time=seq(0.1, 24, 0.1)
n_boot=200

FFPKTmax=function(X) {
  X=as.numeric(X)
  params <- c(F = X[1], KA = X[2], KE =  X[4], V =X[4])
  as.numeric(time[which.max(as.numeric(FFPK(params,time=time,dose=dose)))])
}

FFPKCmax=function(X) {
  X=as.numeric(X)
  params <- c(F = X[1], KA = X[2], KE =  X[4], V =X[4])
  max(as.numeric(FFPK(params,time=time,dose=dose)))
}

FFPKauc=function(X) {
  X=as.numeric(X)
  params <- c(F = X[1], KA = X[2], KE =  X[4], V =X[4])
  sum(as.numeric(FFPK(params,time=time,dose=dose)))
}

sa <- soboljansen(model = NULL, X1, X2, nboot = n_boot, conf = 0.95)

 <- matrix(NA, nrow = nrow(sa$X), ncol = 3)
SimRes[,1]=apply(sa$X,1,FFPKTmax)
SimRes[,2]=apply(sa$X,1,FFPKCmax)
SimRes[,3]=apply(sa$X,1,FFPKauc)

colnames(SimRes) <- c("Tmax","Cmax","AUC")
tell(x = sa, y = SimRes, nboot = n_boot, conf = 0.95)

#################################################################################################################

par_sim=c("Tmax","Cmax","AUC")
n=length(par_sim)
np=length(c("F","KA","KE","V"))

FOI = TI = TI.borninf= TI.bornsup= matrix(NA, nrow = np, ncol = n)
rownames(FOI)= rownames(TI)= rownames(TI.borninf) = rownames(TI.bornsup) =c("F","Kabs","Ke","Vd")
par(mfrow=c(2,2), las=2, cex=0.7)

for(i in 1:n){
     print(i)
     tell(x = sa, y = SimRes[,i], nboot = n_boot, conf = 0.95)
  
     FOI[,i]       = sa$S[,1]    #First order indices
     TI[,i]        = sa$T[,1]    #Total indices
     TI.borninf[,i] = sa$T[,4]   #Lower CL total indices 
     TI.bornsup[,i] = sa$T[,5]   #Upper CL total indices
    
     plot(sa, main=as.character(par_sim[i])) 
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


########################################################################################
  





#############################################################################
# example with PKsim

p <- list(CL = 0.52, V  = 2.01*0.75, KA = 3.72)

pk <- new_ode_model("pk_1cmt_oral")
r1 <- new_regimen(amt = 100*0.75,
                  n = 1,
                  interval = 24)

dat <- sim_ode(
  ode = pk,
  parameters = p,
  regimen = r1,
  t_obs = seq(0, 24, 0.01)
)

PKPDmisc::nca(dat[dat$comp=="obs",]$t, dat[dat$comp=="obs",]$y, dose=100, last_times = c(3, 4, 5), digits = 2)
#############################################################################


#############################################################################
# install.packages("rsconnect")
library(PKPDsim)
library(PKPDsimShiny)
library(shinythemes)


p <- list(CL = 1, V  = 10, KA = 0.5)
pk1 <- new_ode_model("pk_1cmt_oral")
r1 <- new_regimen(amt = 100,
n = 5,
interval = 12)
n <- 100
 omega <- c(0.1,
0.05, 0.1)
#simulation in a Shiny app
PKPDsimShiny::sim_ode_shiny (
 ode = pk1,
 par = p,
 n_ind = 100,
 regimen = r1)
shiny::runApp("shiny-pkpd")
rsconnect::setAccountInfo(name='alfcrisci', token='B92FF3D98C78A375B3156FE4F27C8FF9', secret='ACX1o0PvYRHpgtylngX9+c00tYoUuyYsthDpVE4y')
rsconnect::deployApp("shiny-pkpd",appName = "mychifonecomp",appTitle = "mychif_oneC")

#shinytheme("united")
nca_res=PK::nca(conc, time, n.tail=3, dose=100, method=c("z", "boott"), 
                conf.level=0.95, design="complete",nsample=1000, d.conc)

d.conc=data.frame(conc=C, time=t,subject=1)
conc.obj <- PKNCA::PKNCAconc(d.conc, conc~time|subject)
plot(conc.obj)
