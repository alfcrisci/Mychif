library(foreach)
library(doSNOW)
library(parallel)

# Use the detectCores() function to find the number of cores in system
#no_cores <- detectCores()



multi_tool_blood <- function(par_in,
                       species, 
                       regime, 
                       route, 
                       E_dose, 
                       E_start, 
                       E_end, 
                       E_int,  
                       n_out, 
                       t_start, 
                       t_end, 
                       t_A, 
                       chem_fix,
                       timestep=3) {
 
  ############################################################################################################################
  # General data of integration from minute -1 data.
  require(parallel)
  require(foreach)
  cl <- parallel::makeCluster(detectCores()-1)
  registerDoParallel(cl)
  registerDoSEQ()
  
  timeint=floor(60/timestep) # time resolution
  timescale=floor(60*60/timestep) # steps of integration
  
  ############################################################################################################################
  # General physiology
  
  BW <- par_in$BW            # body weight (kg)
  CO <- par_in$CO/timeint         # cardiac output (L/3sec)
  
  # source data are in min -1.
  
  Qexhale <- (0.499 * (BW^0.81))/timeint        # Ventilation mammalian (L/3sec)
 
  if (species=="chicken") { Qexhale <- ((284 * (BW^0.77))/1000)/timeint}    # Ventilation avian (L/3sec)
  
  # all compartments in one data frame

  comp <- colnames(par_in[grep('fBW',colnames(par_in))])
  comp <- comp[order(comp)]
  phys <- data.frame(comp=c(substr(comp,start=5,stop=nchar(comp)),"ven","art","lumen","milk","egg","urine","air","feces","metab","feed","iv"),fBW=0, fCO=0)
  phys$comp <- as.character(phys$comp)
  phys$fBW <- c(t(par_in[comp]),rep(0,11))
  comp <- gsub('fBW','fCO',comp)
  phys$fCO <- c(t(par_in[comp]),rep(0,11))
  phys$fBW[phys$comp=="ven"] <- 2/3 * phys$fBW[phys$comp=="blood"] # two third of total blood volume assumed venous
  phys$fBW[phys$comp=="art"] <- 1/3 * phys$fBW[phys$comp=="blood"] # one third of total blood volume assumed venous
  
  # normalisation of fractions per tissue

  PCtemp <- 0
  PCsum <- NA
  PCsumnew <- NA
  PCcol <- c(grep('_nl',colnames(par_in)),grep('_pl',colnames(par_in)),grep('_pr',colnames(par_in)),grep('_w',colnames(par_in)))
  fPC <- par_in[PCcol[order(PCcol)]]
  
  # sum of mean fractions per tissue

  for (i in 9:ncol(fPC)) {
    PCtemp <- PCtemp +fPC[,i]
    if ((i-8)%%4==0) {
      PCsum[(i-8)/4] <- PCtemp    
      PCtemp <- 0
    }
  }
  
  # new fractions

  fPC[c(9:56)] <- fPC[c(9:56)]*c(par_in$nl,par_in$pl,par_in$pr,par_in$w)
  for (i in 9:ncol(fPC)) {
    PCtemp <- PCtemp +fPC[,i]
    if ((i-8)%%4==0) {
      PCsumnew[(i-8)/4] <- PCtemp    
      PCtemp <- 0
    }
  }
  for (i in 9:ncol(fPC)) {
    fPC[,i] <- fPC[,i]*(PCsum[ceiling((i-8)/4)]/PCsumnew[ceiling((i-8)/4)])
  }
 
  # partition coefficients
  
 for (i in 1: nrow(phys)) {
    if (phys$comp[i]=="air") {
      #solubility in mol/L instead of mg/L
      S <- (0.001*par_in$S) / par_in$MW 
      Temp <- 298       # Temperature 298 K = 25 degC
      R <- 8.314        # Gas constant (J/mol/K)
      phys$PC[i] <- (par_in$Pv / S) / (R * Temp)  #Kaw (-)

    } else if (phys$comp[i]=="stomach"|
               phys$comp[i]=="abomasum"|
               phys$comp[i]=="omasum"|
               phys$comp[i]=="rumen"|
               phys$comp[i]=="reticulum"|
               phys$comp[i]=="gizzard"|
               phys$comp[i]=="crop") { 

      # partitioning coefficient tissue-water per tissue, based on QSAR from Hendriks et al (2005) (-)
      # tissue composition assumed the same as intestine for these compartments

      phys$PC[i] <- fPC$int_nl * fPC$intestine_nl * par_in$Kow ^ fPC$exp_nl + 
        fPC$int_pl * fPC$intestine_pl * par_in$Kow ^ fPC$exp_pl  + 
        fPC$int_pr * fPC$intestine_pr * par_in$Kow ^ fPC$exp_pr  +
        fPC$int_w * fPC$intestine_w * par_in$Kow ^ fPC$exp_w
      
    } else if (phys$fCO[i]==0|is.na(phys$fCO[i])) {

   # partitioning coefficient tissue-water per tissue, based on QSAR from Hendriks et al (2005) (-)
   # tissue composition same as blood for venous, arterial blood, and 'external' compartments

      phys$PC[i] <- fPC$int_nl * fPC$blood_nl * par_in$Kow ^ fPC$exp_nl + 
        fPC$int_pl * fPC$blood_pl * par_in$Kow ^ fPC$exp_pl  + 
        fPC$int_pr * fPC$blood_pr * par_in$Kow ^ fPC$exp_pr +
        fPC$int_w * fPC$blood_w * par_in$Kow ^ fPC$exp_w
    } else {
      
    # partitioning coefficient tissue-water per tissue, based on QSAR from Hendriks et al (2005) (-)

      phys$PC[i] <- fPC$int_nl * fPC[,paste(phys$comp[i],"nl",sep="_")] * par_in$Kow ^ fPC$exp_nl + 
        fPC$int_pl * fPC[,paste(phys$comp[i],"pl",sep="_")] * par_in$Kow ^ fPC$exp_pl  + 
        fPC$int_pr * fPC[,paste(phys$comp[i],"pr",sep="_")] * par_in$Kow ^ fPC$exp_pr  +
        fPC$int_w * fPC[,paste(phys$comp[i],"w",sep="_")] * par_in$Kow ^ fPC$exp_w
    }
  }
  
  phys$V <- phys$fBW*BW/sum(phys$fBW)                       # volumes (L), based on normalized weight fractions and overall assumed density of 1 L/kg
  phys$Q <- phys$fCO*CO/sum(phys$fCO,na.rm=TRUE)            # blood flows (L/sec), based on normalized CO fractions
  phys$PC <- phys$PC/phys$PC[phys$comp=="blood"]            # tissue-blood partitioning coefficients (-)
  phys <- subset(phys,comp!="blood")                        # remove total blood as compartment
  
  # rates and flows (not chemical-dependent)

  phys$Qingest <- ifelse(phys$comp=="lumen",par_in$Qingest,0)                          # ingestion rate (kgfeed/d)
  phys$Qingest <- phys$Qingest/(24*timescale)                                          # transform feeding rates to kgfeed/10sec
  phys$kgastric <- ifelse(phys$comp=="lumen",par_in$kgastric,0)                        # gastric emptying rate constant (1/min)
  phys$kgastric <- phys$kgastric/timeint                                               # transform gastric emptying rate constant to 1/10sec
  phys$fbile <- ifelse(phys$comp=="lumen",par_in$fbile,0)                              # fraction of metabolites reentering lumen with bile (-)  
  phys$Qmilk <- ifelse(phys$comp=="mamgland"&phys$fCO!=0,par_in$Qmilk,0)/timeint       # milk production rate (L/10sec)
  phys$Qegg <- ifelse(phys$comp=="reprod"&phys$fCO!=0,par_in$Qegg,0)/timeint           # egg production rate (L/10sec)
  
  # toxicokinetic parameters

  phys$kabs <- ifelse(phys$comp=="intestine",par_in$kabs,0)/timeint                                              # absorption rate constant from lumen into intestine (1/10sec)
  phys$Vmax <- ifelse(phys$comp=="liver",par_in$Vmax_tot,NA)/timeint                                             # Vmax (mmol/10sec) (only metabolism in liver assumed)
  phys$Km <- ifelse(phys$comp=="liver",par_in$Km_tot,NA)                                                         # Km (mmol) (only metabolism in liver assumed)
  phys$Cl <- ifelse(phys$comp=="liver",par_in$Cl_hepatic,ifelse(phys$comp=="kidney",par_in$Cl_renal,0))/timeint  # Clearance (L/10sec/kg) for liver and kidney
  phys$Cl <- phys$Cl * BW                                                                                        # transform clearance rates to L/10sec
  phys$fbact <- ifelse(phys$comp=="lumen",par_in$fbact,0)                                                        # fraction of metabolites retransformed into parent compound in lumen (-)
  
  # Create output sheet ----
  


###############################################################################################################################################
# Model simulation 
  
simsteps=as.integer((timescale*t_start):(timescale*t_end))

phys$M <- 0 # mass per compartment (mmol)
phys$C <- 0 # concentration per compartment (mmol/L)


aa=foreach (i=as.integer((timescale*t_start):(timescale*t_end)),.combine = c)  %do%  {
  
  
  # reset dM_in & dM_out for all compartments

  phys$dM_in <- 0
  phys$dM_out <- 0
  
  # ABSORPTION

  # dM/dt as result of feeding

  if (route=="oral") {
      phys$dM_iv <- 0
    
      if (i<(timescale*E_start) | i>(timescale*E_end)) {
          phys$dM_food <- 0 # outside of exposure period
    } else if (regime=="continuous") {
           phys$dM_food <- (E_dose/par_in$MW)*par_in$Qingest # mmol/sec
    } else if (i%%(timescale*E_int)==0) {
      phys$dM_food <- (E_dose/par_in$MW) # one immediate bolus dose (mmol)  (if dose = mg/kgbw, multiply with BW)  
      phys$dM_food[phys$comp!="lumen"] <- 0
    } else {
      phys$dM_food <- 0 # no bolus dose
    }
  }
  
  # dM/dt as result of iv

  if (route=="iv") {
    phys$dM_food <- 0
    if (i<(timescale*E_start)|i>(timescale*E_end)) {  
      phys$dM_iv <- 0
    } else if (i%%(timescale*E_int)==0) {
      phys$dM_iv <- (E_dose/par_in$MW)*BW    # iv dose in mmol (mg/kg / mg/mmol * kg)
      phys$dM_iv[phys$comp!="ven"] <- 0
    } else {
      phys$dM_iv <- 0
    }
  }
  ####################################################################################################
  phys$dM_in <- phys$dM_in + phys$dM_food
  
  phys$dM_out[phys$comp=="feed"] <- phys$dM_food[phys$comp=="lumen"]
  
  phys$dM_in <- phys$dM_in + phys$dM_iv
  
  phys$dM_out[phys$comp=="iv"] <- sum(phys$dM_iv)
  
  # dM/dt as result of absorption over intestinal wall

  phys$dM_abs <- phys$M[phys$comp=="lumen"]*phys$kabs
  phys$dM_in <- phys$dM_in + phys$dM_abs
  phys$dM_out[phys$comp=="lumen"] <- phys$dM_out[phys$comp=="lumen"] + sum(phys$dM_abs)
  
  # DISTRIBUTION

  # dM/dt as result of delivery via arterial blood

  phys$dM_art <- phys$Q*(phys$C[phys$comp=="art"])
  phys$dM_in <- phys$dM_in + phys$dM_art
  phys$dM_out[phys$comp=="art"] <- phys$dM_out[phys$comp=="art"]+sum(phys$dM_art)
  
  # dM/dt as result of delivery to venous blood

  phys$dM_ven <- phys$Q*(phys$C/phys$PC)
  
  # portal vein

  phys$dM_in[phys$comp=="liver"] <- phys$dM_in[phys$comp=="liver"] + phys$dM_ven[phys$comp=="intestine"]
  phys$dM_ven[phys$comp=="liver"] <- phys$dM_ven[phys$comp=="liver"] + phys$Q[phys$comp=="intestine"]*(phys$C[phys$comp=="liver"]/phys$PC[phys$comp=="liver"])
  
  phys$dM_out <- phys$dM_out + phys$dM_ven
  phys$dM_in[phys$comp=="ven"] <- phys$dM_in[phys$comp=="ven"] + sum(phys$dM_ven[phys$comp!="intestine"])
  
  # METABOLISM

  # dM/dt as result of metabolism

  phys$dM_met <- ifelse(!is.na(phys$Vmax)&!is.na(phys$Km),(phys$Vmax*(phys$C/phys$PC))/(phys$Km+phys$C/phys$PC),0)
  phys$dM_met[phys$comp=="liver"] <-ifelse(!is.na(phys$Cl[phys$comp=="liver"]),phys$Cl[phys$comp=="liver"]*phys$C[phys$comp=="liver"],0)
  
  phys$dM_out <- phys$dM_out + phys$dM_met
  
  phys$dM_in <- phys$dM_in + phys$dM_met[phys$comp=="liver"]*phys$fbile*phys$fbact
  phys$dM_in[phys$comp=="metab"] <- phys$dM_in[phys$comp=="metab"] + sum(phys$dM_met) - phys$dM_met[phys$comp=="liver"]*(1-    phys$fbile[phys$comp=="lumen"]*phys$fbact[phys$comp=="lumen"])   
  
  # ELIMINATION

  # dM/dt as result of transport over lung and exhalation

  phys$dM_exh <- ifelse(phys$comp=="ven", CO*phys$C[phys$comp=="ven"],0) # all mass in venous blood is transported
  phys$dM_out <- phys$dM_out + phys$dM_exh
  phys$dM_in[phys$comp=="art"] <- phys$dM_in[phys$comp=="art"] + phys$dM_exh[phys$comp=="ven"] * (CO/(CO+Qexhale*phys$PC[phys$comp=="air"]))
  phys$dM_in[phys$comp=="air"] <- phys$dM_in[phys$comp=="air"] + phys$dM_exh[phys$comp=="ven"] * ((Qexhale*phys$PC[phys$comp=="air"])/(CO+Qexhale*phys$PC[phys$comp=="air"]))
  
  # dM/dt as result of excretion to urine, milk, eggs, feces

  phys$dM_exc <- ifelse(phys$comp=="mamgland",phys$C[phys$comp=="mamgland"]*phys$Qmilk[phys$comp=="mamgland"],
                        ifelse(phys$comp=="kidney",phys$C[phys$comp=="kidney"]*phys$Cl[phys$comp=="kidney"],
                               ifelse(phys$comp=="reprod",phys$C[phys$comp=="reprod"]*phys$Qegg[phys$comp=="reprod"],
                               phys$M*phys$kgastric)))
  
  phys$dM_out <- phys$dM_out + phys$dM_exc
  
  phys$dM_in[phys$comp=="urine"] <- phys$dM_in[phys$comp=="urine"] + phys$dM_exc[phys$comp=="kidney"]
  phys$dM_in[phys$comp=="milk"] <- phys$dM_in[phys$comp=="milk"] + phys$dM_exc[phys$comp=="mamgland"]
  phys$dM_in[phys$comp=="egg"] <- phys$dM_in[phys$comp=="egg"] + phys$dM_exc[phys$comp=="reprod"]
  phys$dM_in[phys$comp=="feces"] <- phys$dM_in[phys$comp=="feces"] + phys$dM_exc[phys$comp=="lumen"]
  
  # Mass and concentration per compartment

  phys$M <- phys$M + phys$dM_in - phys$dM_out  # mmol
  phys$C <- ifelse(phys$V==0,0,phys$M/phys$V)  # mmol/L
  
  # Write to results for t_A
  as.numeric(phys$C[phys$comp=="ven"]*par_in$M)
    
}
parallel::stopCluster(cl)

# End simulation 
###############################################################################################################################################
parallel::stopCluster(cl)

results=data.frame(timesec=simsteps*timestep,blood_c=aa)
return(results)
}


###############################################################################################################################



