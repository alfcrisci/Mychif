#===================================#
# Probabilistic                     #
# multi-compartment model for       #
# veterinary species                #
#===================================#
# Load codes and libraries


library(sensitivity)
library(reshape)
library(ggplot2)
library(compiler)
library(XLConnect)

################################################################################################################################################################
# Set working directory

setwd("/home/alf/Scrivania/lav_michyf/parma_meetings/pbkm_animals")

source("R_code/PBTK_model_core.R")

multi_tool_fast=cmpfun(multi_tool)

source('R_code/Lowry_plots.R')

# cl <- makeCluster(4, type="SOCK") # for 4 cores machine
# registerDoSNOW (cl)
################################################################################################################################################################
# Inputs 
############################################################################################################################
# Loading TK data

TK <- read.csv("data/tk_data/TK_prob.csv",stringsAsFactors=FALSE)                  #toxicokinetic parameters

############################################################################################################################
# Loading physicochemical data

chem <- read.csv("data/chem_data_prob.csv", stringsAsFactors = FALSE)      #physicochemical parameters

############################################################################################################################
# Loading physiology data

fBW <- read.csv("data/physio/fBW_animals_prob.csv",stringsAsFactors = FALSE)      # organ fractions, incl var-distribution parameters
fCO <- read.csv("data/physio/fCO_animals_prob.csv",stringsAsFactors = FALSE)      # blood flow fractions, incl var-distribution parameters
rates <- read.csv("data/physio/rates_animals_prob.csv",stringsAsFactors = FALSE)  # physiological rates (not chemical-dependent)
fPC <- read.csv("data/physio/PC_tissue_prob.csv",stringsAsFactors = FALSE)        # tissue composition, fixed during var-analysis

############################################################################################################################
# Defining species and toxin

idsim="Swine_DON"
species <- "swine"                      # cat/cattle/chicken/sheep/swine
chemical <- "deoxynivalenol"            # for now: chlorpyrifos/melamine/meloxicam/monensin/oxytetracycline/PFOS

#################################################################################################################
# Customize population parameter
# Saint-Cyr

BW_animal=27                            # kg
BW_SD=1                                 # kg standard deviation 


#################################################################################################################
# Paulinski
# BW_animal=49                            # kg
# BW_SD=3.5                               # kg standard deviation

fBW$BW[which(fBW$species==species)]=BW_animal
fBW$BW_sd[which(fBW$species==species)]=BW_SD


############################################################################################################################
# Defining the exposure scenario

regime <- "bolus"                   # exposure regime (bolus/continuous)
route <- "oral"                     # exposure route (oral/iv)
E_dose <-1800                       # 2250 bolus dose 2700/1.2kg Saint Cyr in 2 volte
                                    # exposure dose (mg/kgfeed at continuous; mg at bolus; mg/kgbw at iv)
E_start <- 0                        # start of exposure phase (h)
E_end <- 0.9                          # end of exposure phase (h)
E_int <- 1                          # interval between doses (h); only relevant when regime="bolus" or iv

# estart 0.05 e end 24ore 
############################################################################################################################
# Simulation parameters

A_type <- "VA"                # type of probabilistic analysis ("SA" or "VA")
chem_fix <- TRUE              # fixing the chemical and TK parameters or not (TRUE/FALSE)
n_sim  <- 2                   # number of iterations 1000
n_boot <- 100                 # number of boostrap iterations (for SA) 1000
n_out  <- 2                   # number of compartments to output (blood, total body)



############################################################################################################################
# Simulation timing

t_start <- 0                  # start of simulation (h)
t_end <- 4                  # end of simulation (h)
t_A <- c(seq(0,t_end,by=0.25))     # all time points for analysis (h), only relevant when chem_fix=TRUE


#################################################################################################################################################
# Set up

res_sims=list()

if (chem_fix & A_type=="SA") {

    # Sensitivity analysis with fixed chem parameters and varying physiology  
    # initial list with names of all varying parameters
   
    Names_var <- c("BW",
                 paste("fBW",colnames(fBW[2+3*c(1:19)]),sep="_"),
                 "CO",
                 paste("fCO",colnames(fCO[2+3*c(1:19)]),sep="_"),
                 colnames(rates[2+3*c(0:4)]),
                 "nl","pl","pr","w")
    
    #initial list with median values of all varying parameters 
 
    Medians <- cbind(fBW[fBW$species==species,2+3*c(0:19)],
                     fCO[fCO$species==species,2+3*c(0:19)],
                     rates[rates$species==species,2+3*c(0:4)],
                     data.frame(nl=1,pl=1,pr=1,w=1))
    colnames(Medians) <- Names_var
    Names_var_not <- colnames(Medians[which(Medians==0|is.na(Medians))])
    
    #update distribution parameters for varying parameters

    Medians <- Medians[!(colnames(Medians)%in%Names_var_not)] 

    #exclude all compartments for which median = 0 or NA

    Lower <- Medians - 0.1*Medians
    Upper <- Medians + 0.1*Medians
    
    Names_var <- colnames(Medians)
    NP_var <- length(Names_var) #number of varying parameters
    
    # all fixed parameters (including 'varying parameters' that are 0 or NA)

    Names_fix <- c(colnames(TK[3:8]),
                   colnames(chem[2:5]),
                   colnames(fPC[2:57]),
                   Names_var_not)

    NP_fix <- length(Names_fix) # number of fixed parameters
    
    fix_in <- cbind(TK[TK$species==species&TK$chemical==chemical,colnames(TK)%in%Names_fix],
                    chem[chem$chemical==chemical,colnames(chem)%in%Names_fix],
                    fPC[fPC$species==species,colnames(fPC)%in%Names_fix],
                    fBW[fBW$species==species,colnames(fBW)%in%substr(Names_var_not[grep('fBW',Names_var_not)],start=5,stop=nchar(Names_var_not[grep('fBW',Names_var_not)]))],
                    fCO[fCO$species==species,colnames(fCO)%in%substr(Names_var_not[grep('fCO',Names_var_not)],start=5,stop=nchar(Names_var_not[grep('fCO',Names_var_not)]))],
                    rates[rates$species==species,colnames(rates)%in%Names_var_not])
    colnames(fix_in) <- Names_fix                 
   
 ######################################################
 # Create data frames with random samples


    X1 <- matrix(NA, nrow = n_sim, ncol = NP_var)
    colnames(X1) <- Names_var
    X1 <- as.data.frame(X1)
    
    X2 <- matrix(NA, nrow = n_sim, ncol = NP_var)
    colnames(X2) <- colnames(X1)
    X2 <- as.data.frame(X2)
    
    for(i in 1:NP_var){
      X1[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i]) 
      X2[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
    } 
    
 ######################################################
 # Sobol design
    sa <- soboljansen(model = NULL, X1, X2, nboot = n_boot, conf = 0.95)
    SARes <- sa$X
    
    #create output data frame
    #SimRes <- matrix(NA, nrow = nrow(SARes), ncol = n_out*length(t_A))
    #colnames(SimRes) <- c(paste("blood_t",t_A,sep=""),paste("wholebody_t",t_A,sep=""))
    #SimRes <- as.data.frame(SimRes)  
        
  } else if (A_type=="SA") {

    # Sensitivity analysis with fixed chem parameters and varying physiology
    # initial list with names of all varying parameters

    Names_var <- c(colnames(TK[3:8]),
                   colnames(chem[2:5]))
    
    #initial lower and upper limits for uniform distributions
    Lower <- as.data.frame(cbind(0,0,0,NA,NA,NA,10,0.001,0.00001,0.01))
    Upper <- as.data.frame(cbind(1,1,1,NA,NA,NA,1000,1000,10000000000000,100000))
    colnames(Lower) = colnames(Upper) = Names_var
    
    Names_var_not <- colnames(Lower[which(is.na(Lower))])
    
    Lower <- Lower[!(colnames(Lower)%in%Names_var_not)] #exclude all compartments for which median = 0 or NA
    Upper <- Upper[!(colnames(Upper)%in%Names_var_not)] #exclude all compartments for which median = 0 or NA
    
    Names_var <- colnames(Lower)
    NP_var <- length(Names_var) #number of varying parameters
    
    #list with names of all fixed parameters (including 'varying parameters' that are 0 or NA)

    Names_fix <- c("BW",
                   paste("fBW",colnames(fBW[2+3*c(1:19)]),sep="_"),
                   "CO",
                   paste("fCO",colnames(fCO[2+3*c(1:19)]),sep="_"),
                   colnames(rates[2+3*c(0:4)]),
                   colnames(fPC[2:57]),
                   "nl","pl","pr","w",
                   Names_var_not)
    NP_fix <- length(Names_fix) #number of fixed parameters
    
    fix_in <- cbind(fBW[fBW$species==species,"BW"],
                    fBW[fBW$species==species,colnames(fBW)%in%substr(Names_fix[grep('fBW',Names_fix)],start=5,stop=nchar(Names_fix[grep('fBW',Names_fix)]))],
                    fCO[fCO$species==species,"CO"],
                    fCO[fCO$species==species,colnames(fCO)%in%substr(Names_fix[grep('fCO',Names_fix)],start=5,stop=nchar(Names_fix[grep('fCO',Names_fix)]))],
                    rates[rates$species==species,colnames(rates)%in%Names_fix],
                    fPC[fPC$species==species,colnames(fPC)%in%Names_fix],
                    data.frame(nl=1,pl=1,pr=1,w=1),
                    TK[TK$species==species&TK$chemical==chemical,colnames(TK)%in%Names_var_not],
                    chem[chem$chemical==chemical,colnames(chem)%in%Names_var_not])
    colnames(fix_in) <- Names_fix
    
    #Change upper limit for clearance depending on physiology of species of interest
      #renal clearance = blood flow to kidneys
      Upper[,"Cl_renal"] <- fix_in$fCO_kidney * fix_in$CO / fix_in$BW
      #hepatic clearance = blood flow to liver (incl portal artery)
      Upper[,"Cl_hepatic"] <- (fix_in$fCO_liver+fix_in$fCO_intestine)*fix_in$CO/fix_in$BW
    
    #create data frames with random samples
    X1 <- matrix(NA, nrow = n_sim, ncol = NP_var)
    colnames(X1) <- Names_var
    X1 <- as.data.frame(X1)
    
    X2 <- matrix(NA, nrow = n_sim, ncol = NP_var)
    colnames(X2) <- colnames(X1)
    X2 <- as.data.frame(X2)
    
    for(i in 1:NP_var){
      X1[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i]) 
      X2[,i] <- runif(n_sim, min = Lower[,i], max = Upper[,i])
    } 
    
    #Sobol design
    sa <- soboljansen(model = NULL, X1, X2, nboot = n_boot, conf = 0.95)
    SARes <- sa$X
    
    #create output data frame
    #SimRes <- matrix(NA, nrow = nrow(SARes), ncol = 3)
    #colnames(SimRes) <- c("Cmax","tmax","AUC_24h")
    #SimRes <- as.data.frame(SimRes)  
      
  } else if (A_type=="VA") {

    # variability analysis with fixed chemical and TK parameters

    Names_var <- c("BW",
                   paste("fBW",colnames(fBW[2+3*c(1:19)]),sep="_"),
                   "CO",
                   paste("fCO",colnames(fCO[2+3*c(1:19)]),sep="_"),
                   colnames(rates[2+3*c(0:4)]))
    
    # initial list with mean + sd values of all varying parameters  

    Means <- cbind(fBW[fBW$species==species,2+3*c(0:19)],
                   fCO[fCO$species==species,2+3*c(0:19)],
                   rates[rates$species==species,2+3*c(0:4)])
    SDs <- cbind(fBW[fBW$species==species,3+3*c(0:19)],
                 fCO[fCO$species==species,3+3*c(0:19)],
                 rates[rates$species==species,3+3*c(0:4)])
    Distributions <- cbind(fBW[fBW$species==species,4+3*c(0:19)],
                          fCO[fCO$species==species,4+3*c(0:19)],
                          rates[rates$species==species,4+3*c(0:4)])
    colnames(Means) = colnames(SDs) = colnames(Distributions) = Names_var
    Names_var_not <- colnames(SDs[which(is.na(SDs))])
    
    Means <- Means[!(colnames(Means)%in%Names_var_not)] #exclude all parameters without SD
    SDs <- SDs[!(colnames(SDs)%in%Names_var_not)] #exclude all parameters without SD
    Distributions <- Distributions[!(colnames(Distributions)%in%Names_var_not)] #exclude all parameters without SD
    
    Names_var <- colnames(SDs)
    NP_var <- length(Names_var) #number of varying parameters
    
    #list with names of all fixed parameters (including 'varying parameters' that are 0 or NA)
    Names_fix <- c(colnames(fPC[2:57]),
                   "nl","pl","pr","w",
                   colnames(TK[3:8]),
                   colnames(chem[2:5]),
                   Names_var_not)
    NP_fix <- length(Names_fix) # number of fixed parameters
    
    fix_in <- cbind(fPC[fPC$species==species,colnames(fPC)%in%Names_fix],
                    data.frame(nl=1,pl=1,pr=1,w=1),
                    TK[TK$species==species&TK$chemical==chemical,colnames(TK)%in%Names_fix],
                    chem[chem$chemical==chemical,colnames(chem)%in%Names_fix],
                    fBW[fBW$species==species,colnames(fBW)%in%substr(Names_fix[grep('fBW',Names_fix)],start=5,stop=nchar(Names_fix[grep('fBW',Names_fix)]))],
                    fCO[fCO$species==species,colnames(fCO)%in%substr(Names_fix[grep('fCO',Names_fix)],start=5,stop=nchar(Names_fix[grep('fCO',Names_fix)]))],
                    rates[rates$species==species,colnames(rates)%in%Names_fix])
    colnames(fix_in) <- Names_fix

    # create data frames with random samples
    X1 <- matrix(NA, nrow = n_sim, ncol = NP_var)
    colnames(X1) <- Names_var
    X1 <- as.data.frame(X1)
    SARes <- X1
    
    for(i in 1:NP_var){
      if (Distributions[i] == "N") {
        SARes[,i] <- rnorm(n_sim, mean = Means[,i], sd = SDs[,i])
      } else if (Distributions[i] == "B") {
        alpha <- ((1 - Means[,i]) / SDs[,i]^2 - 1 / Means[,i]) * Means[,i] ^ 2
        beta <- alpha * (1 / Means[,i] - 1)
        SARes[,i] <- rbeta(n_sim, shape1 = alpha, shape2 = beta)
      } else if (Distributions[i] == "LN") {
        CV <- SDs[,i]/Means[,i]
        mlog <- log(Means[,i]/sqrt(1+CV^2)) #mean of log values
        slog <- sqrt(log(1+CV^2)) #sd of log values
        SARes [,i] <- rlnorm(n_sim, meanlog = mlog, sdlog = slog)
      }
    }

    # create output data frame

    SimRes <- matrix(NA, nrow = nrow(SARes), ncol = n_out*length(t_A))
    
    colnames(SimRes) <- c(paste("blood_t",t_A,sep=""),paste("wholebody_t",t_A,sep=""))
    SimRes <- as.data.frame(SimRes)
    write.csv(SimRes,paste0("./Output/SimRes_blood_",idsim,".csv",sep=""),row.names = FALSE)
  }
 
##############################################################################################################################
# Model run for probabilistic analysis (to be done on separate/external server) 

for (k in 1:(nrow(SARes))) {

      var_in_temp <- SARes
      par_in <-suppressWarnings(cbind(var_in_temp,fix_in))
      
      # adjust clearance rates of hypothetical chemicals based on values drawn for physiology

      if (chem_fix & (chemical=="chem1"|chemical=="chem2"|chemical=="chem3")) {
        par_in$Cl_renal <- par_in$Cl_renal*par_in$fCO_kidney*par_in$CO / par_in$BW
        par_in$Cl_hepatic <- par_in$Cl_hepatic*(par_in$fCO_intestine+par_in$fCO_liver)*par_in$CO / par_in$BW
      } 
      write.csv(par_in,file=paste("./temp/par_in_",k,"_",idsim,".csv",sep=""),row.names = FALSE)
    }
    
    SimRes <- matrix(NA, nrow = nrow(par_in), ncol = n_out*length(t_A))
    colnames(SimRes) <- c(paste("blood_t",t_A,sep=""),paste("wholebody_t",t_A,sep=""))
    if (!chem_fix) {
                    SimRes <- matrix(NA, nrow = nrow(par_in), ncol = 3)
                    colnames(SimRes) <- c("Cmax","tmax","AUC_h")
    }
    
    SimRes <- as.data.frame(SimRes)  
    
##########################################################################################
######################################################################################

    
    for (j in 1:nrow(par_in)) {
     
      print(paste0("Running loop: ", j))
      print(paste("Current time: ", Sys.time()))
      
      out_sim<- multi_tool_fast(
        par_in = par_in[j,], # raw index input data frame
        species = species, regime = regime, route = route, 
        E_dose = E_dose, E_start = E_start, E_end = E_end, E_int = E_int, # exposure scenario
        n_out = n_out, t_start = t_start, t_end = t_end, t_A = t_A, 
        chem_fix = chem_fix,timestep = 3)
      
      SimRes[j,]=out_sim$results 
      res_sims[[j]]= out_sim
    }
    
    write.csv(SimRes,paste("./Output/SimRes_",idsim,".csv",sep=""),row.names = FALSE)

    if (chem_fix) {sim_data=t(SimRes)
    data_blood=data.frame(time=t_A,sim_data[grep("blood",rownames(sim_data)),])
    row.names(data_blood)=NULL
    file.remove(paste("./Output/SimRes_",idsim,".xls",sep=""))
    writeWorksheetToFile(paste("./Output/SimRes_",idsim,".xls",sep=""), data_blood, sheet=paste0("blood_",idsim))
    data_body=data.frame(time=t_A,sim_data[grep("whole",rownames(sim_data)),])
    row.names(data_body)=NULL
    writeWorksheetToFile(paste("./Output/SimRes_",idsim,".xls",sep=""), data_body, sheet=paste0("body_",idsim))
    }
    saveRDS(res_sims,paste("./Output/SimRes_",idsim,".rds",sep=""))
##############################################################################################################################

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
    #add legend
    #legend(x="topright", legend=c("P50","min-max","95% interval"),pch=c(16,1,NA),lty=c(NA,NA,1),seg.len=0.3,lwd=c(1,1,2),col=c("black","black","darkgrey"),bty="n",text.col="black",cex=1, y.intersp=0.8,x.intersp=c(1,1,0.9),merge=F) # Adds a legend box to the plot
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
        FOI[,i]       = sa$S[,1]    # First order indices
        TI[,i]        = sa$T[,1]    # Total indices
        TI.borninf[,i] = sa$T[,4]   # Lower CL total indices 
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


