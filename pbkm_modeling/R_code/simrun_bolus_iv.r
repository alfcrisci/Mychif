regime <- "bolus"                   # exposure regime (bolus/continuous)
route <- "iv"                     # exposure route (oral/iv)
                                    # only relevant when regime="bolus" or iv
# Simulation parameters

A_type <- "VA"                # type of probabilistic analysis ("SA" or "VA")
chem_fix <- TRUE              # fixing the chemical and TK parameters or not (TRUE/FALSE)
n_boot <- 1000                # number of boostrap iterations (for SA)
n_out  <- 2                   # number of compartments to output (blood, total body)
t_A <- c(seq(0.05,0.2,by=0.05),seq(0.25,t_end,by=0.25))     # all time points for analysis (h), only relevant when chem_fix=TRUE


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
    
      write.csv(par_in,file=paste("./temp/par_in_",k,"_",idsim,".csv",sep=""),row.names = FALSE)
    }
    
    SimRes <- matrix(NA, nrow = nrow(par_in), ncol = n_out*length(t_A))
    colnames(SimRes) <- c(paste("blood_t",t_A,sep=""),paste("wholebody_t",t_A,sep=""))
    if (!chem_fix) {
                    SimRes <- matrix(NA, nrow = nrow(par_in), ncol = 3)
                    colnames(SimRes) <- c("Cmax","tmax","AUC_h")
    }
    
    SimRes <- as.data.frame(SimRes)  
    
    # par_in
    
    for (j in 1:nrow(par_in)) {
     
      print(paste0("Running loop: ", j))
      print(paste("Current time: ", Sys.time()))
      
      out_sim<- multi_tool_fast(
        par_in = par_in[j,], # raw index input data frame
        species = species, regime = regime, route = route, 
        E_dose = E_dose, E_start = E_start, E_end = E_end, E_int = E_int, # exposure scenario
        n_out = n_out, t_start = t_start, t_end = t_end, t_A = t_A, 
        chem_fix = chem_fix)
      
      SimRes[j,]=out_sim
    }

    data_sim=t(rbind(SimRes,t_A))
    data_blood=data_sim[1:length(t_A),]
    ini=length(t_A)+1
    data_body=data_sim[ini:nrow(data_sim),]
    
    row.names(data_blood)=NULL
    row.names(data_body)=NULL
    
    res_sims_blood=list()
    for ( i in 1:n_sim) {
      t=as.numeric(data_blood[,n_sim+1])
      C=as.numeric(data_blood[,i])
      t=t[which(!is.na(C))]
      C=C[which(!is.na(C))]
      res_sims_blood[[i]]=cbind(t,C)
    }
    res_sims_body=list()
    for ( i in 1:n_sim) {
      t=as.numeric(data_body[,n_sim+1])
      C=as.numeric(data_body[,i])
      t=t[which(!is.na(C))]
      C=C[which(!is.na(C))]
      res_sims_body[[i]]=cbind(t,C)
    }
    