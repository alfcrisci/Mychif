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

if (chem_fix & any(i/timescale==t_A)) {
  
  wholebody_temp<- (sum(phys$M[phys$comp!="feces"&
                                 phys$comp!="urine"&
                                 phys$comp!="lumen"&
                                 phys$comp!="air"&
                                 phys$comp!="metab"&
                                 phys$comp!="feed"&
                                 phys$comp!="milk"&
                                 phys$comp!="egg"])*par_in$MW) / BW
  
  
}
