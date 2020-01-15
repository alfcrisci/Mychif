##########################################################################################################################################################à
# Qt V PC M Qin Qout

############################################################################################################################

setwd("")
chem <- read.csv("data/chem_data_prob.csv", stringsAsFactors = FALSE)             # physicochemical parameters constant are in 1/minutes
TK <- read.csv("data/tk_data/TK_prob.csv",stringsAsFactors=FALSE)                 # toxicokinetic parameters
fBW <- read.csv("data/physio/fBW_animals_prob.csv",stringsAsFactors = FALSE)      # organ fractions, incl var-distribution parameters
fCO <- read.csv("data/physio/fCO_animals_prob.csv",stringsAsFactors = FALSE)      # blood flow fractions, incl var-distribution parameters
rates <- read.csv("data/physio/rates_animals_prob.csv",stringsAsFactors = FALSE)  # physiological rates (not chemical-dependent)
fPC <- read.csv("data/physio/PC_tissue_prob.csv",stringsAsFactors = FALSE)        # tissue composition, fixed during var-analysis

############################################################################################################################

BW <- BW            # bodyweight (kg)
CO <- CO/20         # cardiac output (fro minutes to L/3sec)

# fBW (peso frazione)        fCO (blood flow fractions)           PC           V (volume)          Q

Q_abomasum=Q_adipose=Q_brain=Q_carcass=Q_crop=Q_gizzard=Q_heart=Q_intestine=Q_kidney=0;
Q_liver=Q_lun=Q_mamgland=Q_muscle=Q_omasum=Q_reprod=Q_reticulum=Q_rumen=Q_stomach=Q_lumen=0;
Q_ven=Q_art=0;
Q_milk=Q_egg=Q_urine=Q_air=Q_feces=Q_metab=0




dQ_met <- ifelse(!is.na(phys$Vmax)&!is.na(phys$Km),(phys$Vmax*(phys$C/phys$PC))/(phys$Km+phys$C/phys$PC),0)
dQ_met_liver=phys$C[phys$comp=="liver"]*phys$Cl[phys$comp=="liver"]


phys$dM_in <- phys$dM_in + phys$dM_met[phys$comp=="liver"]*phys$fbile*phys$fbact
phys$dM_met[phys$comp=="liver"]*(1-phys$fbile[phys$comp=="lumen"]*phys$fbact[phys$comp=="lumen"])   

#ELIMINATION

phys$dM_exh <- ifelse(phys$comp=="ven", CO*phys$C[phys$comp=="ven"],0) #all mass in venous blood is transported

phys$dM_out <- phys$dM_out + phys$dM_exh

phys$dM_in[phys$comp=="art"] <- phys$dM_in[phys$comp=="art"] + CO*phys$C[phys$comp=="ven"] * (CO/(CO+Qexhale*phys$PC[phys$comp=="air"]))
phys$dM_in[phys$comp=="air"] <- phys$dM_in[phys$comp=="air"] + CO*phys$C[phys$comp=="ven"] * ((Qexhale*phys$PC[phys$comp=="air"])/(CO+Qexhale*phys$PC[phys$comp=="air"]))


dQ_urine <- phys$C[phys$comp=="kidney"]*phys$Cl[phys$comp=="kidney"]
dQ_milk <- phys$C[phys$comp=="mamgland"]*phys$Qmilk[phys$comp=="mamgland"]
dQ_reprod <- phys$C[phys$comp=="reprod"]*phys$Qegg[phys$comp=="reprod"]
dQ_feces <- phys$M*phys$kgastric

dQ_exh=dQ_urine+dQ_milk+dQ_reprod+dQ_feces
