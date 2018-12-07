#==================================================================
# One-compartment model single dose for veterinary species
# Oral exposure
#=================================================================

# Parameters

B = 0.95       # Bioavailability
D = 5.5        # Dose in mg/kg
V = 0.31       # Volume in L/kg
kel = 0.178    # Elimination rate in hr-1
ka = 0.18      # Absorption rate in hr-1

t <- seq(0, 24, by = 0.5)     #time interval in hours

# Equation one-comaprtment model oral exposure

C <- ((B * D * ka) / (V * (ka-kel))) * ((exp(-kel * t)) - (exp(-ka * t)))

# Creating dataframe 

conc <- data.frame(t, C)

#Plot

plot(t, C, xlab = "Time (h)", ylab = "Concentration (ug/ml)")

#==================================================================
# One-compartment model single dose for veterinary species
# Intravenous
#=================================================================

# Parameters

D = 80           # Dose in mg/kg
V = 0.706        # Volume in L/kg
kel = 0.949      # Elimination rate in hr-1

t <- seq(0, 24, by = 0.5)     #time interval in hours

#Equation one-comaprtment model

C <- (D / V) * (exp(-kel * t))

#Creating dataframe 

conc1 <- data.frame(t, C)

#Plot

plot(t, C, xlab = "Time (h)", ylab = "Concentration (ug/ml)")


