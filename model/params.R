# Create an environment to store the data and params, in order to simplify the environment
mod_data   <- new.env()

# Allocate the epoches to a time vector
epoch_index <- c(rep(1,117),  # Epoch 1, the 2nd wave
                rep(2,53),   # Epoch 2, October lockdown
                          rep(3,28),   # Epoch 3, Christmas opening
                          rep(4,21))   # Epoch 4, January lockdown

mod_data$epoch_01_end <- 118
mod_data$epoch_02_end <- 171
mod_data$epoch_03_end <- 199
mod_data$epoch_04_end <- 219

# Set dimensions for the model, 42 stocks in all
mod_data$NUM_COHORTS <- 6 
mod_data$NUM_STOCKS  <- 7

# Set parameter values - see table 2 in paper
mod_data$latent_time                   <- 3.0
mod_data$subclinical_duration          <- 5.0
mod_data$preclinical_duration          <- 2.1
mod_data$clinical_duration             <- 2.9
mod_data$subclincial_infectiousness    <- 0.5
mod_data$rep_fraction                  <- 0.8

# Transmission parameters
mod_data$Tr        <- c(0.02037858,0.037920761,0.052685521,0.023069804,0.015007166,0.013659096)
names(mod_data$Tr) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")

# Subclinical fraction values
mod_data$scf_0012 <- 1 - 0.44
mod_data$scf_1319 <- 1 - 0.31
mod_data$scf_2024 <- 1 - 0.38
mod_data$scf_2544 <- 1 - 0.44
mod_data$scf_4564 <- 1 - 0.63
mod_data$scf_65P  <- 1 - 0.82
mod_data$SCF <- c(mod_data$scf_0012,
                  mod_data$scf_1319,
                  mod_data$scf_2024,
                  mod_data$scf_2544,
                  mod_data$scf_4564,
                  mod_data$scf_65P)
names(mod_data$SCF) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")

# Set the time vector for the simulation
mod_data$START<-1; 
mod_data$FINISH<-mod_data$epoch_01_end; 
mod_data$STEP<-0.0625;
mod_data$simtime <- seq(mod_data$START, mod_data$FINISH, by=mod_data$STEP)

# Initialise with final states from an earlier calibrated
# IEMAG SEI3R Model of the first wave, up to June 19th 2020

# Initial states for Susceptible (S), equation (1) in paper
mod_data$INIT_S0012 <-	 905290.8497
mod_data$INIT_S1319 <- 	 454873.9279
mod_data$INIT_S2024 <-   298284.5819
mod_data$INIT_S2544 <-  1371881.557
mod_data$INIT_S4564 <-  1214679.199
mod_data$INIT_S65P	<-   692841.885

# Initial states for Exposed (E), equation (2) in paper
mod_data$INIT_E0012	<- 18.33369752
mod_data$INIT_E1319	<- 9.211979781
mod_data$INIT_E2024  <- 6.040776068
mod_data$INIT_E2544  <- 27.78296224
mod_data$INIT_E4564  <- 24.59934398
mod_data$INIT_E65P	  <- 14.0312404

# Initial states for Infected Subclincial (IS), equation (3) in paper
mod_data$INIT_IS0012 <- 	15.03363197
mod_data$INIT_IS1319 <-	7.55382342
mod_data$INIT_IS2024 <-	4.953436376
mod_data$INIT_IS2544 <-	22.78202904
mod_data$INIT_IS4564 <-	20.17146206
mod_data$INIT_IS65P  <-	11.50561713

# Initial states for Infected Preclincial (IP), equation (4) in paper
mod_data$INIT_IP0012 <- 	12.83358826
mod_data$INIT_IP1319 <- 	6.448385847
mod_data$INIT_IP2024 <- 4.228543248
mod_data$INIT_IP2544 <-	19.44807357
mod_data$INIT_IP4564 <- 	17.21954079
mod_data$INIT_IP65P  <- 	9.821868283

# Initial states for Infected Clincial (IC), equation (5) in paper
mod_data$INIT_IC0012 <-	12.83358826
mod_data$INIT_IC1319 <- 	6.448385847
mod_data$INIT_IC2024 <-	4.228543248
mod_data$INIT_IC2544 <-	19.44807357
mod_data$INIT_IC4564 <-	17.21954079
mod_data$INIT_IC65P  <-	9.821868283

# Initial states for Removed (R), equation (6) in paper
mod_data$INIT_R0012 <-	5830.115811
mod_data$INIT_R1319 <-	2929.40957
mod_data$INIT_R2024 <-	1920.96679
mod_data$INIT_R2544 <-	8834.981994
mod_data$INIT_R4564 <-	7822.591386
mod_data$INIT_R65P	 <- 4461.934449

# Create contact matrix, equation (8) in paper. Data calculate via call to socialmixr.
mod_data$POLYMOD_NLAge0012   <- c(10.63768116,	0.823710577,	0.156364331,	3.998045646,	1.372655264,	0.379476843)
mod_data$POLYMOD_NLAge1319   <- c(1.639350163,	10.03846154,	0.723630018,	1.924509662,	1.75181074,	  0.437465485)
mod_data$POLYMOD_NLAge2024   <- c(0.474564246,	1.103511373,	2.416666667,	5.494865795,	3.112879895,	0.469634658)
mod_data$POLYMOD_NLAge2544   <- c(2.638270135,	0.63810849,	  1.19473415,	  7.333333333,	3.609683633,	0.815480045)
mod_data$POLYMOD_NLAge4564   <- c(1.023029168,	0.656019328,	0.764419181,	4.076844658,	4.896551724,	1.230996538)
mod_data$POLYMOD_NLAge65P    <- c(0.495837393,	0.287210759,	0.202188668,	1.614714783,	2.158163243,	2.47826087)
contacts    <- rbind(mod_data$POLYMOD_NLAge0012,
                     mod_data$POLYMOD_NLAge1319,
                     mod_data$POLYMOD_NLAge2024,
                     mod_data$POLYMOD_NLAge2544,
                     mod_data$POLYMOD_NLAge4564,
                     mod_data$POLYMOD_NLAge65P)
colnames(contacts) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")
rownames(contacts) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")

stocks <- c(S0012   = mod_data$INIT_S0012,
            S1319   = mod_data$INIT_S1319,
            S2024   = mod_data$INIT_S2024,
            S2544   = mod_data$INIT_S2544,
            S4564   = mod_data$INIT_S4564,
            S65P    = mod_data$INIT_S65P,
            E0012   = mod_data$INIT_E0012,
            E1319   = mod_data$INIT_E1319,
            E2024   = mod_data$INIT_E2024,
            E2544   = mod_data$INIT_E2544,
            E4564   = mod_data$INIT_E4564,
            E65P    = mod_data$INIT_E65P,
            IS0012  = mod_data$INIT_IS0012,
            IS1319  = mod_data$INIT_IS1319,
            IS2024  = mod_data$INIT_IS2024,
            IS2544  = mod_data$INIT_IS2544,
            IS4564  = mod_data$INIT_IS4564,
            IS65P   = mod_data$INIT_IS65P,
            IP0012  = mod_data$INIT_IP0012,
            IP1319  = mod_data$INIT_IP1319,
            IP2024  = mod_data$INIT_IP2024,
            IP2544  = mod_data$INIT_IP2544,
            IP4564  = mod_data$INIT_IP4564,
            IP65P   = mod_data$INIT_IP65P,
            IC0012  = mod_data$INIT_IC0012,
            IC1319  = mod_data$INIT_IC1319,
            IC2024  = mod_data$INIT_IC2024,
            IC2544  = mod_data$INIT_IC2544,
            IC4564  = mod_data$INIT_IC4564,
            IC65P   = mod_data$INIT_IC65P,
            R0012   = mod_data$INIT_R0012,
            R1319   = mod_data$INIT_R1319,
            R2024   = mod_data$INIT_R2024,
            R2544   = mod_data$INIT_R2544,
            R4564   = mod_data$INIT_R4564,
            R65P    = mod_data$INIT_R65P,
            CI0012  = 0,
            CI1319  = 0,
            CI2024  = 0,
            CI2544  = 0,
            CI4564  = 0,
            CI65P   = 0)

# Create population vector, equation (9) in paper
population        <- c(911180,	457833, 300225, 1380806, 1222581, 697349)
names(population) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")

# Create K matrix, equation (10) in paper
K <- array(dim = c(6,6,4))
rownames(K) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")
colnames(K) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")

# set first epoch values
K[,,1] <- 1

# set second epoch values, based on Stella Powell calibration
K[,,2] <- 0 