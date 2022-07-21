#----------------------------------------------------------------------------------------------
#                     Supporting Model for paper:
#  "An Age-Cohort Simulation Model for Generating COVID-19 Scenarios: 
#               A Study from Ireland’s Pandemic Response"
#
#     Age Cohort Transmission Model 
#     Author: Jim Duggan, School of Computer Science, NUI Galway.
# 
#     R Initial Implementation of the Stella Architect Model
#     Version: 0.1
#     July 2022
#----------------------------------------------------------------------------------------------
#  File name: sei3r.R
#
#  Purpose: This is the callback function accessed by deSolve in order to run the simulation.
#           It implements the main equations documented in the paper
#
#           It is a simplified version of the IEMAG operational age cohort model, and
#           focuses on the transmission structure for the Wuhan variant.
#           Parameters are imported from IEMAG Stella sample runs.
#           See the file params.R for the parameter values, which are stored in the mod_data env
#           Some parameters are in the global environment, specifically:
#               - contacts
#               - K matrix (3 D array)
#               - epoch_index (parallel array indicating epoch for each unit of time)
#----------------------------------------------------------------------------------------------
library(deSolve)
library(dplyr)

sei3r <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{ 
    # (1) Find the model epoch
    epoch<- floor(epoch_index[time])
    
    # (2) Take the stocks vector and convert to a matrix 
    states<-matrix(stocks,nrow=mod_data$NUM_COHORTS,ncol=mod_data$NUM_STOCKS)
    colnames(states) <- c("S","E","IS","IP","IC","R","CI")
    rownames(states) <- c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P")
    
    # (3) For ease of explanation, create separate variables for each stock
    S   <- states[,"S"]
    E   <- states[,"E"]
    IS  <- states[,"IS"]
    IP  <- states[,"IP"]
    IC  <- states[,"IC"]
    R   <- states[,"R"]
    
    # (4) Moderate the contacts with the K matrix for the current epoch
    net_contacts <- contacts * K[,,epoch]
    
    # (5) Create a lambda matrix for the force of infections
    lambda_full <- matrix(rep(0,36),nrow = 6)
    
    # (6) Implement the FOI equation the 6x6 matrix
    for(i in seq_along(c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P"))){
      for(j in seq_along(c("Age0012","Age1319","Age2024","Age2544","Age4564","Age65P"))){
        lambda_full[i,j] <- mod_data$Tr[i] * net_contacts[i,j] * (h*IS[j]+IP[j]+IC[j])/population[j]
      }
    }
    
    # (7) Calculate the lambda FOI for each age cohort, summing over each row
    lambda <- apply(lambda_full,1,sum)
    
    # (8) Update all derivatives
    # Equation (1) - Susceptible - dS_dt
    dS_dt   <- -S * lambda 
    # Equation (2) - Exposed - dE_dt
    dE_dt   <-  S * lambda - sigma * E 
    
    dIS_dt  <-  mod_data$SCF * sigma * E       - gamma1 * IS 
    dIP_dt  <-  (1 - mod_data$SCF) * sigma * E - gamma2 * IP  
    dIC_dt  <-  gamma2 * IP                    - gamma3 * IC  
    dR_dt   <- gamma1 * IS + gamma3 * IC                     
    dCI_dt  <- gamma2 * IP * mod_data$rep_fraction

    return (list(c(dS_dt,   
                   dE_dt,   
                   dIS_dt,
                   dIP_dt,
                   dIC_dt,
                   dR_dt,
                   dCI_dt),
                 CheckSum=sum(states[1:6,1:6]),
                 SCF=mod_data$SCF,
                 Tr=mod_data$Tr,
                 Pop=population,
                 K=K[,,epoch],
                 BaseContacts=contacts,
                 NetContacts=net_contacts,
                 Lambda=lambda,
                 STot=sum(S),
                 ETot=sum(E),
                 IPTot=sum(IP),
                 ICTot=sum(IC),
                 ISTot=sum(IS),
                 RTot=sum(R)))
  })
}

run_sei3r <- function(){
  # Set the model parameters
  auxs <- c(sigma   = 1/mod_data$latent_time, 
            gamma1  = 1/mod_data$subclinical_duration,
            gamma2  = 1/mod_data$preclinical_duration,
            gamma3  = 1/mod_data$clinical_duration,
            h       = mod_data$subclincial_infectiousness)
  
  o<-as_tibble(data.frame(ode(y=stocks, 
                    times=mod_data$simtime, 
                    func = sei3r, 
                    parms=auxs, 
                    method="euler")))
  
  # Return the model values for each unit of time
  filter(o,time %% 1 == 0)
}

