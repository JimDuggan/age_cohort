library(deSolve)
library(dplyr)

sei3r <- function(time, stocks, auxs){
  with(as.list(c(stocks, auxs)),{ 
    #convert the stocks vector to a matrix
    epoch<- floor(mod_data$epoch_index[time])
    # cat("Time = ",time," Epoch = ",epoch,"\n")
    states<-matrix(stocks,nrow=mod_data$NUM_COHORTS,ncol=mod_data$NUM_STOCKS)
    colnames(states) <- c("S","E","IS","IP","IC","R","CI")
    rownames(states) <- c("0012","1319","2024","2544","4564","65P")
    
    S   <- states[,"S"]
    E   <- states[,"E"]
    IS  <- states[,"IS"]
    IP  <- states[,"IP"]
    IC  <- states[,"IC"]
    R   <- states[,"R"]
    
    net_contacts <- contacts * K[,,epoch]
    
    
    lambda <- vector(mode="numeric",length = 6)
    names(lambda) <- c("0012","1319","2024","2544","4564","65P")
  
  
    dS_dt   <- -S * lambda                                   # Equation (1) 
    dE_dt   <-  S * lambda - sigma * E                       # Equation (2) 
    dIS_dt  <-  mod_data$SCF * sigma * E       - gamma1 * IS # Equation (3) 
    dIP_dt  <-  (1 - mod_data$SCF) * sigma * E - gamma2 * IP # Equation (4) 
    dIC_dt  <-  gamma2 * IP                    - gamma3 * IC # Equation (5) 
    dR_dt   <- gamma1 * IS + gamma3 * IC                     # Equation (6) 
    dCI_dt  <- rep(0,6)

    dim(net_contacts) <- c(36,1)
    return (list(c(dS_dt,   # Susceptible
                   dE_dt,   # Exposed
                   dIS_dt,
                   dIP_dt,
                   dIC_dt,
                   dR_dt,
                   dCI_dt),
                 CheckSum=sum(states[1:6,1:6]),
                 NetContacts=net_contacts,
                 Tx=mod_data$
                 Lambda=lambda))
  })
}

run_sei3r <- function(){
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
  
  filter(o,time %% 1 == 0)
}

