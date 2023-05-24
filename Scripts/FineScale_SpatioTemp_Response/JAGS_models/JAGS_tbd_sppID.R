  #'  ----------------------------------------
  #'  Latency model - Previous species detected effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_Time_btwn_Detections_v2.R
  #'  
  #'  Estimate effect of species ID of the most recently detected species on the
  #'  time that elapses before detecting a focal predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_sppID.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      #'  Categorical effect of previous species
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each species - focal predator
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp])
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
