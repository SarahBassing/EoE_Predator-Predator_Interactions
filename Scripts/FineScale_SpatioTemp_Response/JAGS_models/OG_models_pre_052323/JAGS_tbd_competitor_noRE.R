  #'  ----------------------------------------
  #'  Latency model - competitor effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether species of competitor influences TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_competitor_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      #'  Categorical effect of previous competitor
      beta.competitor[1] <- 0
      for(comp in 2:4) {
        beta.competitor[comp] ~ dnorm(0, 0.01)
      }
      
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.competitor[covs[i,2]] 
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor
      for(comp in 1:4) {
        spp.tbd[comp] <- exp(alpha0 + beta.competitor[comp])
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
