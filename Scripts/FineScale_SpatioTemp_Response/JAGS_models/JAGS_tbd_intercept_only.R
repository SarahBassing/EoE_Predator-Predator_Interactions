  #'  ----------------------------------------
  #'  Latency null model - intercept only
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_intercept_only.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0
      }
      
      #'  Derived parameters
      #'  ------------------
      mu.tbd <- exp(alpha0)
      
      } ")
  