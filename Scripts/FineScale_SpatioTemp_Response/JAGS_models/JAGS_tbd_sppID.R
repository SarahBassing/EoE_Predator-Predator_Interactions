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
        
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  Expected observation from fitted model
        y.hat[i] <- tbd_lambda[i]
        #'  Expected observation from simulated data
        y.sim.hat[i] <- tbd_lambda[i]  
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E )
        fit.obs[i] <- pow((y[i]-y.hat[i]),2) / y.hat[i] #(sqrt(y.hat[i])) 
        fit.sim[i] <- pow((y.sim[i]-y.sim.hat[i]),2) / y.sim.hat[i] #(sqrt(y.sim.hat[i])) 
        #https://stackoverflow.com/questions/48024836/posterior-predictive-check-in-jags-dimension-mismatch-error has sqrt(y.hat) in denominator
        #https://www.flutterbys.com.au/stats/tut/tut11.2b.html has y.hat in denominator (no sqrt)
        #Anderson-Darling test as an alternative but not sure how to code this: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm
        
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
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[]) 
      
      } ")
