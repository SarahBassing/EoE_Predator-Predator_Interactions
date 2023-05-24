  #'  ----------------------------------------
  #'  Latency model - Competitor effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_Time_btwn_Detections_v1.R
  #'  
  #'  Estimate effect of species ID of the most recently detected species on the
  #'  time that elapses before detecting a focal predator species.
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_competitor.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      #'  Categorical effect of previous species
      beta.competitor[1] <- 0
      for(comp in 2:4) {
        beta.competitor[comp] ~ dnorm(0, 0.01)
      }
      
      
      #'  Prior for random effect of camera site
      #'  --------------------------------------
      for(j in 1:ncams) {
        alpha[j] ~ dnorm(0, tau.alpha)
      }
      
      #'  Hyperpriors for random effect
      #'  -----------------------------
      sigma ~ dunif(0, 1)
      tau.alpha <- pow(sigma, -2)
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.competitor[covs[i,2]] + alpha[site[i]]
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
