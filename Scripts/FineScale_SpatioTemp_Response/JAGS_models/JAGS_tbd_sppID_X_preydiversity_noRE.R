  #'  ----------------------------------------
  #'  Latency model - species ID x prey diversity effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_Time_btwn_Detections_v2.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether there is an interactive effect between 
  #'  the species of competitor and prey diversity on TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_sppID_X_preydiversity_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)
      
      #'  Categorical effect of previous species
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
      
      #'  Interaction between species and prey diversity
      beta.interaction[1] <- 0
      for(spp in 2:nspp) {
        beta.interaction[spp] ~ dnorm(0, 0.01)
      }

      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.div*covs[i,10] + beta.interaction[covs[i,2]]*covs[i,10]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.div*0 + beta.interaction[spp]*0)
      }
      
      #'  Mean TBD per competitor across range of prey diversity values
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.div[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.div*newcovs[i,5] + beta.interaction[spp]*newcovs[i,5])
        }
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
