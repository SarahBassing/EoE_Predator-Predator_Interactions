  #'  ----------------------------------------
  #'  Latency model - sppID + prey RIA (elk, moose, wtd) effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_Time_btwn_Detections_v2.R
  #'  
  #'  Estimate effect of species ID of the most recently detected species on the
  #'  time that elapses before detecting a focal predator species as well as the
  #'  relative abundance of elk/wtd. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_sppID_elk_wtd_abundance_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
      #'  Categorical effect of previous sppID
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,8] 
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each species at average prey RAI
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0)
      }
      
      #'  Mean TBD per competitor across range of elk/wtd relative abundance values
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.elk[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0)
          spp.tbd.wtd[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,3])
        }
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      # mu.tbd <- mean(spp.tbd[])
      
      } ")
