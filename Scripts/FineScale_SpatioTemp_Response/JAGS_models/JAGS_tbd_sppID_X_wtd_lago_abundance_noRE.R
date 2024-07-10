  #'  ----------------------------------------
  #'  Latency model - species ID * prey RIA (elk & wtd) effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_Time_btwn_Detections_v2.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether there is an interactive effect between
  #'  the species of competitor and the relative abundance of elk/wtd on TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_sppID_X_wtd_lago_abundance_noRE.txt", "
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
      
      #'  Interaction between sppID and prey relative abundance
      beta.interaction.wtd[1] <- 0
      beta.interaction.lago[1] <- 0
      for(spp in 2:nspp) {
        beta.interaction.wtd[spp] ~ dnorm(0, 0.01)
        beta.interaction.lago[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.prey[1]*covs[i,8] + beta.prey[2]*covs[i,9] + 
                          beta.interaction.wtd[covs[i,2]]*covs[i,8] + beta.interaction.lago[covs[i,2]]*covs[i,9] 
        
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E )
        fit.obs[i] <- pow((y[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        fit.sim[i] <- pow((y.sim[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each species ID at average prey RAI
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + 
                             beta.interaction.wtd[spp]*0 + beta.interaction.lago[spp]*0)
      }
      
      #'  Mean TBD per species ID across range of wtd/lagomorph relative abundance values
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.wtd[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*newcovs[i,3] + beta.prey[2]*0 + 
                                     beta.interaction.wtd[spp]*newcovs[i,3] + beta.interaction.lago[spp]*0)
          spp.tbd.lago[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,4] + 
                                     beta.interaction.wtd[spp]*0 + beta.interaction.lago[spp]*newcovs[i,4])
        }
      }
      
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[])
      
      } ")
