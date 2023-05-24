  #'  ----------------------------------------
  #'  Latency model - sppID * prey RIA (elk & wtd) effect, no random effect
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
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_sppID_X_elk_wtd_abundance_noRE.txt", "
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
      beta.interaction.elk[1] <- 0
      beta.interaction.wtd[1] <- 0
      for(spp in 2:nspp) {
        beta.interaction.elk[spp] ~ dnorm(0, 0.01)
        beta.interaction.wtd[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,8] + 
                          beta.interaction.elk[covs[i,2]]*covs[i,5] + beta.interaction.wtd[covs[i,2]]*covs[i,8] 
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each beta.sppID at average prey RAI
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + 
                             beta.interaction.elk[spp]*0 + beta.interaction.wtd[spp]*0)
      }
      
      #'  Mean TBD per beta.sppID across range of elk relative abundance values
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.elk[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0 + 
                                     beta.interaction.elk[spp]*newcovs[i,1] + beta.interaction.wtd[spp]*0)
          spp.tbd.wtd[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,3] + 
                                     beta.interaction.elk[spp]*0 + beta.interaction.wtd[spp]*newcovs[i,3])
        }
      }
      
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
