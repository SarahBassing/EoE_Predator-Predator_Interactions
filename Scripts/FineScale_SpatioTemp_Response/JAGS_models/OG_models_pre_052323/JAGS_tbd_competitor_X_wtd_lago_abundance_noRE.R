  #'  ----------------------------------------
  #'  Latency model - competitor * prey RIA (wtd & lagomorph) effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether there is an interactive effect between
  #'  the species of competitor and the relative abundance of wtd/lagomorph on TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_competitor_X_wtd_lago_abundance_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
      #'  Categorical effect of previous competitor
      beta.competitor[1] <- 0
      for(comp in 2:4) {
        beta.competitor[comp] ~ dnorm(0, 0.01)
      }
      
      #'  Interaction between competitor and prey relative abundance
      beta.interaction.wtd[1] <- 0
      beta.interaction.lago[1] <- 0
      for(comp in 2:4) {
        beta.interaction.wtd[comp] ~ dnorm(0, 0.01)
        beta.interaction.lago[comp] ~ dnorm(0, 0.01)
      }
      
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.competitor[covs[i,2]] + beta.prey[1]*covs[i,8] + beta.prey[2]*covs[i,9] + 
                          beta.interaction.wtd[covs[i,2]]*covs[i,8] + beta.interaction.lago[covs[i,2]]*covs[i,9] 
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor at average prey RAI
      for(comp in 1:4) {
        spp.tbd[comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*0 + 
                             beta.interaction.wtd[comp]*0 + beta.interaction.lago[comp]*0)
      }
      
      #'  Mean TBD per competitor across range of elk relative abundance values
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.wtd[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*newcovs[i,3] + beta.prey[2]*0 + 
                                     beta.interaction.wtd[comp]*newcovs[i,3] + beta.interaction.lago[comp]*0)
        }
      }
      
      #'  Mean TBD per competitor across range of wtd relative abundance values 
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.lago[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,4] + 
                                     beta.interaction.wtd[comp]*0 + beta.interaction.lago[comp]*newcovs[i,4])
        }
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
