  #'  ----------------------------------------
  #'  Latency model - global (prey = elk, moose, wtd), no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether species of competitor, prey diversity, 
  #'  and the relative abundance of elk/moose/wtd interact to influence TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_global_elk_moose_wtd_abundance_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
      #'  Categorical effect of previous competitor
      beta.competitor[1] <- 0
      for(comp in 2:4) {
        beta.competitor[comp] ~ dnorm(0, 0.01)
      }
      
      #'  Interaction between competitor and prey relative abundance
      beta.interaction[1] <- 0
      beta.interaction.elk[1] <- 0
      beta.interaction.moose[1] <- 0
      beta.interaction.wtd[1] <- 0
      for(comp in 2:4) {
        beta.interaction[comp] ~ dnorm(0, 0.01)
        beta.interaction.elk[comp] ~ dnorm(0, 0.01)
        beta.interaction.moose[comp] ~ dnorm(0, 0.01)
        beta.interaction.wtd[comp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.competitor[covs[i,2]] + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,6] + beta.prey[3]*covs[i,8] + beta.div*covs[i,10] +
                          beta.interaction.elk[covs[i,2]]*covs[i,5] + beta.interaction.moose[covs[i,2]]*covs[i,6] + beta.interaction.wtd[covs[i,2]]*covs[i,8] + beta.interaction[covs[i,2]]*covs[i,10]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor at average prey RAI and prey diversity
      for(comp in 1:4) {
        spp.tbd[comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*0 + beta.div*0 +
                             beta.interaction.elk[comp]*0 + beta.interaction.moose[comp]*0 + beta.interaction.wtd[comp]*0 + beta.interaction[comp]*0)
      }
      
      #'  Mean TBD per competitor across range of elk relative abundance values
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.elk[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0 + beta.prey[3]*0 + beta.div*0 +
                                     beta.interaction.elk[comp]*newcovs[i,1] + beta.interaction.moose[comp]*0 + beta.interaction.wtd[comp]*0 + beta.interaction[comp]*0)
        }
      }
      
      #'  Mean TBD per competitor across range of moose relative abundance values 
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.moose[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,2] + beta.prey[3]*0 + beta.div*0 +
                                     beta.interaction.elk[comp]*0 + beta.interaction.moose[comp]*newcovs[i,2] + beta.interaction.wtd[comp]*0 + beta.interaction[comp]*0)
        }
      }
      
      #'  Mean TBD per competitor across range of wtd relative abundance values 
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.wtd[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*newcovs[i,3] + beta.div*0 +
                                     beta.interaction.elk[comp]*0 + beta.interaction.moose[comp]*0 + beta.interaction.wtd[comp]*newcovs[i,3] + beta.interaction[comp]*0)
        }
      }
      
      #'  Mean TBD per competitor across range of prey diversity values 
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.div[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*0 + beta.div*newcovs[i,5] +
                                     beta.interaction.elk[comp]*0 + beta.interaction.moose[comp]*0 + beta.interaction.wtd[comp]*0 + beta.interaction[comp]*newcovs[i,5])
        }
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
