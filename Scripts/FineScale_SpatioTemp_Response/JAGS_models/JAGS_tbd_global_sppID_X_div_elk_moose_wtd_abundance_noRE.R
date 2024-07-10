  #'  ----------------------------------------
  #'  Latency model - global (prey = elk, moose, wtd), no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_Time_btwn_Detections_v2.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether species of competitor, prey diversity, 
  #'  and the relative abundance of elk/moose/wtd interact to influence TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_global_sppID_X_div_elk_moose_wtd_abundance_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
      #'  Categorical effect of previous species
      beta.sppID[1] <- 0
      for(spp in 2:nspp) {
        beta.sppID[spp] ~ dnorm(0, 0.01)
      }
      
      #'  Interaction between competitor and prey relative abundance
      beta.interaction[1] <- 0
      beta.interaction.elk[1] <- 0
      beta.interaction.moose[1] <- 0
      beta.interaction.wtd[1] <- 0
      for(spp in 2:nspp) {
        beta.interaction[spp] ~ dnorm(0, 0.01)
        beta.interaction.elk[spp] ~ dnorm(0, 0.01)
        beta.interaction.moose[spp] ~ dnorm(0, 0.01)
        beta.interaction.wtd[spp] ~ dnorm(0, 0.01)
      }
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.sppID[covs[i,2]] + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,6] + beta.prey[3]*covs[i,8] + beta.div*covs[i,10] +
                          beta.interaction.elk[covs[i,2]]*covs[i,5] + beta.interaction.moose[covs[i,2]]*covs[i,6] + beta.interaction.wtd[covs[i,2]]*covs[i,8] + beta.interaction[covs[i,2]]*covs[i,10]
      
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E )
        fit.obs[i] <- pow((y[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        fit.sim[i] <- pow((y.sim[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor at average prey RAI and prey diversity
      for(spp in 1:nspp) {
        spp.tbd[spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*0 + beta.div*0 +
                             beta.interaction.elk[spp]*0 + beta.interaction.moose[spp]*0 + beta.interaction.wtd[spp]*0 + beta.interaction[spp]*0)
      }
      
      #'  Mean TBD per competitor across range of elk/moose/wtd relative abundance values & prey diversity
      for(i in 1:100){
        for(spp in 1:nspp){
          spp.tbd.elk[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0 + beta.prey[3]*0 + beta.div*0 +
                                     beta.interaction.elk[spp]*newcovs[i,1] + beta.interaction.moose[spp]*0 + beta.interaction.wtd[spp]*0 + beta.interaction[spp]*0)
          spp.tbd.moose[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,2] + beta.prey[3]*0 + beta.div*0 +
                                     beta.interaction.elk[spp]*0 + beta.interaction.moose[spp]*newcovs[i,2] + beta.interaction.wtd[spp]*0 + beta.interaction[spp]*0)
          spp.tbd.wtd[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*newcovs[i,3] + beta.div*0 +
                                     beta.interaction.elk[spp]*0 + beta.interaction.moose[spp]*0 + beta.interaction.wtd[spp]*newcovs[i,3] + beta.interaction[spp]*0)
          spp.tbd.div[i,spp] <- exp(alpha0 + beta.sppID[spp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*0 + beta.div*newcovs[i,5] +
                                     beta.interaction.elk[spp]*0 + beta.interaction.moose[spp]*0 + beta.interaction.wtd[spp]*0 + beta.interaction[spp]*newcovs[i,5])
        }
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[]) 
      
      } ")
