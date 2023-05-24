  #'  ----------------------------------------
  #'  Latency model - competitor + prey RIA (elk, moose, wtd) effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether species of competitor and/or the
  #'  relative abundance of elk/moose/wtd influences TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_competitor_elk_moose_wtd_abundance_noRE.txt", "
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
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.competitor[covs[i,2]] + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,6] + beta.prey[3]*covs[i,8] 
      }
      
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor at average prey RAI
      for(comp in 1:4) {
        spp.tbd[comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*0)
      }
      
      #'  Mean TBD per competitor across range of elk relative abundance values
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.elk[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0 + beta.prey[3]*0)
        }
      }
      
      #'  Mean TBD per competitor across range of moose relative abundance values 
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.moose[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*newcovs[i,2] + beta.prey[3]*0)
        }
      }
      
      #'  Mean TBD per competitor across range of wtd relative abundance values 
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.wtd[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*newcovs[i,3])
        }
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
