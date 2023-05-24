  #'  ----------------------------------------
  #'  Latency model - prey RIA (elk & wtd) effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_Time_btwn_Detections_v2.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether relative abundance of elk/wtd influences TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      
      for(pp in 1:npp) {
        beta.prey[pp] ~ dnorm(0, 0.01)
      }
      
    
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,8]
      }
      
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD at average prey RAI
      mu.tbd <- exp(alpha0 + beta.prey[1]*0 + beta.prey[2]*0)
      
      #'  Mean TBD across range of elk relative abundance values
      for(i in 1:100){
        spp.tbd.elk[i] <- exp(alpha0 + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0)
      }
      
      #'  Mean TBD across range of wtd relative abundance values 
      for(i in 1:100){
        spp.tbd.wtd[i] <- exp(alpha0 + beta.prey[1]*0 + beta.prey[2]*newcovs[i,3])
      }
      
      } ")
