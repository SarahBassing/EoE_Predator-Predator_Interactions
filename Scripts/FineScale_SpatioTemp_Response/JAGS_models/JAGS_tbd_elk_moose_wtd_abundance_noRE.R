  #'  ----------------------------------------
  #'  Latency model - prey RIA (elk, moose, wtd) effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether relative abundance of elk/moose/wtd influences TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_elk_moose_wtd_abundance_noRE.txt", "
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
      
        log(tbd_mu[i]) <- alpha0 + beta.prey[1]*covs[i,5] + beta.prey[2]*covs[i,6] + beta.prey[3]*covs[i,8] 
        
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E )
        fit.obs[i] <- pow((y[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        fit.sim[i] <- pow((y.sim[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
      }
      
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD at average prey RAI
      mu.tbd <- exp(alpha0 + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*0)
      
      #'  Mean TBD across range of elk relative abundance values
      for(i in 1:100){
        spp.tbd.elk[i] <- exp(alpha0 + beta.prey[1]*newcovs[i,1] + beta.prey[2]*0 + beta.prey[3]*0)
      }
      
      #'  Mean TBD across range of moose relative abundance values 
      for(i in 1:100){
        spp.tbd.moose[i] <- exp(alpha0 + beta.prey[1]*0 + beta.prey[2]*newcovs[i,2] + beta.prey[3]*0)
      }
      
      #'  Mean TBD across range of wtd relative abundance values 
      for(i in 1:100){
        spp.tbd.wtd[i] <- exp(alpha0 + beta.prey[1]*0 + beta.prey[2]*0 + beta.prey[3]*newcovs[i,3])
      }
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[]) 
      
      } ")
