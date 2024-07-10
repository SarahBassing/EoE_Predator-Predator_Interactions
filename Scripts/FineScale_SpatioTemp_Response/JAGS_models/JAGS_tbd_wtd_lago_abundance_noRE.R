  #'  ----------------------------------------
  #'  Latency model - prey RIA (wtd & lagomorphs) effect, no random effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether relative abundance of wtd/lagomorphs influences TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_wtd_lago_abundance_noRE.txt", "
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
      
        log(tbd_mu[i]) <- alpha0 + beta.prey[1]*covs[i,8] + beta.prey[2]*covs[i,9] 
        
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E ) where expected value comes from tbd_lambda
        fit.obs[i] <- pow((y[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        fit.sim[i] <- pow((y.sim[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
      }
      
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD at average prey RAI
      mu.tbd <- exp(alpha0 + beta.prey[1]*0 + beta.prey[2]*0)
      
      #'  Mean TBD across range of wtd/lagomorphs relative abundance values
      for(i in 1:100){
        spp.tbd.wtd[i] <- exp(alpha0 + beta.prey[1]*newcovs[i,3] + beta.prey[2]*0)
        spp.tbd.lago[i] <- exp(alpha0 + beta.prey[1]*0 + beta.prey[2]*newcovs[i,4])
      }
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[]) 
      
      } ")
