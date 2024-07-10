#'  ----------------------------------------
#'  Latency model - prey diversity effect, no random effect
#'  ID CRU - Predator Interactions
#'  Sarah Bassing
#'  May 2023
#'  ----------------------------------------
#'  Model sourced to Model_TBD_with_GLMMs.R
#'  
#'  Estimate mean number of minutes elapsed between sequential detections of
#'  different predators and test whether prey diversity influences TBD. 
#'  -----------------------------------------

cat(file = "./Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)

      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.div*covs[i,10]
        
        #'  Goodness-of-fit (Chi-squared test statistic)
        #'  Simulated data from fitted model
        y.sim[i] ~ dexp(tbd_lambda[i])
        #'  GOF (X^2 statistic ---> sum ((o-E)^2)/E )
        fit.obs[i] <- pow((y[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
        fit.sim[i] <- pow((y.sim[i] - tbd_lambda[i]), 2) / tbd_lambda[i]
      }
      
    
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD at average prey diversity value
      mu.tbd <- exp(alpha0 + beta.div*0)
    
      #'  Mean TBD across range of prey diversity values
      for(i in 1:100){
        spp.tbd.div[i] <- exp(alpha0 + beta.div*newcovs[i,5])
      }
      
      #'  For GOF (X^2 statistic)
      chi2.obs <- sum(fit.obs[]) 
      chi2.sim <- sum(fit.sim[]) 
      
      } ")
