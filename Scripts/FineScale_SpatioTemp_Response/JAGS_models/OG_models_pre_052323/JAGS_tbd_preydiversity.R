#'  ----------------------------------------
#'  Latency model - prey diversity effect
#'  ID CRU - Predator Interactions
#'  Sarah Bassing
#'  May 2023
#'  ----------------------------------------
#'  Model sourced to Model_TBD_with_GLMMs.R
#'  
#'  Estimate mean number of minutes elapsed between sequential detections of
#'  different predators and test whether diversity influences TBD. 
#'  -----------------------------------------

cat(file = "./Outputs/Time_btwn_Detections/tbd_preydiversity.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)
      
      #'  Prior for random effect of camera site
      #'  --------------------------------------
      for(j in 1:ncams) {
        alpha[j] ~ dnorm(0, tau.alpha)
      }
      
      #'  Hyperpriors for random effect
      #'  -----------------------------
      sigma ~ dunif(0, 1)
      tau.alpha <- pow(sigma, -2)
    
      
      #'  Define likelihood
      #'  -----------------
      for(i in 1:ntbd){
        y[i] ~ dexp(tbd_lambda[i])
        
        tbd_lambda[i] <- 1/tbd_mu[i]
      
        log(tbd_mu[i]) <- alpha0 + beta.div*covs[i,10] + alpha[site[i]]
      }
    
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD at average prey diversity value
      mu.tbd <- exp(alpha0 + beta.div*0)
    
      #'  Mean TBD across range of prey diversity values
      for(i in 1:100){
        spp.tbd.div[i] <- exp(alpha0 + beta.div*newcovs[i,5])
      }
      
      } ")
