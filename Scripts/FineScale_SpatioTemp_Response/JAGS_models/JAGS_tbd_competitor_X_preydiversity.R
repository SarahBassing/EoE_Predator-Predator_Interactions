  #'  ----------------------------------------
  #'  Latency model - competitor x prey diversity effect
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ----------------------------------------
  #'  Model sourced to Model_TBD_with_GLMMs.R
  #'  
  #'  Estimate mean number of minutes elapsed between sequential detections of
  #'  different predators and test whether there is an interactive effect between 
  #'  the species of competitor and prey diversity on TBD. 
  #'  -----------------------------------------
  
  cat(file = "./Outputs/Time_btwn_Detections/tbd_competitor_X_preydiversity.txt", "
      model{
      
      #'  Priors
      #'  ------
      alpha0 ~ dnorm(0, 0.01)
      beta.div ~ dnorm(0, 0.01)
      
      #'  Categorical effect of previous competitor
      beta.competitor[1] <- 0
      for(comp in 2:4) {
        beta.competitor[comp] ~ dnorm(0, 0.01)
      }
      
      #'  Interaction between competitor and prey diversity
      beta.interaction[1] <- 0
      for(comp in 2:4) {
        beta.interaction[comp] ~ dnorm(0, 0.01)
      }
      
      
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
      
        log(tbd_mu[i]) <- alpha0 + beta.competitor[covs[i,2]] + beta.div*covs[i,10] + beta.interaction[covs[i,2]]*covs[i,10] + alpha[site[i]]
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean TBD of each competitor
      for(comp in 1:4) {
        spp.tbd[comp] <- exp(alpha0 + beta.competitor[comp] + beta.div*0 + beta.interaction[comp]*0)
      }
      
      #'  Mean TBD per competitor across range of prey diversity values
      for(i in 1:100){
        for(comp in 1:4){
          spp.tbd.div[i,comp] <- exp(alpha0 + beta.competitor[comp] + beta.div*newcovs[i,5] + beta.interaction[comp]*newcovs[i,5])
        }
      }
      
      #'  Mean TBD
      #'  Note: this overlooks unequal sample sizes contributing to each spp.tbd
      mu.tbd <- mean(spp.tbd[])
      
      } ")
