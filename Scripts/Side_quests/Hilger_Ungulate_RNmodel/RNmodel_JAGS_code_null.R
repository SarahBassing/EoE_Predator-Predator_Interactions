  #'  -------------------------------
  #'  Royle-Nichols abundance model
  #'  ID CRU - Predator Interactions
  #'  Sarah B. Bassing
  #'  August 2024
  #'  -------------------------------
  #'  RN model to estimate relative abundance from binary detection/non-detection,
  #'  assuming heterogeneous abundance affects detection probability.
  #'  
  #'  Relevant parameters and data:
  #'  lambda: The number of animals available for detection at site i, N[i], is a 
  #'  Poisson-distributed random variable with mean lambda[i]. N is latent.
  #'  
  #'  y: The observed detection/non-detection data for each site i and survey 
  #'  occasion j, y[i,j], is a Bernoulli-distributed random variable with  
  #'  probability p[i,j]. y[i,j] is conditional on N[i].
  #'  
  #'  p & r: The probability of detecting occupancy at site i during occasion j, 
  #'  p[i,j], is a function of the per-individual detection probability, r[i,j], 
  #'  a binomial sampling probability that a particular individual is detected at 
  #'  site i during occasion j, and the number of individuals at site i, N[i].
  #'  -------------------------------
  
  cat(file = './Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt', "
      model{
          
        #'  Define priors
        #'  -------------
        #'  Abundance priors
        lambda ~ dgamma(0.001, 0.001)  
        # beta0 ~ dunif(-10, 10)      # Abundance intercept
        # mean.lambda <- exp(beta0)   # Mean lambda 
          
        #'  Detection priors
        mean.r ~ dunif(0, 1)        # Detection intercept (on probability scale)
        # alpha0 <- logit(mean.r)     # Detection intercept (on logit scale)
          
        #' #'  Categorical effect for camera setup needs multiple a.setup coefficients
        #' a.setup[1] <- 0
        #' for(cam in 2:nsets) {
        #'   a.setup[cam] ~ dnorm(0, 0.001)
        #' }
          

        #'  Define likelihood
        #'  -----------------
        #'  Latent state (abundance)
        for(i in 1:nsites){
          N[i] ~ dpois(lambda)
          # lambda <- mean.lambda 
          
          #'  Log-likelihood of N for WAICj
          log_N[i] <- logdensity.pois(N[i], lambda)
            
          #'  Detection state
          for(j in 1:nsurveys){
            y[i,j] ~ dbern(p[i,j])
            p[i,j] <- 1 - pow((1 - mean.r), N[i])
            # p[i,j] <- 1 - pow((1 - mean.r[i,j]), N[i])
            # logit(r[i,j]) <- alpha0 + a.setup[setup[i]]
            
            #'  Log likelihood of y for WAICj
            loglike.waic[i,j] <- logdensity.bin(y[i,j], p[i,j], N[i])
          }
      
          #'  Joint log-likelihood of N and y for WAICj
          loglike.new[i] <- sum(loglike.waic[i,])+log_N[i]
      
        }
        
        #' #'  Goodness-of-Fit test (code adapted from Mike Meredith: https://github.com/mikemeredith/AHM_code/blob/main/AHM1_ch06/AHM1_06.08.R)
        #' #' Posterior predictive distributions of chi2 discrepancy
        #' for (i in 1:M) {
        #'   for (j in 1:J) {
        #'     C.sim[i,j] ~ dbin(p, N[i]) # Create new data set under model
        #'     e.count[i,j] <- N[i] * p   # Expected datum
        #'     # Chi-square discrepancy for the actual data
        #'     chi2.actual[i,j] <- pow((C[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
        #'     # Chi-square discrepancy for the simulated ('perfect') data
        #'     chi2.sim[i,j] <- pow((C.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
        #'     # Add small value e to denominator to avoid division by zero
        #'   }
        #' }
        #' # Add up individual chi2 values for overall fit statistic
        #' fit.actual <- sum(chi2.actual[,])  # Fit statistic for actual data set
        #' fit.sim <- sum(chi2.sim[,])        # Fit statistic for a fitting model
        #' c.hat <- fit.actual / fit.sim      # c-hat estimate
        #' bpv <- step(fit.sim-fit.actual)    # Bayesian p-value
          
        #'  Derived parameters
        #'  ------------------
        #'  Mean lambda
        mu.lambda <- lambda
        #' lambda <- exp(beta0)
  
        #' #'  Mean per-individual detection probability (r) per camera setup
        #' for(cam in 1:nsets) {
        #'   rSetup[cam] <- 1/(1 + exp(-(alpha0 + a.setup[cam])))
        #' }
        #' #'  per-individual detection probability (r) averaged across all camera setups 
        #' mu.r <- mean(rSetup[])
        mu.r <- mean.r
    
        #'  Mean detection probability (p)
        for(i in 1:nsites) {
          p.occasion[i] <- mean(p[i,])
        }
        mean.p <- mean(p.occasion[])
          
      }
      ")
