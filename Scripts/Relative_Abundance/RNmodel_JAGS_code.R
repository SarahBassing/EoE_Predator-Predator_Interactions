  #'  -------------------------------
  #'  Royle-Nichols abundance model
  #'  ID CRU - Predator Interactions
  #'  Sarah B. Bassing
  #'  November 2023
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
  
  cat(file = './Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt', "
      model{
      
      #'  Define priors
      #'  -------------
      #'  Abundance priors
      beta0 ~ dnorm(0,0.001)
      beta1 ~ dnorm(0,0.001)
      beta2 ~ dnorm(0,0.001)
      beta3 ~ dnorm(0,0.001)
      
      #'  Categorical effect for GMU needs multiple beta4 coefficients
      beta4[1] <- 0
      for(gmu in 2:ngmu) {
        beta4[gmu] ~ dnorm(0, 0.001)
      }
      
      #'  Detection priors
      alpha0 ~ dnorm(0,0.001)
      alpha1 ~ dnorm(0,0.001)
      
      #'  Categorical effect for camera setup needs multiple alpha2 coefficients
      alpha2[1] <- 0
      for(cam in 2:2) {
        alpha2[cam] ~ dnorm(0,0.001)
      }
      
      #'  Define likelihood
      #'  -----------------
      #'  Latent state (abundance)
      for(i in 1:nsite) {
        log.lambda[i] <- beta0 + beta1*PercFor[i] + beta2*Elev[i] + beta3*pow(Elev[i],2) + beta4[GMU[i]]
        log(lambda[i]) <- log.lambda[i] 
        # log(lambda[i]) <- beta0 + beta1*PercFor[i] + beta2*Elev[i] + beta3*pow(Elev[i],2) + beta4[GMU[i]]
        N[i] ~ dpois(lambda[i])
      
      #'  Detection state
        for(j in 1:nsurvey) {
          logit.r[i,j] <- alpha0 + alpha1*nDays[i] + alpha2[Setup[i]]
          logit(r[i,j]) <- logit.r[i,j]
          # logit(r[i,j]) <- alpha0 + alpha1*nDays[i] + alpha2[Setup[i]]
          p[i,j] <- 1 - (1 - r[i,j]) ^ N[i]
          # p[i,j] <- 1 - pow((1 - r[i,j]), N[i])
          y[i,j] ~ dbern(p[i,j])
        }
      }
      
      #' #'  Derived paramters
      #' #'  -----------------
      #' #'  Site-specific abundance
      #' for(i in 1:nsite) {
      #'   siteN[i] <- exp(beta0 + beta1*PercFor[i] + beta2*Elev[i] + beta3*pow(Elev[i],2) + beta4[GMU[i]])
      #' }
      #' 
      #' #'  Mean abundance per GMU
      #' for(gmu in 1:ngmu) {
      #'   gmuN[gmu] <- exp(beta0 + beta4[gmu])
      #' }
      #' 
      #' #'  Mean per-individual detection probability per camera setup
      #' for(cam in 1:2) {
      #'   rSetup[cam] <- 1/(1 + exp(alpha0 + alpha2[cam]))
      #' }
      #' 
      #' #'  Total abundance across camera sites (NOTE this is different than summing across gmuN)
      #' totalN <- sum(N[i:nsite])
      #' 
      #' #'  Total sites and total sites occupied (N > 0)
      #' nSites <- count(N[i:nsite])
      #' occSites <- count(N[i:nsite] > 0)
      #' 
      #' #'  Mean occupancy (psi)
      #' meanpsi <- occSites/nSites
      #' 
      #' #'  Mean detection probability, per-individual detection probability, and lambda
      #' meanp <- mean(p[])
      #' meanr <- mean(r[])
      #' meanlambda <- mean(lambda[])
      
      }
      ")
  