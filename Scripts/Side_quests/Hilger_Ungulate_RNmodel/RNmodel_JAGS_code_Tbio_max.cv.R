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
  
  cat(file = './Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.txt', "
      model{
          
        #'  Define priors
        #'  -------------
        #'  Abundance priors
        beta0 ~ dunif(-10, 10)      # Abundance intercept
        mean.lambda <- exp(beta0)   # Mean lambda for Year 1
          
        #'  Categorical effect for year needs multiple b.year coefficients
        b.year[1] <- 0
        for(yr in 2:nyear) {
          b.year[yr] ~ dnorm(0, 0.001)
        }
          
        #'  Continuous effects for maximum total biomass and CV total biomass
        b.maxTbio ~ dnorm(0, 0.001)
        b.cvTbio ~ dnorm(0, 0.001)
          
        #'  Detection priors
        mean.r ~ dunif(0, 1)        # Detection intercept (on probability scale)
        alpha0 <- logit(mean.r)     # Detection intercept (on logit scale)
          
        #'  Categorical effect for camera setup needs multiple a.setup coefficients
        a.setup[1] <- 0
        for(cam in 2:nsets) {
          a.setup[cam] ~ dnorm(0, 0.001)
        }
          
          
        #'  Define likelihood
        #'  -----------------
        #'  Latent state (abundance)
        for(i in 1:nsites){
          N[i] ~ dpois(lambda[i])
          lambda[i] <- exp(beta0 + b.year[year[i]] + b.maxTbio*max_Tbio[i] + b.cvTbio*cv_Tbio[i])
            
          #'  Log-likelihood of N for WAICj
          log_N[i] <- logdensity.pois(N[i], lambda[i])

          #'  Detection state
          for(j in 1:nsurveys){
            y[i,j] ~ dbern(p[i,j])
            p[i,j] <- 1 - pow((1 - r[i,j]), N[i])
            logit(r[i,j]) <- alpha0 + a.setup[setup[i]]
          
            #'  Log likelihood of y for WAICj
            loglike.waic[i,j] <- logdensity.bin(y[i,j], p[i,j], N[i])
          }
      
          #'  Joint log-likelihood of N and y for WAICj
          loglike.new[i] <- sum(loglike.waic[i,])+log_N[i]
      
        }
          
        #'  Derived parameters
        #'  ------------------
        #'  Mean lambda per year
        for(yr in 1:nyear) {
          lambdaYr[yr] <- exp(beta0 + b.year[yr])
        }
      
        #'  Mean lambda
        mu.lambda <- mean(lambda[])
     
        #'  Mean per-individual detection probability (r) per camera setup
        for(cam in 1:nsets) {
          rSetup[cam] <- 1/(1 + exp(-(alpha0 + a.setup[cam])))
        }
        #'  per-individual detection probability (r) averaged across all camera setups 
        mu.r <- mean(rSetup[])
    
        #'  Mean detection probability (p)
        for(i in 1:nsites) {
          p.occasion[i] <- mean(p[i,])
        }
        mean.p <- mean(p.occasion[])
          
      }
      ")
