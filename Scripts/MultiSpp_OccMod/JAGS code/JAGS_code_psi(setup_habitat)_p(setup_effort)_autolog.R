  #'  ------------------------------------
  #'  Autologistic habitat model, no species interactions
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2023
  #'  ------------------------------------
  #'  Model to test whether predator occurrence is influenced by basic habitat 
  #'  features. Assumes species occur and are detected independently of one another.
  #'  ------------------------------------
  
  cat(file = './Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat)_p(setup_effort)_autolog.txt', "
      model{
          
        #### Define Priors  ####
        #'  ================
        #'  Priors for parameters of interest 
        #'  Intercepts and slopes for linear models associated with each natural parameter
            
        #'  First order occupancy paramters
        for(fo_psi in 1:4){                         
          betaSpp1[fo_psi] ~ dnorm(0, 0.1)
          betaSpp2[fo_psi] ~ dnorm(0, 0.1)
        }
      
        #'  Second order occupancy parameters
        #'  Fix second-order interaction to 0
        betaSpp12 <- 0
        
        #'  First order detection paramters and phi for year 2 autologisitc model   
        for(k in 1:nspec) {
          alpha0[k] ~ dnorm(0, 0.1)
          alpha1[k] ~ dnorm(0, 0.1)
          alpha2[k] ~ dnorm(0, 0.1)
          phi[k] ~ dnorm(0, 0.1)
        }
      
        #'  Second order detection priors (rho)
        #'  Assumes no second-order interactions by setting these to 0
        alphaSpp12 <- 0
        alphaSpp21 <- 0
          
            
        ####  Define Likelihood  ####
        #'  =====================
      
        #'  Latent state and observation processes uses a Categorical distribution 
        #'  since this is essentially a multi-state occupancy model.
            
        for(i in 1:nsites) {
    
          #'  Year 1 latent state model
          #'  -------------------------
          #'  For each site, true occupancy (z) is drawn from a categorical distribution
          #'  with 4 mututally exclusive occupancy probabilities
          
          #'  Natural parameters
          #'  First order covariate order: Intercept[1] + Setup[2] + Elevation[3] + Forest[4]
          f1[i,1] <- betaSpp1[1]*psi_cov[i,1] + betaSpp1[2]*psi_cov[i,2] + betaSpp1[3]*psi_cov[i,3] + betaSpp1[4]*psi_cov[i,4]
          f2[i,1] <- betaSpp2[1]*psi_cov[i,1] + betaSpp2[2]*psi_cov[i,2] + betaSpp2[3]*psi_cov[i,3] + betaSpp2[4]*psi_cov[i,4]
    
          #'  Second order
          #'  No interaction right now
          f12[i,1] <- 0
    
          #'  Psi: probability of latent occupancy state for each site
          #'  1) both present, 2) species 1 present, 3) species 2 present, 4) both absent
          denom[i,1] <- 1 + exp(f1[i,1]) + exp(f2[i,1]) + exp(f1[i,1] + f2[i,1] + f12[i,1])
          psi[1,i,1] <- exp(f1[i,1] + f2[i,1] + f12[i,1]) / denom[i,1]
          psi[2,i,1] <- exp(f1[i,1]) / denom[i,1]
          psi[3,i,1] <- exp(f2[i,1]) / denom[i,1]
          psi[4,i,1] <- 1 / denom[i,1]
      
          #'  Model latent occupancy state (z) as a categorical random variable
          #'  with 4 mutually exclusive occupancy probabilities (psi)
          z[i,1] ~ dcat(psi[,j,1])
    
          #'  Year 1 observation model
          #'  ------------------------
          #'  Loop over each species (k) and sampling occasion (j)
          for(k in 1:nspec) {
            for(j in 1:nsurveys) {
              #'  Detection as a function of camera deployment and sampling effort
              logit(p[k,i,j,1]) <- alpha0[k]*rho_cov[i,j,1,1] + alpha1[k]*rho_cov[i,j,1,3] + alpha3[k]*rho_cov[i,j,1,5]
              #'  Model observations (y) as true occupancy state (z) * detection probability (p)
              #'  Note: Xcat is a matrix of 2^n rows and n columns with 0s & 1s
              y[k,i,j,1] ~ dbern(Xcat[z[i,1],k] * p[k,i,j,1])
            }
          }
      
          #'  Year 2 latent state model
          #'  -------------------------
          #'  For years >1, first order natural parameters have additional parameter 
          #'  with coefficient phi that is multiplied by each site's occupancy
          #'  state from previous year
          for(t in 2:nyear) {
            f1[i,t] <- betaSpp1[1]*psi_cov[i,1] + betaSpp1[2]*psi_cov[i,2] + betaSpp1[3]*psi_cov[i,3] + betaSpp1[4]*psi_cov[i,4] + phi[1]*Xcat[z[i,t-1],1]
            f2[i,t] <- betaSpp2[1]*psi_cov[i,1] + betaSpp2[2]*psi_cov[i,2] + betaSpp2[3]*psi_cov[i,3] + betaSpp2[4]*psi_cov[i,4] + phi[2]*Xcat[z[i,t-1],2]
            f12[i,t] <- 0
    
          #'  Psi
          denom[i,t] <- 1 + exp(f1[i,t]) + exp(f2[i,t]) + exp(f1[i,t] + f2[i,t] + f12[i,t])
          psi[1,i,t] <- exp(f1[i,t] + f2[i,t] + f12[i,t]) / denom[i,t]
          psi[2,i,t] <- exp(f1[i,t]) / denom[i,t]
          psi[3,i,t] <- exp(f2[i,t]) / denom[i,t]
          psi[4,i,t] <- 1 / denom[i,t]
          
          #'  Model latent occupancy based on categories of psi
          z[i,t] ~ dcat(psi[,i,t])
    
          #'  Year 2 observation model
          #'  ------------------------
          for(k in 1:nspec) {
            for(j in 1:nsurveys) {
              logit(p[k,i,j,t]) <- alpha0[k]*rho_cov[i,j,t,1] + alpha1[k]*rho_cov[i,j,t,3] + alpha3[k]*rho_cov[i,j,t,5]
              y[k,i,j,t] ~ dbern(Xcat[z[i,t],k] * p[k,i,j,t])
            }
          }
        }
      }

      ", fill=TRUE)
