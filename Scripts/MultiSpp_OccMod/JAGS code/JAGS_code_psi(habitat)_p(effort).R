  #'  ------------------------------------
  #'  Habitat model
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  ------------------------------------
  #'  Model to test whether predator co-occurrence is non-independent and whether
  #'  basic habitat features influence that relationship.
  #'  ------------------------------------
  
  cat(file = './Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(habitat)_p(effort).txt', "
      model{
      
      ####  Define priors  ####
      #'  -----------------
      #'  Intercepts & slopes associated with each natural parameter
      #'  First order psi intercepts
      mean.psi1 ~ dunif(0, 1)
      mean.psi2 ~ dunif(0, 1)
      betaSpp1 <- logit(mean.psi1)
      betaSpp2 <- logit(mean.psi2)
      
      #'  First order psi slopes
      #'  Keep in mind this will generate priors for all nfirst_order_psi covariates 
      #'  even if don't need them all in the model
      for(fo_psi in 2:nfirst_order_psi) {
        betaSpp1[fo_psi] ~ dnorm(0, 0.001)
        betaSpp2[fo_psi] ~ dnorm(0, 0.001)
      }
      
      #'  Second order psi priors
      for(so_psi in 1:nsecond_order_psi) {
        betaSpp12[so_psi] ~ dnorm(0, 0.001)
      }
      
      #'  First order rho intercepts
      mean.p1 ~ dunif(0, 1)
      mean.p2 ~ dunif(0, 1)
      alphaSpp1 <- logit(mean.p1)
      alphaSpp2 <- logit(mean.p2)
      
      #'  First order rho slopes
      for(fo_rho in 2:nfirst_order_rho) {
        alphaSpp1[fo_rho] ~ dnorm(0, 0.001)
        alphaSpp2[fo_rho] ~ dnorm(0, 0.001)
      }
      
      #'  Second order rho priors
      for(so_rho in 1:nsecond_order_psi) {
        alpha12[so_rho] ~ dnorm(0, 0.001)
      }

      
      ####  Likelihood  ####
      #'  --------------
      #'  1. Basic heirarchical model
      #'  Latent state model
      for(i in 1:nsites) {
        z[i] ~ dcat(lsv[i, (1:ncat)])
      }
      
      #'  Detection model
      for(i in 1:nsites) {
        for(j in 1:nsurveys) {
          y[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]])
        }
      }
      
      #'  2. Define latent state vector & observation matrices
      for(i in 1:nsites) {
        #'  Latent stete probabilities in latent state vector (lsv)
        #'  Probabilities for each state (ncat)
        lsv[i,1] <- 1  # -------------------------------------------- 00
        lsv[i,2] <- exp(psi1[i])  # --------------------------------- 10
        lsv[i,3] <- exp(psi2[i])  # --------------------------------- 01
        lsv[i,4] <- exp(psi12[i]) # --------------------------------- 11
      
        for(j in 1:nsurveys) {
          #'  Detection matrix (i, j, OS = observation state, TS = true state)
          #'  rdm = rho detection matrix; each row sums to 1
          #'  OS along rows, TS, along columns
          #'  True state = 00
          rdm[i, j, 1, 1] <- 1  # ----------------------------------- OS = none
          rdm[i, j, 2, 1] <- 0  # ----------------------------------- OS = Spp1
          rdm[i, j, 3, 1] <- 0  # ----------------------------------- OS = Spp2
          rdm[i, j, 4, 1] <- 0  # ----------------------------------- OS = both
      
          #'  True state = 10
          rdm[i, j, 1, 2] <- 1  # ----------------------------------- OS = none
          rdm[i, j, 2, 2] <- exp(rho1[i, j])  # --------------------- OS = Spp1
          rdm[i, j, 3, 2] <- 0  # ----------------------------------- OS = Spp2
          rdm[i, j, 4, 2] <- 0  # ----------------------------------- OS = both
      
          #'  True state = 01
          rdm[i, j, 1, 3] <- 1  # ----------------------------------- OS = none
          rdm[i, j, 2, 3] <- 0  # ----------------------------------- OS = Spp1
          rdm[i, j, 3, 3] <- exp(rho2[i, j])  # --------------------- OS = Spp2
          rdm[i, j, 4, 3] <- 0  # ----------------------------------- OS = both
      
          #'  True state = 11
          rdm[i, j, 1, 4] <- 1  # ----------------------------------- OS = none
          rdm[i, j, 2, 4] <- exp(rho12[i, j])  # -------------------- OS = Spp1
          rdm[i, j, 3, 4] <- exp(rho21[i, j])  # -------------------- OS = Spp2
          rdm[i, j, 4, 4] <- exp(rho12[i, j] + rho21[i, j])  # ------ OS = both
        }
      
        #'  3. Define linear models for each fundamental parameter that governs the cell probs (lsv & rdm)
        #'  These are my natural parameters (f1, f2, f12)!
        #'  Linear models for the occupancy parameters on the logit scale
        #'  Covariate order: Intercept + Setup + Elevation + Forest
        psi1[i] <- betaSpp1[1]*psi_covs[i,1] + betaSpp1[2]*psi_covs[i,2] + betaSpp1[3]*psi_covs[i,3] + betaSpp1[4]*psi_covs[i,4]
        psi2[i] <- betaSpp2[1]*psi_covs[i,1] + betaSpp2[2]*psi_covs[i,2] + betaSpp2[3]*psi_covs[i,3] + betaSpp2[4]*psi_covs[i,4]
          
        #'  Linear models for species co-occurrence
        psi12[i] <- betaSpp12[1]*psi_inxs_covs[i,1] + betaSpp12[2]*psi_inxs_covs[i,2] + betaSpp12[3]*psi_inxs_covs[i,3] + betaSpp12[4]*psi_inxs_covs[i,4]
      
        #'  Linear models for detection parameters
        for(j in 1:nsurveys) {
          #'  Basesline detection linear predictors: Intercept + Setup + Sampling Effort
          rho1[i, j] <- alphaSpp1[1]*rho_covs[i, j, 1] + alphaSpp1[2]*rho_covs[i, j, 3] + alphaSpp1[3]*rho_covs[i, j, 5]
          rho2[i, j] <- alphaSpp2[1]*rho_covs[i, j, 1] + alphaSpp2[2]*rho_covs[i, j, 3] + alphaSpp2[3]*rho_covs[i, j, 5] 
          #'  Asymetric interaciton between both species (currently just using baseline probabilities but could use rho_inxs_covs)
          rho12[i, j] <- rho1[i, j]
          rho21[i, j] <- rho2[i, j]
        }
      }
      
      ####  Derived parameters  ####
      #'  ----------------------
      #'  TBD
      }
      ")