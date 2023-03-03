  #'  ------------------------------------
  #'  Prey diversity model
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  ------------------------------------
  #'  Model to test whether predator co-occurrence is non-independent and whether
  #'  basic habitat features influence that relationship.
  #'  ------------------------------------
  
  cat(file = './Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity)_p(setup_effort).txt', "
        model{
        
          ##### Define Priors  ####
          #'  =================
          #'  Priors for parameters of interest 
          #'  Intercepts and slopes for linear models associated with each natural parameter
        
          betaSpp1[1] <- logit(mean.psiSpp1)          # fo occupancy intercepts 
          betaSpp2[1] <- logit(mean.psiSpp2)
          mean.psiSpp1[1] ~ dunif(0, 1)               
          mean.psiSpp2[1] ~ dunif(0, 1)
          
          for(fo_psi in 2:9){                         # fo occupancy slopes (nfirst_order_psi should go through 16 but want to exclude most)
            betaSpp1[fo_psi] ~ dnorm(0, 0.1)
            betaSpp2[fo_psi] ~ dnorm(0, 0.1)
          }
    
          #'  Second order psi priors                   # so occupancy intercepts (nsecond_order_psi)
          for(so_psi in 1:9){
            betaSpp12[so_psi] ~ dnorm(0, 0.1)
          }
      
          #'  First order detection priors (rho)
          alphaSpp1[1] <- logit(mean.pSpp1)           # fo detection intercepts 
          alphaSpp2[1] <- logit(mean.pSpp2)
          mean.pSpp1 ~ dunif(0, 1)                    
          mean.pSpp2 ~ dunif(0, 1)
          
          for(fo_rho in 2:3){          # fo detection slopes (nfirst_order_rho)
            alphaSpp1[fo_rho] ~ dnorm(0, 0.1)  
            alphaSpp2[fo_rho] ~ dnorm(0, 0.1)
          }
        
          #'  Second order detection priors (rho)
          #'  none for now
              
          
          ####  Define Likelihood  ####
          #'  =====================
      
          #'  1. Set up basic hierarchical model
          #'  Latent state and observation processes uses a Categorical distribution 
          #'  since this is essentially a multi-state occupancy model.
      
          #'  Latent state model
          #'  ------------------
          #'  For each site, true occupancy (z) is drawn from a categorical distribution
          #'  with 4 mututally exclusive occupancy probabilities
          
          for(i in 1:nsites) {
            z[i] ~ dcat(lsv[i, (1:ncat)])
          }
        
          #'  Observation model
          #'  -----------------
          #'  For each site and survey occasion, the deteciton data are drawn from a
          #'  categorical distribution with 8 latent states (z)
        
          for(i in 1:nsites) {
            for(j in 1:nsurveys) {
              y[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]])
            }
          }
        
          #'  2. Define arrays containing cell probabilities for categorical distributions
            
          for(i in 1:nsites) {
            #'  Probabilities for each latent state (z), held in latent state vector (lsv)
            lsv[i, 1] <- 1                   # Unoccupied
            lsv[i, 2] <- exp(psiSpp1[i])     # Pr(Spp1 present)
            lsv[i, 3] <- exp(psiSpp2[i])     # Pr(Spp2 present)
            lsv[i, 4] <- exp(psiSpp12[i])    # Pr(Spp1 & Spp2 present)
       
            for(j in 1:nsurveys) {
              #'  Probabilities for each detection array, held in rho detection matrix (rdm) 
              #'  where OS = observed state, TS = true state and each row sums to 1. 
              #'  Exponentiating log odds so rdm holds estimates on probability scale.???
              #'  Reminder - this model assumes NO false positives in the data so
              #'  probability is 0 when OS x TS combinations are not possible.
              #'  Mmmk don't freak out over this section!
              #'  Example 1: when only Spp1 is observed and only Spp1 is truly 
              #'  present, the detection probability is rhoSpp1.
              #'  Example 2: when only Spp1 is observed by in reality Spp1 & Spp2
              #'  are truly present, the detection probability is 
              #'  True state = unoccupied (z = 1 --> 000)
              rdm[i, j, 1, 1] <- 1 # ------------------------------------ OS = unoccupied
              rdm[i, j, 2, 1] <- 0 # ------------------------------------ OS = Spp1 present
              rdm[i, j, 3, 1] <- 0 # ------------------------------------ OS = Spp2 present
              rdm[i, j, 4, 1] <- 0 # ------------------------------------ OS = Spp12 present
              #'  True state = Spp1 present (z = 2 --> 100)
              rdm[i, j, 1, 2] <- 1 # ------------------------------------ OS = unoccupied
              rdm[i, j, 2, 2] <- exp(rhoSpp1[i, j]) # ------------------- OS = Spp1 present
              rdm[i, j, 3, 2] <- 0 # ------------------------------------ OS = Spp2 present
              rdm[i, j, 4, 2] <- 0 # ------------------------------------ OS = Spp12 present
              #'  True state = Spp2 present (z = 3 --> 010 )
              rdm[i, j, 1, 3] <- 1 # ------------------------------------ OS = unoccupied
              rdm[i, j, 2, 3] <- 0 # ------------------------------------ OS = Spp1 present
              rdm[i, j, 3, 3] <- exp(rhoSpp2[i, j]) # ------------------- OS = Spp2 present
              rdm[i, j, 4, 3] <- 0 # ------------------------------------ OS = Spp12 present
              #'  True state = Spp1 & Spp2 present (z = 4 --> 110)
              rdm[i, j, 1, 4] <- 1 # ------------------------------------ OS = unoccupied
              rdm[i, j, 2, 4] <- exp(rhoSpp12[i, j]) # ------------------ OS = Spp1 present
              rdm[i, j, 3, 4] <- exp(rhoSpp21[i, j]) # ------------------ OS = Spp2 present
              rdm[i, j, 4, 4] <- exp(rhoSpp12[i, j] + rhoSpp21[i, j]) # - OS = Spp12 present
            }
            
            #'  3. Define linear models for each fundamental parameter that governs the cell probs
            #'  These are my natural parameters (f1, f2, f3, f12, f13, f23, f123)!
            #'  Linear models for the occupancy parameters on the logit scale
            
            #'  ...for states Spp1, Spp2, Spp3
            #'  Covariate order: Intercept + Setup + Elevation + Forest + Elk + Moose + Mule deer + White-tailed deer + Lagomorph
            psiSpp1[i] <- betaSpp1[1]*psi_cov[i,1] + betaSpp1[2]*psi_cov[i,2] + betaSpp1[3]*psi_cov[i,3] + betaSpp1[4]*psi_cov[i,4] + betaSpp1[5]*psi_cov[i,7] + betaSpp1[6]*psi_cov[i,8] + betaSpp1[7]*psi_cov[i,9] + betaSpp1[8]*psi_cov[i,10] + betaSpp1[9]*psi_cov[i,11] 
            psiSpp2[i] <- betaSpp2[1]*psi_cov[i,1] + betaSpp2[2]*psi_cov[i,2] + betaSpp2[3]*psi_cov[i,3] + betaSpp2[4]*psi_cov[i,4] + betaSpp2[5]*psi_cov[i,7] + betaSpp2[6]*psi_cov[i,8] + betaSpp2[7]*psi_cov[i,9] + betaSpp2[8]*psi_cov[i,10] + betaSpp2[9]*psi_cov[i,11] 
        
            #'  ...for state Spp12
            psiSpp12[i] <- betaSpp12[1]*psi_inxs_cov[i,1] + betaSpp12[2]*psi_inxs_cov[i,2] + betaSpp12[3]*psi_inxs_cov[i,3] + betaSpp12[4]*psi_inxs_cov[i,4] + betaSpp12[5]*psi_inxs_cov[i,7] + betaSpp12[6]*psi_inxs_cov[i,8] + betaSpp12[7]*psi_inxs_cov[i,9] + betaSpp12[8]*psi_inxs_cov[i,10] + betaSpp12[9]*psi_inxs_cov[i,11]   
        
            #'  Linear models for the detection parameters on the logit scale
            for(j in 1:nsurveys) {
              #'  Intercept + Setup + Sampling Effort
              rhoSpp1[i, j] <- alphaSpp1[1]*rho_cov[i,j,1] + alphaSpp1[2]*rho_cov[i,j,2] + alphaSpp1[2]*rho_cov[i,j,3]
              rhoSpp2[i, j] <- alphaSpp2[1]*rho_cov[i,j,1] + alphaSpp2[2]*rho_cov[i,j,2] + alphaSpp2[2]*rho_cov[i,j,3]
        
              #'  Asymetric interactirons between all 3 species
              #'  Currently forcing interactions to equal detection probs above
              rhoSpp12[i, j] <- rhoSpp1[i, j]
              rhoSpp21[i, j] <- rhoSpp2[i, j]
            }
          }
        }
        ", fill=TRUE)
  sink()