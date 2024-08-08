  #'  ------------------------------------
  #'  Habitat model, no species interactions
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2023
  #'  ------------------------------------
  #'  Model to test whether predator occurrence is influenced by basic habitat 
  #'  features. Assumes species occur and are detected independently of one another.
  #'  ------------------------------------
  
  cat(file = './Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_GoF.txt', "
      model{
          
        #### Define Priors  ####
        #'  ================
        #'  Priors for parameters of interest 
        #'  Intercepts and slopes for linear models associated with each natural parameter
        
        #'  First order occupancy intercerpts (psi) 
        betaSpp1[1] <- logit(mean.psiSpp1)           
        betaSpp2[1] <- logit(mean.psiSpp2)
        mean.psiSpp1 ~ dunif(0, 1)               
        mean.psiSpp2 ~ dunif(0, 1)
            
        #'  First order occupancy slopes (psi)
        for(fo_psi in 2:6){                         
          betaSpp1[fo_psi] ~ dnorm(0, 0.1)
          betaSpp2[fo_psi] ~ dnorm(0, 0.1)
        }
      
        #'  Second order occupancy intercerpt (psi) 
        #'  Fix second-order interaction to 0
        betaSpp12 <- 0
        
        #'  First order detection intercepts (rho)
        alphaSpp1[1] <- logit(mean.pSpp1)           
        alphaSpp2[1] <- logit(mean.pSpp2)
        mean.pSpp1 ~ dunif(0, 1)                    
        mean.pSpp2 ~ dunif(0, 1)
         
        #'  First order detection slopes (rho)   
        for(fo_rho in 2:3){                         
          alphaSpp1[fo_rho] ~ dnorm(0, 0.1)  
          alphaSpp2[fo_rho] ~ dnorm(0, 0.1)
        }
      
        #'  Second order detection priors (rho)
        #'  Assumes no second-order interactions by setting these to 0
        alphaSpp12 <- 0
        alphaSpp21 <- 0
                
            
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
          # z.sim[i] ~ dcat(lsv[i,(1:ncat)]) # simulated data for GoF, latent state
        }
          
        #'  Observation model
        #'  -----------------
        #'  For each site and survey occasion, the deteciton data are drawn from a
        #'  categorical distribution with 4 latent states (z)
        
        #Indx <- c(1,2,3,4)
          
        for(i in 1:nsites) {
          for(j in 1:nsurveys) {
            y[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]])
            
            #'  Draw a replicate data set under fitted model
            y.sim[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]]) #(1:ncat)
            
            #' #'  Derived parameters for GoF check
            #' y.hat[i,j] <- y[i,j]
            #' y.sim.hat[i,j] <- max(rdm[i, j, , z[i]])
            
            # r.obs[i,j] <- (y[i,j]-y.hat[i,j])/sqrt(y.hat[i,j]*(1-y.hat[i,j]))
            # r.sim[i,j] <- (y.sim[i,j]-y.sim.hat[i,j])/sqrt(y.sim.hat[i,j]*(1-y.sim.hat[i,j]))
            
            #' #'  Grab expected value based on model
            #' # y.hat[i,j] <- Indx[rdm[i, j, , z[i]] == max(rdm[i, j, , z[i]])]
            #' # y.hat[i,j] <- Indx[max(rdm[i, j, (1:ncat), z[i]]) == rdm[i, j, (1:ncat), z[i]]]
            #' y.hat[i,j] <- max(rdm[i, j, (1:ncat), z[i]]) # returns only 1s
            #' # r.obs[i,j] <- (y[i,j]-y.hat[i,j])/sqrt(y.hat[i,j]*(1-y.hat[i,j]))
            #' 
            #' #'  Simulate replicate data set under fitted model for GoF
            #' y.sim[i,j] ~ dcat(rdm[i, j, (1:ncat), z.sim[i]])
            #' # y.sim.hat[i,j] <- Indx[max(rdm[i, j, , z.sim[i]])] 
            #' # # y.sim.hat[i,j] <- Indx[max(rdm[i, j, (1:ncat), z.sim[i]]) == rdm[i, j, (1:ncat), z.sim[i]]] 
            #' # #y.sim.hat[i,j] <- max(rdm[i, j, (1:ncat), z.sim[i]]) # returns only 1s
            #' # r.sim[i,j] <- (y.sim[i,j]-y.sim.hat[i,j])/sqrt(y.sim.hat[i,j]*(1-y.sim.hat[i,j]))
            
            #derived parameters for Goodness-of-Fit check
            y2[i,j] <- y[i,j]
            yrep2[i,j] <- y.sim[i,j]
            
            # seperate species
            y_A[i,j] <- ifelse(y2[i,j]==2 || y2[i,j]==4, 1, 0)
            y_B[i,j] <- ifelse(y2[i,j]==3 || y2[i,j]==4, 1, 0)

            yrep_A[i,j] <- ifelse(yrep2[i,j]==2 || yrep2[i,j]==4, 1, 0)
            yrep_B[i,j] <- ifelse(yrep2[i,j]==3 || yrep2[i,j]==4, 1, 0)
            
          }
        }
        
        #'  Calcualte observed, replicate, expected detection frequencies
        for(i in 1:nsites) {
          
          #'  Det. frequencies for observed and replicated data
          detfreq_A[i] <- sum(y_A[i,])
          detfreqrep_A[i] <- sum(yrep_A[i,])

          detfreq_B[i] <- sum(y_B[i,])
          detfreqrep_B[i] <- sum(yrep_B[i,])

          #'  Separate z by species
          z_A[i] <- ifelse(z[i]==2 || z[i]==4, 1, 0)
          z_B[i] <- ifelse(z[i]==3 || z[i]==4, 1, 0)

          #'  Expected detection frequencies under the model
          for (j in 1:nsurveys){
            tmp_A[i,j] <- z_A[i] * rdm[i, j, 2, 2]
            tmp_B[i,j] <- z_B[i] * rdm[i, j, 3, 3]
          } 
      
          E_A[i] <- sum(tmp_A[i,])     # Expected number of detections for A
          E_B[i] <- sum(tmp_B[i,])     # Expected number of detections for B
      
          #'  Chi-square and Freeman-Tukey discrepancy measures
          #'  ..... for actual data 
          x2_A[i] <- pow((detfreq_A[i] - E_A[i]), 2) / (E_A[i] + 0.0001)
          x2_B[i] <- pow((detfreq_B[i] - E_B[i]), 2) / (E_B[i] + 0.0001)
        
          ft_A[i] <- pow((sqrt(detfreq_A[i]) - sqrt(E_A[i])), 2) 
          ft_B[i] <- pow((sqrt(detfreq_B[i]) - sqrt(E_B[i])), 2)
        
          #'  ..... for replicated data set
          x2rep_A[i] <- pow((detfreqrep_A[i] - E_A[i]), 2) / (E_A[i] + 0.0001)
          x2rep_B[i] <- pow((detfreqrep_B[i] - E_B[i]), 2) / (E_B[i] + 0.0001)
        
          ftrep_A[i] <- pow((sqrt(detfreqrep_A[i]) - sqrt(E_A[i])), 2)
          ftrep_B[i] <- pow((sqrt(detfreqrep_B[i]) - sqrt(E_B[i])), 2)
        } 
        
        #'  Add up overall test statistic and compute fit stat ratio
        chi2.obs_A <- sum(x2_A[])
        chi2.obs_B <- sum(x2_B[])
        
        chi2.sim_A <- sum(x2rep_A[])
        chi2.sim_B <- sum(x2rep_B[])
        
        chi2ratio_A <- chi2.obs_A/chi2.sim_A
        chi2ratio_B <- chi2.obs_B/chi2.sim_B
        
        ft.obs_A <- sum(ft_A[])
        ft.obs_B <- sum(ft_B[])
        
        ft.sim_A <- sum(ftrep_A[])
        ft.sim_B <- sum(ftrep_B[])
        
        ftratio_A <- ft.obs_A/ft.sim_B
        ftratio_B <- ft.obs_A/ft.sim_B
        
        #' #'  GOF Chi2 test statistic
        #' chi2.obs <- sum(r.obs[,]^2)
        #' chi2.sim <- sum(r.sim[,]^2)
          
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
            #'  Exponentiating log odds so rdm holds estimates on probability scale.
            #'  Reminder - this model assumes NO false positives in the data so
            #'  probability is 0 when OS x TS combinations are not possible.
            #'  Mmmk don't freak out over this section!
            #'  Example 1: when only Spp1 is observed and only Spp1 is truly 
            #'  present, the detection probability is rhoSpp1.
            #'  Example 2: when only Spp1 is observed by in reality Spp1 & Spp2
            #'  are truly present, the detection probability is rhoSpp12
            #'  True state = unoccupied (z = 1 --> 00)
            rdm[i, j, 1, 1] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 1] <- 0 # ------------------------------------ OS = Spp1 present
            rdm[i, j, 3, 1] <- 0 # ------------------------------------ OS = Spp2 present
            rdm[i, j, 4, 1] <- 0 # ------------------------------------ OS = Spp12 present
            #'  True state = Spp1 present (z = 2 --> 10)
            rdm[i, j, 1, 2] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 2] <- exp(rhoSpp1[i, j]) # ------------------- OS = Spp1 present
            rdm[i, j, 3, 2] <- 0 # ------------------------------------ OS = Spp2 present
            rdm[i, j, 4, 2] <- 0 # ------------------------------------ OS = Spp12 present
            #'  True state = Spp2 present (z = 3 --> 01)
            rdm[i, j, 1, 3] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 3] <- 0 # ------------------------------------ OS = Spp1 present
            rdm[i, j, 3, 3] <- exp(rhoSpp2[i, j]) # ------------------- OS = Spp2 present
            rdm[i, j, 4, 3] <- 0 # ------------------------------------ OS = Spp12 present
            #'  True state = Spp1 & Spp2 present (z = 4 --> 11)
            rdm[i, j, 1, 4] <- 1 # ------------------------------------ OS = unoccupied
            rdm[i, j, 2, 4] <- exp(rhoSpp1[i, j])  # ------------------ OS = Spp1 present
            rdm[i, j, 3, 4] <- exp(rhoSpp2[i, j])  # ------------------ OS = Spp2 present
            rdm[i, j, 4, 4] <- exp(rhoSpp12[i, j] + rhoSpp21[i, j]) # - OS = Spp12 present
            #' #'  True state = unoccupied (z = 1 --> 00)
            #' rdm[i, j, 1, 1] <- 1 # ------------------------------------ OS = unoccupied
            #' rdm[i, j, 2, 1] <- 0 # ------------------------------------ OS = Spp1 present
            #' rdm[i, j, 3, 1] <- 0 # ------------------------------------ OS = Spp2 present
            #' rdm[i, j, 4, 1] <- 0 # ------------------------------------ OS = Spp12 present
            #' #'  True state = Spp1 present (z = 2 --> 10)
            #' rdm[i, j, 1, 2] <- 1 - exp(rhoSpp1[i, j]) # --------------- OS = unoccupied
            #' rdm[i, j, 2, 2] <- exp(rhoSpp1[i, j]) # ------------------- OS = Spp1 present
            #' rdm[i, j, 3, 2] <- 0 # ------------------------------------ OS = Spp2 present
            #' rdm[i, j, 4, 2] <- 0 # ------------------------------------ OS = Spp12 present
            #' #'  True state = Spp2 present (z = 3 --> 01)
            #' rdm[i, j, 1, 3] <- 1 - exp(rhoSpp2[i, j]) # --------------- OS = unoccupied
            #' rdm[i, j, 2, 3] <- 0 # ------------------------------------ OS = Spp1 present
            #' rdm[i, j, 3, 3] <- exp(rhoSpp2[i, j]) # ------------------- OS = Spp2 present
            #' rdm[i, j, 4, 3] <- 0 # ------------------------------------ OS = Spp12 present
            #' #'  True state = Spp1 & Spp2 present (z = 4 --> 11)
            #' rdm[i, j, 1, 4] <- 1 - exp(rhoSpp12[i, j] + rhoSpp21[i, j]) # - OS = unoccupied
            #' rdm[i, j, 2, 4] <- exp(rhoSpp12[i, j]) # ------------------ OS = Spp1 present
            #' rdm[i, j, 3, 4] <- exp(rhoSpp21[i, j]) # ------------------ OS = Spp2 present
            #' rdm[i, j, 4, 4] <- exp(rhoSpp12[i, j] + rhoSpp21[i, j]) # - OS = Spp12 present
          }
              
          #'  3. Define linear models for each fundamental parameter that governs the cell probs
          #'  These are my natural parameters (f1, f2, f12)!
          #'  Linear models for the occupancy parameters on the logit scale
              
          #'  ...for states Spp1, Spp2
          #'  Covariate order: Intercept[1] + Setup[2] + Year[3] + Forest[4] + Elevation[5] + TRI[17]# + HumanDet[14] + Footprint[16] 
          psiSpp1[i] <- betaSpp1[1]*psi_cov[i,1] + betaSpp1[2]*psi_cov[i,2] + betaSpp1[3]*psi_cov[i,3] + betaSpp1[4]*psi_cov[i,4] + betaSpp1[5]*psi_cov[i,5] + betaSpp1[6]*psi_cov[i,17] # + betaSpp1[7]*psi_cov[i,14] + betaSpp1[6]*psi_cov[i,16]
          psiSpp2[i] <- betaSpp2[1]*psi_cov[i,1] + betaSpp2[2]*psi_cov[i,2] + betaSpp2[3]*psi_cov[i,3] + betaSpp2[4]*psi_cov[i,4] + betaSpp2[5]*psi_cov[i,5] + betaSpp2[6]*psi_cov[i,17] # + betaSpp2[7]*psi_cov[i,14] + betaSpp2[6]*psi_cov[i,16]
          
          #'  ...for state Spp12
          #'  Don't forget - second order parameter set to 0 so no interaction
          psiSpp12[i] <- psiSpp1[i] + psiSpp2[i] + betaSpp12*psi_inxs_cov[i,1]
          
          #'  Baseline linear predictors for detection
          #'  Covariate order: Intercept[1] + Setup[3] + Sampling Effort[5]
          for(j in 1:nsurveys) {
            rhoSpp1[i, j] <- alphaSpp1[1]*rho_cov[i,j,1] + alphaSpp1[2]*rho_cov[i,j,3] + alphaSpp1[3]*rho_cov[i,j,5] 
            rhoSpp2[i, j] <- alphaSpp2[1]*rho_cov[i,j,1] + alphaSpp2[2]*rho_cov[i,j,3] + alphaSpp2[3]*rho_cov[i,j,5] 
          
            #'  Asymetric interactions between both species
            #'  Don't forget - second order parameters set to 0 so no interactions
            rhoSpp12[i, j] <- rhoSpp1[i, j] + alphaSpp12*rho_inxs_cov[i,j,1] 
            rhoSpp21[i, j] <- rhoSpp2[i, j] + alphaSpp21*rho_inxs_cov[i,j,1] 
          }
        }
      }
      ", fill=TRUE)
