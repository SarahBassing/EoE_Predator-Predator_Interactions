  #'  ---------------------------------
  #'  Null model 1: No interactions
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  ---------------------------------
  #'  Included year effects to account for temporal correlation in occurrence 
  #'  related to re-sampling same site each year but otherwise unrelated to 
  #'  environmental/ biotic variables of interest.
  #'  ---------------------------------
  
  cat(file = './Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_altGoF.txt', "
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
        for(fo_psi in 2:2){                         
          betaSpp1[fo_psi] ~ dnorm(0, 0.1)
          betaSpp2[fo_psi] ~ dnorm(0, 0.1)
        }
      
        #'  Second order occupancy intercerpt (psi)
        #'  Fix second-order interaction to 0
        betaSpp12 <- 0
            
        #'  First order detection priors (rho)
        alphaSpp1[1] <- logit(mean.pSpp1)           
        alphaSpp2[1] <- logit(mean.pSpp2)
        mean.pSpp1 ~ dunif(0, 1)                    
        mean.pSpp2 ~ dunif(0, 1)
      
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
        }
          
        #'  Observation model
        #'  -----------------
        #'  For each site and survey occasion, the deteciton data are drawn from a
        #'  categorical distribution with 8 latent states (z)
          
        for(i in 1:nsites) {
          for(j in 1:nsurveys) {
            y[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]])
            
            #' #'  Draw a replicate data set under fitted model
            #' y.sim[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]]) 
            #' 
            #' #'  Derived parameters for Goodness-of-Fit check
            #' y2[i,j] <- y[i,j]
            #' yrep2[i,j] <- y.sim[i,j]
            #' 
            #' #' Seperate out by species
            #' y_A[i,j] <- ifelse(y2[i,j]==2 || y2[i,j]==4, 1, 0)
            #' y_B[i,j] <- ifelse(y2[i,j]==3 || y2[i,j]==4, 1, 0)
            #' 
            #' yrep_A[i,j] <- ifelse(yrep2[i,j]==2 || yrep2[i,j]==4, 1, 0)
            #' yrep_B[i,j] <- ifelse(yrep2[i,j]==3 || yrep2[i,j]==4, 1, 0)
          }
        }
        
        #' #'  Calcualte observed, replicate, expected detection frequencies
        #' for(i in 1:nsites) {
        #'   
        #'   #'  Det. frequencies for observed and replicated data
        #'   detfreq_A[i] <- sum(y_A[i,])
        #'   detfreqrep_A[i] <- sum(yrep_A[i,])
        #' 
        #'   detfreq_B[i] <- sum(y_B[i,])
        #'   detfreqrep_B[i] <- sum(yrep_B[i,])
        #' 
        #'   #'  Separate z by species
        #'   z_A[i] <- ifelse(z[i]==2 || z[i]==4, 1, 0)
        #'   z_B[i] <- ifelse(z[i]==3 || z[i]==4, 1, 0)
        #' 
        #'   #'  Expected detection frequencies under the model
        #'   for (j in 1:nsurveys){
        #'     tmp_A[i,j] <- z_A[i] * rdm[i, j, 2, 2]
        #'     tmp_B[i,j] <- z_B[i] * rdm[i, j, 3, 3]
        #'   } 
        #' 
        #'   E_A[i] <- sum(tmp_A[i,])     # Expected number of detections for A
        #'   E_B[i] <- sum(tmp_B[i,])     # Expected number of detections for B
        #' 
        #'   #'  Chi-square and Freeman-Tukey discrepancy measures
        #'   #'  ..... for actual data 
        #'   x2_A[i] <- pow((detfreq_A[i] - E_A[i]), 2) / (E_A[i] + 0.0001)
        #'   x2_B[i] <- pow((detfreq_B[i] - E_B[i]), 2) / (E_B[i] + 0.0001)
        #' 
        #'   ft_A[i] <- pow((sqrt(detfreq_A[i]) - sqrt(E_A[i])), 2) 
        #'   ft_B[i] <- pow((sqrt(detfreq_B[i]) - sqrt(E_B[i])), 2)
        #' 
        #'   #'  ..... for replicated data set
        #'   x2rep_A[i] <- pow((detfreqrep_A[i] - E_A[i]), 2) / (E_A[i] + 0.0001)
        #'   x2rep_B[i] <- pow((detfreqrep_B[i] - E_B[i]), 2) / (E_B[i] + 0.0001)
        #' 
        #'   ftrep_A[i] <- pow((sqrt(detfreqrep_A[i]) - sqrt(E_A[i])), 2)
        #'   ftrep_B[i] <- pow((sqrt(detfreqrep_B[i]) - sqrt(E_B[i])), 2)
        #' } 
        #' 
        #' #'  Add up overall test statistic and compute fit stat ratio
        #' chi2.obs_A <- sum(x2_A[])
        #' chi2.obs_B <- sum(x2_B[])
        #' 
        #' chi2.sim_A <- sum(x2rep_A[])
        #' chi2.sim_B <- sum(x2rep_B[])
        #' 
        #' chi2ratio_A <- chi2.obs_A/chi2.sim_A
        #' chi2ratio_B <- chi2.obs_B/chi2.sim_B
        #' 
        #' ft.obs_A <- sum(ft_A[])
        #' ft.obs_B <- sum(ft_B[])
        #' 
        #' ft.sim_A <- sum(ftrep_A[])
        #' ft.sim_B <- sum(ftrep_B[])
        #' 
        #' ftratio_A <- ft.obs_A/ft.sim_A
        #' ftratio_B <- ft.obs_B/ft.sim_B
          
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
            #'  are truly present, the detection probability is rhoSpp12.
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
            rdm[i, j, 2, 4] <- exp(rhoSpp12[i, j]) # ------------------ OS = Spp1 present
            rdm[i, j, 3, 4] <- exp(rhoSpp21[i, j]) # ------------------ OS = Spp2 present
            rdm[i, j, 4, 4] <- exp(rhoSpp12[i, j] + rhoSpp21[i, j]) # - OS = Spp12 present
          }
              
          #'  3. Define linear models for each fundamental parameter that governs the cell probs
          #'  These are my natural parameters (f1, f2, f12)!
          #'  Linear models for the occupancy parameters on the logit scale
              
          #'  ...for states Spp1, Spp2
          #'  Covariate order: Intercept[1] + Year[5]
          psiSpp1[i] <- betaSpp1[1]*psi_cov[i,1] + betaSpp1[2]*psi_cov[i,5] 
          psiSpp2[i] <- betaSpp2[1]*psi_cov[i,1] + betaSpp2[2]*psi_cov[i,5] 
          
          #'  ...for state Spp12
          #'  Don't forget - second order parameter set to 0 so no interaction
          psiSpp12[i] <- psiSpp1[i] + psiSpp2[i] + betaSpp12*psi_inxs_cov[i,1]
          
          #'  Linear models for the detection parameters on the logit scale
          for(j in 1:nsurveys) {
            #'  Intercept 
            rhoSpp1[i, j] <- alphaSpp1[1]*rho_cov[i,j,1] 
            rhoSpp2[i, j] <- alphaSpp2[1]*rho_cov[i,j,1] 
          
            #'  Asymetric interactions between both species
            #'  Fixing to be same as species-sepcific detection probability
            rhoSpp12[i, j] <- rhoSpp1[i, j] 
            rhoSpp21[i, j] <- rhoSpp2[i, j] 
          }
        }
        
        #'  For Goodness-of-Fit test
        for(i in 1:nsites) {
          for(j in 1:nsurveys) {
            #'  Draw a replicate data set under fitted model
            y.sim[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]]) 
            
            #'  Expected detection probability based each underlying occurrence state
            #'  Need to normalized rdms because dcat can take values that sum to >1
            #'  rdm is a vector of unnormalized non-negative probability *weights* of length N
            #'  dcat pmf --> P(Y = y) = pi[y]/sum[pi[i, 1:N]] (pg. 50 JAGS manual)
            #'  True state = unoccupied (z = 1 --> 00)
            p[i, j, 1, 1] <- rdm[i, j, 1, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(unoccupied)
            p[i, j, 2, 1] <- rdm[i, j, 2, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(Spp1 detected)
            p[i, j, 3, 1] <- rdm[i, j, 3, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(Spp2 detected)
            p[i, j, 4, 1] <- rdm[i, j, 4, 1]/sum(rdm[i, j, 1:4, 1]) # -- Pr(Spp12 detected)
            #'  True state = Spp1 present (z = 2 --> 10)
            p[i, j, 1, 2] <- rdm[i, j, 1, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(unoccupied)
            p[i, j, 2, 2] <- rdm[i, j, 2, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(Spp1 detected)
            p[i, j, 3, 2] <- rdm[i, j, 3, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(Spp2 detected)
            p[i, j, 4, 2] <- rdm[i, j, 4, 2]/sum(rdm[i, j, 1:4, 2]) # -- Pr(Spp12 detected)
            #'  True state = Spp2 present (z = 3 --> 01)
            p[i, j, 1, 3] <- rdm[i, j, 1, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(unoccupied)
            p[i, j, 2, 3] <- rdm[i, j, 2, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(Spp1 detected)
            p[i, j, 3, 3] <- rdm[i, j, 3, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(Spp2 detected)
            p[i, j, 4, 3] <- rdm[i, j, 4, 3]/sum(rdm[i, j, 1:4, 3]) # -- Pr(Spp12 detected)
            #'  True state = Spp1 & Spp2 present (z = 4 --> 11)
            p[i, j, 1, 4] <- rdm[i, j, 1, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(unoccupied)
            p[i, j, 2, 4] <- rdm[i, j, 2, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(Spp1 detected)
            p[i, j, 3, 4] <- rdm[i, j, 3, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(Spp2 detected)
            p[i, j, 4, 4] <- rdm[i, j, 4, 4]/sum(rdm[i, j, 1:4, 4]) # -- Pr(Spp12 detected)
      
            y.hat[i,j,1:4] <- p[i, j, 1:4, z[i]]
            #'  Rank index position of p in descending order (i.e., identify index position of largest to smallest value of p) 
            #'  Then order data by rank (i.e, order index position based on ranking from highest to lowest value of p)
            y.hat.index[i, j, 1:4] <- order(-rank(p[i, j, 1:4, z[i]]))
            #'  Retain index position of the highest value of p
            y.hat.maxindex[i, j] <- y.hat.index[i, j, 1]

            #'  Calculate Chi-squared test statistic for observed & simulated data sets
            x2[i,j] <- pow((y[i,j] - y.hat.maxindex[i,j]), 2) / (y.hat.maxindex[i,j] + 0.0001)
            x2.sim[i,j] <- pow((y.sim[i,j] - y.hat.maxindex[i,j]), 2) / (y.hat.maxindex[i,j] + 0.0001)
            
            #'  Calculate Freeman-Tukey test statistic for observed & simulated data sets
            ft[i,j] <- pow((sqrt(y[i,j]) - sqrt(y.hat.maxindex[i,j])), 2) 
            ft.sim[i,j] <- pow((sqrt(y.sim[i,j]) - sqrt(y.hat.maxindex[i,j])), 2)
          }
      
        #'  Sum across surveys
        x2.obs[i] <- sum(x2[i,])
        x2.sims[i] <- sum(x2.sim[i,])
        
        ft.obs[i] <- sum(ft[i,])
        ft.sims[i] <- sum(ft.sim[i,])
  
        }
      
        #'  Sum across sites for final Chi-squared & Freeman-Tukey test statistics
        chi2.obs <- sum(x2.obs[])
        chi2.sim <- sum(x2.sims[])
        
        FT.obs <- sum(ft.obs[])
        FT.sims <- sum(ft.sims[])
      }
      ")