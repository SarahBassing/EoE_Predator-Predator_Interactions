  #'  --------------------------------------------------------
  #'  JAGS model: top-down, interference SEM with cross lags
  #'  
  #'  Model description: 
  #'    Structural equation model testing hypothesis that top-down processes and 
  #'    interference competition among predators most strongly determine the summer 
  #'    relative density indices of species in Northern Idaho's medium- and large-
  #'    bodied wildlife community. Model includes cross lags where the relative 
  #'    density of a species in the current time step [t] is affected by the 
  #'    relative density of itself and other species from the previous time step [t-1].
  #'  
  #'  Parameters:
  #'    beta.int: intercept for each regression
  #'    beta.wolf: effect of wolf relative density index from previous time step
  #'    beta.lion: effect of mountain lion relative density index from previous time step
  #'    beta.bear: effect of black bear relative density index from previous time step
  #'    beta.coy: effect of coyote relative density index from previous time step
  #'    beta.elk: effect of elk relative density index from previous time step
  #'    beta.moose: effect of moose relative density index from previous time step
  #'    beta.wtd: effect of white-tailed deer relative density index from previous time step
  #'    beta.harvest: effect of wolf harvest from previous time step
  #'    beta.forest: effect of forest disturbance from previous time step
  #'    sigma.cluster: random effect for cluster (accounting for repeat measures across time step)
  #'  
  #'  Indices:
  #'    k: number of species (nSpp, 1:7), where
  #'       1 = wolf, 2 = mountain lion, 3 = black bear, 4 = coyote, 5 = elk, 6 = moose, 7 = white-tailed deer
  #'    mod: number of regressions (1 for t1, 2 for t>1)
  #'    t: number of time steps (nTimesetp, 1:3)
  #'    w: number of betas for lag wolf effect (nWolf, 1:7), where
  #'       1 = 1L auto-regressive term in wolf regression, >1 = wolf effect on other species
  #'    l: number of betas for lag mountain lion effect (nLion, 1:4, see w for details)
  #'    b: number of betas for lag black bear effect (nBear, 1:5, see w for details)
  #'    c: number of betas for lag coyote effect (nCoy, 1:2, see w for details)
  #'    e: number of betas for lag elk effect (nElk, 1, see w for details)
  #'    m: number of betas for lag moose effect (nMoose, 1, see w for details)
  #'    d: number of betas for lag white-tailed deer effect (nDeer, 1, see w for details)
  #'    h: number of betas for lag harvest effect (nharvest, 1)
  #'    f: number of betas for lag forest effect (nforest, 0)
  #'    cl: number of clusters (nCluster, 1:24), where
  #'       each cluster represents the area over which species-specific RDI, harvest,
  #'       and forest variables were generated from
  #'  --------------------------------------------------------
  
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_crosslag.txt', "
        model{
      
      #'  Define priors
      #'  -------------
      #'  Intercept for each species-specific regression 
      #'  Need two per species because regressions differ between t1 and t>1
      for(k in 1:nSpp) {
        for(mod in 1:2) {
          beta.int[k,mod] ~ dnorm(0, 0.1)  # look more into precision value... I think FSP review commented on this?
        }
      }
      
      #'  Standard deviation for each timestep and species-specific regression  
      for(k in 1:nSpp) {
        for(t in 1:nTimestep) {
          sigma.spp[k,t] ~ dunif(0, 10)
          tau.spp[k,t] <- 1 / pow(sigma.spp[k,t], 2)
        }
      }
      
      ####  DO I REALLY WANT SIGMA.SPP TO VARY WITH TIMESTEP TOO????
      ####  Currently a different estiamted SD for year 1, 2, and 3 regressions
      ####  but do I need that to vary by year if the observation data come from
      ####  the overall same process and their just 3 realizations of that process? 
      ####  But the regressions differ each year based on different parameters
      ####  (t1 vs t>1) and different observations each year so the expected RDI
      ####  (mu.spp) differs each year and so should have a different est. SD, right?
      
      
      #'  Slopes representing species-specific, human, and landscape lag effects
      for(w in 1:nWolf) {
        beta.wolf[w] ~ dnorm(0, 0.1)  # look more into precision value... I think FSP review commented on this?
      }
      for(l in 1:nLion) {
        beta.lion[l] ~ dnorm(0, 0.1)  # 0.01 originally
      }
      for(b in 1:nBear) {
        beta.bear[b] ~ dnorm(0, 0.1)
      }
      for(c in 1:nCoy) {
        beta.coy[c] ~ dnorm(0, 0.1)
      }
      for(e in 1:nElk) {
        beta.elk[e] ~ dnorm(0, 0.1)
      }
      for(m in 1:nMoose) {
        beta.moose[m] ~ dnorm(0, 0.1)
      }
      for(d in 1:nDeer) {
        beta.wtd[d] ~ dnorm(0, 0.1)
      }
      for(h in 1:nharvest) {
        beta.harvest[h] ~ dnorm(0, 0.1)
      }
      for(f in 1:nforest) {
        beta.forest[f] ~ dnorm(0, 0.1)
      }
      
      #'  Standard deviation for among-cluster random intercept term for each species 
      for(k in 1:nSpp) {
        for(cl in 1:nCluster) {
          sigma.cluster[k,cl] ~ dnorm(0, 0.1) T(0,)     #dunif(0, 10)       
          tau.cluster[k,cl] <- 1 / pow(sigma.cluster[k,cl], 2)    
        }
      }
      
      #### DO I WANT A RANDOM EFFECT FOR EACH SPP AND CLUSTER OR JUST ONE PER CLUSTER??? 
      #### I mainly want to account for repeat measure of clusters across years, right?
      #### Do I expect the impact of repeart measures vary by species?
      
      
      
      #'  Likelihood
      #'  ----------
      for(i in 1:nCluster) {
        #'  If t = 1
        wolf[i,1] ~ dnorm(mu.wolf[i,1], tau.spp[1,1])      # Note mu.wolf[i,1] & beta.wolf[1] are different from those in later timesteps
        mu.wolf[i,1] <- beta.int[1,1] + tau.cluster[1,i]   # Double check indexing on tau.cluster[spp,cluster]
          
        lion[i,1] ~ dnorm(mu.lion[i,1], tau.spp[2,1])
        mu.lion[i,1] <- beta.int[2,1] + tau.cluster[2,i]
        
        bear[i,1] ~ dnorm(mu.bear[i,1], tau.spp[3,1])
        mu.bear[i,1] <- beta.int[3,1] + tau.cluster[3,i]
        
        coy[i,1] ~ dnorm(mu.coy[i,1], tau.spp[4,1])
        mu.coy[i,1] <- beta.int[4,1] + tau.cluster[4,i]
        
        elk[i,1] ~ dnorm(mu.elk[i,1], tau.spp[5,1])
        mu.elk[i,1] <- beta.int[5,1] + tau.cluster[5,i]
        
        moose[i,1] ~ dnorm(mu.moose[i,1], tau.spp[6,1])
        mu.moose[i,1] <- beta.int[6,1] + tau.cluster[6,i]
        
        wtd[i,1] ~ dnorm(mu.wtd[i,1], tau.spp[7,1])
        mu.wtd[i,1] <- beta.int[7,1] + tau.cluster[7,i]
          
          
        #'  If t > 1              
        #'  NOTE: spp-specific betas are constant for t>1 (i.e., the effect of 
        #'  lag wolf rdi on lion rdi is assumed to be the same across years)
        for(t in 2:nTimestep) {
          wolf[i,t] ~ dnorm(mu.wolf[i,t], tau.spp[1,t])
          mu.wolf[i,t] <- beta.int[1,2] + beta.wolf[1] * wolf[i,t-1] + beta.harvest[1] * harv[i,t-1] + tau.cluster[1,i]

          lion[i,t] ~ dnorm(mu.lion[i,t], tau.spp[2,t])
          mu.lion[i,t] <- beta.int[2,2] + beta.lion[1] * lion[i,t-1] + beta.wolf[2] * wolf[i,t-1] + beta.bear[2] * bear[i,t-1] + tau.cluster[2,i]

          bear[i,t] ~ dnorm(mu.bear[i,t], tau.spp[3,t])
          mu.bear[i,t] <- beta.int[3,2] + beta.bear[1] * bear[i,t-1] + beta.wolf[3] * wolf[i,t-1] + tau.cluster[3,i]
          
          coy[i,t] ~ dnorm(mu.coy[i,t], tau.spp[4,t])
          mu.coy[i,t] <- beta.int[4,2] + beta.coy[1] * coy[i,t-1] + beta.wolf[4] * wolf[i,t-1] + beta.lion[2] * lion[i,t-1] + beta.bear[3] * bear[i,t-1] + tau.cluster[4,i]
        
          elk[i,t] ~ dnorm(mu.elk[i,t], tau.spp[5,t])
          mu.elk[i,t] <- beta.int[5,2] + beta.elk[1] * elk[i,t-1] + beta.wolf[5] * wolf[i,t-1] + beta.lion[3] * lion[i,t-1] + beta.bear[4] * bear[i,t-1] + tau.cluster[5,i]
          
          moose[i,t] ~ dnorm(mu.moose[i,t], tau.spp[6,t])
          mu.moose[i,t] <- beta.int[6,2] + beta.moose[1] * moose[i,t-1] + beta.wolf[6] * wolf[i,t-1] + tau.cluster[6,i]
          
          wtd[i,t] ~ dnorm(mu.wtd[i,t], tau.spp[7,t])
          mu.wtd[i,t] <- beta.int[7,2] + beta.wtd[1] * wtd[i,t-1] + beta.wolf[7] * wolf[i,t-1] + beta.lion[4] * lion[i,t-1] + beta.bear[5] * bear[i,t-1] + beta.coy[2] * coy[i,t-1] + tau.cluster[7,i]
        
        }
        
      }
      
      
      #'  Derived parameters
      #'  ------------------
      
      
      }")