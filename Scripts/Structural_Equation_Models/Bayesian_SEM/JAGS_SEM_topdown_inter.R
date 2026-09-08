  #'  --------------------------------------------------------
  #'  JAGS model: top-down, interference SEM
  #'  
  #'  Model description: 
  #'    Structural equation model testing hypothesis that top-down processes and 
  #'    interference competition among predators most strongly determine the summer 
  #'    relative density indices of species in Northern Idaho's medium- and large-
  #'    bodied wildlife community. Model includes 1 year lag effect where the relative 
  #'    density of a species in the current time step [t] is affected by the relative 
  #'    density of itself, other species, and harvest from the previous time step [t-1].
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
  #'    beta.wolfharvest: effect of wolf harvest from previous time step
  #'    beta.lionharvest: effect of lion harvest from previous time step
  #'    beta.bearharvest: effect of bear harvest from previous time step
  #'    beta.elkharvest: effect of elk harvest from previous time step
  #'    beta.mooseharvest: effect of moose harvest from previous time step
  #'    beta.deerharvest: effect of deer harvest from previous time step
  #'    beta.wsi: effect of GMU-wide winter severity  from previous time step
  #'    beta.forest: effect of proportion of forest disturbance from previous 20 years
  #'    sigma.cluster: random effect for cluster (accounting for repeat measures across time step)
  #'  
  #'  Indices:
  #'    k: number of species (nSpp, 1:7), where
  #'       1 = wolf, 2 = mountain lion, 3 = black bear, 4 = coyote, 5 = elk, 6 = moose, 7 = white-tailed deer
  #'    t: number of time steps (nTimesetp, 1:3)
  #'    w: number of betas for 1-yr lag wolf effect (nWolf, 1:7), where
  #'       1 = 1L auto-regressive term in wolf regression, 2:7 = wolf effect on other species
  #'    l: number of betas for 1-yr lag mountain lion effect (nLion, 1:4, see w for details)
  #'    b: number of betas for 1-yr lag black bear effect (nBear, 1:2, see w for details)
  #'    c: number of betas for 1-yr lag coyote effect (nCoy, 1, see w for details)
  #'    e: number of betas for 1-yr lag elk effect (nElk, 1, see w for details)
  #'    m: number of betas for 1-yr lag moose effect (nMoose, 1, see w for details)
  #'    d: number of betas for 1-yr lag white-tailed deer effect (nDeer, 1, see w for details)
  #'    h: number of betas for 1-yr lag harvest effect (nharvest, 1, one per species except coy)
  #'    f: number of betas for 1-yr lag forest effect (nforest, 0)
  #'    s: number of betas for 1-yr lag winter severity effect (nWSI, 0)
  #'    i: number of clusters (nCluster, 1:23), where
  #'       each cluster represents the area over which species-specific RDI, harvest,
  #'       forest and WSI variables were generated from
  #'  --------------------------------------------------------
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter.txt', "
      model{
      
      #'  Define priors
      #'  -------------
      #'  Priors for intercepts
      #'  Use more informed prior for lion intercept
      beta.int[1] ~ dnorm(0, 1) # poor convergence with weaker priors
      beta.int.tmin1[1] ~ dnorm(0, 0.01)
      
      #'  Intercept priors for all other species (note, intercept prior for wolf is beta.int[2])
      for(k in 2:nSpp) {
        beta.int[k] ~ dnorm(0, 0.01) 
        beta.int.tmin1[k] ~ dnorm(0, 0.01)
      }
      
      #'  Priors for species lag effects
      #'  As a reminder: precision = 0.01 --> sqrt(0.01^-1) --> SD = 10
      for(l in 1:nLion) {
        beta.lion[l] ~ dnorm(0, 1)  # poor convergence with weaker priors
      }
      for(w in 1:nWolf) {
        beta.wolf[w] ~ dnorm(0, 0.1)  # using more informed prior to improve convergence
      }
      for(b in 1:nBear) {
        beta.bear[b] ~ dnorm(0, 0.01)
      }
      for(c in 1:nCoy) {
        beta.coy[c] ~ dnorm(0, 0.01)
      }
      for(e in 1:nElk) {
        beta.elk[e] ~ dnorm(0, 0.01)
      }
      for(m in 1:nMoose) {
        beta.moose[m] ~ dnorm(0, 0.01)
      }
      for(d in 1:nWtd) {
        beta.wtd[d] ~ dnorm(0, 0.01)
      }

      #'  Priors for anthropogenic & landscape effects
      for(h in 1:nharvest) {
        beta.harvest[h] ~ dnorm(0, 0.01)
      }
      
      #'  Half-Cauchy prior for latent variable SD
      #'  Define scale paramter (represents median of distribution)
      #'  Expresses prior belief of typical SD for latent variable
      scale.aux ~ dgamma(0.5, 0.5)
      scale ~ dnorm(0, 1 / (pow(2, 2) * scale.aux)) T(0,)
      for(k in 1:nSpp) {
        aux[k] ~ dgamma(0.5, 0.5)
        sigma.spp[k] ~ dnorm(0, 1 / (pow(scale, 2) * aux[k])) T(0,)
        tau.spp[k] <- 1 / pow(sigma.spp[k], 2)
      }
      
      #'  Likelihood
      #'  ----------
      #'  Measurement error from RN models for each species and cluster-level RDI
      #'  Posterior summaries (mean & sigma) treated as noisy observations [data] 
      #'  conditional on cluster-level latent parameter (truth) 
      for(i in 1:nSites) {
        for(y in 1:nYear) {
          lion.hat[i,y]  ~ dnorm(lion.latent[i,y],  lion.tau_hat[i,y])
          wolf.hat[i,y]  ~ dnorm(wolf.latent[i,y],  wolf.tau_hat[i,y])
          bear.hat[i,y]  ~ dnorm(bear.latent[i,y],  bear.tau_hat[i,y])
          coy.hat[i,y]   ~ dnorm(coy.latent[i,y],   coy.tau_hat[i,y])
          elk.hat[i,y]   ~ dnorm(elk.latent[i,y],   elk.tau_hat[i,y])
          moose.hat[i,y] ~ dnorm(moose.latent[i,y], moose.tau_hat[i,y])
          wtd.hat[i,y]   ~ dnorm(wtd.latent[i,y],   wtd.tau_hat[i,y])

          lion.tau_hat[i,y]  <- 1 / pow(lion.sigma_hat[i,y], 2)
          wolf.tau_hat[i,y]  <- 1 / pow(wolf.sigma_hat[i,y], 2)
          bear.tau_hat[i,y]  <- 1 / pow(bear.sigma_hat[i,y], 2)
          coy.tau_hat[i,y]   <- 1 / pow(coy.sigma_hat[i,y], 2)
          elk.tau_hat[i,y]   <- 1 / pow(elk.sigma_hat[i,y], 2)
          moose.tau_hat[i,y] <- 1 / pow(moose.sigma_hat[i,y], 2)
          wtd.tau_hat[i,y]   <- 1 / pow(wtd.sigma_hat[i,y], 2)
        }
      }
      
      #'  Ecological process model
      #'  Latent cluster-level RDIs (spp.true[i,y]) govern RN posterior summaries 
      #'  (spp.hat[i,y] and spp.sigma_hat[i,y]) and are in turn drawn from a normal 
      #'  distribution whose mean is defined by a species-specific autoregressive 
      #'  term, the RDIs of other species RDIs, and other variables. Year 1 has 
      #'  no prior year, so it gets an intercept-only baseline model.
      for(i in 1:nSites) {

        #'  Year 1: baseline latent states, no previous information available
        lion.latent[i,1]  ~ dnorm(beta.int.tmin1[1], tau.spp[1])
        wolf.latent[i,1]  ~ dnorm(beta.int.tmin1[2], tau.spp[2])
        bear.latent[i,1]  ~ dnorm(beta.int.tmin1[3], tau.spp[3])
        coy.latent[i,1]   ~ dnorm(beta.int.tmin1[4], tau.spp[4])
        elk.latent[i,1]   ~ dnorm(beta.int.tmin1[5], tau.spp[5])
        moose.latent[i,1] ~ dnorm(beta.int.tmin1[6], tau.spp[6])
        wtd.latent[i,1]   ~ dnorm(beta.int.tmin1[7], tau.spp[7])

        #'  Years 2-4: process model driven by the same previous-year latent
        #'  nodes used as outcomes above
        for(y in 2:nYear) {

          lion.latent[i,y] ~ dnorm(mu.lion[i,y], tau.spp[1])
          mu.lion[i,y] <- beta.int[1] + beta.harvest[1] * lionHarv[i,y-1] + beta.wolf[2] * wolf.latent[i,y-1] # + beta.lion[1] * lion.latent[i,y-1] 
             
          wolf.latent[i,y] ~ dnorm(mu.wolf[i,y], tau.spp[2])
          mu.wolf[i,y] <- beta.int[2] + beta.wolf[1] * wolf.latent[i,y-1] + beta.harvest[2] * wolfHarv[i,y-1] 

          bear.latent[i,y] ~ dnorm(mu.bear[i,y], tau.spp[3])
          mu.bear[i,y] <- beta.int[3] + beta.bear[1] * bear.latent[i,y-1] + beta.harvest[3] * bearHarv[i,y-1] + beta.wolf[3] * wolf.latent[i,y-1]
          
          coy.latent[i,y] ~ dnorm(mu.coy[i,y], tau.spp[4])
          mu.coy[i,y] <- beta.int[4] + beta.coy[1] * coy.latent[i,y-1] + beta.wolf[4] * wolf.latent[i,y-1] + beta.lion[2] * lion.latent[i,y-1]

          elk.latent[i,y] ~ dnorm(mu.elk[i,y], tau.spp[5])
          mu.elk[i,y] <- beta.int[5] + beta.elk[1] * elk.latent[i,y-1] + beta.wolf[5] * wolf.latent[i,y-1] + beta.lion[3] * lion.latent[i,y-1]  
          #  Removed bear - assuming wolves and lions are primary predators, bears are incidental and only affect neonates

          moose.latent[i,y] ~ dnorm(mu.moose[i,y], tau.spp[6])
          mu.moose[i,y] <- beta.int[6] + beta.moose[1] * moose.latent[i,y-1] + beta.wolf[6] * wolf.latent[i,y-1] 

          wtd.latent[i,y] ~ dnorm(mu.wtd[i,y], tau.spp[7])
          mu.wtd[i,y] <- beta.int[7] + beta.wtd[1] * wtd.latent[i,y-1] + beta.lion[4] * lion.latent[i,y-1] 
          #  Removed wolf, bear, and coy - assuming lions are primary predator, all others or incidental
      
        }
      }
      
      
      #'  Derived parameters
      #'  ------------------
      #'  d-Separation...
      
      #'  Total and indirect effects...
      
    }")
  
 