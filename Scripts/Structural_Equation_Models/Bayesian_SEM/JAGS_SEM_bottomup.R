  #'  --------------------------------------------------------
  #'  JAGS model: bottom-up SEM
  #'  
  #'  Model description: 
  #'    Structural equation model testing hypothesis that top-down processes most 
  #'    strongly determine the summer relative density indices of species in Northern 
  #'    Idaho's medium- and large-bodied wildlife community. Model includes 1 year 
  #'    lag effect where the relative density of a species in the current time 
  #'    step [t] is affected by the relative density of itself, other species, 
  #'    and harvest from the previous time step [t-1].
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
  #'    w: number of betas for 1-yr lag wolf effect (nWolf, 1:4), where
  #'       1 = 1L auto-regressive term in wolf regression, 2:4 = wolf effect on other species
  #'    l: number of betas for 1-yr lag mountain lion effect (nLion, 1, see w for details)
  #'    b: number of betas for 1-yr lag black bear effect (nBear, 1, see w for details)
  #'    c: number of betas for 1-yr lag coyote effect (nCoy, 1, see w for details)
  #'    e: number of betas for 1-yr lag elk effect (nElk, 1:5, see w for details)
  #'    m: number of betas for 1-yr lag moose effect (nMoose, 1:2, see w for details)
  #'    d: number of betas for 1-yr lag white-tailed deer effect (nDeer, 1:3, see w for details)
  #'    h: number of betas for 1-yr lag harvest effect (nharvest, 0)
  #'    f: number of betas for 1-yr lag forest effect (nforest, 1:4)
  #'    s: number of betas for 1-yr lag winter severity effect (nWSI, 1:3)
  #'    i: number of clusters (nCluster, 1:23), where
  #'       each cluster represents the area over which species-specific RDI, harvest,
  #'       forest and WSI variables were generated from
  #'  --------------------------------------------------------
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_bottomup.txt', "
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
      # for(h in 1:nharvest) {
      #   beta.harvest[h] ~ dnorm(0, 0.01)
      # }
      for(s in 1:nWSI) {
        beta.wsi[s] ~ dnorm(0, 0.01)
      }
      for(f in 1:nforest) {
        beta.forest[f] ~ dnorm(0, 0.01)
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
      
      #' #'  Alternate prior for latent variable SD
      #' #'  Define scale paramter (represents median of distribution)
      #' #'  Expresses prior belief of typical SD for latent variable
      #' scale <- 1  #  keep scale of standardized RDIs in mind
      #' for(k in 1:nSpp) {
      #'   #'  Draw random value from gamma distribution
      #'   aux[k] ~ dgamma(0.5, 0.5)
      #'   #'  Draw random value from normal distribution using squared scale * aux
      #'   #'  value as the variance, 1/var = precision. Then truncating at zero to
      #'   #'  create a Half-Cauchy distribution.
      #'   #'  time step t regression
      #'   sigma.spp[k] ~ dnorm(0, 1 / (pow(scale, 2) * aux[k])) T(0,)
      #'   tau.spp[k] <- 1 / pow(sigma.spp[k], 2)
      #'   #'  time step t-1 lantent variables (i.e., random effect explaining
      #'   #'  variation around the true ecological state)
      #'   sigma.spp.tmin1[k] ~ dnorm(0, 1 / (pow(scale, 2) * aux[k])) T(0,)
      #'   tau.spp.tmin1[k] <- 1 / pow(sigma.spp.tmin1[k], 2)
      #' }
      
      
      
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
          mu.lion[i,y] <- beta.int[1] + beta.elk[2] * elk.latent[i,y-1] + beta.wtd[2] * wtd.latent[i,y-1] # + beta.lion[1] * lion.latent[i,y-1]
             
          wolf.latent[i,y] ~ dnorm(mu.wolf[i,y], tau.spp[2])
          mu.wolf[i,y] <- beta.int[2] + beta.wolf[1] * wolf.latent[i,y-1] + beta.elk[3] * elk.latent[i,y-1] + beta.moose[2] * moose.latent[i,y-1]

          bear.latent[i,y] ~ dnorm(mu.bear[i,y], tau.spp[3])
          mu.bear[i,y] <- beta.int[3] + beta.bear[1] * bear.latent[i,y-1] + beta.elk[4] * elk.latent[i,y-1] + beta.forest[4] * forest[i,y-1]
          
          coy.latent[i,y] ~ dnorm(mu.coy[i,y], tau.spp[4])
          mu.coy[i,y] <- beta.int[4] + beta.coy[1] * coy.latent[i,y-1] + beta.wtd[3] * wtd.latent[i,y-1] 
          #'  Dropped elk - assuming only neonates have any affect on coyotes and only for short period of time

          elk.latent[i,y] ~ dnorm(mu.elk[i,y], tau.spp[5])
          mu.elk[i,y] <- beta.int[5] + beta.elk[1] * elk.latent[i,y-1] + beta.forest[1] * forest[i,y-1] + beta.wsi[1] * wsi[i,y-1]

          moose.latent[i,y] ~ dnorm(mu.moose[i,y], tau.spp[6])
          mu.moose[i,y] <- beta.int[6] + beta.moose[1] * moose.latent[i,y-1] + beta.forest[2] * forest[i,y-1] + beta.wsi[2] * wsi[i,y-1]

          wtd.latent[i,y] ~ dnorm(mu.wtd[i,y], tau.spp[7])
          mu.wtd[i,y] <- beta.int[7] + beta.wtd[1] * wtd.latent[i,y-1] + beta.forest[3] * forest[i,y-1] + beta.wsi[3] * wsi[i,y-1] 
      
        }
      }
      
      
      
      #'  Derived parameters
      #'  ------------------
      #'  d-Separation...
      
      #'  Total and indirect effects...
      
      }")
  
  # for(i in 1:nCluster) {
  #   wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
  #   mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.elk[2] * elk.tmin1[i] + beta.moose[2] * moose.tmin1[i] 
  #   
  #   wolf.tmin1[i] ~ dnorm(mu.wolf.tmin1[i], tau.spp.tmin1[1])
  #   mu.wolf.tmin1[i] <- beta.int.tmin1[2]
  #   
  #   lion.t[i] ~ dnorm(mu.lion.t[i], tau.spp[2])  
  #   mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.elk[3] * elk.tmin1[i] + beta.wtd[2] * wtd.tmin1[i]
  #   
  #   lion.tmin1[i] ~ dnorm(mu.lion.tmin1[i], tau.spp.tmin1[2])
  #   mu.lion.tmin1[i] <- beta.int.tmin1[1] 
  #   
  #   bear.t[i] ~ dnorm(mu.bear.t[i], tau.spp[3]) 
  #   mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.elk[4] * elk.tmin1[i] + beta.forest[1] * forest.tmin1[i]
  #   
  #   bear.tmin1[i] ~ dnorm(mu.bear.tmin1[i], tau.spp.tmin1[3]) 
  #   mu.bear.tmin1[i] <- beta.int.tmin1[3] 
  #   
  #   coy.t[i] ~ dnorm(mu.coy.t[i], tau.spp[4])
  #   mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wtd[3] * wtd.tmin1[i]
  #   #  Removed beta.elk[5] * elk.tmin1[i] + 
  #   #  Assuming only neonates have any affect on coyotes and only for short period of time
  #   
  #   coy.tmin1[i] ~ dnorm(mu.coy.tmin1[i], tau.spp.tmin1[4])
  #   mu.coy.tmin1[i] <- beta.int.tmin1[4] 
  #   
  #   elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
  #   mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.forest[2] * forest.tmin1[i] + beta.wsi[1] * wsi.tmin1[i]
  #   
  #   elk.tmin1[i] ~ dnorm(mu.elk.tmin1[i], tau.spp.tmin1[5])
  #   mu.elk.tmin1[i] <- beta.int.tmin1[5]
  #   
  #   moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
  #   mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.forest[3] * forest.tmin1[i] + beta.wsi[2] * wsi.tmin1[i]
  #   
  #   moose.tmin1[i] ~ dnorm(mu.moose.tmin1[i], tau.spp.tmin1[6])
  #   mu.moose.tmin1[i] <- beta.int.tmin1[6]
  #   
  #   wtd.t[i] ~ dnorm(mu.wtd.t[i], tau.spp[7])
  #   mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.forest[4] * forest.tmin1[i] + beta.wsi[3] * wsi.tmin1[i]
  #   
  #   wtd.tmin1[i] ~ dnorm(mu.wtd.tmin1[i], tau.spp.tmin1[7])
  #   mu.wtd.tmin1[i] <- beta.int.tmin1[7] 
  #   
  # }