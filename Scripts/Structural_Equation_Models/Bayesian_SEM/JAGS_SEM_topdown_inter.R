#'  --------------------------------------------------------
#'  JAGS model: top-down, interference SEM
#'  
#'  Model description: 
#'    Structural equation model testing hypothesis that top-down processes and 
#'    interference competition among predators most strongly determine the summer 
#'    relative density indices of species in Northern Idaho's medium- and large-
#'    bodied wildlife community. Model includes 1 year lag effect where the relative 
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
#'    t: number of time steps (nTimesetp, 1:3)
#'    w: number of betas for 1-yr lag wolf effect (nWolf, 1:7), where
#'       1 = 1L auto-regressive term in wolf regression, >1 = wolf effect on other species
#'    l: number of betas for 1-yr lag mountain lion effect (nLion, 1:4, see w for details)
#'    b: number of betas for 1-yr lag black bear effect (nBear, 1:5, see w for details)
#'    c: number of betas for 1-yr lag coyote effect (nCoy, 1:2, see w for details)
#'    e: number of betas for 1-yr lag elk effect (nElk, 1, see w for details)
#'    m: number of betas for 1-yr lag moose effect (nMoose, 1, see w for details)
#'    d: number of betas for 1-yr lag white-tailed deer effect (nDeer, 1, see w for details)
#'    h: number of betas for 1-yr lag harvest effect (nharvest, 1)
#'    f: number of betas for 1-yr lag forest effect (nforest, 0)
#'    i: number of clusters (nCluster, 1:24), where
#'       each cluster represents the area over which species-specific RDI, harvest,
#'       and forest variables were generated from
#'  --------------------------------------------------------
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter.txt', "
      model{
      
      #'  Define priors
      #'  -------------
      #'  Priors for intercepts
      for(k in 1:nSpp) {
        beta.int[k] ~ dnorm(0, 0.01)  
      }
      
      #'  Priors for species lag effects
      for(w in 1:nWolf) {
        beta.wolf[w] ~ dnorm(0, 0.01)  # precision = 0.01 --> sqrt(0.01^-1) --> SD = 10
      }
      for(l in 1:nLion) {
        beta.lion[l] ~ dnorm(0, 0.01)
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
      # for(f in 1:nforest) {
      #   beta.forest[f] ~ dnorm(0, 0.01)
      # }
      
      #'  SD prior for each regression
      for(k in 1:nSpp) {
        sigma.spp[k] ~ dunif(0, 10)
        tau.spp[k] <- 1 / pow(sigma.spp[k], 2)
      }
      
      #'  SD prior for among-cluster random intercept term for each species 
      for(k in 1:nSpp) {
        for(i in 1:nCluster) {
          sigma.cluster[k,i] ~ dnorm(0, 0.01) T(0,)      #dunif(0, 10)      
          tau.cluster[k,i] <- 1 / pow(sigma.cluster[k,i], 2)    
        }
      }
      #'  dnorm(0, 0.01) T(0,) some seem to hitting boundary at 0 and a few not converging well 
      #'  dunif(0, 10) most are just returning the prior, a few have poor convergence
      #'  dnorm seems to be the better prior if sticking with a species - cluster random effect
      
      #'  Likelihood
      #'  ----------
      for(i in 1:nCluster) {
        wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
        mu.wolf.t[i] <- beta.int[1] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * harvest.tmin1[i] + tau.cluster[1,i] 
        
        lion.t[i] ~ dnorm(mu.lion.t[i], tau.spp[2])
        mu.lion.t[i] <- beta.int[2] + beta.lion[1] * lion.tmin1[i] + beta.wolf[2] * wolf.tmin1[i] + beta.bear[2] * bear.tmin1[i] + tau.cluster[2,i]

        bear.t[i] ~ dnorm(mu.bear.t[i], tau.spp[3])
        mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[3] * wolf.tmin1[i] + tau.cluster[3,i]

        coy.t[i] ~ dnorm(mu.coy.t[i], tau.spp[4])
        mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[4] * wolf.tmin1[i] + beta.lion[2] * lion.tmin1[i] + beta.bear[3] * bear.tmin1[i] + tau.cluster[4,i]

        elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
        mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[5] * wolf.tmin1[i] + beta.lion[3] * lion.tmin1[i] + beta.bear[4] * bear.tmin1[i] + tau.cluster[5,i]

        moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
        mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[6] * wolf.tmin1[i] + tau.cluster[6,i]

        wtd.t[i] ~ dnorm(mu.wtd.t[i], tau.spp[7])
        mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[7] * wolf.tmin1[i] + beta.lion[4] * lion.tmin1[i] + beta.bear[5] * bear.tmin1[i] + beta.coy[2] * coy.tmin1[i] + tau.cluster[7,i]

      }
      
      
      #'  Derived parameters
      #'  ------------------
      
      
      }")