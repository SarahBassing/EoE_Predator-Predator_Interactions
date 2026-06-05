  #'  --------------------------------------------------------
  #'  JAGS model: top-down, interference SEM - no random effect for cluster
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
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_no_clusterRE.txt', "
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
      # for(s in 1:nWSI) {
      #   beta.wsi[s] ~ dnorm(0, 0.01)
      # }
      # for(f in 1:nforest) {
      #   beta.forest[f] ~ dnorm(0, 0.01)
      # }
      
      #'  SD prior for each regression and latent variables          
      for(k in 1:nSpp) {
        #'  time step t regression
        sigma.spp[k] ~ dunif(0, 10)
        tau.spp[k] <- 1 / pow(sigma.spp[k], 2)
        
        #'  time step t-1 lantent variables (i.e., random effect explaining 
        #'  variation around the true ecological state)
        sigma.spp.tmin1[k] ~ dunif(0, 10)
        tau.spp.tmin1[k] <- 1 / pow(sigma.spp.tmin1[k], 2)
      }
      
      
      #'  Likelihood
      #'  ----------
      #'  Measurement error from RN models for each species and cluster-level RDI
      #'  Posterior summaries (mean & sigma) treated as noisy observations [data] 
      #'  conditional on cluster-level latent parameter (truth) 
      for(i in 1:nCluster) {
        #'  RN model posterior means (spp.t_hat) arise from latent true RDI (spp.t).
        #'  RDI estimates from RN model are a noisy observation of true RDI,
        #'  governed by the true RDI and observed variability. True RDI (spp.t) are
        #'  determined by hypothesized ecological factors in likelihood.
        wolf.t_hat[i] ~ dnorm(wolf.t[i], wolf.t.tau_hat[i])
        wolf.tmin1_hat[i] ~ dnorm(wolf.tmin1[i], wolf.tmin1.tau_hat[i])
        lion.t_hat[i] ~ dnorm(lion.t[i], lion.t.tau_hat[i])
        lion.tmin1_hat[i] ~ dnorm(lion.tmin1[i], lion.tmin1.tau_hat[i])
        bear.t_hat[i] ~ dnorm(bear.t[i], bear.t.tau_hat[i])
        bear.tmin1_hat[i] ~ dnorm(bear.tmin1[i], bear.tmin1.tau_hat[i])
        coy.t_hat[i] ~ dnorm(coy.t[i], coy.t.tau_hat[i])
        coy.tmin1_hat[i] ~ dnorm(coy.tmin1[i], coy.tmin1.tau_hat[i])
        elk.t_hat[i] ~ dnorm(elk.t[i], elk.t.tau_hat[i])
        elk.tmin1_hat[i] ~ dnorm(elk.tmin1[i], elk.tmin1.tau_hat[i])
        moose.t_hat[i] ~ dnorm(moose.t[i], moose.t.tau_hat[i])
        moose.tmin1_hat[i] ~ dnorm(moose.tmin1[i], moose.tmin1.tau_hat[i])
        wtd.t_hat[i] ~ dnorm(wtd.t[i], wtd.t.tau_hat[i])
        wtd.tmin1_hat[i] ~ dnorm(wtd.tmin1[i], wtd.tmin1.tau_hat[i])
        
        #'  RN model posterior SD (spp.t.sigma_hat) used to calculate spp.t.tau_hat
        #'  spp.t.tau_hat are known constraints from the RN model posteriors that
        #'  inform how variable (noisy) the observed RDI can be given the latent truth
        wolf.t.tau_hat[i] <- 1 / pow(wolf.t.sigma_hat[i], 2)
        wolf.tmin1.tau_hat[i] <- 1 / pow(wolf.tmin1.sigma_hat[i], 2)
        lion.t.tau_hat[i] <- 1 / pow(lion.t.sigma_hat[i], 2)
        lion.tmin1.tau_hat[i] <- 1 / pow(lion.tmin1.sigma_hat[i], 2)
        bear.t.tau_hat[i] <- 1 / pow(bear.t.sigma_hat[i], 2)
        bear.tmin1.tau_hat[i] <- 1 / pow(bear.tmin1.sigma_hat[i], 2)
        coy.t.tau_hat[i] <- 1 / pow(coy.t.sigma_hat[i], 2)
        coy.tmin1.tau_hat[i] <- 1 / pow(coy.tmin1.sigma_hat[i], 2)
        elk.t.tau_hat[i] <- 1 / pow(elk.t.sigma_hat[i], 2)
        elk.tmin1.tau_hat[i] <- 1 / pow(elk.tmin1.sigma_hat[i], 2)
        moose.t.tau_hat[i] <- 1 / pow(moose.t.sigma_hat[i], 2)
        moose.tmin1.tau_hat[i] <- 1 / pow(moose.tmin1.sigma_hat[i], 2)
        wtd.t.tau_hat[i] <- 1 / pow(wtd.t.sigma_hat[i], 2)
        wtd.tmin1.tau_hat[i] <- 1 / pow(wtd.tmin1.sigma_hat[i], 2)
      }
      
      #'  Ecological process model
      #'  Latent cluster-level RDIs (spp.t) govern RN posterior summaries (spp.t_hat 
      #'  and spp.t.sigma_hat) and are in turn drawn from a normal distaribution
      #'  whose mean is defined by a species-specific autoregressive term, the 
      #'  RDIs of other species RDIs, and other variables.
      #'  A random effect is also included for repeat measures at the cluster-level.
      for(i in 1:nCluster) {
        wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
        mu.wolf.t[i] <- beta.int[1] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * harvest.tmin1[i]
        
        wolf.tmin1[i] ~ dnorm(mu.wolf.tmin1[i], tau.spp.tmin1[1])
        mu.wolf.tmin1[i] <- beta.int.tmin1[1]
        
        lion.t[i] ~ dnorm(mu.lion.t[i], tau.spp[2])
        mu.lion.t[i] <- beta.int[2] + beta.lion[1] * lion.tmin1[i] + beta.wolf[2] * wolf.tmin1[i] + beta.bear[2] * bear.tmin1[i]
        
        lion.tmin1[i] ~ dnorm(mu.lion.tmin1[i], tau.spp.tmin1[2])
        mu.lion.tmin1[i] <- beta.int.tmin1[2]

        bear.t[i] ~ dnorm(mu.bear.t[i], tau.spp[3])
        mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[3] * wolf.tmin1[i] 
        
        bear.tmin1[i] ~ dnorm(mu.bear.tmin1[i], tau.spp.tmin1[3]) 
        mu.bear.tmin1[i] <- beta.int.tmin1[3] 

        coy.t[i] ~ dnorm(mu.coy.t[i], tau.spp[4])
        mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[4] * wolf.tmin1[i] + beta.lion[2] * lion.tmin1[i] + beta.bear[3] * bear.tmin1[i] 
        
        coy.tmin1[i] ~ dnorm(mu.coy.tmin1[i], tau.spp.tmin1[4])
        mu.coy.tmin1[i] <- beta.int.tmin1[4] 

        elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
        mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[5] * wolf.tmin1[i] + beta.lion[3] * lion.tmin1[i] + beta.bear[4] * bear.tmin1[i]
        
        elk.tmin1[i] ~ dnorm(mu.elk.tmin1[i], tau.spp.tmin1[5])
        mu.elk.tmin1[i] <- beta.int.tmin1[5]

        moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
        mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[6] * wolf.tmin1[i] 
        
        moose.tmin1[i] ~ dnorm(mu.moose.tmin1[i], tau.spp.tmin1[6])
        mu.moose.tmin1[i] <- beta.int.tmin1[6]

        wtd.t[i] ~ dnorm(mu.wtd.t[i], tau.spp[7])
        mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[7] * wolf.tmin1[i] + beta.lion[4] * lion.tmin1[i] + beta.bear[5] * bear.tmin1[i] + beta.coy[2] * coy.tmin1[i] 
        
        wtd.tmin1[i] ~ dnorm(mu.wtd.tmin1[i], tau.spp.tmin1[7])
        mu.wtd.tmin1[i] <- beta.int.tmin1[7] 
        
      }
      
      
      #'  Derived parameters
      #'  ------------------
      
      
      }")