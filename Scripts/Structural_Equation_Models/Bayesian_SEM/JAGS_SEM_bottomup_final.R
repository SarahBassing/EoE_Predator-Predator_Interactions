  #'  --------------------------------------------------------
  #'  JAGS model: bottom-up exploitative competition SEM finalized after d-Sep
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
  #'    w: number of betas for 1-yr lag wolf effect (nWolf), where
  #'       1 = 1L auto-regressive term in wolf regression, 2:4 = wolf effect on other species
  #'    l: number of betas for 1-yr lag mountain lion effect (nLion, see w for details)
  #'    b: number of betas for 1-yr lag black bear effect (nBear, see w for details)
  #'    c: number of betas for 1-yr lag coyote effect (nCoy, see w for details)
  #'    e: number of betas for 1-yr lag elk effect (nElk, see w for details)
  #'    m: number of betas for 1-yr lag moose effect (nMoose, see w for details)
  #'    d: number of betas for 1-yr lag white-tailed deer effect (nDeer, see w for details)
  #'    f: number of betas for 1-yr lag forest effect (nforest, 1:4)
  #'    s: number of betas for 1-yr lag winter severity effect (nWSI, 1:3)
  #'    i: number of clusters (nCluster, 1:23), where
  #'       each cluster represents the area over which species-specific RDI, harvest,
  #'       forest and WSI variables were generated from
  #'  --------------------------------------------------------
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_final.txt', "
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
          mu.lion[i,y] <- beta.int[1] + beta.elk[2] * elk.latent[i,y-1] + beta.elk[3] * elk.latent[i,y] + beta.wtd[2] * wtd.latent[i,y-1] + beta.wtd[3] * wtd.latent[i,y] 
             
          wolf.latent[i,y] ~ dnorm(mu.wolf[i,y], tau.spp[2])
          mu.wolf[i,y] <- beta.int[2] + beta.wolf[1] * wolf.latent[i,y-1] + beta.elk[4] * elk.latent[i,y-1] + beta.elk[5] * elk.latent[i,y] + beta.moose[2] * moose.latent[i,y-1] + beta.moose[3] * moose.latent[i,y] + beta.wtd[4] * wtd.latent[i,y] 

          bear.latent[i,y] ~ dnorm(mu.bear[i,y], tau.spp[3])
          mu.bear[i,y] <- beta.int[3] + beta.bear[1] * bear.latent[i,y-1] + beta.elk[6] * elk.latent[i,y-1] + beta.forest[4] * forest[i,y-1] + beta.wtd[5] * wtd.latent[i,y]
          
          coy.latent[i,y] ~ dnorm(mu.coy[i,y], tau.spp[4])
          mu.coy[i,y] <- beta.int[4] + beta.coy[1] * coy.latent[i,y-1] + beta.wtd[6] * wtd.latent[i,y-1] + beta.wtd[7] * wtd.latent[i,y]
          
          elk.latent[i,y] ~ dnorm(mu.elk[i,y], tau.spp[5])
          mu.elk[i,y] <- beta.int[5] + beta.elk[1] * elk.latent[i,y-1] + beta.forest[1] * forest[i,y-1] + beta.wsi[1] * wsi[i,y-1] + beta.bear[2] * bear.latent[i,y-1]

          moose.latent[i,y] ~ dnorm(mu.moose[i,y], tau.spp[6])
          mu.moose[i,y] <- beta.int[6] + beta.moose[1] * moose.latent[i,y-1] + beta.forest[2] * forest[i,y-1] + beta.wsi[2] * wsi[i,y-1]

          wtd.latent[i,y] ~ dnorm(mu.wtd[i,y], tau.spp[7])
          mu.wtd[i,y] <- beta.int[7] + beta.wtd[1] * wtd.latent[i,y-1] + beta.forest[3] * forest[i,y-1] + beta.wsi[3] * wsi[i,y-1] + beta.coy[2] * coy.latent[i,y-1] + beta.bear[3] * bear.latent[i,y-1]
      
        }
      }
      
      
      #'  Derived parameters: Indirect Effects
      #'  ------------------------------------
      #'  Exploitation competition where predator.tmin1 --> prey.t --> predator.t
      #'  Via elk.........................................
      indirect.bear.elk.lion <- beta.bear[2] * beta.elk[3]     # effect of bear.tmin1 on elk.t --> effect of elk.t on lion.t (contemporaneous)
      indirect.bear.elk.wolf <- beta.bear[2] * beta.elk[5]     # effect of bear.tmin1 on elk.t --> effect of elk.t on wolf.t (contemporaneous)
      
      #'  Via white-tailed deer...........................
      indirect.coy.wtd.lion <- beta.coy[2] * beta.wtd[3]       # effect of coy.tmin1 on wtd.t --> effect of wtd.t on lion.t (contemporaneous)
      indirect.coy.wtd.wolf <- beta.coy[2] * beta.wtd[4]       # effect of coy.tmin1 on wtd.t --> effect of wtd.t on wolf.t (contemporaneous)
      indirect.coy.wtd.bear <- beta.coy[2] * beta.wtd[5]       # effect of coy.tmin1 on wtd.t --> effect of wtd.t on bear.t (contemporaneous)
      indirect.coy.wtd.coy <- beta.coy[2] * beta.wtd[7]        # effect of coy.tmin1 on wtd.t --> effect of wtd.t on coy.t (contemporaneous, coy depleting coy food over time)
      indirect.bear.wtd.lion <- beta.bear[3] * beta.wtd[3]     # effect of bear.tmin1 on wtd.t --> effect of wtd.t on lion.t (contemporaneous)
      indirect.bear.wtd.wolf <- beta.bear[3] * beta.wtd[4]     # effect of bear.tmin1 on wtd.t --> effect of wtd.t on wolf.t (contemporaneous)
      indirect.bear.wtd.bear <- beta.bear[3] * beta.wtd[5]     # effect of bear.tmin1 on wtd.t --> effect of wtd.t on bear.t (contemporaneous)  NOTE: different mechanism but same outcome as indirect.wtd.bear.wtd
      indirect.bear.wtd.coy <- beta.bear[3] * beta.wtd[7]      # effect of bear.tmin1 on wtd.t --> effect of wtd.t on coy.t (contemporaneous)
      
      #'  Exploitation competition where predator.tmin2 --> prey.tmin1 --> predator.t
      #'  Via elk...........................................
      indirect.bear.elklag.lion <- beta.bear[2] * beta.elk[2]  # effect of bear.tmin2 on elk.tmin1 --> effect of elk.tmin1 on lion.t
      indirect.bear.elklag.wolf <- beta.bear[2] * beta.elk[4]  # effect of bear.tmin2 on elk.tmin1 --> effect of elk.tmin1 on wolf.t
      indirect.bear.elklag.bear <- beta.bear[2] * beta.elk[6]  # effect of bear.tmin2 on elk.tmin1 --> effect of elk.tmin1 on bear.t  (bear depleting bear food over time) NOTE: different mechanism but same outcome as indirect.elk.bearlag.elk
      
      #'  Via white-tailed deer..............................
      indirect.coy.wtdlag.lion <- beta.coy[2] * beta.wtd[2]    # effect of coy.tmin2 on wtd.tmin1 --> effect of wtd.tmin1 on lion.t
      indirect.coy.wtdlag.coy <- beta.coy[2] * beta.wtd[6]     # effect of coy.tmin2 on wtd.tmin1 --> effect of wtd.tmin1 on coy.t  (coy depleting coy food over time)
      indirect.bear.wtdlag.lion <- beta.bear[3] * beta.wtd[2]  # effect of bear.tmin2 on wtd.tmin1 --> effect of wtd.tmin1 on lion.t
      indirect.bear.wtdlag.coy <- beta.bear[3] * beta.wtd[6]   # effect of bear.tmin2 on wtd.tmin1 --> effect of wtd.tmin1 on coy.t  
      
      #'  2-yr lag effect on same species where predator effect on prey is extended via prey's AR1 term where predator.tmin2 --> prey.tmin1 --> prey.t
      #'  Via prey's AR1 term.............................
      indirect.bear.elk.self <- beta.bear[2] * beta.elk[1]     # effect of bear.tmin2 on elk.tmin1 --> effect of elk.tmin1 on elk.t   
      indirect.bear.wtd.self <- beta.bear[3] * beta.wtd[1]     # effect of bear.tmin2 on wtd.tmin1 --> effect of wtd.tmin1 on wtd.t
      indirect.coy.wtd.self <- beta.coy[2] * beta.wtd[1]       # effect of coy.tmin2 on wtd.tmin1 --> effect of wtd.tmin1 on wtd.t
      
      #'  2-yr lag effect on same species where prey effect on predator is extended via predator's AR1 term where prey.tmin2 --> predator.tmin1 --> predator.t
      #'  Via predator's AR1 term.........................
      indirect.elk.wolf.self <- beta.elk[4] * beta.wolf[1]     # effect of elk.tmin2 on wolf.tmin1 --> effect of wolf.tmin1 on wolf.t
      indirect.moose.wolf.self <- beta.moose[2] * beta.wolf[1] # effect of moose.tmin2 on wolf.tmin1 --> effect of wolf.tmin1 on wolf.t
      indirect.wtd.wolf.self <- beta.wtd[4] * beta.wolf[1]     # effect of wtd.tmin1 on wolf.tmin1 (contemporaneous) --> effect of wolf.tmin1 on wolf.t
      indirect.elk.bear.self <- beta.elk[6] * beta.bear[1]     # effect of elk.tmin2 on bear.tmin1 --> effect of bear.tmin1 on bear.t
      indirect.wtd.bear.self <- beta.wtd[5] * beta.bear[1]     # effect of wtd.tmin1 on bear.tmin1 (contemporaneous) --> effect of bear.tmin1 on bear.t
      indirect.wtd.coy.self <- beta.wtd[6] * beta.coy[1]       # effect of wtd.tmin2 on coy.tmin1 --> effect of coy.tmin1 on coy.t
      
      #'  Apparent competition where prey.A.tmin2 --> predator.tmin1 --> prey.B.t
      #'  Via black bear.............................
      indirect.elk.bearlag.elk <- beta.elk[6] * beta.bear[2]   # effect of elk.tmin2 on bear.tmin1 --> effect of bear.tmin1 on elk.t (elk affecting itself via effect on bears) NOTE: different mechanism but same outcome as indirect.bear.elklag.bear
      indirect.elk.bearlag.wtd <- beta.elk[6] * beta.bear[3]   # effect of elk.tmin2 on bear.tmin1 --> effect of bear.tmin1 on wtd.t
      indirect.wtd.bear.elk <- beta.wtd[5] * beta.bear[2]      # effect of wtd.tmin1 on bear.tmin1 (contemporaneous) --> effect of bear.tmin1 on elk.t
      # indirect.wtd.bear.wtd <- beta.wtd[5] * beta.bear[3]      # effect of wtd.tmin1 on bear.tmin1 (contemporaneous) --> effect of wtd.tmin1 on wtd.t (wtd affecting itself via effect on bears) NOTE: different mechanism but same outcome as indirect.bear.wtd.bear
      
      #'  Derived parameters: Sum Effects
      #'  -------------------------------
      #'  Summed lagged and contemporaneous direct effects, ignoring paths through other species or its own persistence
      sum.elk.lion <- beta.elk[2] + beta.elk[3]         # effect of elk.tmin1 on lion.t AND effect of elk.t on lion.t
      sum.elk.wolf <- beta.elk[4] + beta.elk[5]         # effect of elk.tmin1 on wolf.t And effect of elk.t on wolf.t
      sum.moose.wolf <- beta.moose[2] + beta.moose[3]   # effect of moose.tmin1 on wolf.t AND effect of moose.t on wolf.t
      sum.wtd.lion <- beta.wtd[2] + beta.wtd[3]         # effect of wtd.tmin1 on lion.t AND effect of wtd.t on lion.t
      sum.wtd.coy <- beta.wtd[6] + beta.wtd[7]          # effect of wtd.tmin1 on coy.t AND effect of wtd.t on coy.t
      
      #'  Derived parameters: Total Effects
      #'  ---------------------------------
      #'  Total effect of spp.A.tmin1 --> spp.B.t, summing over every direct path where A[t-1] --> B[t] and indirect paths
      #'  Total effect of bear on lion through effects on elk and wtd 
      #'  bear[t-1] --> elk[t] --> lion[t]  
      #'  bear[t-1] --> wtd.[t] --> lion[t]
      total.bear.lion <- beta.bear[2] * beta.elk[3] + beta.bear[3] * beta.wtd[3] 
      #'  Total effect of bear on wolf through effects on elk and wtd
      #'  bear[t-1] --> elk[t] --> wolf[t]  
      #'  bear[t-1] --> wtd[t] --> wolf[t]
      total.bear.wolf <- beta.bear[2] * beta.elk[5] + beta.bear[3] * beta.wtd[4] 
      #'  Total effect of elk on lion through its effect on itself 
      #'  elk[t-1] --> lion[t]  
      #'  elk[t-1] --> elk[t] --> lion[t]
      total.elk.lion <- beta.elk[2] + beta.elk[1] * beta.elk[3]         
      #'  Total effect of elk on wolf through its effect on itself and wolf effect on itself 
      #'  elk[t-1] --> wolf[t]  
      #'  elk[t-1] --> elk[t] --> wolf[t]  
      #'  elk[t-1] --> wolf[t-1] --> wolf[t]
      total.elk.wolf <- beta.elk[4] + beta.elk[1] * beta.elk[5] + beta.elk[5] * beta.wolf[1]  
      #'  Total effect of moose on wolf through its effect on itself and wolf effect on itself
      #'  same as above but swap moose for elk
      total.moose.wolf <- beta.moose[2] + beta.moose[1] * beta.moose[3] + beta.moose[3] * beta.wolf[1] 
      #'  Total effect of wtd on lion through effects on itself, bears, elk, and coy
      #'  wtd[t-1] --> lion[t]
      #'  wtd[t-1] --> wtd[t] --> lion[t]
      #'  wtd[t-1] --> bear[t-1] --> elk[t] --> lion[t]
      #'  wtd[t-1] --> bear[t-1] --> wtd[t] --> lion[t]
      #'  wtd[t-1] --> coy[t-1] --> wtd[t] --> lion[t]
      total.wtd.lion <- beta.wtd[2] + beta.wtd[1] * beta.wtd[3] + beta.wtd[5] * beta.bear[2] * beta.elk[3] + beta.wtd[5] * beta.bear[3] * beta.wtd[3] + beta.wtd[7] * beta.coy[2] * beta.wtd[3]
      #'  Total effect of wtd on wolf through effects on itself, bear, elk, and coy
      #'  wtd[t-1] --> wtd[t] --> wolf[t]
      #'  wtd[t-1] --> wolf[t-1] --> wolf[t]
      #'  wtd[t-1] --> bear[t-1] --> elk[t] --> wolf[t]
      #'  wtd[t-1] --> bear[t-1] --> wtd[t] --> wolf[t]
      #'  wtd[t-1] --> coy[t-1] --> wtd[t] --> wolf[t]
      total.wtd.wolf <- beta.wtd[1] * beta.wtd[4] + beta.wtd[4] * beta.wolf[1] + beta.wtd[5] * beta.bear[2] * beta.elk[5] + beta.wtd[5] * beta.bear[3] * beta.wtd[4] + beta.wtd[7] * beta.coy[2] * beta.wtd[4]
      #'  Total effect of wtd on bear through effects on itself, coy, and bear effect on itself
      #'  wtd[t-1] --> wtd[t] --> bear[t]
      #'  wtd[t-1] --> bear[t-1] --> bear[t]
      #'  wtd[t-1] --> bear[t-1] --> wtd[t] --> bear[t]
      #'  wtd[t-1] --> coy[t-1] --> wtd[t] --> bear[t]
      total.wtd.bear <- beta.wtd[1] * beta.wtd[5] + beta.wtd[5] * beta.bear[1] + beta.wtd[5] * beta.bear[3] * beta.wtd[5] + beta.wtd[7] * beta.coy[2] * beta.wtd[5]
      #'  Total effect of wtd on coy through effects on itself, bear, and coy effect on itself
      #'  wtd[t-1] --> coy[t]
      #'  wtd[t-1] --> wtd[t] --> coy[t]
      #'  wtd[t-1] --> coy[t-1] --> coy[t]
      #'  wtd[t-1] --> bear[t-1] --> wtd[t] --> coy[t]
      #'  wtd[t-1] --> coy[t-1] --> wtd[t] --> coy[t]
      total.wtd.coy <- beta.wtd[6] + beta.wtd[1] * beta.wtd[7] + beta.wtd[7] * beta.coy[1] + beta.wtd[5] * beta.bear[3] * beta.wtd[7] + beta.wtd[7] * beta.coy[2] * beta.wtd[7]
      #' #'  Total effect of bear on itself through AR1 and effect on wtd 
      #' #'  bear[t-1] --> bear[t]  
      #' #'  bear[t-1] --> wtd[t] --> bear[t]
      #' total.bear.self <- beta.bear[1] + beta.bear[3] * beta.wtd[5] 
      #' #'  Total effect of coy on itself through AR1 and effect on wtd  
      #' #'  coy[t-1] --> coy[t]  
      #' #'  coy[t-1] --> wtd[t] --> coy[t]
      #' total.coy.self <- beta.coy[1] + beta.coy[2] * beta.wtd[7]   
      #' #'  Total effect of wtd on itself through effects on bear and coy
      #' # total.wtd.wtd # ignoring 
      
      
      
      }")

