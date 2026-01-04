  #'  JAGS model (top-down, interference model)
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
        wolf.t_1[w] ~ dnorm(0, 0.01)  # precision = 0.01 --> sqrt(0.01^-1) --> SD = 10
      }
      for(l in 1:nLion) {
        lion.t_1[l] ~ dnorm(0, 0.01)
      }
      for(b in 1:nBear) {
        bear.t_1[b] ~ dnorm(0, 0.01)
      }
      for(c in 1:nCoy) {
        coy.t_1[c] ~ dnorm(0, 0.01)
      }
      for(e in 1:nElk) {
        elk.t_1[e] ~ dnorm(0, 0.01)
      }
      for(m in 1:nMoose) {
        moose.t_1[m] ~ dnorm(0, 0.01)
      }
      for(d in 1:nWtd) {
        wtd.t_1[d] ~ dnorm(0, 0.01)
      }

      #'  Priors for anthropogenic & landscape effects
      for(h in 1:nharvest) {
        harvest.t_1[h] ~ dnorm(0, 0.01)
      }
      for(f in 1:nforest) {
        forest.t_1[f] ~ dnorm(0, 0.01)
      }
      
      #'  SD prior for each regression
      for(k in 1:nSpp) {
        sigma.spp[k] ~ dunif(0, 10)
        tau.spp[k] <- 1 / pow(sigma.spp[k], 2)
      }
      
      #'  SD prior for among-cluster random intercept term for each species 
      for(k in 1:nSpp) {
        for(cl in 1:nCluster) {
          sigma.cluster[k,cl] ~ dnorm(0, 0.01) T(0,)      #dunif(0, 10)      
          tau.cluster[k,cl] <- 1 / pow(sigma.cluster[k,cl], 2)    
        }
      }
      #'  dnorm(0, 0.01) T(0,) some seem to hitting boundary at 0 and a few not converging well 
      #'  dunif(0, 10) most are just returning the prior, a few have poor convergence
      #'  dnorm seems to be the better prior if sticking with a species - cluster random effect
      
      #'  Likelihood
      #'  ----------
      for(i in 1:nCluster) {
        wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
        mu.wolf.t[i] <- beta.int[1] + wolf.t_1[1] * wolf.tmin1[i] + harvest.t_1[1] * harvest.tmin1[i] + tau.cluster[1,i] 
        
        lion.t[i] ~ dnorm(mu.lion.t[i], tau.spp[2])
        mu.lion.t[i] <- beta.int[2] + lion.t_1[1] * lion.tmin1[i] + wolf.t_1[2] * wolf.tmin1[i] + bear.t_1[2] * bear.tmin1[i] + tau.cluster[2,i]

        bear.t[i] ~ dnorm(mu.bear.t[i], tau.spp[3])
        mu.bear.t[i] <- beta.int[3] + bear.t_1[1] * bear.tmin1[i] + wolf.t_1[3] * wolf.tmin1[i] + tau.cluster[3,i]

        coy.t[i] ~ dnorm(mu.coy.t[i], tau.spp[4])
        mu.coy.t[i] <- beta.int[4] + coy.t_1[1] * coy.tmin1[i] + wolf.t_1[4] * wolf.tmin1[i] + lion.t_1[2] * lion.tmin1[i] + bear.t_1[3] * bear.tmin1[i] + tau.cluster[4,i]

        elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
        mu.elk.t[i] <- beta.int[5] + elk.t_1[1] * elk.tmin1[i] + wolf.t_1[5] * wolf.tmin1[i] + lion.t_1[3] * lion.tmin1[i] + bear.t_1[4] * bear.tmin1[i] + tau.cluster[5,i]

        moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
        mu.moose.t[i] <- beta.int[5] + moose.t_1[1] * moose.tmin1[i] + wolf.t_1[6] * wolf.tmin1[i] + tau.cluster[6,i]

        wtd.t[i] ~ dnorm(mu.wtd.t[i], tau.spp[7])
        mu.wtd.t[i] <- beta.int[7] + wtd.t_1[1] * wtd.tmin1[i] + wolf.t_1[7] * wolf.tmin1[i] + lion.t_1[4] * lion.tmin1[i] + bear.t_1[5] * bear.tmin1[i] + coy.t_1[2] * coy.tmin1[i] + tau.cluster[7,i]

      }
      
      
      #'  Derived parameters
      #'  ------------------
      
      
      }")