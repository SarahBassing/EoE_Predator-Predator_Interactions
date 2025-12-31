  #'  JAGS model (top-down, interference model with cross lags)
  #'  Species indices: 1 = wolf, 2 = mountain lion, 3 = black bear, 4 = coyote,
  #'  5 = elk, 6 = moose, 7 = white-tailed deer
  cat(file = './Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_crosslag.txt', "
        model{
      
      #'  Define priors
      #'  -------------
      #'  Intercept priors for each species-specific regression (2 priors per species 
      #'  b/c models vary between timestep 1 and later timesteps)
      
      for(k in 1:nSpp) {
        for(mod in 1:2) {
          beta.int[k,mod] ~ dnorm(0, 0.01)  # look more into precision value... I think FSP review commented on this?
        }
      }
      
      #'  SD prior for each timestep and species-specific regression
      for(k in 1:nSpp) {
        for(t in 1:nTimestep) {
          sigma.spp[k,t] ~ dunif(0, 10)
          tau.spp[k,t] <- 1 / pow(sigma.spp[k,t], 2)
        }
      }
      
      
      #'  Priors for species lag effects
      
      for(w in 1:nWolf) {
        beta.wolf[w] ~ dnorm(0, 0.01)  # look more into precision value... I think FSP review commented on this?
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
      for(d in 1:nDeer) {
        beta.wtd[d] ~ dnorm(0, 0.01)
      }

      #'  Priors for anthropogenic & landscape effects
      for(h in 1:nharvest) {
        beta.harvest[h] ~ dnorm(0, 0.01)
      }
      # for(f in 1:nforest) {
      #   forest[f] ~ dnorm(0, 0.01)
      # }
      

      #'  SD priors for among-cluster random intercept term for each species -- DO I WANT A RANDOM EFFECT FOR EACH SPP AND CLUSTER OR JUST ONE PER CLUSTER??? Mainly accounting for repeat measure of clusters across years, right?
      for(k in 1:nSpp) {
        for(cl in 1:nCluster) {
          sigma.cluster[k,cl] ~ dunif(0, 10) #dnorm(0, 0.01) T(0,)            
          tau.cluster[k,cl] <- 1 / pow(sigma.cluster[k,cl], 2)    
        }
      }
      
      
      
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
          
          
        #'  If t > 1
        for(t in 2:nTimestep) {
          wolf[i,t] ~ dnorm(mu.wolf[i,t], tau.spp[1,t])
          mu.wolf[i,t] <- beta.int[1,2] + beta.wolf[1] * wolf[i,t-1] + beta.harvest[1] * harv[i,t-1] + tau.cluster[1,i]

          lion[i,t] ~ dnorm(mu.lion[i,t], tau.spp[2,t])
          mu.lion[i,t] <- beta.int[2,2] + beta.lion[1] * lion[i,t-1] + beta.wolf[2] * wolf[i,t-1] + beta.bear[2] * bear[i,t-1] + tau.cluster[2,i]

          bear[i,t] ~ dnorm(mu.bear[i,t], tau.spp[3,t])
          mu.bear[i,t] <- beta.int[3,2] + beta.bear[1] * bear[i,t-1] + beta.wolf[3] * wolf[i,t-1]
        
        }
        
      }
      
      
      #'  Derived parameters
      #'  ------------------
      
      
      }")