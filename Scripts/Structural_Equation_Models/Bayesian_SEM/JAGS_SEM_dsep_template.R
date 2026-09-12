  #'  --------------------------------------
  #'  JAGS template for d-separation test
  #'  August 2026
  #'  --------------------------------------
  #'  Keeping all priors consistent with original SEMs. Template arises in the 
  #'  likelihood, where the regression for each species and time step is coded
  #'  using a placeholder (percent sign s) that is updated with each iteration for each
  #'  d-separation test.
  #'  --------------------------------------

  model_template <- "
      model {
      
      #'  Define priors
      #'  -------------
      #'  Priors for intercepts
      #'  Use more informed prior for lion intercept
      beta.int[1] ~ dnorm(0, 1) # poor convergence with weaker priors
      beta.int.tmin1[1] ~ dnorm(0, 0.01)
      
      #'  Intercept priors for all other species 
      for(k in 2:nSpp) {
        beta.int[k] ~ dnorm(0, 0.01) 
        beta.int.tmin1[k] ~ dnorm(0, 0.01)
      }
      
      #'  Priors for species lag effects
      #'  As a reminder: precision = 0.01 --> sqrt(0.01^-1) --> SD = 10
      for(l in 1:nLion) {
        beta.lion[l] ~ dnorm(0, 1)  
      }
      for(w in 1:nWolf) {
        beta.wolf[w] ~ dnorm(0, 0.1)  
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
      for(f in 1:nforest) {
        beta.forest[f] ~ dnorm(0, 0.01)
      }
      for(s in 1:nWSI) {
        beta.wsi[s] ~ dnorm(0, 0.01)
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
          mu.lion[i,y] <- %s
          
          wolf.latent[i,y] ~ dnorm(mu.wolf[i,y], tau.spp[2])
          mu.wolf[i,y] <- %s
          
          bear.latent[i,y] ~ dnorm(mu.bear[i,y], tau.spp[3])
          mu.bear[i,y] <- %s
          
          coy.latent[i,y] ~ dnorm(mu.coy[i,y], tau.spp[4])
          mu.coy[i,y] <- %s
          
          elk.latent[i,y] ~ dnorm(mu.elk[i,y], tau.spp[5])
          mu.elk[i,y] <- %s
          
          moose.latent[i,y] ~ dnorm(mu.moose[i,y], tau.spp[6])
          mu.moose[i,y] <- %s
          
          wtd.latent[i,y] ~ dnorm(mu.wtd[i,y], tau.spp[7])
          mu.wtd[i,y] <- %s
          
        }
      }
      
    }"
