#'  --------------------------------------------------------
#'  JAGS model: top-down SEM d-separation
#'  "beta.harvest[6]" and "mu.moose.t" are independent given "beta.wolf[3]" and "beta.moose[1]"
#'  --------------------------------------------------------
cat(file = './Outputs/SEM/JAGS_out/d_Sep/JAGS_SEM_topdown_dsep.txt', "
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
      for(w in 1:nWolf) {
        beta.wolf[w] ~ dnorm(0, 0.1)  
      }
      for(l in 1:nLion) {
        beta.lion[l] ~ dnorm(0, 1)  
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

      #'  Priors for anthropogenic effects
      for(h in 1:nharvest) {
        beta.harvest[h] ~ dnorm(0, 0.01)
      }
      #' Half-Cauchy prior for latent variable SD
      scale <- 1  
      for(k in 1:nSpp) {
        aux[k] ~ dgamma(0.5, 0.5)
        sigma.spp[k] ~ dnorm(0, 1 / (pow(scale, 2) * aux[k])) T(0,)
        tau.spp[k] <- 1 / pow(sigma.spp[k], 2)
        sigma.spp.tmin1[k] ~ dnorm(0, 1 / (pow(scale, 2) * aux[k])) T(0,)
        tau.spp.tmin1[k] <- 1 / pow(sigma.spp.tmin1[k], 2)
      }
      
      #'  Likelihood
      #'  ----------
      for(i in 1:nCluster) {
        
        #'  bs_topdown_skinny[[8]]: beta.coy[2] and mu.elk.t are independent given beta.harvest[4], beta.bear[2], beta.lion[2], beta.wolf[2], and beta.elk[1]
        #'  Noisy observation data arises from latent RDI
        wolf.tmin1_hat[i] ~ dnorm(wolf.tmin1[i], wolf.tmin1.tau_hat[i])
        lion.tmin1_hat[i] ~ dnorm(lion.tmin1[i], lion.tmin1.tau_hat[i])
        bear.tmin1_hat[i] ~ dnorm(bear.tmin1[i], bear.tmin1.tau_hat[i])
        coy.tmin1_hat[i] ~ dnorm(coy.tmin1[i], coy.tmin1.tau_hat[i])
        elk.t_hat[i] ~ dnorm(elk.t[i], elk.t.tau_hat[i])
        elk.tmin1_hat[i] ~ dnorm(elk.tmin1[i], elk.tmin1.tau_hat[i])
        
        #'  Precision based on SD of noisy observations
        wolf.tmin1.tau_hat[i] <- 1 / pow(wolf.tmin1.sigma_hat[i], 2)
        lion.tmin1.tau_hat[i] <- 1 / pow(lion.tmin1.sigma_hat[i], 2)
        bear.tmin1.tau_hat[i] <- 1 / pow(bear.tmin1.sigma_hat[i], 2)
        coy.tmin1.tau_hat[i] <- 1 / pow(coy.tmin1.sigma_hat[i], 2)
        elk.t.tau_hat[i] <- 1 / pow(elk.t.sigma_hat[i], 2)
        elk.tmin1.tau_hat[i] <- 1 / pow(elk.tmin1.sigma_hat[i], 2)
        
        #'  Latent RDI
        wolf.tmin1[i] ~ dnorm(mu.wolf.tmin1[i], tau.spp.tmin1[1])
        mu.wolf.tmin1[i] <- beta.int.tmin1[2]
        lion.tmin1[i] ~ dnorm(mu.lion.tmin1[i], tau.spp.tmin1[2])
        mu.lion.tmin1[i] <- beta.int.tmin1[1]
        bear.tmin1[i] ~ dnorm(mu.bear.tmin1[i], tau.spp.tmin1[3]) 
        mu.bear.tmin1[i] <- beta.int.tmin1[3]
        coy.tmin1[i] ~ dnorm(mu.coy.tmin1[i], tau.spp.tmin1[4])
        mu.coy.tmin1[i] <- beta.int.tmin1[4]
        elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
        mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[2] * wolf.tmin1[i] + beta.lion[2] * lion.tmin1[i] + beta.bear[2] * bear.tmin1[i] + beta.harvest[4] * elkHarv.tmin1[i] + beta.coy[2] * coy.tmin1[i]        elk.tmin1[i] ~ dnorm(mu.elk.tmin1[i], tau.spp.tmin1[5])
        mu.elk.tmin1[i] <- beta.int.tmin1[5]
        
        #' #'  bs_topdown[[7]]: beta.coy[2] and mu.moose.t are independent given beta.wolf[3] and beta.moose[1]
        #' #'  Noisy observation data arises from latent RDI
        #' wolf.tmin1_hat[i] ~ dnorm(wolf.tmin1[i], wolf.tmin1.tau_hat[i])
        #' coy.tmin1_hat[i] ~ dnorm(coy.tmin1[i], coy.tmin1.tau_hat[i])
        #' moose.t_hat[i] ~ dnorm(moose.t[i], moose.t.tau_hat[i])
        #' moose.tmin1_hat[i] ~ dnorm(moose.tmin1[i], moose.tmin1.tau_hat[i])
        #' 
        #' #'  Precision based on SD of noisy observations
        #' wolf.tmin1.tau_hat[i] <- 1 / pow(wolf.tmin1.sigma_hat[i], 2)
        #' coy.tmin1.tau_hat[i] <- 1 / pow(coy.tmin1.sigma_hat[i], 2)
        #' moose.t.tau_hat[i] <- 1 / pow(moose.t.sigma_hat[i], 2)
        #' moose.tmin1.tau_hat[i] <- 1 / pow(moose.tmin1.sigma_hat[i], 2)
        #' 
        #' #'  Latent RDI
        #' wolf.tmin1[i] ~ dnorm(mu.wolf.tmin1[i], tau.spp.tmin1[1])
        #' mu.wolf.tmin1[i] <- beta.int.tmin1[2]
        #' coy.tmin1[i] ~ dnorm(mu.coy.tmin1[i], tau.spp.tmin1[4])
        #' mu.coy.tmin1[i] <- beta.int.tmin1[4]
        #' moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[3] * wolf.tmin1[i] + beta.coy[2] * coy.tmin1[i] 
        #' moose.tmin1[i] ~ dnorm(mu.moose.tmin1[i], tau.spp.tmin1[6])
        #' mu.moose.tmin1[i] <- beta.int.tmin1[6]
        
      }
      
      
      
      
      }")