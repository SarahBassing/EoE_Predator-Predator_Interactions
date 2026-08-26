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
      
      #' #' Half-Cauchy prior for latent variable SD
      #' scale <- 1  
      #' for(k in 1:nSpp) {
      #'   aux[k] ~ dgamma(0.5, 0.5)
      #'   sigma.spp[k] ~ dnorm(0, 1 / (pow(scale, 2) * aux[k])) T(0,)
      #'   tau.spp[k] <- 1 / pow(sigma.spp[k], 2)
      #'   sigma.spp.tmin1[k] ~ dnorm(0, 1 / (pow(scale, 2) * aux[k])) T(0,)
      #'   tau.spp.tmin1[k] <- 1 / pow(sigma.spp.tmin1[k], 2)
      #' }
      
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
  
  
  #' #'  Likelihood
  #' #'  ----------
  #' for(i in 1:nCluster) {
  #'   #'  RN model posterior means (spp.t_hat) arise from latent true RDI (spp.t).
  #'   wolf.t_hat[i] ~ dnorm(wolf.t[i], wolf.t.tau_hat[i])
  #'   wolf.tmin1_hat[i] ~ dnorm(wolf.tmin1[i], wolf.tmin1.tau_hat[i])
  #'   lion.t_hat[i] ~ dnorm(lion.t[i], lion.t.tau_hat[i])
  #'   lion.tmin1_hat[i] ~ dnorm(lion.tmin1[i], lion.tmin1.tau_hat[i])
  #'   bear.t_hat[i] ~ dnorm(bear.t[i], bear.t.tau_hat[i])
  #'   bear.tmin1_hat[i] ~ dnorm(bear.tmin1[i], bear.tmin1.tau_hat[i])
  #'   coy.t_hat[i] ~ dnorm(coy.t[i], coy.t.tau_hat[i])
  #'   coy.tmin1_hat[i] ~ dnorm(coy.tmin1[i], coy.tmin1.tau_hat[i])
  #'   elk.t_hat[i] ~ dnorm(elk.t[i], elk.t.tau_hat[i])
  #'   elk.tmin1_hat[i] ~ dnorm(elk.tmin1[i], elk.tmin1.tau_hat[i])
  #'   moose.t_hat[i] ~ dnorm(moose.t[i], moose.t.tau_hat[i])
  #'   moose.tmin1_hat[i] ~ dnorm(moose.tmin1[i], moose.tmin1.tau_hat[i])
  #'   wtd.t_hat[i] ~ dnorm(wtd.t[i], wtd.t.tau_hat[i])
  #'   wtd.tmin1_hat[i] ~ dnorm(wtd.tmin1[i], wtd.tmin1.tau_hat[i])
  #'   
  #'   #'  RN model posterior SD (spp.t.sigma_hat) used to calculate spp.t.tau_hat
  #'   wolf.t.tau_hat[i] <- 1 / pow(wolf.t.sigma_hat[i], 2)
  #'   wolf.tmin1.tau_hat[i] <- 1 / pow(wolf.tmin1.sigma_hat[i], 2)
  #'   lion.t.tau_hat[i] <- 1 / pow(lion.t.sigma_hat[i], 2)
  #'   lion.tmin1.tau_hat[i] <- 1 / pow(lion.tmin1.sigma_hat[i], 2)
  #'   bear.t.tau_hat[i] <- 1 / pow(bear.t.sigma_hat[i], 2)
  #'   bear.tmin1.tau_hat[i] <- 1 / pow(bear.tmin1.sigma_hat[i], 2)
  #'   coy.t.tau_hat[i] <- 1 / pow(coy.t.sigma_hat[i], 2)
  #'   coy.tmin1.tau_hat[i] <- 1 / pow(coy.tmin1.sigma_hat[i], 2)
  #'   elk.t.tau_hat[i] <- 1 / pow(elk.t.sigma_hat[i], 2)
  #'   elk.tmin1.tau_hat[i] <- 1 / pow(elk.tmin1.sigma_hat[i], 2)
  #'   moose.t.tau_hat[i] <- 1 / pow(moose.t.sigma_hat[i], 2)
  #'   moose.tmin1.tau_hat[i] <- 1 / pow(moose.tmin1.sigma_hat[i], 2)
  #'   wtd.t.tau_hat[i] <- 1 / pow(wtd.t.sigma_hat[i], 2)
  #'   wtd.tmin1.tau_hat[i] <- 1 / pow(wtd.tmin1.sigma_hat[i], 2)
  #'   
  #'   #'  Ecological process model
  #'   #'  Regression #1
  #'   lion.t[i] ~ dnorm(mu.lion.t[i], tau.spp[2])  
  #'   mu.lion.t[i] <- %s  
  #'   
  #'   #'  Regression #2
  #'   lion.tmin1[i] ~ dnorm(mu.lion.tmin1[i], tau.spp.tmin1[2])
  #'   mu.lion.tmin1[i] <- %s 
  #'   
  #'   #'  Regression #3
  #'   wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
  #'   mu.wolf.t[i] <- %s  
  #'   
  #'   #'  Regression #4
  #'   wolf.tmin1[i] ~ dnorm(mu.wolf.tmin1[i], tau.spp.tmin1[1])
  #'   mu.wolf.tmin1[i] <- %s 
  #'   
  #'   #'  Regression #5
  #'   bear.t[i] ~ dnorm(mu.bear.t[i], tau.spp[3]) 
  #'   mu.bear.t[i] <- %s  
  #'   
  #'   #'  Regression #6
  #'   bear.tmin1[i] ~ dnorm(mu.bear.tmin1[i], tau.spp.tmin1[3]) 
  #'   mu.bear.tmin1[i] <- %s 
  #'   
  #'   #'  Regression #7
  #'   coy.t[i] ~ dnorm(mu.coy.t[i], tau.spp[4])
  #'   mu.coy.t[i] <- %s 
  #'   
  #'   #'  Regression #8
  #'   coy.tmin1[i] ~ dnorm(mu.coy.tmin1[i], tau.spp.tmin1[4])
  #'   mu.coy.tmin1[i] <- %s  
  #'   
  #'   #'  Regression #9
  #'   elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
  #'   mu.elk.t[i] <- %s 
  #'   
  #'   #'  Regression #10
  #'   elk.tmin1[i] ~ dnorm(mu.elk.tmin1[i], tau.spp.tmin1[5])
  #'   mu.elk.tmin1[i] <- %s 
  #'   
  #'   #'  Regression #11
  #'   moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
  #'   mu.moose.t[i] <- %s  
  #'   
  #'   #'  Regression #12
  #'   moose.tmin1[i] ~ dnorm(mu.moose.tmin1[i], tau.spp.tmin1[6])
  #'   mu.moose.tmin1[i] <- %s 
  #'   
  #'   #'  Regression #13
  #'   wtd.t[i] ~ dnorm(mu.wtd.t[i], tau.spp[7])
  #'   mu.wtd.t[i] <- %s    
  #'   
  #'   #'  Regression #14
  #'   wtd.tmin1[i] ~ dnorm(mu.wtd.tmin1[i], tau.spp.tmin1[7])
  #'   mu.wtd.tmin1[i] <- %s 
  #' }