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
      
        #'  RN model posterior means (spp.t_hat) arise from latent true RDI (spp.t).
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
        
        #'  Ecological process model
        #'  All null models except the one conditional independence claim of interest
        #'  ---------------------------------------
        #'  NOTE ON USING THIS SCRIPT: For each claim, comment out the null model for
        #'  that species (e.g., mu.wolf.t[i] <- beta.int[2]) and all other claims
        #'  included in this script so only the claim of interest for a single 
        #'  species is run (e.g., mu.wolf.t[i] <- beta.int[2] + ... + beta.harvest[2] * deerHarv.tmin1[i])
        #'  Keep all other null models for remaining species and time lags.
        #'  ---------------------------------------
        
        wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
        # mu.wolf.t[i] <- beta.int[2] 
        #' #' bs_topdown_skinny[[6]] # deerHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.harvest[2] * deerHarv.tmin1[i]
        #' #' bs_topdown_skinny[[12]] # wtd.tmin1 and wolf.t are independent given wolfHarv.tmin and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        #' #' bs_topdown_skinny[[18]] # moose.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        #' #' bs_topdown_skinny[[24]] # elkHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1 
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.harvest[2] * elkHarv.tmin1[i]
        #' #' bs_topdown_skinny[[30]] # elk.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        #' #' bs_topdown_skinny[[47]] # coy.t and wolf.t are independent given coy.tmin1, wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.coy[2] * coy.t[i]
        #' #' bs_topdown_skinny[[52]] # bearHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i]
        #' #' bs_topdown_skinny[[55]] # bear.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        #' #' bs_topdown_skinny[[64]] # bear.t and wolf.t are independent given bearHarv.tmin1, bear.tmin1, wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i] + beta.bear[2] * bear.t[i]
        #' #' bs_topdown_skinny[[68]] # lionHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.harvest[2] * lionHarv.tmin1[i]
        #' #' bs_topdown_skinny[[70]] # lion.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.lion[1] + lion.tmin1[i]
        #' #' bs_topdown_skinny[[76]] # lion.t and wolf.t are independent given lionHarv.tmin1, lion.tmin1, wolfHarv.tmin1 and wolf.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.lion[1] + lion.tmin1[i] + beta.harvest[2] * lionHarv.tmin1[i] + beta.lion[2] * lion.t[i]
        #' #' bs_topdown_skinny[[82]] # wtd.t and wolf.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1, wolf.tmin1 and wolfHarv.tmin1
        #' mu.wolf.t[i] <- beta.int[2]  + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.lion[1] + lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.harvest[2] * deerHarv.tmin1[i] + beta.wtd[2] * wtd.t[i]
        #' #' bs_topdown_skinny[[84]] # moose.t and wolf.t are independent given moose.tmin1, wolf.tmin1 and wolfHarv.tmin1
        #' mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i] + beta.moose[2] * moose.t[i]
        #' bs_topdown_skinny[[85]] # elk.t and wolf.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1, wolf.tmin1 and wolfHarv.tmin1
        mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i] + beta.harvest[2] * elkHarv.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.elk[2] * elk.t[i]

        wolf.tmin1[i] ~ dnorm(mu.wolf.tmin1[i], tau.spp.tmin1[1])
        mu.wolf.tmin1[i] <- beta.int.tmin1[2]
        
        lion.t[i] ~ dnorm(mu.lion.t[i], tau.spp[2])  
        mu.lion.t[i] <- beta.int[1] 
        #' #' bs_topdown_skinny[[3]]: deerHarv.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.harvest[2] * deerHarv.tmin1[i]
        #' #' bs_topdown_skinny[[9]] # wtd.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        #' #' bs_topdown_skinny[[15]] # moose.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        #' #' bs_topdown_skinny[[21]] # elkHarv.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.harvest[2] * elkHarv.tmin1[i]
        #' #' bs_topdown_skinny[[27]] # elk.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1 
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        #' #' bs_topdown_skinny[[32]] # coy.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        #' #' bs_topdown_skinny[[41]] # coy.t and lion.t are independent given lionHarv.tmin1, lion.tmin1 and coy.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.coy[2] * coy.t[i]
        #' #' bs_topdown_skinny[[48]] # bearHarv.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i]
        #' #' bs_topdown_skinny[[53]] # bear.tmin1 and lion.t are independent given lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        #' #' bs_topdown_skinny[[58]] # bear.t and lion.t are independent given bearHarv.tmin1, bear.tmin1, lionHarv.tmin1 and lion.tmin1
        #' mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i] + beta.bear[2] * bear.t[i]
        
        lion.tmin1[i] ~ dnorm(mu.lion.tmin1[i], tau.spp.tmin1[2])
        mu.lion.tmin1[i] <- beta.int.tmin1[1] 
        
        bear.t[i] ~ dnorm(mu.bear.t[i], tau.spp[3]) 
        mu.bear.t[i] <- beta.int[3] 
        #' #'  bs_topdown_skinny[[2]]: deerHarv.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
        #' mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.harvest[2] * deerHarv.tmin1[i]
        #' #' bs_topdown_skinny[[8]] # wtd.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
        #' mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        #' #' bs_topdown_skinny[[14]] # moose.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
        #' mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        #' #' bs_topdown_skinny[[20]] # elkHarv.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
        #' mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.harvest[2] * elkHarv.tmin1[i]
        #' #' bs_topdown_skinny[[26]] # elk.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1 
        #' mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        #' #' bs_topdown_skinny[[31]] # coy.tmin1 and bear.t are independent given bearHarv.tmin1 and bear.tmin1
        #' mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        #' #' bs_topdown_skinny[[38]] # coy.t and bear.t are independent given coy.tmin1, bearHarv.tmin1 and bear.tmin1
        #' mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.coy[2] * coy.t[i]
        
        bear.tmin1[i] ~ dnorm(mu.bear.tmin1[i], tau.spp.tmin1[3]) 
        mu.bear.tmin1[i] <- beta.int.tmin1[3] 

        coy.t[i] ~ dnorm(mu.coy.t[i], tau.spp[4])
        mu.coy.t[i] <- beta.int[4] 
        #' #'  bs_topdown_skinny[[1]]: deerHarv.tmin1 and coy.t are independent given coy.tmin1
        #' mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i]
        #' #' bs_topdown_skinny[[7]] # wtd.tmin1 and coy.t are independent given coy.tmin1
        #' mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        #' #' bs_topdown_skinny[[13]] # moose.tmin1 and coy.t are independent given coy.tmin1
        #' mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        #' #' bs_topdown_skinny[[19]] # elkHarv.tmin1 and coy.t are independent given coy.tmin1
        #' mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i]
        #' #' bs_topdown_skinny[[25]] # elk.tmin1 and coy.t are independent given coy.tmin1 
        #' mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.elk[1] * elk.tmin1[i]
  
        coy.tmin1[i] ~ dnorm(mu.coy.tmin1[i], tau.spp.tmin1[4])
        mu.coy.tmin1[i] <- beta.int.tmin1[4] 
        
        elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
        mu.elk.t[i] <- beta.int[5] 
        #'  # bs_topdown_skinny[[5]] # deerHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.harvest[2] * deerHarv.tmin1[i] 
        #' #' bs_topdown_skinny[[11]] # wtd.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        #' #' bs_topdown_skinny[[17]] # moose.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        #' #' bs_topdown_skinny[[34]] # coy.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        #' #' bs_topdown_skinny[[46]] # coy.t and elk.t are independent given coy.tmin1, elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.coy[2] * coy.t[i]
        #' #' bs_topdown_skinny[[51]] # bearHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i]
        #' bs_topdown_skinny[[63]] # bear.t and elk.t are independent given bearHarv.tmin1, bear.tmin1, elkHarv.tmin1, elk.tmin1, lion.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i] + beta.bear[2] * bear.t[i]
        #' #' bs_topdown_skinny[[67]] # lionHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.harvest[2] * lionHarv.tmin1[i]
        #' #' bs_topdown_skinny[[75]] # lion.t and elk.t are independent given lionHarv.tmin1, lion.tmin1, elkHarv.tmin1, elk.tmin1, bear.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.harvest[2] * lionHarv.tmin1[i] + beta.lion[2] * lion.t[i]
        #' #' bs_topdown_skinny[[79]] # wolfHarv.tmin1 and elk.t are independent given elkHarv.tmin1, elk.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.harvest[2] * wolfHarv.tmin1[i]
        #' #' bs_topdown_skinny[[81]] # wtd.t and elk.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1, elk.tmin1, elkHarv.tmin1, and wolf.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.harvest[2] * deerHarv.tmin1[i] + beta.wtd[2] * wtd.t[i]
        #' #' bs_topdown_skinny[[83]] # moose.t and elk.t are independent given moose.tmin1, wolf.tmin1, elkHarv.tmin1, elk.tmin1, bear.tmin1 and lion.tmin1
        #' mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.elk[1] * elk.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i] + beta.moose[2] * moose.t[i]
        
        elk.tmin1[i] ~ dnorm(mu.elk.tmin1[i], tau.spp.tmin1[5])
        mu.elk.tmin1[i] <- beta.int.tmin1[5]
        
        moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
        mu.moose.t[i] <- beta.int[6] 
        #' #' bs_topdown_skinny[[4]] # deerHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i]
        #' #'  bs_topdown_skinny[[10]] # wtd.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        #' #' bs_topdown_skinny[[23]] # elkHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1 
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * elkHarv.tmin1[i]
        #' #' bs_topdown_skinny[[29]] # elk.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        #' #' bs_topdown_skinny[[33]] # coy.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        #' #' bs_topdown_skinny[[45]] # coy.t and moose.t are independent given coy.tmin1, moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.coy[2] * coy.t[i]
        #' #' bs_topdown_skinny[[50]] # bearHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i]
        #' #' bs_topdown_skinny[[54]] # bear.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        #' #' bs_topdown_skinny[[62]] # bear.t and moose.t are independent given bearHarv.tmin1, bear.tmin1, moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.harvest[1] * bear.tmin1[i] + beta.bear[2] * bear.t[i]
        #' #' bs_topdown_skinny[[66]] # lionHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i]
        #' #'  bs_topdown_skinny[[69]] # lion.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i]
        #' #' bs_topdown_skinny[[74]] # lion.t and moose.t are independent given lionHarv.tmin1, lion.tmin1, moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.lion[2] * lion.t[i]
        #' #' bs_topdown_skinny[[78]] # wolfHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i]
        #' #' bs_topdown_skinny[[80]] # wtd.t and moose.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1, moose.tmin1 and wolf.tmin1
        #' mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.wtd[2] * wtd.t[i]

        moose.tmin1[i] ~ dnorm(mu.moose.tmin1[i], tau.spp.tmin1[6])
        mu.moose.tmin1[i] <- beta.int.tmin1[6]
        
        wtd.t[i] ~ dnorm(mu.wtd.t[i], tau.spp[7])
        mu.wtd.t[i] <- beta.int[7] 
        #' #' bs_topdown_skinny[[16]] # moose.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i] 
        #' #' bs_topdown_skinny[[22]] # elkHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1 
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.harvest[2] * elkHarv.tmin1[i] 
        #' #' bs_topdown_skinny[[28]] # elk.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        #' #' bs_topdown_skinny[[44]] # coy.t and wtd.t are independent given coy.tmin1, deerHarv.tmin1, wtd.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.coy[2] * coy.t[i]
        #' #' bs_topdown_skinny[[49]] # bearHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i]
        #' #' bs_topdown_skinny[[61]] # bear.t and wtd.t are independent given bearHarv.tmin1, bear.tmin1, deerHarv.tmin1, wtd.tmin1, coy.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i] + beta.bear[2] * bear.t[i] 
        #' #' bs_topdown_skinny[[65]] # lionHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.harvest[2] * lionHarv.tmin1[i]
        #' #' bs_topdown_skinny[[73]] # lion.t and wtd.t are independent given lionHarv.tmin1, lion.tmin1, deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.lion[2] * lion.t[i]
        #' #' bs_topdown_skinny[[77]] # wolfHarv.tmin1 and wtd.t are independent given deerHarv.tmin1, wtd.tmin1, coy.tmin1, bear.tmin1, lion.tmin1 and wolf.tmin1
        #' mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i] + beta.coy[1] * coy.tmin1[i] + beta.harvest[1] * deerHarv.tmin1[i] + beta.harvest[2] * wolfHarv.tmin1[i]
        
        wtd.tmin1[i] ~ dnorm(mu.wtd.tmin1[i], tau.spp.tmin1[7])
        mu.wtd.tmin1[i] <- beta.int.tmin1[7] 
        
      }
      
      }")