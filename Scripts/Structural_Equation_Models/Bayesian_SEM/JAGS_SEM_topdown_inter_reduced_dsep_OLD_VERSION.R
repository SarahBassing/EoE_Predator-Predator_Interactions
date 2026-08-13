  #'  --------------------------------------------------------
  #'  JAGS model: top-down interference SEM d-separation
  #'  e.g., wtd.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
  #'  --------------------------------------------------------
  cat(file = './Outputs/SEM/JAGS_out/d_Sep/JAGS_SEM_topdown_inter_reduced_dsep.txt', "
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

      #'  Priors for anthropogenic & landscape effects
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
        wolf.t[i] ~ dnorm(mu.wolf.t[i], tau.spp[1])
        mu.wolf.t[i] <- beta.int[2] 
        # bs_topdown_inter_skinny[[6]] # wtd.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        # bs_topdown_inter_skinny[[12]] # moose.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.moose[1] * wolf.tmin1[i]
        # bs_topdown_inter_skinny[[18]] # elk.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        # bs_topdown_inter_skinny[[24]] # coy.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        # bs_topdown_inter_skinny[[30]] # bearHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i]
        # bs_topdown_inter_skinny[[36]] # bear.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        # bs_topdown_inter_skinny[[42]] # lionHarv.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.harvest[2] * lionHarv.tmin1[i]
        # bs_topdown_inter_skinny[[45]] # lion.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.lion[1] * lion.tmin1[i]
        # bs_topdown_inter_skinny[[53]] # wtd.t and wolf.t are independent given wtd.tmin1, lion.tmin1, wolfHarv.tmin1 and wolf.tmin1
        # mu.wolf.t[i] <- beta.int[2] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] + wolfHarv.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.wtd[2] * wtd.t[i]
        
        wolf.tmin1[i] ~ dnorm(mu.wolf.tmin1[i], tau.spp.tmin1[1])
        mu.wolf.tmin1[i] <- beta.int.tmin1[2]
        
        lion.t[i] ~ dnorm(mu.lion.t[i], tau.spp[2])  
        mu.lion.t[i] <- beta.int[1] 
        # bs_topdown_inter_skinny[[5]] # wtd.tmin1 and lion.t are independent given lionHarv.tmin1, lion.tmin1 and wolf.tmin1
        # mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        # bs_topdown_inter_skinny[[11]] # moose.tmin1 and lion.t are independent given lionHarv.tmin1, lion.tmin1 and wolf.tmin1
        # mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        # bs_topdown_inter_skinny[[17]] # elk.tmin1 and lion.t are independent given lionHarv.tmin1, lion.tmin1 and wolf.tmin1
        # mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        # bs_topdown_inter_skinny[[23]] # coy.tmin1 and lion.t are independent given lionHarv.tmin1, lion.tmin1 and wolf.tmin1
        # mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        # bs_topdown_inter_skinny[[29]] # bearHarv.tmin1, lion.t are independent given lionHarv.tmin1, lion.tmin1 and wolf.tmin1
        # mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.harvest[2] * bearHarv.tmin1[i]
        # bs_topdown_inter_skinny[[35]] # bear.tmin1 and lion.t are independent given lionHarv.tmin1, lion.tmin1 and wolf.tmin1
        # mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        # bs_topdown_inter_skinny[[52]] # wtd.t and lion.t are independent given wtd.tmin1, lion.tmin1, lionHarv.tmin1 and wolf.tmin1
        # mu.lion.t[i] <- beta.int[1] + beta.lion[1] * lion.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.wtd[2] * wtd.t[i]
        
        lion.tmin1[i] ~ dnorm(mu.lion.tmin1[i], tau.spp.tmin1[2])
        mu.lion.tmin1[i] <- beta.int.tmin1[1] 
        
        bear.t[i] ~ dnorm(mu.bear.t[i], tau.spp[3]) 
        mu.bear.t[i] <- beta.int[3] 
        # bs_topdown_inter_skinny[[4]] # wtd.tmin1 and bear.t are independent given bearHarv.tmin1, bear.tmin1 and wolf.tmin1
        # mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        # bs_topdown_inter_skinny[[10]] # moose.tmin1 and bear.t are independent given bearHarv.tmin1, bear.tmin1 and wolf.tmin1
        # mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        # bs_topdown_inter_skinny[[16]] # elk.tmin1 and bear.t are independent given bearHarv.tmin1, bear.tmin1 and wolf.tmin1
        # mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        # bs_topdown_inter_skinny[[22]] # coy.tmin1 and bear.t are independent given bearHarv.tmin1, bear.tmin1 and wolf.tmin1
        # mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        # bs_topdown_inter_skinny[[41]] # lionHarv.tmin1 and bear.t are independent given bearHarv.tmin1, bear.tmin1 and wolf.tmin1
        # mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.harvest[2] * lionHarv.tmin1[i]
        # bs_topdown_inter_skinny[[44]] # lion.tmin1 and bear.t are independent given bearHarv.tmin1, bear.tmin1 and wolf.tmin1
        # mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.lion[1] * lion.tmin1[i]
        # bs_topdown_inter_skinny[[51]] # wtd.t and bear.t are independent given wtd.tmin1, lion.tmin1, bearHarv.tmin1, bear.tmin1 and wolf.tmin1
        # mu.bear.t[i] <- beta.int[3] + beta.bear[1] * bear.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.wtd[2] * wtd.t[i]
        
        bear.tmin1[i] ~ dnorm(mu.bear.tmin1[i], tau.spp.tmin1[3]) 
        mu.bear.tmin1[i] <- beta.int.tmin1[3] 

        coy.t[i] ~ dnorm(mu.coy.t[i], tau.spp[4])
        mu.coy.t[i] <- beta.int[4] 
        # bs_topdown_inter_skinny[[3]] # wtd.tmin1 and coy.t are independent given coy.tmin1, lion.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        # bs_topdown_inter_skinny[[9]] # moose.tmin1 and coy.t are independent given coy.tmin1, lion.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        # bs_topdown_inter_skinny[[15]] # elk.tmin1 and coy.t are independent given coy.tmin1, lion.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        # bs_topdown_inter_skinny[[28]] # bearHarv.tmin1 and coy.t are independent given coy.tmin1, lion.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i]
        # bs_topdown_inter_skinny[[34]] # bear.tmin1 and coy.t are independent given coy.tmin1, lion.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        # bs_topdown_inter_skinny[[40]] #lionHarv.tmin1 and coy.t are independent given coy.tmin1, lion.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i]
        # bs_topdown_inter_skinny[[50]] # wtd.t and coy.t are independent given wtd.tmin1, lion.tmin1, coy.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] * beta.wtd[2] * wtd.t[i]
        # bs_topdown_inter_skinny[[56]] # wolfHarv.tmin1 and coy.t are independnet given coy.tmin1, lion.tmin1 and wolf.tmin1
        # mu.coy.t[i] <- beta.int[4] + beta.coy[1] * coy.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i]
        
        coy.tmin1[i] ~ dnorm(mu.coy.tmin1[i], tau.spp.tmin1[4])
        mu.coy.tmin1[i] <- beta.int.tmin1[4] 
        
        elk.t[i] ~ dnorm(mu.elk.t[i], tau.spp[5])
        # mu.elk.t[i] <- beta.int[5] 
        # bs_topdown_inter_skinny[[2]] # wtd.tmin1 and elk.t are independent given elk.tmin1, lion.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        # bs_topdown_inter_skinny[[8]] # moose.tmin1 and elk.t are independent given elk.tmin1, lion.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        # bs_topdown_inter_skinny[[21]] # coy.tmin1 and elk.t are independent given elk.tmin1, lion.tmin1 and wolf.tmin1
        mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        # bs_topdown_inter_skinny[[27]] # bearHarv.tmin1 and elk.t are independent given elk.tmin1, lion.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i]
        # bs_topdown_inter_skinny[[33]] # bear.tmin1 and elk.t are independent given elk.tmin1, lion.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        # bs_topdown_inter_skinny[[39]] # lionHarv.tmin1 and elk.t and independent given elk.tmin1, lion.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i]
        # bs_topdown_inter_skinny[[49]] # wtd.t and elk.t are independent given wtd.tmin1, lion.tmin1, elk.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.wtd[2] * wtd.t[i]
        # bs_topdown_inter_skinny[[55]] # wolfHarv.tmin1 and elk.t are independent given elk.tmin1, lion.tmin1 and wolf.tmin1
        # mu.elk.t[i] <- beta.int[5] + beta.elk[1] * elk.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] * beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i]
        
        elk.tmin1[i] ~ dnorm(mu.elk.tmin1[i], tau.spp.tmin1[5])
        mu.elk.tmin1[i] <- beta.int.tmin1[5]
        
        moose.t[i] ~ dnorm(mu.moose.t[i], tau.spp[6])
        mu.moose.t[i] <- beta.int[6] 
        # bs_topdown_inter_skinny[[1]] # wtd.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.wtd[1] * wtd.tmin1[i]
        # bs_topdown_inter_skinny[[14]] # elk.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        # bs_topdown_inter_skinny[[20]] # coy.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        # bs_topdown_inter_skinny[[26]] # bearHarv.tmin1, moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i]
        # bs_topdown_inter_skinny[[32]] # bear.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        # bs_topdown_inter_skinny[[38]] # lionHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i]
        # bs_topdown_inter_skinny[[43]] # lion.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i]
        # bs_topdown_inter_skinny[[48]] # wtd.t and moose.t are independent given wtd.tmin1, lion.tmin1, moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.wtd[1] * wtd.tmin1[i] + beta.wtd[2] * wtd.t[i]
        # bs_topdown_inter_skinny[[54]] # wolfHarv.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
        # mu.moose.t[i] <- beta.int[6] + beta.moose[1] * moose.tmin1[i] + beta.wolf[1] * wolf.tmin1[i] + beta.harvest[1] * wolfHarv.tmin1[i]
        
        moose.tmin1[i] ~ dnorm(mu.moose.tmin1[i], tau.spp.tmin1[6])
        mu.moose.tmin1[i] <- beta.int.tmin1[6]
        
        wtd.t[i] ~ dnorm(mu.wtd.t[i], tau.spp[7])
        mu.wtd.t[i] <- beta.int[7]   
        # bs_topdown_inter_skinny[[7]] # moose.tmin1 and wtd.t are independent given wtd.tmin1 and lion.tmin1
        # mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.moose[1] * moose.tmin1[i]
        # bs_topdown_inter_skinny[[13]] # elk.tmin1 and wtd.t are independent given wtd.tmin1 and lion.tmin1
        # mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.elk[1] * elk.tmin1[i]
        # bs_topdown_inter_skinny[[19]] # coy.tmin1 and wtd.t are independent given wtd.tmin1 and lion.tmin1
        # mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.coy[1] * coy.tmin1[i]
        # bs_topdown_inter_skinny[[25]] # bearHarv.tmin1 amd wtd.t are independent given wtd.tmin1 and lion.tmin1
        # mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * bearHarv.tmin1[i]
        # bs_topdown_inter_skinny[[31]] # bear.tmin1, wtd.t are independent given wtd.tmin1 and lion.tmin1
        # mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.bear[1] * bear.tmin1[i]
        # bs_topdown_inter_skinny[[37]] # lionHarv.tmin1 and wtd.t are independent given wtd.tmin1 and lion.tmin1
        # mu.wtd.t[i] <- beta.int[7] + beta.wtd[1] * wtd.tmin1[i] + beta.lion[1] * lion.tmin1[i] + beta.harvest[1] * lionHarv.tmin1[i]
        
        wtd.tmin1[i] ~ dnorm(mu.wtd.tmin1[i], tau.spp.tmin1[7])
        mu.wtd.tmin1[i] <- beta.int.tmin1[7] 
      }
      
    }")