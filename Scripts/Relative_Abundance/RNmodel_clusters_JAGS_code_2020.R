  #'  -------------------------------
  #'  Royle-Nichols abundance model
  #'  ID CRU - Predator Interactions
  #'  Sarah B. Bassing
  #'  November 2023
  #'  -------------------------------
  #'  RN model to estimate relative abundance from binary detection/non-detection,
  #'  assuming heterogeneous abundance affects detection probability.
  #'  
  #'  Relevant parameters and data:
  #'  lambda: The number of animals available for detection at site i, N[i], is a 
  #'  Poisson-distributed random variable with mean lambda[i]. N is latent.
  #'  
  #'  y: The observed detection/non-detection data for each site i and survey 
  #'  occasion j, y[i,j], is a Bernoulli-distributed random variable with  
  #'  probability p[i,j]. y[i,j] is conditional on N[i].
  #'  
  #'  p & r: The probability of detecting occupancy at site i during occasion j, 
  #'  p[i,j], is a function of the per-individual detection probability, r[i,j], 
  #'  a binomial sampling probability that a particular individual is detected at 
  #'  site i during occasion j, and the number of individuals at site i, N[i].
  #'  -------------------------------
  
  cat(file = './Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt', "
      model{
      
      #'  Define priors
      #'  -------------
      #'  Abundance priors
      beta0 ~ dunif(-3, 3)        # Abundance intercept
      mean.lambda <- exp(beta0)   # Mean lambda for GMU10A
      # beta1 ~ dnorm(0, 0.01) 
      # beta2 ~ dnorm(0, 0.01) 
      # beta3 ~ dnorm(0, 0.01) 
      
      #'  Categorical effect for GMU needs multiple beta4 coefficients
      beta4[1] <- 0
      for(gmu in 2:ngmu) {
        beta4[gmu] ~ dnorm(0, 0.01) 
      }
      
      #'  Detection priors
      mean.r ~ dunif(0, 1)        # Detection intercept (on probability scale)
      alpha0 <- logit(mean.r)     # Detection intercept (on logit scale)
      #alpha1 ~ dnorm(0, 0.01) 
      
      #'  Categorical effect for camera setup needs multiple alpha2 coefficients
      alpha2[1] <- 0
      for(cam in 2:nsets) {
        alpha2[cam] ~ dnorm(0,0.01)  
      }
      
      
      #'  Define likelihood
      #'  -----------------
      #'  Latent state (abundance)
      for(i in 1:nsites){
        N[i] ~ dpois(lambda[i])
        lambda[i] <- exp(beta0 + beta4[gmu[i]])
        
        #'  Detection state
        for(j in 1:nsurveys){
          y[i,j] ~ dbern(p[i,j])
          p[i,j] <- 1 - pow((1 - r[i,j]), N[i])
          logit(r[i,j]) <- alpha0 + alpha2[setup[i]]  
        }
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean lambda per GMU
      for(gmu in 1:ngmu) {
        lambdaGMU[gmu] <- exp(beta0 + beta4[gmu])
      }
      
      #'  Mean lambda averaged across GMUs
      mu.lambda <- mean(lambdaGMU[])
      
      #' #'  Predicted site-level abundance per GMU 
      #' for(i in 1:ncams1) {
      #'   # Ngmu1[i] <- exp(beta0 + beta1 * forest[i] + beta2 * elev[i] + beta3 * pow(elev[i],2) + beta4[gmu[1]])
      #'   Ngmu1[i] <- exp(beta0 + beta4[gmu[1]])
      #' }
      #' for(i in 1:ncams2) {
      #'   # Ngmu2[i] <- exp(beta0 + beta1 * forest[i] + beta2 * elev[i] + beta3 * pow(elev[i],2) + beta4[gmu[2]])
      #'   Ngmu2[i] <- exp(beta0 + beta4[gmu[2]])
      #' }
      
      #'  Total abundance across camera sites
      totalN <- sum(N[])
      
      #'  Mean density per cluster
      for(i in 1:nsites) {
        N.cl1[i] <- N[i] * cluster_matrix[i,1] # GMU10A cams
        N.cl2[i] <- N[i] * cluster_matrix[i,2] 
        N.cl3[i] <- N[i] * cluster_matrix[i,3] 
        N.cl4[i] <- N[i] * cluster_matrix[i,4] 
        N.cl5[i] <- N[i] * cluster_matrix[i,5] 
        N.cl6[i] <- N[i] * cluster_matrix[i,6] 
        N.cl7[i] <- N[i] * cluster_matrix[i,7] 
        N.cl8[i] <- N[i] * cluster_matrix[i,8] 
        N.cl9[i] <- N[i] * cluster_matrix[i,9] 
        # N.cl10[i] <- N[i] * cluster_matrix[i,10] # GMU1 cams
        # N.cl11[i] <- N[i] * cluster_matrix[i,11] 
        # N.cl12[i] <- N[i] * cluster_matrix[i,12] 
        # N.cl13[i] <- N[i] * cluster_matrix[i,13] 
        # N.cl14[i] <- N[i] * cluster_matrix[i,14] 
        # N.cl15[i] <- N[i] * cluster_matrix[i,15] 
        # N.cl16[i] <- N[i] * cluster_matrix[i,16] 
        # N.cl17[i] <- N[i] * cluster_matrix[i,17] 
        # N.cl18[i] <- N[i] * cluster_matrix[i,18] 
        # N.cl19[i] <- N[i] * cluster_matrix[i,19] 
        # N.cl20[i] <- N[i] * cluster_matrix[i,20] 
        N.cl21[i] <- N[i] * cluster_matrix[i,21] # GMU6 cams
        N.cl22[i] <- N[i] * cluster_matrix[i,22] 
        N.cl23[i] <- N[i] * cluster_matrix[i,23] 
        N.cl24[i] <- N[i] * cluster_matrix[i,24] 
      }
      
      #'  Estimate relative abundance index by summing N per cluster, dividing
      #'  the sum of camera cluster indicators (i.e., number of cameras in a cluster) 
      #'  to estimate the average abundance/camera, then dividing by the area of 
      #'  each cluster to estimate relative density index (per sq-km) and finally 
      #'  multiplying by 100 so RDI per 100 sq-km
      rdi.cl1 <- ((sum(N.cl1) / sum(cluster_matrix[,1])) / cluster_area[1]) * 100 # GMU10A clusters
      rdi.cl2 <- ((sum(N.cl2) / sum(cluster_matrix[,2])) / cluster_area[2]) * 100
      rdi.cl3 <- ((sum(N.cl3) / sum(cluster_matrix[,3])) / cluster_area[3]) * 100
      rdi.cl4 <- ((sum(N.cl4) / sum(cluster_matrix[,4])) / cluster_area[4]) * 100
      rdi.cl5 <- ((sum(N.cl5) / sum(cluster_matrix[,5])) / cluster_area[5]) * 100
      rdi.cl6 <- ((sum(N.cl6) / sum(cluster_matrix[,6])) / cluster_area[6]) * 100
      rdi.cl7 <- ((sum(N.cl7) / sum(cluster_matrix[,7])) / cluster_area[7]) * 100
      rdi.cl8 <- ((sum(N.cl8) / sum(cluster_matrix[,8])) / cluster_area[8]) * 100
      rdi.cl9 <- ((sum(N.cl9) / sum(cluster_matrix[,9])) / cluster_area[9]) * 100
      # rdi.cl10 <- ((sum(N.cl10) / sum(cluster_matrix[,10])) / cluster_area[10]) * 100 # GMU 1 clusters
      # rdi.cl11 <- ((sum(N.cl11) / sum(cluster_matrix[,11])) / cluster_area[11]) * 100
      # rdi.cl12 <- ((sum(N.cl12) / sum(cluster_matrix[,12])) / cluster_area[12]) * 100
      # rdi.cl13 <- ((sum(N.cl13) / sum(cluster_matrix[,13])) / cluster_area[13]) * 100
      # rdi.cl14 <- ((sum(N.cl14) / sum(cluster_matrix[,14])) / cluster_area[14]) * 100
      # rdi.cl15 <- ((sum(N.cl15) / sum(cluster_matrix[,15])) / cluster_area[15]) * 100
      # rdi.cl16 <- ((sum(N.cl16) / sum(cluster_matrix[,16])) / cluster_area[16]) * 100
      # rdi.cl17 <- ((sum(N.cl17) / sum(cluster_matrix[,17])) / cluster_area[17]) * 100
      # rdi.cl18 <- ((sum(N.cl18) / sum(cluster_matrix[,18])) / cluster_area[18]) * 100
      # rdi.cl19 <- ((sum(N.cl19) / sum(cluster_matrix[,19])) / cluster_area[19]) * 100
      # rdi.cl20 <- ((sum(N.cl20) / sum(cluster_matrix[,20])) / cluster_area[20]) * 100
      rdi.cl21 <- ((sum(N.cl21) / sum(cluster_matrix[,21])) / cluster_area[21]) * 100 # GMU6 clusters
      rdi.cl22 <- ((sum(N.cl22) / sum(cluster_matrix[,22])) / cluster_area[22]) * 100
      rdi.cl23 <- ((sum(N.cl23) / sum(cluster_matrix[,23])) / cluster_area[23]) * 100
      rdi.cl24 <- ((sum(N.cl24) / sum(cluster_matrix[,24])) / cluster_area[24]) * 100
       
      
      # totalN.gmu10a <- sum(Ngmu1[])
      # densitykm2.gmu10a <- totalN.gmu10a/area1
      # density100km2.gmu10a <- densitykm2.gmu10a * 100
      # 
      # totalN.gmu6 <- sum(Ngmu2[])
      # densitykm2.gmu6 <- totalN.gmu6/area2
      # density100km2.gmu6 <- densitykm2.gmu6 * 100

      #' #'  Total sites occupied (N > 0)
      #' for(i in 1:nsites) {
      #'   occupied[i] <- ifelse(N[i] > 0, 1, 0)
      #' }
      #' occSites <- sum(occupied[])
      #' 
      #' #'  Mean occupancy probability
      #' mean.psi <- 1 - exp(-mu.lambda)
 
      #'  Mean per-individual detection probability (r) per camera setup
      for(cam in 1:nsets) {
        rSetup[cam] <- 1/(1 + exp(-(alpha0 + alpha2[cam])))
      }
      #'  per-individual detection probability (r) averaged across all camera setups 
      mu.r <- mean(rSetup[])

      #'  Mean detection probability (p)
      for(i in 1:nsites) {
        p.occasion[i] <- mean(p[i,])
      }
      mean.p <- mean(p.occasion[])
      
    }
  ")
