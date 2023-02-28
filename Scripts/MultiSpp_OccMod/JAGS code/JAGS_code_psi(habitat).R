  #'  ------------------------------------
  #'  Habitat model
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  ------------------------------------
  #'  Model to test whether predator co-occurrence is non-independent and whether
  #'  basic habitat features influence that relationship.
  #'  ------------------------------------
  
  cat(file = './Outputs/MultiSpp_OccMod_Outputs/JAGS_output/psi(habitat)_p(effort).txt', "
      model{
      
      ####  Define priors  ####
      #'  -----------------
      #'  Intercepts & slopes associated with each natural parameter
      #'  First order psi intercepts
      mean.psi1 ~ dunif(0, 1)
      mean.psi2 ~ dunif(0, 1)
      betaSpp1 <- logit(mean.psi1)
      betaSpp2 <- logit(mean.psi2)
      
      #'  First order psi slopes
      for(fo_psi in 2:nfirst_order_psi) {
        betaSpp1[fo_psi] ~ dnorm(0, 0.001)
        betaSpp2[fo_psi] ~ dnorm(0, 0.001)
      }
      
      #'  Second order psi priors
      for(so_psi in 1:nsecond_order_psi) {
        betaSpp12[so_psi] ~ dnorm(0, 0.001)
      }
      
      #'  First order rho intercepts
      mean.p1 ~ dunif(0, 1)
      mean.p2 ~ dunif(0, 1)
      alphaSpp1 <- logit(mean.p1)
      alphaSpp2 <- logit(mean.p2)
      
      #'  First order rho slopes
      for(fo_rho in 2:nfirst_order_rho) {
        alphaSpp1[fo_rho] ~ dnorm(0, 0.001)
        alphaSpp2[fo_rho] ~ dnorm(0, 0.001)
      }
      
      #'  Second order rho priors
      for(so_rho in 1:nsecond_order_psi) {
        alpha12[so_rho] ~ dnorm(0, 0.001)
      }

      
      ####  Likelihood  ####
      #'  --------------
      #'  1. Basic heirarchical model
      #'  Latent state model
      for(i in 1:nsites) {
        z[i] ~ dcat(lsv[i, (1:ncat)])
      }
      
      #'  Detection model
      for(i in 1:nsites) {
        for(j in 1:nsurveys) {
          y[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]])
        }
      }
      
      #'  2. Define latent state vector & observation matrices
      for(i in 1:nsites) {
        #'  Latent stete probabilities in latent state vector (lsv)
        #'  Probabilities for each state (ncat)
        lsv[i,1] <- 1
        lsv[i,2] <- exp(psi1[i])
        lsv[i,3] <- exp(psi2[i])
        lsv[i,4] <- exp(psi12[i])
      
        for(j in 1:nsurveys) {
          #'  Detection matrix (OS = observation state, TS = true state)
          #'  rdm = rho detection matrix; each row sums to 1
          #'  OS along rows, TS, along columns
          #'  True state = 00
          
          #'  True state = 10
      
          #'  True state = 01
      
          #'  True state = 11
        }
      }
      
      
      
      
      
      
      
      
      
      
      }
      ")