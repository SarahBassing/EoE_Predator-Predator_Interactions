  #'  -------------------------------------
  #'  Bayesian Structural Equation Models
  #'  Sarah Bassing
  #'  December 2025
  #'  -------------------------------------
  #'  Source data formatting script and run structural equation models (SEM) in
  #'  a Bayesian framework to test hypotheses about how predator-prey and predator- 
  #'  predator interactions influence wildlife populations in northern Idaho. This
  #'  formulation relies on original data structure (stacking 2020/2021 and 
  #'  2021/2022 for Yr1 --> Yr2 effect).
  #'  Species index order: 
  #'    1 = wolf; 2 = cougar; 3 = black bear; 4 = coyote; 5 = elk; 6 = moose; 7 = wtd
  #'  -------------------------------------
  
  #'  Clean workspace
  rm(list = ls())

  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  Run script that formats covariate data
  source("./Scripts/Structural_Equation_Models/Format_covariate_data_for_SEMs.R") 
  
  #'  Run script that formats density data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_density_data_for_SEMs.R")
  
  #'  Take a quick look
  head(density_wide_1YrLag_20s_22s)
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  # dat_final <- density_wide_1YrLag_20s_22s %>%
  #   mutate(GMU_cluster = paste0(GMU, "_", ClusterID)) %>%
  #   group_by(GMU_cluster) %>%
  #   mutate(uniqueCluster = cur_group_id()) %>%
  #   ungroup() %>%
  #   relocate(uniqueCluster, .after = ClusterID) %>%
  #   dplyr::select(-GMU_cluster) %>%
  #   arrange(uniqueCluster, timestep)
  
  dat_final <- density_wide_1YrLag_20s_22s %>%
    mutate(obs = seq(1:nrow(.)),
           GMU_cluster = paste0(GMU, "_", ClusterID)) %>%
    group_by(GMU_cluster) %>%
    mutate(uniqueCluster = cur_group_id()) %>%
    ungroup() #%>%
    # arrange(obs)
  
  head(dat_final)
    
  
  #'  ------------------------
  #####  Setup data for JAGS  #####
  #'  ------------------------
  #'  Bundle data for JAGS
  bundle_dat <- function(dat, nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor) {
    #'  Bundle data for JAGS
    bundled <- list(nWolf = nwolf,
                    nLion = nlion, 
                    nBear = nbear, 
                    nCoy = ncoy, 
                    nElk = nelk,
                    nMoose = nmoose, 
                    nWtd = nwtd, 
                    nharvest = nharv,
                    nforest = nfor,
                    nCluster = as.numeric(length(unique(dat$uniqueCluster))), 
                    nSpp = 7,
                    #'  Standardize and set to numeric for each variable
                    wolf.t = as.numeric(scale(dat$wolf.T)),
                    wolf.tmin1 = as.numeric(scale(dat$wolf.Tminus1)),
                    lion.t = as.numeric(scale(dat$mountain_lion.T)),
                    lion.tmin1 = as.numeric(scale(dat$mountain_lion.Tminus1)),
                    bear.t = as.numeric(scale(dat$bear_black.T)), 
                    bear.tmin1 = as.numeric(scale(dat$bear_black.Tminus1)), 
                    coy.t = as.numeric(scale(dat$coyote.T)), 
                    coy.tmin1 = as.numeric(scale(dat$coyote.Tminus1)), 
                    elk.t = as.numeric(scale(dat$elk.T)), 
                    elk.tmin1 = as.numeric(scale(dat$elk.Tminus1)), 
                    moose.t = as.numeric(scale(dat$moose.T)), 
                    moose.tmin1 = as.numeric(scale(dat$moose.Tminus1)), 
                    wtd.t = as.numeric(scale(dat$whitetailed_deer.T)), 
                    wtd.tmin1 = as.numeric(scale(dat$whitetailed_deer.Tminus1)), 
                    harvest.t = as.numeric(scale(dat$annual_harvest.T)), 
                    harvest.tmin1 = as.numeric(scale(dat$annual_harvest.Tminus1)), 
                    forest.t = as.numeric(scale(dat$DisturbedForest_last20Yrs.T)), 
                    forest.tmin1 = as.numeric(scale(dat$DisturbedForest_last20Yrs.Tminus1)))
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle_topinter <- bundle_dat(dat_final, nwolf = 7, nlion = 4, nbear = 5, 
                                          ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
  data_JAGS_bundle_topexploit <- bundle_dat(dat_final, nwolf = 7, nlion = 4, nbear = 5, 
                                            ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
  data_JAGS_bundle_bottominter <- bundle_dat(dat_final, nwolf = 7, nlion = 4, nbear = 5, 
                                             ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
  data_JAGS_bundle_bottomexploit <- bundle_dat(dat_final, nwolf = 7, nlion = 4, nbear = 5, 
                                               ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 1, nfor = 1)
                                 
  
  # save(data_JAGS_bundle_topinter, file = "./Data/Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_topinter.RData")
  # save(data_JAGS_bundle_topexploit, file = "./Data/Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_topexploit.RData")
  # save(data_JAGS_bundle_bottominter, file = "./Data/Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_bottominter.RData")
  # save(data_JAGS_bundle_bottomexploit, file = "./Data/Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_bottomexploit.RData")
  
  #'  Generate initial values for each parameter (random node)
  generate_inits <- function(nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor) {
    
    #'  Generate random values for each species-specific beta (nwolf, nlion, etc.
    #'  based on number of species-specific betas to be estimated)
    list(
      wolf.t_1 = runif(nwolf, -0.5, 0.5),  # consider -1, 1
      lion.t_1 = runif(nlion, -0.5, 0.5),
      bear.t_1 = runif(nbear, -0.5, 0.5),
      coy.t_1 = runif(ncoy, -0.5, 0.5),
      elk.t_1 = runif(nelk, -0.5, 0.5),
      moose.t_1 = runif(nmoose, -0.5, 0.5),
      wtd.t_1 = runif(nwtd, -0.5, 0.5),
      harvest.t_1 = runif(nharv, -0.5, 0.5),
      forest.t_1 = runif(nfor, -0.5, 0.5)#,
      #' #'  Fix random number generator and seed for every run of this function
      #' .RNG.name = "base::Wichmann-Hill",
      #' .RNG.seed = 182  
      #' #'  Setting RNG seed leads to the same random number stream during
      #' #'  adaptation and sampling, causing different inits to be rapidly erased
      #' #'  and all chains to follow the same deterministic path. This does not 
      #' #'  appear to happen with the cross lag model though.
    )
  }
  #'  Define number of chains
  num.chains <- 3
  #'  Create empty lists
  initsList_topinter <- initsList_topexploit <- initsList_bottominter <- initsList_bottomexploit <- vector('list', num.chains)
  #'  Setting seed for reproducibility
  set.seed(9983)
  #'  Loop through generate_inits function 3 times (1 for each chain) 
  for(i in 1:num.chains){
    initsList_topinter[[i]] <- generate_inits(nwolf = 7, nlion = 4, nbear = 5, ncoy = 2, nelk = 1, 
                                     nmoose = 1, nwtd = 1, nharv = 1, nfor = 0)
    initsList_topexploit[[i]] <- generate_inits(nwolf = 4, nlion = 3, nbear = 3, ncoy = 2, nelk = 1, 
                                              nmoose = 1, nwtd = 1, nharv = 1, nfor = 0)
    initsList_bottominter[[i]] <- generate_inits(nwolf = 4, nlion = 2, nbear = 3, ncoy = 1, nelk = 4, 
                                              nmoose = 2, nwtd = 5, nharv = 0, nfor = 1)
    initsList_bottomexploit[[i]] <- generate_inits(nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, 
                                              nmoose = 2, nwtd = 5, nharv = 0, nfor = 1)
  }
  
  
  #'  Parameters monitored
  params <- c("beta.int", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk", 
              "beta.moose", "beta.wtd", "beta.harvest", "beta.forest", "sigma.spp", "sigma.cluster") 
  
  #'  MCMC settings
  nc <- 3
  ni <- 50000
  nb <- 10000
  nt <- 10
  na <- 5000
  
  #'  Fit Bayesian SEMs
  #'  ------------------
  ####  Top-down, interference model  ####
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter.R")
  start.time = Sys.time()
  SEM_topdown_inter <- jags(data_JAGS_bundle, params, inits = initsList_topinter, 
                            "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_inter$summary)
  which(SEM_topdown_inter$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_topdown_inter$samples)
  save(SEM_topdown_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_inter_", Sys.Date(), ".RData"))
  
  ####  Top-down, exploitative model  ####
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_exploit.R")
  start.time = Sys.time()
  SEM_topdown_exploit <- jags(data_JAGS_bundle, params, inits = initsList_topexploit, 
                              "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_exploit.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_exploit$summary)
  which(SEM_topdown_exploit$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_topdown_exploit$samples)
  save(SEM_topdown_exploit, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_exploit_", Sys.Date(), ".RData"))
  
  ####  Bottom-up, interference model  ####
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_inter.R")
  start.time = Sys.time()
  SEM_bottomup_inter <- jags(data_JAGS_bundle, params, inits = initsList_bottominter, 
                             "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_inter.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_inter$summary)
  which(SEM_bottomup_inter$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_inter$samples)
  save(SEM_bottomup_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_inter_", Sys.Date(), ".RData"))
  
  ####  Bottom-up, exploitative model  ####
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_exploit.R")
  start.time = Sys.time()
  SEM_bottomup_exploit <- jags(data_JAGS_bundle, params, inits = initsList_bottomexploit, 
                               "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_exploit.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_exploit$summary)
  which(SEM_bottomup_exploit$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_exploit$samples)
  save(SEM_bottomup_exploit, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_exploit_", Sys.Date(), ".RData"))
  
  
  
  
  
  
  