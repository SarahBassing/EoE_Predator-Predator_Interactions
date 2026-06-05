  #'  -------------------------------------
  #'  Bayesian Structural Equation Models
  #'  Sarah Bassing
  #'  December 2025
  #'  -------------------------------------
  #'  Source data formatting script and run structural equation models (SEM) in
  #'  a Bayesian framework to test hypotheses about how predator-prey and predator- 
  #'  predator interactions influence wildlife populations in northern Idaho. This
  #'  formulation relies on stacked data structure (stacking 2020/2021 and 
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
  source("./Scripts/Structural_Equation_Models/Format_spatial_covariates_for_SEMs.R") 
  
  #'  Run script that formats density data for SEMs
  # source("./Scripts/Structural_Equation_Models/Format_density_data_for_SEMs.R")
  source("./Scripts/Structural_Equation_Models/Format_RNmodel_Posteriors_for_SEM.R")
  
  #' #'  Take a quick look
  #' head(density_wide_1YrLag_20s_22s)
  
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
  
  # dat_final <- density_wide_1YrLag_20s_22s %>%
  #   mutate(obs = seq(1:nrow(.)),
  #          GMU_cluster = paste0(GMU, "_", ClusterID)) %>%
  #   group_by(GMU_cluster) %>%
  #   mutate(uniqueCluster = cur_group_id()) %>%
  #   ungroup() #%>%
  #   # arrange(obs)
  # 
  # head(dat_final)
  
  
  #'  ------------------------
  ####  Setup data for JAGS  ####
  #'  ------------------------
  #'  Bundle data for JAGS
  bundle_dat <- function(dat, covs, nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor, nwsi) {
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
                    nWSI = nwsi,
                    nCluster = as.numeric(length(unique(dat$wolf$time_t$cluster))),
                    nSpp = 7,
                    #'  Standardized posterior means and SD for each species and time step
                    wolf.t_hat = dat$wolf$time_t$posterior_mu_z,
                    wolf.t.sigma_hat = dat$wolf$time_t$posterior_sd_z,
                    wolf.tmin1_hat = dat$wolf$time_tmin1$posterior_mu_z,
                    wolf.tmin1.sigma_hat = dat$wolf$time_tmin1$posterior_sd_z,
                    lion.t_hat = dat$lion$time_t$posterior_mu_z,
                    lion.t.sigma_hat = dat$lion$time_t$posterior_sd_z,
                    lion.tmin1_hat = dat$lion$time_tmin1$posterior_mu_z,
                    lion.tmin1.sigma_hat = dat$lion$time_tmin1$posterior_sd_z,
                    bear.t_hat = dat$bear$time_t$posterior_mu_z,
                    bear.t.sigma_hat = dat$bear$time_t$posterior_sd_z,
                    bear.tmin1_hat = dat$bear$time_tmin1$posterior_mu_z,
                    bear.tmin1.sigma_hat = dat$bear$time_tmin1$posterior_sd_z,
                    coy.t_hat = dat$coy$time_t$posterior_mu_z,
                    coy.t.sigma_hat = dat$coy$time_t$posterior_sd_z, 
                    coy.tmin1_hat = dat$coy$time_tmin1$posterior_mu_z,
                    coy.tmin1.sigma_hat = dat$coy$time_tmin1$posterior_sd_z,
                    elk.t_hat = dat$elk$time_t$posterior_mu_z,
                    elk.t.sigma_hat = dat$elk$time_t$posterior_sd_z,
                    elk.tmin1_hat = dat$elk$time_tmin1$posterior_mu_z,
                    elk.tmin1.sigma_hat = dat$elk$time_tmin1$posterior_sd_z,
                    moose.t_hat = dat$moose$time_t$posterior_mu_z,
                    moose.t.sigma_hat = dat$moose$time_t$posterior_sd_z,
                    moose.tmin1_hat = dat$moose$time_tmin1$posterior_mu_z,
                    moose.tmin1.sigma_hat = dat$moose$time_tmin1$posterior_sd_z,
                    wtd.t_hat = dat$wtd$time_t$posterior_mu_z,
                    wtd.t.sigma_hat = dat$wtd$time_t$posterior_sd_z,
                    wtd.tmin1_hat = dat$wtd$time_tmin1$posterior_mu_z,
                    wtd.tmin1.sigma_hat = dat$wtd$time_tmin1$posterior_sd_z,
                    #'  Standardized harvest and habitat variables for 1 year time lag
                    wolfHarv.tmin1 = covs$wolf_harv$time_tmin1$wolf_harv_z,
                    lionHarv.tmin1 = covs$lion_harv$time_tmin1$lion_harv_z,
                    bearHarv.tmin1 = covs$bear_harv$time_tmin1$bear_harv_z,
                    elkHarv.tmin1 = covs$elk_harv$time_tmin1$elk_harv_z,
                    mooseHarv.tmin1 = covs$moose_harv$time_tmin1$moose_harv_z,
                    deerHarv.tmin1 = covs$deer_harv$time_tmin1$deer_harv_z,
                    wsi.tmin1 = covs$wsi$time_tmin1$wsi_z,
                    forest.tmin1 = covs$prop_disturbed$time_tmin1$prop_disturb_z,
                    road.tmin1 = covs$road_density$time_tmin1$road_density_z,
                    public.tmin1 = covs$public_land$time_tmin1$public_land_z)
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle_topinter <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 7, nlion = 4, nbear = 5, 
                                          ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 6, nfor = 1, nwsi = 1)
  data_JAGS_bundle_topexploit <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 7, nlion = 4, nbear = 5, 
                                            ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 6, nfor = 1, nwsi = 1)
  data_JAGS_bundle_bottominter <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 7, nlion = 4, nbear = 5, 
                                             ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 0, nfor = 1, nwsi = 1)
  data_JAGS_bundle_bottomexploit <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 7, nlion = 4, nbear = 5, 
                                               ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 0, nfor = 1, nwsi = 1)
                                 
  # save(data_JAGS_bundle_topinter, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_topinter.RData")
  # save(data_JAGS_bundle_topexploit, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_topexploit.RData")
  # save(data_JAGS_bundle_bottominter, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_bottominter.RData")
  # save(data_JAGS_bundle_bottomexploit, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_bottomexploit.RData")
  
  #'  Generate initial values for each parameter (random node)
  generate_inits <- function(nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor, nwsi) {
    
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
      wolfharvest.t_1 = runif(nharv, -0.5, 0.5),
      lionharvest.t_1 = runif(nharv, -0.5, 0.5),
      bearharvest.t_1 = runif(nharv, -0.5, 0.5),
      elkharvest.t_1 = runif(nharv, -0.5, 0.5),
      mooseharvest.t_1 = runif(nharv, -0.5, 0.5),
      deerharvest.t_1 = runif(nharv, -0.5, 0.5),
      wsi.t_1 = runif(nwsi, -0.5, 0.5),
      forest.t_1 = runif(nfor, -0.5, 0.5),
      road.t_1 = runif(nwsi, -0.5, 0.5),
      public.t_1 = runif(nwsi, -0.5, 0.5) # doesn't matter that using nwsi or other b/c all = 1
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
  initsList_topinter <- initsList_topinter_reduced <- initsList_topexploit <- initsList_bottominter <- initsList_bottomexploit <- vector('list', num.chains)
  #'  Setting seed for reproducibility
  set.seed(9983)
  #'  Loop through generate_inits function 3 times (1 for each chain) 
  for(i in 1:num.chains){
    initsList_topinter[[i]] <- generate_inits(nwolf = 7, nlion = 4, nbear = 5, ncoy = 2, nelk = 1, 
                                     nmoose = 1, nwtd = 1, nharv = 6, nfor = 0, nwsi = 0)
    initsList_topinter_reduced[[i]] <- generate_inits(nwolf = 6, nlion = 4, nbear = 1, ncoy = 1, nelk = 1, 
                                              nmoose = 1, nwtd = 1, nharv = 3, nfor = 0, nwsi = 0)
    initsList_topexploit[[i]] <- generate_inits(nwolf = 4, nlion = 3, nbear = 3, ncoy = 2, nelk = 1, 
                                              nmoose = 1, nwtd = 1, nharv = 6, nfor = 0, nwsi = 0)
    initsList_bottominter[[i]] <- generate_inits(nwolf = 4, nlion = 2, nbear = 3, ncoy = 1, nelk = 4, 
                                              nmoose = 2, nwtd = 5, nharv = 0, nfor = 1, nwsi = 1)
    initsList_bottomexploit[[i]] <- generate_inits(nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, 
                                              nmoose = 2, nwtd = 5, nharv = 0, nfor = 1, nwsi = 1)
  }
  
  #'  Parameters monitored
  params <- c("beta.int", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk", 
              "beta.moose", "beta.wtd", "beta.harvest", "beta.wsi","beta.forest", 
              "sigma.spp", "sigma.cluster", "cluster.randeff") 
  
  #'  MCMC settings
  nc <- 3
  ni <- 100000
  nb <- 50000
  nt <- 10
  na <- 5000
  
  
  #'  ---------------------------------
  ####  Call JAGS & Fit Bayesian SEMs  ####
  #'  ---------------------------------
  #####  Top-down, interference model  #####
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter.R")
  start.time = Sys.time()
  SEM_topdown_inter <- jags(data_JAGS_bundle_topinter, inits = initsList_topinter, params, 
                            "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_inter$summary)  
  which(SEM_topdown_inter$summary[,"Rhat"] > 1.1)    
  mcmcplot(SEM_topdown_inter$samples)
  save(SEM_topdown_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_inter_", Sys.Date(), ".RData"))
  
  #'  Notes: 
  #'  Black bear intercept & auto-regressive term struggle to converge w <75000 iterations
  #'  sigma.clusters struggle to converge w <100000 iterations & 50000 burnin
  #'  sigma.clusters seem large compared to sigma.spp and betas (mostly ranging 7 - 9 vs 0.1 - 0.9)
  
  #####  Top-down, interference model  #####
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter_reduced.R")
  start.time = Sys.time()
  SEM_topdown_inter <- jags(data_JAGS_bundle_topinter, inits = initsList_topinter_reduced, params, 
                            "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_reduced.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_inter$summary)  
  which(SEM_topdown_inter$summary[,"Rhat"] > 1.1)    
  mcmcplot(SEM_topdown_inter$samples)
  save(SEM_topdown_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_inter_", Sys.Date(), ".RData"))
  
  
  #'  Same model but no random effect for cluster
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter_no_clusterRE.R")
  start.time = Sys.time()
  SEM_topdown_inter <- jags(data_JAGS_bundle_topinter, inits = initsList_topinter, params,
                            "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_no_clusterRE.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni,
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_inter$summary)
  which(SEM_topdown_inter$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_topdown_inter$samples)
  save(SEM_topdown_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_inter_no_clusterRE_", Sys.Date(), ".RData"))
  

  #####  Top-down, exploitative model  #####  
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_exploit.R")
  start.time = Sys.time()
  SEM_topdown_exploit <- jags(data_JAGS_bundle_topexploit, inits = initsList_topexploit, params, 
                              "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_exploit.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_exploit$summary)
  which(SEM_topdown_exploit$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_topdown_exploit$samples)
  save(SEM_topdown_exploit, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_exploit_", Sys.Date(), ".RData"))
  
  #####  Bottom-up, interference model  #####  
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_inter.R")
  start.time = Sys.time()
  SEM_bottomup_inter <- jags(data_JAGS_bundle_bottominter, inits = initsList_bottominter, params, 
                             "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_inter.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_inter$summary)
  which(SEM_bottomup_inter$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_inter$samples)
  save(SEM_bottomup_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_inter_", Sys.Date(), ".RData"))
  
  #####  Bottom-up, exploitative model  #####  
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_exploit.R")
  start.time = Sys.time()
  SEM_bottomup_exploit <- jags(data_JAGS_bundle_bottomexploit, inits = initsList_bottomexploit, params, 
                               "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_exploit.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_exploit$summary)
  which(SEM_bottomup_exploit$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_exploit$samples)
  save(SEM_bottomup_exploit, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_exploit_", Sys.Date(), ".RData"))
  
  
  
  
  
  
  