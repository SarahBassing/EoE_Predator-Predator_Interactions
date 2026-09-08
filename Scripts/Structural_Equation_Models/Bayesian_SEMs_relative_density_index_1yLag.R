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
  
  #'  Run script that formats RDI posteriors, bundles for JAGS, and draws inits
  source("./Scripts/Structural_Equation_Models/Format_RNmodel_Posteriors_for_SEM.R")
  
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
  
  
  #' #'  ------------------------
  #' ####  Setup data for JAGS  ####
  #' #'  ------------------------
  #' #'  Bundle data for JAGS
  #' bundle_dat <- function(dat, covs, nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor, nwsi) {
  #'   #'  Bundle data for JAGS
  #'   bundled <- list(nWolf = nwolf,
  #'                   nLion = nlion, 
  #'                   nBear = nbear, 
  #'                   nCoy = ncoy, 
  #'                   nElk = nelk,
  #'                   nMoose = nmoose, 
  #'                   nWtd = nwtd, 
  #'                   nharvest = nharv,
  #'                   nforest = nfor,
  #'                   nWSI = nwsi,
  #'                   nCluster = as.numeric(length(unique(dat$wolf$time_t$cluster))),
  #'                   nSpp = 7,
  #'                   #'  Standardized posterior means and SD for each species and time step
  #'                   wolf.t_hat = dat$wolf$time_t$posterior_mu_z,
  #'                   wolf.t.sigma_hat = dat$wolf$time_t$posterior_sd_z,
  #'                   wolf.tmin1_hat = dat$wolf$time_tmin1$posterior_mu_z,
  #'                   wolf.tmin1.sigma_hat = dat$wolf$time_tmin1$posterior_sd_z,
  #'                   lion.t_hat = dat$lion$time_t$posterior_mu_z,
  #'                   lion.t.sigma_hat = dat$lion$time_t$posterior_sd_z,
  #'                   lion.tmin1_hat = dat$lion$time_tmin1$posterior_mu_z,
  #'                   lion.tmin1.sigma_hat = dat$lion$time_tmin1$posterior_sd_z,
  #'                   bear.t_hat = dat$bear$time_t$posterior_mu_z,
  #'                   bear.t.sigma_hat = dat$bear$time_t$posterior_sd_z,
  #'                   bear.tmin1_hat = dat$bear$time_tmin1$posterior_mu_z,
  #'                   bear.tmin1.sigma_hat = dat$bear$time_tmin1$posterior_sd_z,
  #'                   coy.t_hat = dat$coy$time_t$posterior_mu_z,
  #'                   coy.t.sigma_hat = dat$coy$time_t$posterior_sd_z, 
  #'                   coy.tmin1_hat = dat$coy$time_tmin1$posterior_mu_z,
  #'                   coy.tmin1.sigma_hat = dat$coy$time_tmin1$posterior_sd_z,
  #'                   elk.t_hat = dat$elk$time_t$posterior_mu_z,
  #'                   elk.t.sigma_hat = dat$elk$time_t$posterior_sd_z,
  #'                   elk.tmin1_hat = dat$elk$time_tmin1$posterior_mu_z,
  #'                   elk.tmin1.sigma_hat = dat$elk$time_tmin1$posterior_sd_z,
  #'                   moose.t_hat = dat$moose$time_t$posterior_mu_z,
  #'                   moose.t.sigma_hat = dat$moose$time_t$posterior_sd_z,
  #'                   moose.tmin1_hat = dat$moose$time_tmin1$posterior_mu_z,
  #'                   moose.tmin1.sigma_hat = dat$moose$time_tmin1$posterior_sd_z,
  #'                   wtd.t_hat = dat$wtd$time_t$posterior_mu_z,
  #'                   wtd.t.sigma_hat = dat$wtd$time_t$posterior_sd_z,
  #'                   wtd.tmin1_hat = dat$wtd$time_tmin1$posterior_mu_z,
  #'                   wtd.tmin1.sigma_hat = dat$wtd$time_tmin1$posterior_sd_z,
  #'                   #'  Standardized harvest and habitat variables for 1 year time lag
  #'                   wolfHarv.tmin1 = covs$wolf_harv$time_tmin1$wolf_harv_z,
  #'                   lionHarv.tmin1 = covs$lion_harv$time_tmin1$lion_harv_z,
  #'                   bearHarv.tmin1 = covs$bear_harv$time_tmin1$bear_harv_z,
  #'                   elkHarv.tmin1 = covs$elk_harv$time_tmin1$elk_harv_z,
  #'                   mooseHarv.tmin1 = covs$moose_harv$time_tmin1$moose_harv_z,
  #'                   deerHarv.tmin1 = covs$deer_harv$time_tmin1$deer_harv_z,
  #'                   wsi.tmin1 = covs$wsi$time_tmin1$wsi_z,
  #'                   forest.tmin1 = covs$prop_disturbed$time_tmin1$prop_disturb_z,
  #'                   road.tmin1 = covs$road_density$time_tmin1$road_density_z,
  #'                   public.tmin1 = covs$public_land$time_tmin1$public_land_z)
  #'   str(bundled)
  #'   return(bundled)
  #' }
  #' data_JAGS_bundle_ar1 <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 1, nlion = 1, nbear = 1, 
  #'                                    ncoy = 1, nelk = 1, nmoose = 1, nwtd = 1, nharv = 0, nfor = 0, nwsi = 0)
  #' data_JAGS_bundle_top <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 4, nlion = 3, nbear = 3, 
  #'                                         ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 6, nfor = 0, nwsi = 0)
  #' data_JAGS_bundle_bottom <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 1, nlion = 1, nbear = 1, 
  #'                                       ncoy = 1, nelk = 5, nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  #' data_JAGS_bundle_topinter <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 7, nlion = 4, nbear = 5, 
  #'                                         ncoy = 2, nelk = 1, nmoose = 1, nwtd = 1, nharv = 6, nfor = 1, nwsi = 1)
  #' data_JAGS_bundle_topinter_reduced <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 6, nlion = 4, nbear = 2, 
  #'                                         ncoy = 1, nelk = 1, nmoose = 1, nwtd = 1, nharv = 3, nfor = 1, nwsi = 1)
  #' # data_JAGS_bundle_topexploit <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 4, nlion = 3, nbear = 1, 
  #' #                                           ncoy = 1, nelk = 1, nmoose = 1, nwtd = 1, nharv = 3, nfor = 1, nwsi = 1)
  #' data_JAGS_bundle_bottominter <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 4, nlion = 1, nbear = 1, 
  #'                                            ncoy = 1, nelk = 5, nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  #' data_JAGS_bundle_bottominter_reduced <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 4, nlion = 1, nbear = 1, 
  #'                                                    ncoy = 1, nelk = 2, nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  #' data_JAGS_bundle_topbottom <- bundle_dat(post_summaries, covs = covs_ztransformed, nwolf = 7, nlion = 4, nbear = 2, 
  #'                                            ncoy = 1, nelk = 3, nmoose = 2, nwtd = 2, nharv = 5, nfor = 4, nwsi = 3)
  #' 
  #'                                
  #' # save(data_JAGS_bundle_topinter, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_topinter.RData")
  #' # save(data_JAGS_bundle_topexploit, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_topexploit.RData")
  #' # save(data_JAGS_bundle_bottominter, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_bottominter.RData")
  #' # save(data_JAGS_bundle_bottomexploit, file = "./Outputs/SEM/JAGS_data_bundle/data_JAGS_bundle_bottomexploit.RData")
  #' 
  #' #'  Generate initial values for each parameter (random node)
  #' generate_inits <- function(nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, nfor, nwsi) {
  #'   
  #'   #'  Generate random values for each species-specific beta (nwolf, nlion, etc.
  #'   #'  based on number of species-specific betas to be estimated)
  #'   list(
  #'     wolf.t_1 = runif(nwolf, -0.5, 0.5),  # consider -1, 1
  #'     lion.t_1 = runif(nlion, -0.5, 0.5),
  #'     bear.t_1 = runif(nbear, -0.5, 0.5),
  #'     coy.t_1 = runif(ncoy, -0.5, 0.5),
  #'     elk.t_1 = runif(nelk, -0.5, 0.5),
  #'     moose.t_1 = runif(nmoose, -0.5, 0.5),
  #'     wtd.t_1 = runif(nwtd, -0.5, 0.5),
  #'     wolfharvest.t_1 = runif(nharv, -0.5, 0.5),
  #'     lionharvest.t_1 = runif(nharv, -0.5, 0.5),
  #'     bearharvest.t_1 = runif(nharv, -0.5, 0.5),
  #'     elkharvest.t_1 = runif(nharv, -0.5, 0.5),
  #'     mooseharvest.t_1 = runif(nharv, -0.5, 0.5),
  #'     deerharvest.t_1 = runif(nharv, -0.5, 0.5),
  #'     wsi.t_1 = runif(nwsi, -0.5, 0.5),
  #'     forest.t_1 = runif(nfor, -0.5, 0.5),
  #'     road.t_1 = runif(nwsi, -0.5, 0.5),
  #'     public.t_1 = runif(nwsi, -0.5, 0.5) # doesn't matter that using nwsi or other b/c all = 1
  #'     #' #'  Fix random number generator and seed for every run of this function
  #'     #' .RNG.name = "base::Wichmann-Hill",
  #'     #' .RNG.seed = 182  
  #'     #' #'  Setting RNG seed leads to the same random number stream during
  #'     #' #'  adaptation and sampling, causing different inits to be rapidly erased
  #'     #' #'  and all chains to follow the same deterministic path. This does not 
  #'     #' #'  appear to happen with the cross lag model though.
  #'   )
  #' }
  #' #'  Define number of chains
  #' num.chains <- 3
  #' #'  Create empty lists
  #' initsList_ar1 <- initsList_topdown <- initsList_bottomup <- vector('list', num.chains) 
  #' initsList_topinter <- initsList_topexploit <- initsList_bottominter <- initsList_topbottom <- vector('list', num.chains)
  #' initsList_topinter_reduced <- initsList_bottominter_reduced <- vector('list', num.chains)
  #' #'  Setting seed for reproducibility
  #' set.seed(9983)
  #' #'  Loop through generate_inits function 3 times (1 for each chain) 
  #' for(i in 1:num.chains){
  #'   initsList_ar1[[i]] <- generate_inits(nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 1, 
  #'                                    nmoose = 1, nwtd = 1, nharv = 0, nfor = 0, nwsi = 0)
  #'   initsList_topdown[[i]] <- generate_inits(nwolf = 4, nlion = 3, nbear = 3, ncoy = 2, nelk = 1, 
  #'                                            nmoose = 1, nwtd = 1, nharv = 6, nfor = 0, nwsi = 0)
  #'   initsList_bottomup[[i]] <- generate_inits(nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 5, 
  #'                                             nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  #'   initsList_topinter[[i]] <- generate_inits(nwolf = 7, nlion = 4, nbear = 5, ncoy = 2, nelk = 1, 
  #'                                             nmoose = 1, nwtd = 1, nharv = 6, nfor = 0, nwsi = 0)
  #'   initsList_topinter_reduced[[i]] <- generate_inits(nwolf = 6, nlion = 4, nbear = 2, ncoy = 1, nelk = 1, 
  #'                                             nmoose = 1, nwtd = 1, nharv = 1, nfor = 0, nwsi = 0)
  #'   # initsList_topexploit[[i]] <- generate_inits(nwolf = 4, nlion = 3, nbear = 1, ncoy = 1, nelk = 1, 
  #'   #                                           nmoose = 1, nwtd = 1, nharv = 3, nfor = 0, nwsi = 0)
  #'   initsList_bottominter[[i]] <- generate_inits(nwolf = 4, nlion = 1, nbear = 1, ncoy = 1, nelk = 5, 
  #'                                                nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  #'   initsList_bottominter_reduced[[i]] <- generate_inits(nwolf = 4, nlion = 1, nbear = 1, ncoy = 1, nelk = 2, 
  #'                                                        nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  #'   initsList_topbottom[[i]] <- generate_inits(nwolf = 7, nlion = 4, nbear = 2, ncoy = 1, nelk = 3, 
  #'                                              nmoose = 2, nwtd = 2, nharv = 5, nfor = 4, nwsi = 3)
  #' }
  
  #'  Parameters monitored
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.lion", "beta.bear", 
              "beta.coy", "beta.elk", "beta.moose", "beta.wtd", "beta.harvest", 
              "beta.wsi","beta.forest",  "sigma.spp", "lion.latent", "wolf.latent", 
              "bear.latent", "coy.latent", "elk.latent", "moose.latent", "wtd.latent") # "sigma.spp.tmin1", "sigma.cluster", "cluster.randeff" 
   
  #'  MCMC settings
  nc <- 3
  ni <- 100000
  nb <- 50000
  nt <- 10
  na <- 5000
  
  
  #'  ---------------------------------
  ####  Call JAGS & Fit Bayesian SEMs  ####
  #'  ---------------------------------
  #####  Top-down model  #####
  #'  Call bundle_data function (Format_RNmodel_Posteriors_for_SEM.R) to bundle
  #'  input data for JAGS
  data_JAGS_bundle_topdown <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                         dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                         covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                         covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                         nwolf = 3, nlion = 2, nbear = 1, ncoy = 1, nelk = 1, 
                                         nmoose = 1, nwtd = 1, nharv = 5, nfor = 0, nwsi = 0)
                                          
  #'  Define number of chains for initial values to be parallelized
  num.chains <- 3
  #'  Create empty fector to hold initial values
  initsList_topdown <- vector('list', num.chains) 
  #'  Call generate_inits function (in Format_RNmodel_Posteriors_for_SEM.R) and draw
  #'  random starting values
  for(i in 1:num.chains) {
    initsList_topdown[[i]] <- generate_inits(nwolf = 3, nlion = 2, nbear = 1, ncoy = 1, nelk = 1, nmoose = 1, 
                                             nwtd = 1, nharv = 5, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  #'  Source top-down SEM script
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown.R")
  #'  Fit model, review output and assess convergence
  start.time = Sys.time()
  SEM_topdown <- jagsUI::jags(data_JAGS_bundle_topdown, inits = initsList_topdown, params, 
                              "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown$summary)  
  which(SEM_topdown$summary[,"Rhat"] < 0.9)
  which(SEM_topdown$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_topdown$samples)
  #'  Save model output
  save(SEM_topdown, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_", Sys.Date(), ".RData"))
  
  #####  Top-down, exploitation d-Sep updated  #####
  #'  Update JAGS inputs
  data_JAGS_bundle_top_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                           dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                           covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                           covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                           nwolf = 4, nlion = 2, nbear = 2, ncoy = 2, nelk = 3, 
                                           nmoose = 2, nwtd = 7, nharv = 5, nfor = 0, nwsi = 0)
  initsList_topdown_final <- vector('list', num.chains)
  for(i in 1:num.chains){
    initsList_topdown_final[[i]] <- generate_inits(nwolf = 4, nlion = 2, nbear = 2, ncoy = 2, nelk = 3, nmoose = 2, 
                                                   nwtd = 7, nharv = 5, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_final.R")
  
  #'  Updated list of parameters to follow (includes derived parameters for indirect effects now)
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk",
              "beta.moose", "beta.wtd", "beta.harvest", "sigma.spp", "indirect.lion.wtd.lion", 
              "indirect.lion.wtd.wolf", "indirect.lion.wtd.bear", "indirect.wolf.wtd.lion", 
              "indirect.wolf.wtd.wolf", "indirect.wolf.wtd.bear", "indirect.coy.wtd.lion",  
              "indirect.coy.wtd.wolf",  "indirect.coy.wtd.bear", "indirect.wolf.elk.lion", 
              "indirect.wolf.elk.wolf", "indirect.lion.elk.lion", "indirect.lion.elk.wolf",
              "indirect.bear.elk.wolf", "indirect.bear.elk.lion", "indirect.wolf.moose.wolf",
              "indirect.lion.wtdlag.lion", "indirect.lion.wtdlag.wolf", "indirect.lion.wtdlag.coy",
              "indirect.wolf.wtdlag.lion", "indirect.wolf.wtdlag.wolf", "indirect.wolf.wtdlag.coy",
              "indirect.coy.wtdlag.lion",  "indirect.coy.wtdlag.wolf",  "indirect.coy.wtdlag.coy",
              "indirect.wolf.elk.self", "indirect.lion.elk.self", "indirect.bear.elk.self",
              "indirect.wolf.moose.self", "indirect.lion.wtd.self", "indirect.wolf.wtd.self",
              "indirect.coy.wtd.self", "indirect.elk.lion.elk", "indirect.elk.lion.wtd",
              "indirect.wtd.lion.elk", "indirect.wtd.lion.wtd", "indirect.moose.wolf.elk", 
              "indirect.moose.wolf.moose", "indirect.moose.wolf.wtd", "indirect.elk.wolf.elk",   
              "indirect.elk.wolf.moose",   "indirect.elk.wolf.wtd", "indirect.wtd.wolf.elk",   
              "indirect.wtd.wolf.moose",   "indirect.wtd.wolf.wtd", "indirect.wtd.lion.elk.v2", 
              "indirect.wtd.lion.wtd.v2", "indirect.wtd.wolf.elk.v2", "indirect.wtd.wolf.moose.v2", 
              "indirect.wtd.wolf.wtd.v2", "indirect.wtd.bear.elk.v2")
  
  start.time = Sys.time()
  SEM_topdown_final <- jagsUI::jags(data_JAGS_bundle_top_final, inits = initsList_topdown_final, params, 
                      "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_final.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_final$summary)  
  which(SEM_topdown_final$summary[,"Rhat"] < 0.9)
  which(SEM_topdown_final$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_topdown_final$samples)
  save(SEM_topdown_final, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_exploitation_final_", Sys.Date(), ".RData"))
  
  
  #####  Top-down, interference model  #####
  data_JAGS_bundle_topdown_inter <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                         dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                         covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                         covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                         nwolf = 6, nlion = 4, nbear = 1, ncoy = 1, nelk = 1, 
                                         nmoose = 1, nwtd = 1, nharv = 3, nfor = 0, nwsi = 0)
  num.chains <- 3
  initsList_topdown_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_topdown_inter[[i]] <- generate_inits(nwolf = 6, nlion = 4, nbear = 1, ncoy = 1, nelk = 1, nmoose = 1, 
                                             nwtd = 1, nharv = 3, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter.R")
  start.time = Sys.time()
  SEM_topdown_inter <- jagsUI::jags(data_JAGS_bundle_topdown_inter, inits = initsList_topdown_inter, params, 
                                            "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter.txt",
                                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_inter$summary)  
  which(SEM_topdown_inter$summary[,"Rhat"] < 0.9)    
  which(SEM_topdown_inter$summary[,"Rhat"] > 1.1)    
  mcmcplot(SEM_topdown_inter$samples)
  save(SEM_topdown_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_inter_", Sys.Date(), ".RData"))
  
  #####  Top-down, interference, d-Sep updated  #####
  data_JAGS_bundle_topdown_inter_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                     dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                                     covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                     covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                     nwolf = 6, nlion = 3, nbear = 3, ncoy = 2, nelk = 1, 
                                                     nmoose = 1, nwtd = 1, nharv = 3, nfor = 0, nwsi = 0)
                                               
  num.chains <- 3
  initsList_topdown_inter_final <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_topdown_inter_final[[i]] <- generate_inits(nwolf = 6, nlion = 3, nbear = 3, ncoy = 2, nelk = 1, nmoose = 1, 
                                                         nwtd = 1, nharv = 3, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_inter_final.R")
  
  #'  Updated list of parameters to follow (includes derived parameters for indirect effects now)
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk",
              "beta.moose", "beta.wtd", "beta.harvest", "sigma.spp", "indirect.wolf.lion.elk", 
              "indirect.wolf.lion.wtd", "indirect.wolf.lion.coy", "indirect.wolf.bear.wtd",
              "indirect.wolf.coy.wtd", "indirect.lion.coy.wtd", "indirect.wolf.bear.wolf",
              "indirect.bear.wolf.bear", "indirect.bear.wolf.lion", "indirect.bear.wolf.coy",
              "indirect.wolf.bear.self", "indirect.wolf.coy.self", "indirect.lion.coy.self", 
              "indirect.bear.wolf.self", "indirect.wolf.elk.self", "indirect.lion.elk.self", 
              "indirect.wolf.moose.self", "indirect.lion.wtd.self", "indirect.coy.wtd.self", 
              "indirect.bear.wtd.self")
  
  start.time = Sys.time()
  SEM_topdown_inter_final <- jagsUI::jags(data_JAGS_bundle_topdown_inter_final, inits = initsList_topdown_inter_final, params, 
                                          "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_final.txt",
                                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                          n.burnin = nb, parallel = TRUE)
                                    
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_topdown_inter_final$summary)  
  which(SEM_topdown_inter_final$summary[,"Rhat"] < 0.9)    
  which(SEM_topdown_inter_final$summary[,"Rhat"] > 1.1)    
  mcmcplot(SEM_topdown_inter_final$samples)
  save(SEM_topdown_inter_final, file = paste0("./Outputs/SEM/JAGS_out/SEM_topdown_inter_final_", Sys.Date(), ".RData"))
  
  
  
  #####  Bottom-up model  ##### 
  data_JAGS_bundle_bottomup <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                        dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                        covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                        covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                        nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, 
                                        nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  num.chains <- 3
  initsList_bottomup <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup[[i]] <- generate_inits(nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, nmoose = 2, 
                                              nwtd = 3, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup.R")
  start.time = Sys.time()
  SEM_bottomup <- jagsUI::jags(data_JAGS_bundle_bottomup, inits = initsList_bottomup, params, 
                               "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup$summary)
  which(SEM_bottomup$summary[,"Rhat"] < 0.9)
  which(SEM_bottomup$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup$samples)
  save(SEM_bottomup, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_", Sys.Date(), ".RData"))
  
  
  #####  Bottom-up, exploitation d-Sep updated  ##### 
  data_JAGS_bundle_bottomup_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                          dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                          covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                          covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                          nwolf = 1, nlion = 1, nbear = 3, ncoy = 2, nelk = 6, 
                                          nmoose = 3, nwtd = 7, nharv = 0, nfor = 4, nwsi = 3)
  num.chains <- 3
  initsList_bottomup_final <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup_final[[i]] <- generate_inits(nwolf = 1, nlion = 1, nbear = 3, ncoy = 2, nelk = 6, nmoose = 3, 
                                              nwtd = 7, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_final.R")
  
  #'  Updated list of parameters to follow (includes derived parameters for indirect effects now)
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.bear", "beta.coy",
              "beta.elk", "beta.moose", "beta.wtd", "beta.forest", "beta.wsi", "sigma.spp",
              "indirect.bear.elk.lion", "indirect.bear.elk.wolf", "indirect.coy.wtd.lion", 
              "indirect.coy.wtd.wolf", "indirect.coy.wtd.bear", "indirect.coy.wtd.coy",
              "indirect.bear.wtd.lion", "indirect.bear.wtd.wolf", "indirect.bear.wtd.bear", 
              "indirect.bear.wtd.coy", "indirect.bear.elklag.lion", "indirect.bear.elklag.wolf", 
              "indirect.bear.elklag.bear", "indirect.coy.wtdlag.lion", "indirect.coy.wtdlag.coy",
              "indirect.bear.wtdlag.lion", "indirect.bear.wtdlag.coy", "indirect.bear.elk.self", 
              "indirect.bear.wtd.self", "indirect.coy.wtd.self", "indirect.elk.bear.elk", 
              "indirect.elk.bear.wtd")
  
  start.time = Sys.time()
  SEM_bottomup_final <- jagsUI::jags(data_JAGS_bundle_bottomup_final, inits = initsList_bottomup_final, params, 
                               "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_final.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_final$summary)
  which(SEM_bottomup_final$summary[,"Rhat"] < 0.9)
  which(SEM_bottomup_final$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_final$samples)
  save(SEM_bottomup_final, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_exploitation_final", Sys.Date(), ".RData"))
  
  
  #####  Bottom-up, interference model  #####
  data_JAGS_bundle_bottom_inter <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                     dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                     covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                     covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                     nwolf = 3, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, 
                                     nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  num.chains <- 3
  initsList_bottom_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottom_inter[[i]] <- generate_inits(nwolf = 3, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, nmoose = 2, 
                                                    nwtd = 3, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_inter.R")
  start.time = Sys.time()
  SEM_bottomup_inter <- jagsUI::jags(data_JAGS_bundle_bottom_inter, inits = initsList_bottom_inter, params,
                                     "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_inter.txt",
                                     n.adapt = na, n.chains = nc, n.thin = nt, 
                                     n.iter = ni, n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_inter$summary)
  which(SEM_bottomup_inter$summary[,"Rhat"] < 0.9)
  which(SEM_bottomup_inter$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_inter$samples)
  save(SEM_bottomup_inter, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_inter_", Sys.Date(), ".RData"))
  
  
  #####  Bottom-up, interference d-Sep updated  #####
  data_JAGS_bundle_bottom_inter_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                    dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                                    covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                    covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                    nwolf = 3, nlion = 0, nbear = 5, ncoy = 1, nelk = 6, 
                                                    nmoose = 3, nwtd = 6, nharv = 0, nfor = 4, nwsi = 3)
                                             
  num.chains <- 3
  initsList_bottom_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottom_inter[[i]] <- generate_inits(nwolf = 3, nlion = 0, nbear = 5, ncoy = 1, nelk = 6, nmoose = 3, 
                                                 nwtd = 6, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_inter_final.R")
  
  #'  Updated list of parameters to follow (includes derived parameters for indirect effects now)
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.bear", "beta.coy",
              "beta.elk", "beta.moose", "beta.wtd", "beta.forest", "beta.wsi", "sigma.spp",
              "indirect.elk.wolf.bear", "indirect.moose.wolf.bear.v1", "indirect.moose.wolf.bear.v2",
              "indirect.elk.wolf.coy",  "indirect.moose.wolf.coy.v1",  "indirect.moose.wolf.coy.v2",
              "indirect.elk.bear.lion", "indirect.wtd.bear.lion", "indirect.elk.bear.wolf.v1", 
              "indirect.wtd.bear.wolf.v1", "indirect.elk.bear.wolf.v2", "indirect.wtd.bear.wolf.v2",
              "indirect.elk.bear.coy", "indirect.wtd.bear.coy", "indirect.wolf.bear.lion",
              "indirect.wolf.bear.wolf.v1", "indirect.wolf.bear.wolf.v2", "indirect.wolf.bear.coy",
              "indirect.bear.wolf.bear.v1", "indirect.bear.wolf.bear.v2", "indirect.bear.wolf.coy.v1",  
              "indirect.bear.wolf.coy.v2", "indirect.elk.wolf.self", "indirect.moose.wolf.self",
              "indirect.elk.bear.self", "indirect.wtd.bear.self", "indirect.wtd.coy.self",
              "indirect.bear.wolf.self", "indirect.wolf.bear.self", "indirect.wolf.coy.self",  
              "indirect.bear.coy.self")
  
  start.time = Sys.time()
  SEM_bottomup_inter_final <- jagsUI::jags(data_JAGS_bundle_bottom_inter_final, inits = initsList_bottom_inter, params,
                                           "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_inter_final.txt",
                                           n.adapt = na, n.chains = nc, n.thin = nt, 
                                           n.iter = ni, n.burnin = nb, parallel = TRUE)
                                     
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_inter_final$summary)
  which(SEM_bottomup_inter_final$summary[,"Rhat"] < 0.9)
  which(SEM_bottomup_inter_final$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_inter_final$samples)
  save(SEM_bottomup_inter_final, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_inter_final_", Sys.Date(), ".RData"))
  
  
  
  
  
  
  
  