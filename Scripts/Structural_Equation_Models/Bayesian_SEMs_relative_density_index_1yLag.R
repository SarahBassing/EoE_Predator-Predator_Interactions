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
  data_JAGS_bundle_topdown_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                               dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                               covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                               covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                               nwolf = 4, nlion = 2, nbear = 2, ncoy = 2, nelk = 3, 
                                               nmoose = 2, nwtd = 7, nharv = 5, nfor = 0, nwsi = 0)
  num.chains <- 3                                         
  initsList_topdown_final <- vector('list', num.chains)
  for(i in 1:num.chains){
    initsList_topdown_final[[i]] <- generate_inits(nwolf = 4, nlion = 2, nbear = 2, ncoy = 2, nelk = 3, nmoose = 2, 
                                                   nwtd = 7, nharv = 5, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_topdown_final.R")
  
  #'  Latent parameters to monitor if needed (necessary when rerunning this model for final d-Sep and Fisher's C)
  latent_params <- c("lion.latent", "wolf.latent", "bear.latent", "coy.latent", "elk.latent", "moose.latent", "wtd.latent")
  
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
              "indirect.wtd.wolf.wtd.v2", "indirect.wtd.bear.elk.v2",latent_params) 
  
  start.time = Sys.time()
  SEM_topdown_final <- jagsUI::jags(data_JAGS_bundle_topdown_final, inits = initsList_topdown_final, params, 
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
  
  #'  Latent parameters to monitor if needed (necessary when rerunning this model for final d-Sep and Fisher's C)
  latent_params <- c("lion.latent", "wolf.latent", "bear.latent", "coy.latent", "elk.latent", "moose.latent", "wtd.latent")
  
  #'  Updated list of parameters to follow (includes derived parameters for indirect effects now)
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk",
              "beta.moose", "beta.wtd", "beta.harvest", "sigma.spp", "indirect.wolf.lion.elk", 
              "indirect.wolf.lion.wtd", "indirect.wolf.lion.coy", "indirect.wolf.bear.wtd",
              "indirect.wolf.coy.wtd", "indirect.lion.coy.wtd", "indirect.wolf.bear.wolf",
              "indirect.bear.wolf.bear", "indirect.bear.wolf.lion", "indirect.bear.wolf.coy",
              "indirect.wolf.bear.self", "indirect.wolf.coy.self", "indirect.lion.coy.self", 
              "indirect.bear.wolf.self", "indirect.wolf.elk.self", "indirect.lion.elk.self", 
              "indirect.wolf.moose.self", "indirect.lion.wtd.self", "indirect.coy.wtd.self", 
              "indirect.bear.wtd.self", latent_params) 
  
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
  
  #'  Latent parameters to monitor if needed (necessary when rerunning this model for final d-Sep and Fisher's C)
  latent_params <- c("lion.latent", "wolf.latent", "bear.latent", "coy.latent", "elk.latent", "moose.latent", "wtd.latent")
  
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
              "indirect.elk.bear.wtd", latent_params) 
  
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
  save(SEM_bottomup_final, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_exploitation_final_", Sys.Date(), ".RData"))
  
  
  #####  Bottom-up, interference model  #####
  data_JAGS_bundle_bottomup_inter <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                     dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                     covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                     covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                     nwolf = 3, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, 
                                     nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  num.chains <- 3
  initsList_bottomup_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup_inter[[i]] <- generate_inits(nwolf = 3, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, nmoose = 2, 
                                                    nwtd = 3, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_inter.R")
  start.time = Sys.time()
  SEM_bottomup_inter <- jagsUI::jags(data_JAGS_bundle_bottomup_inter, inits = initsList_bottomup_inter, params,
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
  data_JAGS_bundle_bottomup_inter_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                    dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                                    covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                    covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                    nwolf = 3, nlion = 0, nbear = 5, ncoy = 1, nelk = 6, 
                                                    nmoose = 3, nwtd = 6, nharv = 0, nfor = 4, nwsi = 3)
                                             
  num.chains <- 3
  initsList_bottomup_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup_inter[[i]] <- generate_inits(nwolf = 3, nlion = 0, nbear = 5, ncoy = 1, nelk = 6, nmoose = 3, 
                                                    nwtd = 6, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  }
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_bottomup_inter_final.R")
  
  #'  Latent parameters to monitor if needed (necessary when rerunning this model for final d-Sep and Fisher's C)
  latent_params <- c("lion.latent", "wolf.latent", "bear.latent", "coy.latent", "elk.latent", "moose.latent", "wtd.latent")
  
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
              "indirect.bear.coy.self", latent_params) 
  
  start.time = Sys.time()
  SEM_bottomup_inter_final <- jagsUI::jags(data_JAGS_bundle_bottomup_inter_final, inits = initsList_bottomup_inter, params,
                                           "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_inter_final.txt",
                                           n.adapt = na, n.chains = nc, n.thin = nt, 
                                           n.iter = ni, n.burnin = nb, parallel = TRUE)
                                     
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(SEM_bottomup_inter_final$summary)
  which(SEM_bottomup_inter_final$summary[,"Rhat"] < 0.9)
  which(SEM_bottomup_inter_final$summary[,"Rhat"] > 1.1)
  mcmcplot(SEM_bottomup_inter_final$samples)
  save(SEM_bottomup_inter_final, file = paste0("./Outputs/SEM/JAGS_out/SEM_bottomup_inter_final_", Sys.Date(), ".RData"))
  
  
  #'  ---------------------------
  ####  LOSO CV model selection  ####
  #'  ---------------------------
  #'  Withhold one site (all species across years) per fold and refit the final
  #'  model with the training data (all sites not withheld). Then compare refit
  #'  posterior to original estimate and calculate the likelihood of the original
  #'  value if each draw of the refit model's estimate was correct.
  #'  ---------------------------
  #'  install.packages("loo")
  library(loo)
  library(future.apply)
  plan(multisession, workers = parallel::detectCores() - 1)
  
  #'  Function to leave-one-site-out (LOSO) cross validation (CV)
  run_loso_cv <- function(mod, dat, spp_names, nSites, nYear, n.chains,
                          n.adapt, n.burnin, n.iter, n.thin, out_dir) {
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    #'  Store unmodified spp.hat and spp.sigma_hat arrays 
    true_hat <- setNames(lapply(spp_names, function(sp) dat[[paste0(sp, ".hat")]]), spp_names)
    true_sigma_hat <- setNames(lapply(spp_names, function(sp) dat[[paste0(sp, ".sigma_hat")]]), spp_names)
    
    #'  Build a list of parameter names for JAGS to monitor
    latent_params <- paste0(spp_names, ".latent")
    
    # site_results <- vector("list", nSites) 
    #' #'  Each loop holds one site out as test data while retaining the other sites as training data
    #' for(s in 1:nSites) {
    #'   message(sprintf("LOSO fold %d/%d (site %d held out)", s, nSites, s))
    
    site_results <- future_lapply(1:nSites, function(s) {
      #'  Copy the full data set
      fold_data <- dat
      
      #'  Then create this fold's testing data by setting site s's .hat to NA for 
      #'  every species and year
      #'  Note: covariate data from site s remain untouched - allows us to test
      #'  "given what we know about this site's habitat/harvest, does the model 
      #'  correctly predict each species' RDI?
      for(sp in spp_names) {
        fold_data[[paste0(sp, ".hat")]][s,] <- NA
      }
      
      #'  Refit the model from scratch with the modified data
      fit <- jagsUI::jags(data = fold_data, inits = NULL, parameters.to.save = latent_params,
                          model.file = mod, n.chains = n.chains, n.adapt = n.adapt,
                          n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                          parallel = FALSE, verbose = FALSE)
      
      #'  Compute posterior predictive log density of TRUE observed value for 
      #'  each species/year at the held-out site
      #'  log(mean_over_draws( p(y_true | draw's latent state))) 
      #'  via log-sum-exp for numerical stability
      site_loglik <- 0
      # site_detail <- list()
      #'  Grab the refit model's full posterior for the held-out site for a given species and year
      for(sp in spp_names) {
        latent_node <- fit$sims.list[[paste0(sp, ".latent")]]  # [n_draws, nSites, nYear]
        
        for(y in 1:nYear) {
          #'  Check what the real "observed" value was from the saved copy of the OG data
          y_true <- true_hat[[sp]][s,y]
          #'  If the original data was missing for this site and year, skip (nothing to score)
          if(is.na(y_true)) next              
          
          #'  Compare every posterior draw from refit model to the true "observed" value
          #'  and ask how likely was the true observation if this particular draw's
          #'  latent-state estimate was correct
          #'  Grab the true "observed" sigma
          sigma_true <- true_sigma_hat[[sp]][s,y]
          #'  Grab every posterior draw of the refit model
          draws <- latent_node[,s,y]
          
          #'  Compute the log-likelihood per draw using the known measure of uncertainty 
          #'  from RN model as the spread. This returns the density of how each possible
          #'  version of reality (draw) in the posterior predicted what actually happened
          logdens_draws <- dnorm(y_true, mean = draws, sd = sigma_true, log = TRUE)
          #'  Generate an overall predictive score - 
          #'  Calculate the posterior predictive log density (lpd) by subtracting the max value,
          #'  exponentiating, averaging, and then adding the max back on the log scale. This
          #'  keeps everything numerically stable.
          m <- max(logdens_draws)
          site_loglik <- site_loglik + m + log(mean(exp(logdens_draws - m)))  # log-sum-exp
          
          #' #' Add predictive score to running total for log-likelihood of whole 
          #' #' held-out site accumulating across every species and score-able year
          #' site_loglik <- site_loglik + lpd
          #' site_detail[[paste0(sp, "_yr", y)]] <- lpd
        }
      }
      
      #'  Save score for every species/year for site s
      # site_results[[s]] <- list(site = s, elpd = site_loglik, detail = site_detail)
      # saveRDS(site_results[[s]], file.path(out_dir, sprintf("site_%03d.rds", s)))
      saveRDS(list(site = s, elpd = site_loglik), file.path (out_dir, sprintf("site_%03d.rds", s)))
      list(site = s, elpd = site_loglik)
      
    }, future.seed = TRUE)
    
    #'  Sum every site's score into an overall score for the entire model: expected log predictive density (elpd)
    #'  Higher elpd means indicate the model, on average, made better honest predictions
    #'  about sites it did not see. se_elpd estimates how much that total could plausibly
    #'  vary, treating each score as on independent data point
    elpd_per_site <- sapply(site_results, function(x) x$elpd)
    
    list(elpd_per_site = elpd_per_site,
         total_elpd = sum(elpd_per_site),
         se_elpd = sd(elpd_per_site) * sqrt(nSites), # treats sites as the exchangeable unit, similar to loo's SE
         site_results = site_results)
    }
  
  #'  Run LOSO CV for each model
  start.time = Sys.time()
  loso_topdown_exploit <- run_loso_cv(mod = "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_final.txt", dat = data_JAGS_bundle_topdown_final,
                                      spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd"), nSites = 23, 
                                      nYear = 4, n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt,
                                      out_dir = "./Outputs/SEM/LOSO_CV/topdown_exploit")
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(loso_topdown_exploit, file = "./Outputs/SEM/LOSO_CV/LOSO_CV_topdown_exploit.RData")
  
  start.time = Sys.time()
  loso_topdown_inter <- run_loso_cv(mod = "./Outputs/SEM/JAGS_out/JAGS_SEM_topdown_inter_final.txt", dat = data_JAGS_bundle_topdown_inter_final,
                                      spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd"), nSites = 23,  
                                      nYear = 4, n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt,
                                      out_dir = "./Outputs/SEM/LOSO_CV/topdown_inter")
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(loso_topdown_inter, file = "./Outputs/SEM/LOSO_CV/LOSO_CV_topdown_inter.RData")
  
  start.time = Sys.time()
  loso_bottomup_exploit <- run_loso_cv(mod = "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_final.txt", dat = data_JAGS_bundle_bottomup_final,
                                      spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd"), nSites = 23,  
                                      nYear = 4, n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt,
                                      out_dir = )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(loso_bottomup_exploit, file = "./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_exploit.RData")
  
  start.time = Sys.time()
  loso_bottomup_inter <- run_loso_cv(mod = "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_inter_final.txt", dat = data_JAGS_bundle_bottomup_inter_final,
                                      spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd"), nSites = 23,  
                                      nYear = 4, n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt,
                                      out_dir = "./Outputs/SEM/LOSO_CV/bottomup_inter")
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(loso_bottomup_inter, file = "./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_inter.RData")
  
  #'  Read in LOSO_CV results
  load("./Outputs/SEM/LOSO_CV/LOSO_CV_topdown_exploit.RData")
  load("./Outputs/SEM/LOSO_CV/LOSO_CV_topdown_inter.RData")
  load("./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_exploit.RData")
  load("./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_inter.RData")
  
  #'  Compare LOSO CV across models
  compare_loso <- data.frame(model = c("topdown_exploitative", "topdown_interference", 
                                       "bottomup_exploitative", "bottomup_interference"),
                             elpd = c(loso_topdown_exploit$total_elpd, loso_topdown_inter$total_elpd, 
                                      loso_bottomup_exploit$total_elpd, loso_bottomup_inter$total_elpd),
                             se = c(loso_topdown_exploit$se_elpd, loso_topdown_inter$se_elpd, 
                                    loso_bottomup_exploit$se_elpd, loso_bottomup_inter$se_elpd))
  
  #'  Calculate difference between largest elpd ("best" supported model) and other models
  compare_loso$elpd_diff <- compare_loso$elpd - max(compare_loso$elpd)
  compare_loso$SEx2 <- compare_loso$se * 2
  compare_loso[order(-compare_loso$elpd), ]
  #' Differences in elpd are worth treating as meaningful when |elpd_diff| is roughly >=2x its se
  #' In other words, the models are distinguishable when |elpd_diff| is meaningfully larger than its se
  
  #'  Paired comparisons to evaluate site-by-site differences between models 
  #'  Some sites were likely harder to estimate than others and this problem likely
  #'  arose for each model so comparing paired models helps differentiate them
  compare_paired_mods <- function(elpd_per_site_A, elpd_per_site_B) {
    #'  Site-level differences between two models (A - B)
    diff_vec <- elpd_per_site_A - elpd_per_site_B
    elpd_diff <- sum(diff_vec)
    se_diff <- sd(diff_vec) * sqrt(length(diff_vec))
    list(elpd_diff = elpd_diff, se_diff = se_diff, ratio = abs(elpd_diff) / se_diff) # >= ~2 suggests meaningful difference
  }
  compare_paired_mods(loso_bottomup_exploit$elpd_per_site, loso_topdown_exploit$elpd_per_site)
  compare_paired_mods(loso_bottomup_exploit$elpd_per_site, loso_topdown_inter$elpd_per_site)   # suggests a meaningful difference from bottomup_exploit
  compare_paired_mods(loso_bottomup_exploit$elpd_per_site, loso_bottomup_inter$elpd_per_site)
  
  #'  bottomup_exploit and topdown_exploit are essentially indistinguishable 
  #'  bottomup_inter is not meaningfully different from bottomup_exploit
  #'  topdown_inter is distinctly different and less supported than bottomup_exploit
  #'  topdown_exploit and bottomup_inter are statistically indistinguishable alternatives 
  #'  to bottomup_exploit given your data
  
  
  #' ------------------
  #####  LOO and WAIC  #####
  #' ------------------
  #' #'  Extract log-likelihood samples
  #' build_loglik_matrix <- function(mod, dat, spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd")) {
  #'   #'  Create empty list to hold matrix
  #'   mat_list <- list()
  #'   
  #'   for(i in spp_names) {
  #'     #'  extract log-likelihood samples that match spp name
  #'     loglik_node <- mod$sims.list[[paste0("loglik.", i)]]   # dims [n_draws, nSites, nYear]
  #'     
  #'     #'  Extract "raw" observation data that match species name
  #'     hat_array <- dat[[paste0(i, ".hat")]]  # dims [nSites, nYear], NA where missing values
  #'     #'  Create vector that returns TRUE for observations that are not NA
  #'     keep <- !is.na(as.vector(hat_array))
  #'     
  #'     #'  Create matrix with "raw" observed data, excluding NAs
  #'     #'  Count number of draws in loglike_node
  #'     n_draws <- dim(loglik_node)[1]
  #'     #'  Flatten 3D array of loglik_node into 2D matrix where draws stay as rows
  #'     flat <- matrix(loglik_node, nrow = n_draws)
  #'     #'  Remove observations that were flagged as FALSE (NAs in input data)
  #'     flat <- flat[, keep, drop = FALSE]
  #'     
  #'     colnames(flat) <- paste0(i, ".", which(keep))
  #'     mat_list[[i]] <- flat
  #'   }
  #'   
  #'   #'  Return [n_draws x total_real_observations] matrix
  #'   do.call(cbind, mat_list)
  #'   
  #' }
  #' 
  #' #'  Leave-one-out (LOO)
  #' run_loo <- function(mod, dat) {
  #'   #'  Call build_loglik_matrix() function to format log-likelihood data
  #'   loglik_matrix <- build_loglik_matrix(mod, dat)
  #'   
  #'   #'  Generate LOO and WAIC values
  #'   loo_out <- loo(loglik_matrix)
  #'   waic_out <- waic(loglik_matrix)
  #'   
  #'   #'  Check Pareto-k diagnostics for LOO
  #'   #'  NOTE: values > 0.7 flag observations where importance-sampling approximation
  #'   #'  is unreliable (Vehtari, Gelman & Gabry 2017). A few flagged points is normal 
  #'   #'  but a large proportion suggests PSIS-LOO is struggling and results should
  #'   #'  be interpreted with caution.
  #'   n_bad_k <- sum(loo_out$diagnostics$pareto_k > 0.7)
  #'   if(n_bad_k > 0) {
  #'     message(sprintf(
  #'       "%d / %d observations have Pareto k> 0.7 - PSIS-LOO approx. may be unreliable for these. Check loo_result$diagnostics$pareto_k.",
  #'       n_bad_k, length(loo_out$diagnostics$pareto_k)
  #'     ))
  #'   }
  #'   
  #'   loo_list <- list(loo = loo_out, waic = waic_out, loglik_matrix = loglik_matrix)
  #'   # loo_only <- loo_out
  #'   return(loo_list)
  #'   
  #' }
  #' 
  #' #'  Run loo function for each model
  #' loo_topdown_exploit <- run_loo(SEM_topdown_final, data_JAGS_bundle_topdown_final)
  #' loo_topdown_inter <- run_loo(SEM_topdown_inter_final, data_JAGS_bundle_topdown_inter_final)
  #' loo_bottomup_exploit <- run_loo(SEM_bottomup_final, data_JAGS_bundle_bottomup_final)
  #' loo_bottomup_inter <- run_loo(SEM_bottomup_inter_final, data_JAGS_bundle_bottom_inter_final)
  #' 
  #' #'  Compare loo (and WAIC) across models
  #' #'  NOTE: loo_compare ranks models by expected log predictive density (ELPD). 
  #' #'  Top row is the best supported model and elpd_diff / se_diff indicates how 
  #' #'  many standard errors separate each model from the "top" model. Typically, 
  #' #'  a |elpd_diff| less than x2 its se_diff is not clearly different from top model.
  #' loo_compare(loo_topdown_exploit[[1]], loo_topdown_inter[[1]], loo_bottomup_exploit[[1]], loo_bottomup_inter[[1]]) # loo
  #' loo_compare(loo_topdown_exploit[[2]], loo_topdown_inter[[2]], loo_bottomup_exploit[[2]], loo_bottomup_inter[[2]]) # waic
  
  
  