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
  
  
  # #####  Bottom-up, Top-down, Interference & Exploitative Hybrid  #####
  # data_JAGS_bundle_hybrid <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
  #                                         dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
  #                                         covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
  #                                         covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
  #                                         nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, 
  #                                         nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  # num.chains <- 3
  # initsList_hybrid <- vector('list', num.chains) 
  # for(i in 1:num.chains) {
  #   initsList_hybrid[[i]] <- generate_inits(nwolf = 1, nlion = 1, nbear = 1, ncoy = 1, nelk = 4, nmoose = 2, 
  #                                             nwtd = 3, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  # }
  # source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_posthoc_hybrid.R")
  # start.time = Sys.time()
  # SEM_hybrid <- jagsUI::jags(data_JAGS_bundle_hybrid, inits = initsList_hybrid, params, 
  #                              "./Outputs/SEM/JAGS_out/JAGS_SEM_posthoc_hybrid.txt",
  #                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
  #                              n.burnin = nb, parallel = TRUE)
  # end.time <- Sys.time(); (run.time <- end.time - start.time)
  # print(SEM_hybrid$summary)
  # which(SEM_hybrid$summary[,"Rhat"] < 0.9)
  # which(SEM_hybrid$summary[,"Rhat"] > 1.1)
  # mcmcplot(SEM_hybrid$samples)
  # save(SEM_hybrid, file = paste0("./Outputs/SEM/JAGS_out/SEM_posthoc_hybrid_", Sys.Date(), ".RData"))
  

  
  
  
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
    
    site_results <- future_lapply(1:nSites, function(s) {
      #'  Copy the full data set
      fold_data <- dat
      
      #'  Then create this fold's testing data by setting the withheld site's spp.hat 
      #'  to NA for every species and year
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
      
      #'  Compute posterior predictive log density of TRUE observed value for each 
      #'  species/year at the held-out site: log(mean_over_draws( p(y_true | draw's latent state))) 
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
          # site_detail[[paste0(sp, "_yr", y)]] <- lpd
        }
      }
      
      #'  Save score for every species/year for site s
      saveRDS(list(site = s, elpd = site_loglik), file.path (out_dir, sprintf("site_%03d.rds", s))) #, detail = site_detail
      list(site = s, elpd = site_loglik)
      
    }, future.seed = TRUE)
    
    #'  Sum every site's score into an overall score for the entire model: expected log predictive density (elpd)
    #'  Higher elpd means indicate the model, on average, made better honest predictions
    #'  about sites it did not see. se_elpd estimates how much that total could plausibly
    #'  vary, treating each score as on independent data point
    elpd_per_site <- sapply(site_results, function(x) x$elpd)
    
    list(elpd_per_site = elpd_per_site,
         total_elpd = sum(elpd_per_site),
         se_elpd = sd(elpd_per_site) * sqrt(nSites), # treats sites as the exchangeable unit
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
                                      out_dir = "./Outputs/SEM/LOSO_CV/bottomup_exploit")
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(loso_bottomup_exploit, file = "./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_exploit.RData")
  
  start.time = Sys.time()
  loso_bottomup_inter <- run_loso_cv(mod = "./Outputs/SEM/JAGS_out/JAGS_SEM_bottomup_inter_final.txt", dat = data_JAGS_bundle_bottomup_inter_final,
                                      spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd"), nSites = 23,  
                                      nYear = 4, n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt,
                                      out_dir = "./Outputs/SEM/LOSO_CV/bottomup_inter")
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(loso_bottomup_inter, file = "./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_inter.RData")
  
  #' #'  Read in LOSO_CV results
  #' load("./Outputs/SEM/LOSO_CV/LOSO_CV_topdown_exploit.RData")
  #' load("./Outputs/SEM/LOSO_CV/LOSO_CV_topdown_inter.RData")
  #' load("./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_exploit.RData")
  #' load("./Outputs/SEM/LOSO_CV/LOSO_CV_bottomup_inter.RData")
  
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
  
  #'  LOSO CV results table
  loso_tbl <- compare_loso %>%
    transmute(Model = model, ELPD = elpd, SE = se, elpd_diff = elpd_diff) %>%
    mutate(Model = ifelse(Model == "topdown_exploitative", "Top-down, exploitative", Model),
           Model = ifelse(Model == "topdown_interference", "Top-down, interference", Model),
           Model = ifelse(Model == "bottomup_exploitative", "Bottom-up, exploitative", Model),
           Model = ifelse(Model == "bottomup_interference", "Bottom-up, interference", Model)) %>%
    arrange(-elpd_diff)
  # names(loso_tbl)[names(loso_tbl) == "elpd_diff"] <- paste0("\u0394", "ELPD")  #' "ΔELPD"
  print(loso_tbl)
  write_csv(loso_tbl, "./Outputs/SEM/Tables_for_publication/LOSO_CV_results.csv")
  
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
  #'  How does each model compare to top-model?
  # (compare_topdown <- compare_paired_mods(loso_bottomup_inter$elpd_per_site, loso_topdown_exploit$elpd_per_site))
  # (compare_bottomup <- compare_paired_mods(loso_bottomup_inter$elpd_per_site, loso_bottomup_exploit$elpd_per_site))
  # (compare_topdown_inter <- compare_paired_mods(loso_bottomup_inter$elpd_per_site, loso_topdown_inter$elpd_per_site))
  (compare_bottomup <- compare_paired_mods(loso_topdown_exploit$elpd_per_site, loso_bottomup_exploit$elpd_per_site))
  (compare_bottomup_inter <- compare_paired_mods(loso_topdown_exploit$elpd_per_site, loso_bottomup_inter$elpd_per_site))
  (compare_topdown_inter <- compare_paired_mods(loso_topdown_exploit$elpd_per_site, loso_topdown_inter$elpd_per_site))
  
  #'  Create table of pairwise results
  # loso_pair <- c("Top-down, exploitative", "Bottom-up, exploitative", "Top-down, interference")
  # loso_compare_tbl <- data.frame(Model = loso_pair, 
  #                                elpd_diff = c(compare_topdown$elpd_diff, compare_bottomup$elpd_diff, compare_topdown_inter$elpd_diff),
  #                                se_diff = c(compare_topdown$se_diff, compare_bottomup$se_diff, compare_topdown_inter$se_diff),
  #                                Ratio = c(compare_topdown$ratio, compare_bottomup$ratio, compare_topdown_inter$ratio)) %>%
  #   arrange(Ratio)
  loso_pair <- c("Bottom-up, exploitative", "Bottom-up, interference", "Top-down, interference")
  loso_compare_tbl <- data.frame(Model = loso_pair, 
                                 elpd_diff = c(compare_bottomup$elpd_diff, compare_bottomup_inter$elpd_diff, compare_topdown_inter$elpd_diff),
                                 se_diff = c(compare_bottomup$se_diff, compare_bottomup_inter$se_diff, compare_topdown_inter$se_diff),
                                 Ratio = c(compare_bottomup$ratio, compare_bottomup_inter$ratio, compare_topdown_inter$ratio)) %>%
    arrange(Ratio)
  # names(loso_compare_tbl)[names(loso_compare_tbl) == "elpd_diff"] <- paste0("\u0394", "ELPD")  #' "ΔELPD"
  # names(loso_compare_tbl)[names(loso_compare_tbl) == "se_diff"] <- paste0("\u0394", "SE")  #' "ΔELPD"
  print(loso_compare_tbl)
  # write_csv(loso_compare_tbl, "./Outputs/SEM/Tables_for_publication/LOSO_CV_compared2TopDown_results.csv")
  
  loso_full_tbl <- full_join(loso_tbl, loso_compare_tbl, by = c("Model")) %>%
    dplyr::select(-elpd_diff.x) %>%
    transmute(Model = Model,
           ELPD = round(ELPD, 2),
           SE = round(SE, 2),
           elpd_diff = round(elpd_diff.y, 2),
           se_diff = round(se_diff, 2),
           Ratio = round(Ratio, 2))
  names(loso_full_tbl)[names(loso_full_tbl) == "elpd_diff"] <- paste0("\u0394", "ELPD")  #' "ΔELPD"
  names(loso_full_tbl)[names(loso_full_tbl) == "se_diff"] <- paste0("SE(", "\u0394", "ELPD)")  #' "ΔELPD"
  print(loso_full_tbl)
  write_csv(loso_full_tbl, "./Outputs/SEM/Tables_for_publication/LOSO_CV_full_results.csv")
  
  #'  bottomup_inter and topdown_exploit are essentially indistinguishable 
  #'  bottomup_exploit is not meaningfully different from bottomup_inter
  #'  topdown_inter is distinctly different and less supported than bottomup_inter
  #'  topdown_exploit and bottomup_exploit are statistically indistinguishable alternatives 
  #'  to bottomup_inter given the data
  
  #'  -----------------------------
  #####  Stacking model averaging  #####
  #'  -----------------------------
  #'  Generate log predictive densities (LPD) matrix (dims = nSite x model)
  # lpd_matrix <- cbind(
  #   bottomup_interference = loso_bottomup_inter$elpd_per_site,
  #   topdown_exploitative   = loso_topdown_exploit$elpd_per_site,
  #   bottomup_exploitative  = loso_bottomup_exploit$elpd_per_site
  # )
  lpd_matrix <- cbind(
    topdown_exploitative   = loso_topdown_exploit$elpd_per_site,
    bottomup_exploitative  = loso_bottomup_exploit$elpd_per_site,
    bottomup_interference = loso_bottomup_inter$elpd_per_site
  )
 
  #'  Compute stacking weights using the loo::stacking_weights function
  stack_wgts <- stacking_weights(lpd_point = lpd_matrix) 
  print(stack_wgts)
  #'  Compute pseudo-Bayesian Model Averaging (BMA) for comparison 
  #'  Yao et al. (2018) suggest this is second best approach
  pbma_wgts <- pseudobma_weights(lpd_point = lpd_matrix, BB = TRUE)
  print(pbma_wgts)
  
  #'  Calculate model-averaged (stacked) log predictive density for each observation
  #'  using the log-sum-exp method
  wgtd_logsumexp <- function(lpd_row, w) {
    #'  Identify maximum lpd
    m <- max(lpd_row)
    #'  Weight the exponentiated lpd (that was scaled by the max lpd) by the model weights,
    #'  then sum, log and add the max lpd back on 
    wgtd_site_loglik <- m + log(sum(w * exp(lpd_row - m)))
    return(wgtd_site_loglik)
  }
  stacked_lpd <- apply(lpd_matrix, 1, wgtd_logsumexp, w = stack_wgts) 
  stacked_total_elpd <- sum(stacked_lpd)
  
  #'  Print ELPD for stacked vs original best model
  print(stacked_total_elpd)
  print(max(colSums(lpd_matrix)))
  
  #'  Compare stacked model's performance 
  stacked_elpd_per_site <- stacked_lpd
  se_stacked <- sd(stacked_elpd_per_site) * sqrt(length(stacked_elpd_per_site))
  
  updated_loso_tbl <- data.frame(Model = c("Model-averaged (stacked)",
                                           "Top-down, exploitative", 
                                           "Bottom-up, exploitative", 
                                           "Bottom-up, interference"),
                                           # "Bottom-up, interference",
                                           # "Top-down, exploitative", 
                                           # "Bottom-up, exploitative"),
                                 ELPD = c(sum(stacked_elpd_per_site),
                                          loso_topdown_exploit$total_elpd,
                                          loso_bottomup_exploit$total_elpd,
                                          loso_bottomup_inter$total_elpd),
                                         # loso_bottomup_inter$total_elpd,
                                         # loso_topdown_exploit$total_elpd,
                                         # loso_bottomup_exploit$total_elpd),
                                 SE = c(se_stacked,
                                        loso_topdown_exploit$se_elpd,
                                        loso_bottomup_exploit$se_elpd,
                                        loso_bottomup_inter$se_elpd)) %>%
                                        # loso_bottomup_inter$se_elpd,
                                        # loso_topdown_exploit$se_elpd,
                                        # loso_bottomup_exploit$se_elpd)) %>%
    mutate(ELPD = round(ELPD, 2),
           SE = round(SE, 2))
  updated_loso_tbl <- updated_loso_tbl[order(-updated_loso_tbl$ELPD),]
  updated_loso_tbl$deltaELPD <- updated_loso_tbl$ELPD - updated_loso_tbl$ELPD[1]
  names(updated_loso_tbl)[names(updated_loso_tbl) == "deltaELPD"] <- paste0("\u0394", "ELPD")  #' "ΔELPD"
  rownames(updated_loso_tbl) <- NULL
  mod_wgts <- c(NA, round(stack_wgts[1], 2), round(stack_wgts[2], 2), round(stack_wgts[3], 2))
  updated_loso_tbl$Weights <- mod_wgts
  print(updated_loso_tbl)
  write.csv(updated_loso_tbl, "./Outputs/SEM/Tables_for_publication/LOSO_CV_results_modelAvg.csv")
  
  #'  Pairwise comparison between stacked vs best-supported individual model 
  #'  to determine if averaging made a large enough improvement to matter.
  #'  Ratio >= ~2 suggests stacked model meaningfully outperforms original best model.
  best_individual <- colnames(lpd_matrix)[which.max(colSums(lpd_matrix))]
  best_individual_per_site <- lpd_matrix[, best_individual]
  
  compare_paired_mods(stacked_elpd_per_site, best_individual_per_site)
  
  
  #'  ----------------------------------------------
  #####  In-sample expected log predictive density  #####
  #'  ----------------------------------------------
  compute_insample_elpd <- function(mod, dat, spp_names, nSites, nYear) {
    #'  Store unmodified spp.hat and spp.sigma_hat arrays 
    true_hat <- setNames(lapply(spp_names, function(sp) dat[[paste0(sp, ".hat")]]), spp_names)
    true_sigma_hat <- setNames(lapply(spp_names, function(sp) dat[[paste0(sp, ".sigma_hat")]]), spp_names)
    
    #'  Grab site number
    insample_elpd_per_site <- numeric(nSites)
    
    for (s in 1:nSites) {
      #'  Set loglikelihood for each site to 0 to start
      site_loglik <- 0
      for (sp in spp_names) {
        #'  Grab the latent posterior for each species
        latent_node <- mod$sims.list[[paste0(sp, ".latent")]]
        for (y in 1:nYear) {
          #'  Grab the unmodified spp.hat and spp.sigma for each species and year
          y_true <- true_hat[[sp]][s, y]
          if (is.na(y_true)) next
          sigma_true <- true_sigma_hat[[sp]][s, y]
          #'  Grab every posterior draw of the original model
          draws <- latent_node[, s, y]
          #'  Compute the log-likelihood of the observed spp.hat per draw using 
          #'  the known latent posterior and measure of uncertainty from RN model 
          #'  as the spread. This returns the density of how each draw in the posterior 
          #'  predicted what actually happened
          logdens_draws <- dnorm(y_true, mean = draws, sd = sigma_true, log = TRUE)
          #'  Generate an overall predictive score
          m <- max(logdens_draws)
          site_loglik <- site_loglik + m + log(mean(exp(logdens_draws - m)))
        }
      }
      insample_elpd_per_site[s] <- site_loglik
    }
    insample_elpd_per_site
  }
  
  #'  Compute in-sample ELPD for the top model as indicated by Fisher's C AICc and LOSO-CV
  insample_bottomup <- compute_insample_elpd(SEM_bottomup_final, dat = data_JAGS_bundle_bottomup_final, 
                                    spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd"), 
                                    nSites = 23, nYear = 4)  # top model based on Fisher's AICc
  insample_bottomup_inter <- compute_insample_elpd(SEM_bottomup_inter_final, dat = data_JAGS_bundle_bottomup_inter_final, 
                                    spp_names = c("lion", "wolf", "bear", "coy", "elk", "moose", "wtd"), 
                                    nSites = 23, nYear = 4) # top model based on LOSO-CV
  
  gap_per_site_bottomup <- insample_bottomup - loso_bottomup_exploit$elpd_per_site
  cat(sprintf("In-sample total ELPD: %.2f\n", sum(gap_per_site_bottomup)))
  cat(sprintf("Out-of-sample (LOSO) total ELPD: %.2f\n", loso_bottomup_exploit$total_elpd))
  cat(sprintf("Overfitting gap: %.2f\n", sum(gap_per_site_bottomup) - loso_bottomup_exploit$total_elpd))
  
  
  gap_per_site_bottomup_inter <- insample_bottomup_inter - loso_bottomup_inter$elpd_per_site
  cat(sprintf("In-sample total ELPD: %.2f\n", sum(insample_bottomup_inter)))
  cat(sprintf("Out-of-sample (LOSO) total ELPD: %.2f\n", loso_bottomup_inter$total_elpd))
  cat(sprintf("Overfitting gap: %.2f\n", sum(insample_bottomup_inter) - loso_bottomup_inter$total_elpd))
  
  
  #'  -----------------------
  ####  Model Result Tables  ####
  #'  -----------------------
  #####  Coefficients and 95% CRI  #####
  #'  -----------------------------
  #'  Results table of direct effects for each SEM
  #'  Load model outputs
  load("./Outputs/SEM/JAGS_out/SEM_topdown_exploitation_final_2026-10-05.RData") #2026-09-13
  load("./Outputs/SEM/JAGS_out/SEM_topdown_inter_final_2026-10-05.RData")
  load("./Outputs/SEM/JAGS_out/SEM_bottomup_exploitation_final_2026-10-05.RData")
  load("./Outputs/SEM/JAGS_out/SEM_bottomup_inter_final_2026-10-05.RData")
  
  grab_coefs <- function(mod_out) {
    coef_mean <- mod_out$mean 
    coef_ll <- mod_out$q2.5
    coef_ul <- mod_out$q97.5
    coef_overlap0 <- mod_out$overlap0
    coef_list <- list(coef_mean, coef_ll, coef_ul, coef_overlap0)
    names(coef_list) <- c("coef_mean", "coef_ll", "coef_ul", "coef_overlap0")
    
    #'  Retain just the beta coefficients (no indirect effects, spp.latent, sigma, etc.)
    for(i in 1:length(coef_list)) {
      coef_list[[i]] <- coef_list[[i]] %>% 
        subset(., startsWith(names(.), "beta")) 
    }
    return(coef_list)
  }
  topdown_exploit_list <- grab_coefs(SEM_topdown_final)
  topdown_inter_list <- grab_coefs(SEM_topdown_inter_final)
  bottomup_exploit_list <- grab_coefs(SEM_bottomup_final)
  bottomup_inter_list <- grab_coefs(SEM_bottomup_inter_final)
  
  #'  Function to extract and reformat all coefficients for a given parameter into a data.frame
  coef_tbl <- function(coef_list, focal_param) {
    #'  Grab mean, 95% CRI and 0 overlap indicator per focal parameter
    param_mu <- as.data.frame(coef_list$coef_mean[focal_param])
    param_ll <- as.data.frame(coef_list$coef_ll[focal_param])
    param_ul <- as.data.frame(coef_list$coef_ul[focal_param])
    param_overlap0 <- as.data.frame(coef_list$coef_overlap0[focal_param])
    
    #'  Bring all together into a single data frame
    param.df <- data.frame(param_mu, param_ll, param_ul, param_overlap0)
    
    #'  Add a column indicating the parameter name and its index
    param.df$param <- sprintf(paste0(focal_param, "[%d]"), as.numeric(rownames(param.df)))
    
    #'  Name and reorganize columns
    names(param.df) <- c("mean", "ll", "ul", "overlap0", "param")
    param.df <- param.df %>% relocate(param, .before = mean) %>%
      mutate(mean = round(mean, 2),
             ll = round(ll, 2),
             ul = round(ul, 2),
             overlap0 = ifelse(overlap0 == TRUE, "T", "F"))
    
    return(param.df)
  }
  #'  List parameter names to reformat into a single table
  focal_param_topdown <- list("beta.int.tmin1", "beta.int", "beta.lion", "beta.wolf", "beta.bear", 
                              "beta.coy", "beta.elk", "beta.moose", "beta.wtd", "beta.harvest")
  focal_param_bottomup <- list("beta.int.tmin1", "beta.int", "beta.wolf", "beta.bear", #"beta.lion", 
                               "beta.coy", "beta.elk", "beta.moose", "beta.wtd", "beta.forest", "beta.wsi")
  #'  Call function and convert to a single table per model
  topdown_exploit_coefs <- lapply(focal_param_topdown, coef_tbl, coef_list = topdown_exploit_list) %>%
    bind_rows()
  topdown_inter_coefs <- lapply(focal_param_topdown, coef_tbl, coef_list = topdown_inter_list) %>%
    bind_rows()
  bottomup_exploit_coefs <- lapply(focal_param_bottomup, coef_tbl, coef_list = bottomup_exploit_list) %>%
    bind_rows()
  bottomup_inter_coefs <- lapply(focal_param_bottomup, coef_tbl, coef_list = bottomup_inter_list) %>%
    bind_rows()
  
  #'  Order coefficients by species-specific regression
  topdown_exploit_coef_order <- c("beta.int.tmin1[1]", "beta.int[1]", "beta.harvest[1]", "beta.elk[2]", "beta.wtd[2]", "beta.wtd[3]", 
                                  "beta.int.tmin1[2]", "beta.int[2]", "beta.wolf[1]", "beta.harvest[2]", "beta.moose[2]", "beta.elk[3]", "beta.wtd[4]", "beta.wtd[5]", 
                                  "beta.int.tmin1[3]", "beta.int[3]", "beta.bear[1]", "beta.harvest[3]", "beta.wtd[6]",
                                  "beta.int.tmin1[4]", "beta.int[4]", "beta.coy[1]", "beta.wtd[7]", 
                                  "beta.int.tmin1[5]", "beta.int[5]", "beta.elk[1]", "beta.wolf[2]", "beta.lion[1]", "beta.harvest[4]", "beta.bear[2]", 
                                  "beta.int.tmin1[6]", "beta.int[6]", "beta.moose[1]", "beta.wolf[3]", 
                                  "beta.int.tmin1[7]", "beta.int[7]", "beta.wtd[1]", "beta.lion[2]", "beta.harvest[5]", "beta.wolf[4]", "beta.coy[2]")
  topdown_exploit_data_order <- c("Intercept t-1", "Intercept t", "Lion harvest t-1", "Elk t-1", "Deer t-1", "Deer t", 
                                  "Intercept t-1", "Intercept t", "Wolf t-1", "Wolf harvest t-1", "Moose t-1", "Elk t-1", "Deer t-1", "Deer t", 
                                  "Intercept t-1", "Intercept t", "Bear t-1", "Bear harvest t-1", "Deer t",
                                  "Intercept t-1", "Intercept t", "Coy t-1", "Deer t-1", 
                                  "Intercept t-1", "Intercept t", "Elk t-1", "Wolf t-1", "Lion t-1", "Elk harvest t-1", "Bear t-1", 
                                  "Intercept t-1", "Intercept t", "Moose t-1", "Wolf t-1", 
                                  "Intercept t-1", "Intercept t", "Deer t-1", "Lion t-1", "Deer harvest t-1", "Wolf t-1", "Coy t-1")
  sub_reg_spp_topdown_ex <- c(rep("Mountain lion", 6), rep("Wolf", 8), rep("Black bear", 5), rep("Coyote", 4), rep("Elk", 7), rep("Moose", 4), rep("White-tailed deer", 7))
  topdown_exploit_coefs_tbl <- topdown_exploit_coefs %>% arrange(factor(param, levels = topdown_exploit_coef_order)) %>%
    mutate(Model = "Top-down, Exploitative", 
           Species_submodel = sub_reg_spp_topdown_ex, 
           #'  Name parameters based on something more meaningful than their indexing
           Parameter = topdown_exploit_data_order) %>% 
    relocate(Model, .before = param) %>% relocate(Species_submodel, .before = param) %>%
    relocate(Parameter, .before = param)
  
  topdown_inter_coef_order <- c("beta.int.tmin1[1]", "beta.int[1]", "beta.harvest[1]", "beta.wolf[2]", 
                                "beta.int.tmin1[2]", "beta.int[2]", "beta.wolf[1]", "beta.harvest[2]", "beta.bear[2]", 
                                "beta.int.tmin1[3]", "beta.int[3]", "beta.bear[1]", "beta.harvest[3]", "beta.wolf[3]", 
                                "beta.int.tmin1[4]", "beta.int[4]", "beta.coy[1]", "beta.wolf[4]", "beta.lion[1]", 
                                "beta.int.tmin1[5]", "beta.int[5]", "beta.elk[1]", "beta.wolf[5]", "beta.lion[2]", 
                                "beta.int.tmin1[6]", "beta.int[6]", "beta.moose[1]", "beta.wolf[6]",
                                "beta.int.tmin1[7]", "beta.int[7]", "beta.wtd[1]", "beta.lion[3]", "beta.coy[2]", "beta.bear[3]")
  topdown_inter_data_order <- c("Intercept t-1", "Intercept t", "Lion harvest t-1", "Wolf t-1", 
                                "Intercept t-1", "Intercept t", "Wolf t-1", "Wolf harvest t-1", "Bear t", 
                                "Intercept t-1", "Intercept t", "Bear t-1", "Bear harvest t-1", "Wolf t-1",  
                                "Intercept t-1", "Intercept t", "Coy t-1", "Wolf t-1", "Lion t-1",  
                                "Intercept t-1", "Intercept t", "Elk t-1", "Wolf t-1", "Lion t-1", 
                                "Intercept t-1", "Intercept t", "Moose t-1", "Wolf t-1", 
                                "Intercept t-1", "Intercept t", "Deer t-1", "Lion t-1", "Coy t-1", "Bear t-1")
  sub_reg_spp_topdown_int <- c(rep("Mountain lion", 4), rep("Wolf", 5), rep("Black bear", 5), rep("Coyote", 5), rep("Elk", 5), rep("Moose", 4), rep("White-tailed deer", 6))
  topdown_inter_coefs_tbl <- topdown_inter_coefs %>% arrange(factor(param, levels = topdown_inter_coef_order)) %>%
    mutate(Model = "Top-down, Interference", 
           Species_submodel = sub_reg_spp_topdown_int,
           Parameter = topdown_inter_data_order) %>% 
    relocate(Model, .before = param) %>% relocate(Species_submodel, .before = param) %>%
    relocate(Parameter, .before = param)
  
  bottomup_exploit_coefs_order <- c("beta.int.tmin1[1]", "beta.int[1]", "beta.elk[2]", "beta.elk[3]", "beta.wtd[2]", "beta.wtd[3]", 
                                    "beta.int.tmin1[2]", "beta.int[2]", "beta.wolf[1]", "beta.elk[4]", "beta.elk[5]", "beta.moose[2]", "beta.moose[3]", "beta.wtd[4]",
                                    "beta.int.tmin1[3]", "beta.int[3]", "beta.bear[1]", "beta.elk[6]", "beta.forest[4]", "beta.wtd[5]", 
                                    "beta.int.tmin1[4]", "beta.int[4]", "beta.coy[1]", "beta.wtd[6]", "beta.wtd[7]", 
                                    "beta.int.tmin1[5]", "beta.int[5]", "beta.elk[1]", "beta.forest[1]", "beta.wsi[1]", "beta.bear[2]", 
                                    "beta.int.tmin1[6]", "beta.int[6]", "beta.moose[1]", "beta.forest[2]", "beta.wsi[2]",
                                    "beta.int.tmin1[7]", "beta.int[7]", "beta.wtd[1]", "beta.forest[3]", "beta.wsi[3]", "beta.coy[2]", "beta.bear[3]")
  bottomup_exploit_data_order <- c("Intercept t-1", "Intercept t", "Elk t-1", "Elk t", "Deer t-1", "Deer t", 
                                   "Intercept t-1", "Intercept t", "Wolf t-1", "Elk t-1", "Elk t", "Moose t-1", "Moose t", "Deer t", 
                                   "Intercept t-1", "Intercept t", "Bear t-1", "Elk t-1", "Forest t-1", "Deer t", 
                                   "Intercept t-1", "Intercept t", "Coy t-1", "Deer t-1", "Deer t", 
                                   "Intercept t-1", "Intercept t", "Elk t-1", "Forest t-1", "WSI t-1", "Bear t-1", 
                                   "Intercept t-1", "Intercept t", "Moose t-1", "Forest t-1", "WSI t-1",
                                   "Intercept t-1", "Intercept t", "Deer t-1", "Forest t-1", "WSI t-1", "Coy t-1", "Bear t-1")
  sub_reg_spp_bottomup_ex <- c(rep("Mountain lion", 6), rep("Wolf", 8), rep("Black bear", 6), rep("Coyote", 5), rep("Elk", 6), rep("Moose", 5), rep("White-tailed deer", 7))
  bottomup_exploit_coefs_tbl <- bottomup_exploit_coefs %>% arrange(factor(param, levels = bottomup_exploit_coefs_order)) %>%
    mutate(Model = "Bottom-up, Exploitative", 
           Species_submodel = sub_reg_spp_bottomup_ex,
           Parameter = bottomup_exploit_data_order) %>% 
    relocate(Model, .before = param) %>% relocate(Species_submodel, .before = param) %>%
    relocate(Parameter, .before = param)
  
  bottomup_inter_coefs_order <- c("beta.int.tmin1[1]", "beta.int[1]", "beta.elk[2]", "beta.elk[3]", "beta.wtd[2]", "beta.wtd[3]", "beta.bear[2]", 
                                  "beta.int.tmin1[2]", "beta.int[2]", "beta.wolf[1]", "beta.elk[4]", "beta.elk[5]", "beta.moose[2]", "beta.moose[3]", "beta.bear[3]", "beta.bear[4]", 
                                  "beta.int.tmin1[3]", "beta.int[3]", "beta.bear[1]", "beta.elk[6]", "beta.forest[4]", "beta.wolf[2]", "beta.wtd[4]", 
                                  "beta.int.tmin1[4]", "beta.int[4]", "beta.coy[1]", "beta.wtd[5]", "beta.wtd[6]", "beta.wolf[3]", "beta.bear[5]", 
                                  "beta.int.tmin1[5]", "beta.int[5]", "beta.elk[1]", "beta.forest[1]", "beta.wsi[1]",
                                  "beta.int.tmin1[6]", "beta.int[6]", "beta.moose[1]", "beta.forest[2]", "beta.wsi[2]",
                                  "beta.int.tmin1[7]", "beta.int[7]", "beta.wtd[1]", "beta.forest[3]", "beta.wsi[3]")
  bottomup_inter_data_order <- c("Intercept t-1", "Intercept t", "Elk t-1", "Elk t", "Deer t-1", "Deer t", "Bear t-1", 
                                 "Intercept t-1", "Intercept t", "Wolf t-1", "Elk t-1", "Elk t", "Moose t-1", "Moose t", "Bear t-1", "Bear t", 
                                 "Intercept t-1", "Intercept t", "Bear t-1", "Elk t-1", "Forest t-1", "Wolf t-1", "Deer t", 
                                 "Intercept t-1", "Intercept t", "Coy t-1", "Deer t-1", "Deer t", "Wolf t-1", "Bear t-1", 
                                 "Intercept t-1", "Intercept t", "Elk t-1", "Forest t-1", "WSI t-1", 
                                 "Intercept t-1", "Intercept t", "Moose t-1", "Forest t-1", "WSI t-1", 
                                 "Intercept t-1", "Intercept t", "Deer t-1", "Forest t-1", "WSI t-1")
  sub_reg_spp_bottomup_int <- c(rep("Mountain lion", 7), rep("Wolf", 9), rep("Black bear", 7), rep("Coyote", 7), rep("Elk", 5), rep("Moose", 5), rep("White-tailed deer", 5))
  bottomup_inter_coefs_tbl <- bottomup_inter_coefs %>% arrange(factor(param, levels = bottomup_inter_coefs_order)) %>%
    mutate(Model = "Bottom-up, Interference", 
           Species_submodel = sub_reg_spp_bottomup_int,
           Parameter = bottomup_inter_data_order) %>% 
    relocate(Model, .before = param) %>% relocate(Species_submodel, .before = param) %>%
    relocate(Parameter, .before = param)
  
  SEM_final_coefs_tbl <- bind_rows(bottomup_exploit_coefs_tbl, bottomup_inter_coefs_tbl, topdown_inter_coefs_tbl, topdown_exploit_coefs_tbl)
  write_csv(SEM_final_coefs_tbl, "./Outputs/SEM/Tables_for_publication/SEM_coefficients_table_2026-10-05.csv")
  
  #'  ----------------------------
  #####  Coefficient comparisons  #####
  #'  ----------------------------
  #'  Results table to quickly compare similar and dissimilar relationships across SEMs
  SEM_final_coefs_tbl_skinny <- SEM_final_coefs_tbl %>%
    filter(Parameter != "Intercept t-1") %>%
    filter(Parameter != "Intercept t") %>%
    mutate(mean_overlap0 = paste(mean, overlap0)) %>%
    dplyr::select(c(Model, Species_submodel, Parameter, mean_overlap0)) %>%
    pivot_wider(names_from = Species_submodel, values_from = mean_overlap0)
  
  SEM_final_coefs_tbl_skinny_symbols <- SEM_final_coefs_tbl %>%
    filter(Parameter != "Intercept t-1") %>%
    filter(Parameter != "Intercept t") %>%
    mutate(mean_direction = ifelse(mean >= 0, " + ", mean),
           mean_direction = ifelse(mean < 0, " - ", mean_direction), 
           mean_direction = ifelse(is.na(mean_direction), "NA", mean_direction)) %>% 
    dplyr::select(c(Model, Species_submodel, Parameter, mean_direction)) %>%
    pivot_wider(names_from = Species_submodel, values_from = mean_direction)
  
  write_csv(SEM_final_coefs_tbl_skinny, "./Outputs/SEM/Tables_for_publication/SEM_direct_effects_table.csv")
  write_csv(SEM_final_coefs_tbl_skinny_symbols, "./Outputs/SEM/Tables_for_publication/SEM_direct_effect_directions_table.csv")
  
  