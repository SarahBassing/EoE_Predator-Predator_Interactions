  #'  ---------------------------------------
  #'  Bayesian multispecies occupancy model
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  ---------------------------------------
  #'  Script sources data formatting scripts and JAGS code to run single-season
  #'  multispecies occupancy models. Code adapted from Kery & Royle AHM2 book 
  #'  (Ch. 8.2.3). Code runs 2-spp models using detection/non-detection data for 
  #'  5 predator species: black bears, bobcats, coyotes, mountain lions, and wolves
  #'  detected during summers 2020 & 2021 (July - mid Sept) via camera traps in 
  #'  northern Idaho. Co-occurrence models include 11 7-day sampling occasions.
  #'  Competing models test whether co-occurrence is non-independent and whether
  #'  predator occurrence and co-occurrence are influenced by habitat, prey, and/or
  #'  anthropogenic factors.
  #'  
  #'  Detection histories are generated with Detection_histories_for_occmod.R. 
  #'  Covariate data are formatted with Format_data_2spp_occmod_for_JAGS.R
  #'  ---------------------------------------
  
  #'  Clean work space and load libraries
  rm(list = ls())
  
  # library(R2jags)
  library(jagsUI)
  library(abind)
  library(mcmcplots)
  # library(loo)
  library(tidyverse)
  
  #'  Source detection histories and covariate data
  source("./Scripts/MultiSpp_OccMod/Detection_histories_for_occmod.R")
  source("./Scripts/MultiSpp_OccMod/Format_data_2spp_occmod_for_JAGS.R")
  
  
  #'  ----------------------------
  ####  Pair detection histories  ####
  #'  ----------------------------
  #'  Combine annual detection histories
  all_detections <- function(dh1, dh2) {
    dh1 <- dh1[[1]]
    dh2 <- dh2[[1]]
    
    #'  Rename columns in each detection history
    newcols <- c("occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11")
    colnames(dh1) <- newcols; colnames(dh2) <- newcols
    
    #'  Bind year 1 and 2 data together 
    #'  Remember that annual detection and covariate data are each arranged by  
    #'  location then stacked!
    dh <- rbind(dh1, dh2) 
    
    return(dh)
  }
  DH_bear <- all_detections(DH_eoe20s_predators[[1]], DH_eoe21s_predators[[1]])
  DH_bob <- all_detections(DH_eoe20s_predators[[2]], DH_eoe21s_predators[[2]])
  DH_coy <- all_detections(DH_eoe20s_predators[[3]], DH_eoe21s_predators[[3]])
  DH_lion <- all_detections(DH_eoe20s_predators[[4]], DH_eoe21s_predators[[4]])
  DH_wolf <- all_detections(DH_eoe20s_predators[[5]], DH_eoe21s_predators[[5]])
  
  #'  Combine species detection histories into an array (site x survey x species)
  detection_array <- function(spp1, spp2, name1, name2) {
    #'  List detection histories
    spp12_DH <- list(spp1, spp2)
    #'  Name lists based on specific species pairing
    names(spp12_DH) <- c(name1, name2)
    #'  Format list into a 3D array
    spp12_array <- abind(spp12_DH, along = 3)
    
    return(spp12_array)
  }
  wolf_bear_DH <- detection_array(spp1 = DH_wolf, spp2 = DH_bear, name1 = "wolf", name2 = "bear")
  wolf_coy_DH <- detection_array(spp1 = DH_wolf, spp2 = DH_coy, name1 = "wolf", name2 = "coyote")
  wolf_lion_DH <- detection_array(spp1 = DH_wolf, spp2 = DH_lion, name1 = "wolf", name2 = "lion")
  lion_bear_DH <- detection_array(spp1 = DH_lion, spp2 = DH_bear, name1 = "lion", name2 = "bear")
  lion_bob_DH <- detection_array(spp1 = DH_lion, spp2 = DH_bob, name1 = "lion", name2 = "bobcat")
  coy_bob_DH <- detection_array(spp1 = DH_coy, spp2 = DH_bob, name1 = "coyote", name2 = "bobcat")
  
  #'  List 2-species detection arrays together for faster formatting below
  DH_array_list <- list(wolf_bear_DH, wolf_coy_DH, wolf_lion_DH, lion_bear_DH, lion_bob_DH, coy_bob_DH)
  
  
  #'  -------------------------
  ####  Format data for JAGS  ####
  #'  -------------------------
  #'  Define survey dimensions
  nsites <- dim(wolf_bear_DH)[1]
  nsurveys <- dim(wolf_bear_DH)[2]
  nspecies <- dim(wolf_bear_DH)[3]
  #'  Number of possible community states (species interactions): 00, 10, 01, 11
  ncat <- 2^nspecies
  
  #####  Format detection histories  ####
  #'  --------------------------------
  #'  Function to convert 3D detection array into 2D multi-species detection history
  observation_state <- function(array_list) {
    #'  Merge species-specific observations into single 2-species observation state
    #'  (e.g., 00, 01, NANA) for each site and survey occasion
    ycat <- apply(array_list, c(1,2), paste, collapse = "")
    #'  Reclassify each observation state
    ycat[ycat == "00"] <- 1      # Neither species detected
    ycat[ycat == "10"] <- 2      # Only predator 1 (spp1) detected
    ycat[ycat == "01"] <- 3      # Only predator 2 (spp2) detected
    ycat[ycat == "11"] <- 4      # Both predator species detected
    ycat[ycat =="NANA"] <- NA    # Not sampled, no data
    #'  Make each column numeric so JAGS knows to handle it as a response variable
    ycat <- apply(ycat, 2, as.numeric)
    return(ycat)
  }
  multi_spp_DH_list <- lapply(DH_array_list, observation_state)
  names(multi_spp_DH_list) <- c("wolf_bear_DH", "wolf_coy_DH", "wolf_lion_DH", 
                                "lion_bear_DH", "lion_bob_DH", "coy_bob_DH")
  
  #####  Format covariate data  ####
  #'  ---------------------------
  #'  Format site-level covariates for detection sub-model
  det_covs <- stations_eoe20s21s %>%
    mutate(CameraFacing = as.factor(CameraFacing),
           Setup = ifelse(Setup == "ungulate", 0, 1),
           Height = as.numeric(Height))
  table(det_covs[,"CameraFacing"])
  table(det_covs[,"Setup"])
  
  ######  First order occupancy (psi|no second spp)  ####
  psi_cov <- matrix(NA, ncol = 16, nrow = nsites)
  psi_cov[,1] <- 1
  psi_cov[,2] <- det_covs$Setup
  psi_cov[,3] <- stations_eoe20s21s$Elev
  psi_cov[,4] <- stations_eoe20s21s$PercForest
  psi_cov[,5] <- stations_eoe20s21s$Nsmall_deer  #DomPrey
  psi_cov[,6] <- stations_eoe20s21s$Nbig_deer #SppDiversity
  psi_cov[,7] <- stations_eoe20s21s$Nelk
  psi_cov[,8] <- stations_eoe20s21s$Nmoose
  psi_cov[,9] <- stations_eoe20s21s$Nmd
  psi_cov[,10] <- stations_eoe20s21s$Nwtd
  psi_cov[,11] <- stations_eoe20s21s$Nlagomorph
  psi_cov[,12] <- stations_eoe20s21s$Dist2Burbs
  psi_cov[,13] <- stations_eoe20s21s$logNearestRd
  psi_cov[,14] <- stations_eoe20s21s$Nhuman
  psi_cov[,15] <- stations_eoe20s21s$Nlivestock
  psi_cov[,16] <- stations_eoe20s21s$NewLocationID
  head(psi_cov)

  ######  Second order occupancy (psi): 2-way interactions  ####
  psi_inxs_cov <- matrix(NA, ncol = 16, nrow = nsites)
  psi_inxs_cov[,1] <- 1
  psi_inxs_cov[,2] <- det_covs$Setup
  psi_inxs_cov[,3] <- stations_eoe20s21s$Elev
  psi_inxs_cov[,4] <- stations_eoe20s21s$PercForest
  psi_inxs_cov[,5] <- stations_eoe20s21s$Nsmall_deer  #DomPrey
  psi_inxs_cov[,6] <- stations_eoe20s21s$Nbig_deer #SppDiversity
  psi_inxs_cov[,7] <- stations_eoe20s21s$Nelk
  psi_inxs_cov[,8] <- stations_eoe20s21s$Nmoose
  psi_inxs_cov[,9] <- stations_eoe20s21s$Nmd
  psi_inxs_cov[,10] <- stations_eoe20s21s$Nwtd
  psi_inxs_cov[,11] <- stations_eoe20s21s$Nlagomorph
  psi_inxs_cov[,12] <- stations_eoe20s21s$Dist2Burbs
  psi_inxs_cov[,13] <- stations_eoe20s21s$logNearestRd
  psi_inxs_cov[,14] <- stations_eoe20s21s$Nhuman
  psi_inxs_cov[,15] <- stations_eoe20s21s$Nlivestock
  psi_inxs_cov[,16] <- stations_eoe20s21s$NewLocationID
  head(psi_inxs_cov)
  
  ######  First order detection (rho|no second spp)  ####
  rho_cov <- array(NA, dim = c(nsites, nsurveys, 5)) # last digit is number of covariates + intercept
  rho_cov[,,1] <- 1
  rho_cov[,,2] <- det_covs$CameraFacing
  rho_cov[,,3] <- det_covs$Setup
  rho_cov[,,4] <- det_covs$Height
  rho_cov[,,5] <- effort_eoe20s21s
  head(rho_cov)
  
  ######  Second order detection (rho): 2-way interactions  ####
  rho_inxs_cov <- array(NA, dim = c(nsites, nsurveys, 5))
  rho_inxs_cov[,,1] <- 1
  rho_inxs_cov[,,2] <- det_covs$CameraFacing
  rho_inxs_cov[,,3] <- det_covs$Setup
  rho_inxs_cov[,,4] <- det_covs$Height
  rho_inxs_cov[,,5] <- effort_eoe20s21s
  head(rho_inxs_cov)
  
  
  #'  -----------------------------------
  ####  Data and MCMC settings for JAGS  ####
  #'  -----------------------------------
  #####  Bundle detection history and covariate data  ####
  #'  -------------------------------------------------
  #'  Function to bundle detection and covariate data
  bundle_data <- function(obs_array, psi_covs, psi_inxs, rho_covs, rho_inxs, 
                          sites, surveys, psi_1order, psi_2order, rho_1order, 
                          rho_2order, ncats, uniquesites) {
    #'  list all pieces of data together
    bundled <- list(y = obs_array, psi_cov = psi_covs, psi_inxs_cov = psi_inxs,
                    rho_cov = rho_covs, rho_inxs_cov = rho_inxs, nsites = sites,
                    nsurveys = surveys, nfirst_order_psi = ncol(psi_1order), 
                    nsecond_order_psi = ncol(psi_2order), 
                    nfirst_order_rho = dim(rho_1order)[3], 
                    nsecond_order_rho = dim(rho_2order)[3], ncat = ncats, #rho_2order
                    uniquesites = as.numeric(factor(uniquesites), levels = uniquesites))
    #'  Summarize to make sure it looks right
    str(bundled)
    return(bundled)
  }
  bundled_pred_list <- lapply(multi_spp_DH_list, bundle_data, psi_covs = psi_cov, 
                               psi_inxs = psi_inxs_cov, rho_covs = rho_cov, 
                               rho_inxs = rho_inxs_cov, sites = nsites, 
                               surveys = nsurveys, psi_1order = psi_cov, 
                               psi_2order = psi_inxs_cov, rho_1order = rho_cov, 
                               rho_2order = rho_inxs_cov, ncats = ncat, 
                               uniquesites = unique(psi_cov[,16])) 
  # save(bundled_pred_list, file = "./Data/MultiSpp_OccMod_Outputs/bundled_predator_data_list.RData")
  
  
  #####  Initial values for model  ####
  #'  -----------------------------
  #'  Naive occupancy for each species at each site (site x spp matrix)
  initial_z <- function(bundled_dh) {
    #'  Naive occupancy for each species at each site (site x spp matrix)
    zinit <- apply(bundled_dh, c(1, 3), sum, na.rm = TRUE)
    zinit[zinit > 1] <- 1
    #'  Collapse 2-species detection state into 4 categories
    zcat <- apply(zinit, 1, paste, collapse = "")
    zcat[zcat == "00"] <- 1
    zcat[zcat == "10"] <- 2
    zcat[zcat == "01"] <- 3
    zcat[zcat == "11"] <- 4
    #'  Make z numeric again
    zcat <- as.numeric(zcat)
    
    return(zcat)
  }
  zinits <- lapply(DH_array_list, initial_z)
  names(zinits) <- c("wolf_bear_zcat", "wolf_coy_zcat", "wolf_lion_zcat", 
                     "lion_bear_zcat", "lion_bob_zcat", "coy_bob_zcat")
  
  #####  Parameters monitored  ####
  #'  -------------------------
  params <- c("betaSpp1", "betaSpp2", "alphaSpp1", "alphaSpp2", "betaSpp12", 
              "alphaSpp12", "alphaSpp21","sigmaSpp1", "sigmaSpp2",
              "mean.psiSpp1", "mean.psiSpp2", "mean.pSpp1", "mean.pSpp2", "z") 
  
  #####  MCMC settings  ####
  #'  ------------------
  nc <- 3
  ni <- 100000 
  nb <- 50000 
  nt <- 5
  na <- 5000 
  
  #'  -------------------
  ####  RUN JAGS MODELS  ####
  #'  -------------------
  #'  For each predator pairing:
  #'    1. set initial values with correct detection data,
  #'    2. source and run each model in JAGS
  #'    3. visually inspect traceplots
  #'    4. review model summary and any parameters that didn't converge well
  #'    5. save results
  #'    6. model selection
  
  ####  Wolf-Lion Models  ####
  #'  ----------------------
  inits.wolf.lion <- function(){list(z = zinits[[3]])}
  
  #'  Null 1: psi(.), p(.); no interactions, includes random effect for site
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(rx)_p(.).R")
  # start.time = Sys.time()
  # wolf.lion.null1 <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
  #                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(rx)_p(.).txt",
  #                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  # end.time <- Sys.time(); (run.time <- end.time - start.time)
  # print(wolf.lion.null1$summary)
  # print(wolf.lion.null1$DIC)
  # which(wolf.lion.null1$summary[,"Rhat"] > 1.1)
  # mcmcplot(wolf.lion.null1$samples)
  # save(wolf.lion.null1, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(rx)_p(.).RData")
  start.time = Sys.time()
  wolf.lion.null1 <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(rx)_p(.).txt",
                          n.chains = nc, n.iter = ni, n.burnin = bn, n.thin = nt,
                          DIC = TRUE, jags.module = c("glm", "dic"))
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  wolf.lion.null1.sum <- wolf.lion.null1$BUGSoutput$summary
  print(wolf.lion.null1.sum, dig = 3)
  mcmcplot(wolf.lion.null1.sum)
  save(wolf.lion.null1, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(rx)_p(.).RData")
  
  #'  Null 2: psi(.), psix(.), p(.), px(.), includes random effect for site
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(rx)_px(.).R")
  start.time = Sys.time()
  wolf.lion.null2 <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(rx)_px(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null2$summary)
  print(wolf.lion.null2$DIC)
  which(wolf.lion.null2$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null2$samples)
  save(wolf.lion.null2, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(.rx)_px(.).RData")
  
  #'  Habitat model: psi & psix = setup, elevation, forest; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_habitat_rx)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.hab <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_habitat_rx)_px(setup_effort).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.hab$summary)
  print(wolf.lion.hab$DIC)
  which(wolf.lion.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.hab$samples)
  save(wolf.lion.hab, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(setup_habitat_rx)_px(setup_effort).RData")
  
  #'  Prey groups model: psi & psix = setup, elevation, forest, small deer, big deer, livestock; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_wolflion.R")
  start.time = Sys.time()
  wolf.lion.preygroup <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_wolflion.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preygroup$summary)
  print(wolf.lion.preygroup$DIC)
  which(wolf.lion.preygroup$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preygroup$samples)
  save(wolf.lion.preygroup, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(setup_preygroup_rx)_px(setup_effort)_wolflion.RData")
  
  #'  Prey diversity model: psi & psix = setup, elevation, forest, elk, moose, md, wtd, livestock; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_preydiversity_rx)_px(setup_effort)_wolflion_wolflion.R")
  start.time = Sys.time()
  wolf.lion.preydiversity <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_preydiversity_rx)_px(setup_effort)_wolflion.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preydiversity$summary)
  print(wolf.lion.preydiversity$DIC)
  which(wolf.lion.preydiversity$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preydiversity$samples)
  save(wolf.lion.preydiversity, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(setup_preydiversity_rx)_px(setup_effort)_wolflion.RData")
  
  #'  Anthropogenic model: psi & psix = setup, elevation, forest, Dist2Burbs, LogNearestRd, human; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_anthro_rx)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.anthro <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_anthro_rx)_px(setup_effort).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.anthro$summary)
  print(wolf.lion.anthro$DIC)
  which(wolf.lion.anthro$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.anthro$samples)
  save(wolf.lion.anthro, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(setup_anthro_rx)_px(setup_effort).RData")
  
  #'  Global model #1: psi & psix = setup, elevation, forest, small deer, big deer, livestock, lagomorph, Dist2Burbs, LogNearestRd, human; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(global1)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.global1 <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(global1)_px(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.global1$summary)
  print(wolf.lion.global1$DIC)
  which(wolf.lion.global1$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.global1$samples)
  save(wolf.lion.global1, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(global1)_px(setup_effort).RData")
  
  #'  Global model #2: psi & psix = setup, elevation, forest, elk, moose, md, wtd, livestock, lagomorph, Dist2Burbs, LogNearestRd, human; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(global2)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.global2 <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(global2)_px(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.global2$summary)
  print(wolf.lion.global2$DIC)
  which(wolf.lion.global2$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.global2$samples)
  save(wolf.lion.global2, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psix(global2)_px(setup_effort).RData")
  
  
  
  ####  Wolf-Coyote Models  ####
  #'  ----------------------
  inits.wolf.coy <- function(){list(z = zinits[[2]])}
  
  #'  Habitat model: psi & psix = setup, elevation, forest; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_habitat_rx)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.hab <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                                 "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_habitat_rx)_px(setup_effort).txt",
                                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab$summary)
  print(wofl.coy.hab$DIC)
  which(wolf.coy.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab$samples)
  save(wolf.coy.hab, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psix(setup_habitat_rx)_px(setup_effort).RData")
  
  #'  Prey groups model: psi & psix = setup, elevation, forest, small deer, big deer, livestock, lagomorph; px = setup, effort
  # source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp2.R")
  # psi & psix = setup, elevation, forest, dominant prey spp, species diversity, lagomorph; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp12.R")
  start.time = Sys.time()
  wolf.coy.preygroup <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                             # "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp2.txt",
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp12.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preygroup$summary)
  print(wolf.coy.preygroup$DIC)
  which(wolf.coy.preygroup$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preygroup$samples)
  save(wolf.coy.preygroup, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp2.RData")
  
  #'  Prey diversity model: psi & psix = setup, elevation, forest, elk, moose, md, wtd, livestock, lagomorph; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_preydiversity_rx)_px(setup_effort)_bunnySpp2.R")
  start.time = Sys.time()
  wolf.coy.preydiversity <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params, 
                                 "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_preydiversity_rx)_px(setup_effort)_bunnySpp2.txt",
                                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preydiversity$summary)
  print(wolf.coy.preydiversity$DIC)
  which(wolf.coy.preydiversity$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preydiversity$samples)
  save(wolf.coy.preydiversity, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psix(setup_preydiversity_rx)_px(setup_effort)_bunnySpp2.RData")
  
  #'  Anthropogenic model: psi & psix = setup, elevation, forest, Dist2Burbs, LogNearestRd, human; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_anthro_rx)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.anthro <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_anthro_rx)_px(setup_effort).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.anthro$summary)
  print(wolf.coy.anthro$DIC)
  which(wolf.coy.anthro$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.anthro$samples)
  save(wolf.coy.anthro, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psix(setup_anthro_rx)_px(setup_effort).RData")
  
  #'  Global model #1: psi & psix = setup, elevation, forest, Dist2Burbs, LogNearestRd, human, livestock, small deer, big deer, lagomorph; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(global1)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.global1 <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(global1)_px(setup_effort).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.global1$summary)
  print(wolf.coy.global1$DIC)
  which(wolf.coy.global1$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.global1$samples)
  save(wolf.coy.global1, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psix(global1)_px(setup_effort).RData")
  
  #'  Global model #2: psi & psix = setup, elevation, forest, Dist2Burbs, LogNearestRd, human, livestock, elk, moose, md, wtd, lagomorph; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(global2)_px(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.global2 <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(global2)_px(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.global2$summary)
  print(wolf.coy.global2$DIC)
  which(wolf.coy.global2$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.global2$samples)
  save(wolf.coy.global2, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psix(global2)_px(setup_effort).RData")
  
  
  
  ####  Coyote-Bobcat Models  ####
  #'  ------------------------
  inits.coy.bob <- function(){list(z = zinits[[6]])}
  
  #'  Habitat model: psi = setup, elevation, forest; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_habitat_rx)_px(setup_effort).R")
  start.time = Sys.time()
  coy.bob.hab <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_habitat_rx)_px(setup_effort).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab$summary)
  which(coy.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab$samples)
  save(coy.bob.hab, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psix(setup_habitat_rx)_px(setup_effort).RData") 
  
  #'  Prey groups model: psi & psix = setup, elevation, forest, small deer, big deer, lagomorph; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp12.R")
  start.time = Sys.time()
  coy.bob.preygroup <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp12.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preygroup$summary)
  which(coy.bob.preygroup$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preygroup$samples)
  save(coy.bob.preygroup, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psix(setup_preygroup_rx)_px(setup_effort)_bunnySpp12.RData")
  
  #'  Prey diversity model: psi & psix = setup, elevation, forest, elk, moose, md, wtd, lagomorph; p = setup, effort, and interaction term
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_preydiversity_rx)_px(setup_effort)_bunnySpp12.R")
  start.time = Sys.time()
  coy.bob.preydiversity <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_preydiversity_rx)_px(setup_effort)_bunnySpp12.txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preydiversity$summary)
  which(coy.bob.preydiversity$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preydiversity$samples)
  save(coy.bob.preydiversity, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psix(setup_preydiversity_rx)_px(setup_effort)_bunnySpp12.RData")
  
  #'  Anthropogenic model: psi & psix = setup, elevation, forest, Dist2Burbs, LogNearestRd, human; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(setup_anthro_rx)_px(setup_effort).R")
  start.time = Sys.time()
  coy.bob.anthro <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(setup_anthro_rx)_px(setup_effort).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.anthro$summary)
  which(coy.bob.anthro$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.anthro$samples)
  save(coy.bob.anthro, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psix(setup_anthro_rx)_px(setup_effort).RData")
  
  #'  Global model #1: psi & psix = setup, elevation, forest, Dist2Burbs, LogNearestRd, human, livestock, small deer, big deer, lagomorph; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(global1)_px(setup_effort).R")
  start.time = Sys.time()
  coy.bob.global1 <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(global1)_px(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.global1$summary)
  which(coy.bob.global1$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.global1$samples)
  save(coy.bob.global1, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psix(global1)_px(setup_effort).RData")
  
  #'  Global model #2: psi & psix = setup, elevation, forest, Dist2Burbs, LogNearestRd, human, livestock, elk, moose, md, wtd, lagomorph; px = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psix(global2)_px(setup_effort).R")
  start.time = Sys.time()
  coy.bob.global2 <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psix(global2)_px(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.global2$summary)
  which(coy.bob.global2$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.global2$samples)
  save(coy.bob.global2, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psix(global2)_px(setup_effort).RData")
  