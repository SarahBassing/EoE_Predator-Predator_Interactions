  #'  ---------------------------------------
  #'  Bayesian multispecies occupancy model - separate years
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
  
  library(jagsUI)
  library(abind)
  library(mcmcplots)
  library(tidyverse)
  
  #' #'  Source detection histories and covariate data
  #' source("./Scripts/MultiSpp_OccMod/Detection_histories_for_occmod.R")
  #' source("./Scripts/MultiSpp_OccMod/Format_data_2spp_occmod_for_JAGS.R")
  
  #'  Load covariate and detection history data
  load("./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_sepYr_img.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")
  
  
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
    
    dh <- list(dh1, dh2)
    
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
  wolf_bear_DH_20s <- detection_array(spp1 = DH_wolf[[1]], spp2 = DH_bear[[1]], name1 = "wolf", name2 = "bear")
  wolf_bear_DH_21s <- detection_array(spp1 = DH_wolf[[2]], spp2 = DH_bear[[2]], name1 = "wolf", name2 = "bear")
  wolf_coy_DH_20s <- detection_array(spp1 = DH_wolf[[1]], spp2 = DH_coy[[1]], name1 = "wolf", name2 = "coyote")
  wolf_coy_DH_21s <- detection_array(spp1 = DH_wolf[[2]], spp2 = DH_coy[[2]], name1 = "wolf", name2 = "coyote")
  wolf_lion_DH_20s <- detection_array(spp1 = DH_wolf[[1]], spp2 = DH_lion[[1]], name1 = "wolf", name2 = "lion")
  wolf_lion_DH_21s <- detection_array(spp1 = DH_wolf[[2]], spp2 = DH_lion[[2]], name1 = "wolf", name2 = "lion")
  lion_bear_DH_20s <- detection_array(spp1 = DH_lion[[1]], spp2 = DH_bear[[1]], name1 = "lion", name2 = "bear")
  lion_bear_DH_21s <- detection_array(spp1 = DH_lion[[2]], spp2 = DH_bear[[2]], name1 = "lion", name2 = "bear")
  lion_bob_DH_20s <- detection_array(spp1 = DH_lion[[1]], spp2 = DH_bob[[1]], name1 = "lion", name2 = "bobcat")
  lion_bob_DH_21s <- detection_array(spp1 = DH_lion[[2]], spp2 = DH_bob[[2]], name1 = "lion", name2 = "bobcat")
  coy_bob_DH_20s <- detection_array(spp1 = DH_coy[[1]], spp2 = DH_bob[[1]], name1 = "coyote", name2 = "bobcat")
  coy_bob_DH_21s <- detection_array(spp1 = DH_coy[[2]], spp2 = DH_bob[[2]], name1 = "coyote", name2 = "bobcat")
  
  #'  List 2-species detection arrays together for faster formatting below
  DH_array_list_20s <- list(wolf_bear_DH_20s, wolf_coy_DH_20s, wolf_lion_DH_20s, lion_bear_DH_20s, lion_bob_DH_20s, coy_bob_DH_20s)
  DH_array_list_21s <- list(wolf_bear_DH_21s, wolf_coy_DH_21s, wolf_lion_DH_21s, lion_bear_DH_21s, lion_bob_DH_21s, coy_bob_DH_21s)
  
  
  #'  -------------------------
  ####  Format data for JAGS  ####
  #'  -------------------------
  #'  Define survey dimensions
  nsites_20s <- dim(wolf_bear_DH_20s)[1]
  nsites_21s <- dim(wolf_bear_DH_21s)[1]
  nsurveys <- dim(wolf_bear_DH_20s)[2]
  nspecies <- dim(wolf_bear_DH_20s)[3]
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
  multi_spp_DH_20s_list <- lapply(DH_array_list_20s, observation_state)
  multi_spp_DH_21s_list <- lapply(DH_array_list_21s, observation_state)
  spp_names <- c("wolf_bear_DH", "wolf_coy_DH", "wolf_lion_DH", "lion_bear_DH", "lion_bob_DH", "coy_bob_DH")
  names(multi_spp_DH_20s_list) <- spp_names
  names(multi_spp_DH_21s_list) <- spp_names
  
  #####  Format covariate data  ####
  #'  ---------------------------
  #'  Format site-level covariates for detection sub-model
  format_det_covs <- function(stations) {
    det_covs <- stations %>%
      mutate(CameraFacing = as.factor(CameraFacing),
             Setup = ifelse(Setup == "ungulate", 0, 1),
             Height = as.numeric(Height))
    print(table(det_covs[,"CameraFacing"]))
    print(table(det_covs[,"Setup"]))
    return(det_covs)
  }
  det_covs <- lapply(stations_eoe20s21s, format_det_covs)
 
  #'  Format first and second-order covariates for each year
  format_natural_param_covs <- function(nsites, nsurveys, det, stations, effort) {
    
    ######  First order occupancy (psi|no second spp)  ####
    psi_cov <- matrix(NA, ncol = 16, nrow = nsites)
    psi_cov[,1] <- 1
    psi_cov[,2] <- det$Setup
    psi_cov[,3] <- stations$Elev
    psi_cov[,4] <- stations$PercForest
    psi_cov[,5] <- factor(stations$DomPrey, levels = c("elk", "whitetaileddeer", "other"))
    psi_cov[,6] <- stations$SppDiversity
    psi_cov[,7] <- stations$Nelk
    psi_cov[,8] <- stations$Nmoose
    psi_cov[,9] <- stations$Nmd
    psi_cov[,10] <- stations$Nwtd
    psi_cov[,11] <- stations$Nlagomorph
    psi_cov[,12] <- stations$Dist2Burbs
    psi_cov[,13] <- stations$logNearestRd
    psi_cov[,14] <- stations$Nhuman
    psi_cov[,15] <- stations$Nlivestock
    psi_cov[,16] <- stations$NewLocationID
    print(head(psi_cov))
    
    ######  Second order occupancy (psi): 2-way interactions  ####
    psi_inxs_cov <- matrix(NA, ncol = 16, nrow = nsites)
    psi_inxs_cov[,1] <- 1
    psi_inxs_cov[,2] <- det$Setup
    psi_inxs_cov[,3] <- stations$Elev
    psi_inxs_cov[,4] <- stations$PercForest
    psi_inxs_cov[,5] <- factor(stations$DomPrey, levels = c("elk", "whitetaileddeer", "other"))
    psi_inxs_cov[,6] <- stations$SppDiversity
    psi_inxs_cov[,7] <- stations$Nelk
    psi_inxs_cov[,8] <- stations$Nmoose
    psi_inxs_cov[,9] <- stations$Nmd
    psi_inxs_cov[,10] <- stations$Nwtd
    psi_inxs_cov[,11] <- stations$Nlagomorph
    psi_inxs_cov[,12] <- stations$Dist2Burbs
    psi_inxs_cov[,13] <- stations$logNearestRd
    psi_inxs_cov[,14] <- stations$Nhuman
    psi_inxs_cov[,15] <- stations$Nlivestock
    psi_inxs_cov[,16] <- stations$NewLocationID
    print(head(psi_inxs_cov))
    
    ######  First order detection (rho|no second spp)  ####
    rho_cov <- array(NA, dim = c(nsites, nsurveys, 5)) # last digit is number of covariates + intercept
    rho_cov[,,1] <- 1
    rho_cov[,,2] <- det$CameraFacing
    rho_cov[,,3] <- det$Setup
    rho_cov[,,4] <- det$Height
    rho_cov[,,5] <- effort
    print(head(rho_cov))
    
    ######  Second order detection (rho): 2-way interactions  ####
    rho_inxs_cov <- array(NA, dim = c(nsites, nsurveys, 5))
    rho_inxs_cov[,,1] <- 1
    rho_inxs_cov[,,2] <- det$CameraFacing
    rho_inxs_cov[,,3] <- det$Setup
    rho_inxs_cov[,,4] <- det$Height
    rho_inxs_cov[,,5] <- effort
    print(head(rho_inxs_cov))
    
    #'  List matrices together
    cov_matrix_list <- list(psi_cov, psi_inxs_cov, rho_cov, rho_inxs_cov)
    return(cov_matrix_list)
  }
  cov_matrices_20s <- format_natural_param_covs(nsites = nsites_20s, nsurveys = nsurveys, 
                                                det = det_covs[[1]], stations = stations_eoe20s21s[[1]], 
                                                effort = effort_eoe20s21s[[1]])
  cov_matrices_21s <- format_natural_param_covs(nsites = nsites_21s, nsurveys = nsurveys, 
                                                det = det_covs[[2]], stations = stations_eoe20s21s[[2]], 
                                                effort = effort_eoe20s21s[[2]])
  
  
  #'  -----------------------------------
  ####  Data and MCMC settings for JAGS  ####
  #'  -----------------------------------
  #####  Bundle detection history and covariate data  ####
  #'  -------------------------------------------------
  #'  Function to bundle detection and covariate data
  bundle_data <- function(obs_array, psi_covs, psi_inxs, rho_covs, rho_inxs, 
                          sites, surveys, psi_1order, psi_2order, rho_1order, 
                          rho_2order, ncats, nspecies, uniquesites) {
    #'  list all pieces of data together
    bundled <- list(y = obs_array, psi_cov = psi_covs, psi_inxs_cov = psi_inxs,
                    rho_cov = rho_covs, rho_inxs_cov = rho_inxs, nsites = sites,
                    nsurveys = surveys, nfirst_order_psi = ncol(psi_1order), 
                    nsecond_order_psi = ncol(psi_2order), 
                    nfirst_order_rho = dim(rho_1order)[3], 
                    nsecond_order_rho = dim(rho_2order)[3], ncat = ncats, 
                    nspec = nspecies, uniquesites = as.numeric(factor(uniquesites), levels = uniquesites)) #rho_2order
    #'  Summarize to make sure it looks right
    str(bundled)
    return(bundled)
  }
  bundled_pred_list_20s <- lapply(multi_spp_DH_20s_list, bundle_data, psi_covs = cov_matrices_20s[[1]], 
                              psi_inxs = cov_matrices_20s[[2]], rho_covs = cov_matrices_20s[[3]], 
                              rho_inxs = cov_matrices_20s[[4]], sites = nsites_20s, 
                              surveys = nsurveys, psi_1order = cov_matrices_20s[[1]], 
                              psi_2order = cov_matrices_20s[[2]], rho_1order = cov_matrices_20s[[3]], 
                              rho_2order = cov_matrices_20s[[4]], ncats = ncat, 
                              nspecies = nspecies, uniquesites = unique(cov_matrices_20s[[1]][,16])) 
  bundled_pred_list_21s <- lapply(multi_spp_DH_21s_list, bundle_data, psi_covs = cov_matrices_21s[[1]], 
                              psi_inxs = cov_matrices_21s[[2]], rho_covs = cov_matrices_21s[[3]], 
                              rho_inxs = cov_matrices_21s[[4]], sites = nsites_21s, 
                              surveys = nsurveys, psi_1order = cov_matrices_21s[[1]], 
                              psi_2order = cov_matrices_21s[[2]], rho_1order = cov_matrices_21s[[3]], 
                              rho_2order = cov_matrices_21s[[4]], ncats = ncat, 
                              nspecies = nspecies, uniquesites = unique(cov_matrices_21s[[1]][,16])) 
  # save(bundled_pred_list_20s, file = "./Data/MultiSpp_OccMod_Outputs/bundled_predator_data_list_20s.RData")
  # save(bundled_pred_list_21s, file = "./Data/MultiSpp_OccMod_Outputs/bundled_predator_data_list_21s.RData")
  
  
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
  zinits_20s <- lapply(DH_array_list_20s, initial_z)
  zinits_21s <- lapply(DH_array_list_21s, initial_z)
  zcat_names <- c("wolf_bear_zcat", "wolf_coy_zcat", "wolf_lion_zcat", 
                  "lion_bear_zcat", "lion_bob_zcat", "coy_bob_zcat")
  names(zinits_20s) <- zcat_names
  names(zinits_21s) <- zcat_names
  
  #####  Parameters monitored  ####
  #'  -------------------------
  params <- c("betaSpp1", "betaSpp2", "alphaSpp1", "alphaSpp2", "betaSpp12", 
              "alphaSpp12", "alphaSpp21", "mean.psiSpp1", "mean.psiSpp2", 
              "mean.pSpp1", "mean.pSpp2", "z") 
  
  #####  MCMC settings  ####
  #'  ------------------
  # nc <- 3
  # ni <- 100000
  # nb <- 50000  #75000
  # nt <- 10 #5
  # na <- 5000
  nc <- 3
  ni <- 20000
  nb <- 10000
  nt <- 5
  na <- 1000
  
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
  
  
  ####  Wolf-Bear Models  ####
  #'  ----------------------
  #' #' Provide somewhat informed starting values for intercepts & sigmas
  #' inits.wolf.bear <- function(){list(z = zinits[[1]], mean.psiSpp1 = runif(1,0,0.2),
  #'                                    mean.psiSpp2 = runif(1,0.5,0.9), mean.pSpp1 = runif(1,0,0.2),
  #'                                    mean.pSpp2 = runif(1,0,0.3), sigmaSpp1 = runif(1,1.1,1.9), sigmaSpp2 = runif(1,1,2))}
  #' #'  Attempting to identify if there are 2+ local optima by providing sigma1 different initial values
  #' inits.wolf.bear_sigma1.98 <- function(){list(z = zinits[[1]], mean.psiSpp1 = runif(1,0,0.2),
  #'                                    mean.psiSpp2 = runif(1,0.5,0.9), mean.pSpp1 = runif(1,0,0.2),
  #'                                    mean.pSpp2 = runif(1,0,0.3), sigmaSpp1 = 1.98, sigmaSpp2 = runif(1,1,2))}
  #' inits.wolf.bear_sigma1.55 <- function(){list(z = zinits[[1]], mean.psiSpp1 = runif(1,0,0.2),
  #'                                    mean.psiSpp2 = runif(1,0.5,0.9), mean.pSpp1 = runif(1,0,0.2),
  #'                                    mean.pSpp2 = runif(1,0,0.3), sigmaSpp1 = 1.55, sigmaSpp2 = runif(1,1,2))}
  inits.wolf.bear_20s <- function(){list(z = zinits_20s[[1]])}
  inits.wolf.bear_21s <- function(){list(z = zinits_21s[[1]])}
  
  #####  Null model  ####
  #'  psi = random effect
  run_wolfbear_nullmod <- function(bundled, inits) {
    source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(.)_p(.).R")
    start.time = Sys.time()
    wolf.bear.null <- jags(bundled, inits = inits, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(.)_p(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
    end.time <- Sys.time(); print(run.time <- end.time - start.time)
    print(wolf.bear.null$summary)
    print(wolf.bear.null$DIC)
    which(wolf.bear.null$summary[,"Rhat"] > 1.1)
    mcmcplot(wolf.bear.null$samples)
    return(wolf.bear.null)
  }
  wolf.bear.null.20s <- run_wolfbear_nullmod(bundled = bundled_pred_list_20s[[1]], inits = inits.wolf.bear_20s)
  save(wolf.bear.null.20s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear20s_psi(.)_p(.).RData")
  wolf.bear.null.21s <- run_wolfbear_nullmod(bundled = bundled_pred_list_21s[[1]], inits = inits.wolf.bear_21s)
  save(wolf.bear.null.21s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear21s_psi(.)_p(.).RData")
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  run_wolfbear_habmod <- function(bundled, inits) {
    source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat)_p(setup_effort).R")
    start.time = Sys.time()
    wolf.bear.hab <- jags(bundled, inits, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat)_p(setup_effort).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
    end.time <- Sys.time(); print(run.time <- end.time - start.time)
    print(wolf.bear.hab$summary)
    print(wolf.bear.hab$DIC)
    which(wolf.bear.hab$summary[,"Rhat"] > 1.1)
    mcmcplot(wolf.bear.hab$samples)
    return(wolf.bear.hab)
  }
  wolf.bear.hab.20s <- run_wolfbear_habmod(bundled = bundled_pred_list_20s[[1]], inits = inits.wolf.bear_20s)
  save(wolf.bear.hab.20s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear20s_psi(setup_habitat)_p(setup_effort).RData")
  wolf.bear.hab.21s <- run_wolfbear_habmod(bundled = bundled_pred_list_21s[[1]], inits = inits.wolf.bear_21s)
  save(wolf.bear.hab.21s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear21s_psi(setup_habitat)_p(setup_effort).RData")
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, elevation, forest, elk, moose, wtd; p = setup, effort  
  run_wolfbear_npreymod <- function(bundled, inits) {
    source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund)_p(setup_effort)_wolfbearlion.R")
    start.time = Sys.time()
    wolf.bear.preyabund <- jags(bundled, inits = inits, params,
                                "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund)_p(setup_effort)_wolfbearlion.txt",
                                n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
    end.time <- Sys.time(); print(run.time <- end.time - start.time)
    print(wolf.bear.preyabund$summary)
    print(wolf.bear.preyabund$DIC)
    which(wolf.bear.preyabund$summary[,"Rhat"] > 1.1)
    mcmcplot(wolf.bear.preyabund$samples)
    return(wolf.bear.preyabund)
  }
  wolf.bear.preyabund.20s <- run_wolfbear_npreymod(bundled = bundled_pred_list_20s[[1]], inits = inits.wolf.bear_20s)
  save(wolf.bear.preyabund.20s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear20s_psi(setup_preyabund)_p(setup_effort).RData")
  wolf.bear.preyabund.21s <- run_wolfbear_npreymod(bundled = bundled_pred_list_21s[[1]], inits = inits.wolf.bear_21s)
  save(wolf.bear.preyabund.21s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear21s_psi(setup_preyabund)_p(setup_effort).RData")
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, elevation, forest, spp diversity; p = setup, effort  
  run_wolfbear_pdivmod <- function(bundled, inits) {
    source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity)_p(setup_effort).R")
    start.time = Sys.time()
    wolf.bear.preydiv <- jags(bundled, inits = inits, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity)_p(setup_effort).txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
    end.time <- Sys.time(); print(run.time <- end.time - start.time)
    print(wolf.bear.preydiv$summary)
    print(wolf.bear.preydiv$DIC)
    which(wolf.bear.preydiv$summary[,"Rhat"] > 1.1)
    mcmcplot(wolf.bear.preydiv$samples)
    return(wolf.bear.preydiv)
  }
  wolf.bear.preydiv.20s <- run_wolfbear_pdivmod(bundled = bundled_pred_list_20s[[1]], inits = inits.wolf.bear_20s)
  save(wolf.bear.preydiv.20s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear20s_psi(setup_preydiversity)_p(setup_effort).RData")
  wolf.bear.preydiv.21s <- run_wolfbear_pdivmod(bundled = bundled_pred_list_21s[[1]], inits = inits.wolf.bear_21s)
  save(wolf.bear.preydiv.21s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear21s_psi(setup_preydiversity)_p(setup_effort).RData")
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, elevation, forest; psix(.); p = setup, effort  
  run_wolfbear_habxmod <- function(bundled, inits) {
    source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat)_psix(.)_p(setup_effort).R")
    start.time = Sys.time()
    wolf.bear.habx <- jags(bundled, inits = inits, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat)_psix(.)_p(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
    end.time <- Sys.time(); print(run.time <- end.time - start.time)
    print(wolf.bear.habx$summary)
    print(wolf.bear.habx$DIC)
    which(wolf.bear.habx$summary[,"Rhat"] > 1.1)
    mcmcplot(wolf.bear.habx$samples)
    return(wolf.bear.habx)
  }
  wolf.bear.habx.20s <- run_wolfbear_habxmod(bundled = bundled_pred_list_20s[[1]], inits = inits.wolf.bear_20s)
  save(wolf.bear.habx.20s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear20s_psi(setup_habitat)_psix(.)_p(setup_effort).RData")
  wolf.bear.habx.21s <- run_wolfbear_habxmod(bundled = bundled_pred_list_21s[[1]], inits = inits.wolf.bear_21s)
  save(wolf.bear.habx.21s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear21s_psi(setup_habitat)_psix(.)_p(setup_effort).RData")
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, elevation, forest; psix(elk, moose, wtd); p = setup, effort  
  run_wolfbear_npreyxmod <- function(bundled, inits) {
    source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
    start.time = Sys.time()
    wolf.bear.preyabundx <- jags(bundled, inits = inits, params,
                                 "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
    end.time <- Sys.time(); print(run.time <- end.time - start.time)
    print(wolf.bear.preyabundx$summary)
    print(wolf.bear.preyabundx$DIC)
    which(wolf.bear.preyabundx$summary[,"Rhat"] > 1.1)
    mcmcplot(wolf.bear.preyabundx$samples)
    return(wolf.bear.preyabundx)
  }
  wolf.bear.preyabundx.20s <- run_wolfbear_npreyxmod(bundled = bundled_pred_list_20s[[1]], inits = inits.wolf.bear_20s)
  save(wolf.bear.preyabundx.20s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear20s_psi(setup_habitat)_psix(preyabund)_p(setup_effort).RData")
  wolf.bear.preyabundx.21s <- run_wolfbear_npreyxmod(bundled = bundled_pred_list_21s[[1]], inits = inits.wolf.bear_21s)
  save(wolf.bear.preyabundx.21s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear21s_psi(setup_habitat)_psix(preyabund)_p(setup_effort).RData")
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, elevation, forest; psix(spp diversity); p = setup, effort  
  run_wolfbear_pdivxmod <- function(bundled, inits) {
    source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat)_psix(preydiversity)_p(setup_effort).R")
    start.time = Sys.time()
    wolf.bear.preydivx <- jags(bundled, inits = inits, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat)_psix(preydiversity)_p(setup_effort).txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
    end.time <- Sys.time(); print(run.time <- end.time - start.time)
    print(wolf.bear.preydivx$summary)
    print(wolf.bear.preydivx$DIC)
    which(wolf.bear.preydivx$summary[,"Rhat"] > 1.1)
    mcmcplot(wolf.bear.preydivx$samples)
    return(wolf.bear.preydivx)
  }
  wolf.bear.preydivx.20s <- run_wolfbear_pdivxmod(bundled = bundled_pred_list_20s[[1]], inits = inits.wolf.bear_20s)
  save(wolf.bear.preydivx.20s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear20s_psi(setup_habitat)_psix(preydiversity)_p(setup_effort).RData")
  wolf.bear.preydivx.21s <- run_wolfbear_pdivxmod(bundled = bundled_pred_list_21s[[1]], inits = inits.wolf.bear_21s)
  save(wolf.bear.preydivx.21s, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear21s_psi(setup_habitat)_psix(preydiversity)_p(setup_effort).RData")
  
  
  #'  ----------------------
  ####  Wolf-Coyote Models  ####
  #'  ----------------------
  #' Provide somewhat informed starting values for each parameter
  inits.wolf.coy <- function(){list(z = zinits[[2]], mean.psiSpp1 = runif(1,0,0.1),
                                    mean.psiSpp2 = runif(1,0.4,0.7), mean.pSpp1 = runif(1,0,0.2),
                                    mean.pSpp2 = runif(1,0.3,0.5), sigmaSpp1 = runif(1,2,3), 
                                    sigmaSpp2 = runif(1,4.5,5.5))}
  
  inits.wolf.coy <- function(){list(z = zinits[[2]], mean.psiSpp1 = runif(1),
                                    mean.psiSpp2 = runif(1), mean.pSpp1 = runif(1),
                                    mean.pSpp2 = runif(1), sigmaSpp1 = runif(1), 
                                    sigmaSpp2 = runif(1))}
  
  #' #' Provide somewhat informed starting values for each parameter based on alternative parameterization
  #' inits.wolf.coy <- function(){list(z = zinits[[2]], betaSpp1 = runif(1,-4,-3), betaSpp1[2] = runif(1,2,5), 
  #'                                   betaSpp1[3] = runif(1,-1,1), betaSpp1[4] = runif(1,-1,1),
  #'                                   betaSpp2[1] = runif(1,0,-2), betaSpp2[2] = runif(1,2,5),
  #'                                   betaSpp2[3] = runif(1,-1,1), betaSpp2[4] = runif(1,-1,1),
  #'                                   alphaSpp1 = runif(3,-2,2), alphaSpp2 = runif(3,-2,2),
  #'                                   sigmaSpp1 = runif(1,1,2), sigmaSpp2 = runif(1,2,3))}
  
  #####  Null model  ####
  #'  psi = random effect
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(rx)_p(.).R")
  start.time = Sys.time()
  wolf.coy.null <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(rx)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.null$summary)
  print(wolf.coy.null$DIC)
  which(wolf.coy.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.null$samples)
  save(wolf.coy.null, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(rx)_p(.).RData")
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.hab <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab$summary)
  print(wolf.coy.hab$DIC)
  which(wolf.coy.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab$samples)
  save(wolf.coy.hab, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_rx)_p(setup_effort).RData")
  
  #####  Habitat no inxs model  ####
  #'  psi = setup, elevation, forest; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.hab <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab$summary)
  print(wolf.coy.hab$DIC)
  which(wolf.coy.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab$samples)
  save(wolf.coy.hab, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort).RData")
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, elevation, forest, elk, moose, wtd, livestock; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabund <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_wolfcoy.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabund$summary)
  print(wolf.coy.preyabund$DIC)
  which(wolf.coy.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabund$samples)
  save(wolf.coy.preyabund, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preyabund_rx)_p(setup_effort).RData")
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, elevation, forest, spp diversity; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.preydiv <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preydiv$summary)
  print(wolf.coy.preydiv$DIC)
  which(wolf.coy.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preydiv$samples)
  save(wolf.coy.preydiv, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preydiversity_rx)_p(setup_effort).RData")
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, elevation, forest; psix(.); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.habx <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.habx$summary)
  print(wolf.coy.habx$DIC)
  which(wolf.coy.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.habx$samples)
  save(wolf.coy.habx, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_rx)_psix(.)_p(setup_effort).RData")
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, elevation, forest; psix(elk, moose, wtd, lagomorph, livestock); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabundx <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_wolfcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabundx$summary)
  print(wolf.coy.preyabundx$DIC)
  which(wolf.coy.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabundx$samples)
  save(wolf.coy.preyabundx, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort).RData")
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, elevation, forest; psix(spp diversity); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.preydivx <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preydivx$summary)
  print(wolf.coy.preydivx$DIC)
  which(wolf.coy.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preydivx$samples)
  save(wolf.coy.preydivx, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).RData")
  
  
  #'  ---------------------
  ####  Wolf-Lion Models  ####
  #'  ---------------------
  inits.wolf.lion <- function(){list(z = zinits[[3]], mean.psiSpp1 = runif(1,0,0.3), # 0.08 vs 0.19 local optima for mean.psiSpp1
                                     mean.psiSpp2 = runif(1,0.1,0.4), mean.pSpp1 = runif(1,0,0.2), # 0.16 vs 0.28 local optima for mean.psiSpp2
                                     mean.pSpp2 = runif(1,0.1,0.2), sigmaSpp1 = runif(1,2,3), sigmaSpp2 = runif(1,4.5,5.5))}
  # 0.25 vs 2.5 local optima for sigmaSpp1; returns prior for sigmaSpp2
  
  #####  Null model  ####
  #'  psi = random effect
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(rx)_p(.).R")
  start.time = Sys.time()
  wolf.lion.null <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(rx)_p(.).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null$summary)
  print(wolf.lion.null$DIC)
  which(wolf.lion.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null$samples)
  save(wolf.lion.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(rx)_p(.)_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.hab <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.hab$summary)
  print(wolf.lion.hab$DIC)
  which(wolf.lion.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.hab$samples)
  save(wolf.lion.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, elevation, forest, elk, moose, wtd, livestock; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabund <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_wolfbearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabund$summary)
  print(wolf.lion.preyabund$DIC)
  which(wolf.lion.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabund$samples)
  save(wolf.lion.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preyabund_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, elevation, forest, spp diversity; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.preydiv <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preydiv$summary)
  print(wolf.lion.preydiv$DIC)
  which(wolf.lion.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preydiv$samples)
  save(wolf.lion.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preydiversity_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, elevation, forest; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.habx <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.habx$summary)
  print(wolf.lion.habx$DIC)
  which(wolf.lion.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.habx$samples)
  save(wolf.lion.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_rx)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, elevation, forest; psix(elk, moose, wtd, livestock); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabundx <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabundx$summary)
  print(wolf.lion.preyabundx$DIC)
  which(wolf.lion.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabundx$samples)
  save(wolf.lion.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, elevation, forest; psix(spp diversity); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.preydivx <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preydivx$summary)
  print(wolf.lion.preydivx$DIC)
  which(wolf.lion.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preydivx$samples)
  save(wolf.lion.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  
  #'  --------------------
  ####  Lion-Bear Models  ####
  #'  --------------------
  inits.lion.bear <- function(){list(z = zinits[[4]], mean.psiSpp1 = runif(1,0.1,0.3),
                                     mean.psiSpp2 = runif(1,0.6,0.8), mean.pSpp1 = runif(1,0.05,0.15), 
                                     mean.pSpp2 = runif(1,0.1,0.2), sigmaSpp1 = runif(1,1,1.5), sigmaSpp2 = runif(1,1,3))}
  inits.lion.bear <- function(){list(z = zinits[[4]], mean.psiSpp1 = runif(1),
                                     mean.psiSpp2 = runif(1), mean.pSpp1 = runif(1), 
                                     mean.pSpp2 = runif(1), sigmaSpp1 = runif(1), sigmaSpp2 = runif(1))}
  
  
  #####  Null model  ####
  #'  psi = random effect
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(rx)_p(.).R")
  start.time = Sys.time()
  lion.bear.null <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(rx)_p(.).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null$summary)
  print(lion.bear.null$DIC)
  which(lion.bear.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null$samples)
  save(lion.bear.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(rx)_p(.)_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.hab <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.hab$summary)
  print(lion.bear.hab$DIC)
  which(lion.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.hab$samples)
  save(lion.bear.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat)_p(setup_effort_rx).R")
  start.time = Sys.time()
  lion.bear.hab <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat)_p(setup_effort_rx).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.hab$summary)
  print(lion.bear.hab$DIC)
  which(lion.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.hab$samples)
  save(lion.bear.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat)_p(setup_effort_rx)_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.hab <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.hab$summary)
  print(lion.bear.hab$DIC)
  which(lion.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.hab$samples)
  save(lion.bear.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, elevation, forest, elk, wtd; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabund <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_bearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabund$summary)
  print(lion.bear.preyabund$DIC)
  which(lion.bear.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabund$samples)
  save(lion.bear.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preyabund_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, elevation, forest, spp diversity; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.preydiv <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preydiv$summary)
  print(lion.bear.preydiv$DIC)
  which(lion.bear.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preydiv$samples)
  save(lion.bear.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preydiversity_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, elevation, forest; psix(.); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.habx <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.habx$summary)
  print(lion.bear.habx$DIC)
  which(lion.bear.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.habx$samples)
  save(lion.bear.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_rx)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, elevation, forest; psix(elk, wtd); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabundx <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_bearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabundx$summary)
  print(lion.bear.preyabundx$DIC)
  which(lion.bear.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabundx$samples)
  save(lion.bear.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, elevation, forest; psix(spp diversity); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.preydivx <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preydivx$summary)
  print(lion.bear.preydivx$DIC)
  which(lion.bear.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preydivx$samples)
  save(lion.bear.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  
  #'  ----------------------
  ####  Lion-Bobcat Models  ####
  #'  ----------------------
  inits.lion.bob <- function(){list(z = zinits[[5]], mean.psiSpp1 = runif(1,0.1,0.4), # 0.15 vs 0.27 local optima for mean.psiSpp1
                                    mean.psiSpp2 = runif(1,0,0.3), mean.pSpp1 = runif(1,0,0.2), # 0.12 vs 0.26 local optima for mean.psiSpp2
                                    mean.pSpp2 = runif(1,0.1,0.2), sigmaSpp1 = runif(1,1.5,2.5), sigmaSpp2 = runif(1,1,3))}
  # 0.25 vs 2.5 local optima for sigmaSpp1; returns prior for sigmaSpp2
  
  #####  Null model  ####
  #'  psi = random effect
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(rx)_p(.).R")
  start.time = Sys.time()
  lion.bob.null <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(rx)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null$summary)
  print(lion.bob.null$DIC)
  which(lion.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null$samples)
  save(lion.bob.null, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(rx)_p(.).RData")
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.hab <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.hab$summary)
  print(lion.bob.hab$DIC)
  which(lion.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.hab$samples)
  save(lion.bob.hab, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_rx)_p(setup_effort).RData")
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, elevation, forest, elk, wtd, lagomorphs; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.preyabund <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_lionbob.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabund$summary)
  print(lion.bob.preyabund$DIC)
  which(lion.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabund$samples)
  save(lion.bob.preyabund, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preyabund_rx)_p(setup_effort).RData")
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, elevation, forest, spp diversity; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.preydiv <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preydiv$summary)
  print(lion.bob.preydiv$DIC)
  which(lion.bob.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preydiv$samples)
  save(lion.bob.preydiv, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preydiversity_rx)_p(setup_effort).RData")
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, elevation, forest; psix(.); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.habx <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.habx$summary)
  print(lion.bob.habx$DIC)
  which(lion.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.habx$samples)
  save(lion.bob.habx, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_rx)_psix(.)_p(setup_effort).RData")
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, elevation, forest; psix(elk, wtd, lagomorphs); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.preyabundx <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_lionbob.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabundx$summary)
  print(lion.bob.preyabundx$DIC)
  which(lion.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabundx$samples)
  save(lion.bob.preyabundx, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort).RData")
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, elevation, forest; psix(spp diversity); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.preydivx <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preydivx$summary)
  print(lion.bob.preydivx$DIC)
  which(lion.bob.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preydivx$samples)
  save(lion.bob.preydivx, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).RData")
  
  
  #'  ------------------------
  ####  Coyote-Bobcat Models  ####
  #'  ------------------------
  # inits.coy.bob <- function(){list(z = zinits[[6]], mean.psiSpp1 = runif(1,0.4,0.6),
  #                                  mean.psiSpp2 = runif(1,0,0.3), mean.pSpp1 = runif(1,0.4,0.5), 
  #                                  mean.pSpp2 = runif(1,0.1,0.2), sigmaSpp1 = runif(1,1,2), sigmaSpp2 = runif(1,1,3))}
  inits.coy.bob <- function(){list(z = zinits[[6]], mean.psiSpp1 = runif(1),
                                   mean.psiSpp2 = runif(1), mean.pSpp1 = runif(1), 
                                   mean.pSpp2 = runif(1), sigmaSpp1 = runif(1), sigmaSpp2 = runif(1))}
  
  #####  Null model  ####
  #'  psi = random effect
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(rx)_p(.).R")
  start.time = Sys.time()
  coy.bob.null <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(rx)_p(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.null$summary)
  print(coy.bob.null$DIC)
  which(coy.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.null$samples)
  save(coy.bob.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(rx)_p(.)_", Sys.Date(), ".RData"))
  
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_p(setup_effort)_v2.R")
  start.time = Sys.time()
  coy.bob.hab <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_p(setup_effort)_v2.txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab$summary)
  print(coy.bob.hab$DIC)
  which(coy.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab$samples)
  save(coy.bob.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_rx)_p(setup_effort)_v2_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, elevation, forest; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.hab <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_p(setup_effort).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab$summary)
  print(coy.bob.hab$DIC)
  which(coy.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab$samples)
  save(coy.bob.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, elevation, forest, wtd, lagomorphs; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabund <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_rx)_p(setup_effort)_coybob.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabund$summary)
  print(coy.bob.preyabund$DIC)
  which(coy.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabund$samples)
  save(coy.bob.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preyabund_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, elevation, forest, spp diversity; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.preydiv <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_rx)_p(setup_effort).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preydiv$summary)
  print(coy.bob.preydiv$DIC)
  which(coy.bob.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preydiv$samples)
  save(coy.bob.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preydiversity_rx)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, elevation, forest; psix(.); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.habx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(.)_p(setup_effort).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx$summary)
  print(coy.bob.habx$DIC)
  which(coy.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.habx$samples)
  save(coy.bob.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_rx)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, elevation, forest; psix(wtd, lagomorphs); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabundx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_coybob.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabundx$summary)
  print(coy.bob.preyabundx$DIC)
  which(coy.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabundx$samples)
  save(coy.bob.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_rx)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, elevation, forest; psix(spp diversity); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.preydivx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preydivx$summary)
  print(coy.bob.preydivx$DIC)
  which(coy.bob.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preydivx$samples)
  save(coy.bob.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_rx)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  
