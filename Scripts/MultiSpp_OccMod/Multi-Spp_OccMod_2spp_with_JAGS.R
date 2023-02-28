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
  names(multi_spp_DH_list) <- c("wolf_bear_HD", "wolf_coy_DH", "wolf_lion_DH", 
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
  psi_covs <- matrix(NA, ncol = 16, nrow = nsites)
  psi_covs[,1] <- 1
  psi_covs[,2] <- det_covs$Setup
  psi_covs[,3] <- stations_eoe20s21s$Elev
  psi_covs[,4] <- stations_eoe20s21s$PercForest
  psi_covs[,5] <- stations_eoe20s21s$Nsmall_deer
  psi_covs[,6] <- stations_eoe20s21s$Nbig_deer
  psi_covs[,7] <- stations_eoe20s21s$Nelk
  psi_covs[,8] <- stations_eoe20s21s$Nmoose
  psi_covs[,9] <- stations_eoe20s21s$Nmd
  psi_covs[,10] <- stations_eoe20s21s$Nwtd
  psi_covs[,11] <- stations_eoe20s21s$Nlagomorph
  psi_covs[,12] <- stations_eoe20s21s$Dist2Burbs
  psi_covs[,13] <- stations_eoe20s21s$logNearestRd
  psi_covs[,14] <- stations_eoe20s21s$Nhuman
  psi_covs[,15] <- stations_eoe20s21s$Nlivestock
  psi_covs[,16] <- stations_eoe20s21s$NewLocationID
  head(psi_covs)

  ######  Second order occupancy (psi): 2-way interactions  ####
  psi_inxs_covs <- matrix(NA, ncol = 16, nrow = nsites)
  psi_inxs_covs[,1] <- 1
  psi_inxs_covs[,2] <- det_covs$Setup
  psi_inxs_covs[,3] <- stations_eoe20s21s$Elev
  psi_inxs_covs[,4] <- stations_eoe20s21s$PercForest
  psi_inxs_covs[,5] <- stations_eoe20s21s$Nsmall_deer
  psi_inxs_covs[,6] <- stations_eoe20s21s$Nbig_deer
  psi_inxs_covs[,7] <- stations_eoe20s21s$Nelk
  psi_inxs_covs[,8] <- stations_eoe20s21s$Nmoose
  psi_inxs_covs[,9] <- stations_eoe20s21s$Nmd
  psi_inxs_covs[,10] <- stations_eoe20s21s$Nwtd
  psi_inxs_covs[,11] <- stations_eoe20s21s$Nlagomorph
  psi_inxs_covs[,12] <- stations_eoe20s21s$Dist2Burbs
  psi_inxs_covs[,13] <- stations_eoe20s21s$logNearestRd
  psi_inxs_covs[,14] <- stations_eoe20s21s$Nhuman
  psi_inxs_covs[,15] <- stations_eoe20s21s$Nlivestock
  psi_inxs_covs[,16] <- stations_eoe20s21s$NewLocationID
  head(psi_inxs_covs)
  
  ######  First order detection (rho|no second spp)  ####
  rho_covs <- array(NA, dim = c(nsites, nsurveys, 5)) # last digit is number of covariates + intercept
  rho_covs[,,1] <- 1
  rho_covs[,,2] <- det_covs$CameraFacing
  rho_covs[,,3] <- det_covs$Setup
  rho_covs[,,4] <- det_covs$Height
  rho_covs[,,5] <- effort_eoe20s21s
  head(rho_covs)
  
  ######  Second order detection (rho): 2-way interactions  ####
  rho_inxs_covs <- array(NA, dim = c(nsites, nsurveys, 5))
  rho_inxs_covs[,,1] <- 1
  rho_inxs_covs[,,2] <- det_covs$CameraFacing
  rho_inxs_covs[,,3] <- det_covs$Setup
  rho_inxs_covs[,,4] <- det_covs$Height
  rho_inxs_covs[,,5] <- effort_eoe20s21s
  head(rho_inxs_covs)
  
  
  #'  -----------------------------------
  ####  Data and MCMC settings for JAGS  ####
  #'  -----------------------------------
  #####  Bundle detection history and covariate data  ####
  #'  -------------------------------------------------
  #'  Function to bundle detection and covariate data
  bundle_data <- function(obs_array, psi_covs, psi_inxs, rho_covs, rho_inxs, 
                          sites, surveys, ncats) {
    #'  List all data streams together
    bundled <- list(y = obs_array, psi_covs = psi_covs, psi_inxs_covs = psi_inxs,
                    rho_covs = rho_covs, rho_inxs_cov = rho_inxs, nsites = sites,
                    nsurveys = surveys, nfirst_order_psi = ncol(psi_covs),
                    nsecond_order_psi = ncol(psi_inxs), nfirst_order_rho = dim(rho_covs)[3],
                    nsecond_order_rho = dim(rho_inxs)[3], ncat = ncats)
    #'  Summarize to make sure it looks right
    str(bundled)
    return(bundled)
  }
  #'  Run list of 2-species detection histories through function
  bundled_pred_list <- lapply(multi_spp_DH_list, bundle_data, psi_covs = psi_covs,
                              psi_inxs = psi_inxs_covs, rho_covs = rho_covs, 
                              rho_inxs = rho_inxs_covs, sites = nsites, 
                              surveys = nsurveys, ncats = ncat)
  
  
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
  params <- c("betaSpp1", "betaSpp2", "betaSpp12", "alphaSpp1", "alphaSpp2", "alphaSpp12", 
              "mean.psiSpp1", "mean.psiSpp2", "mean.pSpp1", "mean.pSpp2", "z")
  
  #####  MCMC settings  ####
  #'  ------------------
  nc <- 3
  ni <- 50000
  nb <- 20000
  nt <- 3
  na <- 10000
  
  
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
  
  ####  Coyote-Bobcat Models  ####
  #'  ------------------------
  inits.coy.bob <- function(){list(z = zinits[[6]])}
  
  #'  Habitat model: psi = setup, elevation, forest; p = effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(habitat).R")
  start.time = Sys.time()
  coy.bob.hab <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(habitat)_p(effort).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab$summary)
  which(coy.bob.hab.summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab$samples)
  save(coy.bob.hab, file = "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(habitat)_p(effort).RData")
  
  