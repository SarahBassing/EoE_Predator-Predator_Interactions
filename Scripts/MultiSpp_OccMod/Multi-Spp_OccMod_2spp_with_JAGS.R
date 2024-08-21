  #'  ---------------------------------------
  #'  Bayesian multi-species occupancy model
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
  
  #'  Load covariate and detection history data
  load("./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_img_updated_072924.RData") # NOTE: includes Covariates_EoE_..._updated_070824.RData which includes TRI, updated sppDiversity, and 100m radius PercForest 
  #' load("./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_img.RData")
  #' 
  #' #'  Remove problematic row in 2020 detection history (GMU6_U_23 only operational for 0.5 days)
  #' drop_row <- function(det_hist, rm_rows) {
  #'   short_dh <- det_hist[[1]]
  #'   short_effort <- det_hist[[2]]
  #'   
  #'   short_dh <- short_dh[-rm_rows,]
  #'   short_effort <- short_effort[-rm_rows,]
  #'   
  #'   dh_list <- list(short_dh, short_effort)
  #' }
  #' DH_eoe20s_predators <- lapply(DH_eoe20s_predators, drop_row, rm_rows = 393)
  #' 
  #' stations_eoe20s21s <- stations_eoe20s21s %>%
  #'   filter(NewLocationID != "GMU6_U_23" | Season != "Smr20")
  #' effort_eoe20s21s <- effort_eoe20s21s[-393,]
  
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
  bear_coy_DH <- detection_array(spp1 = DH_bear, spp2 = DH_coy, name1 = "bear", name2 = "coyote")
  
  #'  List 2-species detection arrays together for faster formatting below
  DH_array_list <- list(wolf_bear_DH, wolf_coy_DH, wolf_lion_DH, lion_bear_DH, lion_bob_DH, coy_bob_DH, bear_coy_DH)
  
  
  #'  -------------------------
  ####  Format data for JAGS  ####
  #'  -------------------------
  #'  Define survey dimensions
  nsites <- dim(wolf_bear_DH)[1]
  nsurveys <- dim(wolf_bear_DH)[2]
  nspecies <- dim(wolf_bear_DH)[3]
  #'  Number of possible community states (species interactions): 00, 10, 01, 11
  ncat <- 2^nspecies
  #'  Number of sampling years
  nyears <- 2
  
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
                                "lion_bear_DH", "lion_bob_DH", "coy_bob_DH", "bear_coy_DH")
  
  #####  Format covariate data  ####
  #'  ---------------------------
  #'  Format site-level covariates for detection sub-model
  det_covs <- stations_eoe20s21s %>%
    mutate(CameraFacing = as.factor(CameraFacing),
           Setup = ifelse(Setup == "ungulate", 0, 1),
           Height = as.numeric(Height))
  table(det_covs[,"CameraFacing"])
  table(det_covs[,"Setup"])
  
  stations_eoe20s21s <- mutate(stations_eoe20s21s, Year = ifelse(Season == "Smr20", 0, 1))
  
  ######  First order occupancy (psi|no second spp)  ####
  psi_cov <- matrix(NA, ncol = 17, nrow = nsites)
  psi_cov[,1] <- 1
  psi_cov[,2] <- det_covs$Setup
  psi_cov[,3] <- stations_eoe20s21s$Year 
  psi_cov[,4] <- stations_eoe20s21s$PercForest
  psi_cov[,5] <- stations_eoe20s21s$Elev
  psi_cov[,6] <- stations_eoe20s21s$SppDiversity
  psi_cov[,7] <- stations_eoe20s21s$Nelk
  psi_cov[,8] <- stations_eoe20s21s$Nmoose
  psi_cov[,9] <- stations_eoe20s21s$Nmd
  psi_cov[,10] <- stations_eoe20s21s$Nwtd
  psi_cov[,11] <- stations_eoe20s21s$Nlagomorph
  psi_cov[,12] <- stations_eoe20s21s$Dist2Burbs
  psi_cov[,13] <- stations_eoe20s21s$logNearestRd
  psi_cov[,14] <- stations_eoe20s21s$Nhuman
  psi_cov[,15] <- stations_eoe20s21s$Nlivestock
  psi_cov[,16] <- stations_eoe20s21s$Footprint #NewLocationID
  psi_cov[,17] <- stations_eoe20s21s$TRI
  head(psi_cov)

  ######  Second order occupancy (psi): 2-way interactions  ####
  psi_inxs_cov <- matrix(NA, ncol = 17, nrow = nsites)
  psi_inxs_cov[,1] <- 1
  psi_inxs_cov[,2] <- det_covs$Setup
  psi_inxs_cov[,3] <- stations_eoe20s21s$Year 
  psi_inxs_cov[,4] <- stations_eoe20s21s$PercForest
  psi_inxs_cov[,5] <- stations_eoe20s21s$Elev
  psi_inxs_cov[,6] <- stations_eoe20s21s$SppDiversity
  psi_inxs_cov[,7] <- stations_eoe20s21s$Nelk
  psi_inxs_cov[,8] <- stations_eoe20s21s$Nmoose
  psi_inxs_cov[,9] <- stations_eoe20s21s$Nmd
  psi_inxs_cov[,10] <- stations_eoe20s21s$Nwtd
  psi_inxs_cov[,11] <- stations_eoe20s21s$Nlagomorph
  psi_inxs_cov[,12] <- stations_eoe20s21s$Dist2Burbs
  psi_inxs_cov[,13] <- stations_eoe20s21s$logNearestRd
  psi_inxs_cov[,14] <- stations_eoe20s21s$Nhuman
  psi_inxs_cov[,15] <- stations_eoe20s21s$Nlivestock
  psi_inxs_cov[,16] <- stations_eoe20s21s$Footprint #NewLocationID
  psi_inxs_cov[,17] <- stations_eoe20s21s$TRI
  head(psi_inxs_cov)
  
  ######  First order detection (rho|no second spp)  ####
  rho_cov <- array(NA, dim = c(nsites, nsurveys, 6)) # last digit is number of covariates + intercept
  rho_cov[,,1] <- 1
  rho_cov[,,2] <- det_covs$CameraFacing
  rho_cov[,,3] <- det_covs$Setup
  rho_cov[,,4] <- det_covs$Height
  rho_cov[,,5] <- effort_eoe20s21s
  rho_cov[,,6] <- stations_eoe20s21s$Year
  head(rho_cov)
  
  ######  Second order detection (rho): 2-way interactions  ####
  rho_inxs_cov <- array(NA, dim = c(nsites, nsurveys, 6))
  rho_inxs_cov[,,1] <- 1
  rho_inxs_cov[,,2] <- det_covs$CameraFacing
  rho_inxs_cov[,,3] <- det_covs$Setup
  rho_inxs_cov[,,4] <- det_covs$Height
  rho_inxs_cov[,,5] <- effort_eoe20s21s
  rho_inxs_cov[,,6] <- stations_eoe20s21s$Year
  head(rho_inxs_cov)
  
  
  #'  -----------------------------------
  ####  Data and MCMC settings for JAGS  ####
  #'  -----------------------------------
  #####  Bundle detection history and covariate data  ####
  #'  -------------------------------------------------
  #'  Function to bundle detection and covariate data
  bundle_data <- function(obs_array, psi_covs, psi_inxs, rho_covs, rho_inxs, 
                          sites, surveys, psi_1order, psi_2order, rho_1order, 
                          rho_2order, ncats, nspecies, nyears, uniquesites) {
    #'  list all pieces of data together
    bundled <- list(y = obs_array, psi_cov = psi_covs, psi_inxs_cov = psi_inxs,
                    rho_cov = rho_covs, rho_inxs_cov = rho_inxs, nsites = sites,
                    nsurveys = surveys, nfirst_order_psi = ncol(psi_1order), 
                    nsecond_order_psi = ncol(psi_2order), 
                    nfirst_order_rho = dim(rho_1order)[3], 
                    nsecond_order_rho = dim(rho_2order)[3], ncat = ncats, 
                    nspec = nspecies, nyear = nyears, Indx = c(1, 2, 3, 4),
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
                               nspecies = nspecies, nyears = nyears, uniquesites = unique(psi_cov[,16])) 
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
                     "lion_bear_zcat", "lion_bob_zcat", "coy_bob_zcat", "bear_coy_zcat")
  
  #####  Parameters monitored  ####
  #'  -------------------------
  params <- c("betaSpp1", "betaSpp2", "alphaSpp1", "alphaSpp2", "betaSpp12", 
              "alphaSpp12", "alphaSpp21", "mean.psiSpp1", "mean.psiSpp2", 
              "mean.pSpp1", "mean.pSpp2", "z", 
              "y.hat", "y.hat.index", "y.hat.maxindex", "x2", "x2.sim", "chi2.obs", "chi2.sim")
              # "y2", "y_A", "y_B", "yrep2", "yrep_A", "yrep_B",
              # "pA", "pB", "pAB", "pSpp1", "pSpp2",
              # "detfreq_A", "detfreq_B", "detfreqrep_A", "detfreqrep_B",
              # "tmp_A", "tmp_B", "E_A", "E_B",
              # "x2_A", "x2_B", "x2rep_A", "x2rep_B",
              # "chi2.obs_A", "chi2.obs_B", "chi2.sim_A", "chi2.sim_B",
              # "ft_A", "ft_B", "ftrep_A", "ftrep_B", "ft.obs_A", "ft.obs_B", "ft.sim_A", "ft.sim_B",
              # # "d_A", "d_B", "d2_A", "d2_B", "dnew_A", "dnew_B", "dnew2_A", "dnew2_B",
              # # "dsum_A", "dsum_B", "dnewsum_A", "dnewsum_B",
              # "chi2ratio_A", "chi2ratio_B", "ftratio_A", "ftratio_B")
              # # "y.hat", "y.sim", "y.sim.hat", "chi2.obs", "chi2.sim", ) #"z.sim",
  
  #####  MCMC settings  ####
  #'  ------------------
  nc <- 3
  # ni <- 50000 
  nb <- 15000 
  nt <- 10
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
  inits.wolf.bear <- function(){list(z = zinits[[1]])}    
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.bear.null <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.null$summary)
  print(wolf.bear.null$DIC)
  which(wolf.bear.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.null$samples)
  save(wolf.bear.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_", Sys.Date(), ".RData"))
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_GoF.R")
  start.time = Sys.time()
  wolf.bear.null <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_GoF.txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.null$summary)
  print(wolf.bear.null$DIC)
  which(wolf.bear.null$summary[,"Rhat"] > 1.1)
  (wolf.bear.null_X2pB.wolf <- mean(wolf.bear.null$sims.list$chi2.sim_A > wolf.bear.null$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  (wolf.bear.null_X2pB.bear <- mean(wolf.bear.null$sims.list$chi2.sim_B > wolf.bear.null$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  (wolf.bear.null_FTpB.wolf <- mean(wolf.bear.null$sims.list$ft.sim_A > wolf.bear.null$sims.list$ft.obs_A)) # Bayesian p-value GOF
  (wolf.bear.null_FTpB.bear <- mean(wolf.bear.null$sims.list$ft.sim_B > wolf.bear.null$sims.list$ft.obs_B)) # Bayesian p-value GOF
  mean(wolf.bear.null$sims.list$chi2ratio_A); mean(wolf.bear.null$sims.list$chi2ratio_B); mean(wolf.bear.null$sims.list$ftratio_A); mean(wolf.bear.null$sims.list$ftratio_B)
  # mcmcplot(wolf.bear.null$samples)
  save(wolf.bear.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_GoF_", Sys.Date(), ".RData"))
  
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.hab <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.hab$summary)
  print(wolf.bear.hab$DIC)
  which(wolf.bear.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.hab$samples)
  save(wolf.bear.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_GoF.R")
  start.time = Sys.time()
  wolf.bear.hab <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_GoF.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.hab$summary)
  print(wolf.bear.hab$DIC)
  which(wolf.bear.hab$summary[,"Rhat"] > 1.1)
  (wolf.bear.hab_X2pB.wolf <- mean(wolf.bear.hab$sims.list$chi2.sim_A > wolf.bear.hab$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  (wolf.bear.hab_X2pB.bear <- mean(wolf.bear.hab$sims.list$chi2.sim_B > wolf.bear.hab$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  (wolf.bear.hab_FTpB.wolf <- mean(wolf.bear.hab$sims.list$ft.sim_A > wolf.bear.hab$sims.list$ft.obs_A)) # Bayesian p-value GOF
  (wolf.bear.hab_FTpB.bear <- mean(wolf.bear.hab$sims.list$ft.sim_B > wolf.bear.hab$sims.list$ft.obs_B)) # Bayesian p-value GOF
  mean(wolf.bear.hab$sims.list$chi2ratio_A); mean(wolf.bear.hab$sims.list$chi2ratio_B); mean(wolf.bear.hab$sims.list$ftratio_A); mean(wolf.bear.hab$sims.list$ftratio_B)
  # mcmcplot(wolf.bear.hab$samples)
  save(wolf.bear.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_GoF_", Sys.Date(), ".RData"))
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, elk, moose, wtd; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.preyabund <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preyabund$summary)
  print(wolf.bear.preyabund$DIC)
  which(wolf.bear.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preyabund$samples)
  save(wolf.bear.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preyabund_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #' #####  Prey diversity no inxs model  #### 
  #' #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).R")
  #' start.time = Sys.time()
  #' wolf.bear.preydiv <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
  #'                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).txt",
  #'                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(wolf.bear.preydiv$summary)
  #' print(wolf.bear.preydiv$DIC)
  #' which(wolf.bear.preydiv$summary[,"Rhat"] > 1.1)
  #' mcmcplot(wolf.bear.preydiv$samples)
  #' save(wolf.bear.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.habx <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.habx$summary)
  print(wolf.bear.habx$DIC)
  which(wolf.bear.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.habx$samples)
  save(wolf.bear.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ habitat inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(forest, elevation, tri); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.bear.hab2x <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.hab2x$summary)
  print(wolf.bear.hab2x$DIC)
  which(wolf.bear.hab2x$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.hab2x$samples)
  save(wolf.bear.hab2x, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, moose, wtd); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.preyabundx <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.preyabundx$summary)
  print(wolf.bear.preyabundx$DIC)
  which(wolf.bear.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.preyabundx$samples)
  save(wolf.bear.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Full habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd, lagomorphs, forest, elevation, tri); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.bear.hab.preyabundx <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                                 "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_wolfbearlion.txt",
                                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.hab.preyabundx$summary)
  print(wolf.bear.hab.preyabundx$DIC)
  which(wolf.bear.hab.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.hab.preyabundx$samples)
  save(wolf.bear.hab.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(preyabund_habitat)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  
  #' #####  Habitat w/ prey diversity inx model  #### 
  #' #'  psi = setup, year, forest, elevation, tri; psix(spp diversity); p = setup, effort  
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  #' start.time = Sys.time()
  #' wolf.bear.preydivx <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
  #'                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
  #'                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(wolf.bear.preydivx$summary)
  #' print(wolf.bear.preydivx$DIC)
  #' which(wolf.bear.preydivx$summary[,"Rhat"] > 1.1)
  #' mcmcplot(wolf.bear.preydivx$samples)
  #' save(wolf.bear.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  #' 
  #' #####  Global model  #### 
  #' #'  psi = setup, year, forest, elevation, tri; psix(elk, moose, wtd, spp diversity); p = setup, effort  
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.R")
  #' start.time = Sys.time()
  #' wolf.bear.global <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
  #'                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.txt",
  #'                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(wolf.bear.global$summary)
  #' print(wolf.bear.global$DIC)
  #' which(wolf.bear.global$summary[,"Rhat"] > 1.1)
  #' mcmcplot(wolf.bear.global$samples)
  #' save(wolf.bear.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #' #####  Top model w/ intx on detection model v1  #### 
  #' #'  Parameterization tests whether detection of one predator affects detection of the other
  #' #'  Top model:  Prey diversity, no intx on psi
  #' #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort; px(.)
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort)_px(.).R")
  #' start.time = Sys.time()
  #' wolf.bear.preydiv.px <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
  #'                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort)_px(.).txt",
  #'                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(wolf.bear.preydiv.px$summary)
  #' print(wolf.bear.preydiv.px$DIC)
  #' which(wolf.bear.preydiv.px$summary[,"Rhat"] > 1.1)
  #' mcmcplot(wolf.bear.preydiv.px$samples)
  #' save(wolf.bear.preydiv.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_px(.)_", Sys.Date(), ".RData"))
  #' 
  #' #####  Top model w/ intx on detection model v2  #### 
  #' #'  Parameterization tests whether presence of one predator affects detection of the other
  #' #'  Top model:  Prey diversity, no intx on psi
  #' #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort; px(psi) 
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort)_px(psi).R")
  #' start.time = Sys.time()
  #' wolf.bear.preydiv.px2 <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
  #'                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort)_px(psi).txt",
  #'                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(wolf.bear.preydiv.px2$summary)
  #' print(wolf.bear.preydiv.px2$DIC)
  #' which(wolf.bear.preydiv.px2$summary[,"Rhat"] > 1.1)
  #' mcmcplot(wolf.bear.preydiv.px2$samples)
  #' save(wolf.bear.preydiv.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_px(psi)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v1  #### 
  #'  Top model: null
  #'  psi = year; p(.); px(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  wolf.bear.null.px <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(.).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.null.px$summary)
  print(wolf.bear.null.px$DIC)
  which(wolf.bear.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.null.px$samples)
  save(wolf.bear.null.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_px(.)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model: null
  #'  psi = year; p(.); px(psi)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(psi).R")
  start.time = Sys.time()
  wolf.bear.null.px2 <- jags(bundled_pred_list[[1]], inits = inits.wolf.bear, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(psi).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.bear.null.px2$summary)
  print(wolf.bear.null.px2$DIC)
  which(wolf.bear.null.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.bear.null.px2$samples)
  save(wolf.bear.null.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_px(psi)_", Sys.Date(), ".RData"))
  
  
  #'  ----------------------
  ####  Wolf-Coyote Models  ####
  #'  ----------------------
  inits.wolf.coy <- function(){list(z = zinits[[2]])} 
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.coy.null <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.null$summary)
  print(wolf.coy.null$DIC)
  which(wolf.coy.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.null$samples)
  save(wolf.coy.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(yr)_p(.)_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, tri; p = setup, effort  
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
  save(wolf.coy.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_altGoF.R")
  start.time = Sys.time()
  wolf.coy.hab <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_altGoF.txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab$summary)
  print(wolf.coy.hab$DIC)
  which(wolf.coy.hab$summary[,"Rhat"] > 1.1)
  (wolf.coy.hab_X2 <- mean(wolf.coy.hab$sims.list$chi2.sim > wolf.coy.hab$sims.list$chi2.obs)) # Bayesian p-value GOF
  # (wolf.coy.hab_X2pB.wolf <- mean(wolf.coy.hab$sims.list$chi2.sim_A > wolf.coy.hab$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  # (wolf.coy.hab_X2pB.coy <- mean(wolf.coy.hab$sims.list$chi2.sim_B > wolf.coy.hab$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  # (wolf.coy.hab_FTpB.wolf <- mean(wolf.coy.hab$sims.list$ft.sim_A > wolf.coy.hab$sims.list$ft.obs_A)) # Bayesian p-value GOF
  # (wolf.coy.hab_FTpB.coy <- mean(wolf.coy.hab$sims.list$ft.sim_B > wolf.coy.hab$sims.list$ft.obs_B)) # Bayesian p-value GOF
  # mean(wolf.coy.hab$sims.list$chi2ratio_A); mean(wolf.coy.hab$sims.list$chi2ratio_B); mean(wolf.coy.hab$sims.list$ftratio_A); mean(wolf.coy.hab$sims.list$ftratio_B)
  # mcmcplot(wolf.coy.hab$samples)
  save(wolf.coy.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_GoF_", Sys.Date(), ".RData"))
  
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, elk, moose, wtd; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabund <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_wolfcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabund$summary)
  print(wolf.coy.preyabund$DIC)
  which(wolf.coy.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabund$samples)
  save(wolf.coy.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preyabund_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.preydiv <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preydiv$summary)
  print(wolf.coy.preydiv$DIC)
  which(wolf.coy.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preydiv$samples)
  save(wolf.coy.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_preydiversity_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.habx <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.habx$summary)
  print(wolf.coy.habx$DIC)
  which(wolf.coy.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.habx$samples)
  save(wolf.coy.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, moose, wtd, lagomorph); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.preyabundx <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfcoy.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preyabundx$summary)
  print(wolf.coy.preyabundx$DIC)
  which(wolf.coy.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preyabundx$samples)
  save(wolf.coy.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(spp diversity); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.coy.preydivx <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.preydivx$summary)
  print(wolf.coy.preydivx$DIC)
  which(wolf.coy.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.preydivx$samples)
  save(wolf.coy.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))

  #####  Global model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, moose, wtd, lagomorph, spp diversity); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_wolfcoy.R")
  start.time = Sys.time()
  wolf.coy.global <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_wolfcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.global$summary)
  print(wolf.coy.global$DIC)
  which(wolf.coy.global$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.global$samples)
  save(wolf.coy.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v1  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort; px(.) 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_px(.).R")
  start.time = Sys.time()
  wolf.coy.hab.px <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_px(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab.px$summary)
  print(wolf.coy.hab.px$DIC)
  which(wolf.coy.hab.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab.px$samples)
  save(wolf.coy.hab.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_px(.)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort; px(psi) 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_px(psi).R")
  start.time = Sys.time()
  wolf.coy.hab.px2 <- jags(bundled_pred_list[[2]], inits = inits.wolf.coy, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_px(psi).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.coy.hab.px2$summary)
  print(wolf.coy.hab.px2$DIC)
  which(wolf.coy.hab.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.coy.hab.px2$samples)
  save(wolf.coy.hab.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_px(psi)_", Sys.Date(), ".RData"))
  
  
  
  #'  ---------------------
  ####  Wolf-Lion Models  ####
  #'  ---------------------
  inits.wolf.lion <- function(){list(z = zinits[[3]], mean.psiSpp1 = runif(1, 0.01, 0.15), 
                                     mean.psiSpp2 = runif(1, 0.2, 0.3),
                                     mean.pSpp1 = runif(1, 0.02, 0.12), mean.pSpp2 =  runif(1, 0.01, 0.06))}  
  ni <- 100000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.).R")
  start.time = Sys.time()
  wolf.lion.null <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null$summary)
  print(wolf.lion.null$DIC)
  which(wolf.lion.null$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null$samples)
  save(wolf.lion.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_", Sys.Date(), ".RData"))

  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_GoF.R")
  start.time = Sys.time()
  wolf.lion.null <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_GoF.txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null$summary)
  print(wolf.lion.null$DIC)
  which(wolf.lion.null$summary[,"Rhat"] > 1.1)
  (wolf.lion.null_X2pB.wolf <- mean(wolf.lion.null$sims.list$chi2.sim_A > wolf.lion.null$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  (wolf.lion.null_X2pB.lion <- mean(wolf.lion.null$sims.list$chi2.sim_B > wolf.lion.null$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  (wolf.lion.null_FTpB.wolf <- mean(wolf.lion.null$sims.list$ft.sim_A > wolf.lion.null$sims.list$ft.obs_A)) # Bayesian p-value GOF
  (wolf.lion.null_FTpB.lion <- mean(wolf.lion.null$sims.list$ft.sim_B > wolf.lion.null$sims.list$ft.obs_B)) # Bayesian p-value GOF
  mean(wolf.lion.null$sims.list$chi2ratio_A); mean(wolf.lion.null$sims.list$chi2ratio_B); mean(wolf.lion.null$sims.list$ftratio_A); mean(wolf.lion.null$sims.list$ftratio_B)
  # mcmcplot(wolf.lion.null$samples)
  save(wolf.lion.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_GoF_", Sys.Date(), ".RData"))
  
    
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.hab <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.hab$summary)
  print(wolf.lion.hab$DIC)
  which(wolf.lion.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.hab$samples)
  save(wolf.lion.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, elk, moose, wtd; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabund <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_wolfbearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabund$summary)
  print(wolf.lion.preyabund$DIC)
  which(wolf.lion.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabund$samples)
  save(wolf.lion.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preyabund_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.preydiv <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preydiv$summary)
  print(wolf.lion.preydiv$DIC)
  which(wolf.lion.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preydiv$samples)
  save(wolf.lion.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_preydiversity_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.habx <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.habx$summary)
  print(wolf.lion.habx$DIC)
  which(wolf.lion.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.habx$samples)
  save(wolf.lion.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, moose, wtd); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.preyabundx <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_wolfbearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preyabundx$summary)
  print(wolf.lion.preyabundx$DIC)
  which(wolf.lion.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preyabundx$samples)
  save(wolf.lion.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(spp diversity); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  wolf.lion.preydivx <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.preydivx$summary)
  print(wolf.lion.preydivx$DIC)
  which(wolf.lion.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.preydivx$samples)
  save(wolf.lion.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
 
  #####  Global model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, moose, wtd, spp diversity); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.R")
  start.time = Sys.time()
  wolf.lion.global <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_wolfbearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.global$summary)
  print(wolf.lion.global$DIC)
  which(wolf.lion.global$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.global$samples)
  save(wolf.lion.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v1  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  wolf.lion.null.px <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(.).txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null.px$summary)
  print(wolf.lion.null.px$DIC)
  which(wolf.lion.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null.px$samples)
  save(wolf.lion.null.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_px(.)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = year; p(.); px(psi)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(psi).R")
  start.time = Sys.time()
  wolf.lion.null.px2 <- jags(bundled_pred_list[[3]], inits = inits.wolf.lion, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(psi).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(wolf.lion.null.px2$summary)
  print(wolf.lion.null.px2$DIC)
  which(wolf.lion.null.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(wolf.lion.null.px2$samples)
  save(wolf.lion.null.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_px(psi)_", Sys.Date(), ".RData"))
  
  
  
  #'  --------------------
  ####  Lion-Bear Models  ####
  #'  --------------------
  inits.lion.bear <- function(){list(z = zinits[[4]], mean.psiSpp1 = runif(1, 0.1, 0.35), 
                                     mean.psiSpp2 = runif(1, 0.6, 0.7), mean.pSpp1 = runif(1, 0.01, 0.06),
                                     mean.pSpp2 = runif(1, 0.15, 0.2))}
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.).R")
  start.time = Sys.time()
  lion.bear.null <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null$summary)
  print(lion.bear.null$DIC)
  which(lion.bear.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null$samples)
  save(lion.bear.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_", Sys.Date(), ".RData"))
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_GoF.R")
  start.time = Sys.time()
  lion.bear.null <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_GoF.txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null$summary)
  print(lion.bear.null$DIC)
  which(lion.bear.null$summary[,"Rhat"] > 1.1)
  (lion.bear.null_X2pB.lion <- mean(lion.bear.null$sims.list$chi2.sim_A > lion.bear.null$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  (lion.bear.null_X2pB.bear <- mean(lion.bear.null$sims.list$chi2.sim_B > lion.bear.null$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  (lion.bear.null_FTpB.lion <- mean(lion.bear.null$sims.list$ft.sim_A > lion.bear.null$sims.list$ft.obs_A)) # Bayesian p-value GOF
  (lion.bear.null_FTpB.bear <- mean(lion.bear.null$sims.list$ft.sim_B > lion.bear.null$sims.list$ft.obs_B)) # Bayesian p-value GOF
  mean(lion.bear.null$sims.list$chi2ratio_A); mean(lion.bear.null$sims.list$chi2ratio_B); mean(lion.bear.null$sims.list$ftratio_A); mean(lion.bear.null$sims.list$ftratio_B)
  # mcmcplot(lion.bear.null$samples)
  save(lion.bear.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_GoF_", Sys.Date(), ".RData"))
  
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort
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
  #'  psi = setup, year, forest, elevation, tri, elk, wtd; p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabund <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_bearlion.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabund$summary)
  print(lion.bear.preyabund$DIC)
  which(lion.bear.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabund$samples)
  save(lion.bear.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preyabund_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.preydiv <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preydiv$summary)
  print(lion.bear.preydiv$DIC)
  which(lion.bear.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preydiv$samples)
  save(lion.bear.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_preydiversity_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.habx <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.habx$summary)
  print(lion.bear.habx$DIC)
  which(lion.bear.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.habx$samples)
  save(lion.bear.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.preyabundx <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearlion.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preyabundx$summary)
  print(lion.bear.preyabundx$DIC)
  which(lion.bear.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preyabundx$samples)
  save(lion.bear.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(spp diversity); p = setup, effort  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bear.preydivx <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.preydivx$summary)
  print(lion.bear.preydivx$DIC)
  which(lion.bear.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.preydivx$samples)
  save(lion.bear.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Global model  ####
  #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd, spp diversity); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_bearlion.R")
  start.time = Sys.time()
  lion.bear.global <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_bearlion.txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.global$summary)
  print(lion.bear.global$DIC)
  which(lion.bear.global$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.global$samples)
  save(lion.bear.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_bearlion_GoF.R")
  start.time = Sys.time()
  lion.bear.global <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params, 
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_bearlion_GoF.txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.global$summary)
  (lion.bear.global.pval <- mean(lion.bear.global$sims.list$chi2.sim > lion.bear.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  print(lion.bear.global$DIC)
  which(lion.bear.global$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.global$samples)
  save(lion.bear.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  
  
  #####  Top model w/ intx on detection model v1  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  lion.bear.null.px <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(.).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null.px$summary)
  print(lion.bear.null.px$DIC)
  which(lion.bear.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null.px$samples)
  save(lion.bear.null.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_px(.)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = year; p(.); px(psi)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(psi).R")
  start.time = Sys.time()
  lion.bear.null.px2 <- jags(bundled_pred_list[[4]], inits = inits.lion.bear, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(psi).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bear.null.px2$summary)
  print(lion.bear.null.px2$DIC)
  which(lion.bear.null.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bear.null.px2$samples)
  save(lion.bear.null.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_px(psi)_", Sys.Date(), ".RData"))
  
  
  
  #'  ----------------------
  ####  Lion-Bobcat Models  ####
  #'  ----------------------
  inits.lion.bob <- function(){list(z = zinits[[5]], mean.psiSpp1 = runif(1, 0.2, 0.35),
                                    mean.psiSpp2 = runif(1, 0.1, 0.3), mean.pSpp1 = runif(1, 0.01, 0.1),
                                    mean.pSpp2 = runif(1, 0.02, 0.1))}
  ni <- 75000

  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.).R")
  start.time = Sys.time()
  lion.bob.null <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null$summary)
  print(lion.bob.null$DIC)
  which(lion.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null$samples)
  save(lion.bob.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_", Sys.Date(), ".RData"))

  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_GoF.R")
  start.time = Sys.time()
  lion.bob.null <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_GoF.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null$summary)
  print(lion.bob.null$DIC)
  which(lion.bob.null$summary[,"Rhat"] > 1.1)
  (lion.bob.null_X2pB.lion <- mean(lion.bob.null$sims.list$chi2.sim_A > lion.bob.null$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  (lion.bob.null_X2pB.bob <- mean(lion.bob.null$sims.list$chi2.sim_B > lion.bob.null$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  (lion.bob.null_FTpB.lion <- mean(lion.bob.null$sims.list$ft.sim_A > lion.bob.null$sims.list$ft.obs_A)) # Bayesian p-value GOF
  (lion.bob.null_FTpB.bob <- mean(lion.bob.null$sims.list$ft.sim_B > lion.bob.null$sims.list$ft.obs_B)) # Bayesian p-value GOF
  mean(lion.bob.null$sims.list$chi2ratio_A); mean(lion.bob.null$sims.list$chi2ratio_B); mean(lion.bob.null$sims.list$ftratio_A); mean(lion.bob.null$sims.list$ftratio_B)
  # mcmcplot(lion.bob.null$samples)
  save(lion.bob.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_GoF_", Sys.Date(), ".RData"))
  
    
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.hab <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.hab$summary)
  print(lion.bob.hab$DIC)
  which(lion.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.hab$samples)
  save(lion.bob.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, elk, wtd, lagomorphs; p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.preyabund <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_lionbob.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabund$summary)
  print(lion.bob.preyabund$DIC)
  which(lion.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabund$samples)
  save(lion.bob.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preyabund_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey diversity no inxs model  #### 
  #'  psi = setup, year, elevation, forest, spp diversity; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.preydiv <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preydiv$summary)
  print(lion.bob.preydiv$DIC)
  which(lion.bob.preydiv$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preydiv$samples)
  save(lion.bob.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_preydiversity_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.habx <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.habx$summary)
  print(lion.bob.habx$DIC)
  which(lion.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.habx$samples)
  save(lion.bob.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd, lagomorphs); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.preyabundx <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_lionbob.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preyabundx$summary)
  print(lion.bob.preyabundx$DIC)
  which(lion.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preyabundx$samples)
  save(lion.bob.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey diversity inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(spp diversity); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  start.time = Sys.time()
  lion.bob.preydivx <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.preydivx$summary)
  print(lion.bob.preydivx$DIC)
  which(lion.bob.preydivx$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.preydivx$samples)
  save(lion.bob.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Global model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd, lagomorphs, spp diversity); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_lionbob.R")
  start.time = Sys.time()
  lion.bob.global <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_lionbob.txt",
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.global$summary)
  print(lion.bob.global$DIC)
  which(lion.bob.global$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.global$samples)
  save(lion.bob.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v1  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(.).R")
  start.time = Sys.time()
  lion.bob.null.px <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(.).txt",
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null.px$summary)
  print(lion.bob.null.px$DIC)
  which(lion.bob.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null.px$samples)
  save(lion.bob.null.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_px(.)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = year; p(.); px(psi)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.)_px(psi).R")
  start.time = Sys.time()
  lion.bob.null.px2 <- jags(bundled_pred_list[[5]], inits = inits.lion.bob, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.)_px(psi).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(lion.bob.null.px2$summary)
  print(lion.bob.null.px2$DIC)
  which(lion.bob.null.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(lion.bob.null.px2$samples)
  save(lion.bob.null.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_px(psi)_", Sys.Date(), ".RData"))
  
  
  
  #'  ------------------------
  ####  Coyote-Bobcat Models  ####
  #'  ------------------------
  inits.coy.bob <- function(){list(z = zinits[[6]], mean.psiSpp1 = runif(1, 0.3, 0.4), 
                                   mean.psiSpp2 = runif(1, 0.03, 0.18),
                                   mean.pSpp1 = runif(1, 0.15, 0.2), mean.pSpp2 = runif(1, 0.03, 0.1))}
  ni <- 100000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.).R")
  start.time = Sys.time()
  coy.bob.null <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.null$summary)
  print(coy.bob.null$DIC)
  which(coy.bob.null$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.null$samples)
  save(coy.bob.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(yr)_p(.)_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort)_GoF.R")
  start.time = Sys.time()
  coy.bob.hab <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                      "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).txt",
                      n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab$summary)
  print(coy.bob.hab$DIC)
  which(coy.bob.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab$samples)
  save(coy.bob.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, wtd, lagomorphs; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabund <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_coybob.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabund$summary)
  print(coy.bob.preyabund$DIC)
  which(coy.bob.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabund$samples)
  save(coy.bob.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preyabund_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #' #####  Prey diversity no inxs model  #### 
  #' #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).R")
  #' start.time = Sys.time()
  #' coy.bob.preydiv <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
  #'                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).txt",
  #'                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(coy.bob.preydiv$summary)
  #' print(coy.bob.preydiv$DIC)
  #' which(coy.bob.preydiv$summary[,"Rhat"] > 1.1)
  #' mcmcplot(coy.bob.preydiv$samples)
  #' save(coy.bob.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_preydiversity_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.habx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx$summary)
  print(coy.bob.habx$DIC)
  which(coy.bob.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.habx$samples)
  save(coy.bob.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_altGoF.R")
  start.time = Sys.time()
  coy.bob.habx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_altGoF.txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx$summary)
  print(coy.bob.habx$DIC)
  which(coy.bob.habx$summary[,"Rhat"] > 1.1)
  (coy.bob.habx_X2 <- mean(coy.bob.habx$sims.list$chi2.sim > coy.bob.habx$sims.list$chi2.obs)) # Bayesian p-value GOF
  # (coy.bob.habx_X2pB.coy <- mean(coy.bob.habx$sims.list$chi2.sim_A > coy.bob.habx$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  # (coy.bob.habx_X2pB.bob <- mean(coy.bob.habx$sims.list$chi2.sim_B > coy.bob.habx$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  # (coy.bob.habx_FTpB.coy <- mean(coy.bob.habx$sims.list$ft.sim_A > coy.bob.habx$sims.list$ft.obs_A)) # Bayesian p-value GOF
  # (coy.bob.habx_FTpB.bob <- mean(coy.bob.habx$sims.list$ft.sim_B > coy.bob.habx$sims.list$ft.obs_B)) # Bayesian p-value GOF
  # # (coy.bob.habx_PpB.coy <- mean(coy.bob.habx$sims.list$dnewsum_A > coy.bob.habx$sims.list$dsum_A)) # Bayesian p-value GOF
  # # (coy.bob.habx_PpB.bob <- mean(coy.bob.habx$sims.list$dnewsum_B > coy.bob.habx$sims.list$dsum_B)) # Bayesian p-value GOF
  # mean(coy.bob.habx$sims.list$chi2ratio_A); mean(coy.bob.habx$sims.list$chi2ratio_B); mean(coy.bob.habx$sims.list$ftratio_A); mean(coy.bob.habx$sims.list$ftratio_B)
  # mcmcplot(coy.bob.habx$samples)
  save(coy.bob.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF_", Sys.Date(), ".RData"))

  
  #####  Habitat w/ habitat inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(forest, elevation, tri); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort).R")
  start.time = Sys.time()
  coy.bob.hab2x <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort).txt",
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab2x$summary)
  print(coy.bob.hab2x$DIC)
  which(coy.bob.hab2x$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab2x$samples)
  save(coy.bob.hab2x, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(wtd, lagomorphs); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.preyabundx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                               "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_coybob.txt",
                               n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.preyabundx$summary)
  print(coy.bob.preyabundx$DIC)
  which(coy.bob.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.preyabundx$samples)
  save(coy.bob.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Full habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(wtd, lagomorphs, forest, elevation, tri); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_coybob.R")
  start.time = Sys.time()
  coy.bob.hab.preyabundx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                                 "./Outputs/MultiSpp_OccMod_Outputs/JAGS_code_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_coybob.txt",
                                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.hab.preyabundx$summary)
  print(coy.bob.hab.preyabundx$DIC)
  which(coy.bob.hab.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.hab.preyabundx$samples)
  save(coy.bob.hab.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  
  #' #####  Habitat w/ prey diversity inx model  #### 
  #' #'  psi = setup, year, forest, elevation, tri; psix(spp diversity); p = setup, effort
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  #' start.time = Sys.time()
  #' coy.bob.preydivx <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
  #'                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
  #'                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(coy.bob.preydivx$summary)
  #' print(coy.bob.preydivx$DIC)
  #' which(coy.bob.preydivx$summary[,"Rhat"] > 1.1)
  #' mcmcplot(coy.bob.preydivx$samples)
  #' save(coy.bob.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  #' 
  #' #####  Global model  #### 
  #' #'  psi = setup, year, forest, elevation, tri; psix(wtd, lagomorphs, spp diversity); p = setup, effort
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_coybob.R")
  #' start.time = Sys.time()
  #' coy.bob.global <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
  #'                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_coybob.txt",
  #'                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(coy.bob.global$summary)
  #' print(coy.bob.global$DIC)
  #' which(coy.bob.global$summary[,"Rhat"] > 1.1)
  #' mcmcplot(coy.bob.global$samples)
  #' save(coy.bob.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  #' 
  #' 
  #' 
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_coybob_GoF.R")
  #' start.time = Sys.time()
  #' coy.bob.global <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
  #'                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_coybob_GoF.txt",
  #'                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(coy.bob.global$summary)
  #' (coy.bob.global.pval <- mean(coy.bob.global$sims.list$chi2.sim > coy.bob.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  #' print(coy.bob.global$DIC)
  #' which(coy.bob.global$summary[,"Rhat"] > 1.1)
  #' mcmcplot(coy.bob.global$samples)
  #' save(coy.bob.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  
  #' #####  Top model w/ intx on detection model v1  #### 
  #' #'  Parameterization tests whether detection of one predator affects detection of the other
  #' #'  Top model:  Global
  #' #'  psi = setup, year, forest, elevation, tri; psix(wtd, lagomorphs, spp diversity); p = setup, effort; px(.)
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_px(.)_coybob.R")
  #' start.time = Sys.time()
  #' coy.bob.global.px <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
  #'                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_px(.)_coybob.txt",
  #'                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(coy.bob.global.px$summary)
  #' print(coy.bob.global.px$DIC)
  #' which(coy.bob.global.px$summary[,"Rhat"] > 1.1)
  #' mcmcplot(coy.bob.global.px$samples)
  #' save(coy.bob.global.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_px(.)_coybob_", Sys.Date(), ".RData"))
  #' 
  #' #####  Top model w/ intx on detection model v2  #### 
  #' #'  Parameterization tests whether presence of one predator affects detection of the other
  #' #'  Top model:  Global
  #' #'  psi = setup, year, forest, elevation, tri; psix(wtd, lagomorphs, spp diversity); p = setup, effort; px(psi)
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_px(psi)_coybob.R")
  #' start.time = Sys.time()
  #' coy.bob.global.px2 <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
  #'                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_px(psi)_coybob.txt",
  #'                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(coy.bob.global.px2$summary)
  #' print(coy.bob.global.px2$DIC)
  #' which(coy.bob.global.px2$summary[,"Rhat"] > 1.1)
  #' mcmcplot(coy.bob.global.px2$samples)
  #' save(coy.bob.global.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_px(psi)_coybob_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v1  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model: Habitat, with interaction
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort; px(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.).R")
  start.time = Sys.time()
  coy.bob.habx.px <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx.px$summary)
  print(coy.bob.habx.px$DIC)
  which(coy.bob.habx.px$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.habx.px$samples)
  save(coy.bob.habx.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_", Sys.Date(), ".RData"))
  
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_GoF.R")
  start.time = Sys.time()
  coy.bob.habx.px <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_GoF.txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx.px$summary)
  (coy.bob.habx.px_X2pB.coy <- mean(coy.bob.habx.px$sims.list$chi2.sim_A > coy.bob.habx.px$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  (coy.bob.habx.px_X2pB.bob <- mean(coy.bob.habx.px$sims.list$chi2.sim_B > coy.bob.habx.px$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  (coy.bob.habx.px_FTpB.coy <- mean(coy.bob.habx.px$sims.list$ft.sim_A > coy.bob.habx.px$sims.list$ft.obs_A)) # Bayesian p-value GOF
  (coy.bob.habx.px_FTpB.bob <- mean(coy.bob.habx.px$sims.list$ft.sim_B > coy.bob.habx.px$sims.list$ft.obs_B)) # Bayesian p-value GOF
  mean(coy.bob.habx.px$sims.list$chi2ratio_A); mean(coy.bob.habx.px$sims.list$chi2ratio_B); mean(coy.bob.habx.px$sims.list$ftratio_A); mean(coy.bob.habx.px$sims.list$ftratio_B)
  
  
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model:  Global
  #'  psi = setup, year, forest, elevation, tri; psix(wtd, lagomorphs, spp diversity); p = setup, effort; px(psi)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(psi).R")
  start.time = Sys.time()
  coy.bob.habx.px2 <- jags(bundled_pred_list[[6]], inits = inits.coy.bob, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(psi).txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(coy.bob.habx.px2$summary)
  print(coy.bob.habx.px2$DIC)
  which(coy.bob.habx.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(coy.bob.habx.px2$samples)
  save(coy.bob.habx.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(psi)_", Sys.Date(), ".RData"))
  
  
  
  #'  ----------------------
  ####  Bear-Coyote Models  ####
  #'  ----------------------
  inits.bear.coy <- function(){list(z = zinits[[7]], mean.psiSpp1 = runif(1, 0, 1),
                                    mean.psiSpp2 = runif(1, 0, 1), mean.pSpp1 = runif(1, 0, 1),
                                    mean.pSpp2 = runif(1, 0, 1))}
  ni <- 75000
  
  #####  Null model  ####
  #'  psi = year; p(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(yr)_p(.).R")
  start.time = Sys.time()
  bear.coy.null <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(yr)_p(.).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.null$summary)
  print(bear.coy.null$DIC)
  which(bear.coy.null$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.null$samples)
  save(bear.coy.null, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(yr)_p(.)_", Sys.Date(), ".RData"))
  
  #####  Habitat no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri; p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).R")
  start.time = Sys.time()
  bear.coy.hab <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                       "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_p(setup_effort).txt",
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.hab$summary)
  print(bear.coy.hab$DIC)
  which(bear.coy.hab$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.hab$samples)
  save(bear.coy.hab, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Prey abundance no inxs model  #### 
  #'  psi = setup, year, forest, elevation, tri, elk, wtd, lagomorphs; p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_bearcoy.R")
  start.time = Sys.time()
  bear.coy.preyabund <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                             "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preyabund_yr)_p(setup_effort)_bearcoy.txt",
                             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.preyabund$summary)
  print(bear.coy.preyabund$DIC)
  which(bear.coy.preyabund$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.preyabund$samples)
  save(bear.coy.preyabund, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_preyabund_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #' #####  Prey diversity no inxs model  #### 
  #' #'  psi = setup, year, forest, elevation, tri, spp diversity; p = setup, effort
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).R")
  #' start.time = Sys.time()
  #' bear.coy.preydiv <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
  #'                          "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_preydiversity_yr)_p(setup_effort).txt",
  #'                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(bear.coy.preydiv$summary)
  #' print(bear.coy.preydiv$DIC)
  #' which(bear.coy.preydiv$summary[,"Rhat"] > 1.1)
  #' mcmcplot(bear.coy.preydiv$samples)
  #' save(bear.coy.preydiv, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_preydiversity_yr)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(.); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).R")
  start.time = Sys.time()
  bear.coy.habx <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.habx$summary)
  print(bear.coy.habx$DIC)
  which(bear.coy.habx$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.habx$samples)
  save(bear.coy.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.R")
  start.time = Sys.time()
  bear.coy.habx <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF.txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.habx$summary)
  print(bear.coy.habx$DIC)
  which(bear.coy.habx$summary[,"Rhat"] > 1.1)
  (bear.coy.habx_X2pB.bear <- mean(bear.coy.habx$sims.list$chi2.sim_A > bear.coy.habx$sims.list$chi2.obs_A)) # Bayesian p-value GOF
  (bear.coy.habx_X2pB.coy <- mean(bear.coy.habx$sims.list$chi2.sim_B > bear.coy.habx$sims.list$chi2.obs_B)) # Bayesian p-value GOF
  (bear.coy.habx_FTpB.bear <- mean(bear.coy.habx$sims.list$ft.sim_A > bear.coy.habx$sims.list$ft.obs_A)) # Bayesian p-value GOF
  (bear.coy.habx_FTpB.coy <- mean(bear.coy.habx$sims.list$ft.sim_B > bear.coy.habx$sims.list$ft.obs_B)) # Bayesian p-value GOF
  mean(bear.coy.habx$sims.list$chi2ratio_A); mean(bear.coy.habx$sims.list$chi2ratio_B); mean(bear.coy.habx$sims.list$ftratio_A); mean(bear.coy.habx$sims.list$ftratio_B)
  # mcmcplot(bear.coy.habx$samples)
  save(bear.coy.habx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_GoF_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ habitat inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(forest, elevation, tri); p = setup, effort
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort).R")
  start.time = Sys.time()
  bear.coy.hab2x <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                        "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort).txt",
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.hab2x$summary)
  print(bear.coy.hab2x$DIC)
  which(bear.coy.hab2x$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.hab2x$samples)
  save(bear.coy.hab2x, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(habitat)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd, lagomorphs); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearcoy.R")
  start.time = Sys.time()
  bear.coy.preyabundx <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_bearcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.preyabundx$summary)
  print(bear.coy.preyabundx$DIC)
  which(bear.coy.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.preyabundx$samples)
  save(bear.coy.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Full habitat w/ prey abundance inx model  #### 
  #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd, lagomorphs, forest, elevation, tri); p = setup, effort 
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_bearcoy.R")
  start.time = Sys.time()
  bear.coy.hab.preyabundx <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                              "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_bearcoy.txt",
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.hab.preyabundx$summary)
  print(bear.coy.hab.preyabundx$DIC)
  which(bear.coy.hab.preyabundx$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.hab.preyabundx$samples)
  save(bear.coy.hab.preyabundx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(habitat_preyabund)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #' #####  Habitat w/ prey diversity inx model  #### 
  #' #'  psi = setup, year, forest, elevation, tri; psix(spp diversity); p = setup, effort
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).R")
  #' start.time = Sys.time()
  #' bear.coy.preydivx <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
  #'                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort).txt",
  #'                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(bear.coy.preydivx$summary)
  #' print(bear.coy.preydivx$DIC)
  #' which(bear.coy.preydivx$summary[,"Rhat"] > 1.1)
  #' mcmcplot(bear.coy.preydivx$samples)
  #' save(bear.coy.preydivx, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(preydiversity)_p(setup_effort)_", Sys.Date(), ".RData"))
  #' 
  #' #####  Global model  #### 
  #' #'  psi = setup, year, forest, elevation, tri; psix(elk, wtd, lagomorphs, spp diversity); p = setup, effort
  #' source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(global)_psix(global)_p(setup_effort)_bearcoy.R")
  #' start.time = Sys.time()
  #' bear.coy.global <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
  #'                         "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(global)_psix(global)_p(setup_effort)_bearcoy.txt",
  #'                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(bear.coy.global$summary)
  #' print(bear.coy.global$DIC)
  #' which(bear.coy.global$summary[,"Rhat"] > 1.1)
  #' mcmcplot(bear.coy.global$samples)
  #' save(bear.coy.global, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(global)_psix(global)_p(setup_effort)_", Sys.Date(), ".RData"))
  
  #####  Top model w/ intx on detection model v1  #### 
  #'  Parameterization tests whether detection of one predator affects detection of the other
  #'  Top model:  null
  #'  psi = year; p(.); px(.)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.).R")   
  start.time = Sys.time()
  bear.coy.null.px <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                           "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.).txt",    
                           n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.null.px$summary)
  print(bear.coy.null.px$DIC)
  which(bear.coy.null.px$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.null.px$samples)
  save(bear.coy.null.px, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_", Sys.Date(), ".RData")) 
  
  #####  Top model w/ intx on detection model v2  #### 
  #'  Parameterization tests whether presence of one predator affects detection of the other
  #'  Top model:  Habitat, no intx on psi
  #'  psi = year; p(.); px(psi)
  source("./Scripts/MultiSpp_OccMod/JAGS code/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(psi).R")   
  start.time = Sys.time()
  bear.coy.null.px2 <- jags(bundled_pred_list[[7]], inits = inits.bear.coy, params,
                            "./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/JAGS_code_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(psi).txt",
                            n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, DIC = TRUE, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(bear.coy.null.px2$summary)
  print(bear.coy.null.px2$DIC)
  which(bear.coy.null.px2$summary[,"Rhat"] > 1.1)
  mcmcplot(bear.coy.null.px2$samples)
  save(bear.coy.null.px2, file = paste0("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(psi)_", Sys.Date(), ".RData")) 
  
  