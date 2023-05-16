  #'  ------------------------------------------
  #'  Permutation test for latency of site use
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ------------------------------------------
  #'  Scramble time-between-detections of different predator species and estimate
  #'  latency to then compare against observed TBD for each predator pairing.
  #'  Script generates 1000 permutations of each data set, then fits an intercept-only 
  #'  model and the top model identified in TBD_Model_Selection_DIC.R for each species. 
  #'  ------------------------------------------
  
  #'  Clear workspace, load libraries & read in data
  rm(list = ls())

  library(jagsUI)
  library(mcmcplots)
  library(AICcmodavg)
  library(tidyverse)
  
  # load("./Data/Time_btwn_Detections/TBD_all_predator_pairs_2023-05-05.RData")
  load("./Data/Time_btwn_Detections/pred_tbd_short_2023-05-15.RData")
  
  #'  ---------------------
  ####  Scramble TBD data  ####
  #'  ---------------------
  #'  Create empty lists to hold re-sampled data sets
  bear_list <- bob_list <- coy_list <- lion_list <- wolf_list <- list()
  
  #'  Re-sample TBD data to create new response variables for permutation test
  resample_tbd <- function(tbd_df, new_list, npermutations) {
    for(i in 1:npermutations) {
      new_tbd <- sample(tbd_df$TimeSinceLastDet, replace = TRUE)
      new_list[[i]] <- cbind(tbd_df, new_tbd)
    }
    return(new_list)
  }
  new_bear_tbd <- resample_tbd(pred_tbd_short[[1]], new_list = bear_list, npermutations = 1000)
  new_bob_tbd <- resample_tbd(pred_tbd_short[[2]], new_list = bob_list, npermutations = 1000)
  new_coy_tbd <- resample_tbd(pred_tbd_short[[3]], new_list = coy_list, npermutations = 1000)
  new_lion_tbd <- resample_tbd(pred_tbd_short[[4]], new_list = lion_list, npermutations = 1000)
  new_wolf_tbd <- resample_tbd(pred_tbd_short[[5]], new_list = wolf_list, npermutations = 1000)
  
  #'  List and save
  pred_tbd_short_resampled <- list(new_bear_tbd, new_bob_tbd, new_coy_tbd, new_lion_tbd, new_wolf_tbd)
  # save(pred_tbd_short_resampled, file = paste0("./Data/Time_btwn_Detections/pred_tbd_short_resampled_", Sys.Date(), ".RData"))
  load("./Data/Time_btwn_Detections/pred_tbd_short_resampled_2023-05-15.RData")
  
  #'  ---------------------------------------
  ####  Set up MCMC settings and run models  ####
  #'  ---------------------------------------
  #'  Function to define and bundle data
  bundle_dat_data <- function(dat, npreyspp, species_order) {
    #'  Number of observations
    ntbd <- nrow(dat)
    #'  Number of unique camera locations
    ncams <- length(unique(dat$NewLocationID))
    #'  Number of primary prey species
    npp <- npreyspp
    #'  Format covariate data
    tbd_dat <- dat %>%
      transmute(cams = as.numeric(factor(NewLocationID), levels = NewLocationID), # must be 1-n (not 0-n) for nested indexing 
                CompetitorID = as.numeric(factor(Previous_Spp), levels = species_order), # must be 1-4 for nested indexing
                GMU = as.numeric(factor(GMU), levels = c("GMU10A", "GMU6", "GMU1")),
                TBD_mins = TimeSinceLastDet,
                TBD_hrs = HoursSinceLastDet,
                TBD_days = DaysSinceLastDet,
                Elev = scale(Elevation__10m2), 
                PercForest = scale(perc_forest),
                SppDiversity = scale(H),
                Nelk = scale(elk_perday),
                Nmoose = scale(moose_perday),
                Nmd = scale(muledeer_perday),
                Nwtd = scale(whitetaileddeer_perday),
                Nlagomorph = scale(lagomorphs_perday),
                Nlivestock = scale(livestock_perday))
    
    #'  Covariate matrix for JAGS
    covs <- matrix(NA, ncol = 10, nrow = ntbd)
    covs[,1] <- tbd_dat$GMU
    covs[,2] <- tbd_dat$CompetitorID
    covs[,3] <- tbd_dat$Elev
    covs[,4] <- tbd_dat$PercForest
    covs[,5] <- tbd_dat$Nelk
    covs[,6] <- tbd_dat$Nmoose
    covs[,7] <- tbd_dat$Nmd
    covs[,8] <- tbd_dat$Nwtd
    covs[,9] <- tbd_dat$Nlagomorph
    covs[,10] <- tbd_dat$SppDiversity
    
    #'  Generate range of covariate values to predict across
    newElk <- seq(from = min(tbd_dat$Nelk), to = max(tbd_dat$Nelk), length.out = 100)
    newMoose <- seq(from = min(tbd_dat$Nmoose), to = max(tbd_dat$Nmoose), length.out = 100)
    newWTD <- seq(from = min(tbd_dat$Nwtd), to = max(tbd_dat$Nwtd), length.out = 100)
    newBunnies <- seq(from = min(tbd_dat$Nlagomorph), to = max(tbd_dat$Nlagomorph), length.out = 100)
    newSppDiv <- seq(from = min(tbd_dat$SppDiversity), to = max(tbd_dat$SppDiversity), length.out = 100)
    newcovs <- as.matrix(cbind(newElk, newMoose, newWTD, newBunnies, newSppDiv))
    
    #'  Number of covariates
    ncovs <- ncol(covs)
    
    #'  Time between detections
    tbd <- tbd_dat$TBD_mins
    
    bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                    npp = npp, site = tbd_dat$cams, newcovs = newcovs)
    return(bundled)
    
  }
  #'  Bundle formatted data for each of 1000 data sets per species
  new_bear_bundled <- lapply(pred_tbd_short_resampled[[1]], bundle_dat_data, npreyspp = 2, species_order = c("coyote", "bobcat", "mountain_lion", "wolf"))
  new_bob_bundled <- lapply(pred_tbd_short_resampled[[2]], bundle_dat_data, npreyspp = 2, species_order = c("coyote", "bear_black", "mountain_lion", "wolf"))
  new_coy_bundled <- lapply(pred_tbd_short_resampled[[3]], bundle_dat_data, npreyspp = 2, species_order = c("bear_black", "bobcat", "mountain_lion", "wolf"))
  new_lion_bundled <- lapply(pred_tbd_short_resampled[[4]], bundle_dat_data, npreyspp = 2, species_order = c("coyote", "bear_black", "bobcat", "wolf"))
  new_wolf_bundled <- lapply(pred_tbd_short_resampled[[5]], bundle_dat_data, npreyspp = 3, species_order = c("coyote", "bear_black", "bobcat", "mountain_lion"))
  
  #'  MCMC settings
  nc <- 3
  ni <- 25000
  nb <- 5000
  nt <- 10
  na <- 1000
  
  #'  Parameters to monitor
  params <- c("alpha0", "beta.competitor", "beta.prey", "beta.div", "beta.interaction", 
              "beta.interaction.elk", "beta.interaction.wtd", "beta.interaction.moose",
              "beta.interaction.lago", "sigma", "mu.tbd", "spp.tbd", "spp.tbd.elk", 
              "spp.tbd.moose", "spp.tbd.wtd", "spp.tbd.lago", "spp.tbd.div")
  

  #'  -------------------------
  ####  Intercept-only models  ####
  #'  -------------------------
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")

  null_mod_permutation <- function(bundled) {
    new.init <- log(aggregate(bundled$y, list(bundled$site), FUN = mean)[,2])
    inits <- function(){list(alpha = new.init)}
    tbd.null <- jags(bundled, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                     inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                     n.adapt = na, parallel = TRUE)
    return(tbd.null)
  } 
  #'  BEAR
  start.time <- Sys.time()
  tbd_null_bear_perm <- lapply(new_bear_bundled, null_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_null_bear_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_null_bear.RData")
  #'  BOBCAT
  start.time <- Sys.time()
  tbd_null_bob_perm <- lapply(new_bob_bundled, null_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_null_bob_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_null_bob.RData")
  #'  COYOTE
  start.time <- Sys.time()
  tbd_null_coy_perm <- lapply(new_coy_bundled, null_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_null_coy_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_null_coy.RData")
  #'  LION
  start.time <- Sys.time()
  tbd_null_lion_perm <- lapply(new_lion_bundled, null_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_null_lion_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_null_lion.RData")
  #'  WOLF
  start.time <- Sys.time()
  tbd_null_wolf_perm <- lapply(new_wolf_bundled, null_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_null_wolf_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_null_wolf.RData")
  
  
  #'  --------------
  ####  Top models  ####
  #'  --------------
  #####  Bear & Lion models (prey abundance)  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_wtd_abundance_noRE.R")
  
  preyabund_elkwtd_mod_permutation <- function(bundled) {
    new.init <- log(aggregate(bundled$y, list(bundled$site), FUN = mean)[,2])
    inits <- function(){list(alpha = new.init)}
    tbd.preyabund <- jags(bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt', 
                     inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                     n.adapt = na, parallel = TRUE)
    return(tbd.preyabund)
  } 
  #'  BEAR
  start.time <- Sys.time()
  tbd_preyabund_bear_perm <- lapply(new_bear_bundled, preyabund_elkwtd_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_preyabund_bear_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_preyabund_bear.RData")
  #'  LION
  start.time <- Sys.time()
  tbd_preyabund_lion_perm <- lapply(new_lion_bundled, preyabund_elkwtd_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_preyabund_lion_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_preyabund_lion.RData")
  
  
  #####  Wolf model (prey abundance)  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_moose_wtd_abundance_noRE.R")
  
  preyabund_elkmoosewtd_mod_permutation <- function(bundled) {
    new.init <- log(aggregate(bundled$y, list(bundled$site), FUN = mean)[,2])
    inits <- function(){list(alpha = new.init)}
    tbd.preyabund <- jags(bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_moose_wtd_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                          n.adapt = na, parallel = TRUE)
    return(tbd.preyabund)
  } 
  start.time <- Sys.time()
  tbd_preyabund_wolf_perm <- lapply(new_bear_bundled, preyabund_elkmoosewtd_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_preyabund_wolf_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_preyabund_wolf.RData")
  
  
  #####  Bobcat model (competitor ID + prey abundance)  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_wtd_lago_abundance_noRE.R")
  
  compID_preyabund_mod_permutation <- function(bundled) {
    new.init <- log(aggregate(bundled$y, list(bundled$site), FUN = mean)[,2])
    inits <- function(){list(alpha = new.init)}
    tbd.preyabund <- jags(bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_wtd_lago_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                          n.adapt = na, parallel = TRUE)
    return(tbd.preyabund)
  } 
  start.time <- Sys.time()
  tbd_compID_preyabund_bob_perm <- lapply(new_bear_bundled, compID_preyabund_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_compID_preyabund_bob_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_competitor_preyRAI_bob.RData")
  
  
  #####  Coyote model (competitor ID * prey abundance)  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_wtd_lago_abundance_noRE.R")
  
  compIDxpreyabund_mod_permutation <- function(bundled) {
    new.init <- log(aggregate(bundled$y, list(bundled$site), FUN = mean)[,2])
    inits <- function(){list(alpha = new.init)}
    tbd.preyabund <- jags(bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_wtd_lago_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                          n.adapt = na, parallel = TRUE)
    return(tbd.preyabund)
  } 
  start.time <- Sys.time()
  tbd_compIDxpreyabund_coy_perm <- lapply(new_bear_bundled, compIDxpreyabund_mod_permutation)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  save(tbd_compIDxpreyabund_coy_perm, file = "./Outputs/Time_btwn_Detections/permutation_tbd_competitor_X_preyRAI_coy.RData")
  
  #'  Load permutation runs
  load("./Outputs/Time_btwn_Detections/permutation_tbd_null_bear.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_null_bob.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_null_coy.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_null_lion.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_null_wolf.RData")
  
  load("./Outputs/Time_btwn_Detections/permutation_tbd_preyabund_bear.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_competitor_preyRAI_bob.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_competitor_X_preyRAI_coy.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_preyabund_lion.RData")
  load("./Outputs/Time_btwn_Detections/permutation_tbd_preyabund_wolf.RData")
  
  
  #'  --------------------
  ####  Permutation test  ####
  #'  --------------------
  #'  Grab posterior means from each permutation
  perm_coefs <- function(mod_out) {
    alpha0 <- as.numeric(mod_out$mean$alpha0) 
    mean_tbd <- as.numeric(mod_out$mean$mu.tbd) 
    mu_tbd <- cbind(alpha0, mean_tbd)
    mu_tbd <- as.data.frame(mu_tbd) 
    return(mu_tbd)
  }
  mu_bear_tbd <- lapply(tbd_null_bear_perm, perm_coefs); null_mu_tbd_bear <- as.data.frame(do.call(rbind, mu_bear_tbd))
  mu_bob_tbd <- lapply(tbd_null_bob_perm, perm_coefs); null_mu_tbd_bob <- as.data.frame(do.call(rbind, mu_bob_tbd))
  mu_coy_tbd <- lapply(tbd_null_coy_perm, perm_coefs); null_mu_tbd_coy <- as.data.frame(do.call(rbind, mu_coy_tbd))
  mu_lion_tbd <- lapply(tbd_null_lion_perm, perm_coefs); null_mu_tbd_lion <- as.data.frame(do.call(rbind, mu_lion_tbd))
  mu_wolf_tbd <- lapply(tbd_null_wolf_perm, perm_coefs); null_mu_tbd_wolf <- as.data.frame(do.call(rbind, mu_wolf_tbd))
  
  top_mu_bear_tbd <- lapply(tbd_preyabund_bear_perm, perm_coefs); top_mu_tbd_bear <- as.data.frame(do.call(rbind, top_mu_bear_tbd))
  top_mu_bob_tbd <- lapply(tbd_compID_preyabund_bob_perm, perm_coefs); top_mu_tbd_bob <- as.data.frame(do.call(rbind, top_mu_bob_tbd))
  top_mu_coy_tbd <- lapply(tbd_compIDxpreyabund_coy_perm, perm_coefs); top_mu_tbd_coy <- as.data.frame(do.call(rbind, top_mu_coy_tbd))
  top_mu_lion_tbd <- lapply(tbd_preyabund_lion_perm, perm_coefs); top_mu_tbd_lion <- as.data.frame(do.call(rbind, top_mu_lion_tbd))
  top_mu_wolf_tbd <- lapply(tbd_preyabund_wolf_perm, perm_coefs); top_mu_tbd_wolf <- as.data.frame(do.call(rbind, top_mu_wolf_tbd))
  
  #'  Load top models
  load("./Outputs/Time_btwn_Detections/tbd.comp.bear_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_X_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.lion_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.wolf_preyRAI.RData")
  
  #'  Load competitor models
  load("./Outputs/Time_btwn_Detections/tbd.comp.bear_competitor_detection.RData")
  # load("./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_detection.RData")
  # load("./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_detection.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.lion_competitor_detection.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.wolf_competitor_detection.RData")
  
  #'  Grab posterior means from original analyses
  perm_coefs_og <- function(mod_out) {
    alpha0 <- as.numeric(mod_out$mean$alpha0)
    mean_tbd <- as.numeric(mod_out$mean$mu.tbd)
    mu_tbd <- cbind(alpha0, mean_tbd)
    mu_tbd <- as.data.frame(mu_tbd)
    colnames(mu_tbd) <- c("alpha0", "mean.tbd")
    return(mu_tbd)
  }
  og_mu_bear_tbd <- perm_coefs_og(tbd.bear.preyabund)
  og_mu_lion_tbd <- perm_coefs_og(tbd.lion.preyabund)
  og_mu_wolf_tbd <- perm_coefs_og(tbd.wolf.preyabund)
  
  perm_coefs_compID <- function(mod_out, spp1tbd, spp2tbd, spp3tbd, spp4tbd) {
    alpha0 <- as.numeric(mod_out$mean$alpha0)
    mean_tbd <- as.numeric(mod_out$mean$mu.tbd)
    spp_tbd1 <- as.numeric(mod_out$mean$spp.tbd[[1]])
    spp_tbd2 <- as.numeric(mod_out$mean$spp.tbd[[2]])
    spp_tbd3 <- as.numeric(mod_out$mean$spp.tbd[[3]])
    spp_tbd4 <- as.numeric(mod_out$mean$spp.tbd[[4]])
    mu_tbd <- cbind(alpha0, mean_tbd, spp_tbd1, spp_tbd2, spp_tbd3, spp_tbd4)
    mu_tbd <- as.data.frame(mu_tbd)
    colnames(mu_tbd) <- c("alpha0", "mean.tbd", spp1tbd, spp2tbd, spp3tbd, spp4tbd)
    return(mu_tbd)
  }
  og_mu_bob_tbd <- perm_coefs_compID(tbd.bob.compID.preyabund, spp1tbd = "coy.tbd", spp2tbd = "bear.tbd", spp3tbd = "lion.tbd", spp4tbd = "wolf.tbd")
  og_mu_coy_tbd <- perm_coefs_compID(tbd.coy.compIDxpreyabund, spp1tbd = "bear.tbd", spp2tbd = "bob.tbd", spp3tbd = "lion.tbd", spp4tbd = "wolf.tbd")
  comp_mu_bear_tbd <- perm_coefs_compID(tbd.bear.compID, spp1tbd = "coy.tbd", spp2tbd = "bob.tbd", spp3tbd = "lion.tbd", spp4tbd = "wolf.tbd")
  comp_mu_lion_tbd <- perm_coefs_compID(tbd.lion.compID, spp1tbd = "coy.tbd", spp2tbd = "bear.tbd", spp3tbd = "bob.tbd", spp4tbd = "wolf.tbd")
  comp_mu_wolf_tbd <- perm_coefs_compID(tbd.wolf.compID, spp1tbd = "coy.tbd", spp2tbd = "bear.tbd", spp3tbd = "bob.tbd", spp4tbd = "lion.tbd")
  
  #'  Compare permutations to observed mean
  
    
  #'  Plot results
  hist(top_mu_tbd_bear$mean_tbd); abline(v = c(mean(top_mu_tbd_bear$mean_tbd), og_mu_bear_tbd$mean.tbd), col = c("red", "blue"), lty = c(2, 1), lwd = c(2, 2))
  
  