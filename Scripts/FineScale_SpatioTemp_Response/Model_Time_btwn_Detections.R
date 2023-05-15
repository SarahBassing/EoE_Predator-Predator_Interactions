  #'  -------------------------------------------------
  #'  Competition & prey effects on latency of site use
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  -------------------------------------------------
  #'  Load up time-between-detection data and run exponential GLMMs in JAGS to
  #'  assess effect of prey availability on fine scale spatiotemporal interactions
  #'  between competing predators.
  #'  -------------------------------------------------
  
  #'  Clear workspace & load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(mcmcplots)
  library(AICcmodavg)
  library(tidyverse)
  
  #'  Read in data
  load("./Data/Time_btwn_Detections/TBD_all_predator_pairs_2023-05-05.RData")
  
  #'  Covariates
  load("./Data/Covariates_extracted/Covariates_EoE_Smr20.RData")
  load("./Data/Covariates_extracted/Covariates_EoE_Smr21.RData")
  eoe_covs_20s$Year <- "Smr20"
  eoe_covs_21s$Year <- "Smr21"
  
  #'  Join tbd & covariate data together
  tbd_w_covs_20s <- left_join(tbd_pred_pairs_all[tbd_pred_pairs_all$Year == "Smr20",], eoe_covs_20s, by = c("NewLocationID", "GMU", "Year")) 
  tbd_w_covs_21s <- left_join(tbd_pred_pairs_all[tbd_pred_pairs_all$Year == "Smr21",], eoe_covs_21s, by = c("NewLocationID", "GMU", "Year")) 
  tbd_w_covs <- rbind(tbd_w_covs_20s, tbd_w_covs_21s)
  
  #'  Filter by species
  bear <- filter(tbd_w_covs, Species == "bear_black") %>% filter(Predator_pair != "bear_black-bear_black")
  bob <- filter(tbd_w_covs, Species == "bobcat") %>% filter(Predator_pair != "bobcat-bobcat")
  coy <- filter(tbd_w_covs, Species == "coyote") %>% filter(Predator_pair != "coyote-coyote")
  lion <- filter(tbd_w_covs, Species == "mountain_lion") %>% filter(Predator_pair != "mountain_lion-mountain_lion")
  wolf <- filter(tbd_w_covs, Species == "wolf") %>% filter(Predator_pair != "wolf-wolf")
  
  #'  List data sets with extreme values removed
  pred_tbd <- list(bear, bob, coy, lion, wolf)
  
  #'  Function to identify potential outliers
  tbd_summary <- function(tbd, spp, quant) {
    #'  Plot frequency of time-between-detections (should look exponential)
    hist(tbd$DaysSinceLastDet, breaks = 50, main = paste("Number of days between detections for\n", spp))
    boxplot(tbd$DaysSinceLastDet, ylab = "Days", main = paste("Number of days between detections for\n", spp))
    
    #'  Review range of TBD values
    print("Quantiles of days between sequential detections of different predators")
    print(quantile(tbd$DaysSinceLastDet))
    #'  Review 90 - 100th quantiles- where are there gaps and outliers in the distribution?
    print("90th, 95th, and 99th quantiles of days between sequentail detections of predators") 
    print(quantile(tbd$DaysSinceLastDet, c(0.9, 0.95, 0.97, 0.99, 1.0)))
    
    #'  Re-plot frequency of time-btwn-detections after removing extreme values
    short_tbd <- filter(tbd, DaysSinceLastDet <= quantile(tbd$DaysSinceLastDet, c(quant)))
    hist(short_tbd$DaysSinceLastDet, breaks = 25, main = paste("Number of days between detections for\n", spp, "up to quantile =", quant))
    boxplot(short_tbd$DaysSinceLastDet, ylab = "Days", main = paste("Number of days between detections for\n", spp, "up to quantile =", quant))
    
    #'  Summary of observations with each predator species
    print("Total TBDs with each predator species")
    print(table(short_tbd$Species))
    #'  Summary of observations in each year
    print("Total TBDs for each year")
    print(table(short_tbd$Year))
    
    #'  Return dataset after removing extreme values
    return(short_tbd)
  }
  bear_short <- tbd_summary(bear, spp = "bear", quant = 0.99)
  bob_short <- tbd_summary(bob, spp = "bob", quant = 0.99)
  coy_short <- tbd_summary(coy, spp = "coy", quant = 0.99)
  lion_short <- tbd_summary(lion, spp = "lion", quant = 0.99)
  wolf_short <- tbd_summary(wolf, spp = "wolf", quant = 0.99)
  
  #'  List data sets with extreme values removed
  pred_tbd_short <- list(bear_short, bob_short, coy_short, lion_short, wolf_short)
  #' #'  Save for permutation test
  #' save(pred_tbd_short, file = paste0("./Data/Time_btwn_Detections/pred_tbd_short_", Sys.Date(), ".RData")) 
  
  #'  Table observations with each competitor to get a feel for sample size
  table(bear_short$Previous_Spp)
  table(bob_short$Previous_Spp)
  table(coy_short$Previous_Spp)
  table(lion_short$Previous_Spp)
  table(wolf_short$Previous_Spp)
  #'  Oof not a lot of observations for some of these species combos
  #'  Use coyote as indicator variable since has the most observations per species
  
  #'  Are there differences by GMU? 
  #'  Remember: GMUs 10A & 6 have 2 yrs data, GMU 1 only has 1 yr data
  gmu_tbd <- function(tbd) {
    tbd <- tbd %>% group_by(GMU) %>%
      summarise(n_tbd = n(),
                mean_tbd = round(mean(HoursSinceLastDet), 2),
                sd_tbd = round(sd(HoursSinceLastDet), 3),
                se_tbd = round(sd_tbd/sqrt(n_tbd))) %>%
      ungroup() %>%
      dplyr::select(-sd_tbd)
    print(tbd)
  }
  gmu_summary <- lapply(pred_tbd_short, gmu_tbd)
  
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
    print(summary(tbd_dat))
    print(head(tbd_dat))
    
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
    print(summary(tbd))
    hist(tbd)
    
    bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                    npp = npp, site = tbd_dat$cams, newcovs = newcovs)
    return(bundled)
    
  }
  #'  Provide specific order for CompetitorID levels - will differ for each species
  #'  Order generally goes black bear, bobcat, coyote, mountain lion, wolf
  bear_bundled <- bundle_dat_data(pred_tbd_short[[1]], npreyspp = 2, species_order = c("coyote", "bobcat", "mountain_lion", "wolf")) 
  bob_bundled <- bundle_dat_data(pred_tbd_short[[2]], npreyspp = 2, species_order = c("coyote", "bear_black", "mountain_lion", "wolf"))
  coy_bundled <- bundle_dat_data(pred_tbd_short[[3]], npreyspp = 2, species_order = c("bear_black", "bobcat", "mountain_lion", "wolf"))
  lion_bundled <- bundle_dat_data(pred_tbd_short[[4]], npreyspp = 2, species_order = c("coyote", "bear_black", "bobcat", "wolf"))
  wolf_bundled <- bundle_dat_data(pred_tbd_short[[5]], npreyspp = 3, species_order = c("coyote", "bear_black", "bobcat", "mountain_lion"))
  
  #' #'  Save for making figures later
  #' save(bear_bundled, file = "./Data/Time_btwn_Detections/bear_bundled.RData")
  #' save(bob_bundled, file = "./Data/Time_btwn_Detections/bob_bundled.RData")
  #' save(coy_bundled, file = "./Data/Time_btwn_Detections/coy_bundled.RData")
  #' save(lion_bundled, file = "./Data/Time_btwn_Detections/lion_bundled.RData")
  #' save(wolf_bundled, file = "./Data/Time_btwn_Detections/wolf_bundled.RData")
  
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
  
  
  #'  Call Tyra, we need our Next Top Model!
  
  
  #'  -----------------
  ####  BEAR Analyses  ####
  #'  -----------------
  #'  Setup initial values
  bear.init <- log(aggregate(bear_bundled$y, list(bear_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bear.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.null <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.null$summary)
  mcmcplot(tbd.bear.null$samples)
  save(tbd.bear.null, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_intercept_only.RData") 
    
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.compID <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.compID$summary)
  mcmcplot(tbd.bear.compID$samples)
  save(tbd.bear.compID, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_competitor_detection.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bobcat [2], lion [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.div <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.div$summary[1:5,])
  mcmcplot(tbd.bear.div$samples)
  save(tbd.bear.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.preyabund <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.preyabund$summary[1:5,])
  mcmcplot(tbd.bear.preyabund$samples)
  save(tbd.bear.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.compID.div <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_preydiversity_noRE.txt', 
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.compID.div$summary[1:15,])
  mcmcplot(tbd.bear.compID.div$samples)
  save(tbd.bear.compID.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_competitor_preydiv.RData") 
  
  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.compIDxdiv <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_preydiversity_noRE.txt', 
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.compIDxdiv$summary[1:15,])
  mcmcplot(tbd.bear.compIDxdiv$samples)
  save(tbd.bear.compIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_competitor_X_preydiv.RData") 
 
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.compID.preyabund <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_elk_wtd_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.compID.preyabund$summary[1:15,])
  mcmcplot(tbd.bear.compID.preyabund$samples)
  save(tbd.bear.compID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_competitor_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bobcat [2], lion [3], wolf [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.compIDxpreyabund <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_elk_wtd_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.compIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.bear.compIDxpreyabund$samples)
  save(tbd.bear.compIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_competitor_X_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bobcat [2], lion [3], wolf [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.global <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_elk_wtd_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.global$summary[1:21,])
  mcmcplot(tbd.bear.global$samples)
  save(tbd.bear.global, file = "./Outputs/Time_btwn_Detections/tbd.comp.bear_global.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bobcat [2], lion [3], wolf [4]
  
  
  #'  -----------------
  ####  BOBCAT Analyses  ####
  #'  -----------------
  #'  Setup initial values
  bob.init <- log(aggregate(bob_bundled$y, list(bob_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bob.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.null <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.null$summary)
  mcmcplot(tbd.bob.null$samples)
  save(tbd.bob.null, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.compID <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.compID$summary)
  mcmcplot(tbd.bob.compID$samples)
  save(tbd.bob.compID, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_detection.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.div <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.div$summary[1:5,])
  mcmcplot(tbd.bob.div$samples)
  save(tbd.bob.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.preyabund <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_wtd_lago_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.preyabund$summary[1:5,])
  mcmcplot(tbd.bob.preyabund$samples)
  save(tbd.bob.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_preyRAI.RData")
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.compID.div <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_preydiversity_noRE.txt', 
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.compID.div$summary[1:15,])
  mcmcplot(tbd.bob.compID.div$samples)
  save(tbd.bob.compID.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_preydiv.RData") 
  
  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.compIDxdiv <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_preydiversity_noRE.txt', 
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.compIDxdiv$summary[1:18,])
  mcmcplot(tbd.bob.compIDxdiv$samples)
  save(tbd.bob.compIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_X_preydiv.RData") # Not converging well, probably over-parameterized
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.compID.preyabund <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_wtd_lago_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.compID.preyabund$summary[1:15,])
  mcmcplot(tbd.bob.compID.preyabund$samples)
  save(tbd.bob.compID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.compIDxpreyabund <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_wtd_lago_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.compIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.bob.compIDxpreyabund$samples)
  save(tbd.bob.compIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_X_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.global <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_wtd_lago_abundance_noRE.txt', 
                                   inits = inits, n.chains = nc, n.iter = ni, 
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.global$summary[1:21,])
  mcmcplot(tbd.bob.global$samples)
  save(tbd.bob.global, file = "./Outputs/Time_btwn_Detections/tbd.comp.bob_global.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  
  #'  -------------------
  ####  COYOTE Analyses  ####
  #'  -------------------
  #'  Setup initial values
  coy.init <- log(aggregate(coy_bundled$y, list(coy_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = coy.init)}
  #'  NOTE: random effect for site excluded from these models owing to model failure when included
  #'  Error in checkForRemoteErrors(val): 3 nodes produced errors; first error: Error in node tbd_lambda[5] Invalid parent values
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.null <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.null$summary)
  mcmcplot(tbd.coy.null$samples)
  save(tbd.coy.null, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.compID <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_noRE.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.compID$summary)
  mcmcplot(tbd.coy.compID$samples)
  save(tbd.coy.compID, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_detection.RData") 
  #'  Keep in mind CompetitorID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.div <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.div$summary[1:5,])
  mcmcplot(tbd.coy.div$samples)
  save(tbd.coy.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.preyabund <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_wtd_lago_abundance_noRE.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.preyabund$summary[1:5,])
  mcmcplot(tbd.coy.preyabund$samples)
  save(tbd.coy.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.compID.div <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_preydiversity_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                             n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.compID.div$summary[1:15,])
  mcmcplot(tbd.coy.compID.div$samples)
  save(tbd.coy.compID.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_preydiv.RData") 
  
  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.compIDxdiv <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_preydiversity_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.compIDxdiv$summary[1:18,])
  mcmcplot(tbd.coy.compIDxdiv$samples)
  save(tbd.coy.compIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_X_preydiv.RData") 
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.compID.preyabund <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_wtd_lago_abundance_noRE.txt', 
                                   inits = inits, n.chains = nc, n.iter = ni, 
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.compID.preyabund$summary[1:15,])
  mcmcplot(tbd.coy.compID.preyabund$samples)
  save(tbd.coy.compID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.compIDxpreyabund <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_wtd_lago_abundance_noRE.txt', 
                                   inits = inits, n.chains = nc, n.iter = ni, 
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.compIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.coy.compIDxpreyabund$samples)
  save(tbd.coy.compIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_X_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.global <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_wtd_lago_abundance_noRE.txt', 
                                   inits = inits, n.chains = nc, n.iter = ni, 
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.global$summary[1:21,])
  mcmcplot(tbd.coy.global$samples)
  save(tbd.coy.global, file = "./Outputs/Time_btwn_Detections/tbd.comp.coy_global.RData") 
  #'  Keep in mind CompetitorID levels are bear [1], bobcat [2], lion [3], wolf [4]

  
  #'  --------------------------
  ####  MOUNTAIN LION Analyses  ####
  #'  --------------------------
  #'  Setup initial values
  lion.init <- log(aggregate(lion_bundled$y, list(lion_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = lion.init)} 
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.null <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.null$summary)
  mcmcplot(tbd.lion.null$samples)
  save(tbd.lion.null, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.compID <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_noRE.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.compID$summary[1:6,])
  mcmcplot(tbd.lion.compID$samples)
  save(tbd.lion.compID, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_competitor_detection.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.div <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.div$summary[1:5,])
  mcmcplot(tbd.lion.div$samples)
  save(tbd.lion.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.preyabund <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.preyabund$summary[1:5,])
  mcmcplot(tbd.lion.preyabund$samples)
  save(tbd.lion.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.compID.div <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_preydiversity_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                             n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.compID.div$summary[1:15,])
  mcmcplot(tbd.lion.compID.div$samples)
  save(tbd.lion.compID.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_competitor_preydiv.RData") 
  
  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.compIDxdiv <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_preydiversity_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.compIDxdiv$summary[1:18,])
  mcmcplot(tbd.lion.compIDxdiv$samples)
  save(tbd.lion.compIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_competitor_X_preydiv.RData") 
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.compID.preyabund <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_elk_wtd_abundance_noRE.txt', 
                                   inits = inits, n.chains = nc, n.iter = ni, 
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.compID.preyabund$summary[1:15,])
  mcmcplot(tbd.lion.compID.preyabund$samples)
  save(tbd.lion.compID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_competitor_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.compIDxpreyabund <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_elk_wtd_abundance_noRE.txt', 
                                   inits = inits, n.chains = nc, n.iter = ni, 
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.compIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.lion.compIDxpreyabund$samples)
  save(tbd.lion.compIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_competitor_X_preyRAI.RData") # Not converging well, probably over-parameterized
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.global <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_elk_wtd_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.global$summary[1:21,])
  mcmcplot(tbd.lion.global$samples)
  save(tbd.lion.global, file = "./Outputs/Time_btwn_Detections/tbd.comp.lion_global.RData") # Not converging well, probably over-parameterized
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  
  #'  -----------------
  ####  WOLF Analyses  ####
  #'  -----------------
  #'  Setup initial values
  wolf.init <- log(aggregate(wolf_bundled$y, list(wolf_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = wolf.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.null <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.null$summary)
  mcmcplot(tbd.wolf.null$samples)
  save(tbd.wolf.null, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.compID <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.compID$summary[1:10,])
  mcmcplot(tbd.wolf.compID$samples)
  save(tbd.wolf.compID, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_competitor_detection.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.div <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.div$summary[1:5,])
  mcmcplot(tbd.wolf.div$samples)
  save(tbd.wolf.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.preyabund <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_moose_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.preyabund$summary[1:5,])
  mcmcplot(tbd.wolf.preyabund$samples)
  save(tbd.wolf.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.compID.div <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_preydiversity_noRE.txt', 
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.compID.div$summary[1:15,])
  mcmcplot(tbd.wolf.compID.div$samples)
  save(tbd.wolf.compID.div, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_competitor_preydiv.RData") 
  
  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.compIDxdiv <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_preydiversity_noRE.txt', 
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.compIDxdiv$summary[1:15,])
  mcmcplot(tbd.wolf.compIDxdiv$samples)
  save(tbd.wolf.compIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_competitor_X_preydiv.RData") 
  
  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.compID.preyabund <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_elk_moose_wtd_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.compID.preyabund$summary[1:15,])
  mcmcplot(tbd.wolf.compID.preyabund$samples)
  save(tbd.wolf.compID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_competitor_preyRAI.RData") 
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_competitor_X_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.compIDxpreyabund <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_competitor_X_elk_moose_wtd_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.compIDxpreyabund$summary[1:21,])
  mcmcplot(tbd.wolf.compIDxpreyabund$samples)
  save(tbd.wolf.compIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_competitor_X_preyRAI.RData") # Not converging well, probably over-parameterized
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.global <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_elk_moose_wtd_abundance_noRE.txt', 
                                    inits = inits, n.chains = nc, n.iter = ni, 
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.global$summary[1:30,])
  mcmcplot(tbd.wolf.global$samples)
  save(tbd.wolf.global, file = "./Outputs/Time_btwn_Detections/tbd.comp.wolf_global.RData") # Not converging well, probably over-parameterized
  #'  Keep in mind CompetitorID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  
  #'  Fin
  #'  Next stop, TBD_Model_Selection_DIC.R for model selection and table formatting
  
  
  