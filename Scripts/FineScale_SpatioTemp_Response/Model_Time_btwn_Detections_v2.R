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
  library(patchwork)
  
  #'  Read in data
  load("./Data/Time_btwn_Detections/TBD_all_spp-pred_pairs_2023-06-05.RData") #TBD_all_spp-pred_pairs_2023-05-23
  
  #'  Covariates
  load("./Data/Covariates_extracted/Covariates_EoE_Smr20_updated_072924.RData") # updated_070824 with TRI, different PercForest & RAI
  load("./Data/Covariates_extracted/Covariates_EoE_Smr21_updated_072924.RData")
  eoe_covs_20s$Year <- "Smr20"
  eoe_covs_21s$Year <- "Smr21"
  
  #'  Join tbd & covariate data together
  tbd_w_covs_20s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr20",], eoe_covs_20s, by = c("NewLocationID", "GMU", "Year")) 
  tbd_w_covs_21s <- left_join(tbd_spp_pairs_all[tbd_spp_pairs_all$Year == "Smr21",], eoe_covs_21s, by = c("NewLocationID", "GMU", "Year")) 
  tbd_w_covs <- rbind(tbd_w_covs_20s, tbd_w_covs_21s)
  
  #'  Filter by species
  bear <- filter(tbd_w_covs, Focal_predator == "bear_black") %>% filter(Species_pair != "bear_black-bear_black")
  bob <- filter(tbd_w_covs, Focal_predator == "bobcat") %>% filter(Species_pair != "bobcat-bobcat")
  coy <- filter(tbd_w_covs, Focal_predator == "coyote") %>% filter(Species_pair != "coyote-coyote")
  lion <- filter(tbd_w_covs, Focal_predator == "mountain_lion") %>% filter(Species_pair != "mountain_lion-mountain_lion")
  wolf <- filter(tbd_w_covs, Focal_predator == "wolf") %>% filter(Species_pair != "wolf-wolf")
  
  #'  List data sets with extreme values removed
  pred_tbd <- list(bear, bob, coy, lion, wolf)
  
  #'  Function to identify potential outliers and remove extreme observations
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
    hist(short_tbd$DaysSinceLastDet, breaks = 25, main = paste("Number of days between detections for\n", spp, "up to quantile =", quant)); abline(v = 20, col = "red", lty = 2)
    boxplot(short_tbd$DaysSinceLastDet, ylab = "Days", main = paste("Number of days between detections for\n", spp, "up to quantile =", quant)); abline(h = 20, col = "red", lty = 2)
    
    #'  Summary of observations with each predator species
    print("Total TBDs with each predator species")
    print(table(short_tbd$Focal_predator))
    #'  Summary of observations in each year
    print("Total TBDs for each year")
    print(table(short_tbd$Year))
    
    #' #'  Actually just remove any observations over 20 days long since don't expect most 
    #' #'  cues from previous predator to still be detectable beyond then
    #' short_tbd <- filter(short_tbd, DaysSinceLastDet <= 20)
    
    #'  Actually just remove any observations 7-days or longer since co-occurrence model
    #'  considered detections at summer and weekly time scales - want something finer
    #'  than the weekly time scale for this analysis
    short_tbd <- filter(short_tbd, DaysSinceLastDet < 7)
    
    #'  Return dataset after removing extreme values
    return(short_tbd)
  }
  bear_short <- tbd_summary(bear, spp = "bear", quant = 0.99) #0.99 to exclude extreme values (some past 20-day cutoff)
  bob_short <- tbd_summary(bob, spp = "bob", quant = 0.99)
  coy_short <- tbd_summary(coy, spp = "coy", quant = 0.99)
  lion_short <- tbd_summary(lion, spp = "lion", quant = 0.99)
  wolf_short <- tbd_summary(wolf, spp = "wolf", quant = 0.99)
  
  #'  List data sets with extreme values removed
  pred_tbd_short_list <- list(bear_short, bob_short, coy_short, lion_short, wolf_short)
  
  #'  Plot histograms of raw data
  pred_tbd_short_df <- rbind(bear_short, bob_short, coy_short, lion_short, wolf_short) %>%
    mutate(Focal_predator = ifelse(Focal_predator == "bear_black", "Black bear", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "bobcat", "Bobcat", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "coyote", "Coyote", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "mountain_lion", "Mountain lion", Focal_predator),
           Focal_predator = ifelse(Focal_predator == "wolf", "Wolf", Focal_predator))
  #'  Function to create plot for each predator species
  tbd_hist <- function(spp, pred_color) {
    dat <- filter(pred_tbd_short_df, Focal_predator == spp)
    plot_hist <- ggplot(dat, aes(x = HoursSinceLastDet)) + 
      geom_histogram(binwidth = 10, color = "black", fill = pred_color) + 
      theme_bw() + 
      xlab("Wait time (hours)") + 
      ylab("Frequency") +
      ggtitle(spp)
    print(plot_hist)
    return(plot_hist)
  }
  #'  List species and color associated with each one
  spp <- list("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")
  pred_color <- list("#98CAE1", "#A50026", "#DD3D2D", "#FDB366", "#364B9A")
  #'  Apply function to dataset
  histograms <- mapply(tbd_hist, spp, pred_color, SIMPLIFY = FALSE)
  
  #'  Create figure for publication
  tbd_histogram_fig <- histograms[[1]] + histograms[[2]] + theme(axis.title.y = element_blank()) + 
    histograms[[3]] + theme(axis.title.y = element_blank()) + 
    histograms[[4]] + histograms[[5]]  + theme(axis.title.y = element_blank()) + 
    plot_layout(ncol = 3) + plot_annotation(title = "Predator-specific wait times following detection of a different species",
                                            tag_levels = 'a')
  #'  Save
  ggsave("./Outputs/Time_btwn_Detections/Figures/TBD_histograms_rawdata.tiff", tbd_histogram_fig, 
         units = "in", width = 8, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Filter to focal predator species 
  focal_species <- function(tbd_dat) {
    tbd_dat <- tbd_dat %>%
      filter(Previous_Spp == "bear_black" | Previous_Spp == "bobcat" | Previous_Spp == "coyote" | 
               Previous_Spp == "mountain_lion" | Previous_Spp == "wolf")
    return(tbd_dat)
  }
  pred_tbd_short <- lapply(pred_tbd_short_list, focal_species)
  
  #'  Filter to focal predator species and handful of non-target species for comparison
  nontarget_species <- function(tbd_dat) {
    tbd_dat <- tbd_dat %>%
      filter(Previous_Spp == "elk" | Previous_Spp == "moose" | Previous_Spp == "whitetaileddeer" | Previous_Spp == "rabbit_hare") %>% 
      mutate(Previous_Spp = ifelse(Previous_Spp == "rabbit_hare", "lagomorph", Previous_Spp)) 
    return(tbd_dat)
  }
  prey_pred_tbd_short <- lapply(pred_tbd_short_list, nontarget_species)
  
  #' #'  Save for permutation test
  #' save(pred_tbd_short, file = paste0("./Data/Time_btwn_Detections/pred_tbd_short_", Sys.Date(), ".RData"))
  #' save(prey_pred_tbd_short, file = paste0("./Data/Time_btwn_Detections/prey_pred_tbd_short_", Sys.Date(), ".RData"))
 
  #'  Table observations with each competitor to get a feel for sample size
  print("bear"); table(pred_tbd_short[[1]]$Previous_Spp); table(prey_pred_tbd_short[[1]]$Previous_Spp)
  print("bobcat"); table(pred_tbd_short[[2]]$Previous_Spp); table(prey_pred_tbd_short[[2]]$Previous_Spp)
  print("coyote"); table(pred_tbd_short[[3]]$Previous_Spp); table(prey_pred_tbd_short[[3]]$Previous_Spp)
  print("lion"); table(pred_tbd_short[[4]]$Previous_Spp); table(prey_pred_tbd_short[[4]]$Previous_Spp)
  print("wolf"); table(pred_tbd_short[[5]]$Previous_Spp); table(prey_pred_tbd_short[[5]]$Previous_Spp)
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
    #'  Number of unique previously detected species
    nspp <- length(unique(dat$Previous_Spp))
    #'  Number of primary prey species
    npp <- npreyspp
    #'  Format covariate data
    tbd_dat <- dat %>%
      transmute(cams = as.numeric(factor(NewLocationID), levels = NewLocationID), # must be 1-n (not 0-n) for nested indexing 
                SpeciesID = as.numeric(factor(Previous_Spp), levels = species_order), # must be 1-n for nested indexing
                GMU = as.numeric(factor(GMU), levels = c("GMU10A", "GMU6", "GMU1")),
                TBD_mins = TimeSinceLastDet,
                TBD_hrs = HoursSinceLastDet,
                TBD_days = DaysSinceLastDet,
                Elev = scale(Elevation__10m2), 
                TRI = scale(TRI),
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
    covs <- matrix(NA, ncol = 11, nrow = ntbd)
    covs[,1] <- tbd_dat$GMU
    covs[,2] <- tbd_dat$SpeciesID
    covs[,3] <- tbd_dat$Elev
    covs[,4] <- tbd_dat$PercForest
    covs[,5] <- tbd_dat$Nelk
    covs[,6] <- tbd_dat$Nmoose
    covs[,7] <- tbd_dat$Nmd
    covs[,8] <- tbd_dat$Nwtd
    covs[,9] <- tbd_dat$Nlagomorph
    covs[,10] <- tbd_dat$SppDiversity
    covs[,11] <- tbd_dat$TRI
    
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
                    nspp = nspp, npp = npp, site = tbd_dat$cams, newcovs = newcovs)
    return(bundled)
    
  }
  #'  Provide specific order for SpeciesID levels - will differ for each species
  #'  Order generally goes black bear, bobcat, coyote, mountain lion, wolf
  bear_bundled <- bundle_dat_data(pred_tbd_short[[1]], npreyspp = 2, species_order = c("coyote", "bobcat", "mountain_lion", "wolf")) 
  bob_bundled <- bundle_dat_data(pred_tbd_short[[2]], npreyspp = 2, species_order = c("coyote", "bear_black", "mountain_lion", "wolf"))
  coy_bundled <- bundle_dat_data(pred_tbd_short[[3]], npreyspp = 2, species_order = c("bear_black", "bobcat", "mountain_lion", "wolf"))
  lion_bundled <- bundle_dat_data(pred_tbd_short[[4]], npreyspp = 2, species_order = c("coyote", "bear_black", "bobcat", "wolf"))
  wolf_bundled <- bundle_dat_data(pred_tbd_short[[5]], npreyspp = 3, species_order = c("coyote", "bear_black", "bobcat", "mountain_lion"))
  
  #'  Filter to only focal prey species for nontarget analysis
  bear_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[1]][prey_pred_tbd_short[[1]]$Previous_Spp == "elk" | prey_pred_tbd_short[[1]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("elk", "whitetaileddeer")) 
  bob_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[2]][prey_pred_tbd_short[[2]]$Previous_Spp == "lagomorph" | prey_pred_tbd_short[[2]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("lagomorph", "whitetaileddeer"))
  coy_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[3]][prey_pred_tbd_short[[3]]$Previous_Spp == "lagomorph" | prey_pred_tbd_short[[3]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("lagomorph", "whitetaileddeer"))
  lion_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[4]][prey_pred_tbd_short[[4]]$Previous_Spp == "elk" | prey_pred_tbd_short[[4]]$Previous_Spp == "whitetaileddeer",], npreyspp = 2, species_order = c("elk", "whitetaileddeer"))
  wolf_bundled_nontarget <- bundle_dat_data(prey_pred_tbd_short[[5]][prey_pred_tbd_short[[5]]$Previous_Spp != "lagomorph",], npreyspp = 3, species_order = c("elk", "moose", "whitetaileddeer"))
  
  #' #'  Save for making figures later
  #' save(bear_bundled, file = "./Data/Time_btwn_Detections/bear_bundled.RData")
  #' save(bob_bundled, file = "./Data/Time_btwn_Detections/bob_bundled.RData")
  #' save(coy_bundled, file = "./Data/Time_btwn_Detections/coy_bundled.RData")
  #' save(lion_bundled, file = "./Data/Time_btwn_Detections/lion_bundled.RData")
  #' save(wolf_bundled, file = "./Data/Time_btwn_Detections/wolf_bundled.RData")
  #' 
  #' save(bear_bundled_nontarget, file = "./Data/Time_btwn_Detections/bear_bundled_nontarget.RData")
  #' save(bob_bundled_nontarget, file = "./Data/Time_btwn_Detections/bob_bundled_nontarget.RData")
  #' save(coy_bundled_nontarget, file = "./Data/Time_btwn_Detections/coy_bundled_nontarget.RData")
  #' save(lion_bundled_nontarget, file = "./Data/Time_btwn_Detections/lion_bundled_nontarget.RData")
  #' save(wolf_bundled_nontarget, file = "./Data/Time_btwn_Detections/wolf_bundled_nontarget.RData")

  #'  MCMC settings
  nc <- 3
  ni <- 30000
  nb <- 5000
  nt <- 10
  na <- 1000
  
  #'  Parameters to monitor
  params <- c("alpha0", "beta.sppID", "beta.prey", "beta.div", "beta.interaction", 
              "beta.interaction.elk", "beta.interaction.wtd", "beta.interaction.moose",
              "beta.interaction.lago", "mu.tbd", "spp.tbd", "spp.tbd.elk", 
              "spp.tbd.moose", "spp.tbd.wtd", "spp.tbd.lago", "spp.tbd.div",
              "chi2.obs", "chi2.sim") #"sigma", 
  
  
  #'  Call Tyra, we need our Next Top Model!
  
  
  #'  ------------------------------
  ####  COMPETITOR - BEAR Analyses  ####
  #'  ------------------------------
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
  (tbd.bear.null.pval <- mean(tbd.bear.null$sims.list$chi2.sim > tbd.bear.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.null$samples)
  save(tbd.bear.null, file = "./Outputs/Time_btwn_Detections/tbd.bear_intercept_only.RData") 
  
  #####  Competitor model  ####    #DOUBLE CHECK GOF WITH MATT BEFORE IMPLIMENTING THROUGHOUT
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppID <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppID$summary)
  (tbd.bear.sppID.pval <- mean(tbd.bear.sppID$sims.list$chi2.sim > tbd.bear.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.sppID$samples)
  save(tbd.bear.sppID, file = "./Outputs/Time_btwn_Detections/tbd.bear_sppID.RData") 
  #'  Keep in mind SpeciesID levels are coyote[1], bobcat[2], lion[3], wolf[4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.div <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.div$summary[1:5,])
  (tbd.bear.div.pval <- mean(tbd.bear.div$sims.list$chi2.sim > tbd.bear.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.div$samples)
  save(tbd.bear.div, file = "./Outputs/Time_btwn_Detections/tbd.bear_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.preyabund <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.preyabund$summary[1:5,])
  (tbd.bear.preyabund.pval <- mean(tbd.bear.preyabund$sims.list$chi2.sim > tbd.bear.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.preyabund$samples)
  save(tbd.bear.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.bear_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppID.div <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_preydiversity_noRE.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppID.div$summary[1:15,])
  (tbd.bear.sppID.div.pval <- mean(tbd.bear.sppID.div$sims.list$chi2.sim > tbd.bear.sppID.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.sppID.div$samples)
  save(tbd.bear.sppID.div, file = "./Outputs/Time_btwn_Detections/tbd.bear_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppIDxdiv <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_preydiversity_noRE.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppIDxdiv$summary[1:15,])
  (tbd.bear.sppIDxdiv.pval <- mean(tbd.bear.sppIDxdiv$sims.list$chi2.sim > tbd.bear.sppIDxdiv$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.sppIDxdiv$samples)
  save(tbd.bear.sppIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.bear_sppID_X_preydiv.RData")

  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_elk_wtd_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppID.preyabund <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_elk_wtd_abundance_noRE.txt',
                                    inits = inits, n.chains = nc, n.iter = ni,
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppID.preyabund$summary[1:15,])
  (tbd.bear.sppID.preyabund.pval <- mean(tbd.bear.sppID.preyabund$sims.list$chi2.sim > tbd.bear.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.sppID.preyabund$samples)
  save(tbd.bear.sppID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.bear_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bobcat [2], lion [3], wolf [4]

  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_elk_wtd_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bear.sppIDxpreyabund <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_elk_wtd_abundance_noRE.txt',
                                    inits = inits, n.chains = nc, n.iter = ni,
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.sppIDxpreyabund$summary[1:30,])
  (tbd.bear.sppIDxpreyabund.pval <- mean(tbd.bear.sppIDxpreyabund$sims.list$chi2.sim > tbd.bear.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.sppIDxpreyabund$samples)
  save(tbd.bear.sppIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.bear_sppID_X_preyRAI.RData")    # Not converging well: likely over-parameterized
  #'  Keep in mind SpeciesID levels are coyote[1], bobcat[2], lion[3], wolf[4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_sppID_X_div_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bear.global <- jags(bear_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_sppID_X_div_elk_wtd_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bear.global$summary[1:21,])
  (tbd.bear.global.pval <- mean(tbd.bear.global$sims.list$chi2.sim > tbd.bear.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bear.global$samples)
  save(tbd.bear.global, file = "./Outputs/Time_btwn_Detections/tbd.bear_global.RData")   # Not converging well: likely over-parameterized
  #'  Keep in mind SpeciesID levels are coyote[1], bobcat[2], lion[3], wolf[4]
  
  
  #'  --------------------------------
  ####  COMPETITOR - BOBCAT Analyses  ####
  #'  --------------------------------
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
  (tbd.bob.null.pval <- mean(tbd.bob.null$sims.list$chi2.sim > tbd.bob.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.null$samples)
  save(tbd.bob.null, file = "./Outputs/Time_btwn_Detections/tbd.bob_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppID <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppID$summary)
  (tbd.bob.sppID.pval <- mean(tbd.bob.sppID$sims.list$chi2.sim > tbd.bob.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.sppID$samples)
  save(tbd.bob.sppID, file = "./Outputs/Time_btwn_Detections/tbd.bob_sppID.RData") 
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.div <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.div$summary[1:5,])
  (tbd.bob.div.pval <- mean(tbd.bob.div$sims.list$chi2.sim > tbd.bob.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.div$samples)
  save(tbd.bob.div, file = "./Outputs/Time_btwn_Detections/tbd.bob_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.preyabund <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_wtd_lago_abundance_noRE.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.preyabund$summary[1:5,])
  (tbd.bob.preyabund.pval <- mean(tbd.bob.preyabund$sims.list$chi2.sim > tbd.bob.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.preyabund$samples)
  save(tbd.bob.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.bob_preyRAI.RData")
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppID.div <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_preydiversity_noRE.txt',
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                             n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppID.div$summary[1:15,])
  (tbd.bob.sppID.div.pval <- mean(tbd.bob.sppID.div$sims.list$chi2.sim > tbd.bob.sppID.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.sppID.div$samples)
  save(tbd.bob.sppID.div, file = "./Outputs/Time_btwn_Detections/tbd.bob_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppIDxdiv <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_preydiversity_noRE.txt',
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppIDxdiv$summary[1:18,])
  (tbd.bob.sppIDxdiv.pval <- mean(tbd.bob.sppIDxdiv$sims.list$chi2.sim > tbd.bob.sppIDxdiv$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.sppIDxdiv$samples)
  save(tbd.bob.sppIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.bob_sppID_X_preydiv.RData")    # Not converging well: likely over-parameterized

  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_wtd_lago_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppID.preyabund <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_wtd_lago_abundance_noRE.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppID.preyabund$summary[1:15,])
  (tbd.bob.sppID.preyabund.pval <- mean(tbd.bob.sppID.preyabund$sims.list$chi2.sim > tbd.bob.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.sppID.preyabund$samples)
  save(tbd.bob.sppID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.bob_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]

  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_wtd_lago_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.bob.sppIDxpreyabund <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_wtd_lago_abundance_noRE.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.sppIDxpreyabund$summary[1:21,])
  (tbd.bob.sppIDxpreyabund.pval <- mean(tbd.bob.sppIDxpreyabund$sims.list$chi2.sim > tbd.bob.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.sppIDxpreyabund$samples)
  save(tbd.bob.sppIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.bob_sppID_X_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_sppID_X_div_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.bob.global <- jags(bob_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_sppID_X_div_wtd_lago_abundance_noRE.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.bob.global$summary[1:21,])
  (tbd.bob.global.pval <- mean(tbd.bob.global$sims.list$chi2.sim > tbd.bob.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.bob.global$samples)
  save(tbd.bob.global, file = "./Outputs/Time_btwn_Detections/tbd.bob_global.RData")   # Not converging well: likely over-parameterized
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], lion [3], wolf [4]
  
  
  #'  -------------------
  ####  COMPETITOR - COYOTE Analyses  ####
  #'  -------------------
  #'  Setup initial values
  coy.init <- log(aggregate(coy_bundled$y, list(coy_bundled$site), FUN = mean)[,2])
  inits <- function(){list(alpha = coy.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.null <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.null$summary); print(tbd.coy.null$DIC)
  mcmcplot(tbd.coy.null$samples)
  save(tbd.coy.null, file = "./Outputs/Time_btwn_Detections/tbd.coy_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppID <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppID$summary); print(tbd.coy.sppID$DIC)
  (tbd.coy.sppID.pval <- mean(tbd.coy.sppID$sims.list$chi2.sim > tbd.coy.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.sppID$samples)
  save(tbd.coy.sppID, file = "./Outputs/Time_btwn_Detections/tbd.coy_sppID.RData") 
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.div <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.div$summary[1:5,]); print(tbd.coy.div$DIC)
  (tbd.coy.div.pval <- mean(tbd.coy.div$sims.list$chi2.sim > tbd.coy.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.div$samples)
  save(tbd.coy.div, file = "./Outputs/Time_btwn_Detections/tbd.coy_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.preyabund <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_wtd_lago_abundance_noRE.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.preyabund$summary[1:5,]); print(tbd.coy.preyabund$DIC)
  (tbd.coy.preyabund.pval <- mean(tbd.coy.preyabund$sims.list$chi2.sim > tbd.coy.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.preyabund$samples)
  save(tbd.coy.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.coy_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppID.div <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_preydiversity_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                             n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppID.div$summary[1:15,]); print(tbd.coy.sppID.div$DIC)
  (tbd.coy.sppID.div.pval <- mean(tbd.coy.sppID.div$sims.list$chi2.sim > tbd.coy.sppID.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.sppID.div$samples)
  save(tbd.coy.sppID.div, file = "./Outputs/Time_btwn_Detections/tbd.coy_sppID_preydiv.RData") 
  
  #####  Species ID * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppIDxdiv <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_preydiversity_noRE.txt',
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppIDxdiv$summary[1:18,]); print(tbd.coy.sppIDxdiv$DIC)
  (tbd.coy.sppIDxdiv.pval <- mean(tbd.coy.sppIDxdiv$sims.list$chi2.sim > tbd.coy.sppIDxdiv$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.sppIDxdiv$samples)
  save(tbd.coy.sppIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.coy_sppID_X_preydiv.RData")

  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_wtd_lago_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppID.preyabund <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_wtd_lago_abundance_noRE.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppID.preyabund$summary[1:15,]); print(tbd.coy.sppID.preyabund$DIC)
  (tbd.coy.sppID.preyabund.pval <- mean(tbd.coy.sppID.preyabund$sims.list$chi2.sim > tbd.coy.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.sppID.preyabund$samples)
  save(tbd.coy.sppID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.coy_sppID_preyRAI.RData")
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]

  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_wtd_lago_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.coy.sppIDxpreyabund <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_wtd_lago_abundance_noRE.txt',
                                   inits = inits, n.chains = nc, n.iter = ni,
                                   n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.sppIDxpreyabund$summary[1:21,]); print(tbd.coy.sppIDxpreyabund$DIC)
  (tbd.coy.sppIDxpreyabund.pval <- mean(tbd.coy.sppIDxpreyabund$sims.list$chi2.sim > tbd.coy.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.sppIDxpreyabund$samples)
  save(tbd.coy.sppIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.coy_sppID_X_preyRAI.RData")  # chi2.sim has convergence issues
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_sppID_X_div_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.coy.global <- jags(coy_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_sppID_X_div_wtd_lago_abundance_noRE.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.coy.global$summary[1:21,]); print(tbd.coy.global$DIC)
  (tbd.coy.global.pval <- mean(tbd.coy.global$sims.list$chi2.sim > tbd.coy.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.coy.global$samples)
  save(tbd.coy.global, file = "./Outputs/Time_btwn_Detections/tbd.coy_global.RData") # chi2.sim has convergence issues
  #'  Keep in mind SpeciesID levels are bear [1], bobcat [2], lion [3], wolf [4]
  
  
  #'  ---------------------------------------
  ####  COMPETITOR - MOUNTAIN LION Analyses  ####
  #'  ---------------------------------------
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
  print(tbd.lion.null$summary); print(tbd.lion.null$DIC)
  (tbd.lion.null.pval <- mean(tbd.lion.null$sims.list$chi2.sim > tbd.lion.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.null$samples)
  save(tbd.lion.null, file = "./Outputs/Time_btwn_Detections/tbd.lion_intercept_only.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppID <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppID$summary); print(tbd.lion.sppID$DIC)
  (tbd.lion.sppID.pval <- mean(tbd.lion.sppID$sims.list$chi2.sim > tbd.lion.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.sppID$samples)
  save(tbd.lion.sppID, file = "./Outputs/Time_btwn_Detections/tbd.lion_sppID.RData") # Chi2.sim having issues converging but seems to settle on a value
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.div <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.div$summary[1:5,]); print(tbd.lion.div$DIC)
  (tbd.lion.div.pval <- mean(tbd.lion.div$sims.list$chi2.sim > tbd.lion.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.div$samples)
  save(tbd.lion.div, file = "./Outputs/Time_btwn_Detections/tbd.lion_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.preyabund <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.preyabund$summary[1:5,]); print(tbd.lion.preyabund$DIC)
  (tbd.lion.preyabund.pval <- mean(tbd.lion.preyabund$sims.list$chi2.sim > tbd.lion.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.preyabund$samples)
  save(tbd.lion.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.lion_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppID.div <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_preydiversity_noRE.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppID.div$summary[1:15,]); print(tbd.lion.sppID.div$DIC)
  (tbd.lion.sppID.div.pval <- mean(tbd.lion.sppID.div$sims.list$chi2.sim > tbd.lion.sppID.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.sppID.div$samples)
  save(tbd.lion.sppID.div, file = "./Outputs/Time_btwn_Detections/tbd.lion_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppIDxdiv <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_preydiversity_noRE.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppIDxdiv$summary[1:18,]); print(tbd.lion.sppIDxdiv$DIC)
  (tbd.lion.sppIDxdiv.pval <- mean(tbd.lion.sppIDxdiv$sims.list$chi2.sim > tbd.lion.sppIDxdiv$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.sppIDxdiv$samples)
  save(tbd.lion.sppIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.lion_sppID_X_preydiv.RData")    # Not converging well: likely over-parameterized

  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_elk_wtd_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppID.preyabund <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_elk_wtd_abundance_noRE.txt',
                                    inits = inits, n.chains = nc, n.iter = ni,
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppID.preyabund$summary[1:15,]); print(tbd.lion.sppID.preyabund$DIC)
  (tbd.lion.sppID.preyabund.pval <- mean(tbd.lion.sppID.preyabund$sims.list$chi2.sim > tbd.lion.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.sppID.preyabund$samples)
  save(tbd.lion.sppID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.lion_sppID_preyRAI.RData")  # Chi2.sim having issues converging
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]

  #####  Species ID * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_elk_wtd_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.lion.sppIDxpreyabund <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_elk_wtd_abundance_noRE.txt',
                                    inits = inits, n.chains = nc, n.iter = ni,
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.sppIDxpreyabund$summary[1:21,]); print(tbd.lion.sppIDxpreyabund$DIC)
  (tbd.lion.sppIDxpreyabund.pval <- mean(tbd.lion.sppIDxpreyabund$sims.list$chi2.sim > tbd.lion.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.sppIDxpreyabund$samples)
  save(tbd.lion.sppIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.lion_sppID_X_preyRAI.RData")     # Not converging well: likely over-parameterized
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_sppID_X_div_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.lion.global <- jags(lion_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_sppID_X_div_elk_wtd_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.lion.global$summary[1:21,]); print(tbd.lion.global$DIC)
  (tbd.lion.global.pval <- mean(tbd.lion.global$sims.list$chi2.sim > tbd.lion.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.lion.global$samples)
  save(tbd.lion.global, file = "./Outputs/Time_btwn_Detections/tbd.lion_global.RData")      # Not converging well: likely over-parameterized
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], wolf [4]
  
  
  #'  ------------------------------
  ####  COMPETITOR - WOLF Analyses  ####
  #'  ------------------------------
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
  print(tbd.wolf.null$summary); print(tbd.wolf.null$DIC)
  (tbd.wolf.null.pval <- mean(tbd.wolf.null$sims.list$chi2.sim > tbd.wolf.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.null$samples)
  save(tbd.wolf.null, file = "./Outputs/Time_btwn_Detections/tbd.wolf_intercept_only.RData") 
  
  #####  Competitor model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppID <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                          n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppID$summary); print(tbd.wolf.sppID$DIC)
  (tbd.wolf.sppID.pval <- mean(tbd.wolf.sppID$sims.list$chi2.sim > tbd.wolf.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.sppID$samples)
  save(tbd.wolf.sppID, file = "./Outputs/Time_btwn_Detections/tbd.wolf_sppID.RData") # Chi2.sim a little wonky
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.div <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.div$summary[1:5,]); print(tbd.wolf.div$DIC)
  (tbd.wolf.div.pval <- mean(tbd.wolf.div$sims.list$chi2.sim > tbd.wolf.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.div$samples)
  save(tbd.wolf.div, file = "./Outputs/Time_btwn_Detections/tbd.wolf_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.preyabund <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_elk_moose_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.preyabund$summary[1:5,]); print(tbd.wolf.preyabund$DIC)
  (tbd.wolf.preyabund.pval <- mean(tbd.wolf.preyabund$sims.list$chi2.sim > tbd.wolf.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.preyabund$samples)
  save(tbd.wolf.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.wolf_preyRAI.RData") 
  
  #####  Competitor + prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppID.div <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_preydiversity_noRE.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb,
                              n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppID.div$summary[1:15,]); print(tbd.wolf.sppID.div$DIC)
  (tbd.wolf.sppID.div.pval <- mean(tbd.wolf.sppID.div$sims.list$chi2.sim > tbd.wolf.sppID.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.sppID.div$samples)
  save(tbd.wolf.sppID.div, file = "./Outputs/Time_btwn_Detections/tbd.wolf_sppID_preydiv.RData")

  #####  Competitor * prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_preydiversity_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppIDxdiv <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_preydiversity_noRE.txt',
                              inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                              n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppIDxdiv$summary[1:15,]); print(tbd.wolf.sppIDxdiv$DIC)
  (tbd.wolf.sppIDxdiv.pval <- mean(tbd.wolf.sppIDxdiv$sims.list$chi2.sim > tbd.wolf.sppIDxdiv$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.sppIDxdiv$samples)
  save(tbd.wolf.sppIDxdiv, file = "./Outputs/Time_btwn_Detections/tbd.wolf_sppID_X_preydiv.RData")    # Chi2.sim not converging well

  #####  Competitor + prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_elk_moose_wtd_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppID.preyabund <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_elk_moose_wtd_abundance_noRE.txt',
                                    inits = inits, n.chains = nc, n.iter = ni,
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppID.preyabund$summary[1:25,]); print(tbd.wolf.sppID.preyabund$DIC)
  (tbd.wolf.sppID.preyabund.pval <- mean(tbd.wolf.sppID.preyabund$sims.list$chi2.sim > tbd.wolf.sppID.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.sppID.preyabund$samples)
  save(tbd.wolf.sppID.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.wolf_sppID_preyRAI.RData")   # Chi2.sim not converging well
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]

  #####  Competitor * prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID_X_elk_moose_wtd_abundance_noRE.R")

  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.sppIDxpreyabund <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_sppID_X_elk_moose_wtd_abundance_noRE.txt',
                                    inits = inits, n.chains = nc, n.iter = ni,
                                    n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.sppIDxpreyabund$summary[1:21,]); print(tbd.wolf.sppIDxpreyabund$DIC)
  (tbd.wolf.sppIDxpreyabund.pval <- mean(tbd.wolf.sppIDxpreyabund$sims.list$chi2.sim > tbd.wolf.sppIDxpreyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.sppIDxpreyabund$samples)
  save(tbd.wolf.sppIDxpreyabund, file = "./Outputs/Time_btwn_Detections/tbd.wolf_sppID_X_preyRAI.RData")      # Not converging well: likely over-parameterized
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_sppID_X_div_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.wolf.global <- jags(wolf_bundled, params, './Outputs/Time_btwn_Detections/tbd_global_sppID_X_div_elk_moose_wtd_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.wolf.global$summary[1:30,]); print(tbd.wolf.global$DIC)
  (tbd.wolf.global.pval <- mean(tbd.wolf.global$sims.list$chi2.sim > tbd.wolf.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.wolf.global$samples)
  save(tbd.wolf.global, file = "./Outputs/Time_btwn_Detections/tbd.wolf_global.RData")    # Not converging well: likely over-parameterized
  #'  Keep in mind SpeciesID levels are coyote [1], bear [2], bobcat [3], lion [4]
  
  
  #'  ------------------------------
  ####  NON-TARGET - BEAR Analyses  ####
  #'  ------------------------------
  #'  Setup initial values
  bear.nt.init <- log(aggregate(bear_bundled_nontarget$y, list(bear_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bear.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.null <- jags(bear_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.null$summary); print(tbd.nt.bear.null$DIC)
  (tbd.nt.bear.null.pval <- mean(tbd.nt.bear.null$sims.list$chi2.sim > tbd.nt.bear.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bear.null$samples)
  save(tbd.nt.bear.null, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bear_intercept_only.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.sppID <- jags(bear_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.sppID$summary); print(tbd.nt.bear.sppID$DIC)
  (tbd.nt.bear.sppID.pval <- mean(tbd.nt.bear.sppID$sims.list$chi2.sim > tbd.nt.bear.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bear.sppID$samples)
  save(tbd.nt.bear.sppID, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bear_sppID.RData") 
  #'  Keep in mind SppID levels are elk[1], white-tailed deer[2]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.div <- jags(bear_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.div$summary[1:5,]); print(tbd.nt.bear.div$DIC)
  (tbd.nt.bear.div.pval <- mean(tbd.nt.bear.div$sims.list$chi2.sim > tbd.nt.bear.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bear.div$samples)
  save(tbd.nt.bear.div, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bear_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.preyabund <- jags(bear_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.preyabund$summary[1:5,]); print(tbd.nt.bear.preyabund$DIC)
  (tbd.nt.bear.preyabund.pval <- mean(tbd.nt.bear.preyabund$sims.list$chi2.sim > tbd.nt.bear.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bear.preyabund$samples)
  save(tbd.nt.bear.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bear_preyRAI.RData")
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bear.global <- jags(bear_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_global_elk_wtd_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bear.global$summary[1:21,]); print(tbd.nt.bear.global$DIC)
  (tbd.nt.bear.global.pval <- mean(tbd.nt.bear.global$sims.list$chi2.sim > tbd.nt.bear.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bear.global$samples)
  save(tbd.nt.bear.global, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bear_global.RData") 
  #'  Keep in mind SpeciesID levels are elk[1], white-tailed deer[2]
  
  
  #'  --------------------------------
  ####  NON-TARGET - BOBCAT Analyses  ####
  #'  --------------------------------
  #'  Setup initial values
  bob.nt.init <- log(aggregate(bob_bundled_nontarget$y, list(bob_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = bob.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.null <- jags(bob_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.null$summary); print(tbd.nt.bob.null$DIC)
  (tbd.nt.bob.null.pval <- mean(tbd.nt.bob.null$sims.list$chi2.sim > tbd.nt.bob.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bob.null$samples)
  save(tbd.nt.bob.null, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bob_intercept_only.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.sppID <- jags(bob_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                        n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.sppID$summary); print(tbd.nt.bob.sppID$DIC)
  (tbd.nt.bob.sppID.pval <- mean(tbd.nt.bob.sppID$sims.list$chi2.sim > tbd.nt.bob.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bob.sppID$samples)
  save(tbd.nt.bob.sppID, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bob_sppID.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], white-tailed deer[2]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.div <- jags(bob_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.div$summary[1:5,]); print(tbd.nt.bob.div$DIC)
  (tbd.nt.bob.div.pval <- mean(tbd.nt.bob.div$sims.list$chi2.sim > tbd.nt.bob.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bob.div$samples)
  save(tbd.nt.bob.div, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bob_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.preyabund <- jags(bob_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_wtd_lago_abundance_noRE.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.preyabund$summary[1:5,]); print(tbd.nt.bob.preyabund$DIC)
  (tbd.nt.bob.preyabund.pval <- mean(tbd.nt.bob.preyabund$sims.list$chi2.sim > tbd.nt.bob.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bob.preyabund$samples)
  save(tbd.nt.bob.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bob_preyRAI.RData")
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.bob.global <- jags(bob_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_global_wtd_lago_abundance_noRE.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.bob.global$summary[1:21,]); print(tbd.nt.bob.global$DIC)
  (tbd.nt.bob.global.pval <- mean(tbd.nt.bob.global$sims.list$chi2.sim > tbd.nt.bob.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.bob.global$samples)
  save(tbd.nt.bob.global, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.bob_global.RData")  # GoF chi1.sim not converging well
  #'  Keep in mind SpeciesID levels are lagomorph[1], white-tailed deer[2]
  
  
  #'  --------------------------------
  ####  NON-TARGET - COYOTE Analyses  ####
  #'  --------------------------------
  #'  Setup initial values
  coy.nt.init <- log(aggregate(coy_bundled_nontarget$y, list(coy_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = coy.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.null <- jags(coy_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.null$summary); print(tbd.nt.coy.null$DIC)
  (tbd.nt.coy.null.pval <- mean(tbd.nt.coy.null$sims.list$chi2.sim > tbd.nt.coy.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.coy.null$samples)
  save(tbd.nt.coy.null, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.coy_intercept_only.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.sppID <- jags(coy_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                        n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.sppID$summary); print(tbd.nt.coy.sppID$DIC)
  (tbd.nt.coy.sppID.pval <- mean(tbd.nt.coy.sppID$sims.list$chi2.sim > tbd.nt.coy.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.coy.sppID$samples)
  save(tbd.nt.coy.sppID, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.coy_sppID.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], white-tailed deer[2]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.div <- jags(coy_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                      inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                      n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.div$summary[1:5,]); print(tbd.nt.coy.div$DIC)
  (tbd.nt.coy.div.pval <- mean(tbd.nt.coy.div$sims.list$chi2.sim > tbd.nt.coy.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.coy.div$samples)
  save(tbd.nt.coy.div, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.coy_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.preyabund <- jags(coy_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_wtd_lago_abundance_noRE.txt', 
                            inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                            n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.preyabund$summary[1:5,]); print(tbd.nt.coy.preyabund$DIC)
  (tbd.nt.coy.preyabund.pval <- mean(tbd.nt.coy.preyabund$sims.list$chi2.sim > tbd.nt.coy.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.coy.preyabund$samples)
  save(tbd.nt.coy.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.coy_preyRAI.RData") 
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_wtd_lago_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.coy.global <- jags(coy_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_global_wtd_lago_abundance_noRE.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, 
                         n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.coy.global$summary[1:21,]); print(tbd.nt.coy.global$DIC)
  (tbd.nt.coy.global.pval <- mean(tbd.nt.coy.global$sims.list$chi2.sim > tbd.nt.coy.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.coy.global$samples)
  save(tbd.nt.coy.global, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.coy_global.RData") 
  #'  Keep in mind SpeciesID levels are lagomorph[1], white-tailed deer[2]
  
  
  #'  ---------------------------------------
  ####  NON-TARGET - MOUNTAIN LION Analyses  ####
  #'  ---------------------------------------
  #'  Setup initial values
  lion.nt.init <- log(aggregate(lion_bundled_nontarget$y, list(lion_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = lion.nt.init)} 
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.null <- jags(lion_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.null$summary); print(tbd.nt.lion.null$DIC)
  (tbd.nt.lion.null.pval <- mean(tbd.nt.lion.null$sims.list$chi2.sim > tbd.nt.lion.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.lion.null$samples)
  save(tbd.nt.lion.null, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.lion_intercept_only.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.sppID <- jags(lion_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.sppID$summary); print(tbd.nt.lion.sppID$DIC)
  (tbd.nt.lion.sppID.pval <- mean(tbd.nt.lion.sppID$sims.list$chi2.sim > tbd.nt.lion.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.lion.sppID$samples)
  save(tbd.nt.lion.sppID, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.lion_sppID.RData") 
  #'  Keep in mind SpeciesID levels are elk[1], white-tailed deer[2]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.div <- jags(lion_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.div$summary[1:5,]); print(tbd.nt.lion.div$DIC)
  (tbd.nt.lion.div.pval <- mean(tbd.nt.lion.div$sims.list$chi2.sim > tbd.nt.lion.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.lion.div$samples)
  save(tbd.nt.lion.div, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.lion_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.preyabund <- jags(lion_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_elk_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.preyabund$summary[1:5,]); print(tbd.nt.lion.preyabund$DIC)
  (tbd.nt.lion.preyabund.pval <- mean(tbd.nt.lion.preyabund$sims.list$chi2.sim > tbd.nt.lion.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.lion.preyabund$samples)
  save(tbd.nt.lion.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.lion_preyRAI.RData")
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_elk_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.lion.global <- jags(lion_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_global_elk_wtd_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.lion.global$summary[1:21,]); print(tbd.nt.lion.global$DIC)
  (tbd.nt.lion.global.pval <- mean(tbd.nt.lion.global$sims.list$chi2.sim > tbd.nt.lion.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.lion.global$samples)
  save(tbd.nt.lion.global, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.lion_global.RData") 
  #'  Keep in mind SpeciesID levels are elk[1], white-tailed deer[2]
  
  
  #'  ------------------------------
  ####  NON-TARGET - WOLF Analyses  ####
  #'  ------------------------------
  #'  Setup initial values
  wolf.nt.init <- log(aggregate(wolf_bundled_nontarget$y, list(wolf_bundled_nontarget$site), FUN = mean)[,2])
  inits <- function(){list(alpha = wolf.nt.init)}
  
  #####  Null model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_intercept_only.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.null <- jags(wolf_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_intercept_only.txt', 
                        inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                        n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.null$summary); print(tbd.nt.wolf.null$DIC)
  (tbd.nt.wolf.null.pval <- mean(tbd.nt.wolf.null$sims.list$chi2.sim > tbd.nt.wolf.null$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.wolf.null$samples)
  save(tbd.nt.wolf.null, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_intercept_only.RData") 
  
  #####  Species ID model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_sppID.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.sppID <- jags(wolf_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_sppID.txt', 
                         inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, 
                         n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.sppID$summary); print(tbd.nt.wolf.sppID$DIC)
  (tbd.nt.wolf.sppID.pval <- mean(tbd.nt.wolf.sppID$sims.list$chi2.sim > tbd.nt.wolf.sppID$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.wolf.sppID$samples)
  save(tbd.nt.wolf.sppID, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_sppID.RData") 
  #'  Keep in mind SpeciesID levels are elk[1], moose[2], white-tailed deer[3]
  
  #####  Prey diversity model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_preydiversity_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.div <- jags(wolf_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_preydiversity_noRE.txt', 
                       inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                       n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.div$summary[1:5,]); print(tbd.nt.wolf.div$DIC)
  (tbd.nt.wolf.div.pval <- mean(tbd.nt.wolf.div$sims.list$chi2.sim > tbd.nt.wolf.div$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.wolf.div$samples)
  save(tbd.nt.wolf.div, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_preydiversity.RData") 
  
  #####  Prey relative abundance model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.preyabund <- jags(wolf_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_elk_moose_wtd_abundance_noRE.txt', 
                             inits = inits, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt,
                             n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.preyabund$summary[1:5,]); print(tbd.nt.wolf.preyabund$DIC)
  (tbd.nt.wolf.preyabund.pval <- mean(tbd.nt.wolf.preyabund$sims.list$chi2.sim > tbd.nt.wolf.preyabund$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.wolf.preyabund$samples)
  save(tbd.nt.wolf.preyabund, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_preyRAI.RData") 
  
  #####  Global model  ####
  source("./Scripts/FineScale_SpatioTemp_Response/JAGS_models/JAGS_tbd_global_elk_moose_wtd_abundance_noRE.R")
  
  #'  Run model
  start.time <- Sys.time()
  tbd.nt.wolf.global <- jags(wolf_bundled_nontarget, params, './Outputs/Time_btwn_Detections/tbd_global_elk_moose_wtd_abundance_noRE.txt', 
                          inits = inits, n.chains = nc, n.iter = ni, 
                          n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(tbd.nt.wolf.global$summary[1:30,]); print(tbd.nt.wolf.global$DIC)
  (tbd.nt.wolf.global.pval <- mean(tbd.nt.wolf.global$sims.list$chi2.sim > tbd.nt.wolf.global$sims.list$chi2.obs)) # Bayesian p-value GOF
  mcmcplot(tbd.nt.wolf.global$samples)
  save(tbd.nt.wolf.global, file = "./Outputs/Time_btwn_Detections/tbd.nontarget.wolf_global.RData") 
  #'  Keep in mind SpeciesID levels are elk[1], moose[2], white-tailed deer[3]

  
  #'  Fin
  #'  Next stop, TBD_Model_Selection_DIC.R for model selection and table formatting
  
  
