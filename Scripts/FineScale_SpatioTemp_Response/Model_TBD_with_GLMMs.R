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

  #'  Standardize covariates
  z_transform_covs <- function(tbd) {
    z_covs <- tbd %>%
      transmute(NewLocationID = NewLocationID, 
                Species = Species,
                Previous_Spp = as.factor(as.character(Previous_Spp)),
                GMU = factor(GMU, levels = c("GMU10A", "GMU6", "GMU1")),
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
    return(z_covs)
  }
  pred_tbd_zcov <- lapply(pred_tbd_short, z_transform_covs)
  
  
  #'  Call Tyra, we need our Next Top Model!
  
  #'  ---------------------------------------
  ####  Set up MCMC settings and run models  ####
  #'  ---------------------------------------
  #'  MCMC settigns
  nc <- 3
  ni <- 50000
  nt <- 10
  na <- 5000
  
  #'  Function to define and bundle data
  bundle_dat_data <- function(dat) {
    #'  Number of observations
    ntbd <- nrow(dat)
    #'  Number of unique camera locations
    ncams <- length(unique(dat$NewLocationID))
    #'  Format covariate data
    tbd_dat <- dat %>%
      transmute(cams = as.numeric(factor(NewLocationID), levels = NewLocationID), # must be 1 - n (not 0 - n) for nested indexing 
                CompetitorID = as.numeric(factor(Previous_Spp), levels = c("bear_black", "bobcat", "coyote", "mountain_lion", "wolf")), # must be 1-5 for nested indexing
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
    covs <- matrix(NA, ncol = 10, nrow = ntdb)
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
    newElk <- seq(from = min(tbd_dat$Nelk), to = max(tbd_dat$Nelk))
    newMoose <- seq(from = min(tbd_dat$Nmoose), to = max(tbd_dat$Nmoose))
    newWTD <- seq(from = min(tbd_dat$Nwtd), to = max(tbd_dat$Nwtd))
    newBunnies <- seq(from = min(tbd_dat$Nlagomorph), to = max(tbd_dat$Nlagomorph))
    newSppDiv <- seq(from = min(tbd_dat$SppDiversity), to = max(tbd_dat$SppDiversity))
    newcovs <- as.matrix(cbind(newElk, newMoose, newWTD, newBunnies, newSppDiv))
    
    #'  Number of covariates
    ncovs <- ncol(covs)
    
    #'  Time between detections
    tbd <- tbd_dat$TBD_mins
    print(summary(tbd))
    hist(tbd)
    
    bundled <- list(y = tbd, covs = covs, ncams = ncams, ncovs = ncovs, ntbd = ntbd,
                    site = tbd_dat$cams, newcovs = newcovs)
    
  }
  
  
  