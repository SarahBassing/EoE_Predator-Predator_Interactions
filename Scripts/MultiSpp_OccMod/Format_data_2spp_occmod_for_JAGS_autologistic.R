  #'  --------------------------------------------------------------
  #'  Format Data for 2-spp Autologistic Occupancy Models (for JAGS) 
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  --------------------------------------------------------------
  #'  Script to format covariate data for multi-species occupancy models with 
  #'  EoE data. Focusing on black bear, bobcat, coyote, mountain lion, and wolf 
  #'  interactions. Cameras ran summer 2020, winter 2020-2021, and summer 2021 but
  #'  focusing on summer data only. Co-occurrence models test whether predator
  #'  co-occurrence is non-independent and whether their occurrence, co-occurrence,
  #'  and detection are influenced by other variables of interest.
  #'  
  #'  Summer primary period is considered July 1 - Sept 15, equating to 11 1-wk 
  #'  sampling periods. 
  #'  Winter primary period is considered Dec 1 - Feb 1, equating to 9 1-wk sampling periods
  #'   
  #'  
  #'  Encounter histories are generated with the Detection_histories_for_occmod.R
  #'  Covariate data were extracted and formatted with Covariate_Extract.R
  #'  --------------------------------------------
  
  #'  Load libraries
  library(unmarked)
  library(stringr)
  library(tidyverse)
  
  #'  Load and bind annual detection histories
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")
  
  #'  Combine annual detection histories
  all_detections <- function(dh1, dh2) {
    dh1 <- dh1
    dh2 <- dh2
    
    dh1 <- as.data.frame(dh1) %>% rownames_to_column("NewLocationID")
    dh2 <- as.data.frame(dh2) %>% rownames_to_column("NewLocationID")
    
    #'  Join annual data together so each site has a complete detection history
    dh_full <- full_join(dh1, dh2, by = "NewLocationID") %>%
      arrange(NewLocationID)
    
    #'  Split back into year 1 & 2 dh
    dh1 <- dh_full[,1:12]; dh1 <- column_to_rownames(dh1, var = "NewLocationID")
    dh2 <- dh_full[,c(1,13:23)]; dh2 <- column_to_rownames(dh2, var = "NewLocationID")
    
    #'  Rename columns in each detection history
    newcols <- c("occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11")
    colnames(dh1) <- newcols; colnames(dh2) <- newcols
    #'  Convert back to matrices
    dh1 <- as.matrix(dh1); dh2 <- as.matrix(dh2)
    
    #'  Convert annual detection histories into a single array
    dh <- array(data = c(dh1, dh2), dim = c(727, 11, 2), dimnames = list(rownames(dh1), 
                                                                         colnames(dh1), 
                                                                         c("year1", "year2")))
    
    return(dh)
  }
  #'  Generate detection history array
  dh <- all_detections(DH_eoe20s_predators[[1]][[1]], DH_eoe21s_predators[[1]][[1]])
  
  #'  Load camera station and covariate data
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(NewLocationID != "GMU6_U_160" | CameraHeight_M != 1.2) #' Remove duplicate camera where height changed slightly
  
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  #'  Reformat camera deployment data
  format_cam_station <- function(cams, season) {
    cams <- dplyr::select(cams, c("NewLocationID", "Lat", "Long")) %>%
      mutate(GMU = str_extract(NewLocationID, "[^_]+"),
             Season = season) %>%
      left_join(cams_eoe_long[cams_eoe_long$Season == season,], by = c("NewLocationID", "Season")) %>%
      #'  Drop handful of duplicated rows (not sure why this happens with left_join)
      unique() %>%
      #'  Reduce categories to random, road, trail
      mutate(CameraFacing = ifelse(CameraFacing == "road" | CameraFacing == "atv" | CameraFacing == "gravel" | 
                                     CameraFacing == "decommision" | CameraFacing == "decommission", "road", CameraFacing),
             CameraFacing = ifelse(CameraFacing == "hiking" | CameraFacing == "game" | CameraFacing == "other", "trail", CameraFacing),
             CameraFacing = ifelse(is.na(CameraFacing), "trail", CameraFacing)) %>% # predator cameras so either trail or road
      arrange(NewLocationID)
    return(cams)
  }
  cams_eoe20s <- format_cam_station(eoe_probcams_20s, season = "Smr20")
  cams_eoe21s <- format_cam_station(eoe_probcams_21s, season = "Smr21") 
  cams_eoe20s21s <- full_join(cams_eoe20s, cams_eoe21s, by = c("NewLocationID")) %>%
    arrange(NewLocationID) %>%
    mutate(Lat.x = ifelse(is.na(Lat.x), Lat.y, Lat.x),
           Lat.y = ifelse(is.na(Lat.y), Lat.x, Lat.y),
           Long.x = ifelse(is.na(Long.x), Long.y, Long.x),
           Long.y = ifelse(is.na(Long.y), Long.x, Long.y),
           GMU.x = ifelse(is.na(GMU.x), GMU.y, GMU.x),
           GMU.y = ifelse(is.na(GMU.y), GMU.x, GMU.y),
           Gmu.x = ifelse(is.na(Gmu.x), Gmu.y, Gmu.x),
           Gmu.y = ifelse(is.na(Gmu.y), Gmu.x, Gmu.y),
           Season.x = "Smr20",
           Season.y = "Smr21",
           Setup.x = ifelse(is.na(Setup.x), Setup.y, Setup.x),
           Setup.y = ifelse(is.na(Setup.y), Setup.x, Setup.y),
           CameraHeight_M.x = ifelse(is.na(CameraHeight_M.x), CameraHeight_M.y, CameraHeight_M.x),
           CameraHeight_M.y = ifelse(is.na(CameraHeight_M.y), CameraHeight_M.x, CameraHeight_M.y),
           CameraFacing.x = ifelse(is.na(CameraFacing.x), CameraFacing.y, CameraFacing.x),
           CameraFacing.y = ifelse(is.na(CameraFacing.y), CameraFacing.x, CameraFacing.y))
  
  #'  Split back into year 1 & 2
  cams_eoe20s <- cams_eoe20s21s[,1:9]; 
  cams_eoe21s <- cams_eoe20s21s[,c(1,10:17)]
  #'  Rename columns
  newcols <- c("NewLocationID", "Lat", "Long", "GMU", "Season", "Gmu", "Setup", "CameraHeight_M", "CameraFacing")
  colnames(cams_eoe20s) <- newcols; colnames(cams_eoe21s) <- newcols
  
  #'  Load extracted covariate data
  load("./Data/Covariates_extracted/Covariates_EoE_Smr20.RData")
  load("./Data/Covariates_extracted/Covariates_EoE_Smr21.RData")
  
  
  #'  Scale and format site-level covariates
  format_covs <- function(cams_yr1, cams_yr2, covs_yr1, covs_yr2, rm_rows_yr1, rm_rows_yr2) {
    #'  Join camera data with extracted covariate data
    cam_covs1 <- full_join(cams_yr1, covs_yr1) %>% arrange(NewLocationID)
    cam_covs2 <- full_join(cams_yr2, covs_yr2) %>% arrange(NewLocationID)
    # cam_covs1 <- cam_covs1 %>%
    #   mutate(perc_forest = ifelse(is.na(perc_forest), cam_covs2$perc_forest, perc_forest), 
    #          Elevation__10m2 = ifelse(is.na(Elevation__10m2), cam_covs2$Elevation__10m2, Elevation__10m2),
    #          SR = ifelse(is.na(SR), cam_covs2$SR, SR),
    #          H = ifelse(is.na(H), cam_covs2$H, H),
    #          dominantprey = ifelse(is.na(dominantprey), cam_covs2$dominantprey, dominantprey),
    #          elk_perday = ifelse(is.na(elk_perday), cam_covs2$elk_perday, elk_perday),
    #          lagomorphs_perday = ifelse(is.na(lagomorphs_perday), cam_covs2$lagomorphs_perday, lagomorphs_perday),
    #          livestock_perday = ifelse(is.na(livestock_perday), cam_covs2$livestock_perday, livestock_perday),
    #          moose_perday = ifelse(is.na(moose_perday), cam_covs2$moose_perday, moose_perday),
    #          muledeer_perday = ifelse(is.na(muledeer_perday), cam_covs2$muledeer_perday, muledeer_perday),
    #          whitetaileddeer_perday = ifelse(is.na(whitetaileddeer_perday), cam_covs2$whitetaileddeer_perday, whitetaileddeer_perday))
    # cam_covs2 <- cam_covs2 %>%
    #   mutate(perc_forest = ifelse(is.na(perc_forest), cam_covs1$perc_forest, perc_forest), 
    #          Elevation__10m2 = ifelse(is.na(Elevation__10m2), cam_covs1$Elevation__10m2, Elevation__10m2),
    #          SR = ifelse(is.na(SR), cam_covs1$SR, SR),
    #          H = ifelse(is.na(H), cam_covs1$H, H),
    #          dominantprey = ifelse(is.na(dominantprey), cam_covs1$dominantprey, dominantprey),
    #          elk_perday = ifelse(is.na(elk_perday), cam_covs1$elk_perday, elk_perday),
    #          lagomorphs_perday = ifelse(is.na(lagomorphs_perday), cam_covs1$lagomorphs_perday, lagomorphs_perday),
    #          livestock_perday = ifelse(is.na(livestock_perday), cam_covs1$livestock_perday, livestock_perday),
    #          moose_perday = ifelse(is.na(moose_perday), cam_covs1$moose_perday, moose_perday),
    #          muledeer_perday = ifelse(is.na(muledeer_perday), cam_covs1$muledeer_perday, muledeer_perday),
    #          whitetaileddeer_perday = ifelse(is.na(whitetaileddeer_perday), cam_covs1$whitetaileddeer_perday, whitetaileddeer_perday))

    #'  Bind annual covariate data together
    cam_covs <- rbind(cam_covs1, cam_covs2)
    #'  Rename, format, and scale as needed
    formatted <- transmute(cam_covs,
                           NewLocationID = NewLocationID,
                           Season = as.factor(as.character(Season)),
                           GMU = GMU,
                           CameraFacing = as.factor(as.character(CameraFacing)),
                           Setup = as.factor(as.character(Setup)),
                           DomPrey = as.factor(dominantprey),
                           Height = scale(CameraHeight_M),
                           PercForest = scale(perc_forest), 
                           Elev = scale(Elevation__10m2),
                           SppRich = scale(SR),
                           SppDiversity = scale(H),
                           Nelk = scale(elk_perday),    
                           Nlagomorph = scale(lagomorphs_perday),
                           Nlivestock = scale(livestock_perday),
                           Nmoose = scale(moose_perday),
                           Nmd = scale(muledeer_perday),
                           Nwtd = scale(whitetaileddeer_perday)) 
    
    #'  Adjust reference category for CameraFacing factors
    order_camfacing <- c("random", "trail", "road")
    formatted <- formatted %>%
      mutate(CameraFacing = fct_relevel(CameraFacing, order_camfacing))
    
    #'  Split into annual covariate datasets
    scaled_covs1 <- filter(formatted, Season == "Smr20")
    scaled_covs2 <- filter(formatted, Season == "Smr21")
    
    #'  Fill in missing values for cameras that were only deployed/operational in
    #'  one year but not the other
    #'  NOTE: for prey per day, this assumes average prey activity at each site was the same across years!!! 
    #'  NOTE: filling in NAs AFTER standardizing- meaning these duplicates are not included when scaling above
    scaled_covs1 <- scaled_covs1 %>%
      mutate(PercForest = ifelse(is.na(PercForest), scaled_covs2$PercForest, PercForest), 
             Elev = ifelse(is.na(Elev), scaled_covs2$Elev, Elev),
             SppRich = ifelse(is.na(SppRich), scaled_covs2$SppRich, SppRich),
             SppDiversity = ifelse(is.na(SppDiversity), scaled_covs2$SppDiversity, SppDiversity),
             DomPrey = ifelse(is.na(DomPrey), scaled_covs2$DomPrey, DomPrey),
             Nelk = ifelse(is.na(Nelk), scaled_covs2$Nelk, Nelk),
             Nlagomorph = ifelse(is.na(Nlagomorph), scaled_covs2$Nlagomorph, Nlagomorph),
             Nlivestock = ifelse(is.na(Nlivestock), scaled_covs2$Nlivestock, Nlivestock),
             Nmoose = ifelse(is.na(Nmoose), scaled_covs2$Nmoose, Nmoose),
             Nmd = ifelse(is.na(Nmd), scaled_covs2$Nmd, Nmd),
             Nwtd = ifelse(is.na(Nwtd), scaled_covs2$Nwtd, Nwtd),
             Season = droplevels(Season))
    scaled_covs2 <- scaled_covs2 %>%
      mutate(PercForest = ifelse(is.na(PercForest), scaled_covs1$PercForest, PercForest), 
             Elev = ifelse(is.na(Elev), scaled_covs1$Elev, Elev),
             SppRich = ifelse(is.na(SppRich), scaled_covs1$SppRich, SppRich),
             SppDiversity = ifelse(is.na(SppDiversity), scaled_covs1$SppDiversity, SppDiversity),
             DomPrey = ifelse(is.na(DomPrey), scaled_covs1$DomPrey, DomPrey),
             Nelk = ifelse(is.na(Nelk), scaled_covs1$Nelk, Nelk),
             Nlagomorph = ifelse(is.na(Nlagomorph), scaled_covs1$Nlagomorph, Nlagomorph),
             Nlivestock = ifelse(is.na(Nlivestock), scaled_covs1$Nlivestock, Nlivestock),
             Nmoose = ifelse(is.na(Nmoose), scaled_covs1$Nmoose, Nmoose),
             Nmd = ifelse(is.na(Nmd), scaled_covs1$Nmd, Nmd),
             Nwtd = ifelse(is.na(Nwtd), scaled_covs1$Nwtd, Nwtd),
             Season = droplevels(Season))
    
    #'  List scaled data
    scaled_covs_list <- list(scaled_covs1, scaled_covs2)
    
    return(scaled_covs_list)
  }
  stations_eoe20s21s <- format_covs(cams_yr1 = cams_eoe20s, cams_yr2 = cams_eoe21s, 
                                    covs_yr1 = eoe_covs_20s, covs_yr2 = eoe_covs_21s)
  
  #'  Correct the mismatch - remove cameras with NO detection data either year
  #'  but were included in covariate data sets
  eh1 <- as.data.frame(rownames(dh))
  eh1$occ <- dh[,1,1]
  eh1 <- mutate(eh1, occ = ifelse(is.na(occ), 2, occ))
  names(eh1) <- c("NewLocationID", "det")
  covs1 <- stations_eoe20s21s[[1]][,1:3]
  tst1 <- full_join(eh1, covs1, by = "NewLocationID")
  rmv1 <- filter(tst1, is.na(det)) %>% dplyr::select(NewLocationID)
  
  eh2 <- as.data.frame(rownames(dh))
  eh2$occ <- dh[,1,2]
  eh2 <- mutate(eh2, occ = ifelse(is.na(occ), 2, occ))
  names(eh2) <- c("NewLocationID", "det")
  covs2 <- stations_eoe20s21s[[2]][,1:3]
  tst2 <- full_join(eh2, covs2, by = "NewLocationID")
  rmv2 <- filter(tst2, is.na(det)) %>% dplyr::select(NewLocationID)
  
  #'  Thinned covariate data
  stations_eoe20s <- stations_eoe20s21s[[1]][!(stations_eoe20s21s[[1]]$NewLocationID %in% rmv1$NewLocationID),]
  stations_eoe21s <- stations_eoe20s21s[[2]][!(stations_eoe20s21s[[2]]$NewLocationID %in% rmv2$NewLocationID),]
  stations_eoe20s21s <- list(stations_eoe20s, stations_eoe21s)
  
  #'  Save
  save(stations_eoe20s21s, file = "./Data/Covariates_extracted/stations_eoe20s21s_autolog_mod.RData")
  
  #'  Double check things are ordered correctly
  nrow(stations_eoe20s21s[[1]]); nrow(dh)
  stations_eoe20s21s[[1]][82:90,1:4]
  dh[82:90,1:3,1]
  stations_eoe20s21s[[2]][700:708,1:4]
  dh[700:708,1:3,2]
  
  #'  Correlation matrix to check for collinearity among continuous variables
  #'  Note: the species-specific total mortality and area-weighted mortality are 
  #'  highly correlated (1 or -1) so ignore those coefficients
  #'  Warnings are due to variables with no variation (Lion and Wolf mort km2 data)
  corr_matrix <- function(dat, firstcol, lastcol) {
    continuous_variables <- dat[,firstcol:lastcol]
    corr_all <- cor(continuous_variables)
    corr_all <- as.data.frame(round(corr_all, 2))
    print(corr_all)
    return(corr_all)
  }
  cov_corr_matrix <- corr_matrix(stations_eoe20s21s[[1]], firstcol = 6, lastcol = 17)
  cov_corr_matrix <- corr_matrix(stations_eoe20s21s[[2]], firstcol = 6, lastcol = 17)
  
  
  #'  ---------------------------
  ####  Survey-level covariates  ####
  #'  ---------------------------
  #'  Generate sampling effort array using "all_detections" function from above
  effort <- all_detections(DH_eoe20s_predators[[1]][[2]], DH_eoe21s_predators[[1]][[2]])
  effort_stacked <- rbind(effort[,,1], effort[,,2])
  
  #'  Scale survey-level covariates
  #'  Keep in mind NAs are not being filled in here and have no influence on scaling
  scale_srvy_cov <- function(time_covs) {
    #'  Find mean & standard deviation of covariates across all sites & occasions
    mu <- mean(as.matrix(time_covs), na.rm = TRUE)
    sd <- sd(as.matrix(time_covs), na.rm = TRUE)
    
    #'  Z-transform (center observations around mean & scale by 1 SD)
    scaled <- ((time_covs - mu) / sd)
    scaled <- round(scaled, 3)
    
    return(scaled)
  }
  effort_eoe20s21s <- scale_srvy_cov(effort_stacked)
  #'  Convert back to matrices
  effort_eoe20s <- as.matrix(effort_eoe20s21s[1:727,]); effort_eoe21s <- as.matrix(effort_eoe20s21s[728:1454,])
  
  #'  Convert annual detection histories into a single array
  effort_eoe20s21s <- array(data = c(effort_eoe20s, effort_eoe21s), 
                            dim = c(727, 11, 2), dimnames = list(rownames(effort_eoe20s), 
                                                                 colnames(effort_eoe20s), 
                                                                 c("year1", "year2")))
  
  #'  -----------------------
  ####  Covariate mean & SD  ####
  #'  -----------------------
  #'  Save unscaled covariates after specific rows have been excluded (needed for plotting later on)
  unscaled_covs <- function(cams_yr1, cams_yr2, covs_yr1, covs_yr2, rm_rows_yr1, rm_rows_yr2) {
    #'  Join camera data with extracted covariate data
    cam_covs1 <- full_join(cams_yr1, covs_yr1) %>% arrange(NewLocationID)
    cam_covs2 <- full_join(cams_yr2, covs_yr2) %>% arrange(NewLocationID)
    #'  Fill in NAs for cameras that sampled one but not both years
    #'  Assuming prey activity was similar across years at any given site
    cam_covs1 <- cam_covs1 %>%
      mutate(perc_forest = ifelse(is.na(perc_forest), cam_covs2$perc_forest, perc_forest),
             Elevation__10m2 = ifelse(is.na(Elevation__10m2), cam_covs2$Elevation__10m2, Elevation__10m2),
             SR = ifelse(is.na(SR), cam_covs2$SR, SR),
             H = ifelse(is.na(H), cam_covs2$H, H),
             dominantprey = ifelse(is.na(dominantprey), cam_covs2$dominantprey, dominantprey),
             elk_perday = ifelse(is.na(elk_perday), cam_covs2$elk_perday, elk_perday),
             lagomorphs_perday = ifelse(is.na(lagomorphs_perday), cam_covs2$lagomorphs_perday, lagomorphs_perday),
             livestock_perday = ifelse(is.na(livestock_perday), cam_covs2$livestock_perday, livestock_perday),
             moose_perday = ifelse(is.na(moose_perday), cam_covs2$moose_perday, moose_perday),
             muledeer_perday = ifelse(is.na(muledeer_perday), cam_covs2$muledeer_perday, muledeer_perday),
             whitetaileddeer_perday = ifelse(is.na(whitetaileddeer_perday), cam_covs2$whitetaileddeer_perday, whitetaileddeer_perday))
    cam_covs2 <- cam_covs2 %>%
      mutate(perc_forest = ifelse(is.na(perc_forest), cam_covs1$perc_forest, perc_forest),
             Elevation__10m2 = ifelse(is.na(Elevation__10m2), cam_covs1$Elevation__10m2, Elevation__10m2),
             SR = ifelse(is.na(SR), cam_covs1$SR, SR),
             H = ifelse(is.na(H), cam_covs1$H, H),
             dominantprey = ifelse(is.na(dominantprey), cam_covs1$dominantprey, dominantprey),
             elk_perday = ifelse(is.na(elk_perday), cam_covs1$elk_perday, elk_perday),
             lagomorphs_perday = ifelse(is.na(lagomorphs_perday), cam_covs1$lagomorphs_perday, lagomorphs_perday),
             livestock_perday = ifelse(is.na(livestock_perday), cam_covs1$livestock_perday, livestock_perday),
             moose_perday = ifelse(is.na(moose_perday), cam_covs1$moose_perday, moose_perday),
             muledeer_perday = ifelse(is.na(muledeer_perday), cam_covs1$muledeer_perday, muledeer_perday),
             whitetaileddeer_perday = ifelse(is.na(whitetaileddeer_perday), cam_covs1$whitetaileddeer_perday, whitetaileddeer_perday))
    #'  Remove cameras that were never operable but have covariate data
    cam_covs1 <- cam_covs1[!(cam_covs1$NewLocationID %in% rm_rows_yr1$NewLocationID),]
    cam_covs2 <- cam_covs2[!(cam_covs2$NewLocationID %in% rm_rows_yr2$NewLocationID),]
    #'  List un-scaled data
    unscaled_covs_list <- list(cam_covs1, cam_covs2)
    
    return(unscaled_covs_list)
  }
  stations_skinny_eoe20s21s <- unscaled_covs(cams_yr1 = cams_eoe20s, cams_yr2 = cams_eoe21s, 
                                             covs_yr1 = eoe_covs_20s, covs_yr2 = eoe_covs_21s,
                                             rm_rows_yr1 = rmv1, rm_rows_yr2 = rmv2)
  # save(stations_skinny_eoe20s21s, file = "./Data/Covariates_extracted/Covariate_skinny_EoE20s21s_autologic_mod.RData")
  
  #' #'  Save image of entire environment so it can be used with HPC
  #' save.image(file = "./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_autologic_img.RData")
  
  #'  Fin!
  #'  This is all sourced by Multi-Spp_OccMod_2spp_Bayesian.R
  
