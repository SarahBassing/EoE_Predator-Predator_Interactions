  #'  -------------------------------------------------
  #'  Format Data for 2-spp Occupancy Models (for JAGS) 
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  -------------------------------------------------
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
  dh1 <- DH_eoe20s_predators[[1]][[1]]
  dh2 <- DH_eoe21s_predators[[1]][[1]]
  newcols <- c("occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11")
  colnames(dh1) <- newcols; colnames(dh2) <- newcols
  dh <- list(dh1, dh2) 
  
  #'  Load camera station and covariate data
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Target", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(NewLocationID != "GMU6_U_160" | CameraHeight_M != 1.2) #' Remove duplicate camera where height changed slightly
  
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
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
  cams_eoe20w <- format_cam_station(eoe_probcams_20w, season = "Wtr20")
  cams_eoe21s <- format_cam_station(eoe_probcams_21s, season = "Smr21") 
  
  
  #'  Load extracted covariate data
  load("./Data/Covariates_extracted/Covariates_EoE_Smr20.RData")
  load("./Data/Covariates_extracted/Covariates_EoE_Wtr20.RData")
  load("./Data/Covariates_extracted/Covariates_EoE_Smr21.RData")
  # source("./Scripts/Data_Formatting/Covariate_Extract.R")
  
  
  #'  Scale and format site-level covariates
  format_covs <- function(cams, covs, rm_rows) {
    #'  Join camera data with extracted covariate data
    cam_covs <- full_join(cams, covs) %>% arrange(NewLocationID)
    #'  Remove rows where camera was inoperable the entire season - covariates at 
    #'  these sites shouldn't be included when scaling since they don't contribute
    #'  to detection data
    cam_covs <- cam_covs[-rm_rows,]
    #'  Rename, format, and scale as needed
    formatted <- transmute(cam_covs,
                           NewLocationID = as.factor(NewLocationID),
                           Season = as.factor(Season),
                           GMU = as.factor(GMU),
                           CameraFacing = as.factor(CameraFacing),
                           Setup = as.factor(Setup),
                           Target = as.factor(Target),
                           Habitat = as.factor(HabLayer_30m2),
                           DomPrey = as.factor(dominantprey),
                           Height = scale(CameraHeight_M),
                           PercForest = scale(perc_forest), 
                           Elev = scale(Elevation__10m2),
                           Dist2Burbs = scale(Dist2Suburbs),
                           logDist2Burbs = scale(log(Dist2Suburbs+1)),
                           Dist2Rrl = scale(Dist2Rural),
                           NearestRd = scale(dist2rd),
                           logNearestRd = scale(log(dist2rd+1)),
                           MinGroupSize = scale(avg_min_group_size), 
                           SppRich = scale(SR),
                           SppDiversity = scale(H),
                           Nelk = scale(elk_perday),    
                           Nhuman = scale(human_perday),
                           Nmotorized = scale(human_motorized_perday),
                           Nlagomorph = scale(lagomorphs_perday),
                           Nlivestock = scale(livestock_perday),
                           Nmoose = scale(moose_perday),
                           Nmd = scale(muledeer_perday),
                           Nwtd = scale(whitetaileddeer_perday),
                           Nungulate = scale(ungulate_perday),
                           Nbig_deer = scale(big_deer_perday),
                           Nsmall_deer = scale(small_deer_perday),
                           Bear_mort_n = scale(Bear_mort_n), 
                           Bear_mort_km2 = scale(Bear_mort_km2),
                           Bob_mort_n = scale(Bob_mort_n), 
                           Bob_mort_km2 = scale(Bob_mort_km2),
                           Lion_mort_n = scale(Lion_mort_n), 
                           Lion_mort_km2 = scale(Lion_mort_km2),
                           Wolf_mort_n = scale(Wolf_mort_n), 
                           Wolf_mort_km2 = scale(Wolf_mort_km2)) 
    
    #'  Adjust reference category for CameraFacing factors
    order_camfacing <- c("random", "trail", "road")
    formatted <- formatted %>%
      mutate(
        CameraFacing = fct_relevel(CameraFacing, order_camfacing)
      )
    
    #'  Identify which sites have missing data for distance & height covariates
    ung_cams <- formatted[formatted$Setup == "ungualte",]
    pred_cams <- formatted[formatted$Setup == "predator",]
    noHgt_ung <- which(is.na(ung_cams$Height)); print(noHgt_ung)
    noHgt_pred <- which(is.na(pred_cams$Height)); print(noHgt_pred)
    
    return(formatted)
  }
  rm_rows_eoe20s <- c(61, 79, 82, 98, 125, 157, 171, 177, 178, 181, 186, 192, 200, 214, 228, 235, 236, 259, 311, 334, 346, 361, 371, 379, 380, 385, 433, 437, 439, 458, 493)
  rm_rows_eoe21s <- c(6, 106, 112, 116, 127, 145, 147, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682)
  stations_eoe20s <- format_covs(cams = cams_eoe20s, covs = eoe_covs_20s, rm_rows = rm_rows_eoe20s) %>%
    #'  Force values to 0 (mean when scaled) since scale() doesn't work when no variation in data
    mutate(Lion_mort_km2 = 0)
  stations_eoe21s <- format_covs(cams = cams_eoe21s, covs = eoe_covs_21s, rm_rows = rm_rows_eoe21s) %>%
    #'  Force values to 0 (mean when scaled) since scale() doesn't work when no variation in data
    mutate(Lion_mort_km2 = 0)
  stations_eoe20s21s <- list(stations_eoe20s, stations_eoe21s)
  
  #' #'  Save
  #' save(stations_eoe20s21s, file = "./Data/Covariates_extracted/stations_eoe20s21s.RData")

  #'  Double check things are ordered correctly
  nrow(stations_eoe20s21s[[1]]); nrow(dh[[1]])
  stations_eoe20s21s[[1]][82:90,1:4]
  dh[[1]][82:90,1:3]
  stations_eoe20s21s[[2]][600:608,1:4]
  dh[[2]][600:608,1:3]
  
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
  cov_corr_matrix <- corr_matrix(stations_eoe20s21s[[1]], firstcol = 9, lastcol = 38)
  cov_corr_matrix <- corr_matrix(stations_eoe20s21s[[2]], firstcol = 9, lastcol = 38)
  
  
  #'  ---------------------------
  ####  Survey-level covariates  ####
  #'  ---------------------------
  #'  Load survey-level data
  load("Data/Wolf count data/count_eoe20s_wolf.RData")
  load("Data/Wolf count data/count_eoe21s_wolf.RData")
  
  count_eoe20s21s_wolf <- list(count_eoe20s_wolf[[1]], count_eoe21s_wolf[[1]])
  
  #'  Replace NAs in sampling effort with 0 - these sites truly were not surveyed
  #'  during those sampling occasions so survey really is 0
  count_eoe20s_effort <- replace(count_eoe20s_wolf[[2]], is.na(count_eoe20s_wolf[[2]]), 0)
  count_eoe21s_effort <- replace(count_eoe21s_wolf[[2]], is.na(count_eoe21s_wolf[[2]]), 0)
  count_eoe20s21s_effort <- list(count_eoe20s_effort, count_eoe21s_effort)
  nrow(count_eoe20s21s_effort[[1]]);   nrow(count_eoe20s21s_effort[[2]])
  
  #'  Scale survey-level covariates
  scale_srvy_cov <- function(time_covs) {
    #'  Find mean & standard deviation of covariates across all sites & occasions
    mu <- mean(as.matrix(time_covs), na.rm = TRUE)
    sd <- sd(as.matrix(time_covs), na.rm = TRUE)
    
    #'  Z-transform (center observations around mean & scale by 1 SD)
    scaled <- ((time_covs - mu) / sd)
    scaled <- round(scaled, 3)
    
    return(scaled)
  }
  #'  Note: rm_rows have already been removed when generated in Detection_histories_for_occmod.R
  wolf_activity_eoe20s21s <- lapply(count_eoe20s21s_wolf, scale_srvy_cov)
  effort_eoe20s21s <- lapply(count_eoe20s21s_effort, scale_srvy_cov)
  
  #'  Create list of survey level covariates for unmarked 
  srvy_covs_eoe20s <- list(wolf_activity = wolf_activity_eoe20s21s[[1]], effort = effort_eoe20s21s[[1]])
  srvy_covs_eoe21s <- list(wolf_activity = wolf_activity_eoe20s21s[[2]], effort = effort_eoe20s21s[[2]])
  
  
  #'  -----------------------
  ####  Covariate mean & SD  ####
  #'  -----------------------
  #'  Save unscaled covariates after specific rows have been excluded (needed for plotting later on)
  unscaled_covs <- function(cams, covs, rm_rows) {
    #'  Join camera data with extracted covariate data
    cam_covs <- full_join(cams, covs) %>% arrange(NewLocationID)
    #'  Remove rows where camera was inoperable the entire season - covariates at 
    #'  these sites shouldn't be included when scaling since they don't contribute
    #'  to detection data
    cam_covs <- cam_covs[-rm_rows,]
    cam_covs <- transmute(cam_covs,
                          NewLocationID = as.factor(NewLocationID),
                          Season = as.factor(Season),
                          GMU = as.factor(GMU),
                          CameraFacing = as.factor(CameraFacing),
                          Setup = as.factor(Setup),
                          Target = as.factor(Target),
                          DomPrey = as.factor(dominantprey),
                          Habitat = as.factor(HabLayer_30m2),
                          Height = CameraHeight_M,
                          PercForest = perc_forest, 
                          Elev = Elevation__10m2,
                          Dist2Burbs = Dist2Suburbs,
                          logDist2Burbs = log(Dist2Suburbs+1),
                          Dist2Rrl = Dist2Rural,
                          NearestRd = dist2rd,
                          logNearestRd = log(dist2rd+1),
                          MinGroupSize = avg_min_group_size, 
                          SppRich = SR,
                          SppDiversity = H,
                          Nelk = elk_perday,    
                          Nhuman = human_perday,
                          Nmotorized = human_motorized_perday,
                          Nlagomorph = lagomorphs_perday,
                          Nlivestock = livestock_perday,
                          Nmoose = moose_perday,
                          Nmd = muledeer_perday,
                          Nwtd = whitetaileddeer_perday,
                          Nungulate = ungulate_perday,
                          Nbig_deer = big_deer_perday,
                          Nsmall_deer = small_deer_perday,
                          Bear_mort_n = Bear_mort_n, 
                          Bear_mort_km2 = Bear_mort_km2,
                          Bob_mort_n = Bob_mort_n, 
                          Bob_mort_km2 = Bob_mort_km2,
                          Lion_mort_n = Lion_mort_n, 
                          Lion_mort_km2 = Lion_mort_km2,
                          Wolf_mort_n = Wolf_mort_n, 
                          Wolf_mort_km2 = Wolf_mort_km2) %>%
      return(cam_covs)
  }
  stations_skinny_eoe20s <- unscaled_covs(cams = cams_eoe20s, covs = eoe_covs_20s, rm_rows = rm_rows_eoe20s)
  stations_skinny_eoe21s <- unscaled_covs(cams = cams_eoe21s, covs = eoe_covs_21s, rm_rows = rm_rows_eoe21s)
  stations_skinny_eoe20s21s <- list(stations_skinny_eoe20s, stations_skinny_eoe21s)
  # save(stations_skinny_eoe20s21s, file = "./Data/Covariates_extracted/Covariate_skinny_EoE20s21s_sepYr.RData")
  
  stations_skinny_eoe20s <- unscaled_covs(cams_eoe20s, covs = eoe_covs_20s, rm_rows = rm_rows_eoe20s)
  stations_skinny_eoe21s <- unscaled_covs(cams_eoe21s, covs = eoe_covs_21s, rm_rows = rm_rows_eoe21s)
  
  # save(stations_skinny_eoe20s, file = "./Data/Covariates_extracted/Covariate_skinny_EoE20s_sepYr.RData")
  # save(stations_skinny_eoe21s, file = "./Data/Covariates_extracted/Covariate_skinny_EoE21s_sepYr.RData")
  
  
  #' #'  Save image of entire environment so it can be used with HPC
  #' save.image(file = "./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_sepYr_img.RData")
  
  #'  Fin!
  #'  This is all sourced by Multi-Spp_OccMod_2spp_Bayesian.R
  
  
  
  #' #'  Double check DH and covs line up - same number of observations?!
  #' #'  Function to combine detection and covariate data in an unmarked data frame
  #' #'  Not needed for JAGS but good way to double check these different datasets
  #' #'  are consitent
  #' umf_setup <- function(dh_spp1, dh_spp2, listnames, sitecovs, srvycovs, plotit = T) {
  #' #'  List detection histories of interacting species
  #' spp12_DH <- list(spp1 = dh_spp1, spp2 = dh_spp2)
  #' #'  Rename lists by species
  #' names(spp12_DH) <- listnames
  #' #'  Create array with detection histories, site-level and survey-level covariates
  #' spp12_UMF <- unmarkedFrameOccuMulti(y = spp12_DH,
  #'                                     siteCovs = sitecovs,
  #'                                     obsCovs = srvycovs,
  #'                                     maxOrder = 2)
  #' #'  Plot detection histories
  #' print(plot(spp12_UMF))
  #' #'  Summarize data
  #' summary(spp12_UMF)
  #' return(spp12_UMF)
  #' }
  #' testing_umf <- umf_setup(dh_spp1 = DH_eoe20s_predators[[5]][[1]], dh_spp2 = DH_eoe20s_predators[[1]][[1]], 
  #'                                listnames = c("wolf", "bear"), sitecovs = stations_eoe20s, srvycovs = srvy_covs_eoe20s)
  
  