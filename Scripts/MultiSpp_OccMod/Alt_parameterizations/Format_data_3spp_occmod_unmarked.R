  #'  --------------------------------------------
  #'  Format Data for 3-spp Occupancy Models (unmarked) 
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  --------------------------------------------
  #'  Script to create unmarked dataframes and format covariate data for multi-
  #'  species occupancy models with EoE data. Focusing on black bear, bobcat,
  #'  coyote, mountain lion, and wolf interactions. Cameras ran summer 2020, winter
  #'  2020-2021, and summer 2021. Co-occurrence models test whether predator
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
  
  #'  Load detection histories
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20w_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")
  
  #'  Load camera station and covariate data
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Target", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(NewLocationID != "GMU6_U_160" | CameraHeight_M != 1.2) #' Remove duplicate camera where height changed slightly
  
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
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
    cam_covs <- full_join(cams, covs)
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
                           Height = scale(CameraHeight_M),
                           PercForest = scale(perc_forest), 
                           Elev = scale(Elevation__10m2),
                           Dist2Burbs = scale(Dist2Suburbs),
                           logDist2Burbs = scale(log(Dist2Suburbs+1)),
                           Dist2Rrl = scale(Dist2Rural),
                           NearestRd = scale(dist2rd),
                           logNearestRd = scale(log(dist2rd+1)),
                           MinGroupSize = scale(avg_min_group_size), 
                           Nelk = scale(elk),
                           Nhuman = scale(human),
                           Nhumanplus = scale(human_plus),
                           Nmotorized = scale(human_motorized),
                           Nlagomorph = scale(lagomorphs),
                           Nlivestock = scale(livestock),
                           Nmoose = scale(moose),
                           Nmd = scale(muledeer),
                           Nwtd = scale(whitetaileddeer),
                           Nungulate = scale(ungulate),
                           Bear_mort_n = scale(Bear_mort_n), 
                           Bear_mort_km2 = scale(Bear_mort_km2),
                           Bob_mort_n = scale(Bob_mort_n), 
                           Bob_mort_km2 = scale(Bob_mort_km2),
                           Lion_mort_n = scale(Lion_mort_n), 
                           Lion_mort_km2 = scale(Lion_mort_km2),
                           Wolf_mort_n = scale(Wolf_mort_n), 
                           Wolf_mort_km2 = scale(Wolf_mort_km2)
    ) %>%
      arrange(NewLocationID) #NECESSARY TO MATCH DH's CAMERALOCATION ORDER
    
    #'  Adjust reference category for CameraFacing factors
    order_camfacing <- c("random", "trail", "road")
    formatted <- formatted %>%
      mutate(
        CameraFacing = fct_relevel(CameraFacing, order_camfacing)
      )
    
    #'  Identify which sites have missing data for distance & height covariates
    ung_cams <- formatted[formatted$Setup == "ungualte",]
    pred_cams <- formatted[formatted$Setup == "predator",]
    # noDist <- which(is.na(formatted$Distance)); print(noDist)
    noHgt_ung <- which(is.na(ung_cams$Height)); print(noHgt_ung)
    noHgt_pred <- which(is.na(pred_cams$Height)); print(noHgt_pred)
    
    #'  Replace missing values with mean once covariates are z-transformed
    # formatted$Distance[is.na(formatted$Distance),] <- 0
    # formatted$Height[is.na(formatted$Height),] <- 0   
    
    return(formatted)
  }
  rm_rows_eoe20s <- c(61, 79, 82, 98, 125, 157, 171, 177, 178, 181, 186, 192, 200, 214, 228, 235, 236, 259, 311, 334, 346, 361, 371, 379, 380, 385, 433, 437, 439, 458, 493)
  stations_eoe20s <- format_covs(cams_eoe20s, covs = eoe_covs_20s, rm_rows = rm_rows_eoe20s) %>%
    #'  Force values to 0 (mean when scaled) since scale() doesn't work when no variation in data
    mutate(Lion_mort_km2 = 0)
  
  rm_rows_eoe20w <- c(7, 16, 32, 43, 123, 138, 195, 215, 227, 242, 252, 268)
  stations_eoe20w <- format_covs(cams_eoe20w, covs = eoe_covs_20w, rm_rows = rm_rows_eoe20w) %>%
    mutate(Wolf_mort_km2 = 0)
  
  rm_rows_eoe21s <- c(6, 106, 112, 116, 127, 145, 147, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682)
  stations_eoe21s <- format_covs(cams_eoe21s, covs = eoe_covs_21s, rm_rows = rm_rows_eoe21s) %>%
    mutate(Lion_mort_km2 = 0)
  
  #'  Double check things are ordered correctly
  stations_eoe21s[82:90,1:4]
  DH_eoe21s_predators[[1]][[1]][82:90,1:3]

  
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
  camera_station_list <- list(stations_eoe20s, stations_eoe20w, stations_eoe21s)
  cov_corr_matrix <- lapply(camera_station_list, corr_matrix, firstcol = 8, lastcol = 36)
  #'  Dist2Suburbs & Dist2Rural moderately correlated (r ranged 0.58 - 0.64 depending on year)
  #'  Dist2Suburbs & mortality correlated for most species and years (r values included 0.62, -0.87, etc.)
  #'  Nhuman & Nhumanplus highly correlate (duh), and Nmotorized correlated with both in winter (0.59 - 0.6)
  
  #'  ---------------------------
  ####  Survey-level covariates  ####
  #'  ---------------------------
  #'  Load survey-level data
  load("Data/Wolf count data/count_eoe20s_wolf.RData")
  load("Data/Wolf count data/count_eoe20w_wolf.RData")
  load("Data/Wolf count data/count_eoe21s_wolf.RData")
  
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
  wolf_activity_eoe20s <- scale_srvy_cov(count_eoe20s_wolf[[1]])
  wolf_activity_eoe20w <- scale_srvy_cov(count_eoe20w_wolf[[1]])
  wolf_activity_eoe21s <- scale_srvy_cov(count_eoe21s_wolf[[1]])
  
  #'  Create list of survey level covariates for unmarked 
  srvy_covs_eoe20s <- list(wolf_activity = wolf_activity_eoe20s)
  srvy_covs_eoe20w <- list(wolf_activity = wolf_activity_eoe20w)
  srvy_covs_eoe21s <- list(wolf_activity = wolf_activity_eoe21s)
  
  
  #'  ---------------------------
  ####  Setup data for unmarked  ####
  #'  ---------------------------
  #'  Multi-species unmarkedDF --> unmarkedFrameOccuMulti (pg 151 of unmarked manual)
  #'  List relevant detection histories and create unmarked data frame for each 
  #'  multispecies occupancy model. Currently running 2-species occupancy models.
  #'  Define maximum interaction order with maxOrder. Defaults to all possible
  #'  interactions if not defined. 
  #'  Detection histories have to already be generated for this to work. 
  #'  List order in summer:
  #'  Black bear [1], Bobcat [2], Coyote [3], Mountain lion [4], Wolf [5],
  #'  Detection history [1] and sampling effort [2] listed per species
  #'  List order in winter:
  #'  Bobcat [1], Coyote [2], Mountain lion [3], Wolf [4]
  #'  ---------------------------
  #'  Function to create unmarked dataframes for multispecies occupancy models 
  #'  for each pairwise combination and season
  umf_setup_3wayX <- function(dh_spp1, dh_spp2, dh_spp3, listnames, sitecovs, srvycovs, plotit = T) {
    
    #'  List detection histories of interacting species
    spp123_DH <- list(spp1 = dh_spp1, spp2 = dh_spp2, spp3 = dh_spp3)
    #'  Rename lists by species
    names(spp123_DH) <- listnames
    
    #'  Create array with detection histories, site-level and survey-level covariates
    spp123_UMF <- unmarkedFrameOccuMulti(y = spp123_DH,
                                        siteCovs = sitecovs,
                                        obsCovs = srvycovs,
                                        maxOrder = 2)
    
    #'  Plot detection histories
    print(plot(spp123_UMF))
    #'  Summarize data
    summary(spp123_UMF)
    
    return(spp123_UMF)
  }
  
  ####  Wolf-Bear-Lion UMFs  ####
  wolf_bear_lion_20s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20s_predators[[5]][[1]], dh_spp2 = DH_eoe20s_predators[[1]][[1]], 
                                            dh_spp3 = DH_eoe20s_predators[[4]][[1]], listnames = c("wolf", "bear", "lion"), 
                                            sitecovs = stations_eoe20s, srvycovs = srvy_covs_eoe20s)
  wolf_bear_lion_21s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe21s_predators[[5]][[1]], dh_spp2 = DH_eoe21s_predators[[1]][[1]], 
                                            dh_spp3 = DH_eoe21s_predators[[4]][[1]], listnames = c("wolf", "bear", "lion"), 
                                            sitecovs = stations_eoe21s, srvycovs = srvy_covs_eoe21s)
  
  
  ####  Wolf-Bobcat-Lion UMFs  ####
  wolf_bob_lion_20s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20s_predators[[5]][[1]], dh_spp2 = DH_eoe20s_predators[[2]][[1]], 
                                            dh_spp3 = DH_eoe20s_predators[[4]][[1]], listnames = c("wolf", "bobcat", "lion"), 
                                            sitecovs = stations_eoe20s, srvycovs = srvy_covs_eoe20s)
  wolf_bob_lion_20w_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20w_predators[[4]][[1]], dh_spp2 = DH_eoe20w_predators[[1]][[1]], 
                                           dh_spp3 = DH_eoe20w_predators[[3]][[1]], listnames = c("wolf", "bobcat", "lion"), 
                                           sitecovs = stations_eoe20w, srvycovs = srvy_covs_eoe20w)
  wolf_bob_lion_21s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe21s_predators[[5]][[1]], dh_spp2 = DH_eoe21s_predators[[2]][[1]], 
                                            dh_spp3 = DH_eoe21s_predators[[4]][[1]], listnames = c("wolf", "bobcat", "lion"), 
                                            sitecovs = stations_eoe21s, srvycovs = srvy_covs_eoe21s)
  
  
  ####  Wolf-Coyote-Lion UMFs  ####
  wolf_coy_lion_20s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20s_predators[[5]][[1]], dh_spp2 = DH_eoe20s_predators[[3]][[1]], 
                                           dh_spp3 = DH_eoe20s_predators[[4]][[1]], listnames = c("wolf", "coyote", "lion"), 
                                           sitecovs = stations_eoe20s, srvycovs = srvy_covs_eoe20s)
  wolf_coy_lion_20w_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20w_predators[[4]][[1]], dh_spp2 = DH_eoe20w_predators[[2]][[1]], 
                                           dh_spp3 = DH_eoe20w_predators[[3]][[1]], listnames = c("wolf", "coyote", "lion"), 
                                           sitecovs = stations_eoe20w, srvycovs = srvy_covs_eoe20w)
  wolf_coy_lion_21s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe21s_predators[[5]][[1]], dh_spp2 = DH_eoe21s_predators[[3]][[1]], 
                                           dh_spp3 = DH_eoe21s_predators[[4]][[1]], listnames = c("wolf", "coyote", "lion"), 
                                           sitecovs = stations_eoe21s, srvycovs = srvy_covs_eoe21s)
  
  
  ####  Wolf-Bobcat-Coyote UMFs  ####
  wolf_bob_coy_20s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20s_predators[[5]][[1]], dh_spp2 = DH_eoe20s_predators[[2]][[1]], 
                                           dh_spp3 = DH_eoe20s_predators[[3]][[1]], listnames = c("wolf", "bobcat", "coyote"), 
                                           sitecovs = stations_eoe20s, srvycovs = srvy_covs_eoe20s)
  wolf_bob_coy_20w_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20w_predators[[4]][[1]], dh_spp2 = DH_eoe20w_predators[[1]][[1]], 
                                           dh_spp3 = DH_eoe20w_predators[[2]][[1]], listnames = c("wolf", "bobcat", "coyote"), 
                                           sitecovs = stations_eoe20w, srvycovs = srvy_covs_eoe20w)
  wolf_bob_coy_21s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe21s_predators[[5]][[1]], dh_spp2 = DH_eoe21s_predators[[2]][[1]], 
                                           dh_spp3 = DH_eoe21s_predators[[3]][[1]], listnames = c("wolf", "bobcat", "coyote"), 
                                           sitecovs = stations_eoe21s, srvycovs = srvy_covs_eoe21s)
  
  
  ####  Bear-bobcat-Coyote UMFs  ####
  bear_bob_coy_20s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe20s_predators[[1]][[1]], dh_spp2 = DH_eoe20s_predators[[2]][[1]], 
                                            dh_spp3 = DH_eoe20s_predators[[3]][[1]], listnames = c("bear", "bobcat", "coyote"), 
                                            sitecovs = stations_eoe20s, srvycovs = srvy_covs_eoe20s)
  bear_bob_coy_21s_umf <- umf_setup_3wayX(dh_spp1 = DH_eoe21s_predators[[1]][[1]], dh_spp2 = DH_eoe21s_predators[[2]][[1]], 
                                            dh_spp3 = DH_eoe21s_predators[[3]][[1]], listnames = c("bear", "bobcat", "coyote"), 
                                            sitecovs = stations_eoe21s, srvycovs = srvy_covs_eoe21s)
  