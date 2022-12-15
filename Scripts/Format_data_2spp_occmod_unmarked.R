  #'  --------------------------------------------
  #'  Format Data for 2-spp Occupancy Models (unmarked) 
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  --------------------------------------------
  #'  Script to create unmarked dataframes and format covariate data for multi-
  #'  species occupancy models with EoE data. Focusing on black bear, bobcat,
  #'  coyote, mountain lion, and wolf interactions. Cameras ran summer 2020, winter
  #'  2020-2021, and winter 2021. Co-occurrence models test whether predator
  #'  co-occurrence is non-independent and whether their occurrence, co-occurrence,
  #'  and detection are influenced by other variables of interest.
  #'  
  #'  Summer primary period is considered July 1 - Spet. 15, equating to 11 1-wk 
  #'  sampling periods. 
  #'  Winter primary period is considered ....
  #'   
  #'  
  #'  Encounter histories are generated with the Detection_histories_for_occmod.R
  #'  --------------------------------------------
  
  #'  Load libraries
  library(unmarked)
  library(stringr)
  library(tidyverse)
  
  #'  Load detection histories
  load("./Data/Detection_Histories/DH_eoe21s_predators.RData")
  
  
  #'  Load camera station and covariate data
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Target", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(Season == "Smr21")
  
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  cams_eoe21s <- dplyr::select(eoe_probcams_21s, c("NewLocationID", "Lat", "Long")) %>%
    mutate(GMU = str_extract(NewLocationID, "[^_]+"),
           Season = "Smr21") %>%
    left_join(cams_eoe_long, by = c("NewLocationID", "Season")) %>%
    #'  Drop handful of duplicated rows (not sure why this happens with left_join)
    unique() %>%
    mutate(CameraFacing = ifelse(CameraFacing == "road" | CameraFacing == "atv" | CameraFacing == "gravel" | CameraFacing == "decommision", "road", CameraFacing),
           CameraFacing = ifelse(CameraFacing == "hiking" | CameraFacing == "game" | CameraFacing == "other", "trail", CameraFacing), 
           CameraFacing = ifelse(is.na(CameraFacing), "trail", CameraFacing)) %>% # predator cameras so either trail or road
    arrange(NewLocationID)
  
  #'  Load extracted covariate data
  load("./Data/Covariates_EoE_Smr21.RData")
  # source("./Scripts/Covariate_Extract.R")
  
  #'  Scale and format site-level covariates
  format_covs <- function(cams, covs) {
    #'  Join camera data with extracted covariate data
    cam_covs <- full_join(cams, covs)
    #'  Rename, format, and scale as needed
    formatted <- transmute(cam_covs,
                           NewLocationID = as.factor(NewLocationID),
                           #Year = as.factor(Year),
                           Season = as.factor(Season),
                           GMU = as.factor(GMU),
                           #Distance = scale(Distance_Focal_Point),
                           Height = scale(CameraHeight_M),
                           CameraFacing = as.factor(CameraFacing),
                           Setup = as.factor(Setup),
                           Target = as.factor(Target),
                           perc_forest = scale(perc_forest), 
                           min_group_size = scale(avg_min_group_size)
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
    formatted$Height[is.na(formatted$Height),] <- 0   #########THIS NEEDS TO BE FIXED! WE HAVE THE HEIGHT FOR MOST OF THESE FROM PREVIOUS YEARS
    
    return(formatted)
  }
  stations_eoe21s <- format_covs(cams_eoe21s, covs = eoe_covs_21s)
  #'  Remove rows when camera was inoperable the entire season
  stations_eoe21s <- stations_eoe21s[-c(6, 106, 112, 116, 127, 145, 178, 194, 195, 260, 267, 296, 343, 355, 365, 409, 417, 419, 423, 430, 450, 510, 530, 577, 578, 580, 588, 621, 627, 647, 652, 682),]
  
  #'  Double check things are ordered correctly
  stations_eoe21s[82:90,]
  DH_eoe21s_predators[[1]][[1]][82:90,1:3]
  
  
  ####  EVENTUALLY CHECK FOR COLLINEARITY WITH MORE CONTINUOUS VARIABLES  ####
  cor(stations_eoe21s$perc_forest, stations_eoe21s$min_group_size, use = "complete.obs")
  
  
  #'  ---------------------------
  ####  Survey-level covariates  ####
  #'  ---------------------------
  #'  Load survey-level data
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
  wolf_activity <- scale_srvy_cov(count_eoe21s_wolf[[1]])
  
  #'  Create list of survey level covariates
  srvy_covs <- list(wolf_activity = wolf_activity)
  
  #'  ---------------------------
  ####  Setup data for unmarked  ####
  #'  ---------------------------
  #'  Multi-species unmarkedDF --> unmarkedFrameOccuMulti (pg 151 of unmarked manual)
  #'  List relevant detection histories and create unmarked data frame for each 
  #'  multi-species occupancy model. Currently running 3-species occupancy models.
  #'  Define maximum interaction order with maxOrder. Defaults to all possible
  #'  interactions if not defined. 
  #'  Detection histories have to already be generated for this to work.
  #'  ---------------------------
  
  ####  BEAR-LION UMF  ####
  bear_lion_smr21_DH <- list(bear = DH_eoe21s_predators[[1]][[1]], lion = DH_eoe21s_predators[[4]][[1]])
  bear_lion_smr21_UMF <- unmarkedFrameOccuMulti(y = bear_lion_smr21_DH,
                                           siteCovs = stations_eoe21s,
                                           obsCovs = srvy_covs,
                                           maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_lion_smr21_UMF)
  #'  Review covariates
  summary(bear_lion_smr21_UMF)
  #'  Look at natural parameters f design matrix
  bear_lion_smr21_UMF@fDesign
  
  
  ####  COY-LION UMF  ####
  coy_lion_smr21_DH <- list(coy = DH_eoe21s_predators[[3]][[1]], lion = DH_eoe21s_predators[[4]][[1]])
  coy_lion_UMF <- unmarkedFrameOccuMulti(y = coy_lion_smr21_DH,
                                              siteCovs = stations_eoe21s,
                                              obsCovs = srvy_covs,
                                              maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_lion_UMF)
  #'  Review covariates
  summary(coy_lion_UMF)
  
  
  ####  BOB-LION UMF  ####
  bob_lion_smr21_DH <- list(bob = DH_eoe21s_predators[[2]][[1]], lion = DH_eoe21s_predators[[4]][[1]])
  bob_lion_UMF <- unmarkedFrameOccuMulti(y = bob_lion_smr21_DH,
                                              siteCovs = stations_eoe21s,
                                              obsCovs = srvy_covs,
                                              maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_lion_UMF)
  #'  Review covariates
  summary(bob_lion_UMF)
  
  
  
