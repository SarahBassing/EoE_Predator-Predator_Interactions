  #'  -------------------------------
  #'  Royle-Nichols Abundance Model
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  November 2023
  #'  -------------------------------
  #'  Script to generate detection histories and run RN models to estimate relative
  #'  abundance at each camera site.
  #'  
  #'  Camera operations table generated in Detection_data_cleaning.R
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(camtrapR)
  library(chron)
  library(lubridate)
  library(unmarked)
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  ------------------------------
  ####  Load and format image data  ####
  #'  ------------------------------
  #'  Load detection data
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_NewLocationID_2023-09-26.RData") 
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_NewLocationID_2023-09-26.RData") 
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID_2023-09-26.RData") 
  
  #'  Sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe22s_sequential_probimgs.RData")
  
  #'  Problem cameras
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe22s_problem_cams.RData")
  
  #'  Filter detection data to time period of interest and remove problem images
  thin_detections <- function(dets, seqprobs, start_date, end_date) {
    #'  Remove images with known problems (e.g., viewshed obscured)
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    
    #'  Filter images
    clean_dets <- skinny_dets %>%
      filter(Species != "none") %>%
      filter(Vehicle != "TRUE") %>%
      filter(OpState != "maintenance") %>%
      #'  Add count = 1 for species missing count data (mainly humans, rabbit_hare, cattle_cow)
      mutate(Count = ifelse(Count == 0, 1, Count),
             Count = ifelse(is.na(Count), 1, Count),
             Date = as.Date(Date, format = "%Y-%m-%d")) %>% 
      #'  Filter to images to desired date range
      filter(Date >= start_date & Date <= end_date) %>%
      #'  Remove single mislabeled image from 2022 data (will not affect detection histories)
      filter(NewLocationID != "GMU10A_U_50" | Date != "2022-08-27" | Species != "human") %>%
      #'  Remove observations that can't be linked to a camera with coordinates
      filter(!is.na(NewLocationID)) %>%
      arrange(NewLocationID, posix_date_time) %>%
      dplyr::select(c("NewLocationID", "posix_date_time", "Date", "Time", "Species", "Count", "OpState")) 
    
    return(clean_dets)
  }
  df_all_20s <- thin_detections(eoe20s_allM, seqprobs = eoe_seqprob_20s, start_date = "2020-07-01", end_date = "2020-09-15") 
  df_all_21s <- thin_detections(eoe21s_allM, seqprobs = eoe_seqprob_21s, start_date = "2021-07-01", end_date = "2021-09-15")
  df_all_22s <- thin_detections(eoe22s_allM, seqprobs = eoe_seqprob_22s, start_date = "2022-07-01", end_date = "2022-09-15") 

  #'  ------------------------------------------
  #####  Generate independent detection events  #####
  #'  ------------------------------------------
  #'  Create a column grouping images into "independent" detection events based 
  #'  on defined amount of time that must elapse between images of the same species 
  #'  at the same camera site. Filter to the first image of each event and then
  #'  split into species-specific detection events.
  #'  Currently using 5 min as interval (5*60 = 300 seconds)
  #'  -----------------------------------------
  unique_detections <- function(dets, elapsed_time) {
    det_events <- dets %>%
      arrange(NewLocationID, posix_date_time) %>%
      #'  Flag images of same species at same camera as being a different detection event
      #'  when time since last image of that species is greater than defined time interval
      group_by(NewLocationID, Species) %>%
      mutate(det_events = cumsum(c(1, diff(posix_date_time) > elapsed_time))) %>% # units in seconds! 
      ungroup() %>%
      #'  Retain only the first image from each unique detection event
      group_by(NewLocationID, Species, det_events) %>%
      slice(1L) %>%
      ungroup()
    
    return(det_events)
  }
  npp20s_det_events <- unique_detections(df_all_20s, elapsed_time = 300)
  npp21s_det_events <- unique_detections(df_all_21s, elapsed_time = 300)
  npp22s_det_events <- unique_detections(df_all_22s, elapsed_time = 300)
  
  #'  ---------------------------
  #####  Camera operation table  #####
  #'  --------------------------- 
  #'  Create a matrix with each camera & dates it deployed
  #'  1 = operating; 0 = not operating but deployed; NA = not deployed
  camera_operation_tbl <- function(cams) {
    cams <- arrange(cams, NewLocationID)
    #'  Add one day to retrieval date when setup & retrieval dates are the same
    #'  Necessary for cameraOperation function below... so annoying
    same_startend <- cams[cams$Setup_date == cams$Retrieval_date,]
    same_startend <- mutate(same_startend, Retrieval_date = as.Date(Retrieval_date) + 1)
    same_startend <- mutate(same_startend, Problem1_to = as.Date(Problem1_to) +1)
    cams[match(same_startend$NewLocationID, cams$NewLocationID),] <- same_startend
    
    #'  Make sure camera data are organized by NewLocationID
    cams <- arrange(cams, NewLocationID)
    
    #'  Create camera operation table
    camop_problem <- cameraOperation(CTtable = cams,
                                     stationCol = "NewLocationID",
                                     setupCol = "Setup_date",
                                     retrievalCol = "Retrieval_date",
                                     hasProblems = TRUE,
                                     dateFormat = "%Y-%m-%d", 
                                     writecsv = FALSE) 
    
    return(camop_problem)
  }
  npp20s_probs <- camera_operation_tbl(eoe_probcams_20s) 
  npp21s_probs <- camera_operation_tbl(eoe_probcams_21s)
  npp22s_probs <- camera_operation_tbl(eoe_probcams_22s)
  
  #'  ------------------------
  #####  Detection histories  #####
  #'  ------------------------
  #'  Create species-specific detection histories and sampling effort matrix
  #'  FYI: July 1 - Sept 15 = 11 7-day sampling periods and 77 1-day sampling periods
  DH <- function(dets, cam_probs, spp, start_date, y, rm_rows, oc) {
    det_hist <- detectionHistory(recordTable = dets,
                                 camOp = cam_probs,
                                 stationCol = "NewLocationID",
                                 # speciesCol = "Species",
                                 recordDateTimeCol = "posix_date_time",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 species = spp,
                                 occasionLength = 1,
                                 day1 = start_date, 
                                 datesAsOccasionNames = FALSE,
                                 # occasionStartTime = 12, # starts at noon
                                 timeZone = "America/Edmonton",
                                 output = y,
                                 includeEffort = TRUE,
                                 scaleEffort = FALSE,
                                 # writecsv = TRUE,
                                 outDir = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories")
    
    #'  Reduce detection histories to sampling occasions of interest (drop extra
    #'  occasions after focal period of interest)
    short_dh <- det_hist[[1]][,1:oc]
    short_effort <- det_hist[[2]][,1:oc]
    
    #'  Remove rows where camera was inoperable during the entire sampling season
    short_dh <- as.data.frame(short_dh)
    short_effort <- as.data.frame(short_effort)
    short_dh <- filter(short_dh, rowSums(is.na(short_dh)) != ncol(short_dh))
    short_effort <- filter(short_effort, rowSums(is.na(short_effort)) != ncol(short_effort))
    
    dh_list <- list(short_dh, short_effort)
    
    return(dh_list)
  }
  #'  Create season-specific detection histories for species listed below
  spp_smr <- list("bear_black", "bobcat", "coyote", "mountain_lion", "wolf", "elk", "moose", "muledeer", "whitetaileddeer", "rabbit_hare")
  DHeff_npp20s_RNmod <- lapply(spp_smr, DH, dets = npp20s_det_events, cam_probs = npp20s_probs, start_date = "2020-07-01", y = "binary", rm_rows = rm_rows_npp20s, oc = 77) 
  DHeff_npp21s_RNmod <- lapply(spp_smr, DH, dets = npp21s_det_events, cam_probs = npp21s_probs, start_date = "2021-07-01", y = "binary", rm_rows = rm_rows_npp21s, oc = 77)
  DHeff_npp22s_RNmod <- lapply(spp_smr, DH, dets = npp22s_det_events, cam_probs = npp22s_probs, start_date = "2022-07-01", y = "binary", rm_rows = rm_rows_npp22s, oc = 77)
  
  #'  Remove sampling effort from species-specific DH_npp lists
  strip_list <- function(dh) {
    #'  Keep only the detection history per species
    dh_only <- dh[[1]]
    return(dh_only)
  }
  DH_npp20s_RNmod <- lapply(DHeff_npp20s_RNmod, strip_list)
  DH_npp21s_RNmod <- lapply(DHeff_npp21s_RNmod, strip_list)
  DH_npp22s_RNmod <- lapply(DHeff_npp22s_RNmod, strip_list)
  
  #'  --------------------
  #####  Sampling effort  #####
  #'  --------------------
  #'  Number of days per sampling occasion each camera was operational
  #' Second list generated by camtrapR's detection histories function
  sampling_effort <- function(dh) {
    camdays <- dh[[1]][[2]]
    effort <- as.data.frame(camdays) %>%
      rownames_to_column(var = "NewLocationID") %>%
      #'  Count total operational days and what that is in hours
      mutate(ndays = rowSums(.[-1], na.rm = T),
             nhrs = ndays*24)
    return(effort)
  }
  effort_20s_RNmod <- sampling_effort(DHeff_npp20s_RNmod)
  effort_21s_RNmod <- sampling_effort(DHeff_npp21s_RNmod)
  effort_22s_RNmod <- sampling_effort(DHeff_npp22s_RNmod)
  
  #' #'  ------------------------
  #' #####  SAVE detection data  #####
  #' #'  ------------------------
  #' #'  Seasonal detection histories
  #' save(DH_npp20s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp20s_RNmod.RData")
  #' save(DH_npp21s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp21s_RNmod.RData")
  #' save(DH_npp22s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp22s_RNmod.RData")
  #' #' Sampling effort
  #' save(effort_20s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp20s_RNmod.RData")
  #' save(effort_21s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp21s_RNmod.RData")
  #' save(effort_22s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp22s_RNmod.RData")
  #' 
  #' #'  Save for publication
  #' DetectionHist_Smr20 <- list(DH_npp20s_RNmod[[1]][[1]], DH_npp20s_RNmod[[2]][[1]], DH_npp20s_RNmod[[3]][[1]], DH_npp20s_RNmod[[4]][[1]], DH_npp20s_RNmod[[5]][[1]])
  #' DetectionHist_Smr21 <- list(DH_npp21s_RNmod[[1]][[1]], DH_npp21s_RNmod[[2]][[1]], DH_npp21s_RNmod[[3]][[1]], DH_npp21s_RNmod[[4]][[1]], DH_npp21s_RNmod[[5]][[1]])
  #' DetectionHist_Smr22 <- list(DH_npp22s_RNmod[[1]][[1]], DH_npp22s_RNmod[[2]][[1]], DH_npp22s_RNmod[[3]][[1]], DH_npp22s_RNmod[[4]][[1]], DH_npp22s_RNmod[[5]][[1]])
  #' SamplingEffort_Smr20 <- rn_effort_20s
  #' SamplingEffort_Smr21 <- rn_effort_21s
  #' SamplingEffort_Smr22 <- rn_effort_22s
  #' ####  UPDATE WORKING DIRECTORY EVENTUALLY
  #' # save(DetectionHist_Smr20, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/DetectionHist_smr20_RNmod.RData")
  #' # save(DetectionHist_Smr21, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/DetectionHist_smr21_RNmod.RData")
  #' # save(DetectionHist_Smr22, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/DetectionHist_smr22_RNmod.RData")
  #' # save(SamplingEffort_Smr20, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/SamplingEffort_smr20_RNmod.RData")
  #' # save(SamplingEffort_Smr21, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/SamplingEffort_smr21_RNmod.RData")
  #' # save(SamplingEffort_Smr22, file = "C:/Users/sbassing/OneDrive - University of Idaho/Repositories/Idaho_Predator_CoOccurrence/Data/SamplingEffort_smr22_RNmod.RData")

  #'  ---------------------------------
  #####  Format site-level covariates  #####
  #'  ---------------------------------
  #'  Load camera site covariates
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long_Smr2020-2022.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(NewLocationID != "GMU6_U_160" | CameraHeight_M != 1.2) %>%
    filter(NewLocationID != "GMU1_P_27" | CameraFacing != "atv")
  
  load("./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022.RData")
  
  #'  Reformat camera deployment data
  format_cam_station <- function(cams, season, habitat_covs) {
    cams <- dplyr::select(cams, c("NewLocationID", "Lat", "Long")) %>%
      mutate(GMU = str_extract(NewLocationID, "[^_]+"),
             Season = season) %>%
      left_join(cams_eoe_long[cams_eoe_long$Season == season,], by = c("NewLocationID", "Season")) %>%
      #'  Drop handful of duplicated rows (not sure why this happens with left_join)
      unique() %>%
      #'  Drop the few instances where camera height changed but site is the same
      group_by(NewLocationID) %>%
      slice(1L) %>%
      ungroup() %>%
      #'  Reduce categories to random, road, trail
      mutate(CameraFacing = ifelse(CameraFacing == "road" | CameraFacing == "atv" | CameraFacing == "gravel" | 
                                     CameraFacing == "decommision" | CameraFacing == "decommission", "road", CameraFacing),
             CameraFacing = ifelse(CameraFacing == "hiking" | CameraFacing == "game" | CameraFacing == "other", "trail", CameraFacing),
             CameraFacing = ifelse(is.na(CameraFacing), "trail", CameraFacing)) %>% # predator cameras so either trail or road
      arrange(NewLocationID) %>%
      left_join(habitat_covs, by = c("NewLocationID", "GMU"))
    return(cams)
  }
  cams_eoe20s <- format_cam_station(eoe_probcams_20s, season = "Smr20", habitat_covs = covariate_list[[1]])
  cams_eoe21s <- format_cam_station(eoe_probcams_21s, season = "Smr21", habitat_covs = covariate_list[[2]])
  cams_eoe22s <- format_cam_station(eoe_probcams_22s, season = "Smr22", habitat_covs = covariate_list[[3]])

  #'  Scale and format site-level covariates
  format_covs <- function(cams, dets, effort) {   
    #'  Remove rows where camera was inoperable the entire season - covariates at 
    #'  these sites shouldn't be included when scaling since they don't contribute
    #'  to detection data
    camsites <- row.names(dets)
    cam_covs <- cams[(cams$NewLocationID %in% camsites),]
    #'  Grab total number of days each camera was operable during season
    cam_operable <- effort %>% dplyr::select(c(NewLocationID, ndays))
    #'  Add total number of operable camera days to camera covariate dataset
    cam_covs <- left_join(cam_covs, cam_operable, by = "NewLocationID")
    #'  Rename, format, and scale as needed
    formatted <- cam_covs %>%
      mutate(GMUs = ifelse(GMU == "GMU10A", 1, 2),  # GMU10A represents the intercept!
             GMUs = ifelse(GMU == "GMU1", 3, GMUs),
             Setup = ifelse(Setup == "ungulate", 1, 2)) %>%  #'  Ungulate (random) cameras represent the intercept!
      transmute(NewLocationID = as.factor(NewLocationID),
                Season = as.factor(Season),
                GMU = as.factor(GMUs),
                Setup = as.factor(Setup),
                nDays = scale(ndays),
                PercFor = scale(perc_forest),
                Elev = scale(Elevation__10m2)) %>%
      #'  Arrange by camera location -- NECESSARY TO MATCH DH's CAMERA LOCATION ORDER
      arrange(NewLocationID) 
    
    return(formatted)
  }
  stations_npp20s <- format_covs(cams_eoe20s, dets = DHeff_npp20s_RNmod[[1]][[1]], effort = effort_20s_RNmod) 
  stations_npp21s <- format_covs(cams_eoe21s, dets = DHeff_npp21s_RNmod[[1]][[1]], effort = effort_21s_RNmod) 
  stations_npp22s <- format_covs(cams_eoe22s, dets = DHeff_npp22s_RNmod[[1]][[1]], effort = effort_22s_RNmod) 
  
  #'  Double check things are ordered correctly!!!!
  stations_npp22s[82:90,1:4]
  DH_npp22s_RNmod[[1]][82:90,1:3]
  
  #'  Double check percent forested habitat and elevation aren't too correlated
  cor(stations_npp20s$PercFor, stations_npp20s$Elev)
  cor(stations_npp21s$PercFor, stations_npp21s$Elev)
  cor(stations_npp22s$PercFor, stations_npp22s$Elev)
  
  #'  -----------------------------------
  #####  Format survey-level covariates  #####
  #'  -----------------------------------
  #'  Scale survey-level covariates (really only relevant if using 7-day sampling occasions)
  scale_srvy_cov <- function(srvy_covs) {
    #'  Reduce to just number of days per sampling occasion camera was operating
    skinny_srvy_covs <- srvy_covs %>% dplyr::select(-c(NewLocationID, ndays, nhrs))
    
    #'  Replace NAs with 0 - these sites truly not surveyed so effort is 0
    srvy_covs <- replace(skinny_srvy_covs, is.na(skinny_srvy_covs), 0)
  
    #'  Find mean & standard deviation of covariates across all sites & occasions
    mu <- mean(as.matrix(srvy_covs), na.rm = TRUE)
    sd <- sd(as.matrix(srvy_covs), na.rm = TRUE)
    
    #'  Z-transform (center observations around mean & scale by 1 SD)
    scaled <- ((srvy_covs - mu) / sd)
    scaled <- round(scaled, 3)
    
    scaled <- as.matrix(scaled)
    
    return(scaled)
  }
  effort_list <- list(effort_20s_RNmod, effort_21s_RNmod, effort_22s_RNmod)  
  zeffort_RNmod <- lapply(effort_list, scale_srvy_cov)
  
  #'  Create list of survey level covariates for unmarked (one list per year)
  srvy_covs_20s <- list(effort = zeffort_RNmod[[1]])
  srvy_covs_21s <- list(effort = zeffort_RNmod[[2]])
  srvy_covs_22s <- list(effort = zeffort_RNmod[[3]])
  
  
  #'  -----------------------
  ####  Royle-Nichols model  ####
  #'  -----------------------
  #'  Run Royle-Nichols model using a Bayesian statistical framework in JAGS and
  #'  a likelihood framework with unmarked
  
  #'  ------------------------
  #####  Setup data for JAGS  #####
  #'  ------------------------
  #'  Bundle detection histories and covariates for each species and year
  bundle_dat <- function(dh, nsite, nsurvey, cov, effort) {
    #'  Convert detection history to matrix
    dh <- as.matrix(dh)
    dimnames(dh) <- NULL
    #'  Count number of sites per GMU
    ncams_perGMU <- cov %>%
      group_by(GMU) %>%
      summarise(nsites = n()) %>%
      ungroup()
    #'  Split up covariates by GMU
    covs_GMU10A <- filter(cov, GMU == 1)
    covs_GMU6 <- filter(cov, GMU == 2)
    covs_GMU1 <- filter(cov, GMU == 3)
    #'  Bundle data for JAGS
    bundled <- list(y = dh, 
                    nsites = dim(dh)[1], 
                    nsurveys = dim(dh)[2], 
                    ngmu = max(as.numeric(cov$GMU)),
                    nsets = max(as.numeric(cov$Setup)),
                    ncams1 = as.numeric(ncams_perGMU[1,2]), # GMU10A
                    ncams2 = as.numeric(ncams_perGMU[2,2]), # GMU6
                    ncams3 = as.numeric(ifelse(is.na(ncams_perGMU[3,2]), 0, ncams_perGMU[3,2])), #GMU1
                    gmu = as.numeric(cov$GMU), 
                    setup = as.numeric(cov$Setup),
                    forest = as.numeric(cov$PercFor), 
                    elev = as.numeric(cov$Elev), 
                    ndays = as.numeric(cov$nDays), 
                    seffort = as.numeric(effort),
                    #'  GMU-specific covariates for predicting N
                    forest1 = as.numeric(covs_GMU10A$PercFor),
                    forest2 = as.numeric(covs_GMU6$PercFor),
                    forest3 = as.numeric(covs_GMU1$PercFor),
                    elev1 = as.numeric(covs_GMU10A$Elev),
                    elev2 = as.numeric(covs_GMU6$Elev),
                    elev3 = as.numeric(covs_GMU1$Elev),
                    #'  Area of each (km2)
                    area1 = as.numeric(8527.31),
                    area2 = as.numeric(5905.44),
                    area3 = as.numeric(14648.92))
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle_20s <- lapply(DH_npp20s_RNmod, bundle_dat, cov = stations_npp20s, effort = zeffort_RNmod[[1]])
  data_JAGS_bundle_21s <- lapply(DH_npp21s_RNmod, bundle_dat, cov = stations_npp21s, effort = zeffort_RNmod[[2]])
  data_JAGS_bundle_22s <- lapply(DH_npp22s_RNmod, bundle_dat, cov = stations_npp22s, effort = zeffort_RNmod[[3]])
  
  save(data_JAGS_bundle_20s, file = "./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_20s.RData")
  save(data_JAGS_bundle_21s, file = "./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_21s.RData")
  save(data_JAGS_bundle_22s, file = "./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_22s.RData")
  
  #'  Initial values
  #'  Using naive occupancy as a starting point for local abundance
  initial_n <- function(dh) {
    #'  Max value per row
    ninit <- apply(dh, 1, max, na.rm = TRUE)
    ninit <- as.vector(ninit)
    return(ninit)
  }
  #'  Apply function per species for each year
  ninit_20s <- lapply(DH_npp20s_RNmod, initial_n)
  ninit_21s <- lapply(DH_npp21s_RNmod, initial_n)
  ninit_22s <- lapply(DH_npp22s_RNmod, initial_n)
  
  #'  Parameters monitored
  params <- c("beta0", "beta1", "beta2", "beta3", "beta4", #"mean.lambda", 
              "alpha0", "alpha1", "alpha2", "rSetup", "mu.r", "mean.p", #"mean.r", 
              "mu.lambda", "totalN", "occSites", "mean.psi", 
              "totalN.gmu10a", "densitykm2.gmu10a", "density100km2.gmu10a", 
              "totalN.gmu6", "densitykm2.gmu6", "density100km2.gmu6", 
              "totalN.gmu1", "densitykm2.gmu1", "density100km2.gmu1","N")
  #'  NOTE about mean vs mu lambda and r: 
  #'  mean.lambda = the intercept, i.e., mean lambda for GMU10A 
  #'  mean.r = the intercept, i.e., per-individual detection probability at random sites
  #'  mu.lambda = lambda averaged across all GMUs
  #'  mu.r = per-individual detection probability averaged across all sites 
  
  #'  MCMC settings
  nc <- 3
  ni <- 50000
  nb <- 10000
  nt <- 10
  na <- 5000
  
  #'  ---------------------------
  #####  Run RN model with JAGS  #####
  #'  ---------------------------
  ######  2020  Analyses  ######
  source("./Scripts/Relative_Abundance/RNmodel_JAGS_code_2020.R")
  
  ######  Black bear  ######
  start.time = Sys.time()
  inits_bear20s <- function(){list(N = ninit_20s[[1]])}
  RN_bear_20s <- jags(data_JAGS_bundle_20s[[1]], inits = inits_bear20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_20s$summary)
  which(RN_bear_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_20s$samples)
  save(RN_bear_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_bear_20s_", Sys.Date(), ".RData"))

  ######  Bobcat  ######
  start.time = Sys.time()
  inits_bob20s <- function(){list(N = ninit_20s[[2]])}
  RN_bob_20s <- jags(data_JAGS_bundle_20s[[2]], inits = inits_bob20s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bob_20s$summary)
  which(RN_bob_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bob_20s$samples)
  save(RN_bob_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_bob_20s_", Sys.Date(), ".RData")) 
  
  ######  Coyote  ######
  start.time = Sys.time()
  inits_coy20s <- function(){list(N = ninit_20s[[3]])}
  RN_coy_20s <- jags(data_JAGS_bundle_20s[[3]], inits = inits_coy20s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_20s$summary)
  which(RN_coy_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_20s$samples)
  save(RN_coy_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_coy_20s_", Sys.Date(), ".RData"))
  
  ######  Mountain lion  ######
  start.time = Sys.time()
  inits_lion20s <- function(){list(N = ninit_20s[[4]])}
  RN_lion_20s <- jags(data_JAGS_bundle_20s[[4]], inits = inits_lion20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_20s$summary)
  which(RN_lion_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_20s$samples)
  save(RN_lion_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_lion_20s_", Sys.Date(), ".RData")) 
  
  ######  Wolf  ######
  start.time = Sys.time()
  inits_wolf20s <- function(){list(N = ninit_20s[[5]])}
  RN_wolf_20s <- jags(data_JAGS_bundle_20s[[5]], inits = inits_wolf20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_20s$summary)
  which(RN_wolf_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_20s$samples)
  save(RN_wolf_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_wolf_20s_", Sys.Date(), ".RData")) 
  
  
  ######  2021 & 2022  Analyses  ######
  source("./Scripts/Relative_Abundance/RNmodel_JAGS_code.R")
  
  ######  Black bear  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_bear21s <- function(){list(N = ninit_21s[[1]])}
  RN_bear_21s <- jags(data_JAGS_bundle_21s[[1]], inits = inits_bear21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_21s$summary)
  which(RN_bear_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_21s$samples)
  save(RN_bear_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_bear_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_bear22s <- function(){list(N = ninit_22s[[1]])}
  RN_bear_22s <- jags(data_JAGS_bundle_22s[[1]], inits = inits_bear22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_22s$summary)
  which(RN_bear_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_22s$samples)
  save(RN_bear_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_bear_22s_", Sys.Date(), ".RData"))
  
  ######  Bobcat  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_bob21s <- function(){list(N = ninit_21s[[2]])}
  RN_bob_21s <- jags(data_JAGS_bundle_21s[[2]], inits = inits_bob21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bob_21s$summary)
  which(RN_bob_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bob_21s$samples)
  save(RN_bob_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_bob_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_bob22s <- function(){list(N = ninit_22s[[2]])}
  RN_bob_22s <- jags(data_JAGS_bundle_22s[[2]], inits = inits_bob22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bob_22s$summary)
  which(RN_bob_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bob_22s$samples)
  save(RN_bob_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_bob_22s_", Sys.Date(), ".RData"))
  
  ######  Coyote  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_coy21s <- function(){list(N = ninit_21s[[3]])}
  RN_coy_21s <- jags(data_JAGS_bundle_21s[[3]], inits = inits_coy21s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_21s$summary)
  which(RN_coy_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_21s$samples)
  save(RN_coy_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_coy_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_coy22s <- function(){list(N = ninit_22s[[3]])}
  RN_coy_22s <- jags(data_JAGS_bundle_22s[[3]], inits = inits_coy22s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_22s$summary)
  which(RN_coy_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_22s$samples)
  save(RN_coy_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_coy_22s_", Sys.Date(), ".RData"))
  
  ######  Mountain lion  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_lion21s <- function(){list(N = ninit_21s[[4]])}
  RN_lion_21s <- jags(data_JAGS_bundle_21s[[4]], inits = inits_lion21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_21s$summary)
  which(RN_lion_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_21s$samples)
  save(RN_lion_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_lion_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_lion22s <- function(){list(N = ninit_22s[[4]])}
  RN_lion_22s <- jags(data_JAGS_bundle_22s[[4]], inits = inits_lion22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_22s$summary)
  which(RN_lion_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_22s$samples)
  save(RN_lion_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_lion_22s_", Sys.Date(), ".RData"))
  
  ######  Wolf  ######
  ni_wolf <- 75000
  #'  Summer 2021
  start.time = Sys.time()
  inits_wolf21s <- function(){list(N = ninit_21s[[5]])}
  RN_wolf_21s <- jags(data_JAGS_bundle_21s[[5]], inits = inits_wolf21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_21s$summary)
  which(RN_wolf_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_21s$samples)
  save(RN_wolf_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_wolf_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_wolf22s <- function(){list(N = ninit_22s[[5]])}
  RN_wolf_22s <- jags(data_JAGS_bundle_22s[[5]], inits = inits_wolf22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_22s$summary)
  which(RN_wolf_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_22s$samples)
  save(RN_wolf_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_wolf_22s_", Sys.Date(), ".RData"))
  
  
    
  #'  ----------------------------
  #####  Setup data for unmarked  #####
  #'  ----------------------------
  #'  Create list of survey level covariates for unmarked (one list per year)
  srvy_covs_20s <- list(effort = zeffort_RNmod[[1]])
  srvy_covs_21s <- list(effort = zeffort_RNmod[[2]])
  srvy_covs_22s <- list(effort = zeffort_RNmod[[3]])
  
  #'  Function to create unmarked dataframes for Royle-Nichols occupancy models
  #'  Single-season, single-species occupancy unmarkedDF --> unmarkedFrameOccu 
  umf_setup <- function(dh, listnames, sitecovs, srvycovs, plotit = T) {
    
    #'  Create array with detection histories, site-level and survey-level covariates
    UMF <- unmarkedFrameOccu(y = dh, siteCovs = sitecovs, obsCovs = srvycovs)
    
    #'  Plot detection histories
    print(plot(UMF))
    #'  Summarize data
    summary(UMF)
    
    return(UMF)
  }
  bear_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[1]][[1]], sitecovs = stations_npp20s, srvycovs = srvy_covs_20s)
  bear_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[1]][[1]], sitecovs = stations_npp21s, srvycovs = srvy_covs_21s)
  bear_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[1]][[1]], sitecovs = stations_npp22s, srvycovs = srvy_covs_22s)
  
  bob_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[2]][[1]], sitecovs = stations_npp20s, srvycovs = srvy_covs_20s)
  bob_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[2]][[1]], sitecovs = stations_npp21s, srvycovs = srvy_covs_21s)
  bob_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[2]][[1]], sitecovs = stations_npp22s, srvycovs = srvy_covs_22s)
  
  coy_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[3]][[1]], sitecovs = stations_npp20s, srvycovs = srvy_covs_20s)
  coy_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[3]][[1]], sitecovs = stations_npp21s, srvycovs = srvy_covs_21s)
  coy_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[3]][[1]], sitecovs = stations_npp22s, srvycovs = srvy_covs_22s)
  
  lion_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[4]][[1]], sitecovs = stations_npp20s, srvycovs = srvy_covs_20s)
  lion_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[4]][[1]], sitecovs = stations_npp21s, srvycovs = srvy_covs_21s)
  lion_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[4]][[1]], sitecovs = stations_npp22s, srvycovs = srvy_covs_22s)
  
  wolf_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[5]][[1]], sitecovs = stations_npp20s, srvycovs = srvy_covs_20s)
  wolf_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[5]][[1]], sitecovs = stations_npp21s, srvycovs = srvy_covs_21s)
  wolf_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[5]][[1]], sitecovs = stations_npp22s, srvycovs = srvy_covs_22s)

  #'  -------------------------------
  #####  Run RN model with unmarked  #####  
  #'  -------------------------------
  #'  Fit Royle-Nichles model (RN model) to detection/non-detection data to estimate
  #'  relative abundance and detection probability per species and season
  #'  ~ covariates affecting per-individual detection probability ~ covariates affecting abundance
  #'  Use "effort" covariate if sampling occasion = 7 days, "nDays" if sampling occasion = 1 day
  #'  Including all variables (even if not significant for all spp/years) to make as predictive as possible
  #'  Not using model selection here b/c testing competing models with SEM stage of analysis,
  #'  this is simply to derive a predicted representation of variation in predator abundance
  (bear20s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), bear_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (bear21s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), bear_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (bear22s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), bear_22s_umf, se = TRUE, engine = "C", threads = 4))
  bear_rn_list <- list(bear20s.rn, bear21s.rn, bear22s.rn)
  
  (bob20s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), bob_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (bob21s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), bob_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (bob22s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), bob_22s_umf, se = TRUE, engine = "C", threads = 4))
  bob_rn_list <- list(bob20s.rn, bob21s.rn, bob22s.rn)
  
  (coy20s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), coy_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (coy21s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), coy_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (coy22s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), coy_22s_umf, se = TRUE, engine = "C", threads = 4))
  coy_rn_list <- list(coy20s.rn, coy21s.rn, coy22s.rn)
  
  (lion20s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), lion_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (lion21s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), lion_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (lion22s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), lion_22s_umf, se = TRUE, engine = "C", threads = 4)) 
  lion_rn_list <- list(lion20s.rn, lion21s.rn, lion22s.rn)
  
  (wolf20s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), wolf_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (wolf21s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), wolf_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (wolf22s.rn <- occuRN(~ Setup + nDays ~ GMU + PercFor + Elev + I(Elev^2), wolf_22s_umf, se = TRUE, engine = "C", threads = 4))
  wolf_rn_list <- list(wolf20s.rn, wolf21s.rn, wolf22s.rn)
  
  #'  Empirical Bayes estimates of abundance at each site
  bear20s_re <- ranef(bear20s.rn, K = 50)
  plot(bear20s_re)
  
  # Estimates of conditional abundance distribution at each site
  (re <- ranef(bear20s.rn))
  # Best Unbiased Predictors
  bup(re, stat="mean") # Posterior mean
  bup(re, stat="mode") # Posterior mode
  confint(re, level=0.95) # 90% CI
  plot(re, subset=site %in% c(1:10), layout=c(5, 2), xlim=c(-1,20))
  
  