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
  df_all_20s <- thin_detections(eoe20s_allM, seqprobs = eoe_seqprob_20s, start_date = "2020-06-01", end_date = "2020-09-15") 
  df_all_21s <- thin_detections(eoe21s_allM, seqprobs = eoe_seqprob_21s, start_date = "2021-06-01", end_date = "2021-09-15")
  df_all_22s <- thin_detections(eoe22s_allM, seqprobs = eoe_seqprob_22s, start_date = "2022-06-01", end_date = "2022-09-15") 

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
  
  #'  Save
  save(npp20s_det_events, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/Unique_detection_events_20s.RData")
  save(npp21s_det_events, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/Unique_detection_events_21s.RData")
  save(npp22s_det_events, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/Unique_detection_events_22s.RData")
    
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
  #'  FYI: June 1 - Sept 15 = 107 1-day sampling periods / July 1 - Sept 15 = 77 1-day sampling periods
  DH <- function(dets, cam_probs, spp, start_date, y, rm_rows, oc) {
    dets <- dets %>% arrange(NewLocationID, posix_date_time)
    det_hist <- detectionHistory(recordTable = dets,
                                 camOp = cam_probs,
                                 stationCol = "NewLocationID",
                                 recordDateTimeCol = "posix_date_time",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 species = spp,
                                 occasionLength = 1,
                                 day1 = start_date, 
                                 datesAsOccasionNames = FALSE,
                                 timeZone = "America/Edmonton",
                                 output = y,
                                 includeEffort = TRUE,
                                 scaleEffort = FALSE,
                                 outDir = "./Data/MultiSpp_OccMod_Outputs/Detection_Histories")
    
    #'  Reduce detection histories to sampling occasions of interest (drop extra
    #'  occasions after focal period of interest)
    short_dh <- det_hist[[1]][,1:oc]
    short_effort <- det_hist[[2]][,1:oc]
    
    #' #'  Remove rows where camera was inoperable during the entire sampling season
    #' short_dh <- as.data.frame(short_dh)
    #' short_effort <- as.data.frame(short_effort)
    #' short_dh <- filter(short_dh, rowSums(is.na(short_dh)) != ncol(short_dh))
    #' short_effort <- filter(short_effort, rowSums(is.na(short_effort)) != ncol(short_effort))
    
    #'  Remove rows where camera was inoperable during the entire sampling season
    short_dh <- as.data.frame(short_dh)
    short_effort <- as.data.frame(short_effort)  %>%
      mutate(NewLocationID = row.names(.)) %>% relocate(NewLocationID, .before = "o1")
    short_dh <- filter(short_dh, rowSums(is.na(short_dh)) != ncol(short_dh)) %>%
      mutate(NewLocationID = row.names(.)) %>% relocate(NewLocationID, .before = "o1")
    short_effort <- short_effort[short_effort$NewLocationID %in% short_dh$NewLocationID,]
    
    dh_list <- list(short_dh, short_effort)
    
    return(dh_list)
  }
  #'  Create season-specific detection histories for species listed below
  spp_smr <- list("bear_black", "bobcat", "coyote", "mountain_lion", "wolf", "elk", "moose", "muledeer", "whitetaileddeer", "rabbit_hare")
  DHeff_npp20s_RNmod <- lapply(spp_smr, DH, dets = npp20s_det_events, cam_probs = npp20s_probs, start_date = "2020-06-01", y = "binary", oc = 107) #rm_rows = rm_rows_npp20s, 
  DHeff_npp21s_RNmod <- lapply(spp_smr, DH, dets = npp21s_det_events, cam_probs = npp21s_probs, start_date = "2021-06-01", y = "binary", oc = 107) #rm_rows = rm_rows_npp21s, 
  DHeff_npp22s_RNmod <- lapply(spp_smr, DH, dets = npp22s_det_events, cam_probs = npp22s_probs, start_date = "2022-06-01", y = "binary", oc = 107) #rm_rows = rm_rows_npp22s, 
  
  #'  Remove sampling effort from species-specific DH_npp lists
  strip_list <- function(dh) {
    #'  Keep only the detection history per species
    dh_only <- dh[[1]]
    dh_only <- dh_only %>% dplyr::select(-NewLocationID)
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
      dplyr::select(-NewLocationID) %>%
      rownames_to_column(var = "NewLocationID") %>%
      #'  Count total operational days and what that is in hours
      mutate(ndays = rowSums(.[-1], na.rm = T),
             nhrs = ndays*24)
    return(effort)
  }
  effort_20s_RNmod <- sampling_effort(DHeff_npp20s_RNmod)
  effort_21s_RNmod <- sampling_effort(DHeff_npp21s_RNmod)
  effort_22s_RNmod <- sampling_effort(DHeff_npp22s_RNmod)
  
  #'  ------------------------
  #####  SAVE detection data  #####
  #'  ------------------------
  #'  Seasonal detection histories
  save(DH_npp20s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp20s_RNmod.RData")
  save(DH_npp21s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp21s_RNmod.RData")
  save(DH_npp22s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp22s_RNmod.RData")
  #' Sampling effort
  save(effort_20s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp20s_RNmod.RData")
  save(effort_21s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp21s_RNmod.RData")
  save(effort_22s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp22s_RNmod.RData")


  #'  ---------------------------------
  #####  Format site-level covariates  #####
  #'  ---------------------------------
  #'  Load camera site covariates
  cams_eoe_long <- read.csv("./Data/IDFG camera data/cams_eoe_long_Smr2020-2022.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Season", "CameraHeight_M", "CameraFacing")) %>%
    filter(NewLocationID != "GMU6_U_160" | CameraHeight_M != 1.2) %>%
    filter(NewLocationID != "GMU1_P_27" | CameraFacing != "atv") %>%
    filter(Season != "Wtr20") %>%
    arrange(NewLocationID, Season)
  
  load("./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022_Feb2024.RData")
  
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
    cam_covs <- full_join(cam_covs, cam_operable, by = "NewLocationID")
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
                Elev = scale(Elevation__10m2)) #%>%
      #' #'  Arrange by camera location -- NECESSARY TO MATCH DH's CAMERA LOCATION ORDER
      #' arrange(NewLocationID) 
    
    return(formatted)
  }
  stations_npp20s <- format_covs(cams_eoe20s, dets = DH_npp20s_RNmod[[1]], effort = effort_20s_RNmod) 
  stations_npp21s <- format_covs(cams_eoe21s, dets = DH_npp21s_RNmod[[1]], effort = effort_21s_RNmod) 
  stations_npp22s <- format_covs(cams_eoe22s, dets = DH_npp22s_RNmod[[1]], effort = effort_22s_RNmod) 
  
  #'  Double check things are ordered correctly!!!!
  stations_npp20s[82:90,1:4]; DH_npp20s_RNmod[[1]][82:90,1:3]; nrow(stations_npp20s); nrow(DH_npp20s_RNmod[[1]])
  stations_npp21s[82:90,1:4]; DH_npp21s_RNmod[[1]][82:90,1:3]; nrow(stations_npp21s); nrow(DH_npp21s_RNmod[[1]])
  stations_npp22s[82:90,1:4]; DH_npp22s_RNmod[[1]][82:90,1:3]; nrow(stations_npp22s); nrow(DH_npp22s_RNmod[[1]])
  
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
              "alpha0", "alpha2", "rSetup", "mu.r", "mean.p", #"alpha1" #"mean.r", 
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
  ni_wolf <-  100000 
  start.time = Sys.time()
  inits_wolf20s <- function(){list(N = ninit_20s[[5]])}
  RN_wolf_20s <- jags(data_JAGS_bundle_20s[[5]], inits = inits_wolf20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_20s$summary)
  which(RN_wolf_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_20s$samples)
  save(RN_wolf_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_wolf_20s_", Sys.Date(), ".RData")) 
  
  ######  Elk  ######
  start.time = Sys.time()
  inits_elk20s <- function(){list(N = ninit_20s[[6]])}
  RN_elk_20s <- jags(data_JAGS_bundle_20s[[6]], inits = inits_elk20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_20s$summary)
  which(RN_elk_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_20s$samples)
  save(RN_elk_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_elk_20s_", Sys.Date(), ".RData"))
  
  ######  Moose  ######
  start.time = Sys.time()
  inits_moose20s <- function(){list(N = ninit_20s[[7]])}
  RN_moose_20s <- jags(data_JAGS_bundle_20s[[7]], inits = inits_moose20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_20s$summary)
  which(RN_moose_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_20s$samples)
  save(RN_moose_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_moose_20s_", Sys.Date(), ".RData"))
  
  ######  White-tailed Deer  ######
  ni_wtd <- 100000
  start.time = Sys.time()
  inits_wtd20s <- function(){list(N = ninit_20s[[9]])}
  RN_wtd_20s <- jags(data_JAGS_bundle_20s[[9]], inits = inits_wtd20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_20s$summary)
  which(RN_wtd_20s$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_wtd_20s$samples)
  save(RN_wtd_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_wtd_20s_", Sys.Date(), ".RData"))
  
  ######  Lagomorphs  ######
  ni_lag <- 100000
  start.time = Sys.time()
  inits_lago20s <- function(){list(N = ninit_20s[[10]])}
  RN_lago_20s <- jags(data_JAGS_bundle_20s[[10]], inits = inits_lago20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_lag, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lago_20s$summary)
  which(RN_lago_20s$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_lago_20s$samples)
  save(RN_lago_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_lagomorphs_20s_", Sys.Date(), ".RData"))
  
  
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
  ni_wolf <- 100000 
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
  
  ######  Elk  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_elk21s <- function(){list(N = ninit_21s[[6]])}
  RN_elk_21s <- jags(data_JAGS_bundle_21s[[6]], inits = inits_elk21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_21s$summary)
  which(RN_elk_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_21s$samples)
  save(RN_elk_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_elk_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_elk22s <- function(){list(N = ninit_22s[[6]])}
  RN_elk_22s <- jags(data_JAGS_bundle_22s[[6]], inits = inits_elk22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_22s$summary)
  which(RN_elk_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_22s$samples)
  save(RN_elk_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_elk_22s_", Sys.Date(), ".RData"))
  
  ######  Moose  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_moose21s <- function(){list(N = ninit_21s[[7]])}
  RN_moose_21s <- jags(data_JAGS_bundle_21s[[7]], inits = inits_moose21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_21s$summary)
  which(RN_moose_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_21s$samples)
  save(RN_moose_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_moose_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_moose22s <- function(){list(N = ninit_22s[[7]])}
  RN_moose_22s <- jags(data_JAGS_bundle_22s[[7]], inits = inits_moose22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_22s$summary)
  which(RN_moose_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_22s$samples)
  save(RN_moose_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_moose_22s_", Sys.Date(), ".RData"))
  
  ######  White-tailed Deer  ######
  ni_wtd <- 100000
  #'  Summer 2021
  start.time = Sys.time()
  inits_wtd21s <- function(){list(N = ninit_21s[[9]])}
  RN_wtd_21s <- jags(data_JAGS_bundle_21s[[9]], inits = inits_wtd21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_21s$summary)
  which(RN_wtd_21s$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_wtd_21s$samples)
  save(RN_wtd_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_wtd_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_wtd22s <- function(){list(N = ninit_22s[[9]])}
  RN_wtd_22s <- jags(data_JAGS_bundle_22s[[9]], inits = inits_wtd22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_22s$summary)
  which(RN_wtd_22s$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_wtd_22s$samples)
  save(RN_wtd_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_wtd_22s_", Sys.Date(), ".RData"))
  
  ######  Lagomorphs  ######
  ni_lag <- 100000
  #'  Summer 2021
  start.time = Sys.time()
  inits_lago21s <- function(){list(N = ninit_21s[[10]])}
  RN_lago_21s <- jags(data_JAGS_bundle_21s[[10]], inits = inits_lago21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_lag, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lago_21s$summary)
  which(RN_lago_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lago_21s$samples)
  save(RN_lago_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_lagomorphs_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_lago22s <- function(){list(N = ninit_22s[[10]])}
  RN_lago_22s <- jags(data_JAGS_bundle_22s[[10]], inits = inits_lago22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_lag, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lago_22s$summary)
  which(RN_lago_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lago_22s$samples)
  save(RN_lago_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_lagomorphs_22s_", Sys.Date(), ".RData"))
  
  
  #'  --------------------------
  ####  Summarize density data  ####
  #'  --------------------------
  #'  Load model outputs
  filenames <- list.files("./Outputs/Relative_Abundance/RN_model/JAGS_out", pattern="*.RData", full.names=TRUE)
  lapply(filenames, load, environment())
  
  #'  Format GMU specific mean & sd
  gmu_densities <- function(mod, spp, yr) {
    #'  Pull out relevant summary stats 
    gmu10a_mean <- round(mod$mean$density100km2.gmu10a, 2)
    gmu10a_sd <- round(mod$sd$density100km2.gmu10a, 3)
    gmu6_mean <- round(mod$mean$density100km2.gmu6, 2)
    gmu6_sd <- round(mod$sd$density100km2.gmu6, 3)
    gmu1_mean <- round(mod$mean$density100km2.gmu1, 2)
    gmu1_sd <- round(mod$sd$density100km2.gmu1, 3)
    
    #'  Create data frame
    mod_summary <- cbind(spp, yr, gmu1_mean, gmu1_sd, gmu6_mean, gmu6_sd, gmu10a_mean, gmu10a_sd)
    mod_summary <- as.data.frame(mod_summary) %>%
      mutate(gmu1_mean = ifelse(yr == "2020", NA, gmu1_mean),
             gmu1_sd = ifelse(yr == "2020", NA, gmu1_sd))
    names(mod_summary) <- c("Species", "Year", "GMU1 avg. density", "GMU1 SD", 
                            "GMU6 avg. density", "GMU6 SD", "GMU10A avg. density", "GMU10A SD")
    
    return(mod_summary)
  }
  bear_density_20s <- gmu_densities(RN_bear_20s, spp = "black bear", yr = 2020)
  bear_density_21s <- gmu_densities(RN_bear_21s, spp = "black bear", yr = 2021)
  bear_density_22s <- gmu_densities(RN_bear_22s, spp = "black bear", yr = 2022)
  bob_density_20s <- gmu_densities(RN_bob_20s, spp = "bobcat", yr = 2020)
  bob_density_21s <- gmu_densities(RN_bob_21s, spp = "bobcat", yr = 2021)
  bob_density_22s <- gmu_densities(RN_bob_22s, spp = "bobcat", yr = 2022)
  coy_density_20s <- gmu_densities(RN_coy_20s, spp = "coyote", yr = 2020)
  coy_density_21s <- gmu_densities(RN_coy_21s, spp = "coyote", yr = 2021)
  coy_density_22s <- gmu_densities(RN_coy_22s, spp = "coyote", yr = 2022)
  lion_density_20s <- gmu_densities(RN_lion_20s, spp = "mountain lion", yr = 2020)
  lion_density_21s <- gmu_densities(RN_lion_21s, spp = "mountain lion", yr = 2021)
  lion_density_22s <- gmu_densities(RN_lion_22s, spp = "mountain lion", yr = 2022)
  wolf_density_20s <- gmu_densities(RN_wolf_20s, spp = "wolf", yr = 2020)
  wolf_density_21s <- gmu_densities(RN_wolf_21s, spp = "wolf", yr = 2021)
  wolf_density_22s <- gmu_densities(RN_wolf_22s, spp = "wolf", yr = 2022)
  
  elk_density_20s <- gmu_densities(RN_elk_20s, spp = "elk", yr = 2020)
  elk_density_21s <- gmu_densities(RN_elk_21s, spp = "elk", yr = 2021)
  elk_density_22s <- gmu_densities(RN_elk_22s, spp = "elk", yr = 2022)
  lago_density_20s <- gmu_densities(RN_lago_20s, spp = "lagomorphs", yr = 2020)
  lago_density_21s <- gmu_densities(RN_lago_21s, spp = "lagomorphs", yr = 2021)
  lago_density_22s <- gmu_densities(RN_lago_22s, spp = "lagomorphs", yr = 2022)
  moose_density_20s <- gmu_densities(RN_moose_20s, spp = "moose", yr = 2020)
  moose_density_21s <- gmu_densities(RN_moose_21s, spp = "moose", yr = 2021)
  moose_density_22s <- gmu_densities(RN_moose_22s, spp = "moose", yr = 2022)
  wtd_density_20s <- gmu_densities(RN_wtd_20s, spp = "white-tailed deer", yr = 2020)
  wtd_density_21s <- gmu_densities(RN_wtd_21s, spp = "white-tailed deer", yr = 2021)
  wtd_density_22s <- gmu_densities(RN_wtd_22s, spp = "white-tailed deer", yr = 2022)
  
  #'  Merge all together
  RN_density <- rbind(bear_density_20s, bear_density_21s, bear_density_22s, bob_density_20s, bob_density_21s, bob_density_22s,
                      coy_density_20s, coy_density_21s, coy_density_22s, lion_density_20s, lion_density_21s, lion_density_22s,
                      wolf_density_20s, wolf_density_21s, wolf_density_22s, elk_density_20s, elk_density_21s, elk_density_22s,
                      lago_density_20s, lago_density_21s, lago_density_22s, moose_density_20s, moose_density_21s, moose_density_22s,
                      wtd_density_20s, wtd_density_21s, wtd_density_22s)
  #'  GMU1
  RN_density_gmu1 <- dplyr::select(RN_density, c("Species", "Year", "GMU1 avg. density", "GMU1 SD")) %>%
    rename("Mean density" = "GMU1 avg. density") %>%
    rename("SD" = "GMU1 SD") %>%
    filter(Year != "2020") %>%
    arrange(Year, `Mean density`)
  RN_density_gmu6 <- dplyr::select(RN_density, c("Species", "Year", "GMU6 avg. density", "GMU6 SD")) %>%
    rename("Mean density" = "GMU6 avg. density") %>%
    rename("SD" = "GMU6 SD") %>%
    arrange(Year, `Mean density`)
  RN_density_gmu10a <- dplyr::select(RN_density, c("Species", "Year", "GMU10A avg. density", "GMU10A SD")) %>%
    rename("Mean density" = "GMU10A avg. density") %>%
    rename("SD" = "GMU10A SD") %>%
    arrange(Year, `Mean density`)
  
  write_csv(RN_density_gmu1, "./Outputs/Relative_Abundance/RN_model/RN_density_gmu1.csv")
  write_csv(RN_density_gmu6, "./Outputs/Relative_Abundance/RN_model/RN_density_gmu6.csv")
  write_csv(RN_density_gmu10a, "./Outputs/Relative_Abundance/RN_model/RN_density_gmu10a.csv")
  
  #'  List model outputs together 
  rn_bear_list <- list(RN_bear_20s, RN_bear_21s, RN_bear_22s)
  rn_bob_list <- list(RN_bob_20s, RN_bob_21s, RN_bob_22s)
  rn_coy_list <- list(RN_coy_20s, RN_coy_21s, RN_coy_22s)
  rn_lion_list <- list(RN_lion_20s, RN_lion_21s, RN_lion_22s)
  rn_wolf_list <- list(RN_wolf_20s, RN_wolf_21s, RN_wolf_22s)
  
  rn_elk_list <- list(RN_elk_20s, RN_elk_21s, RN_elk_22s)
  rn_lago_list <- list(RN_lago_20s, RN_lago_21s, RN_lago_22s)
  rn_moose_list <- list(RN_moose_20s, RN_moose_21s, RN_moose_22s)
  rn_wtd_list <- list(RN_wtd_20s, RN_wtd_21s, RN_wtd_22s)
  
  #'  Load saved detection histories
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp20s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp21s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp22s_RNmod.RData")
  
  #'  List one detection history per year (doesn't matter which species it relates to)
  dh_list <- list(DH_npp20s_RNmod[[1]], DH_npp21s_RNmod[[1]], DH_npp22s_RNmod[[1]])
  
  #'  Save estimated N per site
  estimated_N <- function(mod, dh, spp) {
    #'  Grab estimated N and SD per site
    RN.n <- mod$mean$N
    RN.sd <- mod$sd$N
    #'  Grab camera location
    # dh <- dh[[1]]
    locs <- rownames(dh)
    #'  Merge and format into single data frame with corresponding N & SD per site
    out <- cbind(locs, RN.n, RN.sd)
    RN_est <- as.data.frame(out) %>%
      mutate(NewLocationID = locs, 
             Species = spp,
             RN.n = as.numeric(RN.n),
             RN.sd = as.numeric(RN.sd)) %>%
      separate(locs, c("GMU", "Setup", "site"), sep = "_") %>%
      mutate(CellID = paste0(GMU, "_", site)) %>%
      dplyr::select(-site) %>%
      relocate(NewLocationID, .before = "GMU") %>%
      relocate(Species, .after = "Setup") %>%
      relocate(CellID, .after = "NewLocationID")
    return(RN_est)
  }
  rn_bear_out <- mapply(estimated_N, rn_bear_list, dh = dh_list, spp = "bear_black", SIMPLIFY = FALSE)
  rn_bob_out <- mapply(estimated_N, rn_bob_list, dh = dh_list, spp = "bobcat", SIMPLIFY = FALSE)
  rn_coy_out <- mapply(estimated_N, rn_coy_list, dh = dh_list, spp = "coyote", SIMPLIFY = FALSE)
  rn_lion_out <- mapply(estimated_N, rn_lion_list, dh = dh_list, spp = "mountain_lion", SIMPLIFY = FALSE)
  rn_wolf_out <- mapply(estimated_N, rn_wolf_list, dh = dh_list, spp = "wolf", SIMPLIFY = FALSE)
  
  rn_elk_out <- mapply(estimated_N, rn_elk_list, dh = dh_list, spp = "elk", SIMPLIFY = FALSE)
  rn_lago_out <- mapply(estimated_N, rn_lago_list, dh = dh_list, spp = "lagomorphs", SIMPLIFY = FALSE)
  rn_moose_out <- mapply(estimated_N, rn_moose_list, dh = dh_list, spp = "moose", SIMPLIFY = FALSE)
  rn_wtd_out <- mapply(estimated_N, rn_wtd_list, dh = dh_list, spp = "whitetailed_deer", SIMPLIFY = FALSE)
  
  rn_2020 <- rbind(rn_bear_out[[1]], rn_bob_out[[1]], rn_coy_out[[1]], rn_lion_out[[1]], rn_wolf_out[[1]],
                   rn_elk_out[[1]], rn_lago_out[[1]], rn_moose_out[[1]], rn_wtd_out[[1]]) %>%
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr20") %>%
    relocate(season, .after = Setup)
  rn_2021 <- rbind(rn_bear_out[[2]], rn_bob_out[[2]], rn_coy_out[[2]], rn_lion_out[[2]], rn_wolf_out[[2]],
                   rn_elk_out[[2]], rn_lago_out[[2]], rn_moose_out[[2]], rn_wtd_out[[2]]) %>%
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr21") %>%
    relocate(season, .after = Setup)
  rn_2022 <- rbind(rn_bear_out[[3]], rn_bob_out[[3]], rn_coy_out[[3]], rn_lion_out[[3]], rn_wolf_out[[3]],
                   rn_elk_out[[3]], rn_lago_out[[3]], rn_moose_out[[3]], rn_wtd_out[[3]]) %>%
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr22") %>%
    relocate(season, .after = Setup)
 
  RN_abundance <- list(rn_2020, rn_2021, rn_2022)
  save(RN_abundance, file = "./Outputs/Relative_Abundance/RN_model/RN_abundance.RData")
  

  #####  Visualize local abundance data  #####
  #'  -----------------------------------
  #'  Map relative density data per species, study area and year
  library(sf)
  library(ggplot2)
  library(patchwork)
  
  #'  Load spatial data
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  cams_20s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_20s_wgs84.shp")
  cams_21s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_21s_wgs84.shp")
  cams_22s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_22s_wgs84.shp")
  bigwater <- st_read("./Shapefiles/National Hydrology Database Idaho State/Idaho_waterbodies_1km2.shp")
  bigwater_wgs84 <- st_transform(bigwater, crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
    filter(gnis_name == "Priest Lake" | gnis_name == "Upper Priest Lake" | 
             gnis_name == "Lake Pend Oreille" | #gnis_name == "Cabinet Gorge Reservoir" | 
             gnis_name == "Chatcolet Lake" | gnis_name == "Dworshak Reservoir") %>%
    dplyr::select(gnis_name)
  priestrivers <- st_read("./Shapefiles/National Hydrology Database Idaho State/PriestRiver_flowline.shp")
  pendoreille <- st_read("./Shapefiles/National Hydrology Database Idaho State/Pendoreille_flowline.shp")
  kootenairiver <- st_read("./Shapefiles/National Hydrology Database Idaho State/KootenaiRiver_flowline.shp")
  clarkfork <- st_read("./Shapefiles/National Hydrology Database Idaho State/ClarkForkRiver_flowline.shp") 
  clearwater <- st_read("./Shapefiles/National Hydrology Database Idaho State/NorthForkClearwater_flowline.shp")
  rivers <- bind_rows(priestrivers, pendoreille, kootenairiver, clarkfork, clearwater) 
  rivers_clip <- st_intersection(rivers, eoe_gmu_wgs84)
  
  #'  List camera spatial data
  cam_list <- list(cams_20s_wgs84, cams_21s_wgs84, cams_22s_wgs84)
  
  #'  Append RN abundance estimates to spatial data
  spatial_rn <- function(rn, spp, cams) {
    #'  Filter data to single species
    single_spp_rn <- rn %>%
      filter(Species == spp) %>%
      #'  Rename camera location column to match spatial data
      rename("NwLctID" = "NewLocationID") %>%
      mutate(RN.n.rounded = round(RN.n, 0))
    
    #'  Join spatial data with rn data
    rn_shp <- full_join(cams, single_spp_rn, by = "NwLctID") %>%
      filter(!is.na(Species)) %>%
      #'  Change camera location column back to something less awkward
      rename("NewLocationID" = "NwLctID")
    
    return(rn_shp)
  }
  spatial_rn_bear <- mapply(rn = rn_bear_out, spatial_rn, spp = "bear_black", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_bob <- mapply(rn = rn_bob_out, spatial_rn, spp = "bobcat", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_coy <- mapply(rn = rn_coy_out, spatial_rn, spp = "coyote", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_lion <- mapply(rn = rn_lion_out, spatial_rn, spp = "mountain_lion", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_wolf <- mapply(rn = rn_wolf_out, spatial_rn, spp = "wolf", cams = cam_list, SIMPLIFY = FALSE)
  
  spatial_rn_elk <- mapply(rn = rn_elk_out, spatial_rn, spp = "elk", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_lago <- mapply(rn = rn_lago_out, spatial_rn, spp = "lagomorphs", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_moose <- mapply(rn = rn_moose_out, spatial_rn, spp = "moose", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_wtd <- mapply(rn = rn_wtd_out, spatial_rn, spp = "whitetailed_deer", cams = cam_list, SIMPLIFY = FALSE)
  
  #'  List spatial RN abundance data and save
  spatial_RN_list <- list(spatial_rn_bear, spatial_rn_bob, spatial_rn_coy, spatial_rn_lion, spatial_rn_wolf,
                          spatial_rn_elk, spatial_rn_lago, spatial_rn_moose, spatial_rn_wtd)
  save(spatial_RN_list, file = "./Shapefiles/IDFG spatial data/Camera_locations/spatial_RN_list.RData")
  
  #'  Merge wolf results into single large spatial dataframe and save as a shapefile
  spatial_rn_wolf_locs <- bind_rows(spatial_rn_wolf[[1]], spatial_rn_wolf[[2]], spatial_rn_wolf[[3]])
  st_write(spatial_rn_wolf_locs, "./Shapefiles/IDFG spatial data/Camera_locations/spatial_rn_wolf_locs.shp")
  
  year_list <- list("2020", "2021", "2022")
  
  #'  Function to map RN abundance
  map_rn <- function(sf_rn, yr, spp) {
    #'  Filter spatial RN data by study area
    sf_rn_gmu1 <- sf_rn[sf_rn$GMU == "GMU1",]
    sf_rn_gmu6 <- sf_rn[sf_rn$GMU == "GMU6",]
    sf_rn_gmu10a <- sf_rn[sf_rn$GMU == "GMU10A",]
    
    size_breaks <- c(0, 1, 2, 3, 5, 7, 9, 12)
    
    #'  GMU 1 plot
    gmu1_rn <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "1",], fill = NA) +
      geom_sf(data = sf_rn_gmu1, aes(size = RN.n.rounded), shape  = 21, 
              col = "darkred", fill = "darkred", alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,12)) +
      labs(size = "Estimated \nlocal abundance", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  GMU 6 plot
    gmu6_rn <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "6",], fill = NA) +
      geom_sf(data = sf_rn_gmu6, aes(size = RN.n.rounded), shape  = 21, 
              col = "darkgreen", fill = "darkgreen", alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,12)) +
      labs(size = "Estimated \nlocal abundance", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  GMU 10A plot
    gmu10a_rn <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "10A",], fill = NA) +
      geom_sf(data = sf_rn_gmu10a, aes(size = RN.n.rounded), shape = 21, 
              col = "darkblue", fill = "darkblue", alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,12)) +
      labs(size = "Estimated \nlocal abundance", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle(yr)
    
    #'  Plot each map
    print(gmu1_rn); print(gmu6_rn); print(gmu10a_rn)
    
    #'  List GMU local abundance estimate maps together
    gmu_maps <- list(gmu1_rn, gmu6_rn, gmu10a_rn)
    
    return(gmu_maps)
  }
  rn_maps_bear <- mapply(map_rn, sf_rn = spatial_rn_bear, yr = year_list, spp = "Black bear", SIMPLIFY = FALSE)
  rn_maps_bob <- mapply(map_rn, sf_rn = spatial_rn_bob, yr = year_list, spp = "Bobcat", SIMPLIFY = FALSE)
  rn_maps_coy <- mapply(map_rn, sf_rn = spatial_rn_coy, yr = year_list, spp = "Coyote", SIMPLIFY = FALSE)
  rn_maps_lion <- mapply(map_rn, sf_rn = spatial_rn_lion, yr = year_list, spp = "Mountain lion", SIMPLIFY = FALSE)
  rn_maps_wolf <- mapply(map_rn, sf_rn = spatial_rn_wolf, yr = year_list, spp = "wolf", SIMPLIFY = FALSE)
  
  rn_maps_elk <- mapply(map_rn, sf_rn = spatial_rn_elk, yr = year_list, spp = "Elk", SIMPLIFY = FALSE)
  rn_maps_lago <- mapply(map_rn, sf_rn = spatial_rn_lago, yr = year_list, spp = "Lagomorphs", SIMPLIFY = FALSE)
  rn_maps_moose <- mapply(map_rn, sf_rn = spatial_rn_moose, yr = year_list, spp = "Moose", SIMPLIFY = FALSE)
  rn_maps_wtd <- mapply(map_rn, sf_rn = spatial_rn_wtd, yr = year_list, spp = "White-tailed deer", SIMPLIFY = FALSE)
  
  #'  Plot same species and study area back to back across years
  gmu_by_yr_plots <- function(fig, spp) {
    #'  Note: list order is [[i]][[j]] 
    #'  where i = 2020, 2021, or 2022 and j = GMU1, GMU6, or GMU10a 
    #'  GMU 1 plots
    gmu1_patch <- fig[[1]][[1]] + fig[[2]][[1]] + fig[[3]][[1]] +
      plot_annotation(paste("GMU 1", spp, "relative local abundance (RN model)"))
    
    #'  GMU 6 plots
    gmu6_patch <- fig[[1]][[2]] + fig[[2]][[2]] + fig[[3]][[2]] +
      plot_annotation(paste("GMU 6", spp, "relative local abundance (RN model)"))
    
    #'  GMU 10A plots
    gmu10a_patch <- fig[[1]][[3]] + fig[[2]][[3]] + fig[[3]][[3]] +
      plot_annotation(paste("GMU 10A", spp, "relative local abundance (RN model)"))
    
    #'  Print figure panels
    print(gmu1_patch); print(gmu6_patch); print(gmu10a_patch)
    
    #'  List
    gmu_patchwork_list <- list(gmu1_patch, gmu6_patch, gmu10a_patch)
    
    return(gmu_patchwork_list)
  }
  rn_gmu_bear <- gmu_by_yr_plots(rn_maps_bear, spp = "black bear")
  rn_gmu_bob <- gmu_by_yr_plots(rn_maps_bob, spp = "bobcat")
  rn_gmu_coy <- gmu_by_yr_plots(rn_maps_coy, spp = "coyote")
  rn_gmu_lion <- gmu_by_yr_plots(rn_maps_lion, spp = "mountain lion")
  rn_gmu_wolf <- gmu_by_yr_plots(rn_maps_wolf, spp = "wolf")
  
  rn_gmu_elk <- gmu_by_yr_plots(rn_maps_elk, spp = "elk")
  rn_gmu_lago <- gmu_by_yr_plots(rn_maps_lago, spp = "lagomorphs")
  rn_gmu_moose <- gmu_by_yr_plots(rn_maps_moose, spp = "moose")
  rn_gmu_wtd <- gmu_by_yr_plots(rn_maps_wtd, spp = "white-tailed deer")
  
  
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_bear.tiff", rn_gmu_bear[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_bear.tiff", rn_gmu_bear[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_bear.tiff", rn_gmu_bear[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_bob.tiff", rn_gmu_bob[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_bob.tiff", rn_gmu_bob[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_bob.tiff", rn_gmu_bob[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_coy.tiff", rn_gmu_coy[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_coy.tiff", rn_gmu_coy[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_coy.tiff", rn_gmu_coy[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_lion.tiff", rn_gmu_lion[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_lion.tiff", rn_gmu_lion[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_lion.tiff", rn_gmu_lion[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_wolf.tiff", rn_gmu_wolf[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_wolf.tiff", rn_gmu_wolf[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_wolf.tiff", rn_gmu_wolf[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_elk.tiff", rn_gmu_elk[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_elk.tiff", rn_gmu_elk[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_elk.tiff", rn_gmu_elk[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_lago.tiff", rn_gmu_lago[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_lago.tiff", rn_gmu_lago[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_lago.tiff", rn_gmu_lago[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_moose.tiff", rn_gmu_moose[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_moose.tiff", rn_gmu_moose[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_moose.tiff", rn_gmu_moose[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu1_wtd.tiff", rn_gmu_wtd[[1]],
         units = "in", width = 13, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu6_wtd.tiff", rn_gmu_wtd[[2]],
         units = "in", width = 15, height = 4, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_gmu10A_wtd.tiff", rn_gmu_wtd[[3]],
         units = "in", width = 15, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  
  
  
  
  ######  Alternative figure arrangement  ######
  #'  ------------------------------------
  #'  Add year to each dataframe and unlist into one single large dataframe per species
  add_yr <- function(dat, yr) {
    dat$Year <- yr
    return(dat)
  }
  rn_bear_all <- mapply(add_yr, dat = spatial_rn_bear, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_bob_all <- mapply(add_yr, dat = spatial_rn_bob, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_coy_all <- mapply(add_yr, dat = spatial_rn_coy, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_lion_all <- mapply(add_yr, dat = spatial_rn_lion, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wolf_all <- mapply(add_yr, dat = spatial_rn_wolf, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  rn_elk_all <- mapply(add_yr, dat = spatial_rn_elk, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_lago_all <- mapply(add_yr, dat = spatial_rn_lago, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_moose_all <- mapply(add_yr, dat = spatial_rn_moose, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wtd_all <- mapply(add_yr, dat = spatial_rn_wtd, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  #'  Make one giant faceted plot where rows represent GMU and columns represent years 
  #'  to ensure that the dot sizes are all consistent for at least a single species
  library(ggh4x)
  map_rn_v2 <- function(sf_rn, spp) {
    #'  Define size of circles
    size_breaks <- c(0, 1, 2, 3, 5, 7, 9, 12)
    
    sf_rn <- mutate(sf_rn, GMU = factor(GMU, levels = c("GMU1", "GMU6", "GMU10A")))
    pal <- c("darkcyan", "lightcoral", "darkgoldenrod3")
    
    #'  Create figure
    spp_rn <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) + 
      geom_sf(data = rivers_clip, color = "lightskyblue3") + 
      geom_sf(data = bigwater_wgs84, fill = "lightskyblue2") +
      geom_sf(data = sf_rn, aes(size = RN.n.rounded, colour = GMU, fill = GMU), shape = 21, alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,12)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      labs(size = "Predicted \nRAI", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 18)) +
      labs(title = paste("Predicted relative abundance indicies for", spp)) + 
      # facet_grid(rows = vars(Species), cols = vars(Year)) +
      # force_panelsizes(rows = 2, cols = 1, TRUE)
      facet_wrap(~Year) 
      
    #'  Plot each map
    print(spp_rn)
    
    return(spp_rn)
  }
  rn_maps_bear <- map_rn_v2(rn_bear_all, spp = "black bears")
  rn_maps_bob <- map_rn_v2(rn_bob_all, spp = "bobcats")
  rn_maps_coy <- map_rn_v2(rn_coy_all, spp = "coyotes")
  rn_maps_lion <- map_rn_v2(rn_lion_all, spp = "mountain lions")
  rn_maps_wolf <- map_rn_v2(rn_wolf_all, spp = "wolves")
  rn_maps_elk <- map_rn_v2(rn_elk_all, spp = "elk")
  rn_maps_lago <- map_rn_v2(rn_lago_all, spp = "lagomorphs")
  rn_maps_moose <- map_rn_v2(rn_moose_all, spp = "moose")
  rn_maps_wtd <- map_rn_v2(rn_wtd_all, spp = "white-tailed deer")
  
  #'  Figures with multiple species together (use facet_grid option in function)
  ungulates <- bind_rows(rn_elk_all, rn_moose_all, rn_wtd_all) %>%
    mutate(Species = ifelse(Species == "elk", "Elk", Species),
           Species = ifelse(Species == "moose", "Moose", Species),
           Species = ifelse(Species == "whitetailed_deer", "White-tailed deer", Species))
  predators_nowolf <- bind_rows(rn_bear_all, rn_coy_all, rn_lion_all) %>%
    mutate(Species = ifelse(Species == "bear_black", "Black bear", Species),
           Species = ifelse(Species == "coyote", "Coyote", Species),
           Species = ifelse(Species == "mountain_lion", "Mountain lion", Species))
  rn_maps_ungulates <- map_rn_v2(ungulates, spp = "ungulates")
  rn_maps_predators_nowolf <- map_rn_v2(predators_nowolf, spp = "predators")
  
  
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_blackbear.tiff", rn_maps_bear,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_bobcat.tiff", rn_maps_bob,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_coyote.tiff", rn_maps_coy,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_lion.tiff", rn_maps_lion,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_wolf.tiff", rn_maps_wolf,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_elk.tiff", rn_maps_elk,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_lagomorphs.tiff", rn_maps_lago,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_moose.tiff", rn_maps_moose,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_wtd.tiff", rn_maps_wtd,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_ungualtes.tiff", rn_maps_ungulates,
         units = "in", width = 20, height = 18, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_predators_nowolf.tiff", rn_maps_predators_nowolf,
         units = "in", width = 20, height = 18, dpi = 400, device = "tiff", compression = "lzw")
  

  
  
  #'  -----------------------
  ####  LIKELIHOOD APPROACH  ####
  #'  -----------------------
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
    UMF <- unmarkedFrameOccu(y = dh, siteCovs = sitecovs)#, obsCovs = srvycovs)
    
    #'  Plot detection histories
    print(plot(UMF))
    #'  Summarize data
    summary(UMF)
    
    return(UMF)
  }
  bear_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[1]][[1]], sitecovs = stations_npp20s)#, srvycovs = srvy_covs_20s)
  bear_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[1]][[1]], sitecovs = stations_npp21s)#, srvycovs = srvy_covs_21s)
  bear_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[1]][[1]], sitecovs = stations_npp22s)#, srvycovs = srvy_covs_22s)
  
  bob_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[2]][[1]], sitecovs = stations_npp20s)#, srvycovs = srvy_covs_20s)
  bob_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[2]][[1]], sitecovs = stations_npp21s)#, srvycovs = srvy_covs_21s)
  bob_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[2]][[1]], sitecovs = stations_npp22s)#, srvycovs = srvy_covs_22s)
  
  coy_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[3]][[1]], sitecovs = stations_npp20s)#, srvycovs = srvy_covs_20s)
  coy_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[3]][[1]], sitecovs = stations_npp21s)#, srvycovs = srvy_covs_21s)
  coy_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[3]][[1]], sitecovs = stations_npp22s)#, srvycovs = srvy_covs_22s)
  
  lion_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[4]][[1]], sitecovs = stations_npp20s)#, srvycovs = srvy_covs_20s)
  lion_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[4]][[1]], sitecovs = stations_npp21s)#, srvycovs = srvy_covs_21s)
  lion_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[4]][[1]], sitecovs = stations_npp22s)#, srvycovs = srvy_covs_22s)
  
  wolf_20s_umf <- umf_setup(dh = DH_npp20s_RNmod[[5]][[1]], sitecovs = stations_npp20s)#, srvycovs = srvy_covs_20s)
  wolf_21s_umf <- umf_setup(dh = DH_npp21s_RNmod[[5]][[1]], sitecovs = stations_npp21s)#, srvycovs = srvy_covs_21s)
  wolf_22s_umf <- umf_setup(dh = DH_npp22s_RNmod[[5]][[1]], sitecovs = stations_npp22s)#, srvycovs = srvy_covs_22s)

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
  (bear20s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), bear_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (bear21s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), bear_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (bear22s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), bear_22s_umf, se = TRUE, engine = "C", threads = 4))
  bear_rn_list <- list(bear20s.rn, bear21s.rn, bear22s.rn)
  
  (bob20s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), bob_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (bob21s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), bob_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (bob22s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), bob_22s_umf, se = TRUE, engine = "C", threads = 4))
  bob_rn_list <- list(bob20s.rn, bob21s.rn, bob22s.rn)
  
  (coy20s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), coy_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (coy21s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), coy_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (coy22s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), coy_22s_umf, se = TRUE, engine = "C", threads = 4))
  coy_rn_list <- list(coy20s.rn, coy21s.rn, coy22s.rn)
  
  (lion20s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), lion_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (lion21s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), lion_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (lion22s.rn <- occuRN(~ Setu + nDays ~ GMU + PercFor + Elev + I(Elev^2), lion_22s_umf, se = TRUE, engine = "C", threads = 4)) 
  lion_rn_list <- list(lion20s.rn, lion21s.rn, lion22s.rn)
  
  (wolf20s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), wolf_20s_umf, se = TRUE, engine = "C", threads = 4)) #effort
  (wolf21s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), wolf_21s_umf, se = TRUE, engine = "C", threads = 4)) 
  (wolf22s.rn <- occuRN(~ Setup ~ GMU + PercFor + Elev + I(Elev^2), wolf_22s_umf, se = TRUE, engine = "C", threads = 4))
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
  
  