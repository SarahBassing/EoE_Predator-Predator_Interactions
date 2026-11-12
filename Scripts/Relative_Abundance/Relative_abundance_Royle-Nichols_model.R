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
  #'  
  #'  2020 - 2022 detection history list order: 
  #'  [[1]] bear_black, [[2]] bobcat, [[3]] coyote, [[4]] mountain_lion, [[5]] wolf, 
  #'  [[6]] elk, [[7]] moose, [[8]] muledeer, [[9]] whitetaileddeer, [[10]] rabbit_hare
  #'  
  #'  2023 detection history list order: 
  #'  [[1]] bear_black, [[2]] coyote, [[3]] mountain_lion, [[4]] wolf, 
  #'  [[5]] elk, [[6]] moose, [[7]] whitetaileddeer
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
  library(sf)
  
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

  #'  Load 2023 detection histories prepared by M. Parsons IDFG (NOTE THE INDEXING DIFFERS FROM OTHER YEARS)
  # load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp20s_RNmod.RData")
  # load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp21s_RNmod.RData")
  # load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp22s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/Unique_detection_events_23s.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp23s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/SamplingEffort_npp2023.RData")

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
    #' #'  Grab total number of days each camera was operable during season
    #' cam_operable <- effort %>% dplyr::select(c(NewLocationID, ndays))
    #' #'  Add total number of operable camera days to camera covariate dataset
    #' cam_covs <- full_join(cam_covs, cam_operable, by = "NewLocationID")
    #'  Rename, format, and scale as needed
    formatted <- cam_covs %>%
      mutate(GMUs = ifelse(GMU == "GMU10A", 1, 2),  # GMU10A represents the intercept!
             GMUs = ifelse(GMU == "GMU1", 3, GMUs),
             Setup = ifelse(Setup == "ungulate", 1, 2)) %>%  #'  Ungulate (random) cameras represent the intercept!
      transmute(NewLocationID = as.factor(NewLocationID),
                Season = as.factor(Season),
                GMU = as.factor(GMUs),
                Setup = as.factor(Setup))#,
                # nDays = scale(ndays),
                # PercFor = scale(perc_forest),
                # Elev = scale(Elevation__10m2)) #%>%
      #' #'  Arrange by camera location -- NECESSARY TO MATCH DH's CAMERA LOCATION ORDER
      #' arrange(NewLocationID) 
    
    return(formatted)
  }
  stations_npp20s <- format_covs(cams_eoe20s, dets = DH_npp20s_RNmod[[1]]) #, effort = effort_20s_RNmod
  stations_npp21s <- format_covs(cams_eoe21s, dets = DH_npp21s_RNmod[[1]]) #, effort = effort_21s_RNmod
  stations_npp22s <- format_covs(cams_eoe22s, dets = DH_npp22s_RNmod[[1]]) #, effort = effort_22s_RNmod
  #'  2023 camera stations
  load("./Data/IDFG camera data/IDFG 2023 detection data/Detection_Histories_RNmodel/stations_smr2023.RData")
  stations_npp23s <- stations_smr2023
  # save(stations_npp23s, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/stations_npp23s.RData")
  
  
  #'  Double check things are ordered correctly!!!!
  stations_npp20s[82:90,1:4]; DH_npp20s_RNmod[[1]][82:90,1:3]; nrow(stations_npp20s); nrow(DH_npp20s_RNmod[[1]])
  stations_npp21s[82:90,1:4]; DH_npp21s_RNmod[[1]][82:90,1:3]; nrow(stations_npp21s); nrow(DH_npp21s_RNmod[[1]])
  stations_npp22s[82:90,1:4]; DH_npp22s_RNmod[[1]][82:90,1:3]; nrow(stations_npp22s); nrow(DH_npp22s_RNmod[[1]])
  stations_npp23s[82:90,1:4]; DH_npp23s_RNmod[[1]][82:90,1:3]; nrow(stations_npp23s); nrow(DH_npp23s_RNmod[[1]])
  
  #' #'  Double check percent forested habitat and elevation aren't too correlated
  #' cor(stations_npp20s$PercFor, stations_npp20s$Elev)
  #' cor(stations_npp21s$PercFor, stations_npp21s$Elev)
  #' cor(stations_npp22s$PercFor, stations_npp22s$Elev)
  
  #' #'  -----------------------------------
  #' #####  Format survey-level covariates  #####
  #' #'  -----------------------------------
  #' #'  Scale survey-level covariates (really only relevant if using 7-day sampling occasions)
  #' scale_srvy_cov <- function(srvy_covs) {
  #'   #'  Reduce to just number of days per sampling occasion camera was operating
  #'   skinny_srvy_covs <- srvy_covs %>% dplyr::select(-c(NewLocationID, ndays, nhrs))
  #'   
  #'   #'  Replace NAs with 0 - these sites truly not surveyed so effort is 0
  #'   srvy_covs <- replace(skinny_srvy_covs, is.na(skinny_srvy_covs), 0)
  #' 
  #'   #'  Find mean & standard deviation of covariates across all sites & occasions
  #'   mu <- mean(as.matrix(srvy_covs), na.rm = TRUE)
  #'   sd <- sd(as.matrix(srvy_covs), na.rm = TRUE)
  #'   
  #'   #'  Z-transform (center observations around mean & scale by 1 SD)
  #'   scaled <- ((srvy_covs - mu) / sd)
  #'   scaled <- round(scaled, 3)
  #'   
  #'   scaled <- as.matrix(scaled)
  #'   
  #'   return(scaled)
  #' }
  #' effort_list <- list(effort_20s_RNmod, effort_21s_RNmod, effort_22s_RNmod)  
  #' zeffort_RNmod <- lapply(effort_list, scale_srvy_cov)
  #' 
  #' #'  Create list of survey level covariates for unmarked (one list per year)
  #' srvy_covs_20s <- list(effort = zeffort_RNmod[[1]])
  #' srvy_covs_21s <- list(effort = zeffort_RNmod[[2]])
  #' srvy_covs_22s <- list(effort = zeffort_RNmod[[3]])
  
  #'  -----------------------
  ####  Append cluster info  ####
  #'  -----------------------
  #'  FYI: This section is a little backwards because this depends on the wolf RN 
  #'  model results estimated below. But once those are fit, clusters are identified 
  #'  with Camera_clustering.R script. Once each camera is assigned a cluster, the 
  #'  cluster info is appended to the RN model input data here and the RN models 
  #'  are run for all species (and re-run for wolves) so the expected relative 
  #'  density index can be estimated for each cluster.
  #'  ----------------------
  #'  Function to clean up camera cluster spatial data
  reformat_clusters <- function(clusters) {
    skinny_clusters <- clusters %>%
      dplyr::select(c(NwLctID, CellID, GMU, Setup, Clustrs, are_km2, geometry)) %>%
      rename(NewLocationID = NwLctID) %>%
      rename(ClusterID = Clustrs) %>%
      rename(area_km2 = are_km2)
    return(skinny_clusters)
  }
  
  #'  Read in camera cluster shapefiles and clean up with reformat_clusters() function
  clusters_gmu1 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu1_05.11.26.shp") %>% #06.30.25
    reformat_clusters(.)
  clusters_gmu6 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu6_05.11.26.shp") %>%
    reformat_clusters(.)
  clusters_gmu10a <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu10a_05.11.26.shp") %>%
    reformat_clusters(.)
  
  #'  Merge spatial camera cluster data together and reduce to one observation per camera
  clusters_all <- bind_rows(clusters_gmu1, clusters_gmu6, clusters_gmu10a) %>%
    mutate(GMU_cluster = paste0(GMU, "_", ClusterID)) %>%
    group_by(NewLocationID) %>%
    slice(1L) %>%
    ungroup() %>%
    as.data.frame(.) %>%
    dplyr::select(-c(CellID, GMU, Setup, geometry)) %>%
    group_by(GMU_cluster) %>%
    mutate(uniqueCluster = cur_group_id()) %>%
    ungroup()
  #'  How many unique clusters are there?
  max(clusters_all$uniqueCluster)
  
  #'  Join camera stations with camera cluster info
  cluster_cams <- function(stations, clusters) {
    clustered_cams <- stations %>%
      left_join(clusters, by = "NewLocationID") %>%
      #'  Count number of unique cameras per cluster
      group_by(GMU_cluster) %>%
      mutate(ncams = n()) %>%
      ungroup()
    return(clustered_cams)
  }
  stations <- list(stations_npp20s, stations_npp21s, stations_npp22s, stations_npp23s)
  clustered_cams <- lapply(stations, cluster_cams, clusters = clusters_all)
  
  #'  Create indicator matrix for each camera per cluster
  cluster_indicator_matrix <- function(clust_cams) {
    cluster_indicator <- clust_cams %>%
      #'  Focus on each camera and which unique cluster it was assigned to
      dplyr::select(c("NewLocationID", "uniqueCluster")) %>%
      #'  If row (camera) was assigned to uniqueCluster X, mark as 1, else, mark as 0
      mutate(cluster.1 = ifelse(uniqueCluster == 1, 1, 0),    #'  GMU10A
             cluster.2 = ifelse(uniqueCluster == 2, 1, 0),
             cluster.3 = ifelse(uniqueCluster == 3, 1, 0),
             cluster.4 = ifelse(uniqueCluster == 4, 1, 0),
             cluster.5 = ifelse(uniqueCluster == 5, 1, 0),
             cluster.6 = ifelse(uniqueCluster == 6, 1, 0),
             cluster.7 = ifelse(uniqueCluster == 7, 1, 0), 
             cluster.8 = ifelse(uniqueCluster == 8, 1, 0),
             cluster.9 = ifelse(uniqueCluster == 9, 1, 0),    #'  GMU1
             cluster.10 = ifelse(uniqueCluster == 10, 1, 0),
             cluster.11 = ifelse(uniqueCluster == 11, 1, 0), 
             cluster.12 = ifelse(uniqueCluster == 12, 1, 0),
             cluster.13 = ifelse(uniqueCluster == 13, 1, 0),
             cluster.14 = ifelse(uniqueCluster == 14, 1, 0),
             cluster.15 = ifelse(uniqueCluster == 15, 1, 0),
             cluster.16 = ifelse(uniqueCluster == 16, 1, 0),
             cluster.17 = ifelse(uniqueCluster == 17, 1, 0),
             cluster.18 = ifelse(uniqueCluster == 18, 1, 0),
             cluster.19 = ifelse(uniqueCluster == 19, 1, 0),   #'  GMU6
             cluster.20 = ifelse(uniqueCluster == 20, 1, 0),
             cluster.21 = ifelse(uniqueCluster == 21, 1, 0),
             cluster.22 = ifelse(uniqueCluster == 22, 1, 0),
             cluster.23 = ifelse(uniqueCluster == 23, 1, 0)) %>%#,
             # cluster.24 = ifelse(uniqueCluster == 24, 1, 0)) %>%
      dplyr::select(-c(NewLocationID, uniqueCluster)) %>%
      as.matrix(.)
    return(cluster_indicator)
  }
  cluster_index <- lapply(clustered_cams, cluster_indicator_matrix)
  head(cluster_index[[1]])
  head(cluster_index[[2]])
  head(cluster_index[[3]])
  head(cluster_index[[4]])
  
  #'  -----------------------
  ####  Royle-Nichols model  ####
  #'  -----------------------
  #'  Run Royle-Nichols model using a Bayesian statistical framework in JAGS and
  #'  a likelihood framework with unmarked
  
  #'  ------------------------
  #####  Setup data for JAGS  #####
  #'  ------------------------
  #'  Bundle detection histories and covariates for each species and year
  bundle_dat <- function(dh, nsite, nsurvey, cov, effort, cluster_indicator) {
    #'  Convert detection history to matrix
    dh <- as.matrix(dh)
    dimnames(dh) <- NULL
    #' #'  Count number of sites per GMU
    #' ncams_perGMU <- cov %>%
    #'   group_by(GMU) %>%
    #'   summarise(nsites = n()) %>%
    #'   ungroup()
    
    #'  Bundle data for JAGS
    bundled <- list(y = dh, 
                    nsites = dim(dh)[1], 
                    nsurveys = dim(dh)[2], 
                    ngmu = max(as.numeric(cov$GMU)),
                    nsets = max(as.numeric(cov$Setup)),
                    # ncams1 = as.numeric(ncams_perGMU[1,2]), # GMU10A
                    # ncams2 = as.numeric(ncams_perGMU[2,2]), # GMU6
                    # ncams3 = as.numeric(ifelse(is.na(ncams_perGMU[3,2]), 0, ncams_perGMU[3,2])), #GMU1
                    gmu = as.numeric(cov$GMU), 
                    setup = as.numeric(cov$Setup),
                    # forest = as.numeric(cov$PercFor), 
                    # elev = as.numeric(cov$Elev), 
                    # ndays = as.numeric(cov$nDays), 
                    # seffort = as.numeric(effort),
                    #'  Cluster-specific data for estimating relative density/cluster
                    ncams = as.numeric(cov$ncams),
                    cluster_area = as.numeric(cov$area_km2),
                    cluster_matrix = cluster_indicator)
              
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle_20s <- lapply(DH_npp20s_RNmod, bundle_dat, cov = clustered_cams[[1]], cluster_indicator = cluster_index[[1]]) #cov = stations_npp20s, effort = zeffort_RNmod[[1]]
  data_JAGS_bundle_21s <- lapply(DH_npp21s_RNmod, bundle_dat, cov = clustered_cams[[2]], cluster_indicator = cluster_index[[2]]) #cov = stations_npp21s, effort = zeffort_RNmod[[2]]
  data_JAGS_bundle_22s <- lapply(DH_npp22s_RNmod, bundle_dat, cov = clustered_cams[[3]], cluster_indicator = cluster_index[[3]]) #cov = stations_npp22s, effort = zeffort_RNmod[[3]]
  data_JAGS_bundle_23s <- lapply(DH_npp23s_RNmod, bundle_dat, cov = clustered_cams[[3]], cluster_indicator = cluster_index[[3]]) #cov = stations_npp23s, effort = zeffort_RNmod[[3]]
  
  #'  Save
  save(data_JAGS_bundle_20s, file = paste0("./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_20s_", Sys.Date(), ".RData")) #12.26.25
  save(data_JAGS_bundle_21s, file = paste0("./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_21s_", Sys.Date(), ".RData"))
  save(data_JAGS_bundle_22s, file = paste0("./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_22s_", Sys.Date(), ".RData"))
  save(data_JAGS_bundle_23s, file = paste0("./Data/Relative abundance data/RAI Phase 2/data_JAGS_bundle_23s_", Sys.Date(), ".RData"))
  
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
  ninit_23s <- lapply(DH_npp23s_RNmod, initial_n)
  
  #'  Parameters monitored
  params <- c("beta0", "beta4", #"beta1", "beta2", "beta3", #"mean.lambda", 
              "alpha0", "alpha2", "rSetup", "mu.r", "mean.p", #"alpha1" #"mean.r", 
              "mu.lambda", "totalN", "rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", 
              "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9", "rdi.cl10", 
              "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", "rdi.cl15", "rdi.cl16", 
              "rdi.cl17", "rdi.cl18", "rdi.cl19", "rdi.cl20", "rdi.cl21", "rdi.cl22", 
              "rdi.cl23", # "rdi.cl24", #"occSites", "mean.psi", 
              # "totalN.gmu10a", "densitykm2.gmu10a", "density100km2.gmu10a", 
              # "totalN.gmu6", "densitykm2.gmu6", "density100km2.gmu6", 
              # "totalN.gmu1", "densitykm2.gmu1", "density100km2.gmu1",
              "N")
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
  source("./Scripts/Relative_Abundance/RNmodel_clusters_JAGS_code_2020.R")
  
  ######  Black bear  ######
  start.time = Sys.time()
  inits_bear20s <- function(){list(N = ninit_20s[[1]])}
  RN_bear_20s <- jags(data_JAGS_bundle_20s[[1]], inits = inits_bear20s, params,
                      # "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_20s$summary)
  which(RN_bear_20s$summary[,"Rhat"] < 0.9)
  which(RN_bear_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_20s$samples)
  save(RN_bear_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_bear_20s_", Sys.Date(), ".RData"))

  # ######  Bobcat  ######
  # start.time = Sys.time()
  # inits_bob20s <- function(){list(N = ninit_20s[[2]])}
  # RN_bob_20s <- jags(data_JAGS_bundle_20s[[2]], inits = inits_bob20s, params,
  #                    "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
  #                    n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
  #                    n.burnin = nb, parallel = TRUE)
  # end.time <- Sys.time(); (run.time <- end.time - start.time)
  # print(RN_bob_20s$summary)
  # which(RN_bob_20s$summary[,"Rhat"] > 1.1)
  # mcmcplot(RN_bob_20s$samples)
  # save(RN_bob_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_bob_20s_", Sys.Date(), ".RData")) 
  
  ######  Coyote  ######
  ni_coy <-  100000
  start.time = Sys.time()
  inits_coy20s <- function(){list(N = ninit_20s[[3]])}
  RN_coy_20s <- jags(data_JAGS_bundle_20s[[3]], inits = inits_coy20s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_coy, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_20s$summary)
  which(RN_coy_20s$summary[,"Rhat"] < 0.9)
  which(RN_coy_20s$summary[,"Rhat"] > 1.1)   # a few r-hat >1.1
  mcmcplot(RN_coy_20s$samples)
  save(RN_coy_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_coy_20s_", Sys.Date(), ".RData"))
  
  ######  Mountain lion  ######
  start.time = Sys.time()
  inits_lion20s <- function(){list(N = ninit_20s[[4]])}
  RN_lion_20s <- jags(data_JAGS_bundle_20s[[4]], inits = inits_lion20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_20s$summary)
  which(RN_lion_20s$summary[,"Rhat"] < 0.9)
  which(RN_lion_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_20s$samples)
  save(RN_lion_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_lion_20s_", Sys.Date(), ".RData")) 
  
  ######  Wolf  ######
  ni_wolf <-  100000 
  start.time = Sys.time()
  inits_wolf20s <- function(){list(N = ninit_20s[[5]])}
  RN_wolf_20s <- jags(data_JAGS_bundle_20s[[5]], inits = inits_wolf20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_20s$summary)
  which(RN_wolf_20s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_20s$samples)
  save(RN_wolf_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wolf_20s_", Sys.Date(), ".RData")) 
  
  ######  Elk  ######
  start.time = Sys.time()
  inits_elk20s <- function(){list(N = ninit_20s[[6]])}
  RN_elk_20s <- jags(data_JAGS_bundle_20s[[6]], inits = inits_elk20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_20s$summary)
  which(RN_elk_20s$summary[,"Rhat"] < 0.9)
  which(RN_elk_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_20s$samples)
  save(RN_elk_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_elk_20s_", Sys.Date(), ".RData"))
  
  ######  Moose  ######
  start.time = Sys.time()
  inits_moose20s <- function(){list(N = ninit_20s[[7]])}
  RN_moose_20s <- jags(data_JAGS_bundle_20s[[7]], inits = inits_moose20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_20s$summary)
  which(RN_moose_20s$summary[,"Rhat"] < 0.9)
  which(RN_moose_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_20s$samples)
  save(RN_moose_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_moose_20s_", Sys.Date(), ".RData"))
  
  ######  White-tailed Deer  ######
  ni_wtd <- 100000
  start.time = Sys.time()
  inits_wtd20s <- function(){list(N = ninit_20s[[9]])}
  RN_wtd_20s <- jags(data_JAGS_bundle_20s[[9]], inits = inits_wtd20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_20s$summary)
  which(RN_wtd_20s$summary[,"Rhat"] < 0.9)
  which(RN_wtd_20s$summary[,"Rhat"] > 1.1) #' several r-hat >1.1
  mcmcplot(RN_wtd_20s$samples)
  save(RN_wtd_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wtd_20s_", Sys.Date(), ".RData"))
  
  # ######  Lagomorphs  ######
  # ni_lag <- 100000
  # start.time = Sys.time()
  # inits_lago20s <- function(){list(N = ninit_20s[[10]])}
  # RN_lago_20s <- jags(data_JAGS_bundle_20s[[10]], inits = inits_lago20s, params,
  #                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters_2020.txt",
  #                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_lag, 
  #                     n.burnin = nb, parallel = TRUE)
  # end.time <- Sys.time(); (run.time <- end.time - start.time)
  # print(RN_lago_20s$summary)
  # which(RN_lago_20s$summary[,"Rhat"] < 0.9)
  # which(RN_lago_20s$summary[,"Rhat"] > 1.1) 
  # mcmcplot(RN_lago_20s$samples)
  # save(RN_lago_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_lagomorphs_20s_", Sys.Date(), ".RData"))
  
  
  ######  2021 - 2023  Analyses  ######
  source("./Scripts/Relative_Abundance/RNmodel_clusters_JAGS_code.R")
  
  ######  Black bear  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_bear21s <- function(){list(N = ninit_21s[[1]])}
  RN_bear_21s <- jags(data_JAGS_bundle_21s[[1]], inits = inits_bear21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_21s$summary)
  which(RN_bear_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_21s$samples)
  save(RN_bear_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_bear_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_bear22s <- function(){list(N = ninit_22s[[1]])}
  RN_bear_22s <- jags(data_JAGS_bundle_22s[[1]], inits = inits_bear22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_22s$summary)
  which(RN_bear_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_22s$samples)
  save(RN_bear_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_bear_22s_", Sys.Date(), ".RData"))
  
  #'  Summer 2023
  start.time = Sys.time()
  inits_bear23s <- function(){list(N = ninit_23s[[1]])}  # NOTE same indexing for 2023 bears
  RN_bear_23s <- jags(data_JAGS_bundle_23s[[1]], inits = inits_bear23s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_23s$summary)
  which(RN_bear_23s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_23s$samples)
  save(RN_bear_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_bear_23s_", Sys.Date(), ".RData"))
  
  
  #' ######  Bobcat  ######
  #' #'  Summer 2021
  #' start.time = Sys.time()
  #' inits_bob21s <- function(){list(N = ninit_21s[[2]])}
  #' RN_bob_21s <- jags(data_JAGS_bundle_21s[[2]], inits = inits_bob21s, params,
  #'                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
  #'                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
  #'                     n.burnin = nb, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(RN_bob_21s$summary)
  #' which(RN_bob_21s$summary[,"Rhat"] > 1.1)
  #' mcmcplot(RN_bob_21s$samples)
  #' save(RN_bob_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_bob_21s_", Sys.Date(), ".RData"))
  #' 
  #' #'  Summer 2022
  #' start.time = Sys.time()
  #' inits_bob22s <- function(){list(N = ninit_22s[[2]])}
  #' RN_bob_22s <- jags(data_JAGS_bundle_22s[[2]], inits = inits_bob22s, params,
  #'                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
  #'                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
  #'                     n.burnin = nb, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(RN_bob_22s$summary)
  #' which(RN_bob_22s$summary[,"Rhat"] > 1.1)
  #' mcmcplot(RN_bob_22s$samples)
  #' save(RN_bob_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_bob_22s_", Sys.Date(), ".RData"))

  
  ######  Coyote  ######
  ni_coy <-  100000
  #'  Summer 2021
  start.time = Sys.time()
  inits_coy21s <- function(){list(N = ninit_21s[[3]])}
  RN_coy_21s <- jags(data_JAGS_bundle_21s[[3]], inits = inits_coy21s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_coy, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_21s$summary)
  which(RN_coy_21s$summary[,"Rhat"] < 0.9)
  which(RN_coy_21s$summary[,"Rhat"] > 1.1)  # a few N had r>1.1
  mcmcplot(RN_coy_21s$samples)
  save(RN_coy_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_coy_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_coy22s <- function(){list(N = ninit_22s[[3]])}
  RN_coy_22s <- jags(data_JAGS_bundle_22s[[3]], inits = inits_coy22s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_coy, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_22s$summary)
  which(RN_coy_22s$summary[,"Rhat"] < 0.9)
  which(RN_coy_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_22s$samples)
  save(RN_coy_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_coy_22s_", Sys.Date(), ".RData"))
  
  #'  Summer 2023
  start.time = Sys.time()
  inits_coy23s <- function(){list(N = ninit_23s[[2]])}  #  NOTE the different indexing
  RN_coy_23s <- jags(data_JAGS_bundle_23s[[2]], inits = inits_coy23s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_coy, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_23s$summary)
  which(RN_coy_23s$summary[,"Rhat"] < 0.9)
  which(RN_coy_23s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_23s$samples)
  save(RN_coy_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_coy_23s_", Sys.Date(), ".RData"))
  
  
  ######  Mountain lion  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_lion21s <- function(){list(N = ninit_21s[[4]])}
  RN_lion_21s <- jags(data_JAGS_bundle_21s[[4]], inits = inits_lion21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_21s$summary)
  which(RN_lion_21s$summary[,"Rhat"] < 0.9)
  which(RN_lion_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_21s$samples)
  save(RN_lion_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_lion_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_lion22s <- function(){list(N = ninit_22s[[4]])}
  RN_lion_22s <- jags(data_JAGS_bundle_22s[[4]], inits = inits_lion22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_22s$summary)
  which(RN_lion_22s$summary[,"Rhat"] < 0.9)
  which(RN_lion_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_22s$samples)
  save(RN_lion_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_lion_22s_", Sys.Date(), ".RData"))
  
  #'  Summer 2023
  start.time = Sys.time()
  inits_lion23s <- function(){list(N = ninit_23s[[3]])}  #NOTE the different indexing
  RN_lion_23s <- jags(data_JAGS_bundle_23s[[3]], inits = inits_lion23s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_23s$summary)
  which(RN_lion_23s$summary[,"Rhat"] < 0.9)
  which(RN_lion_23s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_23s$samples)
  save(RN_lion_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_lion_23s_", Sys.Date(), ".RData"))
  
  
  ######  Wolf  ######
  ni_wolf <- 100000 
  #'  Summer 2021
  start.time = Sys.time()
  inits_wolf21s <- function(){list(N = ninit_21s[[5]])}
  RN_wolf_21s <- jags(data_JAGS_bundle_21s[[5]], inits = inits_wolf21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_21s$summary)
  which(RN_wolf_21s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_21s$samples)
  save(RN_wolf_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wolf_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_wolf22s <- function(){list(N = ninit_22s[[5]])}
  RN_wolf_22s <- jags(data_JAGS_bundle_22s[[5]], inits = inits_wolf22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_22s$summary)
  which(RN_wolf_22s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_22s$samples)
  save(RN_wolf_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wolf_22s_", Sys.Date(), ".RData"))

  #'  Summer 2023
  start.time = Sys.time()
  inits_wolf23s <- function(){list(N = ninit_23s[[4]])}  # NOTE the different indexing
  RN_wolf_23s <- jags(data_JAGS_bundle_23s[[4]], inits = inits_wolf23s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_23s$summary)
  which(RN_wolf_23s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_23s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_23s$samples)
  save(RN_wolf_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wolf_23s_", Sys.Date(), ".RData"))
  
    
  ######  Elk  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_elk21s <- function(){list(N = ninit_21s[[6]])}
  RN_elk_21s <- jags(data_JAGS_bundle_21s[[6]], inits = inits_elk21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_21s$summary)
  which(RN_elk_21s$summary[,"Rhat"] < 0.9)  
  which(RN_elk_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_21s$samples)
  save(RN_elk_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_elk_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_elk22s <- function(){list(N = ninit_22s[[6]])}
  RN_elk_22s <- jags(data_JAGS_bundle_22s[[6]], inits = inits_elk22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_22s$summary)
  which(RN_elk_22s$summary[,"Rhat"] < 0.9)  
  which(RN_elk_22s$summary[,"Rhat"] > 1.1)  # a few N had r>1.1
  mcmcplot(RN_elk_22s$samples)
  save(RN_elk_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_elk_22s_", Sys.Date(), ".RData"))
  
  #'  Summer 2023
  start.time = Sys.time()
  inits_elk23s <- function(){list(N = ninit_23s[[5]])}  # NOTE the different indexing
  RN_elk_23s <- jags(data_JAGS_bundle_23s[[5]], inits = inits_elk23s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_23s$summary)
  which(RN_elk_23s$summary[,"Rhat"] < 0.9)  
  which(RN_elk_23s$summary[,"Rhat"] > 1.1)  
  mcmcplot(RN_elk_23s$samples)
  save(RN_elk_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_elk_23s_", Sys.Date(), ".RData"))
  
  
  ######  Moose  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_moose21s <- function(){list(N = ninit_21s[[7]])}
  RN_moose_21s <- jags(data_JAGS_bundle_21s[[7]], inits = inits_moose21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_21s$summary)
  which(RN_moose_21s$summary[,"Rhat"] < 0.9)  
  which(RN_moose_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_21s$samples)
  save(RN_moose_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_moose_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_moose22s <- function(){list(N = ninit_22s[[7]])}
  RN_moose_22s <- jags(data_JAGS_bundle_22s[[7]], inits = inits_moose22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_22s$summary)
  which(RN_moose_22s$summary[,"Rhat"] < 0.9)  
  which(RN_moose_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_22s$samples)
  save(RN_moose_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_moose_22s_", Sys.Date(), ".RData"))
  
  #'  Summer 2023
  start.time = Sys.time()
  inits_moose23s <- function(){list(N = ninit_23s[[6]])}  # NOTE the different indexing
  RN_moose_23s <- jags(data_JAGS_bundle_23s[[6]], inits = inits_moose23s, params,
                       "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                       n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                       n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_23s$summary)
  which(RN_moose_23s$summary[,"Rhat"] < 0.9)  
  which(RN_moose_23s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_23s$samples)
  save(RN_moose_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_moose_23s_", Sys.Date(), ".RData"))
  
  
  ######  White-tailed Deer  ######
  ni_wtd <- 100000
  #'  Summer 2021
  start.time = Sys.time()
  inits_wtd21s <- function(){list(N = ninit_21s[[9]])}
  RN_wtd_21s <- jags(data_JAGS_bundle_21s[[9]], inits = inits_wtd21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_21s$summary)
  which(RN_wtd_21s$summary[,"Rhat"] < 0.9) 
  which(RN_wtd_21s$summary[,"Rhat"] > 1.1) # a handful of N had r>1.1
  mcmcplot(RN_wtd_21s$samples)
  save(RN_wtd_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wtd_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_wtd22s <- function(){list(N = ninit_22s[[9]])}
  RN_wtd_22s <- jags(data_JAGS_bundle_22s[[9]], inits = inits_wtd22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_22s$summary)
  which(RN_wtd_22s$summary[,"Rhat"] < 0.9) 
  which(RN_wtd_22s$summary[,"Rhat"] > 1.1) # a good handful of N had r>1.1
  mcmcplot(RN_wtd_22s$samples)
  save(RN_wtd_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wtd_22s_", Sys.Date(), ".RData"))
  
  #'  Summer 2023
  start.time = Sys.time()
  inits_wtd23s <- function(){list(N = ninit_23s[[7]])}
  RN_wtd_23s <- jags(data_JAGS_bundle_23s[[7]], inits = inits_wtd23s, params,
                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_23s$summary)
  which(RN_wtd_23s$summary[,"Rhat"] < 0.9) 
  which(RN_wtd_23s$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_wtd_23s$samples)
  save(RN_wtd_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_wtd_23s_", Sys.Date(), ".RData"))
  
  
  #' ######  Lagomorphs  ######
  #' ni_lag <- 100000
  #' #'  Summer 2021
  #' start.time = Sys.time()
  #' inits_lago21s <- function(){list(N = ninit_21s[[10]])}
  #' RN_lago_21s <- jags(data_JAGS_bundle_21s[[10]], inits = inits_lago21s, params,
  #'                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
  #'                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_lag, 
  #'                     n.burnin = nb, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(RN_lago_21s$summary)
  #' which(RN_lago_21s$summary[,"Rhat"] > 1.1)
  #' mcmcplot(RN_lago_21s$samples)
  #' save(RN_lago_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_lagomorphs_21s_", Sys.Date(), ".RData"))
  #' 
  #' #'  Summer 2022
  #' start.time = Sys.time()
  #' inits_lago22s <- function(){list(N = ninit_22s[[10]])}
  #' RN_lago_22s <- jags(data_JAGS_bundle_22s[[10]], inits = inits_lago22s, params,
  #'                     "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_clusters.txt",
  #'                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_lag, 
  #'                     n.burnin = nb, parallel = TRUE)
  #' end.time <- Sys.time(); (run.time <- end.time - start.time)
  #' print(RN_lago_22s$summary)
  #' which(RN_lago_22s$summary[,"Rhat"] > 1.1)
  #' mcmcplot(RN_lago_22s$samples)
  #' save(RN_lago_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_clusterRDI_lagomorphs_22s_", Sys.Date(), ".RData"))
  
  
  #'  -------------------------------
  ####  Summarize density estimates  ####
  #'  -------------------------------
  #'  Load model outputs
  filenames <- list.files("./Outputs/Relative_Abundance/RN_model/JAGS_out", pattern="*.RData", full.names=TRUE)
  lapply(filenames, load, environment())
  
  #'  List model outputs together 
  rn_bear_list <- list(RN_bear_20s, RN_bear_21s, RN_bear_22s, RN_bear_23s)
  rn_coy_list <- list(RN_coy_20s, RN_coy_21s, RN_coy_22s, RN_coy_23s)
  rn_lion_list <- list(RN_lion_20s, RN_lion_21s, RN_lion_22s, RN_lion_23s)
  rn_wolf_list <- list(RN_wolf_20s, RN_wolf_21s, RN_wolf_22s, RN_wolf_23s)
  
  rn_elk_list <- list(RN_elk_20s, RN_elk_21s, RN_elk_22s, RN_elk_23s)
  rn_moose_list <- list(RN_moose_20s, RN_moose_21s, RN_moose_22s, RN_moose_23s)
  rn_wtd_list <- list(RN_wtd_20s, RN_wtd_21s, RN_wtd_22s, RN_wtd_23s)
  
  #'  Load saved detection histories
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp20s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp21s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp22s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp23s_RNmod.RData")
  
  #'  List one detection history per year (doesn't matter which species it relates to)
  dh_list <- list(DH_npp20s_RNmod[[1]], DH_npp21s_RNmod[[1]], DH_npp22s_RNmod[[1]], DH_npp23s_RNmod[[1]])
  
  #'  Load camera spatial data
  cams_20s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_20s_wgs84.shp")
  cams_21s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_21s_wgs84.shp")
  cams_22s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_22s_wgs84.shp")
  cams_23s_wgs84 <- st_read("./Shapefiles/IDFG spatial data 2023 additional/cams_23s_wgs84.shp")
  
  #'  List camera spatial data
  cam_list <- list(cams_20s_wgs84, cams_21s_wgs84, cams_22s_wgs84, cams_23s_wgs84)
  
  #####  Relative Abundance Index (N --> RAI)  #####
  #'  -----------------------------------------
  #'  Save estimated N per site
  estimated_N <- function(mod, dh, spp) {
    #'  Grab estimated N and SD per camera site
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
  # rn_bob_out <- mapply(estimated_N, rn_bob_list, dh = dh_list, spp = "bobcat", SIMPLIFY = FALSE)
  rn_coy_out <- mapply(estimated_N, rn_coy_list, dh = dh_list, spp = "coyote", SIMPLIFY = FALSE)
  rn_lion_out <- mapply(estimated_N, rn_lion_list, dh = dh_list, spp = "mountain_lion", SIMPLIFY = FALSE)
  rn_wolf_out <- mapply(estimated_N, rn_wolf_list, dh = dh_list, spp = "wolf", SIMPLIFY = FALSE)
  
  rn_elk_out <- mapply(estimated_N, rn_elk_list, dh = dh_list, spp = "elk", SIMPLIFY = FALSE)
  # rn_lago_out <- mapply(estimated_N, rn_lago_list, dh = dh_list, spp = "lagomorphs", SIMPLIFY = FALSE)
  rn_moose_out <- mapply(estimated_N, rn_moose_list, dh = dh_list, spp = "moose", SIMPLIFY = FALSE)
  rn_wtd_out <- mapply(estimated_N, rn_wtd_list, dh = dh_list, spp = "whitetailed_deer", SIMPLIFY = FALSE)
  
  rn_2020 <- rbind(rn_bear_out[[1]], rn_coy_out[[1]], rn_lion_out[[1]], rn_wolf_out[[1]],
                   rn_elk_out[[1]], rn_moose_out[[1]], rn_wtd_out[[1]]) %>% #rn_bob_out[[1]], rn_lago_out[[1]], 
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr20") %>%
    relocate(season, .after = Setup)
  rn_2021 <- rbind(rn_bear_out[[2]], rn_coy_out[[2]], rn_lion_out[[2]], rn_wolf_out[[2]],
                   rn_elk_out[[2]], rn_moose_out[[2]], rn_wtd_out[[2]]) %>% #rn_bob_out[[2]], rn_lago_out[[2]], 
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr21") %>%
    relocate(season, .after = Setup)
  rn_2022 <- rbind(rn_bear_out[[3]], rn_coy_out[[3]], rn_lion_out[[3]], rn_wolf_out[[3]],
                   rn_elk_out[[3]], rn_moose_out[[3]], rn_wtd_out[[3]]) %>% # rn_bob_out[[3]], rn_lago_out[[3]], 
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr22") %>%
    relocate(season, .after = Setup)
  rn_2023 <- rbind(rn_bear_out[[4]], rn_coy_out[[4]], rn_lion_out[[4]], rn_wolf_out[[4]],
                   rn_elk_out[[4]], rn_moose_out[[4]], rn_wtd_out[[4]]) %>% 
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr22") %>%
    relocate(season, .after = Setup)
 
  RN_abundance <- list(rn_2020, rn_2021, rn_2022, rn_2023)
  save(RN_abundance, file = paste0("./Outputs/Relative_Abundance/RN_model/RN_abundance_", Sys.Date(), ".RData")) #_2026-10-27
  
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
  spatial_rn_coy <- mapply(rn = rn_coy_out, spatial_rn, spp = "coyote", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_lion <- mapply(rn = rn_lion_out, spatial_rn, spp = "mountain_lion", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_wolf <- mapply(rn = rn_wolf_out, spatial_rn, spp = "wolf", cams = cam_list, SIMPLIFY = FALSE)
  
  spatial_rn_elk <- mapply(rn = rn_elk_out, spatial_rn, spp = "elk", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_moose <- mapply(rn = rn_moose_out, spatial_rn, spp = "moose", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_wtd <- mapply(rn = rn_wtd_out, spatial_rn, spp = "whitetailed_deer", cams = cam_list, SIMPLIFY = FALSE)
  
  #'  List spatial RN abundance data and save
  spatial_RN_list <- list(spatial_rn_bear, spatial_rn_coy, spatial_rn_lion, spatial_rn_wolf,
                          spatial_rn_elk, spatial_rn_moose, spatial_rn_wtd) 
  save(spatial_RN_list, file = paste0("./Shapefiles/IDFG spatial data/Camera_locations/spatial_RN_list_2020_2023_", Sys.Date(), ".RData"))
  
  #'  Merge wolf results into single large spatial dataframe and save as a shapefile
  spatial_rn_wolf_locs <- bind_rows(spatial_rn_wolf[[1]], spatial_rn_wolf[[2]], spatial_rn_wolf[[3]], spatial_rn_wolf[[4]])
  st_write(spatial_rn_wolf_locs, "./Shapefiles/IDFG spatial data/Camera_locations/spatial_rn_wolf_locs_2020_2023.shp")
  

  #####  Relative Density Index (RDI)  #####
  #'  ---------------------------------
  #'  Load cluster polygons
  gmu1_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu1_05.11.26.shp") %>% 
    mutate(GMU = "GMU1",
           Cluster = 9:18) %>%
    rename("GMU_clusterID" = "Clusters")
  gmu6_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu6_05.11.26.shp") %>% 
    mutate(GMU = "GMU6",
           Cluster = 19:23) %>%
    rename("GMU_clusterID" = "Clusters")
  gmu10a_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu10a_05.11.26.shp") %>% 
    mutate(GMU = "GMU10A",
           Cluster = 1:8) %>%
    rename("GMU_clusterID" = "Clusters")
  
  #'  Full set of cluster polygons
  gmu_poly <- bind_rows(gmu10a_poly, gmu1_poly, gmu6_poly) %>%
    relocate(GMU_clusterID, .after = Cluster) %>%
    relocate(area_km2, .after = GMU_clusterID) %>%
    mutate(Cluster = as.numeric(Cluster))
  
  #'  Save estimated RDI per cluster
  estimated_RDI <- function(mod, spp) {
    #'  Grab JAGS output names that start with "rdi.cl"
    rdi_names <- grep("^rdi\\.cl[0-9]+$", names(mod$mean), value = TRUE)
    
    #'  Grab estimated RDI per cluster (GMU 10A: cl1 - cl8; GMU 1: cl9 - cl18; GMU 6: cl19 - cl23)
    rdi_mean <- unlist(mod$mean[rdi_names]) 
    rdi_sd <- unlist(mod$sd[rdi_names]) 
    rdi_ll <- unlist(mod$q2.5[rdi_names]) 
    rdi_ul <- unlist(mod$q97.5[rdi_names]) 
    
    #'  Create data frame with mean and 95% CRI for each cluster's estimated RDI
    rdi <- data.frame(rdi = round(rdi_mean, 2), 
                      rdi_sd = round(rdi_sd, 2),
                      rdi_lower = round(rdi_ll, 2), 
                      rdi_upper = round(rdi_ul, 2)) %>%
      #'  Add cluster and species info
      mutate(Cluster = as.numeric(sub("^rdi\\.cl", "", row.names(.))),
             Species = spp) %>%
      relocate(Cluster, .before = rdi) %>%
      relocate(Species, .before = Cluster)
    
    #'  Grab cluster shapefile
    cluster_sf <- gmu_poly %>% dplyr::select(-GMU_clusterID) %>%
      mutate(Cluster = as.numeric(Cluster))
    
    #'  Join spatial data with RDI estimates
    cluster_RDI <- right_join(cluster_sf, rdi, by = "Cluster")
    
    return(cluster_RDI)
  }
  spatial_rdi_bear <- mapply(estimated_RDI, mod = rn_bear_list, spp = "bear_black", SIMPLIFY = FALSE)
  spatial_rdi_coy <- mapply(estimated_RDI, mod = rn_coy_list, spp = "coyote", SIMPLIFY = FALSE)
  spatial_rdi_lion <- mapply(estimated_RDI, mod = rn_lion_list, spp = "mountain_lion", SIMPLIFY = FALSE)
  spatial_rdi_wolf <- mapply(estimated_RDI, mod = rn_wolf_list, spp = "wolf", SIMPLIFY = FALSE)
  spatial_rdi_elk <- mapply(estimated_RDI, mod = rn_elk_list, spp = "elk", SIMPLIFY = FALSE)
  spatial_rdi_moose <- mapply(estimated_RDI, mod = rn_moose_list, spp = "moose", SIMPLIFY = FALSE)
  spatial_rdi_wtd <- mapply(estimated_RDI, mod = rn_wtd_list, spp = "whitetailed_deer", SIMPLIFY = FALSE)
  
  #'  List spatial RN abundance data and save
  spatial_rdi_list <- list(spatial_rdi_bear, spatial_rdi_coy, spatial_rdi_lion, spatial_rdi_wolf,
                          spatial_rdi_elk, spatial_rdi_moose, spatial_rdi_wtd) 
  save(spatial_rdi_list, file = "./Shapefiles/IDFG spatial data/Camera_locations/spatial_RDI_list_2020_2023.RData")
  
  #'  Merge wolf results into single large spatial data frame and save as a shapefile
  spatial_RDI_wolf_locs <- bind_rows(spatial_rdi_wolf[[1]], spatial_rdi_wolf[[2]], spatial_rdi_wolf[[3]], spatial_rdi_wolf[[4]])
  st_write(spatial_RDI_wolf_locs, "./Shapefiles/IDFG spatial data/Camera_locations/spatial_rdi_wolf_locs_2020_2023.shp")
  
  
  #'  Drop spatial attributes and create summary table of RDI input data for SEM
  #'  List years
  year <- c("2020", "2021", "2022", "2023")
  #'  Function to reformat data 
  reformat_rdi <- function(spatial_rdi, spp) {
    #'  Name list by year
    names(spatial_rdi) <- year
    #'  Bind annual lists into single data frame
    nonspatial_rdi_df <- bind_rows(spatial_rdi, .id = "Year") %>%
      as.data.frame(.) %>%
      #'  Remove spatial data
      dplyr::select(-geometry) %>%
      #'  Reformat columns
      transmute(GMU = GMU, 
                Year = Year, 
                Species = spp,
                Cluster = Cluster,
                Area_km2 = round(area_km2, 2),
                #'  Combine mean and SD from cluster-level posterior distribution
                RDI = paste0(rdi, " (", rdi_sd, ")"))
    return(nonspatial_rdi_df)
  }
  rdi_bear_clean <- reformat_rdi(spatial_rdi_bear, spp = "Black bear")
  rdi_coy_clean <- reformat_rdi(spatial_rdi_coy, spp = "Coyote")
  rdi_lion_clean <- reformat_rdi(spatial_rdi_lion, spp = "Mountain lion")
  rdi_wolf_clean <- reformat_rdi(spatial_rdi_wolf, spp = "Wolf")
  rdi_elk_clean <- reformat_rdi(spatial_rdi_elk, spp = "Elk")
  rdi_moose_clean <- reformat_rdi(spatial_rdi_moose, spp = "Moose")
  rdi_wtd_clean <- reformat_rdi(spatial_rdi_wtd, spp = "White-tailed deer")
  
  #'  Bind all together now
  rdi_master_input <- bind_rows(rdi_bear_clean, rdi_coy_clean, rdi_lion_clean,
                                rdi_wolf_clean, rdi_elk_clean, rdi_moose_clean, 
                                rdi_wtd_clean) %>%
    pivot_wider(id_cols = c(GMU, Year, Cluster, Area_km2), names_from = Species, values_from = RDI) %>%
    #'  Arrange by GMU based on N to S arrangement
    arrange(match(GMU, c("GMU1", "GMU6", "GMU10A")))

  write_csv(rdi_master_input, paste0("./Outputs/Relative_Abundance/RN_model/Tables/RDI_master_input_", Sys.Date(), ".csv"))
  
    
  #' #'  Grab estimated RDI and SD per cluster
  #' estimated_N <- function(mod, dh, spp) {
  #'   #'  Grab estimated RDI and SD per cluster
  #'   RDI.cl1 <- mod$mean$rdi.cl1        # GMU 10A
  #'   RDI.cl2 <- mod$mean$rdi.cl2
  #'   RDI.cl3 <- mod$mean$rdi.cl3
  #'   RDI.cl4 <- mod$mean$rdi.cl4
  #'   RDI.cl5 <- mod$mean$rdi.cl5
  #'   RDI.cl6 <- mod$mean$rdi.cl6
  #'   RDI.cl7 <- mod$mean$rdi.cl7
  #'   RDI.cl8 <- mod$mean$rdi.cl8
  #'   RDI.cl9 <- mod$mean$rdi.cl9        # GMU 1
  #'   RDI.cl10 <- mod$mean$rdi.cl10
  #'   RDI.cl11 <- mod$mean$rdi.cl11
  #'   RDI.cl12 <- mod$mean$rdi.cl12
  #'   RDI.cl13 <- mod$mean$rdi.cl13
  #'   RDI.cl14 <- mod$mean$rdi.cl14
  #'   RDI.cl15 <- mod$mean$rdi.cl15
  #'   RDI.cl16 <- mod$mean$rdi.cl16
  #'   RDI.cl17 <- mod$mean$rdi.cl17
  #'   RDI.cl18 <- mod$mean$rdi.cl18
  #'   RDI.cl19 <- mod$mean$rdi.cl19      # GMU 6
  #'   RDI.cl20 <- mod$mean$rdi.cl20
  #'   RDI.cl22 <- mod$mean$rdi.cl21
  #'   RDI.cl23 <- mod$mean$rdi.cl22
  #'   # RDI.cl24 <- mod$mean$rdi.cl23
  #' 
  #'   #'  Merge and format into single data frame with corresponding N & SD per site
  #'   out <- cbind(locs, RN.n, RN.sd)
  #'   RN_est <- as.data.frame(out) %>%
  #'     mutate(NewLocationID = locs,
  #'            Species = spp,
  #'            RN.n = as.numeric(RN.n),
  #'            RN.sd = as.numeric(RN.sd)) %>%
  #'     separate(locs, c("GMU", "Setup", "site"), sep = "_") %>%
  #'     mutate(CellID = paste0(GMU, "_", site)) %>%
  #'     dplyr::select(-site) %>%
  #'     relocate(NewLocationID, .before = "GMU") %>%
  #'     relocate(Species, .after = "Setup") %>%
  #'     relocate(CellID, .after = "NewLocationID")
  #'   return(RN_est)
  #' }
  #' rn_bear_out <- mapply(estimated_N, rn_bear_list, dh = dh_list, spp = "bear_black", SIMPLIFY = FALSE)
  #' 
  #' #'  Append RN abundance estimates to spatial data
  #' spatial_rn <- function(rn, spp, cams) {
  #'   #'  Filter data to single species
  #'   single_spp_rn <- rn %>%
  #'     filter(Species == spp) %>%
  #'     #'  Rename camera location column to match spatial data
  #'     rename("NwLctID" = "NewLocationID") %>%
  #'     mutate(RN.n.rounded = round(RN.n, 0))
  #' 
  #'   #'  Join spatial data with rn data
  #'   rn_shp <- full_join(cams, single_spp_rn, by = "NwLctID") %>%
  #'     filter(!is.na(Species)) %>%
  #'     #'  Change camera location column back to something less awkward
  #'     rename("NewLocationID" = "NwLctID")
  #' 
  #'   return(rn_shp)
  #' }
  #' spatial_rn_bear <- mapply(rn = rn_bear_out, spatial_rn, spp = "bear_black", cams = cam_list, SIMPLIFY = FALSE)
  #' spatial_rn_coy <- mapply(rn = rn_coy_out, spatial_rn, spp = "coyote", cams = cam_list, SIMPLIFY = FALSE)
  #' spatial_rn_lion <- mapply(rn = rn_lion_out, spatial_rn, spp = "mountain_lion", cams = cam_list, SIMPLIFY = FALSE)
  #' spatial_rn_wolf <- mapply(rn = rn_wolf_out, spatial_rn, spp = "wolf", cams = cam_list, SIMPLIFY = FALSE)
  #' 
  #' spatial_rn_elk <- mapply(rn = rn_elk_out, spatial_rn, spp = "elk", cams = cam_list, SIMPLIFY = FALSE)
  #' spatial_rn_moose <- mapply(rn = rn_moose_out, spatial_rn, spp = "moose", cams = cam_list, SIMPLIFY = FALSE)
  #' spatial_rn_wtd <- mapply(rn = rn_wtd_out, spatial_rn, spp = "whitetailed_deer", cams = cam_list, SIMPLIFY = FALSE)
  #' 
  #' #'  List spatial RN abundance data and save
  #' spatial_RN_list <- list(spatial_rn_bear, spatial_rn_coy, spatial_rn_lion, spatial_rn_wolf,
  #'                         spatial_rn_elk, spatial_rn_moose, spatial_rn_wtd) #spatial_rn_bob, spatial_rn_lago,
  #' save(spatial_RN_list, file = "./Shapefiles/IDFG spatial data/Camera_locations/spatial_RN_list_2020_2023.RData")
  #' 
  #' #'  Merge wolf results into single large spatial dataframe and save as a shapefile
  #' spatial_rn_wolf_locs <- bind_rows(spatial_rn_wolf[[1]], spatial_rn_wolf[[2]], spatial_rn_wolf[[3]], spatial_rn_wolf[[4]])
  #' st_write(spatial_rn_wolf_locs, "./Shapefiles/IDFG spatial data/Camera_locations/spatial_rn_wolf_locs_2020_2023.shp")
  #' 

  #####  Visualize RN model outputs  #####
  #'  -------------------------------
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
  cams_23s_wgs84 <- st_read("./Shapefiles/IDFG spatial data 2023 additional/cams_23s_wgs84.shp")
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

  #'  Join survey year with each GMU (note GMU1 was not surveyed in 2020)
  year_surveyed <- rep(c(2020, 2021, 2022, 2023), 3)
  gmu_by_year <- bind_rows(eoe_gmu_wgs84, eoe_gmu_wgs84, eoe_gmu_wgs84, eoe_gmu_wgs84) %>%
    mutate(Year = year_surveyed) %>%
    filter(NAME != "1" | Year != 2020)
  year_surveyed <- rep(c(2020, 2021, 2022, 2023), 384)
  river_by_year <- bind_rows(rivers_clip, rivers_clip, rivers_clip, rivers_clip) %>%
    mutate(Year = year_surveyed) %>%
    filter(Year != 2020)
  year_surveyed <- rep(c(2020, 2021, 2022, 2023), 6)
  bigwater_by_year <- bind_rows(bigwater_wgs84, bigwater_wgs84, bigwater_wgs84, bigwater_wgs84) %>%
    mutate(Year = year_surveyed) %>%
    filter(gnis_name != "Lake Pend Oreille" | Year != 2020) %>%
    filter(gnis_name != "Upper Priest Lake" | Year != 2020) %>%
    filter(gnis_name != "Priest Lake" | Year != 2020)  
  
  ######  Site-specific Relative Abundance  ######
  #'  --------------------------------------
  year_list <- list("2020", "2021", "2022", "2023")
  
  #'  Add year to each data frame and unlist into one single large data frame per species
  add_yr <- function(dat, yr) {
    dat$Year <- yr
    return(dat)
  }
  rn_bear_all <- mapply(add_yr, dat = spatial_rn_bear, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  # rn_bob_all <- mapply(add_yr, dat = spatial_rn_bob, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_coy_all <- mapply(add_yr, dat = spatial_rn_coy, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_lion_all <- mapply(add_yr, dat = spatial_rn_lion, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wolf_all <- mapply(add_yr, dat = spatial_rn_wolf, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  rn_elk_all <- mapply(add_yr, dat = spatial_rn_elk, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  # rn_lago_all <- mapply(add_yr, dat = spatial_rn_lago, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_moose_all <- mapply(add_yr, dat = spatial_rn_moose, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wtd_all <- mapply(add_yr, dat = spatial_rn_wtd, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  

  #'  Make one giant faceted plot where rows represent GMU and columns represent years 
  #'  to ensure that the dot sizes are all consistent for at least a single species
  library(ggh4x)
  map_rn <- function(sf_rn, spp) {
    #'  Define size of circles
    size_breaks <- c(0, 1, 2, 3, 5, 7, 9, 12)
    
    sf_rn <- mutate(sf_rn, GMU = factor(GMU, levels = c("GMU1", "GMU6", "GMU10A")))
    pal <- c("darkcyan", "lightcoral", "darkgoldenrod3")
    
    #'  Create figure
    spp_rn <- ggplot() +
      geom_sf(data = gmu_by_year, fill = NA) + 
      geom_sf(data = river_by_year, color = "lightskyblue3") + 
      geom_sf(data = bigwater_by_year, fill = "lightskyblue2") +
      geom_sf(data = sf_rn, aes(size = RN.n.rounded, colour = GMU, fill = GMU), shape = 21, alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,12)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      labs(size = "Predicted \nRAI", x = "Longitude", y = "Latitude", 
           title = paste(spp, "relative abundance index")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 22)) +
      #'  Legends: opaque, larger symbols (alpha override only affects the legend)
      guides(
        size = guide_legend(order = 1, ncol = 1, 
                            override.aes = list(alpha = 1, fill = "grey70",
                                                colour = "grey20", stroke = 1)),
        colour = guide_legend(order = 2, ncol = 1, override.aes = list(alpha = 1, size = 7, stroke = 1.2)),
        fill = "none"      # colour and fill share one palette, so show only one GMU legend
      ) +
      # facet_grid(rows = vars(Species), cols = vars(Year)) +
      # force_panelsizes(rows = 2, cols = 1, TRUE)
      facet_wrap(~Year, nrow = 1, ncol = 4) 
      
    #'  Plot each map
    print(spp_rn)
    
    return(spp_rn)
  }
  rn_maps_bear <- map_rn(rn_bear_all, spp = "Black bear")
  rn_maps_coy <- map_rn(rn_coy_all, spp = "Coyote")
  rn_maps_lion <- map_rn(rn_lion_all, spp = "Mountain lion")
  rn_maps_wolf <- map_rn(rn_wolf_all, spp = "Wolf")
  rn_maps_elk <- map_rn(rn_elk_all, spp = "Elk")
  rn_maps_moose <- map_rn(rn_moose_all, spp = "Moose")
  rn_maps_wtd <- map_rn(rn_wtd_all, spp = "White-tailed deer")
  
  
  ######  Cluster-specific Relative Density  ######
  #'  ---------------------------------------
  #'  Add year to RDI data and combine into single data frame per species
  rdi_bear_all <- mapply(add_yr, dat = spatial_rdi_bear, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rdi_coy_all <- mapply(add_yr, dat = spatial_rdi_coy, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rdi_lion_all <- mapply(add_yr, dat = spatial_rdi_lion, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rdi_wolf_all <- mapply(add_yr, dat = spatial_rdi_wolf, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rdi_elk_all <- mapply(add_yr, dat = spatial_rdi_elk, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rdi_moose_all <- mapply(add_yr, dat = spatial_rdi_moose, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rdi_wtd_all <- mapply(add_yr, dat = spatial_rdi_wtd, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  #'  Function to map RDI per cluster for each GMU and year
  map_rdi <- function(sf_rdi, n_breaks = 5, spp) {
    #'  Define min and max bounds or data set
    rdi_min <- floor(min(sf_rdi$rdi, na.rm = TRUE) * 100) / 100
    rdi_max <- ceiling(max(sf_rdi$rdi, na.rm = TRUE) * 100) / 100
    
    #'  Define breaks in continuous fill scale 
    size_breaks <- round(seq(rdi_min, rdi_max, length.out = n_breaks), 2)
    print(size_breaks)
    
    #'  Create figure
    spp_rdi <- ggplot() +
      geom_sf(data = sf_rdi, aes(fill = rdi, group = Year)) +
      scale_fill_gradient2(low = "#f5f5f5", mid = "#54a1a1", high = "#1f6f6f", 
                           midpoint = mean(sf_rdi$rdi),
                           limits = c(rdi_min, rdi_max),
                           breaks = size_breaks,
                           labels = function(x) sprintf("%.2f", x)) + 
      geom_sf(data = gmu_by_year, fill = NA) + 
      geom_sf(data = river_by_year, color = "lightskyblue3") + 
      geom_sf(data = bigwater_by_year, fill = "lightskyblue2") +
      labs(fill = "Relative \nDensity \nIndex (RDI)", x = "Longitude", y = "Latitude",
           title = paste(spp, "relative density index")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 22)) +
      facet_wrap(~Year, nrow = 1, ncol = 4) 
    
    return(spp_rdi)
  }
  (rdi_maps_bear <- map_rdi(sf_rdi = rdi_bear_all, spp = "Black bear"))
  (rdi_maps_coy <- map_rdi(sf_rdi = rdi_coy_all, spp = "Coyote"))
  (rdi_maps_lion <- map_rdi(sf_rdi = rdi_lion_all, spp = "Mountain lion"))
  (rdi_maps_wolf <- map_rdi(sf_rdi = rdi_wolf_all, spp = "Wolf"))
  (rdi_maps_elk <- map_rdi(sf_rdi = rdi_elk_all, spp = "Elk"))
  (rdi_maps_moose <- map_rdi(sf_rdi = rdi_moose_all, spp = "Moose"))
  (rdi_maps_wtd <- map_rdi(sf_rdi = rdi_wtd_all, spp = "White-tailed deer"))
  
  
  #####  Merge RAI and RDI maps  #####
  #'  ---------------------------
  (bear_RN_maps <- rn_maps_bear + rdi_maps_bear + plot_annotation(tag_levels = 'A'))
  (coy_RN_maps <- rn_maps_coy + rdi_maps_coy + plot_annotation(tag_levels = 'A'))
  (lion_RN_maps <- rn_maps_lion + rdi_maps_lion + plot_annotation(tag_levels = 'A'))
  (wolf_RN_maps <- rn_maps_wolf + rdi_maps_wolf + plot_annotation(tag_levels = 'A'))
  (elk_RN_maps <- rn_maps_elk + rdi_maps_elk + plot_annotation(tag_levels = 'A'))
  (moose_RN_maps <- rn_maps_moose + rdi_maps_moose + plot_annotation(tag_levels = 'A'))
  (wtd_RN_maps <- rn_maps_wtd + rdi_maps_wtd + plot_annotation(tag_levels = 'A'))
  
  
  #####  Save figures  #####
  #'  -----------------
  #'  Stand alone relative abundance index (RAI)
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_blackbear_2020_2023.tiff", rn_maps_bear,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_coyote_2020_2023.tiff", rn_maps_coy,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_lion_2020_2023.tiff", rn_maps_lion,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_wolf_2020_2023.tiff", rn_maps_wolf,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_elk_2020_2023.tiff", rn_maps_elk,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_moose_2020_2023.tiff", rn_maps_moose,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RN_map_wtd_2020_2023.tiff", rn_maps_wtd,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  
  #'  Stand alone relative density index (RDI)
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_map_blackbear_2020_2023.tiff", rdi_maps_bear,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_map_coyote_2020_2023.tiff", rdi_maps_coy,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_map_lion_2020_2023.tiff", rdi_maps_lion,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_map_wolf_2020_2023.tiff", rdi_maps_wolf,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_map_elk_2020_2023.tiff", rdi_maps_elk,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_map_moose_2020_2023.tiff", rdi_maps_moose,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_map_wtd_2020_2023.tiff", rdi_maps_wtd,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  
  #'  Merged figures (RAI and RDI)
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RAI_RDI_map_blackbear_2020_2023.tiff", bear_RN_maps,
         units = "in", width = 30, height = 10, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RAI_RDI_map_coyote_2020_2023.tiff", coy_RN_maps,
         units = "in", width = 30, height = 10, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RAI_RDI_map_lion_2020_2023.tiff", lion_RN_maps,
         units = "in", width = 30, height = 10, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RAI_RDI_map_wolf_2020_2023.tiff", wolf_RN_maps,
         units = "in", width = 30, height = 10, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RAI_RDI_map_elk_2020_2023.tiff", elk_RN_maps,
         units = "in", width = 30, height = 10, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RAI_RDI_map_moose_2020_2023.tiff", moose_RN_maps,
         units = "in", width = 30, height = 10, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RAI_RDI_map_wtd_2020_2023.tiff", wtd_RN_maps,
         units = "in", width = 30, height = 10, dpi = 600, device = "tiff", compression = "lzw")
  
  
  
  #'  Next step is reformatting posterior RDI estimates for SEM with
  #'  Format_RNmodel_Posterior_for_SEM.R script
  
  