  #'  ---------------------------
  #'  Local abundance of deer and elk in N. Idaho
  #'  Northern Idaho Predator-Prey Project
  #'  Sarah B. Bassing
  #'  August 2024
  #'  ---------------------------
  #'  Script to format photo-capture data of elk and white-tailed deer in northern 
  #'  Idaho and fit Royle-Nichols abundance models to detection histories to test
  #'  the effects of forage on ungulate relative abundance. Fit RN models per
  #'  month, with annual data stacked together.
  #'  
  #'  Input data: Species-specific detection data (June, July, Aug, 2020, 2021, & 2022)
  #'  from cameras and camera site data regarding camera setup and basic habitat variables. 
  #'  ---------------------------
  
  #'  Clear workspace
  rm(list = ls())
  
  #'  Load libraries
  library(camtrapR)
  library(jagsUI)
  library(chron)
  library(sf)
  library(terra)
  library(tidyverse)
  library(mcmcplots)
  
  #'  Load detection data
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe20s_allM_NewLocationID_2023-09-26.RData") 
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe21s_allM_NewLocationID_2023-09-26.RData")
  load("./Data/IDFG camera data/Split datasets/Updated_EoE_datasets/eoe22s_allM_NewLocationID_2023-09-26.RData")
  
  #'  Load sequential problem images
  load("./Data/IDFG camera data/Problem images/eoe20s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe21s_sequential_probimgs.RData")
  load("./Data/IDFG camera data/Problem images/eoe22s_sequential_probimgs.RData")
  
  #'  Load problem cameras
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe22s_problem_cams.RData")
  
  cams_2020_2022 <- list(eoe_probcams_20s[,1:3], eoe_probcams_21s[,1:3], eoe_probcams_22s[,1:3])
  
  #'  Load basic camera site covariate data
  cams_eoe_long <- read.csv("./Data/side_quests/Hilger/cams_eoe_long_Smr2020-2022.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Season")) %>%
    distinct()
  
  #'  Load forage data for elk
  file_paths_elk <- list.files(path = "./Data/side_quests/Hilger/Ungulate density covariates kg_ha/Elk", pattern = "\\.rds", full.names = TRUE)
  file_names_elk <- gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths_elk))
  forage_list_elk <- lapply(file_paths_elk, readRDS)
  names(forage_list_elk) <- c("Elk2020_kg_ha", "Elk2021_kg_ha", "Elk2022_kg_ha")
  
  #'  Load forage data for white-tailed deer
  file_paths_wtd <- list.files(path = "./Data/side_quests/Hilger/Ungulate density covariates kg_ha/WTD", pattern = "\\.rds", full.names = TRUE)
  file_names_wtd <- gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths_wtd))
  forage_list_wtd <- lapply(file_paths_wtd, readRDS)
  names(forage_list_wtd) <- c("WTD2020_kg_ha", "WTD2021_kg_ha", "WTD2022_kg_ha")
  head(forage_list_wtd[[1]])
  
  #'  Read in spatial data and extract
  pforest_100m <- rast("./Shapefiles/National Land Cover Database (NCLD)/PercentForest_100m.tif")
  elev <- rast("./Shapefiles/IDFG spatial data/Elevation__10m2.tif")
  #'  Check projections and define them
  crs(pforest_100m, describe = TRUE, proj = TRUE)
  crs(elev, describe = TRUE, proj = TRUE)
  wgs84 <- crs("+proj=longlat +datum=WGS84 +no_defs")
  (aea <- crs(pforest_100m, proj = TRUE))
  (nad83 <- crs(elev, proj = TRUE))
  
  #'  Make camera location data spatial sf objects
  spatial_locs <- function(locs, proj) {
    locs <- arrange(locs, NewLocationID)
    sf_locs <- st_as_sf(locs, coords = c("Long", "Lat"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      dplyr::select(c("NewLocationID")) %>%
      st_transform(proj)
    return(sf_locs)
  }
  cams_aea <- lapply(cams_2020_2022, spatial_locs, proj = aea)
  cams_nad83 <- lapply(cams_2020_2022, spatial_locs, proj = nad83)
  
  extract_covs <- function(locs_aea, locs_nad83) {
    #'  Extract percent forest cover within 100m of each camera
    perc_forest <- terra::extract(pforest_100m, vect(locs_aea)) %>%
      transmute(ID = ID, perc_forest = focal_sum) %>%
      mutate(perc_forest = round(perc_forest, 2))
    #'  Extract elevation (m) at each camera
    elev <- terra::extract(elev, vect(locs_nad83))
    
    #'  Combine into a single data frame
    covs <- as.data.frame(locs_aea) %>%
      mutate(ID = seq(1:nrow(.))) %>%
      full_join(perc_forest, by = "ID") %>%
      full_join(elev, by = "ID") %>%
      rename("elev" = "Elevation__10m2") %>%
      dplyr::select(-geometry)
    
    return(covs)
  }
  sitecovs_20s <- extract_covs(cams_aea[[1]], cams_nad83[[1]])
  sitecovs_21s <- extract_covs(cams_aea[[2]], cams_nad83[[2]])
  sitecovs_22s <- extract_covs(cams_aea[[3]], cams_nad83[[3]])
  
  #'  Merge sitecovs and forage data per species and year
  join_covs <- function(sitecovs, numnumcovs) {
    covs <- full_join(sitecovs, numnumcovs, by = c("NewLocationID" = "original_id"))
    return(covs)
  }
  covs_elk_20s <- join_covs(sitecovs_20s, forage_list_elk[[1]])
  covs_elk_21s <- join_covs(sitecovs_21s, forage_list_elk[[2]])
  covs_elk_22s <- join_covs(sitecovs_22s, forage_list_elk[[3]])
  covs_wtd_20s <- join_covs(sitecovs_20s, forage_list_wtd[[1]])
  covs_wtd_21s <- join_covs(sitecovs_21s, forage_list_wtd[[2]])
  covs_wtd_22s <- join_covs(sitecovs_22s, forage_list_wtd[[3]])
  
  #'  --------------------------
  ####  Filter and format data  ####
  #'  --------------------------
  #####  1) Filter detection data to time period of interest and remove problem images  #####
  #'  ----------------------------------------------------------------------------------
  thin_detections <- function(dets, seqprobs, start_date, end_date) {
    #'  Remove images with known problems (e.g., viewshed obscured)
    skinny_dets <- dets[!(dets$File %in% seqprobs$File),]
    
    #'  Filter images to focal species and time period of interest
    clean_dets <- skinny_dets %>%
      filter(Species == "elk" | Species == "whitetaileddeer") %>%
      #'  Add count = 1 for species missing count data 
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
      dplyr::select(c("NewLocationID", "posix_date_time", "Date", "Time", "Species", "Count")) 
    
    return(clean_dets)
  }
  df_all_20s <- thin_detections(eoe20s_allM, seqprobs = eoe_seqprob_20s, start_date = "2020-07-01", end_date = "2020-08-31") 
  df_all_21s <- thin_detections(eoe21s_allM, seqprobs = eoe_seqprob_21s, start_date = "2021-07-01", end_date = "2021-08-31")
  df_all_22s <- thin_detections(eoe22s_allM, seqprobs = eoe_seqprob_22s, start_date = "2022-07-01", end_date = "2022-08-31") 
  
  #####  2) Generate independent detection events  #####
  #'  ---------------------------------------------
  #'  Defined as independent when >= 5 min have elapsed between sequential images 
  #'  of the same species (5*60 = 300 seconds) or a new species is detected
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
  eoe20s_det_events <- unique_detections(df_all_20s, elapsed_time = 300)
  eoe21s_det_events <- unique_detections(df_all_21s, elapsed_time = 300)
  eoe22s_det_events <- unique_detections(df_all_22s, elapsed_time = 300)
  
  #####  3) Create a matrix with each camera & dates deployed  #####
  #'  ---------------------------------------------------------
  #'  1 = operating; 0 = not operating but deployed; NA = not deployed
  camera_operation_tbl <- function(cams) {
    cams <- arrange(cams, NewLocationID)
    #'  Add one day to retrieval date when setup & retrieval dates are the same
    #'  Necessary for cameraOperation function below... so annoying
    same_startend <- cams[cams$Setup_date == cams$Retrieval_date,]
    same_startend <- mutate(same_startend, Retrieval_date = as.Date(Retrieval_date) + 1)
    same_startend <- mutate(same_startend, Problem1_to = as.Date(Problem1_to) + 1)
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
  eoe20s_probs <- camera_operation_tbl(eoe_probcams_20s) 
  eoe21s_probs <- camera_operation_tbl(eoe_probcams_21s)
  eoe22s_probs <- camera_operation_tbl(eoe_probcams_22s)  
  
  #####  4) Create species-specific detection histories  #####
  #'  ---------------------------------------------------
  #'  FYI: July 1 - Aug 31 = 62 1-day sampling periods
  DH <- function(dets, cam_probs, spp, start_date, y, oc) {
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
    
    #'  Remove rows where camera was inoperable during the entire sampling season
    short_dh <- as.data.frame(short_dh)
    short_effort <- as.data.frame(short_effort)  %>%
      mutate(NewLocationID = row.names(.)) %>% relocate(NewLocationID, .before = "o1")
    short_dh <- filter(short_dh, rowSums(is.na(short_dh)) != ncol(short_dh)) %>%
      mutate(NewLocationID = row.names(.)) %>% relocate(NewLocationID, .before = "o1")
    short_effort <- short_effort[short_effort$NewLocationID %in% short_dh$NewLocationID,]
    
    #'  Remove rows where cameras was inoperable for >30 days (we require at least 
    #'  31 days of active surveying for site to be included, i.e., <31 NAs)
    short_dh <- filter(short_dh, rowSums(is.na(short_dh)) < 31)
    short_effort <- short_effort[short_effort$NewLocationID %in% short_dh$NewLocationID,]
    
    #'  Drop NewLocationID column
    short_dh <- dplyr::select(short_dh, -NewLocationID)
    short_effort <- dplyr::select(short_effort, -NewLocationID)
    
    dh_list <- list(short_dh, short_effort)
    
    return(dh_list)
  }
  #'  Create season-specific detection histories for species listed below
  spp_smr <- list("elk", "whitetaileddeer")
  DHeff_eoe20s_RNmod <- lapply(spp_smr, DH, dets = eoe20s_det_events, cam_probs = eoe20s_probs, start_date = "2020-07-01", y = "binary", oc = 62) 
  DHeff_eoe21s_RNmod <- lapply(spp_smr, DH, dets = eoe21s_det_events, cam_probs = eoe21s_probs, start_date = "2021-07-01", y = "binary", oc = 62)
  DHeff_eoe22s_RNmod <- lapply(spp_smr, DH, dets = eoe22s_det_events, cam_probs = eoe22s_probs, start_date = "2022-07-01", y = "binary", oc = 62)
  
  #'  Remove sampling effort from species-specific DH_npp lists
  strip_list <- function(dh) {
    #'  Keep only the detection history per species
    dh_only <- dh[[1]]
    return(dh_only)
  }
  DH_eoe20s_RNmod <- lapply(DHeff_eoe20s_RNmod, strip_list)
  DH_eoe21s_RNmod <- lapply(DHeff_eoe21s_RNmod, strip_list)
  DH_eoe22s_RNmod <- lapply(DHeff_eoe22s_RNmod, strip_list) 
  
  #'  Split by month
  split_DH <- function(dh) {
    july <- dh[,1:31]
    august <- dh[,32:62]
    month_list <- list(july, august)
    return(month_list)
  }
  #'  Monthly list order: Elk[[1]], White-tailed Deer[[2]]; July[[1]], August[[2]]
  DH_eoe20s_monthly <- lapply(DH_eoe20s_RNmod, split_DH)
  DH_eoe21s_monthly <- lapply(DH_eoe21s_RNmod, split_DH)
  DH_eoe22s_monthly <- lapply(DH_eoe22s_RNmod, split_DH)
  
  #'  Stack annual detection histories per month and species
  #'  e.g., all July elk detections (2020, 2021, 2022) stacked into a single detection history 
  DH_elk_jul <- rbind(DH_eoe20s_monthly[[1]][[1]], DH_eoe21s_monthly[[1]][[1]], DH_eoe22s_monthly[[1]][[1]])
  DH_elk_aug <- rbind(DH_eoe20s_monthly[[1]][[2]], DH_eoe21s_monthly[[1]][[2]], DH_eoe22s_monthly[[1]][[2]])
  DH_wtd_jul <- rbind(DH_eoe20s_monthly[[2]][[1]], DH_eoe21s_monthly[[2]][[1]], DH_eoe22s_monthly[[2]][[1]])
  DH_wtd_aug <- rbind(DH_eoe20s_monthly[[2]][[2]], DH_eoe21s_monthly[[2]][[2]], DH_eoe22s_monthly[[2]][[2]])
  
  #'  List for easier use
  DH_elk_list <- list(DH_elk_jul, DH_elk_aug)
  DH_wtd_list <- list(DH_wtd_jul, DH_wtd_aug)
  
  #####  5) Format covariate data  #####
  #'  -----------------------------
  #'  Break up covariate data by summer
  sites_20s <- filter(cams_eoe_long, Season == "Smr20")
  sites_21s <- filter(cams_eoe_long, Season == "Smr21")
  sites_22s <- filter(cams_eoe_long, Season == "Smr22") %>%
    distinct()
  
  #'  Scale and format site-level covariates
  format_covs <- function(dets, cams, sitecovs) {   
    #'  Remove rows where camera was inoperable the entire season - covariates at 
    #'  these sites shouldn't be included when scaling since they don't contribute
    #'  to detection data
    camsites <- row.names(dets)
    cam_covs <- cams[(cams$NewLocationID %in% camsites),]
    
    #'  Rename, format, and scale as needed
    formatted <- cam_covs %>%
      mutate(GMUs = ifelse(Gmu == "10A", 1, 2),  # GMU10A represents the intercept!
             GMUs = ifelse(Gmu == "1", 3, GMUs),
             Setup = ifelse(Setup == "ungulate", 1, 2)) %>%  #'  Ungulate (random) cameras represent the intercept!
      transmute(NewLocationID = as.factor(NewLocationID),
                Season = as.factor(Season),
                GMU = as.factor(GMUs),
                Setup = as.factor(Setup)) %>%
      #'  Join with site-level covariates
      left_join(sitecovs, by = join_by(NewLocationID)) %>%
      #'  Arrange by camera location -- NECESSARY TO MATCH DH's CAMERA LOCATION ORDER
      arrange(factor(NewLocationID, levels = camsites))
    
    return(formatted)
  }
  stations_elk_eoe20s <- format_covs(DH_eoe20s_RNmod[[1]], cams = sites_20s, sitecovs = covs_elk_20s) %>% dplyr::select(-ID)
  stations_elk_eoe21s <- format_covs(DH_eoe21s_RNmod[[1]], cams = sites_21s, sitecovs = covs_elk_21s) %>% dplyr::select(-ID)
  stations_elk_eoe22s <- format_covs(DH_eoe22s_RNmod[[1]], cams = sites_22s, sitecovs = covs_elk_22s) %>% dplyr::select(-ID)
  stations_wtd_eoe20s <- format_covs(DH_eoe20s_RNmod[[1]], cams = sites_20s, sitecovs = covs_wtd_20s) %>% dplyr::select(-ID)
  stations_wtd_eoe21s <- format_covs(DH_eoe21s_RNmod[[1]], cams = sites_21s, sitecovs = covs_wtd_21s) %>% dplyr::select(-ID)
  stations_wtd_eoe22s <- format_covs(DH_eoe22s_RNmod[[1]], cams = sites_22s, sitecovs = covs_wtd_22s) %>% dplyr::select(-ID)
  
  #'  Stack together 
  station_stack_elk <- rbind(stations_elk_eoe20s, stations_elk_eoe21s, stations_elk_eoe22s)
  station_stack_wtd <- rbind(stations_wtd_eoe20s, stations_wtd_eoe21s, stations_wtd_eoe22s)
  
  #'  Seasonal stacked data
  station_stack_elk_july <- station_stack_elk %>% dplyr::select(-c("mean_Tbio_august_kg_ha", "max_Tbio_august_kg_ha", "cv_Tbio_august", "mean_HQ_august", "max_HQ_august", "cv_HQ_august")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_july")), ~ str_remove(., "_july")) %>%
    rename_at(vars(matches("july")), ~ str_remove(., "july")) %>%
    rename_at(vars(matches("_elk")), ~ str_remove(., "_elk"))
  station_stack_elk_aug <- station_stack_elk %>% dplyr::select(-c("mean_Tbio_july_kg_ha", "max_Tbio_july_kg_ha", "cv_Tbio_july", "mean_HQ_july", "max_HQ_july", "cv_HQ_july")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_august")), ~ str_remove(., "_august")) %>%
    rename_at(vars(matches("august")), ~ str_remove(., "august")) %>%
    rename_at(vars(matches("_elk")), ~ str_remove(., "_elk"))
  station_elk_list <- list(station_stack_elk_july, station_stack_elk_aug)
  
  station_stack_wtd_july <- station_stack_wtd %>% dplyr::select(-c("mean_Tbio_august_kg_ha", "max_Tbio_august_kg_ha", "cv_Tbio_august", "mean_HQ_august", "max_HQ_august", "cv_HQ_august")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_july")), ~ str_remove(., "_july")) %>%
    rename_at(vars(matches("july")), ~ str_remove(., "july"))%>%
    rename_at(vars(matches("_deer")), ~ str_remove(., "_deer"))
  station_stack_wtd_aug <- station_stack_wtd %>% dplyr::select(-c("mean_Tbio_july_kg_ha", "max_Tbio_july_kg_ha", "cv_Tbio_july", "mean_HQ_july", "max_HQ_july", "cv_HQ_july")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_august")), ~ str_remove(., "_august")) %>%
    rename_at(vars(matches("august")), ~ str_remove(., "august"))%>%
    rename_at(vars(matches("_deer")), ~ str_remove(., "_deer"))
  station_wtd_list <- list(station_stack_wtd_july, station_stack_wtd_aug)
  
  #'  Double check things are ordered correctly!!!!
  stations_elk_eoe20s[82:90,1:9]; DH_eoe20s_RNmod[[1]][82:90,1:3]; nrow(stations_elk_eoe20s); nrow(DH_eoe20s_RNmod[[1]])
  stations_elk_eoe21s[82:90,1:9]; DH_eoe21s_RNmod[[1]][82:90,1:3]; nrow(stations_elk_eoe21s); nrow(DH_eoe21s_RNmod[[1]])
  stations_elk_eoe22s[82:90,1:9]; DH_eoe22s_RNmod[[1]][82:90,1:3]; nrow(stations_elk_eoe22s); nrow(DH_eoe22s_RNmod[[1]])
  
  #'  Correlation matrix to check for collinearity among continuous variables
  corr_matrix <- function(dat, firstcol, lastcol) {
    continuous_variables <- dat[,firstcol:lastcol]
    corr_all <- cor(continuous_variables)
    corr_all <- as.data.frame(round(corr_all, 2))
    return(corr_all)
  }
  corr_matrix(station_stack_elk_july, firstcol = 5, lastcol = 15) # Mean & max Tbio (0.89), mean & max HQ (0.77), and total predicted & total selected (0.77) 
  corr_matrix(station_stack_elk_aug, firstcol = 5, lastcol = 15) # Mean & max Tbio (0.90), mean & max HQ (0.74), and total predicted & total selected (0.77) 
  corr_matrix(station_stack_wtd_july, firstcol = 5, lastcol = 15) # Mean & max Tbio (0.91), mean & max HQ (0.87), cv Tbio & cv HQ (0.69), and total predicted & total selected (0.79) 
  corr_matrix(station_stack_wtd_aug, firstcol = 5, lastcol = 15) # Mean & max Tbio (0.91), mean & max HQ (0.86), cv Tbio & cv HQ (0.69), and total predicted & total selected (0.79) 
  
  #####  Save!  #####
  #'  ----------
  #'  Seasonal detection histories
  save(DH_eoe20s_RNmod, file = "./Data/Side_quests/Hilger/DH_eoe20s_RNmod.RData")
  save(DH_eoe21s_RNmod, file = "./Data/Side_quests/Hilger/DH_eoe21s_RNmod.RData")
  save(DH_eoe22s_RNmod, file = "./Data/Side_quests/Hilger/DH_eoe22s_RNmod.RData")
  
  #'  Stacked detection histories (all years per species and month)
  save(DH_elk_list, file = "./Data/Side_quests/Hilger/DH_elk_RNmod.RData")
  save(DH_wtd_list, file = "./Data/Side_quests/Hilger/DH_wtd_RNmod.RData")
  
  #'  Seasonal camera station data
  save(stations_elk_eoe20s, file = "./Data/Side_quests/Hilger/stations_elk_eoe20s.RData")
  save(stations_elk_eoe21s, file = "./Data/Side_quests/Hilger/stations_elk_eoe21s.RData")
  save(stations_elk_eoe22s, file = "./Data/Side_quests/Hilger/stations_elk_eoe22s.RData")
  save(stations_wtd_eoe20s, file = "./Data/Side_quests/Hilger/stations_wtd_eoe20s.RData")
  save(stations_wtd_eoe21s, file = "./Data/Side_quests/Hilger/stations_wtd_eoe21s.RData")
  save(stations_wtd_eoe22s, file = "./Data/Side_quests/Hilger/stations_wtd_eoe22s.RData")
  
  #'  ---------------------------------
  ####  Royle-Nichols Abundance Model  ####
  #'  ---------------------------------
  #'  Run RN models for each predator species using a Bayesian framework in JAGS
  #'  ------------------------
  #####  Setup data for JAGS  #####
  #'  ------------------------
  #'  Bundle detection histories and covariates for each species and year
  bundle_dat <- function(dh, nsite, nsurvey, cov) {
    #'  Convert detection history to matrix
    dh <- as.matrix(dh)
    dimnames(dh) <- NULL
    #'  Count number of sites per GMU
    ncams_perGMU <- cov %>%
      group_by(GMU) %>%
      summarise(nsites = n()) %>%
      ungroup()
    #'  Bundle data for JAGS
    bundled <- list(y = dh, 
                    nsites = dim(dh)[1], 
                    nsurveys = dim(dh)[2],
                    nyear = max(as.numeric(cov$Season)),
                    ngmu = max(as.numeric(cov$GMU)),
                    nsets = max(as.numeric(cov$Setup)),
                    gmu = as.numeric(cov$GMU), 
                    setup = as.numeric(cov$Setup),
                    year = as.numeric(cov$Season, levels = c("Smr20", "Smr21", "Smr22")),
                    elev = as.numeric(scale(cov$elev)),
                    forest = as.numeric(scale(cov$perc_forest)),
                    mean_Tbio = as.numeric(scale(cov$mean_Tbio_kg_ha)),
                    max_Tbio = as.numeric(scale(cov$max_Tbio_kg_ha)),
                    cv_Tbio = as.numeric(scale(cov$cv_Tbio)),
                    mean_HQ = as.numeric(scale(cov$mean_HQ)),
                    max_HQ = as.numeric(scale(cov$max_HQ)),
                    cv_HQ = as.numeric(scale(cov$cv_HQ)),
                    total_selected = as.numeric(scale(cov$total_selected_species)),
                    total_predicted = as.numeric(scale(cov$total_predicted_species)),
                    prop_selected = as.numeric(scale(cov$prop_selected_species_weighted)))
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle_elk <- mapply(bundle_dat, dh = DH_elk_list, cov = station_elk_list, SIMPLIFY = FALSE)
  data_JAGS_bundle_wtd <- mapply(bundle_dat, dh = DH_wtd_list, cov = station_wtd_list, SIMPLIFY = FALSE)
  
  #'  Initial values
  #'  Using naive occupancy as a starting point for local abundance
  initial_n <- function(dh) {
    #'  Max value per row
    ninit <- apply(dh, 1, max, na.rm = TRUE)
    ninit <- as.vector(ninit)
    return(ninit)
  }
  #'  Apply function per species for each year
  ninit_elk <- lapply(DH_elk_list, initial_n)
  ninit_wtd <- lapply(DH_wtd_list, initial_n)
  
  #'  Parameters monitored
  params <- c("beta0", "b.year", "alpha0", "a.setup", 
              "b.meanTbio", "b.maxTbio", "b.cvTbio", 
              "b.meanHQ", "b.maxHQ", "b.cvHQ", 
              "b.selected", "b.predicted", "b.prop.selected", 
              "rSetup", "mu.r", "mean.p", "N", "loglike.new")
  #'  NOTE about mean vs mu lambda and r: 
  #'  mean.lambda = the intercept based on reference category, i.e., mean lambda for 2020 
  #'  mean.r = the intercept based on reference category, i.e., per-individual detection probability at random sites
  #'  mu.lambda = lambda averaged across all GMUs
  #'  mu.r = per-individual detection probability averaged across all sites 
  
  #'  MCMC settings
  nc <- 3
  ni <- 75000
  nb <- 5000
  nt <- 10
  na <- 5000
  
  #'  Load competing models
  #'  Univariate models to compare mean & max biomass variables, and selected vs predicted community composition
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_meanHQ.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_maxHQ.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_meanTbio.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_maxTbio.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_totalSelected.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_totalPredicted.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_null.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_HQ_mean.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_HQ_max.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_Tbio_mean.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_selected.propSelected.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_predicted.propSelected.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_elk_july_global.R") 
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_elk_aug_global.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_july_global_1.R") 
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_july_global_2.R") 
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_aug_global_1.R") 
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_aug_global_2.R") 
  
  #'  Function to calculate joint-likelihood approach to WAIC (WAICj) based on Gaya & Ketz (2024)
  calc.jointlike <- function(x){
    like <- as.matrix(x$sims.list$loglike.new) #log-likelihood for every iteration and camera
    fbar <- colMeans(exp(like)) #mean likelihood 
    Pw <- sum(apply(like,2,var)) #mean variance in log-likelihood 
    WAIC_ish<- -2*sum(log(fbar))+2*Pw
    return(WAIC_ish)
  }
  
  #'  -----------------------
  #####  Elk July RN models  #####
  #'  -----------------------
  ######  Univariate models ######  # MaxHQ, MaxTbio, TotalSelected have more DIC support
  #'  Mean HQ (DIC = 24540.61)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_meanHQ <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanHQ.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_meanHQ$summary[1:15,])
  (RN_elk_july_meanHQ_WAICj <- calc.jointlike(RN_elk_july_meanHQ))
  print(RN_elk_july_meanHQ$DIC)
  save(RN_elk_july_meanHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_meanHQ_", Sys.Date(), ".RData"))
  
  #'  Max HQ  (DIC = 24493.78) **
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_maxHQ <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxHQ.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_maxHQ$summary[1:15,])
  (RN_elk_july_maxHQ_WAICj <- calc.jointlike(RN_elk_july_maxHQ))
  print(RN_elk_july_maxHQ$DIC)
  save(RN_elk_july_maxHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_maxHQ_", Sys.Date(), ".RData"))
  
  #'  Mean Tbio (DIC = 24546.61)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_meanTbio <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanTbio.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_meanTbio$summary[1:15,])
  (RN_elk_july_meanTbio_WAICj <- calc.jointlike(RN_elk_july_meanTbio))
  print(RN_elk_july_meanTbio$DIC)
  save(RN_elk_july_meanTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_meanTbio_", Sys.Date(), ".RData"))
  
  #'  Max Tbio (DIC = 24544.07) **
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_maxTbio <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxTbio.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_maxTbio$summary[1:15,])
  (RN_elk_july_maxTbio_WAICj <- calc.jointlike(RN_elk_july_maxTbio))
  print(RN_elk_july_maxTbio$DIC)
  save(RN_elk_july_maxTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_maxTbio_", Sys.Date(), ".RData"))
  
  #'  Total Selected (DIC = 24521.24) **
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_selected <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalSelected.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_selected$summary[1:15,])
  (RN_elk_july_selected_WAICj <- calc.jointlike(RN_elk_july_selected))
  print(RN_elk_july_selected$DIC)
  save(RN_elk_july_selected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_selected_", Sys.Date(), ".RData"))
  
  #'  Total Predicted (DIC = 24545.29)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_predicted <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalPredicted.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_predicted$summary[1:15,])
  (RN_elk_july_predicted_WAICj <- calc.jointlike(RN_elk_july_predicted))
  print(RN_elk_july_predicted$DIC)
  save(RN_elk_july_predicted, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_predicted_", Sys.Date(), ".RData"))
  
  ######  Model 1  ######
  #'  Null  (DIC = 24518.87; WAICj = 191348.6)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_null <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                      "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_null$summary)
  (RN_elk_july_null_WAICj <- calc.jointlike(RN_elk_july_null))
  print(RN_elk_july_null$DIC)
  which(RN_elk_july_null$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_null$samples)
  save(RN_elk_july_null, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_null_", Sys.Date(), ".RData"))
  
  ######  High Quality Biomass  ######
  #'  Mean & CV  (DIC = 24537.45; WAICj = 205032.1)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_HQ_mean.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                      "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_mean.cv.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_HQ_mean.cv$summary[1:15,])
  (RN_elk_july_HQ_mean.cv_WAICj <- calc.jointlike(RN_elk_july_HQ_mean.cv))
  print(RN_elk_july_HQ_mean.cv$DIC)
  which(RN_elk_july_HQ_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_HQ_mean.cv$samples)
  save(RN_elk_july_HQ_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_HQ_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV  (DIC = 24535.7; WAICj = 204042.1)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_HQ_max.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                                 "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_max.cv.txt",
                                 n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                 n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_HQ_max.cv$summary[1:15,])
  (RN_elk_july_HQ_max.cv_WAICj <- calc.jointlike(RN_elk_july_HQ_max.cv))
  print(RN_elk_july_HQ_max.cv$DIC)
  which(RN_elk_july_HQ_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_HQ_max.cv$samples)
  save(RN_elk_july_HQ_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_HQ_max.cv_", Sys.Date(), ".RData"))
  
  ######  Total Biomass  ######
  #'  Mean & CV  (DIC = 24547.04; WAICj = 195755.1)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_Tbio_mean.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_mean.cv.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_Tbio_mean.cv$summary[1:15,])
  (RN_elk_july_Tbio_mean.cv_WAICj <- calc.jointlike(RN_elk_july_Tbio_mean.cv))
  print(RN_elk_july_Tbio_mean.cv$DIC)
  which(RN_elk_july_Tbio_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_Tbio_mean.cv$samples)
  save(RN_elk_july_Tbio_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_Tbio_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV   (DIC = 24537.56; WAICj = 192491.1)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_Tbio_max.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                                   "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.txt",
                                   n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                   n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_Tbio_max.cv$summary[1:15,])
  (RN_elk_july_Tbio_max.cv_WAICj <- calc.jointlike(RN_elk_july_Tbio_max.cv))
  print(RN_elk_july_Tbio_max.cv$DIC)
  which(RN_elk_july_Tbio_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_Tbio_max.cv$samples)
  save(RN_elk_july_Tbio_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_Tbio_max.cv_", Sys.Date(), ".RData"))
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected  (DIC = 24521.25; WAICj = 191263)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_selected.propSelected <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_selected.propSelected.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_selected.propSelected$summary[1:15,])
  (RN_elk_july_selected.propSelected_WAICj <- calc.jointlike(RN_elk_july_selected.propSelected))
  print(RN_elk_july_selected.propSelected$DIC)
  which(RN_elk_july_selected.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_selected.propSelected$samples)
  save(RN_elk_july_selected.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_selected.propSelected_", Sys.Date(), ".RData"))
  
  #'  Predicted & Proportion Selected  (DIC = 24559.59; WAICj = 194585.3)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_predicted.propSelected <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_predicted.propSelected.txt",
                                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_predicted.propSelected$summary[1:15,])
  (RN_elk_july_predicted.propSelected_WAICj <- calc.jointlike(RN_elk_july_predicted.propSelected))
  print(RN_elk_july_predicted.propSelected$DIC)
  which(RN_elk_july_predicted.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_predicted.propSelected$samples)
  save(RN_elk_july_predicted.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_predicted.propSelected_", Sys.Date(), ".RData"))
  
  ######  GLOBAL  ######
  #'  Using only most supported non-correlated covariates (based on WAICj from models above)
  #'  Global: Max HQ, cv HQ, Max Tbio, cv Tbio, Selected, PropSelected
  #'  (DIC = 24532.76; WAICj = 205372.6)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_global <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_elk_july_global.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_global$summary[1:15,])
  (RN_elk_july_global_WAICj <- calc.jointlike(RN_elk_july_global))
  print(RN_elk_july_global$DIC)
  which(RN_elk_july_global$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_global$samples)
  save(RN_elk_july_global, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_global_", Sys.Date(), ".RData"))
  
  #'  -------------------------
  #####  Elk August RN models  #####
  #'  -------------------------
  ######  Univariate models ######  # MaxHQ, MaxTbio, TotalSelected have more DIC support
  #'  Mean HQ (DIC = 23693.48)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_meanHQ <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanHQ.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_meanHQ$summary[1:15,])
  (RN_elk_aug_meanHQ_WAICj <- calc.jointlike(RN_elk_aug_meanHQ))
  print(RN_elk_aug_meanHQ$DIC)
  save(RN_elk_aug_meanHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_meanHQ_", Sys.Date(), ".RData"))
  
  #'  Max HQ (DIC = 23690.05) **
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_maxHQ <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxHQ.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_maxHQ$summary[1:15,])
  (RN_elk_aug_maxHQ_WAICj <- calc.jointlike(RN_elk_aug_maxHQ))
  print(RN_elk_aug_maxHQ$DIC)
  save(RN_elk_aug_maxHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_maxHQ_", Sys.Date(), ".RData"))
  
  #'  Mean Tbio (DIC = 23696.5)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_meanTbio <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanTbio.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_meanTbio$summary[1:15,])
  (RN_elk_aug_meanTbio_WAICj <- calc.jointlike(RN_elk_aug_meanTbio))
  print(RN_elk_aug_meanTbio$DIC)
  save(RN_elk_aug_meanTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_meanTbio_", Sys.Date(), ".RData"))
  
  #'  Max Tbio (DIC = 23692.49) **
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_maxTbio <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                              "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxTbio.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_maxTbio$summary[1:15,])
  (RN_elk_aug_maxTbio_WAICj <- calc.jointlike(RN_elk_aug_maxTbio))
  print(RN_elk_aug_maxTbio$DIC)
  save(RN_elk_aug_maxTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_maxTbio_", Sys.Date(), ".RData"))
  
  #'  Total Selected (DIC = 23666.4) **
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_selected <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalSelected.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_selected$summary[1:15,])
  (RN_elk_aug_selected_WAICj <- calc.jointlike(RN_elk_aug_selected))
  print(RN_elk_aug_selected$DIC)
  save(RN_elk_aug_selected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_selected_", Sys.Date(), ".RData"))
  
  #'  Total Predicted (DIC = 23686.04)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_predicted <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalPredicted.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_predicted$summary)
  (RN_elk_aug_predicted_WAICj <- calc.jointlike(RN_elk_aug_predicted))
  print(RN_elk_aug_predicted$DIC)
  save(RN_elk_aug_predicted, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_predicted_", Sys.Date(), ".RData"))
  
  ######  Model 1  ######
  #'  Null  (DIC = 23690.13; WAICj = 110965.5)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_null <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                      "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_null$summary)
  (RN_elk_aug_null_WAICj <- calc.jointlike(RN_elk_aug_null))
  print(RN_elk_aug_null$DIC)
  which(RN_elk_aug_null$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_null$samples)
  save(RN_elk_aug_null, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_null_", Sys.Date(), ".RData"))
  
  ######  High Quality Biomass  ######
  #'  Mean & CV   (DIC = 23685.62; WAICj = 118302.8)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_HQ_mean.cv <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                     "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_mean.cv.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_HQ_mean.cv$summary[1:15,])
  (RN_elk_aug_HQ_mean.cv_WAICj <- calc.jointlike(RN_elk_aug_HQ_mean.cv))
  print(RN_elk_aug_HQ_mean.cv$DIC)
  which(RN_elk_aug_HQ_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_HQ_mean.cv$samples)
  save(RN_elk_aug_HQ_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_HQ_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV  (DIC = 23667.9; WAICj = 117231.4)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_HQ_max.cv <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_max.cv.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_HQ_max.cv$summary[1:15,])
  (RN_elk_aug_HQ_max.cv_WAICj <- calc.jointlike(RN_elk_aug_HQ_max.cv))
  print(RN_elk_aug_HQ_max.cv$DIC)
  which(RN_elk_aug_HQ_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_HQ_max.cv$samples)
  save(RN_elk_aug_HQ_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_HQ_max.cv_", Sys.Date(), ".RData"))
  
  ######  Total Biomass  ######
  #'  Mean & CV   (DIC = 23672.75; WAICj = 113111)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_Tbio_mean.cv <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_mean.cv.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_Tbio_mean.cv$summary[1:15,])
  (RN_elk_aug_Tbio_mean.cv_WAICj <- calc.jointlike(RN_elk_aug_Tbio_mean.cv))
  print(RN_elk_aug_Tbio_mean.cv$DIC)
  which(RN_elk_aug_Tbio_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_Tbio_mean.cv$samples)
  save(RN_elk_aug_Tbio_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_Tbio_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV  (DIC = 23679.73; WAICj = 112069.8)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_Tbio_max.cv <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                                  "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.txt",
                                  n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                  n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_Tbio_max.cv$summary[1:15,])
  (RN_elk_aug_Tbio_max.cv_WAICj <- calc.jointlike(RN_elk_aug_Tbio_max.cv))
  print(RN_elk_aug_Tbio_max.cv$DIC)
  which(RN_elk_aug_Tbio_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_Tbio_max.cv$samples)
  save(RN_elk_aug_Tbio_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_Tbio_max.cv_", Sys.Date(), ".RData"))
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected  (DIC = 23674.87; WAICj = 111953.8)
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_selected.propSelected <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_selected.propSelected.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_selected.propSelected$summary[1:15,])
  (RN_elk_aug_selected.propSelected_WAICj <- calc.jointlike(RN_elk_aug_selected.propSelected))
  print(RN_elk_aug_selected.propSelected$DIC)
  which(RN_elk_aug_selected.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_selected.propSelected$samples)
  save(RN_elk_aug_selected.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_selected.propSelected_", Sys.Date(), ".RData"))
  
  #'  Predicted & Proportion Selected  (DIC = 23693.75; WAICj = 112792.3)
  start.time = Sys.time()  
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_predicted.propSelected <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_predicted.propSelected.txt",
                                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_predicted.propSelected$summary[1:15,])
  (RN_elk_aug_predicted.propSelected_WAICj <- calc.jointlike(RN_elk_aug_predicted.propSelected))
  print(RN_elk_aug_predicted.propSelected$DIC)
  which(RN_elk_aug_predicted.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_predicted.propSelected$samples)
  save(RN_elk_aug_predicted.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_predicted.propSelected_", Sys.Date(), ".RData"))
  
  ######  GLOBAL  ######
  #'  Using only most supported non-correlated covariates (based on WAICj from models above)
  #'  Global: Max HQ, cv HQ, Max Tbio, cv Tbio, selected, PropSelected
  #'  (DIC = 23659.2; WAICj = 118713.5)
  start.time = Sys.time() 
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_global <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_elk_aug_global.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_global$summary[1:15,])
  (RN_elk_aug_global_WAICj <- calc.jointlike(RN_elk_aug_global))
  print(RN_elk_aug_global$DIC)
  which(RN_elk_aug_global$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_global$samples)
  save(RN_elk_aug_global, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_global_", Sys.Date(), ".RData"))
  
  #'  -----------------------
  #####  WTD July RN models  #####
  #'  -----------------------
  ######  Univariate models ######  # MeanHQ, MaxTbio, TotalSelected have more DIC support
  #'  Mean HQ (DIC = 40592.29) **
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_meanHQ <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanHQ.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_meanHQ$summary[1:15,])
  (RN_wtd_july_meanHQ_WAICj <- calc.jointlike(RN_wtd_july_meanHQ))
  print(RN_wtd_july_meanHQ$DIC)
  save(RN_wtd_july_meanHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_meanHQ_", Sys.Date(), ".RData"))
  
  #'  Max HQ  (DIC = 40642.9)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_maxHQ <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxHQ.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_maxHQ$summary[1:15,])
  (RN_wtd_july_maxHQ_WAICj <- calc.jointlike(RN_wtd_july_maxHQ))
  print(RN_wtd_july_maxHQ$DIC)
  save(RN_wtd_july_maxHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_maxHQ_", Sys.Date(), ".RData"))
  
  #'  Mean Tbio (DIC = 40619.1)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_meanTbio <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanTbio.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_meanTbio$summary[1:15,])
  (RN_wtd_july_meanTbio_WAICj <- calc.jointlike(RN_wtd_july_meanTbio))
  print(RN_wtd_july_meanTbio$DIC)
  save(RN_wtd_july_meanTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_meanTbio_", Sys.Date(), ".RData"))
  
  #'  Max Tbio (DIC = 40602.09) **
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_maxTbio <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                              "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxTbio.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_maxTbio$summary[1:15,])
  (RN_wtd_july_maxTbio_WAICj <- calc.jointlike(RN_wtd_july_maxTbio))
  print(RN_wtd_july_maxTbio$DIC)
  save(RN_wtd_july_maxTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_maxTbio_", Sys.Date(), ".RData"))
  
  #'  Total Selected (DIC = 40597.45) **
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_selected <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalSelected.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_selected$summary[1:15,])
  (RN_wtd_july_selected_WAICj <- calc.jointlike(RN_wtd_july_selected))
  print(RN_wtd_july_selected$DIC)
  save(RN_wtd_july_selected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_selected_", Sys.Date(), ".RData"))
  
  #'  Total Predicted (DIC = 40645.89)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_predicted <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalPredicted.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_predicted$summary[1:15,])
  (RN_wtd_july_predicted_WAICj <- calc.jointlike(RN_wtd_july_predicted))
  print(RN_wtd_july_predicted$DIC)
  save(RN_wtd_july_predicted, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_predicted_", Sys.Date(), ".RData"))
  
  ######  Model 1  ######
  #'  Null  (DIC = 40655.86; WAICj = 1537895)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_null <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_null$summary)
  (RN_wtd_july_null_WAICj <- calc.jointlike(RN_wtd_july_null))
  print(RN_wtd_july_null$DIC)
  which(RN_wtd_july_null$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_null$samples)
  save(RN_wtd_july_null, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_null_", Sys.Date(), ".RData"))
  
  ######  High Quality Biomass  ######
  #'  Mean & CV  (DIC = 40558.56; WAICj = 1767083)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_HQ_mean.cv <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_mean.cv.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_HQ_mean.cv$summary[1:15,])
  (RN_wtd_july_HQ_mean.cv_WAICj <- calc.jointlike(RN_wtd_july_HQ_mean.cv))
  print(RN_wtd_july_HQ_mean.cv$DIC)
  which(RN_wtd_july_HQ_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_HQ_mean.cv$samples)
  save(RN_wtd_july_HQ_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_HQ_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV  (DIC = 40584.65; WAICj = 1694535)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_HQ_max.cv <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                                 "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_max.cv.txt",
                                 n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                 n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_HQ_max.cv$summary[1:15,])
  (RN_wtd_july_HQ_max.cv_WAICj <- calc.jointlike(RN_wtd_july_HQ_max.cv))
  print(RN_wtd_july_HQ_max.cv$DIC)
  which(RN_wtd_july_HQ_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_HQ_max.cv$samples)
  save(RN_wtd_july_HQ_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_HQ_max.cv_", Sys.Date(), ".RData"))
  
  ######  Total Biomass  ######
  #'  Mean & CV  (DIC = 40632.43; WAICj = 1568780)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_Tbio_mean.cv <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_mean.cv.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_Tbio_mean.cv$summary[1:15,])
  (RN_wtd_july_Tbio_mean.cv_WAICj <- calc.jointlike(RN_wtd_july_Tbio_mean.cv))
  print(RN_wtd_july_Tbio_mean.cv$DIC)
  which(RN_wtd_july_Tbio_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_Tbio_mean.cv$samples)
  save(RN_wtd_july_Tbio_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_Tbio_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV  (DIC = 40615.6; WAICj = 1602163)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_Tbio_max.cv <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                                   "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.txt",
                                   n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                   n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_Tbio_max.cv$summary[1:15,])
  (RN_wtd_july_Tbio_max.cv_WAICj <- calc.jointlike(RN_wtd_july_Tbio_max.cv))
  print(RN_wtd_july_Tbio_max.cv$DIC)
  which(RN_wtd_july_Tbio_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_Tbio_max.cv$samples)
  save(RN_wtd_july_Tbio_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_Tbio_max.cv_", Sys.Date(), ".RData"))
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected  (DIC = 40558.06; WAICj = 1796101)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_selected.propSelected <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_selected.propSelected.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_selected.propSelected$summary[1:15,])
  (RN_wtd_july_selected.propSelected_WAICj <- calc.jointlike(RN_wtd_july_selected.propSelected))
  print(RN_wtd_july_selected.propSelected$DIC)
  which(RN_wtd_july_selected.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_selected.propSelected$samples)
  save(RN_wtd_july_selected.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_selected.propSelected_", Sys.Date(), ".RData"))
  
  #'  Predicted & Proportion Selected  (DIC = 40520.48; WAICj = 1886611)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_predicted.propSelected <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_predicted.propSelected.txt",
                                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_predicted.propSelected$summary[1:15,])
  (RN_wtd_july_predicted.propSelected_WAICj <- calc.jointlike(RN_wtd_july_predicted.propSelected))
  print(RN_wtd_july_predicted.propSelected$DIC)
  which(RN_wtd_july_predicted.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_predicted.propSelected$samples)
  save(RN_wtd_july_predicted.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_predicted.propSelected_", Sys.Date(), ".RData"))
  
  ######  GLOBAL  ######
  #'  Using only most supported non-correlated covariates
  #'  Global: Max HQ, cv HQ, Mean Tbio, Selected, PropSelected
  #'  cvHQ & cvTbio are highly correlated so running global model with only one cv biomass variable at a time
  #'  (DIC = 40555.07; WAICj = 1969482)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_global1 <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_wtd_july_global_1.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_global1$summary[1:15,])
  (RN_wtd_july_global1_WAICj <- calc.jointlike(RN_wtd_july_global1))
  print(RN_wtd_july_global1$DIC)
  which(RN_wtd_july_global1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_global1$samples)
  save(RN_wtd_july_global1, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_global.cvHQ_", Sys.Date(), ".RData"))
  
  #' (DIC = 40574.9; WAICj = 1873468)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_global2 <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_wtd_july_global_2.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_global2$summary)
  (RN_wtd_july_global2_WAICj <- calc.jointlike(RN_wtd_july_global2))
  print(RN_wtd_july_global2$DIC)
  which(RN_wtd_july_global2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_global2$samples)
  save(RN_wtd_july_global2, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_global.cvTbio_", Sys.Date(), ".RData"))
  
  #'  -------------------------
  #####  WTD August RN models  #####
  #'  -------------------------
  ######  Univariate models ######  # MeanHQ, MaxTbio, TotalSelected have more DIC support
  #'  Mean HQ (DIC = 39910.93) **
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_meanHQ <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanHQ.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_meanHQ$summary[1:15,])
  (RN_wtd_aug_meanHQ_WAICj <- calc.jointlike(RN_wtd_aug_meanHQ))
  print(RN_wtd_aug_meanHQ$DIC)
  save(RN_wtd_aug_meanHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_meanHQ_", Sys.Date(), ".RData"))
  
  #'  Max HQ  (DIC = 39937.76)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_maxHQ <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxHQ.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_maxHQ$summary[1:15,])
  (RN_wtd_aug_maxHQ_WAICj <- calc.jointlike(RN_wtd_aug_maxHQ))
  print(RN_wtd_aug_maxHQ$DIC)
  save(RN_wtd_aug_maxHQ, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_maxHQ_", Sys.Date(), ".RData"))
  
  #'  Mean Tbio (DIC = 39927.51)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_meanTbio <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanTbio.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_meanTbio$summary[1:15,])
  (RN_wtd_aug_meanTbio_WAICj <- calc.jointlike(RN_wtd_aug_meanTbio))
  print(RN_wtd_aug_meanTbio$DIC)
  save(RN_wtd_aug_meanTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_meanTbio_", Sys.Date(), ".RData"))
  
  #'  Max Tbio (DIC = 39915.58) **
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_maxTbio <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                              "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxTbio.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_maxTbio$summary[1:15,])
  (RN_wtd_aug_maxTbio_WAICj <- calc.jointlike(RN_wtd_aug_maxTbio))
  print(RN_wtd_aug_maxTbio$DIC)
  save(RN_wtd_aug_maxTbio, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_maxTbio_", Sys.Date(), ".RData"))
  
  #'  Total Selected (DIC = 39889.78) **
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_selected <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalSelected.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_selected$summary[1:15,])
  (RN_wtd_aug_selected_WAICj <- calc.jointlike(RN_wtd_aug_selected))
  print(RN_wtd_aug_selected$DIC)
  save(RN_wtd_aug_selected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_selected_", Sys.Date(), ".RData"))
  
  #'  Total Predicted (DIC = 39893.59)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_predicted <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalPredicted.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_predicted$summary[1:15,])
  (RN_wtd_aug_predicted_WAICj <- calc.jointlike(RN_wtd_aug_predicted))
  print(RN_wtd_aug_predicted$DIC)
  save(RN_wtd_aug_predicted, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_predicted_", Sys.Date(), ".RData"))
  
  ######  Model 1  ######
  #'  Null  (DIC = 39916.8; WAICj = 1006710)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_null <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_null$summary)
  (RN_wtd_aug_null_WAICj <- calc.jointlike(RN_wtd_aug_null))
  print(RN_wtd_aug_null$DIC)
  which(RN_wtd_aug_null$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_null$samples)
  save(RN_wtd_aug_null, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_null_", Sys.Date(), ".RData"))
  
  ######  High Quality Biomass models  ######
  #'  Mean & CV  (DIC = 39911.87; WAICj = 1156182)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_HQ_mean.cv <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_mean.cv.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_HQ_mean.cv$summary[1:15,])
  (RN_wtd_aug_HQ_mean.cv_WAICj <- calc.jointlike(RN_wtd_aug_HQ_mean.cv))
  print(RN_wtd_aug_HQ_mean.cv$DIC)
  which(RN_wtd_aug_HQ_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_HQ_mean.cv$samples)
  save(RN_wtd_aug_HQ_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_HQ_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV  (DIC = 39933.27; WAICj = 1085812)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_HQ_max.cv <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_max.cv.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_HQ_max.cv$summary[1:15,])
  (RN_wtd_aug_HQ_max.cv_WAICj <- calc.jointlike(RN_wtd_aug_HQ_max.cv))
  print(RN_wtd_aug_HQ_max.cv$DIC)
  which(RN_wtd_aug_HQ_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_HQ_max.cv$samples)
  save(RN_wtd_aug_HQ_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_HQ_max.cv_", Sys.Date(), ".RData"))
  
  ######  Total Biomass  ######
  #'  Mean & CV  (DIC = 39946.98; WAICj = 1015780)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_Tbio_mean.cv <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_mean.cv.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_Tbio_mean.cv$summary[1:15,])
  (RN_wtd_aug_Tbio_mean.cv_WAICj <- calc.jointlike(RN_wtd_aug_Tbio_mean.cv))
  print(RN_wtd_aug_Tbio_mean.cv$DIC)
  which(RN_wtd_aug_Tbio_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_Tbio_mean.cv$samples)
  save(RN_wtd_aug_Tbio_mean.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_Tbio_mean.cv_", Sys.Date(), ".RData"))
  
  #'  Max & CV  (DIC = 39898.43; WAICj = 1041801)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_Tbio_max.cv <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                                  "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.txt",
                                  n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                  n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_Tbio_max.cv$summary[1:15,])
  (RN_wtd_aug_Tbio_max.cv_WAICj <- calc.jointlike(RN_wtd_aug_Tbio_max.cv))
  print(RN_wtd_aug_Tbio_max.cv$DIC)
  which(RN_wtd_aug_Tbio_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_Tbio_max.cv$samples)
  save(RN_wtd_aug_Tbio_max.cv, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_Tbio_max.cv_", Sys.Date(), ".RData"))
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected  (DIC = 39859.71; WAICj = 1189218)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_selected.propSelected <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_selected.propSelected.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_selected.propSelected$summary[1:15,])
  (RN_wtd_aug_selected.propSelected_WAICj <- calc.jointlike(RN_wtd_aug_selected.propSelected))
  print(RN_wtd_aug_selected.propSelected$DIC)
  which(RN_wtd_aug_selected.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_selected.propSelected$samples)
  save(RN_wtd_aug_selected.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_selected.propSelected_", Sys.Date(), ".RData"))
  
  #'  Predicted & Proportion Selected  (DIC = 39814.27; WAICj = 1223039)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_predicted.propSelected <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_predicted.propSelected.txt",
                                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_predicted.propSelected$summary[1:15,])
  (RN_wtd_aug_predicted.propSelected_WAICj <- calc.jointlike(RN_wtd_aug_predicted.propSelected))
  print(RN_wtd_aug_predicted.propSelected$DIC)
  which(RN_wtd_aug_predicted.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_predicted.propSelected$samples)
  save(RN_wtd_aug_predicted.propSelected, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_predicted.propSelected_", Sys.Date(), ".RData"))
  
  ######  GLOBAL  ######
  #'  Using only most supported non-correlated covariates
  #'  Global: Max HQ, cv HQ, Mean Tbio, cv Tbio, selected, PropSelected
  #'  cvHQ & cvTbio are highly correlated so running global model with only one cv biomass variable at a time
  #'  (DIC = 39881.91; WAICj = 1312306)
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_global1 <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_wtd_aug_global_1.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_global1$summary)
  (RN_wtd_aug_global1_WAICj <- calc.jointlike(RN_wtd_aug_global1))
  print(RN_wtd_aug_global1$DIC)
  which(RN_wtd_aug_global1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_global1$samples)
  save(RN_wtd_aug_global1, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_global.cvHQ_", Sys.Date(), ".RData"))
  
  #'  (DIC = 39877.86; WAICj = 1257700)
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_global2 <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                                  "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_wtd_aug_global_2.txt",
                                  n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                  n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_global2$summary)
  (RN_wtd_aug_global2_WAICj <- calc.jointlike(RN_wtd_aug_global2))
  print(RN_wtd_aug_global2$DIC)
  which(RN_wtd_aug_global2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_global2$samples)
  save(RN_wtd_aug_global2, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_global.cvTbio_", Sys.Date(), ".RData"))
  
