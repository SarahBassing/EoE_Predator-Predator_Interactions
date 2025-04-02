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
  file_paths_elk <- list.files(path = "./Data/side_quests/Hilger/Ungulate density covariates kg_ha/Final covariates/Elk", pattern = "\\.RDS", full.names = TRUE)
  file_names_elk <- gsub(pattern = "\\.RDS$", replacement = "", x = basename(file_paths_elk))
  forage_list_elk <- lapply(file_paths_elk, readRDS)
  names(forage_list_elk) <- c("Elk2020_kg_ha", "Elk2021_kg_ha", "Elk2022_kg_ha")
  
  #'  Load forage data for white-tailed deer
  file_paths_wtd <- list.files(path = "./Data/side_quests/Hilger/Ungulate density covariates kg_ha/Final covariates/WTD", pattern = "\\.RDS", full.names = TRUE)
  file_names_wtd <- gsub(pattern = "\\.RDS$", replacement = "", x = basename(file_paths_wtd))
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
    numnumcovs_by_month <- numnumcovs %>%
      mutate(Month = ifelse(Month == "07", "july", "august")) %>%
      pivot_wider(names_from = Month, values_from = c(mean_HQBio_kg, max_HQBio_kg,
                                                      cv_HQBio_kg, mean_TBio_kg, 
                                                      max_TBio_kg, cv_TBio_kg, 
                                                      selected, total, prop_selected))   
    covs <- full_join(sitecovs, numnumcovs_by_month, by = c("NewLocationID" = "original_id"))
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
                                 outDir = "./Data/Side_quests/Hilger")
    
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
  
  #'  Grab effort for each season (will be same for both species)
  grab_effort <- function(dh) {
    #'  Keep only effort
    effort_only <- dh[[2]]
    effort_july <- effort_only[,1:31]
    effort_august <- effort_only[,32:62]
    july_nDays <- rowSums(effort_july, na.rm = TRUE)
    august_nDays <- rowSums(effort_august, na.rm = TRUE)
    effort <- bind_cols(july_nDays, august_nDays) %>% as.data.frame(.)
    names(effort) <- c("July_nDays", "August_nDays")
    return(effort)
  }
  Effort_eoe20s_RNmod <- grab_effort(DHeff_eoe20s_RNmod[[1]])
  Effort_eoe21s_RNmod <- grab_effort(DHeff_eoe21s_RNmod[[1]])
  Effort_eoe22s_RNmod <- grab_effort(DHeff_eoe22s_RNmod[[1]])
  
  #'  Snag rownames per annaul DH for each species (helpful for reporting results later on)
  grab_location_id <- function(dh, yr) {
    locs <- rownames(dh)
    locs <- as.data.frame(locs)
    names(locs) <- "NewLocationID"
    locs$Season <- yr
    return(locs)
  }
  locs_eoe20s <- lapply(DH_eoe20s_RNmod, grab_location_id, yr = "Smr20")#; save(locs_eoe20s, file = "./Data/Side_quests/Hilger/DH_rownames_eoe20s.RData")
  locs_eoe21s <- lapply(DH_eoe21s_RNmod, grab_location_id, yr = "Smr21")#; save(locs_eoe21s, file = "./Data/Side_quests/Hilger/DH_rownames_eoe21s.RData")
  locs_eoe22s <- lapply(DH_eoe22s_RNmod, grab_location_id, yr = "Smr22")#; save(locs_eoe22s, file = "./Data/Side_quests/Hilger/DH_rownames_eoe22s.RData")
  stacked_rownames_elk <- bind_rows(locs_eoe20s[[1]], locs_eoe21s[[1]], locs_eoe22s[[1]])
  stacked_rownames_wtd <- bind_rows(locs_eoe20s[[2]], locs_eoe21s[[2]], locs_eoe22s[[2]])
  # save(stacked_rownames_elk, file = "./Data/Side_quests/Hilger/stacked_rownames_elk_02.03.25.RData")
  # save(stacked_rownames_wtd, file = "./Data/Side_quests/Hilger/stacked_rownames_wtd_02.03.25.RData")
  
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
  station_stack_elk_july <- station_stack_elk %>% dplyr::select(-contains("_august")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_july")), ~ str_remove(., "_july")) 
  station_stack_elk_aug <- station_stack_elk %>% dplyr::select(-contains("_july")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_august")), ~ str_remove(., "_august")) 
  station_elk_list <- list(station_stack_elk_july, station_stack_elk_aug)
  
  station_stack_wtd_july <- station_stack_wtd %>% dplyr::select(-contains("_august")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_july")), ~ str_remove(., "_july")) 
  station_stack_wtd_aug <- station_stack_wtd %>% dplyr::select(-contains("_july")) %>%
    #'  Remove the month identifier in each column so covariate names are identical across data sets
    rename_at(vars(matches("_august")), ~ str_remove(., "_august")) 
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
  corr_matrix(station_stack_elk_july, firstcol = 5, lastcol = 15) # Mean & max HQ (0.66), max HQ & max Tbio (-0.67), cvHQ & cvTbio (0.81), mean & max Tbio (0.85), and total predicted & total selected (0.77) 
  corr_matrix(station_stack_elk_aug, firstcol = 5, lastcol = 15) # Mean & max HQ (0.67), max HQ & max Tbio (-0.66), cvHQ & cvTbio (0.81), mean & max Tbio (0.82), and total predicted & total selected (0.77) 
  corr_matrix(station_stack_wtd_july, firstcol = 5, lastcol = 15) # Mean & max HQ (0.72), cvHQ & cvTbio (0.95), mean & max Tbio (0.87), and total predicted & total selected (0.79) 
  corr_matrix(station_stack_wtd_aug, firstcol = 5, lastcol = 15) # Mean & max HQ (0.72), cvHQ & cvTbio (0.95), mean & max Tbio (0.85), and total predicted & total selected (0.79) 
  
  #####  Save!  #####
  #'  ----------
  #'  Unique detection events
  save(eoe20s_det_events, file = "./Data/Side_quests/Hilger/DetectionEvents_eoe20s.RData")
  save(eoe21s_det_events, file = "./Data/Side_quests/Hilger/DetectionEvents_eoe21s.RData")
  save(eoe22s_det_events, file = "./Data/Side_quests/Hilger/DetectionEvents_eoe22s.RData")
  
  #'  Seasonal detection histories
  save(DH_eoe20s_RNmod, file = "./Data/Side_quests/Hilger/DH_eoe20s_RNmod.RData")
  save(DH_eoe21s_RNmod, file = "./Data/Side_quests/Hilger/DH_eoe21s_RNmod.RData")
  save(DH_eoe22s_RNmod, file = "./Data/Side_quests/Hilger/DH_eoe22s_RNmod.RData")
  
  #'  Save sampling effort
  save(Effort_eoe20s_RNmod, file = "./Data/Side_quests/Hilger/Effort_eoe20s_RNmod.RData")
  save(Effort_eoe21s_RNmod, file = "./Data/Side_quests/Hilger/Effort_eoe21s_RNmod.RData")
  save(Effort_eoe22s_RNmod, file = "./Data/Side_quests/Hilger/Effort_eoe22s_RNmod.RData")
  
  #'  Stacked detection histories (all years per species and month)
  save(DH_elk_list, file = "./Data/Side_quests/Hilger/DH_elk_RNmod.RData")
  save(DH_wtd_list, file = "./Data/Side_quests/Hilger/DH_wtd_RNmod.RData")
  
  #'  Seasonal camera station data
  save(station_elk_list, file = "./Data/Side_quests/Hilger/stacked_stations_elk_july_aug_02.03.25.RData")
  save(station_wtd_list, file = "./Data/Side_quests/Hilger/stacked_stations_wtd_july_aug_02.03.25.RData")
  
  save(stations_elk_eoe20s, file = "./Data/Side_quests/Hilger/stations_elk_eoe20s_02.03.25.RData")
  save(stations_elk_eoe21s, file = "./Data/Side_quests/Hilger/stations_elk_eoe21s_02.03.25.RData")
  save(stations_elk_eoe22s, file = "./Data/Side_quests/Hilger/stations_elk_eoe22s_02.03.25.RData")
  save(stations_wtd_eoe20s, file = "./Data/Side_quests/Hilger/stations_wtd_eoe20s_02.03.25.RData")
  save(stations_wtd_eoe21s, file = "./Data/Side_quests/Hilger/stations_wtd_eoe21s_02.03.25.RData")
  save(stations_wtd_eoe22s, file = "./Data/Side_quests/Hilger/stations_wtd_eoe22s_02.03.25.RData")
  
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
                    mean_Tbio = as.numeric(scale(cov$mean_TBio_kg)),
                    max_Tbio = as.numeric(scale(cov$max_TBio_kg)),
                    cv_Tbio = as.numeric(scale(cov$cv_TBio_kg)),
                    mean_HQ = as.numeric(scale(cov$mean_HQBio_kg)),
                    max_HQ = as.numeric(scale(cov$max_HQBio_kg)),
                    cv_HQ = as.numeric(scale(cov$cv_TBio_kg)),
                    total_selected = as.numeric(scale(cov$selected)),
                    total_predicted = as.numeric(scale(cov$total)),
                    prop_selected = as.numeric(scale(cov$prop_selected)))
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
              "mu.lambda", "lambdaYr", "rSetup", "mu.r", "mean.p", 
              "N", "log_N", "loglike.new")
  #'  NOTE about mean vs mu lambda and r: 
  #'  mean.lambda = the intercept based on reference category, i.e., mean lambda for 2020 
  #'  mean.r = the intercept based on reference category, i.e., per-individual detection probability at random sites
  #'  mu.lambda = lambda averaged across all years (and GMUs)
  #'  mu.r = per-individual detection probability averaged across all sites and years 
  
  #'  MCMC settings
  nc <- 3
  ni <- 75000 
  nb <- 5000
  nt <- 10
  na <- 5000
  
  #'  Load competing models
  #'  Univariate models to compare mean & max biomass variables, and selected vs predicted community composition
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_null.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_meanHQ.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_maxHQ.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_cvHQ.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_meanTbio.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_maxTbio.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_cvTbio.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_totalSelected.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_totalPredicted.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_PropSelected.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_HQ_mean.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_HQ_max.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_Tbio_mean.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_selected.propSelected.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_predicted.propSelected.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_elk_july_global_1.R") 
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_elk_july_global_2.R") 
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_elk_aug_global_1.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_elk_aug_global_2.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_july_global_1.R") 
  # source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_july_global_2.R") 
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_aug_global_1.R") 
  # source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_wtd_aug_global_2.R") 
  
  #'  Function to calculate WAIC and joint-likelihood approach (WAICj) based on Gaya & Ketz (2024)
  calc.jointlike <- function(x){
    #'  log-likelihood for every iteration and camera
    like <- as.matrix(x$sims.list$log_N)
    like_joint <- as.matrix(x$sims.list$loglike.new) 
    #'  mean likelihood 
    fbar <- colMeans(exp(like)) 
    fbar_joint <- colMeans(exp(like_joint))
    #'  mean variance in log-likelihood 
    Pw <- sum(apply(like,2,var)) 
    Pw_joint <- sum(apply(like_joint,2,var)) 
    #'  WAIC
    WAIC<- -2*sum(log(fbar))+2*Pw
    #'  joint likelihood WAIC
    WAIC_joint<- -2*sum(log(fbar_joint))+2*Pw_joint
    #'  List WAIC and WAICj
    WAICs <- list(WAIC, WAIC_joint)
    return(WAICs)
  }
  

  #'  -----------------------
  #####  Elk July RN models  #####
  #'  -----------------------
  ######  Null model  ######
  #'  Camera setup on detection  
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_null <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_null$summary[1:15,])
  (RN_elk_july_null_WAICs <- calc.jointlike(RN_elk_july_null))
  print(RN_elk_july_null$DIC)
  which(RN_elk_july_null$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_null$samples)
  save(RN_elk_july_null, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_null.RData")
  
  ######  Univariate models ######  
  #'  Mean HQ 
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_meanHQ <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanHQ.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_meanHQ$summary[1:15,])
  (RN_elk_july_meanHQ_WAICs <- calc.jointlike(RN_elk_july_meanHQ))
  print(RN_elk_july_meanHQ$DIC)
  which(RN_elk_july_meanHQ$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_meanHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_meanHQ.RData")
  
  #'  Max HQ  
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_maxHQ <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxHQ.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_maxHQ$summary[1:15,])
  (RN_elk_july_maxHQ_WAICs <- calc.jointlike(RN_elk_july_maxHQ))
  print(RN_elk_july_maxHQ$DIC)
  which(RN_elk_july_maxHQ$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_maxHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_maxHQ.RData")
  
  #'  cv HQ   
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_cvHQ <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvHQ.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_cvHQ$summary[1:15,])
  (RN_elk_july_cvHQ_WAICs <- calc.jointlike(RN_elk_july_cvHQ))
  print(RN_elk_july_cvHQ$DIC)
  which(RN_elk_july_cvHQ$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_cvHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_cvHQ.RData")
  
  #'  Mean Tbio 
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_meanTbio <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_meanTbio.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_meanTbio$summary[1:15,])
  (RN_elk_july_meanTbio_WAICs <- calc.jointlike(RN_elk_july_meanTbio))
  print(RN_elk_july_meanTbio$DIC)
  which(RN_elk_july_meanTbio$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_meanTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_meanTbio.RData")
  
  #'  Max Tbio 
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_maxTbio <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_maxTbio.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_maxTbio$summary[1:15,])
  (RN_elk_july_maxTbio_WAICs <- calc.jointlike(RN_elk_july_maxTbio))
  print(RN_elk_july_maxTbio$DIC)
  which(RN_elk_july_maxTbio$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_maxTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_maxTbio.RData")
  
  #'  cv Tbio 
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_cvTbio <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                              "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvTbio.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_cvTbio$summary[1:15,])
  (RN_elk_july_cvTbio_WAICs <- calc.jointlike(RN_elk_july_cvTbio))
  print(RN_elk_july_cvTbio$DIC)
  which(RN_elk_july_cvTbio$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_cvTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_cvTbio.RData")
  
  #'  Total Selected 
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_selected <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalSelected.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_selected$summary[1:15,])
  (RN_elk_july_selected_WAICs <- calc.jointlike(RN_elk_july_selected))
  print(RN_elk_july_selected$DIC)
  which(RN_elk_july_selected$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_selected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_selected.RData")
  
  #'  Total Predicted 
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_predicted <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_totalPredicted.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_predicted$summary[1:15,])
  (RN_elk_july_predicted_WAICs <- calc.jointlike(RN_elk_july_predicted))
  print(RN_elk_july_predicted$DIC)
  which(RN_elk_july_predicted$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_predicted, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_predicted.RData")
  
  #'  Proportion Selected
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_propselect <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_PropSelected.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_propselect$summary[1:15,])
  (RN_elk_july_propselect_WAICs <- calc.jointlike(RN_elk_july_propselect))
  print(RN_elk_july_propselect$DIC)
  which(RN_elk_july_propselect$summary[,"Rhat"] > 1.1)
  save(RN_elk_july_propselect, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_propselected.RData")
  
  ######  High Quality Biomass  ######
  #'  Mean & CV  
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_HQ_mean.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                      "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_mean.cv.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_HQ_mean.cv$summary[1:15,])
  (RN_elk_july_HQ_mean.cv_WAICs <- calc.jointlike(RN_elk_july_HQ_mean.cv))
  print(RN_elk_july_HQ_mean.cv$DIC)
  which(RN_elk_july_HQ_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_HQ_mean.cv$samples)
  save(RN_elk_july_HQ_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_HQ_mean.cv.RData")
  
  #'  Max & CV  
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_HQ_max.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                                 "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_HQ_max.cv.txt",
                                 n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                 n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_HQ_max.cv$summary[1:15,])
  (RN_elk_july_HQ_max.cv_WAICs <- calc.jointlike(RN_elk_july_HQ_max.cv))
  print(RN_elk_july_HQ_max.cv$DIC)
  which(RN_elk_july_HQ_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_HQ_max.cv$samples)
  save(RN_elk_july_HQ_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_HQ_max.cv.RData")
  
  ######  Total Biomass  ######
  #'  Mean & CV  
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_Tbio_mean.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_mean.cv.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_Tbio_mean.cv$summary[1:15,])
  (RN_elk_july_Tbio_mean.cv_WAICs <- calc.jointlike(RN_elk_july_Tbio_mean.cv))
  print(RN_elk_july_Tbio_mean.cv$DIC)
  which(RN_elk_july_Tbio_mean.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_Tbio_mean.cv$samples)
  save(RN_elk_july_Tbio_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_Tbio_mean.cv.RData")
  
  #'  Max & CV   
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_Tbio_max.cv <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                                   "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_Tbio_max.cv.txt",
                                   n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                   n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_Tbio_max.cv$summary[1:15,])
  (RN_elk_july_Tbio_max.cv_WAICs <- calc.jointlike(RN_elk_july_Tbio_max.cv))
  print(RN_elk_july_Tbio_max.cv$DIC)
  which(RN_elk_july_Tbio_max.cv$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_Tbio_max.cv$samples)
  save(RN_elk_july_Tbio_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_Tbio_max.cv.RData")
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected  
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_selected.propSelected <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_selected.propSelected.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_selected.propSelected$summary[1:15,])
  (RN_elk_july_selected.propSelected_WAICs <- calc.jointlike(RN_elk_july_selected.propSelected))
  print(RN_elk_july_selected.propSelected$DIC)
  which(RN_elk_july_selected.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_selected.propSelected$samples)
  save(RN_elk_july_selected.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_selected.propSelected.RData")
  
  #'  Predicted & Proportion Selected  
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_predicted.propSelected <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_predicted.propSelected.txt",
                                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_predicted.propSelected$summary[1:15,])
  (RN_elk_july_predicted.propSelected_WAICs <- calc.jointlike(RN_elk_july_predicted.propSelected))
  print(RN_elk_july_predicted.propSelected$DIC)
  which(RN_elk_july_predicted.propSelected$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_predicted.propSelected$samples)
  save(RN_elk_july_predicted.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_predicted.propSelected.RData")
  
  ######  GLOBAL 1  ######    
  #'  Using only most supported non-correlated covariates 
  #'  Global 1: Max HQ, cv HQ, Mean Tbio, Selected, PropSelected 
  #'  Univariate models indicate: Max HQ has lower WAIC compared to Mean HQ, 
  #'  Mean Tbio is statistically meaningful but Max Tbio is not even though lower WAIC, 
  #'  cv HQ and selected had lower WAIC compared to cv Tbio and predicted (all statistically meaningful)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_global1 <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_elk_july_global_1.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_global1$summary[1:15,])
  (RN_elk_july_global1_WAICs <- calc.jointlike(RN_elk_july_global1))
  print(RN_elk_july_global1$DIC)
  which(RN_elk_july_global1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_global1$samples)
  save(RN_elk_july_global1, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_global_maxHQ_cvHQ_meanTbio_sel_propSel.RData")
  
  ######  GLOBAL 2  ######    
  #'  Using only most supported non-correlated covariates 
  #'  Global 2: Mean HQ, cv HQ, Mean Tbio, Selected, PropSelected 
  #'  Univariate models indicate: Mean HQ and Mean Tbio are not correlated at all even though MaxHQ has lower WAIC, 
  #'  Mean Tbio is statistically meaningful but Max Tbio is not even though lower WAIC, 
  #'  cv HQ and selected had lower WAIC compared to cv Tbio and predicted (all statistically meaningful)
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_global2 <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_elk_july_global_2.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_global2$summary[1:15,])
  (RN_elk_july_global2_WAICs <- calc.jointlike(RN_elk_july_global2))
  print(RN_elk_july_global2$DIC)
  which(RN_elk_july_global2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_global2$samples)
  save(RN_elk_july_global2, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_july_mods/RN_elk_july_global_meanHQ_cvHQ_meanTbio_sel_propSel.RData")
  
  #'  -------------------------
  #####  Elk August RN models  #####
  #'  -------------------------
  ######  Null model  ######
  #'  Camera setup on detection  
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_null <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_null$summary[1:15,])
  (RN_elk_aug_null_WAICj <- calc.jointlike(RN_elk_aug_null))
  print(RN_elk_aug_null$DIC)
  which(RN_elk_aug_null$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_null$samples)
  save(RN_elk_aug_null, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_null.RData")
  
  ######  Univariate models ######  
  #'  Mean HQ 
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
  save(RN_elk_aug_meanHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_meanHQ.RData")
  
  #'  Max HQ 
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
  save(RN_elk_aug_maxHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_maxHQ.RData")
  
  #'  cv HQ  
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_cvHQ <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvHQ.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_cvHQ$summary[1:15,])
  (RN_elk_aug_cvHQ_WAICj <- calc.jointlike(RN_elk_aug_cvHQ))
  print(RN_elk_aug_cvHQ$DIC)
  save(RN_elk_aug_cvHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_cvHQ.RData")
  
  #'  Mean Tbio 
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
  save(RN_elk_aug_meanTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_meanTbio.RData")
  
  #'  Max Tbio 
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
  save(RN_elk_aug_maxTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_maxTbio.RData")
  
  #'  cv Tbio 
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_cvTbio <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvTbio.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_cvTbio$summary[1:15,])
  (RN_elk_aug_cvTbio_WAICj <- calc.jointlike(RN_elk_aug_cvTbio))
  print(RN_elk_aug_cvTbio$DIC)
  save(RN_elk_aug_cvTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_cvTbio.RData")
  
  #'  Total Selected 
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
  save(RN_elk_aug_selected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_selected.RData")
  
  #'  Total Predicted 
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
  save(RN_elk_aug_predicted, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_predicted.RData")
  
  #'  Proportion Selected 
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_propselect <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_PropSelected.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_propselect$summary)
  (RN_elk_aug_propselect_WAICj <- calc.jointlike(RN_elk_aug_propselect))
  print(RN_elk_aug_propselect$DIC)
  save(RN_elk_aug_propselect, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_propselected.RData")
  
  ######  High Quality Biomass  ######
  #'  Mean & CV  
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
  save(RN_elk_aug_HQ_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_HQ_mean.cv.RData")
  
  #'  Max & CV  
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
  save(RN_elk_aug_HQ_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_HQ_max.cv.RData")
  
  ######  Total Biomass  ######
  #'  Mean & CV  
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
  save(RN_elk_aug_Tbio_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_Tbio_mean.cv.RData")
  
  #'  Max & CV  
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
  save(RN_elk_aug_Tbio_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_Tbio_max.cv.RData")
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected  
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
  save(RN_elk_aug_selected.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_selected.propSelected.RData")
  
  #'  Predicted & Proportion Selected  
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
  save(RN_elk_aug_predicted.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_predicted.propSelected.RData")
  
  ######  GLOBAL 1  ######    
  #'  Using only most supported non-correlated covariates 
  #'  Global 1: Mean HQ, Max Tbio, cv Tbio, predicted, PropSelected  
  #'  Univariate models indicate: Max HQ has slightly lower WAIC compared to Mean HQ but effect size the same, 
  #'  Max Tbio and predicted had lower WAIC compered to Mean Tbio and selected  
  #'  cv Tbio had lower WAIC compared to cv HQ though neither are statistically meaningful
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_global1 <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_elk_aug_global_1.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_global1$summary[1:15,])
  (RN_elk_aug_global1_WAICj <- calc.jointlike(RN_elk_aug_global1))
  print(RN_elk_aug_global1$DIC)
  which(RN_elk_aug_global1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_global1$samples)
  save(RN_elk_aug_global1, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_global_meanHQ_maxTbio_cvTbio_pred_propSel.RData")
  
  ######  GLOBAL 2  ######    
  #'  Using only most supported non-correlated covariates 
  #'  Global 2: Mean HQ, Mean Tbio, cv Tbio, predicted, PropSelected  
  #'  Univariate models indicate: Mean HQ and Mean Tbio are not correlated at all even though MaxHQ has slightly lower WAIC
  #'  predicted had lower WAIC compered to selected 
  #'  cv Tbio had lower WAIC compared to cv HQ though neither are statistically meaningful
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_global2 <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_elk_aug_global_2.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_global2$summary[1:15,])
  (RN_elk_aug_global2_WAICj <- calc.jointlike(RN_elk_aug_global2))
  print(RN_elk_aug_global2$DIC)
  which(RN_elk_aug_global2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_global2$samples)
  save(RN_elk_aug_global2, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/Elk_aug_mods/RN_elk_aug_global_meanHQ_meanTbio_cvTbio_pred_propSel.RData")
  
  
  #'  -----------------------
  #####  WTD July RN models  #####
  #'  -----------------------
  ######  Null model  ######
  #'  Camera setup on detection  
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_null <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_null.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_null$summary[1:15,])
  (RN_wtd_july_null_WAICj <- calc.jointlike(RN_wtd_july_null))
  print(RN_wtd_july_null$DIC)
  which(RN_wtd_july_null$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_null$samples)
  save(RN_wtd_july_null, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_null.RData")
  
  ######  Univariate models ######  
  #'  Mean HQ 
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
  save(RN_wtd_july_meanHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_meanHQ.RData")
  
  #'  Max HQ  
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
  save(RN_wtd_july_maxHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_maxHQ.RData")
  
  #'  cv HQ  
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_cvHQ <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                            "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvHQ.txt",
                            n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                            n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_cvHQ$summary[1:15,])
  (RN_wtd_july_cvHQ_WAICj <- calc.jointlike(RN_wtd_july_cvHQ))
  print(RN_wtd_july_cvHQ$DIC)
  save(RN_wtd_july_cvHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_cvHQ.RData")
  
  #'  Mean Tbio 
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
  save(RN_wtd_july_meanTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_meanTbio.RData")
  
  #'  Max Tbio 
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
  save(RN_wtd_july_maxTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_maxTbio.RData")
  
  #'  cv Tbio 
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_cvTbio <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                              "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvTbio.txt",
                              n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                              n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_cvTbio$summary[1:15,])
  (RN_wtd_july_cvTbio_WAICj <- calc.jointlike(RN_wtd_july_cvTbio))
  print(RN_wtd_july_cvTbio$DIC)
  save(RN_wtd_july_cvTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_cvTbio.RData")
  
  #'  Total Selected 
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
  save(RN_wtd_july_selected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_selected.RData")
  
  #'  Total Predicted 
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
  save(RN_wtd_july_predicted, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_predicted.RData")
  
  #'  Proportion Selected 
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_propselect <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                                "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_PropSelected.txt",
                                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                                n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_propselect$summary[1:15,])
  (RN_wtd_july_propselect_WAICj <- calc.jointlike(RN_wtd_july_propselect))
  print(RN_wtd_july_propselect$DIC)
  save(RN_wtd_july_propselect, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_propselected.RData")
    
  ######  High Quality Biomass  ######
  #'  Mean & CV  
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
  save(RN_wtd_july_HQ_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_HQ_mean.cv.RData")
  
  #'  Max & CV  
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
  save(RN_wtd_july_HQ_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_HQ_max.cv.RData")
  
  ######  Total Biomass  ######
  #'  Mean & CV  
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
  save(RN_wtd_july_Tbio_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_Tbio_mean.cv.RData")
  
  #'  Max & CV  
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
  save(RN_wtd_july_Tbio_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_Tbio_max.cv.RData")
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected 
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
  save(RN_wtd_july_selected.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_selected.propSelected.RData")
  
  #'  Predicted & Proportion Selected  
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
  save(RN_wtd_july_predicted.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_predicted.propSelected.RData")
  
  ######  GLOBAL  ######    
  #'  Using only most supported non-correlated covariates (based on lowest WAIC from univariate models above)
  #'  Global: Max HQ, cv HQ, Max Tbio, Selected, PropSelected   
  #'  Removed Mean HQ, & Mean Tbio, cv Tbio from 2.24.25 run because correlated with max HQ, cvHQ, & max Tbio with updated forage variables
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
  save(RN_wtd_july_global1, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_global_1_redo.RData")
  
  # start.time = Sys.time()
  # inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  # RN_wtd_july_global2 <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
  #                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_wtd_july_global_2.txt",
  #                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
  #                          n.burnin = nb, parallel = TRUE)
  # end.time <- Sys.time(); (run.time <- end.time - start.time)
  # print(RN_wtd_july_global2$summary)
  # (RN_wtd_july_global2_WAICj <- calc.jointlike(RN_wtd_july_global2))
  # print(RN_wtd_july_global2$DIC)
  # which(RN_wtd_july_global2$summary[,"Rhat"] > 1.1)
  # mcmcplot(RN_wtd_july_global2$samples)
  # save(RN_wtd_july_global2, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_july_mods/RN_wtd_july_global.cvTbio.RData")
  
  #'  -------------------------
  #####  WTD August RN models  #####
  #'  -------------------------
  ######  Null model  ######
  #'  Camera setup on detection  
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
  save(RN_wtd_aug_null, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_null.RData")
  
  ######  Univariate models ######  
  #'  Mean HQ 
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
  save(RN_wtd_aug_meanHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_meanHQ.RData")
  
  #'  Max HQ  
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
  save(RN_wtd_aug_maxHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_maxHQ.RData")
  
  #'  cv HQ  
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_cvHQ <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvHQ.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_cvHQ$summary[1:15,])
  (RN_wtd_aug_cvHQ_WAICj <- calc.jointlike(RN_wtd_aug_cvHQ))
  print(RN_wtd_aug_cvHQ$DIC)
  save(RN_wtd_aug_cvHQ, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_cvHQ.RData")
  
  #'  Mean Tbio 
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
  save(RN_wtd_aug_meanTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_meanTbio.RData")
  
  #'  Max Tbio 
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
  save(RN_wtd_aug_maxTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_maxTbio.RData")

  #'  cv Tbio  
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_cvTbio <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                             "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_cvTbio.txt",
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_cvTbio$summary[1:15,])
  (RN_wtd_aug_cvTbio_WAICj <- calc.jointlike(RN_wtd_aug_cvTbio))
  print(RN_wtd_aug_cvTbio$DIC)
  save(RN_wtd_aug_cvTbio, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_cvTbio.RData")
  
  #'  Total Selected 
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
  save(RN_wtd_aug_selected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_selected.RData")
  
  #'  Total Predicted 
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
  save(RN_wtd_aug_predicted, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_predicted.RData")
  
  #'  Proportion Selected
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_propselect <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                               "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_PropSelected.txt",
                               n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                               n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_propselect$summary[1:15,])
  (RN_wtd_aug_propselect_WAICj <- calc.jointlike(RN_wtd_aug_propselect))
  print(RN_wtd_aug_propselect$DIC)
  save(RN_wtd_aug_propselect, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_propselected.RData")
  
  ######  High Quality Biomass models  ######
  #'  Mean & CV  
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
  save(RN_wtd_aug_HQ_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_HQ_mean.cv.RData")
  
  #'  Max & CV  
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
  save(RN_wtd_aug_HQ_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_HQ_max.cv.RData")
  
  ######  Total Biomass  ######
  #'  Mean & CV  
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
  save(RN_wtd_aug_Tbio_mean.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_Tbio_mean.cv.RData")
  
  #'  Max & CV  
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
  save(RN_wtd_aug_Tbio_max.cv, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_Tbio_max.cv.RData")
  
  ######  Community Composition  ######
  #'  Selected & Proportion Selected  
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
  save(RN_wtd_aug_selected.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_selected.propSelected.RData")
  
  #'  Predicted & Proportion Selected  
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
  save(RN_wtd_aug_predicted.propSelected, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_predicted.propSelected.RData")
  
  ######  GLOBAL  ######   
  #'  Using only most supported non-correlated covariates
  #'  Global: Mean HQ, Max Tbio, cv Tbio, predicted, PropSelected   
  #'  Max HQ, MeanTbio, cv HQ, selected are highly correlated with Mean HQ, MaxTbio, cvTbio, and predicted 
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
  save(RN_wtd_aug_global1, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_global_1_redo.RData")
  
  # start.time = Sys.time()
  # inits_wtd_July <- function(){list(N = ninit_wtd[[2]])}
  # RN_wtd_aug_global2 <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
  #                                 "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_wtd_aug_global_2.txt",
  #                                 n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
  #                                 n.burnin = nb, parallel = TRUE)
  # end.time <- Sys.time(); (run.time <- end.time - start.time)
  # print(RN_wtd_aug_global2$summary)
  # (RN_wtd_aug_global2_WAICj <- calc.jointlike(RN_wtd_aug_global2))
  # print(RN_wtd_aug_global2$DIC)
  # which(RN_wtd_aug_global2$summary[,"Rhat"] > 1.1)
  # mcmcplot(RN_wtd_aug_global2$samples)
  # save(RN_wtd_aug_global2, file = "./Outputs/Hilger_RNmodel/JAGS_out/Fit_02.24.25/WTD_aug_mods/RN_wtd_aug_global.cvTbio.RData")
  # 
