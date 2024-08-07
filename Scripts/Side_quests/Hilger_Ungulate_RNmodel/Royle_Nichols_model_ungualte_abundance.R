  #'  ---------------------------
  #'  Local abundance of deer and elk in N. Idaho
  #'  Northern Idaho Predator-Prey Project
  #'  Sarah B. Bassing
  #'  August 2024
  #'  ---------------------------
  #'  Script to format photo-capture data of elk and white-tailed deer in northern 
  #'  Idaho and fit Royle-Nichols abundance models to detection histories to test
  #'  the effects of forage on ungualte relative abundance. Fit RN models per
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
      filter(Species == "elk" | Species == "moose" | Species == "whitetaileddeer") %>%
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
  stations_eoe20s <- format_covs(DH_eoe20s_RNmod[[1]], cams = sites_20s, sitecovs = sitecovs_20s) %>%
    dplyr::select(-ID)
  stations_eoe21s <- format_covs(DH_eoe21s_RNmod[[1]], cams = sites_21s, sitecovs = sitecovs_21s) %>%
    dplyr::select(-ID)
  stations_eoe22s <- format_covs(DH_eoe22s_RNmod[[1]], cams = sites_22s, sitecovs = sitecovs_22s) %>%
    dplyr::select(-ID)
  
  #'  Stack together 
  station_stats <- rbind(stations_eoe20s, stations_eoe21s, stations_eoe22s)
  
  #'  Double check things are ordered correctly!!!!
  stations_eoe20s[82:90,]; DH_eoe20s_RNmod[[1]][82:90,1:3]; nrow(stations_eoe20s); nrow(DH_eoe20s_RNmod[[1]])
  stations_eoe21s[82:90,]; DH_eoe21s_RNmod[[1]][82:90,1:3]; nrow(stations_eoe21s); nrow(DH_eoe21s_RNmod[[1]])
  stations_eoe22s[82:90,]; DH_eoe22s_RNmod[[1]][82:90,1:3]; nrow(stations_eoe22s); nrow(DH_eoe22s_RNmod[[1]])
  
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
  save(stations_eoe20s, file = "./Data/Side_quests/Hilger/stations_eoe20s.RData")
  save(stations_eoe21s, file = "./Data/Side_quests/Hilger/stations_eoe21s.RData")
  save(stations_eoe22s, file = "./Data/Side_quests/Hilger/stations_eoe22s.RData")
  
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
                    forest = as.numeric(scale(cov$perc_forest)))
    str(bundled)
    return(bundled)
  }
  list_names <- c("July_DH", "August_DH")
  data_JAGS_bundle_elk <- lapply(DH_elk_list, bundle_dat, cov = station_stats)
  data_JAGS_bundle_wtd <- lapply(DH_wtd_list, bundle_dat, cov = station_stats)
  
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
  params <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "alpha0", "alpha1", 
              "rSetup", "mu.r", "mean.p", "mu.lambda", "totalN", "occSites", "mean.psi", "N")
  #'  NOTE about mean vs mu lambda and r: 
  #'  mean.lambda = the intercept based on reference category, i.e., mean lambda for GMU10A 
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
  #'  mod1: null (with year effect)
  #'  mod2: GMU and basic landscape characteristics
  #'  mod3: GMU, landscape characteristics, and forage availability
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_mod1.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_mod2.R")
  source("./Scripts/Side_quests/Hilger_Ungulate_RNmodel/RNmodel_JAGS_code_mod3.R")
  
  #'  -----------------------
  #####  Elk July RN models  #####
  #'  -----------------------
  ######  Model 1  ######
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_mod1 <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                      "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod1.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_mod1$summary)
  which(RN_elk_july_mod1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_mod1$samples)
  save(RN_elk_july_mod1, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_mod1_", Sys.Date(), ".RData"))
  
  ######  Model 2  ######
  start.time = Sys.time()
  inits_elk_July <- function(){list(N = ninit_elk[[1]])}
  RN_elk_july_mod2 <- jags(data_JAGS_bundle_elk[[1]], inits = inits_elk_July, params,
                      "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod2.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_july_mod2$summary)
  which(RN_elk_july_mod2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_july_mod2$samples)
  save(RN_elk_july_mod2, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_july_mod2_", Sys.Date(), ".RData"))
  
  #'  -------------------------
  #####  Elk August RN models  #####
  #'  -------------------------
  ######  Model 1  ######
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_mod1 <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                      "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod1.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_mod1$summary)
  which(RN_elk_aug_mod1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_mod1$samples)
  save(RN_elk_aug_mod1, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_mod1_", Sys.Date(), ".RData"))
  
  ######  Model 2  ######
  start.time = Sys.time()
  inits_elk_Aug <- function(){list(N = ninit_elk[[2]])}
  RN_elk_aug_mod2 <- jags(data_JAGS_bundle_elk[[2]], inits = inits_elk_Aug, params,
                     "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod2.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_aug_mod2$summary)
  which(RN_elk_aug_mod2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_aug_mod2$samples)
  save(RN_elk_aug_mod2, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_elk_aug_mod2_", Sys.Date(), ".RData"))
  
  #'  -----------------------
  #####  WTD July RN models  #####
  #'  -----------------------
  ######  Model 1  ######
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_mod1 <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod1.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_mod1$summary)
  which(RN_wtd_july_mod1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_mod1$samples)
  save(RN_wtd_july_mod1, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_mod1_", Sys.Date(), ".RData"))
  
  ######  Model 2  ######
  start.time = Sys.time()
  inits_wtd_July <- function(){list(N = ninit_wtd[[1]])}
  RN_wtd_july_mod2 <- jags(data_JAGS_bundle_wtd[[1]], inits = inits_wtd_July, params,
                           "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod2.txt",
                           n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                           n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_july_mod2$summary)
  which(RN_wtd_july_mod2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_july_mod2$samples)
  save(RN_wtd_july_mod2, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_july_mod2_", Sys.Date(), ".RData"))
  
  #'  -------------------------
  #####  WTD August RN models  #####
  #'  -------------------------
  ######  Model 1  ######
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_mod1 <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod1.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_mod1$summary)
  which(RN_wtd_aug_mod1$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_mod1$samples)
  save(RN_wtd_aug_mod1, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_mod1_", Sys.Date(), ".RData"))
  
  ######  Model 2  ######
  start.time = Sys.time()
  inits_wtd_Aug <- function(){list(N = ninit_wtd[[2]])}
  RN_wtd_aug_mod2 <- jags(data_JAGS_bundle_wtd[[2]], inits = inits_wtd_Aug, params,
                          "./Outputs/Hilger_RNmodel/RNmodel_JAGS_code_mod2.txt",
                          n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                          n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_aug_mod2$summary)
  which(RN_wtd_aug_mod2$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wtd_aug_mod2$samples)
  save(RN_wtd_aug_mod2, file = paste0("./Outputs/Hilger_RNmodel/JAGS_out/RN_wtd_aug_mod2_", Sys.Date(), ".RData"))
  
  
  
  
  
