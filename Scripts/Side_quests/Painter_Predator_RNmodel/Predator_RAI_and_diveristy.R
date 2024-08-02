  #'  ---------------------------
  #'  Detection-corrected relative abundance indices for predators in N. Idaho
  #'  Northern Idaho Predator-Prey Project
  #'  Sarah B. Bassing
  #'  July 2024
  #'  ---------------------------
  #'  Script to format photo-capture data of 5 predator species in northern Idaho 
  #'  and fit Royle-Nichols abundance models to detection histories. Predictions
  #'  of "local" abundance at each camera site can be used as a detection-corrected
  #'  index of relative abundance and used to calculate metrics of species diversity.
  #'  
  #'  Input data: Species-specific detection data (June, July, Aug, 2020, 2021, & 2022)
  #'  from camera traps and camera site data regarding camera setup (i.e., predator-style 
  #'  camera setup vs ungulate-style camera setup). 
  #'  ---------------------------
  
  #'  Clear workspace
  rm(list = ls())

  #'  Load libraries
  library(camtrapR)
  library(jagsUI)
  library(chron)
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
  
  #'  Load basic camera site covariate data
  cams_eoe_long <- read.csv("./Data/side_quests/Painter/cams_eoe_long_Smr2020-2022.csv") %>%
    dplyr::select(c("NewLocationID", "Gmu", "Setup", "Season")) %>%
    distinct()
  
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
      filter(Species == "bear_black" | Species == "bobcat" | Species == "coyote" |
               Species == "mountain_lion" | Species == "wolf") %>%
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
  df_all_20s <- thin_detections(eoe20s_allM, seqprobs = eoe_seqprob_20s, start_date = "2020-06-01", end_date = "2020-08-31") 
  df_all_21s <- thin_detections(eoe21s_allM, seqprobs = eoe_seqprob_21s, start_date = "2021-06-01", end_date = "2021-08-31")
  df_all_22s <- thin_detections(eoe22s_allM, seqprobs = eoe_seqprob_22s, start_date = "2022-06-01", end_date = "2022-08-31") 
  
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
  #'  FYI: June 1 - Aug 31 = 92 1-day sampling periods
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
    
    #'  Remove rows where cameras was inoperable for >61 days (we require at least 
    #'  30 days of active surveying for site to be included, i.e., <62 NAs)
    short_dh <- filter(short_dh, rowSums(is.na(short_dh)) < 62)
    short_effort <- short_effort[short_effort$NewLocationID %in% short_dh$NewLocationID,]
    
    #'  Drop NewLocationID column
    short_dh <- dplyr::select(short_dh, -NewLocationID)
    short_effort <- dplyr::select(short_effort, -NewLocationID)
    
    dh_list <- list(short_dh, short_effort)
    
    return(dh_list)
  }
  #'  Create season-specific detection histories for species listed below
  spp_smr <- list("bear_black", "bobcat", "coyote", "mountain_lion", "wolf")
  DHeff_eoe20s_RNmod <- lapply(spp_smr, DH, dets = eoe20s_det_events, cam_probs = eoe20s_probs, start_date = "2020-06-01", y = "binary", oc = 92) 
  DHeff_eoe21s_RNmod <- lapply(spp_smr, DH, dets = eoe21s_det_events, cam_probs = eoe21s_probs, start_date = "2021-06-01", y = "binary", oc = 92)
  DHeff_eoe22s_RNmod <- lapply(spp_smr, DH, dets = eoe22s_det_events, cam_probs = eoe22s_probs, start_date = "2022-06-01", y = "binary", oc = 92)
  
  #'  Remove sampling effort from species-specific DH_npp lists
  strip_list <- function(dh) {
    #'  Keep only the detection history per species
    dh_only <- dh[[1]]
    return(dh_only)
  }
  DH_eoe20s_RNmod <- lapply(DHeff_eoe20s_RNmod, strip_list)
  DH_eoe21s_RNmod <- lapply(DHeff_eoe21s_RNmod, strip_list)
  DH_eoe22s_RNmod <- lapply(DHeff_eoe22s_RNmod, strip_list)  
  
  #####  5) Format covariate data  #####
  #'  -----------------------------
  #'  Break up covariate data by summer
  covs_20s <- filter(cams_eoe_long, Season == "Smr20")
  covs_21s <- filter(cams_eoe_long, Season == "Smr21")
  covs_22s <- filter(cams_eoe_long, Season == "Smr22") %>%
    distinct()
  
  #'  Scale and format site-level covariates
  format_covs <- function(cams, dets) {   
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
      #'  Arrange by camera location -- NECESSARY TO MATCH DH's CAMERA LOCATION ORDER
      arrange(factor(NewLocationID, levels = camsites))
      
    return(formatted)
  }
  stations_eoe20s <- format_covs(covs_20s, dets = DH_eoe20s_RNmod[[1]]) 
  stations_eoe21s <- format_covs(covs_21s, dets = DH_eoe21s_RNmod[[1]]) 
  stations_eoe22s <- format_covs(covs_22s, dets = DH_eoe22s_RNmod[[1]])
  
  #'  Double check things are ordered correctly!!!!
  stations_eoe20s[82:90,1:4]; DH_eoe20s_RNmod[[1]][82:90,1:3]; nrow(stations_eoe20s); nrow(DH_eoe20s_RNmod[[1]])
  stations_eoe21s[82:90,1:4]; DH_eoe21s_RNmod[[1]][82:90,1:3]; nrow(stations_eoe21s); nrow(DH_eoe21s_RNmod[[1]])
  stations_eoe22s[82:90,1:4]; DH_eoe22s_RNmod[[1]][82:90,1:3]; nrow(stations_eoe22s); nrow(DH_eoe22s_RNmod[[1]])

  #####  Save!  #####
  #'  ----------
  #'  Seasonal detection histories
  save(DH_eoe20s_RNmod, file = "./Data/Side_quests/Painter/DH_eoe20s_RNmod.RData")
  save(DH_eoe21s_RNmod, file = "./Data/Side_quests/Painter/DH_eoe21s_RNmod.RData")
  save(DH_eoe22s_RNmod, file = "./Data/Side_quests/Painter/DH_eoe22s_RNmod.RData")
  
  #'  Seasonal camera station data
  save(stations_eoe20s, file = "./Data/Side_quests/Painter/stations_eoe20s.RData")
  save(stations_eoe21s, file = "./Data/Side_quests/Painter/stations_eoe21s.RData")
  save(stations_eoe22s, file = "./Data/Side_quests/Painter/stations_eoe22s.RData")
  
  
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
                    ngmu = max(as.numeric(cov$GMU)),
                    nsets = max(as.numeric(cov$Setup)),
                    gmu = as.numeric(cov$GMU), 
                    setup = as.numeric(cov$Setup))
    str(bundled)
    return(bundled)
  }
  data_JAGS_bundle_20s <- lapply(DH_eoe20s_RNmod, bundle_dat, cov = stations_eoe20s)
  data_JAGS_bundle_21s <- lapply(DH_eoe21s_RNmod, bundle_dat, cov = stations_eoe21s)
  data_JAGS_bundle_22s <- lapply(DH_eoe22s_RNmod, bundle_dat, cov = stations_eoe22s)

  #'  Initial values
  #'  Using naive occupancy as a starting point for local abundance
  initial_n <- function(dh) {
    #'  Max value per row
    ninit <- apply(dh, 1, max, na.rm = TRUE)
    ninit <- as.vector(ninit)
    return(ninit)
  }
  #'  Apply function per species for each year
  ninit_20s <- lapply(DH_eoe20s_RNmod, initial_n)
  ninit_21s <- lapply(DH_eoe21s_RNmod, initial_n)
  ninit_22s <- lapply(DH_eoe22s_RNmod, initial_n)
  
  #'  Parameters monitored
  params <- c("beta0", "beta1", "beta2", "alpha0", "alpha1", 
              "rSetup", "mu.r", "mean.p", "mu.lambda", 
              "totalN", "occSites", "mean.psi", "N")
  #'  NOTE about mean vs mu lambda and r: 
  #'  mean.lambda = the intercept based on reference category, i.e., mean lambda for GMU10A 
  #'  mean.r = the intercept based on reference category, i.e., per-individual detection probability at random sites
  #'  mu.lambda = lambda averaged across all GMUs
  #'  mu.r = per-individual detection probability averaged across all sites 
  
  #'  MCMC settings
  nc <- 3
  ni <- 100000
  nb <- 20000
  nt <- 10
  na <- 5000
  

  #'  --------------------
  ######  2020  Analyses  ######
  #'  --------------------
  source("./Scripts/Side_quests/Painter_Predator_RNmodel/RNmodel_JAGS_code_2020_WTD_FawnProject.R")
  
  #'  BLACK BEAR JUNE-AUG 2020
  start.time = Sys.time()
  inits_bear20s <- function(){list(N = ninit_20s[[1]])}
  RN_bear_20s <- jags(data_JAGS_bundle_20s[[1]], inits = inits_bear20s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_20s$summary)
  which(RN_bear_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_20s$samples)
  save(RN_bear_20s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_bear_20s_", Sys.Date(), ".RData"))
  
  #'  BOBCAT JUNE-AUG  2020
  start.time = Sys.time()
  inits_bob20s <- function(){list(N = ninit_20s[[2]])}
  RN_bob_20s <- jags(data_JAGS_bundle_20s[[2]], inits = inits_bob20s, params,
                     "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2020.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bob_20s$summary)
  which(RN_bob_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bob_20s$samples)
  save(RN_bob_20s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_bob_20s_", Sys.Date(), ".RData")) 
  
  #'  COYOTE JUNE-AUG 2020
  start.time = Sys.time()
  inits_coy20s <- function(){list(N = ninit_20s[[3]])}
  RN_coy_20s <- jags(data_JAGS_bundle_20s[[3]], inits = inits_coy20s, params,
                     "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2020.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_20s$summary)
  which(RN_coy_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_20s$samples)
  save(RN_coy_20s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_coy_20s_", Sys.Date(), ".RData"))
  
  #'  MOUNTAIN LION JUNE-AUG 2020
  start.time = Sys.time()
  inits_lion20s <- function(){list(N = ninit_20s[[4]])}
  RN_lion_20s <- jags(data_JAGS_bundle_20s[[4]], inits = inits_lion20s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_20s$summary)
  which(RN_lion_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_20s$samples)
  save(RN_lion_20s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_lion_20s_", Sys.Date(), ".RData")) 
  
  #'  WOLF JUNE-AUG 2020
  start.time = Sys.time()
  inits_wolf20s <- function(){list(N = ninit_20s[[5]])}
  RN_wolf_20s <- jags(data_JAGS_bundle_20s[[5]], inits = inits_wolf20s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_20s$summary)
  which(RN_wolf_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_20s$samples)
  save(RN_wolf_20s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_wolf_20s_", Sys.Date(), ".RData")) 
  
  #'  ---------------------------
  ######  2021 & 2022  Analyses  ######
  #'  ---------------------------
  source("./Scripts/Side_quests/Painter_Predator_RNmodel/RNmodel_JAGS_code_2021&2022_WTD_FawnProject.R")
  
  #'  BLACK BEAR JUNE-AUG 2021
  start.time = Sys.time()
  inits_bear21s <- function(){list(N = ninit_21s[[1]])}
  RN_bear_21s <- jags(data_JAGS_bundle_21s[[1]], inits = inits_bear21s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_21s$summary)
  which(RN_bear_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_21s$samples)
  save(RN_bear_21s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_bear_21s_", Sys.Date(), ".RData"))
  
  #'  BLACK BEAR JUNE-AUG 2022
  start.time = Sys.time()
  inits_bear22s <- function(){list(N = ninit_22s[[1]])}
  RN_bear_22s <- jags(data_JAGS_bundle_22s[[1]], inits = inits_bear22s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_22s$summary)
  which(RN_bear_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_22s$samples)
  save(RN_bear_22s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_bear_22s_", Sys.Date(), ".RData"))
  
  #'  BOBCAT JUNE-AUG 2021
  start.time = Sys.time()
  inits_bob21s <- function(){list(N = ninit_21s[[2]])}
  RN_bob_21s <- jags(data_JAGS_bundle_21s[[2]], inits = inits_bob21s, params,
                     "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bob_21s$summary)
  which(RN_bob_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bob_21s$samples)
  save(RN_bob_21s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_bob_21s_", Sys.Date(), ".RData"))
  
  #'  BOBCAT JUNE-AUG 2022
  start.time = Sys.time()
  inits_bob22s <- function(){list(N = ninit_22s[[2]])}
  RN_bob_22s <- jags(data_JAGS_bundle_22s[[2]], inits = inits_bob22s, params,
                     "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bob_22s$summary)
  which(RN_bob_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bob_22s$samples)
  save(RN_bob_22s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_bob_22s_", Sys.Date(), ".RData"))
  
  #'  COYOTE JUNE-AUG 2021
  start.time = Sys.time()
  inits_coy21s <- function(){list(N = ninit_21s[[3]])}
  RN_coy_21s <- jags(data_JAGS_bundle_21s[[3]], inits = inits_coy21s, params,
                     "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_21s$summary)
  which(RN_coy_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_21s$samples)
  save(RN_coy_21s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_coy_21s_", Sys.Date(), ".RData"))
  
  #'  COYOTE JUNE-AUG 2022
  start.time = Sys.time()
  inits_coy22s <- function(){list(N = ninit_22s[[3]])}
  RN_coy_22s <- jags(data_JAGS_bundle_22s[[3]], inits = inits_coy22s, params,
                     "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_22s$summary)
  which(RN_coy_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_22s$samples)
  save(RN_coy_22s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_coy_22s_", Sys.Date(), ".RData"))
  
  #'  MOUNTAIN LION JUNE-AUG 2021
  start.time = Sys.time()
  inits_lion21s <- function(){list(N = ninit_21s[[4]])}
  RN_lion_21s <- jags(data_JAGS_bundle_21s[[4]], inits = inits_lion21s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_21s$summary)
  which(RN_lion_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_21s$samples)
  save(RN_lion_21s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_lion_21s_", Sys.Date(), ".RData"))
  
  #'  MOUNTAIN LION JUNE-AUG 2022
  start.time = Sys.time()
  inits_lion22s <- function(){list(N = ninit_22s[[4]])}
  RN_lion_22s <- jags(data_JAGS_bundle_22s[[4]], inits = inits_lion22s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_22s$summary)
  which(RN_lion_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_22s$samples)
  save(RN_lion_22s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_lion_22s_", Sys.Date(), ".RData"))
  
  #'  WOLF JUNE-AUG 2021
  # ni_wolf <- 100000 
  start.time = Sys.time()
  inits_wolf21s <- function(){list(N = ninit_21s[[5]])}
  RN_wolf_21s <- jags(data_JAGS_bundle_21s[[5]], inits = inits_wolf21s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_21s$summary)
  which(RN_wolf_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_21s$samples)
  save(RN_wolf_21s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_wolf_21s_", Sys.Date(), ".RData"))
  
  #'  WOLF JUNE-AUG 2022
  # ni_wolf <- 100000 
  start.time = Sys.time()
  inits_wolf22s <- function(){list(N = ninit_22s[[5]])}
  RN_wolf_22s <- jags(data_JAGS_bundle_22s[[5]], inits = inits_wolf22s, params,
                      "./Outputs/Painter_RNmodel/RNmodel_JAGS_code_2021&2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_22s$summary)
  which(RN_wolf_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_22s$samples)
  save(RN_wolf_22s, file = paste0("./Outputs/Painter_RNmodel/JAGS_out/RN_wolf_22s_", Sys.Date(), ".RData"))
  
  #'  List and save model outputs together 
  rn_bear_list <- list(RN_bear_20s, RN_bear_21s, RN_bear_22s)
  rn_bob_list <- list(RN_bob_20s, RN_bob_21s, RN_bob_22s)
  rn_coy_list <- list(RN_coy_20s, RN_coy_21s, RN_coy_22s)
  rn_lion_list <- list(RN_lion_20s, RN_lion_21s, RN_lion_22s)
  rn_wolf_list <- list(RN_wolf_20s, RN_wolf_21s, RN_wolf_22s)
  
  RN_model_outputs <- list(rn_bear_list, rn_bob_list, rn_coy_list, rn_lion_list, rn_wolf_list)
  
  save(RN_model_outputs, file = "./Outputs/Painter_RNmodel/RN_model_outputs.RData")
  
  
  #'  -----------------------------
  ####  Summarize model estimates  ####
  #'  -----------------------------
  #'  Load saved detection histories
  load("./Data/Side_quests/Painter/DH_eoe20s_RNmod.RData")
  load("./Data/Side_quests/Painter/DH_eoe21s_RNmod.RData")
  load("./Data/Side_quests/Painter/DH_eoe22s_RNmod.RData")
  
  #'  List one detection history per year (doesn't matter which species it relates to)
  dh_list <- list(DH_eoe20s_RNmod[[1]], DH_eoe21s_RNmod[[1]], DH_eoe22s_RNmod[[1]])
  
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
  
  #'  Create giant data frames of site-specific local abundance estimates
  rn_2020 <- rbind(rn_bear_out[[1]], rn_bob_out[[1]], rn_coy_out[[1]], rn_lion_out[[1]], rn_wolf_out[[1]]) %>%
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr20") %>%
    relocate(season, .after = Setup)
  rn_2021 <- rbind(rn_bear_out[[2]], rn_bob_out[[2]], rn_coy_out[[2]], rn_lion_out[[2]], rn_wolf_out[[2]]) %>%
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr21") %>%
    relocate(season, .after = Setup)
  rn_2022 <- rbind(rn_bear_out[[3]], rn_bob_out[[3]], rn_coy_out[[3]], rn_lion_out[[3]], rn_wolf_out[[3]]) %>%
    arrange(NewLocationID, Species) %>%
    mutate(season = "Smr22") %>%
    relocate(season, .after = Setup)
  
  RN_abundance <- list(rn_2020, rn_2021, rn_2022)
  RN_abundance_df <- rbind(rn_2020, rn_2021, rn_2022)
  
  save(RN_abundance, file = "./Outputs/Painter_RNmodel/RN_abundance.RData")
  write_csv(RN_abundance_df, file = "./Outputs/Painter_RNmodel/RN_abundance.csv")
  
  #'  -----------------------------
  ####  Species diversity metrics  ####
  #'  -----------------------------
  #'  Function to calculate species richness and Shannon's diversity index based 
  #'  on species detected at each camera
  species_diversity <- function(RA, season) {
    #'  Species Richness (S)
    #'  --------------------
    #'  Sum number of unique species detected at each camera
    SR <- RA %>% 
      dplyr::select(c(NewLocationID, Setup, Species, RN.n)) %>%
      mutate(Species_det = ifelse(RN.n < 1, 0, 1)) %>%
      group_by(NewLocationID, Setup) %>%
      mutate(SR = sum(Species_det)) %>%
      dplyr::select(-c(Species, RN.n, Species_det)) %>%
      slice(1L) %>%
      ungroup()
    
    #'  Shannon's diversity index (H)
    #'  -----------------------------
    #'  Considers species richness and evenness (abundance of each species)
    #'  https://www.programmingr.com/shannon-diversity-index-the-diversity-function-in-r/
    Shannon <- as.data.frame(RA) %>%
      dplyr::select(c(NewLocationID, Setup, Species, RN.n)) %>%
      #'  Force species with extremely low predicted relative abundance to 0 (RN < 1)
      mutate(RN.n = ifelse(RN.n < 1, 0, RN.n)) %>%
      pivot_wider(names_from = Species, values_from = RN.n) %>%
      as.data.frame(.)
    
    #'  Loop through each camera site to calculate H
    H <- c(NA)
    for(i in 1:nrow(Shannon)) {
      #'  Relative abundance of each species
      n <- c(Shannon[i,3], Shannon[i,4], Shannon[i,5], Shannon[i,6], Shannon[i,7]) 
      #'  Remove species that were not detected (RA = 0)
      n <- n[n != 0]
      #'  Calculate proportion of community each species represents
      N <- sum(n)
      p <- n/N
      #'  Calculate Shannon's diversity index (H)
      H[i] <- -sum(p * log(p))
    }
    Shannon <- cbind(Shannon, H) %>%
      mutate(H = round(H, 5))
    
    #'  Merge SR, H, and most frequently detected wild ungulate species with raw and 
    #'  weighted relative abundance indices
    Spp_diversity <- full_join(SR, Shannon, by = c("NewLocationID", "Setup")) %>%
      relocate(SR, .after = Setup) %>%
      relocate(H, .after = SR) %>%
      mutate(season = season) %>%
      relocate(season, .after = Setup)
    
    return(Spp_diversity)
  }
  spp_diversity_Smr20 <- species_diversity(rn_2020, season = "Smr20") 
  spp_diversity_Smr21 <- species_diversity(rn_2021, season = "Smr21")
  spp_diversity_Smr22 <- species_diversity(rn_2022, season = "Smr22")
  
  spp_diversity <- rbind(spp_diversity_Smr20, spp_diversity_Smr21, spp_diversity_Smr22) %>%
    dplyr::select(c(NewLocationID, Setup, season, SR, H))
  spp_diversity_list <- list(spp_diversity_Smr20, spp_diversity_Smr21, spp_diversity_Smr22)
  
  save(spp_diversity, file = "./Outputs/Painter_RNmodel/spp_diversity.RData")
  write_csv(spp_diversity, file = "./Outputs/Painter_RNmodel/spp_diversity.csv")
  
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
  
  spatial_diversity <- function(div, spp, cams) {
    spatial_div <- div %>%
      #'  Rename camera location column to match spatial data
      rename("NwLctID" = "NewLocationID") 
    
    #'  Join spatial data with diversity data
    H_shp <- full_join(cams, spatial_div, by = "NwLctID") %>%
      filter(!is.na(H)) %>%
      #'  Change camera location column back to something less awkward
      rename("NewLocationID" = "NwLctID")
    
    return(H_shp)
  }
  spatial_spp_diversity <- mapply(div = spp_diversity_list, spatial_diversity, spp = "species_diversity", cams = cam_list, SIMPLIFY = FALSE)
  
  #'  List spatial RN abundance data and save
  spatial_Painter_RNmodel_list <- list(spatial_rn_bear, spatial_rn_bob, spatial_rn_coy, spatial_rn_lion, spatial_rn_wolf, spatial_spp_diversity)
  save(spatial_Painter_RNmodel_list, file = "./Shapefiles/IDFG spatial data/Camera_locations/spatial_Painter_RNmodel_list.RData")
  
  year_list <- list("2020", "2021", "2022")
  
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
  spp_div_all <- mapply(add_yr, dat = spatial_spp_diversity, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  #'  Make one giant faceted plot where rows represent GMU and columns represent years 
  #'  to ensure that the dot sizes are all consistent for at least a single species
  map_rn <- function(sf_rn, spp) {
    #'  Define size of circles
    size_breaks <- c(0, 1, 2, 3, 5, 7, 9, 12)
    
    sf_rn <- mutate(sf_rn, GMU = factor(GMU, levels = c("GMU1", "GMU6", "GMU10A")))
    pal <- c("darkcyan", "lightcoral", "darkgoldenrod3")
    
    #'  Create figure
    spp_rn <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) + 
      geom_sf(data = sf_rn, aes(size = RN.n.rounded, colour = GMU, fill = GMU), shape = 21, alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,12)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      labs(size = "Estimated \nlocal abundance", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 18)) +
      facet_wrap(~Year) + 
      labs(title = paste("Estimated local abundance of", spp, "from RN model, rounded to whole number"))
    
    #'  Plot each map
    print(spp_rn)
    
    return(spp_rn)
  }
  rn_maps_bear <- map_rn(rn_bear_all, spp = "black bear")
  rn_maps_bob <- map_rn(rn_bob_all, spp = "bobcat")
  rn_maps_coy <- map_rn(rn_coy_all, spp = "coyote")
  rn_maps_lion <- map_rn(rn_lion_all, spp = "mountain lion")
  rn_maps_wolf <- map_rn(rn_wolf_all, spp = "wolf")
  
  ggsave("./Outputs/Painter_RNmodel/Figures/RN_map_blackbear.tiff", rn_maps_bear,
         units = "in", width = 13, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Painter_RNmodel/Figures/RN_map_bobcat.tiff", rn_maps_bob,
         units = "in", width = 13, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Painter_RNmodel/Figures/RN_map_coyote.tiff", rn_maps_coy,
         units = "in", width = 13, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Painter_RNmodel/Figures/RN_map_lion.tiff", rn_maps_lion,
         units = "in", width = 13, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Painter_RNmodel/Figures/RN_map_wolf.tiff", rn_maps_wolf,
         units = "in", width = 13, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  
  map_diversity <- function(sf_div, div_metric, div_type, size_breaks) {
    sf_div <- sf_div %>%
      separate(NewLocationID, c("GMU", "Setup", "site"), sep = "_") %>% 
      mutate(GMU = factor(GMU, levels = c("GMU1", "GMU6", "GMU10A"))) %>%
      dplyr::select(-c(Setup, site))
    pal <- c("darkcyan", "lightcoral", "darkgoldenrod3")
    
    #'  Create figure
    spp_div <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) + 
      geom_sf(data = sf_div, aes(size = div_metric, colour = GMU, fill = GMU), shape = 21, alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,5)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      labs(size = "Predator diversity", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 18)) +
      facet_wrap(~Year) + 
      labs(title = paste0("Predator diversity (", div_type, ")"))
    
    #'  Plot each map
    print(spp_div)
    
    return(spp_div)
  }
  diversity_maps_H <- map_diversity(spp_div_all, div_metric = spp_div_all$H, div_type = "Shannon's H", size_breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5))
  diversity_maps_SR <- map_diversity(spp_div_all, div_metric = spp_div_all$SR, div_type = "Species Richness", size_breaks = c(0, 1, 2, 3, 4, 5))
  
  
  ggsave("./Outputs/Painter_RNmodel/Figures/Diversity_map_ShannonsH.tiff", diversity_maps_H,
         units = "in", width = 13, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Painter_RNmodel/Figures/Diversity_map_SpeciesRichness.tiff", diversity_maps_SR,
         units = "in", width = 13, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  