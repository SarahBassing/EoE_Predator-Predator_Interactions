  #'  -------------------------------------------------
  #'  Royle-Nichols Abundance Model for wolf data only
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2026
  #'  ------------------------------------------------
  #'  Script to run RN models to estimate relative abundance at each camera site for wolves.
  #'  
  #'  Detection histories generated in Relative_abundance_Royle-Nochols_model.R
  #'  and camera operations table generated in Detection_data_cleaning.R
  #'  2023 data formatted by M. Parsons using same methods
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
  #'  Detection histories
  #'  2020 - 2022 detection history list order: 
  #'  [[1]] bear_black, [[2]] bobcat, [[3]] coyote, [[4]] mountain_lion, [[5]] wolf, [[6]] elk, [[7]] moose, [[8]] muledeer, [[9]] whitetaileddeer, [[10]] rabbit_hare
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp20s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp21s_RNmod.RData")
  load("./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp22s_RNmod.RData")
  
  #'  2023 detection history list order: 
  #'  [[1]] bear_black, [[2]] coyote, [[3]] mountain_lion, [[4]] wolf, [[5]] elk, [[6]] moose, [[7]] whitetaileddeer
  load("./Data/IDFG camera data/IDFG 2023 detection data/Detection_Histories_RNmodel/DH_smr2023.RData")
  DH_npp23s_RNmod <- DH_smr2023
  # save(DH_npp23s_RNmod, file = "./Data/Relative abundance data/RAI Phase 2/Detection_Histories_RNmodel/DH_npp23s_RNmod.RData")
  
  #'  Problem cameras
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe22s_problem_cams.RData")
  
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
  stations_npp20s <- format_covs(cams_eoe20s, dets = DH_npp20s_RNmod[[1]]) #, effort = effort_20s_RNmod) 
  stations_npp21s <- format_covs(cams_eoe21s, dets = DH_npp21s_RNmod[[1]]) #, effort = effort_21s_RNmod) 
  stations_npp22s <- format_covs(cams_eoe22s, dets = DH_npp22s_RNmod[[1]]) #, effort = effort_22s_RNmod) 
  #'  2023 camera stations
  load("./Data/IDFG camera data/IDFG 2023 detection data/Detection_Histories_RNmodel/stations_smr2023_withcoords.RData")
  stations_npp23s <- stations_smr2023
  # cams_23s <- dplyr::select(stations_smr2023, c("NewLocationID", "Lat", "Long")) %>%
  #   rename("NwLctID" = "NewLocationID")
  # cams_23s_wgs84 <- st_as_sf(cams_23s, coords = c("Long", "Lat"), crs = "+proj=longlat +datum=WGS84 +no_defs")
  # write_sf(cams_23s_wgs84, "./Shapefiles/IDFG spatial data 2023 additional/cams_23s_wgs84.shp")
  
  #'  Double check things are ordered correctly!!!!
  stations_npp20s[82:90,1:4]; DH_npp20s_RNmod[[1]][82:90,1:3]; nrow(stations_npp20s); nrow(DH_npp20s_RNmod[[1]])
  stations_npp21s[82:90,1:4]; DH_npp21s_RNmod[[1]][82:90,1:3]; nrow(stations_npp21s); nrow(DH_npp21s_RNmod[[1]])
  stations_npp22s[82:90,1:4]; DH_npp22s_RNmod[[1]][82:90,1:3]; nrow(stations_npp22s); nrow(DH_npp22s_RNmod[[1]])
  stations_npp23s[82:90,1:4]; DH_npp23s_RNmod[[1]][82:90,1:3]; nrow(stations_npp23s); nrow(DH_npp23s_RNmod[[1]])
  
  #'  -----------------------
  ####  Royle-Nichols model  ####
  #'  -----------------------
  #'  Run Royle-Nichols model using a Bayesian statistical framework in JAGS and
  #'  a likelihood framework with unmarked
  
  #'  ------------------------
  #####  Setup data for JAGS  #####
  #'  ------------------------
  #'  Bundle detection histories and covariates for each species and year
  bundle_dat <- function(dh, cov) {
    #'  Convert detection history to matrix
    dh <- as.matrix(dh)
    dimnames(dh) <- NULL
    #'  Count number of sites per GMU
    ncams_perGMU <- cov %>%
      group_by(GMU) %>%
      summarise(nsites = n()) %>%
      ungroup()
    print(ncams_perGMU)
    
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
  data_JAGS_bundle_20s <- lapply(DH_npp20s_RNmod, bundle_dat, cov = stations_npp20s) 
  data_JAGS_bundle_21s <- lapply(DH_npp21s_RNmod, bundle_dat, cov = stations_npp21s) 
  data_JAGS_bundle_22s <- lapply(DH_npp22s_RNmod, bundle_dat, cov = stations_npp22s) 
  data_JAGS_bundle_23s <- lapply(DH_npp23s_RNmod, bundle_dat, cov = stations_npp23s) 

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
  params <- c("beta0", "beta4", #"beta1", "beta2", "beta3", "beta4", #"mean.lambda", 
              "alpha0", "alpha2", "rSetup", "mu.r", "mean.p", #"alpha1" #"mean.r", 
              "mu.lambda", "totalN", #"occSites", "mean.psi", 
              "N")
  #'  NOTE about mean vs mu lambda and r: 
  #'  mean.lambda = the intercept, i.e., mean lambda for GMU10A 
  #'  mean.r = the intercept, i.e., per-individual detection probability at random sites
  #'  mu.lambda = lambda averaged across all GMUs
  #'  mu.r = per-individual detection probability averaged across all sites 
  
  #'  MCMC settings
  nc <- 3
  ni <- 100000
  nb <- 10000
  nt <- 10
  na <- 5000
  
  #'  -----------------------------------------
  #####  Run RN model with JAGS for wolf data  #####
  #'  -----------------------------------------
  ######  2020  Analyses  ######
  source("./Scripts/Relative_Abundance/RNmodel_JAGS_code_2020.R")  
  
  start.time = Sys.time()
  inits_wolf20s <- function(){list(N = ninit_20s[[5]])}
  RN_wolf_20s <- jags(data_JAGS_bundle_20s[[5]], inits = inits_wolf20s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_20s$summary)
  which(RN_wolf_20s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_20s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_20s$samples)
  save(RN_wolf_20s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_20s_", Sys.Date(), ".RData")) 
  
  
  ######  2021 - 2023  Analyses  ######
  source("./Scripts/Relative_Abundance/RNmodel_JAGS_code.R")
  
  #'  Summer 2021
  start.time = Sys.time()
  inits_wolf21s <- function(){list(N = ninit_21s[[5]])}
  RN_wolf_21s <- jags(data_JAGS_bundle_21s[[5]], inits = inits_wolf21s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_21s$summary)
  which(RN_wolf_21s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_21s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_21s$samples)
  save(RN_wolf_21s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_21s_", Sys.Date(), ".RData"))
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_wolf22s <- function(){list(N = ninit_22s[[5]])}
  RN_wolf_22s <- jags(data_JAGS_bundle_22s[[5]], inits = inits_wolf22s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_22s$summary)
  which(RN_wolf_22s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_22s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_22s$samples)
  save(RN_wolf_22s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_22s_", Sys.Date(), ".RData"))
  
  #'  Summer 2023
  start.time = Sys.time()
  inits_wolf23s <- function(){list(N = ninit_23s[[4]])}   # note the different indexing
  RN_wolf_23s <- jags(data_JAGS_bundle_23s[[4]], inits = inits_wolf23s, params,
                      "./Outputs/Relative_Abundance/RN_model/JAGS_RNmod.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_23s$summary)
  which(RN_wolf_23s$summary[,"Rhat"] < 0.9)
  which(RN_wolf_23s$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_23s$samples)
  save(RN_wolf_23s, file = paste0("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_23s_", Sys.Date(), ".RData"))
  
  
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
  
  #' #'  Load wolf RN model outputs, if needed
  #' load("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_20s_2026-03-24.RData")
  #' load("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_21s_2026-03-24.RData")
  #' load("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_22s_2026-03-25.RData")
  #' load("./Outputs/Relative_Abundance/RN_model/JAGS_out/RN_RAI_wolf_2020_2023/RN_RAI_wolf_23s_2026-03-25.RData")
  
  #'  List RN outputs
  rn_wolf_list <- list(RN_wolf_20s, RN_wolf_21s, RN_wolf_22s, RN_wolf_23s)
    
  #'  List one detection history per year (doesn't matter which species it relates to)
  dh_list <- list(DH_npp20s_RNmod[[1]], DH_npp21s_RNmod[[1]], DH_npp22s_RNmod[[1]], DH_npp23s_RNmod[[1]])
  
  #'  Save estimated N per site
  estimated_N <- function(mod, dh, spp) {
    #'  Grab estimated N and SD per site
    RN.n <- mod$mean$N
    RN.sd <- mod$sd$N
    #'  Grab camera location
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
  rn_wolf_out <- mapply(estimated_N, rn_wolf_list, dh = dh_list, spp = "wolf", SIMPLIFY = FALSE)
  
  #'  List camera spatial data
  cam_list <- list(cams_20s_wgs84, cams_21s_wgs84, cams_22s_wgs84, cams_23s_wgs84)
  
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
  spatial_rn_wolf <- mapply(rn = rn_wolf_out, spatial_rn, spp = "wolf", cams = cam_list, SIMPLIFY = FALSE)
  
  #'  Merge wolf results into single large spatial dataframe and save as a shapefile
  spatial_rn_wolf_locs <- bind_rows(spatial_rn_wolf[[1]], spatial_rn_wolf[[2]], spatial_rn_wolf[[3]], spatial_rn_wolf[[4]])
  st_write(spatial_rn_wolf_locs, "./Shapefiles/IDFG spatial data/Camera_locations/spatial_rn_wolf_locs_2020_2023.shp")
  