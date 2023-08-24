  #'  --------------------------------------
  #'  Extract covariates at camera locations
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  August 2023
  #'  --------------------------------------
  #'  Read in spatial data and extract covariate values at camera locations.
  #'  
  #'  Requires:
  #'  Problem camera data frames from Detection_data_cleaning.R
  #'  Wolf detection and minimum group count data from Detection_histories_for_occmod.R
  #'  Reported mortality data from Covariates_NLCD_manipulations.R
  #'  And spatial layers
  #'  --------------------------------------
  
  #'  Load libraries
  library(sf)
  library(rgeos)
  library(terra)
  library(tidyverse)
  library(ggplot2)
  library(grid)
  library(MASS)
  
  #'  ----------------
  ####  Spatial data  ####
  #'  ----------------
  #'  Load and review
  pforest <- rast("./Shapefiles/National Land Cover Database (NCLD)/PercentForest_500m.tif")
  nlcd <- rast("./Shapefiles/National Land Cover Database (NCLD)/NLCD19_Idaho.tif")
  id <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp")
  elev <- rast("./Shapefiles/IDFG spatial data/Elevation__10m2.tif")
  habclass <- rast("./Shapefiles/IDFG spatial data/HabLayer_30m2.tif")
  dist2suburbs <- rast("./Shapefiles/GEE/HumanSettlement/Dist2Suburbs.tif")
  dist2rural <- rast("./Shapefiles/GEE/HumanSettlement/Dist2Rural.tif")
  #'  IDFG Geodatabase with roads
  idfg_gdb <- "./Shapefiles/IDFG spatial data/IDFG Geodatabase.gdb"
  st_layers(idfg_gdb)
  rds <- sf::st_read(dsn = idfg_gdb, layer = "Road_OpenStreetMap")
  
  #'  Take a closer look  
  crs(pforest, describe = TRUE, proj = TRUE)
  crs(nlcd, describe = TRUE, proj = TRUE)
  crs(elev, describe = TRUE, proj = TRUE)
  crs(habclass, describe = TRUE, proj = TRUE)
  crs(dist2suburbs, describe = TRUE, proj = TRUE)
  st_crs(id)
  st_crs(rds)
  
  res(pforest)
  res(nlcd)
  res(elev)
  res(habclass)
  res(dist2suburbs)
  
  #'  Define projections to use when reprojecting camera locations
  wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
  (aea <- crs(nlcd, proj = TRUE))
  (nad83 <- crs(elev, proj = TRUE)) # Idaho Transverse Mercator NAD83
  (hab_crs <- crs(habclass, proj = TRUE))
  
  #'  Reproject shapefiles
  id_aea <- st_transform(id, aea)
  id_nad83 <- st_transform(id, nad83)
  
  #'  Convert NAs to 0 - happens when pixel was suburban or rural and surrounded 
  #'  by other suburban or rural pixels so no distance to calculated (should be 0)
  dist2suburbs[is.na(dist2suburbs)] <- 0
  dist2rural[is.na(dist2rural)] <- 0
  
  #'  ---------------
  ####  Camera data  ####
  #'  ---------------
  #'  Load camera location data and format
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  #'  List camera locations
  cams_list <- list(eoe_probcams_20s, eoe_probcams_21s)
  
  #'  Make camera location data spatial sf objects
  spatial_locs <- function(locs, proj) {
    locs <- arrange(locs, NewLocationID)
    sf_locs <- st_as_sf(locs, coords = c("Long", "Lat"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      dplyr::select(c("NewLocationID")) %>%
      st_transform(proj)
    return(sf_locs)
  }
  cams_wgs84 <- lapply(cams_list, spatial_locs, proj = wgs84)
  cams_aea <- lapply(cams_list, spatial_locs, proj = aea)
  cams_nad83 <- lapply(cams_list, spatial_locs, proj = nad83)
  cams_hab_crs <- lapply(cams_list, spatial_locs, proj = hab_crs)
  
  #'  Double check these are plotting correctly
  plot(pforest, main = "Camera locations over percent forested habitat, 500m radius")
  plot(id_aea[1], add = TRUE, col = NA)
  plot(cams_aea[[2]], add = TRUE, col = "black", cex = 0.75)
  
  plot(nlcd, main = "Camera locations over National Landcover Classifications")
  plot(id_aea[1], add = TRUE, col = NA)
  plot(cams_aea[[2]], add = TRUE, col = "black", cex = 0.75)
  
  
  #'  ------------------------
  ####  Additional data sets  ####
  #'  ------------------------
  #####  Predator mortality data  ####
  #'  Mortality data provided by IDFG
  load("./Data/IDFG BGMR data/mort_preSmr19.RData")
  load("./Data/IDFG BGMR data/mort_preSmr20.RData")
  load("./Data/IDFG BGMR data/mort_preWtr20.RData")
  load("./Data/IDFG BGMR data/mort_preSmr21.RData")
  
  reformat_mort_dat <- function(mort) {
    #'  Split predator mortality data by species
    mort <- dplyr::select(mort, -gmu_area)
    mort_bear <- filter(mort, Species == "Black Bear") %>%
      rename(Bear_mort_n = total_mortalities,
             Bear_mort_km2 = mortality_km2) %>%
      dplyr::select(-Species)
    mort_bob <- filter(mort, Species == "Bobcat") %>%
      rename(Bob_mort_n = total_mortalities,
             Bob_mort_km2 = mortality_km2) %>%
      dplyr::select(-c(Species, GMU))
    mort_lion <- filter(mort, Species == "Mountain Lion") %>%
      rename(Lion_mort_n = total_mortalities,
             Lion_mort_km2 = mortality_km2) %>%
      dplyr::select(-c(Species, GMU))
    mort_wolf <- filter(mort, Species == "Wolf") %>%
      rename(Wolf_mort_n = total_mortalities,
             Wolf_mort_km2 = mortality_km2) %>%
      dplyr::select(-c(Species, GMU))
    mort_df <- as.data.frame(cbind(mort_bear, mort_bob, mort_lion, mort_wolf))
    return(mort_df)
  }
  mort_Smr19_df <- reformat_mort_dat(mort_preSmr19) %>%
    filter(GMU != "GMU1")
  mort_Smr20_df <- reformat_mort_dat(mort_preSmr20) %>%
    filter(GMU != "GMU1")
  mort_Wtr20_df <- reformat_mort_dat(mort_preWtr20) %>%
    filter(GMU != "GMU1")
  mort_Smr21_df <- reformat_mort_dat(mort_preSmr21)
  
  #'  summarize mortality data
  Season <- c("2018-2019", "2018-2019", "2019-2020", "2019-2020", "2020-2021", "2020-2021", "2020-2021")
  mort_summary <- rbind(mort_Smr19_df, mort_Smr20_df, mort_Smr21_df)
  mort_summary <- cbind(Season, mort_summary)
  names(mort_summary)[names(mort_summary) == "Season"] <- "Harvest season (fall - winter)"
  
  mort_summary_noGMU1 <- filter(mort_summary, GMU != "GMU1") %>%
    dplyr::select(`Harvest season (fall - winter)`, Bear_mort_n, Bob_mort_n, Lion_mort_n, Wolf_mort_n) %>%
    group_by(`Harvest season (fall - winter)`) %>%
    summarise(across(everything(), sum),
              .groups = 'drop') %>%
    ungroup()
  
  #'  ----------------------------------
  ####  COVARIATE EXTRACTION & MERGING  ####
  #'  ----------------------------------
  cov_extract <- function(locs_aea, locs_nad83, locs_hab_crs, min_group_size, mort, relativeN) {
    
    #'  Extract covariate data for each camera site from spatial layers
    perc_forest <- terra::extract(pforest, vect(locs_aea)) %>%
      transmute(ID = ID, perc_forest = focal_sum) %>%
      mutate(perc_forest = round(perc_forest, 3))
    landcover <- terra::extract(nlcd, vect(locs_aea))
    elev <- terra::extract(elev, vect(locs_nad83))
    habitat <- terra::extract(habclass, vect(locs_hab_crs))
    dist2suburbs <- terra::extract(dist2suburbs, vect(locs_nad83))
    dist2rural <- terra::extract(dist2rural, vect(locs_nad83))
    #'  Find nearest road to each camera
    nearestrd <- st_nearest_feature(locs_nad83, rds)
    #'  Calculate distance to nearest road (in meters)
    dist2rd <- as.numeric(st_distance(locs_nad83, rds[nearestrd,], by_element = T))
    
    #'  Join each extracted covariate to the unique camera location data
    covs <- as.data.frame(locs_aea) %>%
      mutate(ID = seq(1:nrow(.))) %>%
      full_join(perc_forest, by = "ID") %>%
      full_join(landcover, by = "ID") %>%
      full_join(elev, by = "ID") %>%
      full_join(habitat, by = "ID") %>%
      full_join(dist2suburbs, by = "ID") %>%
      full_join(dist2rural, by = "ID") %>%
      cbind(dist2rd) %>%
      mutate(GMU = sub("_.*", "", NewLocationID), 
             dist2rd = round(dist2rd, digits = 2)) %>%
      relocate(GMU, .after = NewLocationID) %>%
      full_join(mort, by = "GMU") %>%
      dplyr::select(-c(geometry, ID)) %>%
      arrange(NewLocationID) %>%
      rename("Landcover_30m2" = "NLCD Land Cover Class")
      
    return(covs)
  }
  eoe_covs_20s <- cov_extract(locs_aea = cams_aea[[1]], locs_nad83 = cams_nad83[[1]], locs_hab_crs = cams_hab_crs[[1]], 
                              mort = mort_Smr20_df, relativeN = spp_diversity_Smr20)
  eoe_covs_21s <- cov_extract(locs_aea = cams_aea[[2]], locs_nad83 = cams_nad83[[2]], locs_hab_crs = cams_hab_crs[[2]], 
                              mort = mort_Smr21_df, relativeN = spp_diversity_Smr21)
  

  # ####  Save  ###
  # write.csv(eoe_covs_20s, file = "./Data/Relative abundance data/Covariates_EoE_Smr20.csv")
  # save(eoe_covs_20s, file = "./Data/Relative abundance data/Covariates_EoE_Smr20.RData")
  # 
  # write.csv(eoe_covs_21s, file = "./Data/Relative abundance data/Covariates_EoE_Smr21.csv")
  # save(eoe_covs_21s, file = "./Data/Relative abundance data/Covariates_EoE_Smr21.RData")
