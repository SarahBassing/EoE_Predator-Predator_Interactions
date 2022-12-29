  #'  --------------------------------------
  #'  Extract covariates at camera locations
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
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
  library(raster)
  library(terra)
  library(tidyverse)
  
  #'  ----------------
  ####  Spatial data  ####
  #'  ----------------
  #'  Load and review
  # tree <- rast("./Shapefiles/Tree Cover/TreeCanopy_Idaho.tif")
  pforest <- rast("./Shapefiles/National Land Cover Database (NCLD)/PercentForeest_500m.tif")
  nlcd <- rast("./Shapefiles/National Land Cover Database (NCLD)/NLCD19_Idaho.tif")
  id <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp") 
  
  crs(pforest, describe = TRUE, proj = TRUE)
  crs(nlcd, describe = TRUE, proj = TRUE)
  # projection(id)
  
  res(pforest)
  res(nlcd)
  
  #'  Define projections to use when reprojecting camera locations
  wgs84 <- crs("+proj=longlat +datum=WGS84 +no_defs")
  (aea <- crs(nlcd, proj = TRUE))
  
  #'  Reproject shapefiles
  id_aea <- st_transform(id, aea)

  
  #'  ---------------
  ####  Camera data  ####
  #'  ---------------
  #'  Load camera location data and format
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  #'  Make camera location data spatial sf objects
  spatial_locs <- function(locs, proj) {
    sf_locs <- st_as_sf(locs, coords = c("Long", "Lat"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      dplyr::select(c("NewLocationID")) %>%
      st_transform(proj)
    # print(projection(sf_locs))
    return(sf_locs)
  }
  cams_aea_eoe20s <- spatial_locs(eoe_probcams_20s, proj = aea)
  cams_aea_eoe20w <- spatial_locs(eoe_probcams_20w, proj = aea)
  cams_aea_eoe21s <- spatial_locs(eoe_probcams_21s, proj = aea)
  
  #'  Double check these are plotting correctly
  plot(pforest, main = "Camera locations over percent forested habitat, 500m radius")
  plot(id_aea[1], add = TRUE, col = NA)
  plot(cams_aea_eoe20s, add = TRUE, col = "black", cex = 0.75)
  
  
  #'  ------------------------
  ####  Additional data sets  ####
  #'  ------------------------
  #'  Wolf data generated from camera detections
  load("./Data/Wolf count data/count_eoe20s_wolf.RData")
  load("./Data/Wolf count data/min_group_size_eoe20s.RData") 
  
  load("./Data/Wolf count data/count_eoe20w_wolf.RData")
  load("./Data/Wolf count data/min_group_size_eoe20w.RData") 
  
  load("./Data/Wolf count data/count_eoe21s_wolf.RData")
  load("./Data/Wolf count data/min_group_size_eoe21s.RData")
  
  #'  Mortality data provided by IDFG
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
  mort_Smr20_df <- reformat_mort_dat(mort_preSmr20) %>%
    filter(GMU != "GMU1")
  mort_Wtr20_df <- reformat_mort_dat(mort_preWtr20) %>%
    filter(GMU != "GMU1")
  mort_Smr21_df <- reformat_mort_dat(mort_preSmr21)
  
  #'  ----------------------------------
  ####  COVARIATE EXTRACTION & MERGING  ####
  #'  ----------------------------------
  cov_extract <- function(locs, min_group_size, mort) {
    
    #'  Extract covariate data for each camera site from spatial layers
    perc_forest <- terra::extract(pforest, vect(locs)) %>%
      transmute(ID = ID, perc_forest = focal_sum) %>%
      mutate(perc_forest = round(perc_forest,3))
    landcover <- terra::extract(nlcd, vect(locs))
    
    #'  Join each extracted covariate to the unique camera location data
    covs <- as.data.frame(locs) %>%
      mutate(ID = seq(1:nrow(.))) %>%
      full_join(perc_forest, by = "ID") %>%
      full_join(landcover, by = "ID") %>%
      #'  Join with additional data sets already formatted for each camera site
      full_join(min_group_size, by = "NewLocationID") %>%
      mutate(GMU = sub("_.*", "", NewLocationID)) %>%
      relocate(GMU, .after = NewLocationID) %>%
      full_join(mort, by = "GMU") %>%
      dplyr::select(-c(geometry, ID, Lat, Long))
    
     return(covs)
  }
  eoe_covs_20s <- cov_extract(cams_aea_eoe20s, min_group_size = min_group_size_eoe20s, mort = mort_Smr20_df)
  eoe_covs_20w <- cov_extract(cams_aea_eoe20w, min_group_size = min_group_size_eoe20w, mort = mort_Wtr20_df)
  eoe_covs_21s <- cov_extract(cams_aea_eoe21s, min_group_size = min_group_size_eoe21s, mort = mort_Smr21_df)
  
  
  #'  Save
  write.csv(eoe_covs_20s, file = "./Data/Covariates_EoE_Smr20.csv")
  save(eoe_covs_20s, file = "./Data/Covariates_EoE_Smr20.RData")
  
  write.csv(eoe_covs_20w, file = "./Data/Covariates_EoE_Wtr20.csv")
  save(eoe_covs_20w, file = "./Data/Covariates_EoE_Wtr20.RData")
  
  write.csv(eoe_covs_21s, file = "./Data/Covariates_EoE_Smr21.csv")
  save(eoe_covs_21s, file = "./Data/Covariates_EoE_Smr21.RData")
  
  
  
  
  
  
  
  
  
  
  