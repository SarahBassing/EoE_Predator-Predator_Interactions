  #'  --------------------------------------
  #'  Extract covariates at camera locations
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  --------------------------------------
  #'  Read in spatial data and extract covariate values at camera locations
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
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  #'  Fix a couple of capitalization issues in the GMU part of NewLocationID
  eoe_probcams_21s <- eoe_probcams_21s %>%
    mutate(NewLocationID = toupper(NewLocationID))
  
  #'  Make camera location data spatial sf objects
  spatial_locs <- function(locs, proj) {
    sf_locs <- st_as_sf(locs, coords = c("Long", "Lat"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      dplyr::select(c("NewLocationID")) %>%
      st_transform(proj)
    # print(projection(sf_locs))
    return(sf_locs)
  }
  cams_aea <- spatial_locs(eoe_probcams_21s, proj = aea)
  
  #'  Double check these are plotting correctly
  plot(pforest, main = "Camera locations over percent forested habitat, 500m radius")
  plot(id_aea[1], add = TRUE, col = NA)
  plot(cams_aea, add = TRUE, col = "black", cex = 0.75)
  
  
  #'  ------------------------
  ####  Additional data sets  ####
  #'  ------------------------
  #'  Wolf data generated from camera detections
  load("./Data/Wolf count data/count_eoe21s_wolf.RData")
  load("./Data/Wolf count data/min_group_size_eoe21s.RData") 
  
  
  #'  ----------------------------------
  ####  COVARIATE EXTRACTION & MERGING  ####
  #'  ----------------------------------
  cov_extract <- function(locs) {
    #'  Extract covariate data for each camera site from spatial layers
    perc_forest <- terra::extract(pforest, vect(locs)) %>%
      transmute(ID = ID, perc_forest = focal_sum) %>%
      mutate(perc_forest = round(perc_forest,3))
    landcover <- terra::extract(nlcd, vect(locs))
    #'  Join each extracted covariate to the unique camera location data
    covs <- as.data.frame(cams_aea) %>%
      mutate(ID = seq(1:nrow(.))) %>%
      full_join(perc_forest, by = "ID") %>%
      full_join(landcover, by = "ID") %>%
      #'  Join with additional data sets already formatted for each camera site
      full_join(min_group_size_eoe21s, by = "NewLocationID") %>%
      dplyr::select(-c(geometry, ID, Lat, Long))
    return(covs)
  }
  eoe_covs_21s <- cov_extract(cams_aea)
  
  
  #'  Save
  write.csv(eoe_covs_21s, file = "./Data/Covariates_EoE_Smr21.csv")
  save(eoe_covs_21s, file = "./Data/Covariates_EoE_Smr21.RData")
  
  
  
  
  
  
  
  
  
  
  