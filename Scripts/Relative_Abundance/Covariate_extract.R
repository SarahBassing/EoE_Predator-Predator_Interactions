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
  #'  Spatial data
  #'  --------------------------------------
  
  #'  Load libraries
  library(sf)
  library(rgeos)
  library(terra)
  library(lubridate)
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
  mtbs <- st_read(dsn = "./Shapefiles/MTBS_perimeter_data", layer = "mtbs_perims_DD")
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
  st_crs(mtbs)
  
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
  mtbs_nad83 <- st_transform(mtbs, nad83) # Takes a hot second
  
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
  load("./Data/IDFG camera data/Problem cams/eoe22s_problem_cams.RData")
  
  #'  List camera locations
  cams_list <- list(eoe_probcams_20s, eoe_probcams_21s, eoe_probcams_22s)
  
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
  
  #' #'  Save camera locations as shapefiles
  #' cams_20s_wgs84 <- cams_wgs84[[1]]; cams_21s_wgs84 <- cams_wgs84[[2]]; cams_22s_wgs84 <- cams_wgs84[[3]]
  #' 
  #' wd <- getwd()
  #' st_write(cams_20s_wgs84, dsn = paste0(wd,"/Shapefiles/IDFG spatial data/Camera_locations/cams_20s_wgs84.shp"), layer = "cams_20s_wgs84.shp")
  #' st_write(cams_21s_wgs84, dsn = paste0(wd,"/Shapefiles/IDFG spatial data/Camera_locations/cams_21s_wgs84.shp"), layer = "cams_21s_wgs84.shp")
  #' st_write(cams_22s_wgs84, dsn = paste0(wd,"/Shapefiles/IDFG spatial data/Camera_locations/cams_22s_wgs84.shp"), layer = "cams_22s_wgs84.shp")
  
  #'  Double check these are plotting correctly
  plot(pforest, main = "Camera locations over percent forested habitat, 500m radius")
  plot(id_aea[1], add = TRUE, col = NA)
  plot(cams_aea[[2]], add = TRUE, col = "black", cex = 0.75)
  
  plot(nlcd, main = "Camera locations over National Landcover Classifications")
  plot(id_aea[1], add = TRUE, col = NA)
  plot(cams_aea[[2]], add = TRUE, col = "black", cex = 0.75)
  
  #'  ----------------------------------
  ####  COVARIATE EXTRACTION & MERGING  ####
  #'  ----------------------------------
  cov_extract <- function(locs_aea, locs_nad83, locs_hab_crs) {
    
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
    #'  Intersect overlapping cameras & burn perimeters
    mtbs_simple <- mtbs_nad83 %>% 
      dplyr::select("Ig_Date", "Incid_Name", "Incid_Type") %>%
      mutate(Burn_year = format(as.Date(Ig_Date, format="%Y-%m-%d"),"%Y")) 
    mtbs_cam_inter <- st_intersection(locs_nad83, mtbs_simple)
    #'  Save burn year
    mtbs_burn_yr <- as.data.frame(mtbs_cam_inter) %>%
      dplyr::select(c("NewLocationID", "Burn_year"))
    
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
      full_join(mtbs_burn_yr, by = "NewLocationID") %>%
      relocate(Burn_year, .after = dist2rd) %>%
      mutate(GMU = sub("_.*", "", NewLocationID), 
             dist2rd = round(dist2rd, digits = 2)) %>%
      relocate(GMU, .after = NewLocationID) %>%
      dplyr::select(-c(geometry, ID)) %>%
      arrange(NewLocationID) %>%
      rename("Landcover_30m2" = "NLCD Land Cover Class")
      
    return(covs)
  }
  covs_20s <- cov_extract(locs_aea = cams_aea[[1]], locs_nad83 = cams_nad83[[1]], locs_hab_crs = cams_hab_crs[[1]])
  covs_21s <- cov_extract(locs_aea = cams_aea[[2]], locs_nad83 = cams_nad83[[2]], locs_hab_crs = cams_hab_crs[[2]])
  covs_22s <- cov_extract(locs_aea = cams_aea[[3]], locs_nad83 = cams_nad83[[3]], locs_hab_crs = cams_hab_crs[[3]])
  
  covariate_list <- list(covs_20s, covs_21s, covs_22s)

  # ####  Save  ###
  # save(covariate_list, file = "./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022.RData")
