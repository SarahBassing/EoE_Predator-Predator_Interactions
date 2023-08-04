  #'  -------------------------------
  #'  Data formatting
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  July 2023
  #'  -------------------------------
  #'  Explore how to cluster cameras together
  #'  
  #'  Camera operations table generated in Detection_data_cleaning.R
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(exactextractr)
  library(sf)
  library(terra)
  library(ggplot2)
  library(tidyverse)
  
  #'  ----------------
  ####  Read in data  ####
  #'  ----------------
  #'  Problem cameras
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  cams_list <- list(eoe_probcams_20s, eoe_probcams_21s)
  
  #'  Spatial data
  gmu <- st_read("./Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp")
  eoe_gmu <- gmu[gmu$NAME == "1" | gmu$NAME == "6" | gmu$NAME == "10A",] 
  usa <- st_read("./Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
  id <- usa[usa$NAME == "Idaho",] 
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  id_wgs84 <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Define projections and units
  gmu_proj <- crs(eoe_gmu) #crs("+init=EPSG:2243") #NAD83 / Idaho West (ftUS)
  id_proj <- crs(id)
  wgs84 <- crs(eoe_gmu_wgs84) 
  
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
  cams_nad83 <- lapply(cams_list, spatial_locs, proj = gmu_proj)
  
  units::set_units(st_area(eoe_gmu), km^2) #' area of each GMU in sq-km (order = GMU1, GMU6, GMU10A)
  units::set_units(st_area(id), km^2) 

  #'  How many cells are needed for raster? 
  #'  Calculate length/width of EoE bbox, divide by length/width of pixel, and round to nearest whole number
  bbox <- st_bbox(eoe_gmu)
  pxl_size <- 6400 * 2 # 6400m = 6.4km min. dist. btwn rendezvous sites of adjacent packs (Ausband et al. 2010)
  xcells <- round(((bbox[3] - bbox[1]) / pxl_size), 0)
  ycells <- round(((bbox[4] - bbox[2]) / pxl_size), 0)
  
  #'  Create raster with specified pixel size
  eoe_rast <- rast(ncols = xcells, nrows = ycells, xmin = bbox[1], xmax = bbox[3], ymin = bbox[2], ymax = bbox[4], crs = gmu_proj)#"+init=EPSG:2243" 
  values(eoe_rast) <- 1:ncell(eoe_rast)
  names(eoe_rast) <- "grid_id"
  
  
  ####  Visualize  ####
  #'  Generate spatial vector data from rasters for visualization
  eoe_poly <- as.polygons(eoe_rast, values = TRUE) 
  eoe_poly_wgs84 <- terra::project(eoe_poly, wgs84)
  
  #'  Plot
  plot(eoe_poly)
  plot(eoe_gmu[1], add = TRUE)
  plot(cams_nad83[[2]], pch = 16, col = "black", cex = 0.8, add = TRUE)
  plot(eoe_poly, add = TRUE)
  
  plot(eoe_poly_wgs84)
  plot(eoe_gmu_wgs84[1], add = TRUE)
  plot(cams_wgs84[[2]], pch = 16, col = "black", cex = 0.8, add = TRUE)
  plot(eoe_poly_wgs84, add = TRUE)
  
  #'  Close up of cameras across grid cells in GMU1
  gmu1 <- gmu[gmu$NAME == "1",][1]
  gmu1_rast_poly <- crop(eoe_poly, gmu1)
  plot(gmu1_rast_poly)
  plot(gmu1, add = TRUE)
  plot(cams_nad83[[2]], pch = 16, col = "black", cex = 0.8, add = TRUE)
  plot(gmu1_rast_poly, add = TRUE)
  
  
  ####  Sampling effort per grid cell  ####
  #'  Mask raster so cell values = NA if outside GMUs 1, 6, 10A 
  eoe_mask <- mask(eoe_rast, eoe_gmu)
  
  #'  Extract area of each cell that is contained within each polygon
  overlap_area <- exact_extract(eoe_rast, eoe_gmu, coverage_area = TRUE)
  overlap_area[[1]]$GMU <- "GMU1"; overlap_area[[2]]$GMU <- "GMU6"; overlap_area[[3]]$GMU <- "GMU10A"
  area_m2 <- rbind(overlap_area[[1]], overlap_area[[2]], overlap_area[[3]])
  max_cell_size <- max(area_m2$coverage_area)
  area_m2$prop_cell <- area_m2$coverage_area/max_cell_size
  area_m2 <- relocate(area_m2, prop_cell, .after = coverage_area)
  area_m2 <- rename(area_m2, "coverage_area_m2" = "coverage_area")

  #'  Extract grid info for each camera site
  cams_grid <- list(NA)
  for(l in 1:2) {
    #'  Extract grid cell ID for each camera location per year using raster
    cam_cell_id <- terra::extract(eoe_mask, cams_nad83[[l]])
    cams_grid[[l]] <- cbind(cams_nad83[[l]], cam_cell_id)
    
    #'  Merge cell area with camera cell ID (includes cells with NO cameras)
    cams_grid[[l]] <- right_join(cams_grid[[l]], area_m2, by = c("grid_id" = "value")) %>%
      relocate(ID, .before = NewLocationID) 
  }
  head(cams_grid[[2]])
  
  #'  Summarize sampling effort per grid cell
  #'  Exclude GMU1 from year 1 (summer 2020) data
  cams_grid[[1]] <- filter(cams_grid[[1]], GMU != "GMU1")
  
  cams_per_monitored_grid <- list(NA)
  for(l in 1:2) {
    #'  Count number of cameras per cell EXCLUDING cells with no cameras
    cams_per_monitored_grid[[l]] <- cams_grid[[l]] %>%
      filter(!is.na(NewLocationID)) %>%
      group_by(grid_id) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      st_drop_geometry(.)
  }
  chk_it <- cams_per_monitored_grid[[2]]
  
  #'  Pull out grid cells with no cameras
  no_cams20s <- filter(cams_grid[[1]], is.na(NewLocationID)) %>%
    dplyr::select(grid_id) %>%
    mutate(n = 0) %>% relocate(n, .after = grid_id) %>%
    st_drop_geometry(.)
  no_cams21s <- filter(cams_grid[[2]], is.na(NewLocationID)) %>%
    dplyr::select(grid_id) %>%
    mutate(n = 0) %>% relocate(n, .after = grid_id) %>%
    st_drop_geometry(.)
  no_cams <- list(no_cams20s, no_cams21s)
  
  #'  Combine counts of cameras in monitored cells & cells without cameras
  cams_per_monitored_grid[[1]] <- rbind(cams_per_monitored_grid[[1]], no_cams20s) %>% arrange(grid_id) %>% rename(n_cams = n)
  cams_per_monitored_grid[[2]] <- rbind(cams_per_monitored_grid[[2]], no_cams21s) %>% arrange(grid_id) %>% rename(n_cams = n)
  
  #'  Average sampling effort per cell
  mean_effort <- lapply(cams_per_monitored_grid, function(x) {mean(x$n_cams, na.rm = TRUE)}) %>% unlist(.)
  sd_effort <- lapply(cams_per_monitored_grid, function(x) {sd(x$n_cams, na.rm = TRUE)}) %>% unlist(.)
  effort_per_grid <- data.frame(year = c("2020", "2021"),
                                mean_cams = round(c(mean_effort[1], mean_effort[2]), 3),
                                sd_cams = round(c(sd_effort[1], sd_effort[2]), 3))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  