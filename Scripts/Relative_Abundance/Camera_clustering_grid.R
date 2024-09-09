  #'  -------------------------------
  #'  Camera clustering
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  July 2023
  #'  -------------------------------
  #'  Explore how to group cameras into clusters based on clustering algorithms and grid cells
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
  
  #'  Camera locations with wolf RAI from RN models
  wolf_cams <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/spatial_rn_wolf_locs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs") %>%
    mutate(N_rounded = round(RN_n, 0)) %>%
    filter(Setup == "P") %>%
    #'  Retain the observation across years with the highest RAI
    group_by(NwLctID) %>%
    slice(which.max(RN_n)) %>%
    ungroup() #%>%
    #distinct(geometry, .keep_all = TRUE)
  
  xy <- st_coordinates(wolf_cams)
  
  crs(wolf_cams)
  
  #'  Split wolf cameras by GMU
  wolf_cams_gmu1 <- wolf_cams[wolf_cams$GMU == "GMU1",]
  wolf_cams_gmu6 <- wolf_cams[wolf_cams$GMU == "GMU6",]
  wolf_cams_gmu10a <- wolf_cams[wolf_cams$GMU == "GMU10A",]
  wolf_cams_list <- list(wolf_cams_gmu1, wolf_cams_gmu6, wolf_cams_gmu10a)
  
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
  
  #'  Visualize wolf RN abundance estimates (June - Sept 2020 - 2022)
  ggplot() +
    geom_sf(data = eoe_gmu_wgs84, fill = NA) +
    geom_sf(data = wolf_cams, aes(size = N_rounded), shape  = 21, 
            col = "darkred", fill = "darkred", alpha = 3/10) +
    scale_size_continuous(breaks = c(0, 1, 2, 3, 5, 7, 9, 12), range = c(0,12)) +
    labs(size = "Estimated \nlocal abundance", x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  
  # units::set_units(st_area(eoe_gmu), km^2) #' area of each GMU in sq-km (order = GMU1, GMU6, GMU10A)
  # units::set_units(st_area(id), km^2) 

  #'  -----------------------
  ####  Clustering approach  ####
  #'  -----------------------
  # options(repos = c(
  #   mpadge = 'https://mpadge.r-universe.dev',
  #   CRAN = 'https://cloud.r-project.org'))
  # install.packages ("spatialcluster")
  library(spatialcluster)
  library(dbscan)
  library(pracma)
  
  #'  Calculate non-spatial distance between RN estimates
  N_distance <- function(locs) {
    RN_dist <- dist(locs$RN_n)
    RN_dist <- as.matrix(RN_dist)
    return(RN_dist)
  }
  N_dist_list <- lapply(wolf_cams_list, N_distance)
  
  #'  Grab xy coordinates and calculate distance between each camera 
  xy_distance <- function(locs) {
    xy <- as.matrix(st_coordinates(locs))
    xy_dist <- dist(xy)
    xy_dist <- as.matrix(xy_dist)
    return(xy_dist)
  }
  xy_dist_list <- lapply(wolf_cams_list, xy_distance)
  
  
  kD <- pdist2(xy_dist_list[[2]], xy_dist_list[[2]])
  kDsort <- sort(kD, decreasing = FALSE)
  plot(kDsort)
  
  N_distance_dataset <- function(locs) {
    xy <- as.data.frame(st_coordinates(locs))
    RN_dist <- bind_cols(locs, xy) %>%
      dplyr::select(c(RN_n, X, Y)) %>%
      as.data.frame() %>%
      dplyr::select(-geometry)
    RN_dist <- as.matrix(RN_dist)
    return(RN_dist)
  }
  Nxy_dist_list <- lapply(wolf_cams_list, N_distance_dataset)
  
  clusters_GMU1 <- dbscan::dbscan(Nxy_dist_list[[1]], eps = 0.1, minPts = 4)
  augment(clusters_GMU1, Nxy_dist_list[[1]]) %>%
    ggplot(aes(x = X, y = Y)) + geom_point(aes(color = .cluster, shape = noise), size = 2.5) + scale_shape_manual(values = c(19, 4))
  
  clusters_GMU6 <- dbscan::dbscan(Nxy_dist_list[[2]], eps = 0.07, minPts = 4)
  augment(clusters_GMU6, Nxy_dist_list[[2]]) %>%
    ggplot(aes(x = X, y = Y)) + geom_point(aes(color = .cluster, shape = noise), size = 2.5) + scale_shape_manual(values = c(19, 4))
  
  clusters_GMU10A <- dbscan::dbscan(Nxy_dist_list[[3]], eps = 0.07, minPts = 4)
  augment(clusters_GMU10A, Nxy_dist_list[[3]]) %>%
    ggplot(aes(x = X, y = Y)) + geom_point(aes(color = .cluster, shape = noise), size = 2.5) + scale_shape_manual(values = c(19, 4))
  
  
  
  #'  Merge two distance matrices
  dbscan_cluster <- function(dist1, dist2, eps_buffer, minimumPts, locs) {
    dist <- dist1 + dist2
    clusters <- dbscan::dbscan(dist, eps = eps_buffer, minPts = minimumPts)
    print(tidy(clusters))
    augment(dist, clusters)
    locs$clusters <- clusters$cluster
    map_clusters <- ggplot() +
      geom_sf(data = locs, aes(color = clusters), size = 3) + scale_color_gradient(low = "lightblue", high = "red") + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    print(map_clusters)
    return(locs)
  }
  double_dist_list <- mapply(dbscan_cluster, dist1 = N_dist_list, dist2 = xy_dist_list, 
                             eps_buffer = 800, minimumPts = 2, locs = wolf_cams_list, SIMPLIFY = FALSE)
  
  scl_cluster <- function(locs, dist1, nclusters, linkages) {
    xy <- as.matrix(st_coordinates(locs))
    clusters <- scl_full(xy, dist1, ncl = nclusters, linkage = linkages)
    print(plot(clusters))
    return(locs)
  }
  scl_cluster_list <- mapply(scl_cluster, locs = wolf_cams_list, dist1 = N_dist_list,
                             nclusters = 15, linkages = "single", SIMPLIFY = FALSE)
  
  
  #'  ----------------------
  ####  Grid cell approach  ####
  #'  ----------------------
  #'  How many cells are needed for raster? 
  #'  Calculate length/width of EoE bbox, divide by length/width of pixel, and round to nearest whole number
  bbox <- st_bbox(eoe_gmu)
  pxl_size <- 6400 * 2 # Minimum dist. btwn rendezvous sites of adjacent packs was 6.4km (Ausband et al. 2010)
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  