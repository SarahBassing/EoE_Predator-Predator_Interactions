  #'  -------------------------------
  #'  Camera clustering
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  July 2023
  #'  -------------------------------
  #'  Group cameras into clusters based on clustering algorithms 
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
  library(ClustGeo)
  library(units)
  library(mapview)#; mapviewOptions(fgb = FALSE)
  library(adehabitatHR)
  library(lwgeom)
  
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
  #'  Idaho waterbodies (polygon features)
  bigwater <- st_read("./Shapefiles/National Hydrology Database Idaho State/Idaho_waterbodies_1km2.shp")
  bigwater_wgs84 <- st_transform(bigwater, crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
    filter(gnis_name == "Priest Lake" | gnis_name == "Upper Priest Lake" | 
             gnis_name == "Lake Pend Oreille" | gnis_name == "Cabinet Gorge Reservoir" |
             gnis_name == "Chatcolet Lake" | gnis_name == "Dworshak Reservoir") %>%
    dplyr::select(gnis_name)
  #' #'  Extract and generate sf line objects from big rivers
  #' #'  Idaho flowlines (line features)
  #' bigwater_line0 <- st_read("./Shapefiles/National Hydrology Database Idaho State/Shape/NHDFlowline_0.shp") %>%
  #'   st_zm(drop = TRUE, what = "ZM") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")
  #' bigwater_line1 <- st_read("./Shapefiles/National Hydrology Database Idaho State/Shape/NHDFlowline_1.shp") %>%
  #'   st_zm(drop = TRUE, what = "ZM") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")
  #' #'  Grab and combine flowlines for Upper Priest River and Priest River
  #' upr0 <- bigwater_line0 %>% filter(gnis_name == "Upper Priest River")
  #' upr1 <- bigwater_line1 %>% filter(gnis_name == "Upper Priest River")
  #' upr <- st_union(upr0, upr1) %>% st_cast(., "MULTILINESTRING")
  #' pr0 <- bigwater_line0 %>% filter(gnis_name == "Priest River")
  #' pr1 <- bigwater_line1 %>% filter(gnis_name == "Priest River")
  #' pr <- st_union(pr0, pr1) %>% st_cast(., "MULTILINESTRING")
  #' #'  Join all rivers, convert to a multilinestring, then to a linestring, and drop duplicate sub-geometries so object is a manageable size
  #' priestrivers <- st_union(upr, pr) %>% st_cast(., "MULTILINESTRING") %>% st_cast(priestrivers, "LINESTRING") %>% st_difference(.)
  # northfork0 <- bigwater_line0 %>% filter(gnis_name == "North Fork Clearwater River") %>% st_cast(., "MULTILINESTRING")
  # northfork1 <- bigwater_line1 %>% filter(gnis_name == "North Fork Clearwater River") %>% st_cast(., "MULTILINESTRING")
  # northforkclearwater <- st_union(northfork0, northfork1) %>% st_cast(., "MULTILINESTRING") %>% st_cast(., "LINESTRING") %>% st_difference(.)
  # clarkfork0 <- bigwater_line0 %>% filter(gnis_name == "Clark Fork")
  # clarkfork1 <- bigwater_line1 %>% filter(gnis_name == "Clark Fork")
  # clarkfork <- st_union(clarkfork0, clarkfork1) %>% st_cast(., "MULTILINESTRING") %>% st_cast(., "LINESTRING") %>% st_difference(.)
  priestrivers <- st_read("./Shapefiles/National Hydrology Database Idaho State/PriestRiver_flowline.shp")
  pendoreille <- st_read("./Shapefiles/National Hydrology Database Idaho State/Pendoreille_flowline.shp")
  kootenairiver <- st_read("./Shapefiles/National Hydrology Database Idaho State/KootenaiRiver_flowline.shp")
  clarkfork <- st_read("./Shapefiles/National Hydrology Database Idaho State/ClarkForkRiver_flowline.shp") 
  clearwater <- st_read("./Shapefiles/National Hydrology Database Idaho State/NorthForkClearwater_flowline.shp")
  
  #'  Extend Priest River line so it goes beyond Pend Oreille River (needed for cropping polygons)
  #'  Convert line to points and add row name to data frame
  pr_points <- st_cast(priestrivers, "MULTIPOINT") %>%
    mutate(rowID = row.names(.))
  mapview(pr_points)
  #'  Identify the lowest reach on the priest river and save as series of points 
  last_reach <-  st_cast(priestrivers[129,], "MULTIPOINT") %>% st_cast("POINT")
  mapview(list(last_reach[90,], last_reach[100,]))
  #'  Identify the last point of the last reach
  (last_point <- as.data.frame(st_coordinates(last_reach[100,])))
  #'  Create new point that extends south
  X <- -116.8972; Y <- 48.16000
  new_last_point <- bind_cols(X,Y) %>% as.data.frame(.) 
  names(new_last_point) <- c("X", "Y")
  #'  Merge with last point of last reach
  new_last_point <- rbind(last_point, new_last_point) %>%
    mutate(gns_n_2 = "Priest River")
  #'  Make df spatial
  new_last_point_sf <- st_as_sf(new_last_point, coords = c("X", "Y"), crs = "+proj=longlat +datum=WGS84 +no_defs")
  #'  Convert spatial points into a line
  new_line <- new_last_point_sf %>% group_by(gns_n_2) %>% dplyr::summarize(do_union = FALSE) %>%  
    st_cast("LINESTRING") 
  priestrivers_basic <- dplyr::select(priestrivers, "gns_n_2")
  #'  Merge with larger Priest River line
  longer_priestriver <- rbind(priestrivers_basic, new_line)
  mapview(longer_priestriver)
  
  #'  Extend Kootenai River line so it goes beyond US-Canada border by a smidge
  #'  Convert line to points and add row name to data frame
  mapview::mapview(list(eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 1,], kootenairiver))
  kr_points <- st_cast(kootenairiver, "MULTIPOINT") %>%
    mutate(rowID = row.names(.))
  mapview(kr_points)
  #'  Identify the last reach on Kootenai and save as series of points 
  last_reach <-  st_cast(kootenairiver[51,], "MULTIPOINT") %>% st_cast("POINT")
  #'  Identify the last point of the last reach
  (last_point <- as.data.frame(st_coordinates(last_reach[8,])))
  #'  Create new point that extends south
  X <- -116.5038; Y <- 49.05000
  new_last_point <- bind_cols(X,Y) %>% as.data.frame(.) 
  names(new_last_point) <- c("X", "Y")
  #'  Merge with last point of last reach
  new_last_point <- rbind(last_point, new_last_point) %>%
    mutate(gns_n_1 = "Kootenai River")
  #'  Make df spatial
  new_last_point_sf <- st_as_sf(new_last_point, coords = c("X", "Y"), crs = "+proj=longlat +datum=WGS84 +no_defs")
  #'  Convert spatial points into a line
  new_line <- new_last_point_sf %>% group_by(gns_n_1) %>% dplyr::summarize(do_union = FALSE) %>%  
    st_cast("LINESTRING") 
  kootenai_basic <- dplyr::select(kootenairiver, "gns_n_1")
  #'  Merge with larger Priest River line
  longer_kootenairiver <- rbind(kootenai_basic, new_line)
  mapview(longer_kootenairiver)
  
  #'  Split GMU1 by Upper Priest River, Priest River, and Kootenai River
  gmu1_priestriver_split <- eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 1,] %>% lwgeom::st_split(longer_priestriver) %>%
    st_collection_extract("POLYGON")
  #'  Save just the E & W sides of GMU1 as unique polygons
  gmu1_E <- gmu1_priestriver_split[1,]; plot(gmu1_E[1])
  gmu1_W <- gmu1_priestriver_split[2,]; plot(gmu1_W[1])
  mapview::mapview(list(gmu1_E, gmu1_W))
  gmu1_kootenairiver_split <- gmu1_E %>% lwgeom::st_split(longer_kootenairiver) %>%  
    st_collection_extract("POLYGON")
  gmu1_NE <- gmu1_kootenairiver_split[2,]; plot(gmu1_NE[1])
  gmu1_C <- gmu1_kootenairiver_split[1,]; plot(gmu1_C[1])
  mapview::mapview(list(gmu1_W, gmu1_NE, gmu1_C))
  
  #'  Grab just Pend Oreille Lake and convert to line feature
  lakependoreille <- bigwater_wgs84 %>% filter(gnis_name == "Lake Pend Oreille") %>% st_cast(., "MULTILINESTRING")
  priestlakes <- bigwater_wgs84 %>% filter(gnis_name == "Priest Lake" | gnis_name == "Upper Priest Lake") %>% st_cast(., "MULTILINESTRING")
  
  #'  Split GMU10A by Dwarshak Reservior
  #'  Grab just Dwarshak polygon and convert to a line feature
  dworshak <- bigwater_wgs84 %>% filter(gnis_name == "Dworshak Reservoir") %>% st_cast(., "MULTILINESTRING")
  #'  Split GMU10A using st_split and then collect all polygons into single sf object
  gmu10a_dworshak_split <- eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "10A",] %>% lwgeom::st_split(dworshak) %>%
    st_collection_extract("POLYGON")
  #'  Save just the N & S ends of GMU10A as unique polygons
  gmu10a_N <- gmu10a_dworshak_split[1,]; plot(gmu10a_N[1])
  gmu10a_S <- gmu10a_dworshak_split[53,]; plot(gmu10a_S[1])
  mapview::mapview(list(gmu10a_N, gmu10a_S))
  
  #'  Camera locations with wolf RAI from RN models
  wolf_cams <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/spatial_rn_wolf_locs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs") #%>%
    # filter(Setup == "P") %>%
    #' #'  Retain the observation across years with the highest RAI
    #' group_by(NwLctID) %>%
    #' slice(which.max(RN_n)) %>%
    #' ungroup() #%>%
    #distinct(geometry, .keep_all = TRUE)
  
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
    geom_sf(data = bigwater_wgs84, fill = "lightblue") +
    geom_sf(data = eoe_gmu_wgs84, fill = NA) +
    geom_sf(data = wolf_cams, aes(size = RN_n_rn), shape  = 21, 
            col = "darkred", fill = "darkred", alpha = 3/10) +
    scale_size_continuous(breaks = c(0, 1, 2, 3, 5, 7, 9, 12), range = c(0,12)) +
    labs(size = "Estimated \nlocal abundance", x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  
  #'  --------------------------------------------------------
  ####  Cluster camera sites by geographic & RAI homogeneity  ####
  #'  --------------------------------------------------------
  #'  Select value for defining nearest neighbor using a k-distance graph
  #'  https://www.mathworks.com/help/stats/dbscan-clustering.html#mw_0ad71796-4ae1-4620-a0e5-9ea5dcded2c6
  kD_plot <- function(locs, n) {
    #'  Remove duplicated cameras so each location represented once (avoids lots of 0 dist)
    skinny_locs <- distinct(locs, geometry, .keep_all = TRUE)
    #'  Calculate distance (m) between all locations
    dist.mat <- st_distance(skinny_locs)
    #'  Sort and return the distances between the first n nearest neighbors... 
    #'  excluding the 1 nearest distance which is 0 b/c the closes point is always 
    #'  going to be itself in st_distance
    nn.dist <- apply(dist.mat, 1, function(x) {return(sort(x, partial = 2)[2:n+1])})
    #'  Find average distance between n nearest neighbors
    avg.nn.dist <- apply(nn.dist, 2, function(x) {return(mean(x))})
    print(mean(avg.nn.dist))
    #'  Sort n nearest neighbor distances 
    kDsort <- sort(nn.dist, decreasing = FALSE)
    #'  Plot k-distance graph
    plot(kDsort)
    return(kDsort)
  }
  kDsort_byGMU <- lapply(wolf_cams_list, kD_plot, n = 1)
  
  #'  Use Ward-like hierarchical clustering to group cameras based on geographical distance & RAI
  ward_like_cluster <- function(dat, k, a, nndist, a2) {
    #'  Create dissimilarity matrices from RAI and geographic data
    #'  Grab just RAI estimates (RN_n vs RN_n_rn)
    rai <- dat %>% dplyr::select(RN_n) %>% as.data.frame() %>% dplyr::select(-geometry)
    #'  Calculate distance for RAI across sites
    D0 <- dist(rai)
    #'  Calculate spatial distances between camera sites
    xy_dist <- st_distance(dat)
    #'  Remove units and turn into a dist object
    D1 <- as.dist(drop_units(xy_dist))
    
    #'  Generate clusters by ignoring geographical distances (for demonstration only)
    tree_D0 <- hclustgeo(D0)
    partition_D0tree <- cutree(tree_D0, k)
    cams <- as(dat, "Spatial")
    sp::plot(cams, col = partition_D0tree, pch = 19, main = paste("Clusters ignoring geographical distances, k =", k)) 
    legend("topleft", legend = paste("cluster", 1:k), fill = 1:k, bty = "n")
    
    #'  Generate clusters by ignoring RAI distances (for demonstration only)
    tree_D1 <- hclustgeo(D1)
    partition_D1tree <- cutree(tree_D1, k)
    cams <- as(dat, "Spatial")
    sp::plot(cams, col = partition_D1tree, pch = 19, main = paste("Clusters ignoring RAI distances, k =", k)) 
    legend("topleft", legend = paste("cluster", 1:k), fill = 1:k, bty = "n")
    
    #'  Look for alpha that retains as much homogeneity as possible from D0 and D1 perspective
    #'  alpha = 0 all weight on D0, alpha = 1 all weight on D1
    #'  Must provide a pre-defined number of clusters
    cr <- choicealpha(D0, D1, range.alpha = seq(0, 1, 0.1), K = k, graph = TRUE)
    #'  Identify clusters with both dissimilarity matrices, weighting importance
    #'  of the constraints to increase geographical homogeneity without losing too 
    #'  much homogeneity of RAI within each cluster
    tree_D0D1 <- hclustgeo(D0, D1, alpha = a)
    #'  Split data into clusters
    Partitions <- cutree(tree_D0D1, k)
    #'  Append partitions to original data
    dat$Clusters <- Partitions
    #'  Visualize
    cams <- as(dat, "Spatial")
    sp::plot(cams, col = Partitions, pch = 19, 
             main = paste("Clusters weighted by RAI & geographical distances \nk =", k, "a =", a))
    legend("topleft", legend = paste("cluster", 1:k), fill = 1:k, bty = "n")
    
    #'  Adjust clusters by accounting for spatial adjacency among cameras (incorporate spatial autocorrelation)
    #'  Find all nearest neighbors within set distance (spherical distance in km)
    nearst_cams <- spdep::dnearneigh(sp::coordinates(cams), 0, nndist, longlat = TRUE) 
    #'  Grab camera names
    NwLctID_label <- as.vector(cams$NwLctID)
    #'  Nearest neighbors for an example camera
    print(NwLctID_label[nearst_cams[[10]]]) 
    #'  Build adjacency matrix with binary coding (flag all neighboring cameras for a given cam)
    adjacency_matrix <- spdep::nb2mat(nearst_cams, style = "B") 
    #'  Fill in 1 on the diagonal (each camera is its own neighbor)
    diag(adjacency_matrix) <- 1
    #'  Label matrix rows & columns with corresponding camera name
    colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- NwLctID_label
    #'  Replace D1 distance matrix with adjacency matrix, but inverted (switch 0s and 1)
    D1 <- 1-adjacency_matrix 
    #'  Trun it into a dist object
    D1 <- as.dist(D1)
    #'  Look for alpha that retains homogeneity for D0 & D1
    #'  Focus on alpha for standardized values of Q if D0 and D1 are on extremely different scales
    cr <- choicealpha(D0, D1, range.alpha = seq(0, 1, 0.1), K = k, graph = TRUE)
    cr$Q
    cr$Qnorm
    #'  Identify clusters again, this time with adjacency matrix and new alpha
    tree_D0D1a <- hclustgeo(D0, D1, alpha = a2)
    #'  And visualize
    plot(tree_D0D1a, hang = -1, label = FALSE, xlab = "", sub = "", main = "")
    Partitions_adjacency <- cutree(tree_D0D1a, k)
    #'  Append partitions to original data
    dat$Clusters_adjacency <- Partitions_adjacency
    sp::plot(cams, col = Partitions_adjacency, pch = 19, 
             main = paste("Clusters weighted by RAI & adjacency matrix \nk =", k, "a =", a2))
    legend("topleft", legend = paste("cluster", 1:k), fill = 1:k, bty = "n")
    
    return(dat)
  }
  #'  NOTES about parameter values:
  #'  k based on 2015 total documented and suspected packs per GMU, at peak of
  #'  wolf population growth & likely all usable N Idaho territories occupied 
  #'  Run to identify a & a2 based on Q0Q1-alpha plots, then re-run w/ updated values
  #'  Current nndist = 10km based on results from k-distance plots (used to define 
  #'  nearest neighbors and plots average distances between each point - values
  #'  past the "elbow" or "knee" are usually noise based on DBSCAN clustering method)
  clusters_gmu1 <- ward_like_cluster(wolf_cams_gmu1, k = 9, a = 0.6, nndist = 10, a2 = 0.3) # current best: all cams, k = 9, a = 0.6, nndist = 10, a2 = 0.3 weighted to favor geographic distance a bit more, RAI and geographic distance that does not account for adjacency
  mapview::mapview(clusters_gmu1, zcol = "Clusters")
  clusters_gmu6 <- ward_like_cluster(wolf_cams_gmu6, k = 4, a = 0.3, nndist = 10, a2 = 0.1) # current best: all cams, k = 4, a = 0.3, nndist = 10, a2 = 0.1 weighted to favor geographic distance a bit more, RAI and geographic distance that does not account for adjacency
  mapview::mapview(clusters_gmu6, zcol = "Clusters")
  clusters_gmu10a <- ward_like_cluster(wolf_cams_gmu10a, k = 7, a = 0.5, nndist = 10, a2 = 0.1) # current best: all cams, k = 7, a = 0.4, nndist = 10, a2 = 0.1 weighted to favor geographic distance a bit more, RAI and geographic distance that does not account for adjacency
  mapview::mapview(clusters_gmu10a, zcol = "Clusters")
  
  #'  Split camera data based on large waterbodies that act as barriers
  wolf_cams_gmu1_E <- st_intersection(wolf_cams_gmu1, gmu1_E) %>% dplyr::select(c("NwLctID", "CellID", "GMU", "Setup", "Species", "RN_n", "RN_sd", "RN_n_rn"))
  wolf_cams_gmu1_NE <- st_intersection(wolf_cams_gmu1, gmu1_NE) %>% dplyr::select(c("NwLctID", "CellID", "GMU", "Setup", "Species", "RN_n", "RN_sd", "RN_n_rn"))
  wolf_cams_gmu1_C <- st_intersection(wolf_cams_gmu1, gmu1_C) %>% dplyr::select(c("NwLctID", "CellID", "GMU", "Setup", "Species", "RN_n", "RN_sd", "RN_n_rn"))
  wolf_cams_gmu1_W <- st_intersection(wolf_cams_gmu1, gmu1_W) %>% dplyr::select(c("NwLctID", "CellID", "GMU", "Setup", "Species", "RN_n", "RN_sd", "RN_n_rn"))
  wolf_cams_gmu10a_N <- st_intersection(wolf_cams_gmu10a, gmu10a_N) %>% dplyr::select(c("NwLctID", "CellID", "GMU", "Setup", "Species", "RN_n", "RN_sd", "RN_n_rn"))
  wolf_cams_gmu10a_S <- st_intersection(wolf_cams_gmu10a, gmu10a_S) %>% dplyr::select(c("NwLctID", "CellID", "GMU", "Setup", "Species", "RN_n", "RN_sd", "RN_n_rn")) %>%
    #'  Add "GMU10A_P_59" back in (falls slightly outside of GMU10a boundaries)
    bind_rows(wolf_cams_gmu10a[wolf_cams_gmu10a$NwLctID == "GMU10A_P_59",])
  
  #'  Rerun clustering for cameras split into GMU1 E & W
  #'  Priest River appears to split 2 of the 9 clusters identified in clusters_gmu1, 
  #'  and 1 cluster is entirely on W side so allowing for 3 on W side, and up to 8 on E side
  clusters_gmu1E <- ward_like_cluster(wolf_cams_gmu1_E, k = 8, a = 0.8, nndist = 10, a2 = 0.4) # current best: E cams, k = 8, a = 0.8, nndist = 10, a2 = 0.4
  mapview::mapview(clusters_gmu1E, zcol = "Clusters")
  clusters_gmu1NE <- ward_like_cluster(wolf_cams_gmu1_NE, k = 2, a = 0.9, nndist = 10, a2 = 0.5) # current best: NE cams, k = 2, a = 0.9, nndist = 10, a2 = 0.5
  mapview::mapview(clusters_gmu1NE, zcol = "Clusters")
  clusters_gmu1C <- ward_like_cluster(wolf_cams_gmu1_C, k = 7, a = 0.8, nndist = 10, a2 = 0.6) # current best: C cams, k = 7, a = 0.8, nndist = 10, a2 = 0.6
  mapview::mapview(clusters_gmu1C, zcol = "Clusters")
  clusters_gmu1W <- ward_like_cluster(wolf_cams_gmu1_W, k = 3, a = 0.7, nndist = 10, a2 = 0.3) # current best: W cams, k = 3, a = 0.7, nndist = 10, a2 = 0.3
  mapview::mapview(clusters_gmu1W, zcol = "Clusters")
  mapview::mapview(list(clusters_gmu1W, clusters_gmu1C, clusters_gmu1NE), zcol = "Clusters")
  
  
  #'  Rerun clustering for cameras split into GMU10A N & S
  #'  Dworshak Reservior appears to split 2 of the 7 clusters identified in clusters_gmu10a,
  #'  1 cluster is entirely on N side, 2 clusters split by Dworshak are overlapping on
  #'  N side so allowing for 2 on N side, and up to 7 on S side
  clusters_gmu10aN <- ward_like_cluster(wolf_cams_gmu10a_N, k = 2, a = 0.6, nndist = 10, a2 = 0.4) # current best: N cams, k = 2, a = 0.6, nndist = 10, a2 = 0.4
  mapview::mapview(clusters_gmu10aN, zcol = "Clusters")
  clusters_gmu10aS <- ward_like_cluster(wolf_cams_gmu10a_S, k = 7, a = 0.5, nndist = 10, a2 = 0.2) # current best: S cams, k = 7, a = 0.5, nndist = 10, a2 = 0.2
  mapview::mapview(clusters_gmu10aS, zcol = "Clusters")
  
  #'  ---------------------------------------
  ####  Create polygons around each cluster  ####
  #'  ---------------------------------------
  #'  Function to create convex hull for cameras in each cluster
  starter_hull <- function(dat, buff) {
    df.sf <- dat %>% st_transform(gmu_proj) %>%
      dplyr::select(Clusters)
    #'  Visualize with Mapview package
    print(mapview::mapview(df.sf))
    #'  Group and summarise by cluster, and draw hulls
    hulls <- df.sf %>%
      group_by(Clusters) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_convex_hull()
    #'  Create Minimum Convex Polygons
    dat.sp <- as(df.sf, "Spatial")
    mcps <- mcp(dat.sp, percent = 99)
    mcps.sf <- st_as_sf(mcps) %>%
      rename("Clusters" = "id") 
    #'  Expand boundary of MCPs by average nearest neighbor distance
    buff_mcps <- st_buffer(mcps.sf, dist = buff) %>%
      st_transform(wgs84)
    #'  View result
    # print(mapview::mapview(list(df.sf, mcps.sf), zcol = "Clusters"))
    print(mapview::mapview(list(df.sf, buff_mcps), zcol = "Clusters"))
    return(buff_mcps)
  }
  starter_hulls_gmu1W <- starter_hull(clusters_gmu1W, buff = 3336.067) # buff based on avg. nearest neighbor distance from k-distance plot
  starter_hulls_gmu1E <- starter_hull(clusters_gmu1E, buff = 3336.067)
  starter_hulls_gmu1NE <- starter_hull(clusters_gmu1NE, buff = 3336.067)
  starter_hulls_gmu1C <- starter_hull(clusters_gmu1C, buff = 3336.067)
  starter_hulls_gmu6 <- starter_hull(clusters_gmu6, buff = 1157.688) 
  starter_hulls_gmu10aN <- starter_hull(clusters_gmu10aN, buff = 1237.789) 
  starter_hulls_gmu10aS <- starter_hull(clusters_gmu10aS, buff = 1237.789) 
  
  #'  Function to create density polygons (utility distributions) for cameras in each cluster
  starter_UDs <- function(dat) {
    #'  Grab just the cluster column (and geometries)
    clusters <- dat %>% dplyr::select(Clusters)
    #'  Convert sf object to a sp SpatialPointsDataFrame
    clusters <- as(clusters, "Spatial")
    #'  Create utilization distribution for each cluster of cameras
    #'  h controls the smoothing parameter
    #'  kern defines how to calculate the UD -- using epa to keep the UD tighter (bivnorm creates huge polygons for some of these)
    kud <- kernelUD(clusters, h = "href", kern = "epa")  
    # kud <- kernelUD(clusters, h = "href", kern = "bivnorm")
    #'  Grab 95% - 85% vertices for each cluster UD
    cluster95 <- getverticeshr(kud, 95)
    cluster90 <- getverticeshr(kud, 90)
    cluster85 <- getverticeshr(kud, 85)
    #'  Convert to sf object
    cluster95_sf <- st_as_sf(cluster95) %>% rename("Clusters" = "id")
    cluster90_sf <- st_as_sf(cluster90) %>% rename("Clusters" = "id")
    cluster85_sf <- st_as_sf(cluster85) %>% rename("Clusters" = "id")
    #'  Visualize clusters with mapview
    print(mapview::mapview(list(dat, cluster95_sf), zcol = "Clusters")) # 95UDs contain all points
    # print(mapview::mapview(list(dat, cluster90_sf), zcol = "Clusters")) # 90UDs exclude some cams... don't want that
    # print(mapview::mapview(list(dat, cluster85_sf), zcol = "Clusters")) # excludes some cams... no good
    return(cluster95_sf)
  }
  UDs_gmu1W <- starter_UDs(clusters_gmu1W) 
  UDs_gmu1E <- starter_UDs(clusters_gmu1E) 
  UDs_gmu1NE <- starter_UDs(clusters_gmu1NE) 
  UDs_gmu1C <- starter_UDs(clusters_gmu1C) 
  UDs_gmu6 <- starter_UDs(clusters_gmu6)
  UDs_gmu10aN <- starter_UDs(clusters_gmu10aN) 
  UDs_gmu10aS <- starter_UDs(clusters_gmu10aS) 

  #'  ----------------------------
  #####  Adjust cluster polygons  #####
  #'  ----------------------------
  #'  Can't have overlapping polygons or polygons with shared camera sites
  #'  Create unique & non-overlapping cluster polygons through unions
  #'  Rules for joining intersected polygons:
  #'    1) If there is a mix of cams from different clusters in intersection,  
  #'    assign intersection to polygon w/ larger proportion of cams in the 
  #'    intersection 
  #'    2) If intersection is empty of cameras, assign to cluster polygon with 
  #'    cameras closest to the edge of the intersection
  #"  ----------------------------
  #'  Turn off spherical geometry (st_union doesn't work with globe type geometry like an S2 object)
  sf::sf_use_s2(FALSE)
  
  ######  GMU1 Cluster Polygons  ######
  #'  ---------------------------
  #'  Drop UDs where all points are also contained within other UDs
  UDs_gmu1E <- UDs_gmu1E %>% filter(Clusters != 5)
  UDs_gmu1C <- UDs_gmu1C %>% filter(Clusters != 5)
  mapview::mapview(list(clusters_gmu1E, UDs_gmu1E), zcol = "Clusters")
  mapview::mapview(list(clusters_gmu1NE, UDs_gmu1NE), zcol = "Clusters")
  mapview::mapview(list(clusters_gmu1C, UDs_gmu1C), zcol = "Clusters")
  mapview::mapview(list(clusters_gmu1W, UDs_gmu1W), zcol = "Clusters")  
  
  #'  Intersect overlapping UDs
  UDs_gmu1E_intersect <- st_intersection(UDs_gmu1E) %>%
    mutate(NewClusters = row.names(.)) %>%
    dplyr::select(c(Clusters, NewClusters)) 
  mapview::mapview(list(clusters_gmu1E, UDs_gmu1E_intersect), zcol = "Clusters")
  UDs_gmu1NE_intersect <- st_intersection(UDs_gmu1NE) %>%
    mutate(NewClusters = row.names(.)) %>%
    dplyr::select(c(Clusters, NewClusters)) 
  mapview::mapview(list(clusters_gmu1NE, UDs_gmu1NE_intersect), zcol = "Clusters")
  UDs_gmu1C_intersect <- st_intersection(UDs_gmu1C) %>%
    mutate(NewClusters = row.names(.)) %>%
    dplyr::select(c(Clusters, NewClusters)) 
  mapview::mapview(list(clusters_gmu1C, UDs_gmu1C_intersect), zcol = "Clusters")
  UDs_gmu1W_intersect <- st_intersection(UDs_gmu1W) %>%
    mutate(NewClusters = row.names(.)) %>%
    dplyr::select(c(Clusters, NewClusters)) 
  mapview::mapview(list(clusters_gmu1W, UDs_gmu1W_intersect), zcol = "Clusters")
  
  #'  Snag individual cluster polygons that are stand alone (different ones depending on how GMU1 is split up by rivers)
  ud_gmu1E_c3 <- filter(UDs_gmu1E, Clusters == 3) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON")
  ud_gmu1E_c6 <- filter(UDs_gmu1E, Clusters == 6) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON")
  ud_gmu1NE_c2 <- filter(UDs_gmu1NE, Clusters == 2) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON")
  ud_gmu1C_c2 <- filter(UDs_gmu1C, Clusters == 2) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON")
  ud_gmu1C_c3 <- filter(UDs_gmu1C, Clusters == 3) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON")
  ud_gmu1W_c1 <- filter(UDs_gmu1W, Clusters == 1) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON")
  ud_gmu1W_c3 <- filter(UDs_gmu1W, Clusters == 3) %>% st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON")
  
  #'  Snag portions of remaining polygons that were created with intersection
  ud_gmu1E_c1 <- filter(UDs_gmu1E_intersect, NewClusters == 1)
  ud_gmu1E_c1.2 <- filter(UDs_gmu1E_intersect, NewClusters == 1.2)
  ud_gmu1E_c4 <- filter(UDs_gmu1E_intersect, NewClusters == 4)
  ud_gmu1E_c2 <- filter(UDs_gmu1E_intersect, NewClusters == 2)
  ud_gmu1E_c2.1 <- filter(UDs_gmu1E_intersect, NewClusters == 2.1)
  ud_gmu1E_c7 <- filter(UDs_gmu1E_intersect, NewClusters == 7)
  ud_gmu1NE_c1 <- filter(UDs_gmu1NE_intersect, NewClusters == 1) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu1NE_c1)
  ud_gmu1C_c4 <- filter(UDs_gmu1C_intersect, NewClusters == 4) %>% st_cast(., "GEOMETRYCOLLECTION") %>% 
    st_collection_extract("POLYGON"); mapview(ud_gmu1C_c4)
  ud_gmu1C_c1 <- filter(UDs_gmu1C_intersect, NewClusters == 1 | NewClusters == 1.2) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu1C_c1)
  ud_gmu1C_c6 <- filter(UDs_gmu1C_intersect, NewClusters == 6 | NewClusters == 6.1) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu1C_c6[c(1,4),])
  ud_gmu1C_c7 <- filter(UDs_gmu1C_intersect, NewClusters == 7 | NewClusters == 6.1) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu1C_c7[c(1,5),])
  ud_gmu1W_c2 <- filter(UDs_gmu1W_intersect, NewClusters == 2) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu1W_c2[1,])
  ud_gmu1W_c2 <- ud_gmu1W_c2[1,]
  
  #'  Split GMU1 E amalgamation of Cluster 8 polygons into individual polygons
  ud_gmu1E_c8 <- filter(UDs_gmu1E_intersect, NewClusters == 8 | NewClusters == 7.1) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu1E_c8[c(1,4),])
  #'  Save individual cluster 8 polygons so they can be joined as needed
  ud_gmu1E_c8a <- ud_gmu1E_c8[4,] %>% st_cast("POLYGON"); ud_gmu1E_c8b <- ud_gmu1E_c8[1,] %>% st_cast("POLYGON"); mapview(list(ud_gmu1E_c8a, ud_gmu1E_c8b))
  ud_gmu1E_c8c <- ud_gmu1E_c8[2,]; ud_gmu1E_c8d <- ud_gmu1E_c8[3,]; mapview(list(ud_gmu1E_c8c, ud_gmu1E_c8d))
  
  #'  Create new polygons by joining specific intersections
  ud_gmu1E_c2 <- st_union(ud_gmu1E_c2, ud_gmu1E_c2.1)
  ud_gmu1E_c4 <- st_union(ud_gmu1E_c4, ud_gmu1E_c1.2)
  ud_gmu1E_c8 <- st_union(ud_gmu1E_c8a, ud_gmu1E_c8b) %>% st_cast(., "GEOMETRYCOLLECTION") %>% 
    st_collection_extract("POLYGON") #'  Not exactly sure why this is needed but throws errors without these steps
  ud_gmu1E_c7 <- st_union(ud_gmu1E_c7, ud_gmu1E_c8c) %>% st_union(., ud_gmu1E_c8d)
  ud_gmu1C_c1 <- st_union(ud_gmu1C_c1[1,], ud_gmu1C_c1[2,]) %>% dplyr::select(Clusters); mapview(ud_gmu1C_c1)
  ud_gmu1C_c6 <- st_union(ud_gmu1C_c6[1,], ud_gmu1C_c6[4,]) %>% dplyr::select(Clusters); mapview(ud_gmu1C_c6)
  ud_gmu1C_c7 <- st_union(ud_gmu1C_c7[5,], ud_gmu1C_c7[1,]) %>% st_cast(., "GEOMETRYCOLLECTION") %>% 
    st_collection_extract("POLYGON") %>% dplyr::select(Clusters); mapview(ud_gmu1C_c7)
 
  #'  Remove large waterbodies from to GMU1 UDs
  ud_gmu1E_c1_pr_split <- ud_gmu1E_c1 %>% lwgeom::st_split(pendoreille) %>%  
    st_collection_extract("POLYGON") %>% lwgeom::st_split(longer_priestriver) %>%  
    st_collection_extract("POLYGON")
  ud_gmu1E_c2_pr_split <- ud_gmu1E_c2 %>% lwgeom::st_split(longer_priestriver) %>%
    st_collection_extract("POLYGON")
  ud_gmu1E_c4_pr_split <- ud_gmu1E_c4 %>% lwgeom::st_split(longer_priestriver) %>%  
    st_collection_extract("POLYGON")
  ud_gmu1E_c7_pr_split <- ud_gmu1E_c7 %>% lwgeom::st_split(pendoreille) %>%
    st_collection_extract("POLYGON")
  
  ud_gmu1NE_c1_kr_split <- ud_gmu1NE_c1 %>% lwgeom::st_split(kootenairiver) %>%  
    st_collection_extract("POLYGON") 
  ud_gmu1NE_c2_kr_split <- ud_gmu1NE_c2 %>% lwgeom::st_split(longer_kootenairiver) %>%
    st_collection_extract("POLYGON")
  
  ud_gmu1C_c1_pr_split <- ud_gmu1C_c1 %>% lwgeom::st_split(pendoreille) %>%  
    st_collection_extract("POLYGON") %>% lwgeom::st_split(longer_priestriver) %>%  
    st_collection_extract("POLYGON") %>% lwgeom::st_split(priestlakes) %>%  
    st_collection_extract("POLYGON")
  ud_gmu1C_c2_kr_split <- ud_gmu1C_c2 %>% lwgeom::st_split(kootenairiver) %>%
    st_collection_extract("POLYGON") %>% lwgeom::st_split(longer_priestriver) %>%  
    st_collection_extract("POLYGON") %>% lwgeom::st_split(priestlakes) %>%  
    st_collection_extract("POLYGON")
  ud_gmu1C_c4_kr_split <- ud_gmu1C_c4 %>% lwgeom::st_split(kootenairiver) %>%
    st_collection_extract("POLYGON") %>% lwgeom::st_split(longer_priestriver) %>%  
    st_collection_extract("POLYGON") %>% lwgeom::st_split(priestlakes) %>%  
    st_collection_extract("POLYGON")
  ud_gmu1C_c6_por_split <- ud_gmu1C_c6 %>% lwgeom::st_split(lakependoreille) %>%
    st_collection_extract("POLYGON") %>% lwgeom::st_split(clarkfork) %>%  
    st_collection_extract("POLYGON")
  ud_gmu1C_c7_kr_split <- ud_gmu1C_c7 %>% lwgeom::st_split(kootenairiver) %>%
    st_collection_extract("POLYGON")
  
  ud_gmu1W_c1_pr_split <- ud_gmu1W_c1 %>% lwgeom::st_split(pendoreille) %>%  
    st_collection_extract("POLYGON") %>% lwgeom::st_split(longer_priestriver) %>%  
    st_collection_extract("POLYGON")
  ud_gmu1W_c2_pr_split <- ud_gmu1W_c2 %>% lwgeom::st_split(longer_priestriver) %>%
    st_collection_extract("POLYGON") %>% lwgeom::st_split(priestlakes) %>%  
    st_collection_extract("POLYGON") 
  ud_gmu1W_c3_pr_split <- ud_gmu1W_c3 %>% lwgeom::st_split(longer_priestriver) %>%  
    st_collection_extract("POLYGON") %>% lwgeom::st_split(priestlakes) %>%  
    st_collection_extract("POLYGON")
  
  #'  Save just the E & W sides of GMU1 as unique polygons
  ud_gmu1E_c1_noriver <- ud_gmu1E_c1_pr_split[2,]; plot(ud_gmu1E_c1_noriver[2])
  ud_gmu1E_c2_noriver <- ud_gmu1E_c2_pr_split[2,]; plot(ud_gmu1E_c2_noriver[1])
  ud_gmu1E_c4_noriver <- ud_gmu1E_c4_pr_split[1,]; plot(ud_gmu1E_c4_noriver[1])
  ud_gmu1E_c7_noriver <- ud_gmu1E_c7_pr_split[1,]; plot(ud_gmu1E_c7_noriver[1]) 
  
  ud_gmu1NE_c1_noriver <- ud_gmu1NE_c1_kr_split[1,]; plot(ud_gmu1NE_c1_noriver[2])
  row.names(ud_gmu1NE_c1_noriver) <- 1
  ud_gmu1NE_c2_noriver <- ud_gmu1NE_c2_kr_split[2,]; plot(ud_gmu1NE_c2_noriver[1])
  row.names(ud_gmu1NE_c2_noriver) <- 2
  
  ud_gmu1C_c1_noriver <- ud_gmu1C_c1_pr_split[2,]; plot(ud_gmu1C_c1_noriver[1])  
  row.names(ud_gmu1C_c1_noriver) <- 1
  ud_gmu1C_c2_noriver <- ud_gmu1C_c2_kr_split[3,]; plot(ud_gmu1C_c2_noriver[1])  
  row.names(ud_gmu1C_c2_noriver) <- 2
  ud_gmu1C_c4_noriver <- ud_gmu1C_c4_kr_split[7,]; plot(ud_gmu1C_c4_noriver[1]) 
  row.names(ud_gmu1C_c4_noriver) <- 4
  ud_gmu1C_c6_noriver <- ud_gmu1C_c6_por_split[1,]; plot(ud_gmu1C_c6_noriver[1])
  row.names(ud_gmu1C_c6_noriver) <- 6
  ud_gmu1C_c7_noriver <- ud_gmu1C_c7_kr_split[1,]; plot(ud_gmu1C_c7_noriver[1])
  row.names(ud_gmu1C_c7_noriver) <- 7
  
  ud_gmu1W_c1_noriver <- ud_gmu1W_c1_pr_split[1,]; plot(ud_gmu1W_c1_noriver[1])
  ud_gmu1W_c2_noriver <- ud_gmu1W_c2_pr_split[1,]; plot(ud_gmu1W_c2_noriver[1])
  ud_gmu1W_c3_noriver <- ud_gmu1W_c3_pr_split[1,]; plot(ud_gmu1W_c3_noriver[1])
  
  #'  Spatial data frame of new cluster polygons
  UDs_gmu1E_poly <- bind_rows(ud_gmu1E_c1_noriver, ud_gmu1E_c2_noriver, ud_gmu1E_c3, 
                              ud_gmu1E_c4_noriver, ud_gmu1E_c6, ud_gmu1E_c7, ud_gmu1E_c8) %>%
    dplyr::select(-c(NewClusters, Clusters.1, NewClusters.1, Clusters.2, NewClusters.2, area))
  mapview::mapview(list(clusters_gmu1E, UDs_gmu1E_poly), zcol = "Clusters")
  UDs_gmu1NE_poly <- bind_rows(ud_gmu1NE_c1_noriver, ud_gmu1NE_c2_noriver) %>% 
    dplyr::select(Clusters)
  mapview::mapview(list(clusters_gmu1NE, UDs_gmu1NE_poly), zcol = "Clusters")
  UDs_gmu1C_poly <- bind_rows(ud_gmu1C_c1_noriver, ud_gmu1C_c2_noriver, ud_gmu1C_c3, 
                              ud_gmu1C_c4_noriver, ud_gmu1C_c6_noriver, ud_gmu1C_c7_noriver) %>%
    dplyr::select(-c(NewClusters, area))
  mapview::mapview(list(clusters_gmu1C, UDs_gmu1C_poly), zcol = "Clusters")
  UDs_gmu1W_poly <- bind_rows(ud_gmu1W_c1_noriver, ud_gmu1W_c2_noriver, ud_gmu1W_c3_noriver) %>%
    dplyr::select(-c(NewClusters, area))
  mapview::mapview(list(clusters_gmu1W, UDs_gmu1W_poly), zcol = "Clusters")

  #'  Calculate area (km^2) for each cluster polygon
  gmu1E_poly_area <- st_area(UDs_gmu1E_poly) %>% set_units(., km^2) %>% drop_units(.) %>% as.data.frame()
  names(gmu1E_poly_area) <- "area_km2"
  UDs_gmu1E_poly <- bind_cols(UDs_gmu1E_poly, gmu1E_poly_area) 
  mean(UDs_gmu1E_poly$area_km2)
  
  gmu1NE_poly_area <- st_area(UDs_gmu1NE_poly) %>% set_units(., km^2) %>% drop_units(.) %>% as.data.frame()
  names(gmu1NE_poly_area) <- "area_km2"
  UDs_gmu1NE_poly <- bind_cols(UDs_gmu1NE_poly, gmu1NE_poly_area) 
  mean(UDs_gmu1NE_poly$area_km2)
  
  gmu1C_poly_area <- st_area(UDs_gmu1C_poly) %>% set_units(., km^2) %>% drop_units(.) %>% as.data.frame()
  names(gmu1C_poly_area) <- "area_km2"
  UDs_gmu1C_poly <- bind_cols(UDs_gmu1C_poly, gmu1C_poly_area) 
  mean(UDs_gmu1C_poly$area_km2)
  
  gmu1W_poly_area <- st_area(UDs_gmu1W_poly) %>% set_units(., km^2) %>% drop_units(.) %>% as.data.frame()
  names(gmu1W_poly_area) <- "area_km2"
  UDs_gmu1W_poly <- bind_cols(UDs_gmu1W_poly, gmu1W_poly_area) 
  mean(UDs_gmu1W_poly$area_km2)
  
  ######  GMU6 Cluster Polygons  ######
  #'  ---------------------------
  #'  Intersect overlapping UDs
  UDs_gmu6_intersect <- st_intersection(UDs_gmu6) %>%
    mutate(NewClusters = row.names(.)) %>%
    dplyr::select(c(Clusters, NewClusters)) 
  mapviewOptions(fgb = TRUE)
  mapview::mapview(list(clusters_gmu6, UDs_gmu6_intersect), zcol = "Clusters")
  
  #'  Snag individual cluster polygon that's not changing
  #'  Other 3 cluster polygons overlap portions of polygon 3 but cluster 3 is most 
  #'  distinct from other clusters so want to retain the original polygon
  ud_gmu6_c3 <- filter(UDs_gmu6, Clusters == 3) %>% st_cast("POLYGON")
  
  #'  Create individual polygons remaining portions of other cluster polygons
  ud_gmu6_c1 <- filter(UDs_gmu6_intersect, NewClusters == 1) %>% st_cast("POLYGON")
  ud_gmu6_c2 <- filter(UDs_gmu6_intersect, NewClusters == 2) %>% st_cast("POLYGON")
  #'  Split amalgamation of Cluster 4 polygons into individual polygons
  ud_gmu6_c4 <- filter(UDs_gmu6_intersect, NewClusters == 4 | NewClusters == 1.2)  %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu6_c4[1:2,])
  #'  Join specific polygons from cluster 4
  ud_gmu6_c4 <- st_union(ud_gmu6_c4[2,], ud_gmu6_c4[1,]) %>% st_cast("POLYGON"); mapview(ud_gmu6_c4)

  #'  Spatial data frame of new cluster polygons
  UDs_gmu6_poly <- bind_rows(ud_gmu6_c1, ud_gmu6_c2, ud_gmu6_c3, ud_gmu6_c4) %>%
    dplyr::select(-c(NewClusters, Clusters.1, NewClusters.1, area))
  mapview::mapview(list(clusters_gmu6, UDs_gmu6_poly), zcol = "Clusters")
  
  #'  Calculate area (km^2) for each cluster polygon
  gmu6_poly_area <- st_area(UDs_gmu6_poly) %>% set_units(., km^2) %>% drop_units(.) %>% as.data.frame()
  names(gmu6_poly_area) <- "area_km2"
  UDs_gmu6_poly <- bind_cols(UDs_gmu6_poly, gmu6_poly_area)
  mean(UDs_gmu6_poly$area_km2)
  
  ######  GMU10A Cluster Polygons  ######
  #'  -----------------------------
  mapview::mapview(list(clusters_gmu10aN, UDs_gmu10aN), zcol = "Clusters")
  mapview::mapview(list(clusters_gmu10aS, UDs_gmu10aS), zcol = "Clusters")
    
  #'  Intersect overlapping UDs
  UDs_gmu10aN_intersect <- st_intersection(UDs_gmu10aN) %>%
    mutate(NewClusters = row.names(.)) %>%
    dplyr::select(c(Clusters, NewClusters)) 
  mapview::mapview(list(clusters_gmu10aN, UDs_gmu10aN_intersect), zcol = "Clusters")
  UDs_gmu10aS_intersect <- st_intersection(UDs_gmu10aS) %>%
    mutate(NewClusters = row.names(.)) %>%
    dplyr::select(c(Clusters, NewClusters)) 
  mapview::mapview(list(clusters_gmu10aS, UDs_gmu10aS_intersect), zcol = "Clusters")
  
  #'  Snag individual cluster polygons that are stand alone
  ud_gmu10aN_c1 <- filter(UDs_gmu10aN, Clusters == 1)
  ud_gmu10aS_c4 <- filter(UDs_gmu10aS, Clusters == 4)
  ud_gmu10aS_c6 <- filter(UDs_gmu10aS, Clusters == 6) %>% st_cast("POLYGON") 
  
  #'  Create individual polygons for places where two cluster polygons intersect
  ud_gmu10aN_c2 <- filter(UDs_gmu10aN_intersect, NewClusters == 2)
  ud_gmu10aS_c1 <- filter(UDs_gmu10aS_intersect, NewClusters == 1 | NewClusters == 1.1 | NewClusters == 1.4) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu10aS_c1[c(1,5:6),])
  ud_gmu10aS_c2 <- filter(UDs_gmu10aS_intersect, NewClusters == 2 | NewClusters == 2.1 | NewClusters == 2.2 | NewClusters == 1.2) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu10aS_c2)
  ud_gmu10aS_c3 <- filter(UDs_gmu10aS_intersect, NewClusters == 3 | NewClusters == 3.1 | NewClusters == 1.3) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu10aS_c3)
  ud_gmu10aS_c5 <- filter(UDs_gmu10aS_intersect, NewClusters == 5) %>% dplyr::select(Clusters)
  ud_gmu10aS_c7 <- filter(UDs_gmu10aS_intersect, NewClusters == 7 | NewClusters == 1.5) %>% 
    st_cast("MULTILINESTRING") %>% st_cast("LINESTRING") %>% st_cast("POLYGON"); mapview(ud_gmu10aS_c7)
  
  #'  Create new polygons by joining specific intersections
  ud_gmu10aS_c1 <- st_union(ud_gmu10aS_c1[1,], ud_gmu10aS_c1[5,]) %>% st_union(., ud_gmu10aS_c1[6,]) %>% 
    dplyr::select(Clusters) %>% st_cast(., "GEOMETRYCOLLECTION") %>% 
    st_collection_extract("POLYGON")
  mapview(ud_gmu10aS_c1)
  row.names(ud_gmu10aS_c1) <- 1
  ud_gmu10aS_c2 <- st_union(ud_gmu10aS_c2[1,], ud_gmu10aS_c2[2,]) %>% st_union(., ud_gmu10aS_c2[3,]) %>% 
    st_union(., ud_gmu10aS_c2[4,]) %>% st_union(., ud_gmu10aS_c2[6,]) %>% 
    # st_cast(., "GEOMETRYCOLLECTION") %>% st_collection_extract("POLYGON") %>% 
    dplyr::select(Clusters); mapview(ud_gmu10aS_c2)
  ud_gmu10aS_c3 <- st_union(ud_gmu10aS_c3[3,], ud_gmu10aS_c3[2,]) %>% st_union(., ud_gmu10aS_c3[4,]) %>% 
    st_union(., ud_gmu10aS_c3[1,]) %>% st_cast("POLYGON") %>%
    dplyr::select(Clusters)
  ud_gmu10aS_c3 <- ud_gmu10aS_c3[1,]; mapview(ud_gmu10aS_c3)
  row.names(ud_gmu10aS_c3) <- 3
  ud_gmu10aS_c7 <- st_union(ud_gmu10aS_c7[2,], ud_gmu10aS_c7[1,]) %>% 
    dplyr::select(Clusters); mapview(ud_gmu10aS_c7)
  
  #'  Remove large waterbodies from to GMU10A UDs
  ud_gmu10aN_c2_dworsak_split <- ud_gmu10aN_c2 %>% lwgeom::st_split(dworshak) %>% 
    st_collection_extract("POLYGON") 
  ud_gmu10aN_c2 <- ud_gmu10aN_c2_dworsak_split[1,]; mapview(ud_gmu10aN_c2)
  UDs_gmu10aS_poly_nowater <- bind_rows(ud_gmu10aS_c3, ud_gmu10aS_c4, ud_gmu10aS_c7) %>% 
    st_cast("MULTIPOLYGON") %>% dplyr::select(Clusters) %>% 
    lwgeom::st_split(dworshak) %>% st_collection_extract("POLYGON") %>% 
    lwgeom::st_split(clearwater) %>% st_collection_extract("POLYGON"); mapview(UDs_gmu10aS_poly_nowater)
  UDs_gmu10aS_poly_nowater <- bind_rows(UDs_gmu10aS_poly_nowater[16,], UDs_gmu10aS_poly_nowater[5,], UDs_gmu10aS_poly_nowater[24,])
  row.names(UDs_gmu10aS_poly_nowater) <- c(4, 3, 7)
  mapview(UDs_gmu10aS_poly_nowater, zcol = "Clusters")
  
  #'  Spatial data frame of new cluster polygons
  UDs_gmu10aN_poly <- bind_rows(ud_gmu10aN_c1, ud_gmu10aN_c2) %>%
    dplyr::select(-c(NewClusters, area))
  mapview::mapview(list(clusters_gmu10aN, UDs_gmu10aN_poly), zcol = "Clusters")
  UDs_gmu10aS_poly <- bind_rows(UDs_gmu10aS_poly_nowater, ud_gmu10aS_c1, ud_gmu10aS_c2, ud_gmu10aS_c5, ud_gmu10aS_c6) %>%
    dplyr::select(-area)
  mapview::mapview(list(clusters_gmu10aS, UDs_gmu10aS_poly), zcol = "Clusters")
  
  #'  Calculate area (km^2) for each cluster polygon
  gmu10aN_poly_area <- st_area(UDs_gmu10aN_poly) %>% set_units(., km^2) %>% drop_units(.) %>% as.data.frame()
  names(gmu10aN_poly_area) <- "area_km2"
  UDs_gmu10aN_poly <- bind_cols(UDs_gmu10aN_poly, gmu10aN_poly_area) 
  mean(UDs_gmu10aN_poly$area_km2)
  
  gmu10aS_poly_area <- st_area(UDs_gmu10aS_poly) %>% set_units(., km^2) %>% drop_units(.) %>% as.data.frame()
  names(gmu10aS_poly_area) <- "area_km2"
  UDs_gmu10aS_poly <- bind_cols(UDs_gmu10aS_poly, gmu10aS_poly_area) 
  mean(UDs_gmu10aS_poly$area_km2)
  
  #'  ------------------------------------------------------
  ####  Relabel camera clusters based on adjusted polygons  ####
  #'  ------------------------------------------------------
  #'  Join GMU1 clusters into single data set
  gmu1_poly <- bind_rows(UDs_gmu1NE_poly, UDs_gmu1C_poly, UDs_gmu1W_poly) %>%
    mutate(Clusters = seq(1:nrow(.))); mapview(gmu1_poly, zcol = "Clusters")
  #'  Relabel GMU1 camera clusters
  cam_clusters_gmu1 <- bind_rows(clusters_gmu1NE, clusters_gmu1C, clusters_gmu1W) %>%
    st_intersection(gmu1_poly) %>%
    rename("Clusters_original" = "Clusters") %>%
    rename("Clusters" = "Clusters.1") 
  mapview(list(gmu1_poly, cam_clusters_gmu1), zcol = "Clusters")
  
  #'  Join GMU6 clusters into single data set
  gmu6_poly <- UDs_gmu6_poly %>%
    mutate(Clusters = seq(1:nrow(.))); mapview(gmu6_poly, zcol = "Clusters")
  #'  Relabel GMU6 camera clusters
  cam_clusters_gmu6 <- st_intersection(clusters_gmu6, gmu6_poly) %>%
    rename("Clusters_original" = "Clusters") %>%
    rename("Clusters" = "Clusters.1") 
  #'  Grab polygon areas
  c2_area <- unique(cam_clusters_gmu6$area_km2[cam_clusters_gmu6$Clusters == 2])
  #'  GMU6_U_125 gets left out during intersection 
  #'  (falls slightly outside of Cluster 2 polygon b/c clusters based on 95% KDE)
  GMU6_U_125 <- filter(clusters_gmu6, NwLctID == "GMU6_U_125") %>%
    rename("Clusters_original" = "Clusters") %>%
    mutate(Clusters = Clusters_original,
           area_km2 = c2_area) %>%
    relocate(Clusters, .after = Clusters_adjacency) %>%
    relocate(area_km2, .after = Clusters)
  cam_clusters_gmu6 <- bind_rows(cam_clusters_gmu6, GMU6_U_125)
  mapview(list(gmu6_poly, cam_clusters_gmu6), zcol = "Clusters")
  
  #'  Join GMU10A clusters into single data set
  gmu10a_poly <- bind_rows(UDs_gmu10aN_poly, UDs_gmu10aS_poly) %>%
    mutate(Clusters = seq(1:nrow(.))); mapview(gmu1_poly, zcol = "Clusters")
  #'  Relabel GMU10A camera clusters
  cam_clusters_gmu10a <- bind_rows(clusters_gmu10aN, clusters_gmu10aS) %>%
    st_intersection(gmu10a_poly) %>%
    rename("Clusters_original" = "Clusters") %>%
    rename("Clusters" = "Clusters.1")
  mapview(list(gmu10a_poly, cam_clusters_gmu10a), zcol = "Clusters")
  
  #'  -----------------------------------------------------------
  ####  Calculate Relative Density Index for wolves per cluster  ####
  #'  -----------------------------------------------------------
  #  Calculate density per species per year per cluster
  wolf_density_gmu1 <- cam_clusters_gmu1 %>%
    group_by(Clusters) %>%
    reframe(nWolf = sum(RN_n_rn),
            `Wolf density` = nWolf/area_km2) %>%
    unique() %>%
    mutate(GMU = "GMU1") %>%
    full_join(gmu1_poly, by = "Clusters") %>%
    st_as_sf()
  wolf_density_gmu6 <- cam_clusters_gmu6 %>%
    group_by(Clusters) %>%
    reframe(nWolf = sum(RN_n_rn),
            `Wolf density` = nWolf/area_km2) %>%
    unique() %>%
    mutate(GMU = "GMU6") %>%
    full_join(gmu6_poly, by = "Clusters") %>%
    st_as_sf()
  wolf_density_gmu10a <- cam_clusters_gmu10a %>%
    group_by(Clusters) %>%
    reframe(nWolf = sum(RN_n_rn),
            `Wolf density` = nWolf/area_km2) %>%
    unique() %>%
    mutate(GMU = "GMU10A") %>%
    full_join(gmu10a_poly, by = "Clusters") %>%
    st_as_sf()
  
  wolf_density_tbl <- bind_rows(wolf_density_gmu1, wolf_density_gmu6, wolf_density_gmu10a) %>%
    as.data.frame(.) %>% 
    mutate(area_km2 = round(area_km2, 2),
           `Wolf density` = round(`Wolf density`, 3)) %>%
    relocate(area_km2, .after = nWolf) %>%
    relocate(GMU, .before = Clusters) %>%
    dplyr::select(-geometry)
  write_csv(wolf_density_tbl, "./Outputs/Relative_Abundance/RN_model/Cluster_wolf_RAI_density.csv")
  
  #'  Back to allowing spherical geometry
  sf::sf_use_s2(TRUE)
  
  
  #'  Save cluster data as shapefiles
  st_write(cam_clusters_gmu1, "./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu1.shp")
  st_write(cam_clusters_gmu6, "./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu6.shp")
  st_write(cam_clusters_gmu10a, "./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu10a.shp")
  st_write(gmu1_poly, "./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu1.shp")
  st_write(gmu6_poly, "./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu6.shp")
  st_write(gmu10a_poly, "./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu10a.shp")
  
  
  #'  Visualize with ggplot and save
  clusters_gmu1 <- clusters_gmu1 %>% dplyr::select(Clusters) %>% mutate(Clusters = as.character(Clusters))
  clusters_gmu6 <- clusters_gmu6 %>% dplyr::select(Clusters) %>% mutate(Clusters = as.character(Clusters))
  clusters_gmu10a <- clusters_gmu10a %>% dplyr::select(Clusters) %>% mutate(Clusters = as.character(Clusters))
  gmu1_clusters_og <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 1,], fill = NA) +
    geom_sf(data = UDs_gmu1W, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = UDs_gmu1C, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = UDs_gmu1NE, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = clusters_gmu1, aes(col = Clusters, fill = Clusters), 
            shape  = 21, size = 2.5) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU1 Camera clusters, original")
  gmu1_clusters_final <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 1,], fill = NA) +
    geom_sf(data = gmu1_poly, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = cam_clusters_gmu1, aes(col = Clusters, fill = Clusters), 
            shape  = 21, size = 2.5) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU1 Camera clusters, final")
  (gmu1_wolf_density <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 1,], fill = NA) +
    geom_sf(data = wolf_density_gmu1, aes(fill = `Wolf density`), alpha = 0.2) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU1 relative wolf density index (wolves/km^2)"))
  
  gmu6_clusters_og <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 6,], fill = NA) +
    geom_sf(data = UDs_gmu6, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = clusters_gmu6, aes(col = Clusters, fill = Clusters), 
            shape  = 21, size = 2.5) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU6 Camera clusters, original")
  gmu6_clusters_final <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 6,], fill = NA) +
    geom_sf(data = gmu6_poly, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = cam_clusters_gmu6, aes(col = Clusters, fill = Clusters), 
            shape  = 21, size = 2.5) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU6 Camera clusters, final")
  (gmu6_wolf_density <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == 6,], fill = NA) +
    geom_sf(data = wolf_density_gmu6, aes(fill = `Wolf density`), alpha = 0.2) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU6 relative wolf density index (wolves/km^2)"))
  
  gmu10a_clusters_og <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "10A",], fill = NA) +
    geom_sf(data = UDs_gmu10aN, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = UDs_gmu10aS, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = clusters_gmu10a, aes(col = Clusters, fill = Clusters), 
            shape  = 21, size = 2.5) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU10A Camera clusters, original")
  gmu10a_clusters_final <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "10A",], fill = NA) +
    geom_sf(data = gmu10a_poly, aes(fill = Clusters), alpha = 0.2) +
    geom_sf(data = cam_clusters_gmu10a, aes(col = Clusters, fill = Clusters), 
            shape  = 21, size = 2.5) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU10A Camera clusters, final")
  (gmu10a_wolf_density <- ggplot() +
    geom_sf(data = eoe_gmu_wgs84[eoe_gmu_wgs84$NAME == "10A",], fill = NA) +
    geom_sf(data = wolf_density_gmu10a, aes(fill = `Wolf density`), alpha = 0.2) +
    labs(x = "Longitude", y = "Latitude") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GMU10A relative wolf density index (wolves/km^2)"))
  
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU1_clusters_og.tiff", gmu1_clusters_og,
         units = "in", width = 5, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU1_clusters_final.tiff", gmu1_clusters_final,
         units = "in", width = 5, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU1_wolf_density.tiff", gmu1_wolf_density,
         units = "in", width = 5, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU6_clusters_og.tiff", gmu6_clusters_og,
         units = "in", width = 6, height = 5, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU6_clusters_final.tiff", gmu6_clusters_final,
         units = "in", width = 6, height = 5, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU6_wolf_density.tiff", gmu6_wolf_density,
         units = "in", width = 5, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU10A_clusters_og.tiff", gmu10a_clusters_og,
         units = "in", width = 5, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU10A_clusters_final.tiff", gmu10a_clusters_final,
         units = "in", width = 5, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/GMU10A_wolf_density.tiff", gmu10a_wolf_density,
         units = "in", width = 5, height = 6, dpi = 300, device = "tiff", compression = "lzw")
  
  
  
  #  SEM your heart out
  
  
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  