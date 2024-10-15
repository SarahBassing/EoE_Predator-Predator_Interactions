  #'  --------------------------------
  #'  Prep Data for Structural Equation Models
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2024
  #'  --------------------------------
  #'  Format and explore data before building SEMs
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(piecewiseSEM)
  library(labelled)
  library(lme4)
  library(MASS)
  #install.packages("multcompView", repos="http://R-Forge.R-project.org")
  # library(multcompView)
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
  library(ggplot2)
  
  #'  Load RN model local abundance estimates
  load("./Outputs/Relative_Abundance/RN_model/RN_abundance.RData")
  
  #'  Function to clean up camera cluster spatial data
  reformat_clusters <- function(clusters) {
    skinny_clusters <- clusters %>%
      dplyr::select(c(NwLctID, CellID, GMU, Setup, Clustrs, are_km2, geometry)) %>%
      rename(NewLocationID = NwLctID) %>%
      rename(ClusterID = Clustrs) %>%
      rename(area_km2 = are_km2)
    return(skinny_clusters)
  }
  
  #'  Read in camera cluster shapefiles and clean up with reformat_clusters() function
  clusters_gmu1 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu1.shp") %>% 
    reformat_clusters(.)
  clusters_gmu6 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu6.shp") %>% 
    reformat_clusters(.)
  clusters_gmu10a <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu10a.shp") %>% 
    reformat_clusters(.)
  
  #'  Merge spatial camera cluster data together and reduce to one observation per camera
  clusters_all <- bind_rows(clusters_gmu1, clusters_gmu6, clusters_gmu10a) %>%
    group_by(NewLocationID) %>%
    slice(1L) %>%
    ungroup()
  
  #'  Read in cluster polygon shapefiles
  gmu1_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu1.shp") %>% mutate(GMU = "GMU1")
  gmu6_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu6.shp") %>% mutate(GMU = "GMU6")
  gmu10a_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu10a.shp") %>% mutate(GMU = "GMU10A")
  
  #'  Merge cluster polygons across GMUs
  cluster_poly <- bind_rows(gmu1_poly, gmu6_poly, gmu10a_poly) %>%
    rename(ClusterID = Clusters)
  mapview::mapview(list(cluster_poly, clusters_all), zcol = "ClusterID")
  
  #'  Read in EoE GMUs shapefile
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Read in forested landcover
  perc_forest <- rast("./Shapefiles/National Land Cover Database (NCLD)/PercentForest_100m.tif")
  forest <- rast("./Shapefiles/National Land Cover Database (NCLD)/Forested_landcover.tif")
  notforest <- rast("./Shapefiles/National Land Cover Database (NCLD)/NonForested_landcover.tif")
  forest_proj <- crs(forest)
  
  #' #'  Calculate average percent forested habitat across cameras in each cluster
  #' extract_forest <- function(cams, gmu) {
  #'   #'  Transform and reduce to a single observation per camera (don't need to 
  #'   #'  repeat information across years)
  #'   cams <- st_transform(cams, forest_proj) %>%
  #'     group_by(NewLocationID) %>%
  #'     slice(1L) %>%
  #'     ungroup()
  #'   #'  Extract % forested habitat within 100m radius of each camera site and
  #'   #'  average across cameras per cluster
  #'   perc_forest_100m <- terra::extract(perc_forest, vect(cams)) %>%
  #'     dplyr::select(-ID) %>%
  #'     bind_cols(cams) %>%
  #'     rename("percent_forest" = "focal_sum") %>%
  #'     group_by(ClusterID) %>%
  #'     summarise(percent_forest = mean(percent_forest)) %>%
  #'     ungroup() %>%
  #'     mutate(GMU = gmu)
  #'   return(perc_forest_100m)
  #' }
  #' forest_gmu1 <- extract_forest(clusters_gmu1, gmu = "GMU1")
  #' forest_gmu6 <- extract_forest(clusters_gmu6, gmu = "GMU6")
  #' forest_gmu10a <- extract_forest(clusters_gmu10a, gmu = "GMU10A")
  #' 
  #' perc_forest <- bind_rows(forest_gmu1, forest_gmu6, forest_gmu10a) %>%
  #'   dplyr::select(c(GMU, ClusterID, percent_forest))
  
  #' #'  Vectorize polygons
  #' gmu1_vect <- st_transform(gmu1_poly, forest_proj) %>% vect(.)
  #' gmu6_vect <- st_transform(gmu6_poly, forest_proj) %>% vect(.)
  #' gmu10a_vect <- st_transform(gmu10a_poly, forest_proj) %>% vect(.)

  #'  Zonal statistics: sum amount of forest and non-forested land cover per cluster
  #'  and calculate % forested
  percent_forested <- function(polygon) {
    #'  Vectorize polygons
    vect_poly <- st_transform(polygon, forest_proj) %>% vect(.)
    #'  Snag cluster ID
    clusters <- as.data.frame(polygon) %>% dplyr::select(c(Clusters, GMU))
    #'  Sum pixels classified as forested land cover per cluster
    sum_forest <- zonal(forest, vect_poly, fun = "sum", na.rm = TRUE) %>%
      rename("total_forest" = "NLCD Land Cover Class")
    #'  Sum pixels classified as non-forested land cover per cluster
    sum_nonforest <- zonal(notforest, vect_poly, fun = "sum", na.rm = TRUE) %>%
      rename("total_nonforest" = "NLCD Land Cover Class")
    #'  Combine and calculate percent forested per cluster
    perc_forest <- cbind(sum_forest, sum_nonforest) %>% as.data.frame(.) %>%
      mutate(total_pix = total_forest + total_nonforest,
             percent_forest = total_forest/total_pix) %>%
      bind_cols(clusters) %>%
      rename("ClusterID" = "Clusters")
    return(perc_forest)
  }
  cluster_perc_forest_gmu1 <- percent_forested(gmu1_poly)
  cluster_perc_forest_gmu6 <- percent_forested(gmu6_poly)
  cluster_perc_forest_gmu10a <- percent_forested(gmu10a_poly)

  perc_forest <- bind_rows(cluster_perc_forest_gmu1, cluster_perc_forest_gmu6, cluster_perc_forest_gmu10a) %>%
    dplyr::select(c(GMU, ClusterID, percent_forest))
  
  #'  ----------------------------
  ####  Percent disturbed forest  ####
  #'  ----------------------------
  #'  Hansen's Global Forest Change dataset (area and year of canopy loss per cluster)
  gfc <- read_csv("./Data/GEE outputs/GFC_annual_canopy_loss_area_clusters.csv") %>% 
    dplyr::select(-`.geo`) %>%
    #'  Make canopy loss year easier to interpret
    mutate(CanopyLossYear = Year_add1 + 2001,
           #'  Calculate proportion of canopy loss within each cluster
           CanopyLossProp = CanopyLossArea_sq_m/BufferArea_sq_m,
           CanopyLossArea_sq_km = CanopyLossArea_sq_m/1000000,
           CanopyLossArea_sq_m = round(CanopyLossArea_sq_m, 3),
           CanopyLossProp = round(CanopyLossProp, 5),
           #'  Relabel GEE indexing to corresponding GMUs
           GMU = ifelse(str_detect(`system:index`, "_1_1_"), "GMU1", `system:index`),  # note the order is different with PRISM data
           GMU = ifelse(str_detect(GMU, "_1_2_"), "GMU6", GMU),
           GMU = ifelse(str_detect(GMU, "_2_"), "GMU10A", GMU)) %>%
    #'  Remove unnecessary columns and rows
    dplyr::select(c(GMU, ClusterID, CanopyLossYear, CanopyLossArea_sq_km, CanopyLossArea_sq_m, CanopyLossProp)) %>%
    filter(CanopyLossArea_sq_m > 0) %>% 
    #'  Retain loss year where largest proportion of area was lost (not necessarily 
    #'  the most recent - want to represent the "dominant" habitat type so even if 
    #'  smaller more recent loss occurred area is best represented by larger, older loss)
    group_by(GMU, ClusterID) %>%
    slice_max(order_by = CanopyLossProp, n = 1) %>% 
    ungroup() %>%
    arrange(GMU, ClusterID)
  
  #'  National Landcover Database (frequency of landcover class per cluster)
  nlcd <- read_csv("./Data/GEE outputs/NLCD_frequencies_2019_2021_clusters.csv") %>%
    #'  Relabel GEE indexing to corresponding GMUs
    mutate(GMU = ifelse(str_detect(`system:index`, "1_1_"), "GMU1", `system:index`),  # note the order is different with PRISM data
           GMU = ifelse(str_detect(GMU, "1_2_"), "GMU6", GMU),
           GMU = ifelse(str_detect(GMU, "2_"), "GMU10A", GMU)) %>%
    dplyr::select(c(GMU, Clusters, landcover_2019, landcover_2021)) %>%
    rename("ClusterID" = "Clusters") %>%
    mutate(landcover_2019 = str_replace_all(landcover_2019, "[[{}]]", ""),
           landcover_2021 = str_replace_all(landcover_2021, "[[{}]]", ""))
  nlcd19 <- nlcd %>% dplyr::select(c(GMU, ClusterID, landcover_2019)) %>%
    #'  Remove GMU1 cameras from this dataset
    filter(GMU != "GMU1") %>%
    #'  Separate wonky GEE character string grouping all landcover outputs together 
    #'  (FYI separate does something weird with extra columns so be sure to provide 
    #'  double the max number of possible landcover types). Ignore warning about missing pieces and NAs.
    separate(landcover_2019, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17", "col18", "col19", "col20", "col21", "col22"), sep = "[, ]")
  nlcd21 <- nlcd %>% dplyr::select(c(GMU, ClusterID, landcover_2021)) %>%
    separate(landcover_2021, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17", "col18", "col19", "col20", "col21", "col22"), sep = "[, ]")
  
  #'  Function to reformat NLCD data extracted from GEE
  habitat_frequencies <- function(dat) {
    reformatted <- dat %>%
      #'  Convert to long format to organize columns better
      pivot_longer(!c(GMU, ClusterID), names_to = "col_name", values_to = "landcover_class") %>%
      #'  Split landcover class from number of pixels
      separate(landcover_class, into = c("NLCD_class", "npixels"), sep = "[=]") %>%
      #'  Drop unneeded columns and rows
      dplyr::select(-col_name) %>%
      filter(!is.na(npixels)) %>%
      #'  Reformat npixels 
      mutate(npixels = as.numeric(npixels)) %>%
             # npixels = round(npixels, 2)) %>%
      #'  Calculate proportion of pixels made up by each cover type near camera
      group_by(GMU, ClusterID) %>%
      reframe(NLCD_class = NLCD_class,
              npixels = npixels,
              totalPix = sum(npixels),
              pixel_area_m2 = totalPix * (30*30),
              pixel_area_km2 = pixel_area_m2/1000000,
              percentPix = npixels/totalPix) %>%
      ungroup() %>%
      mutate(percentPix = round(percentPix, 3)) %>%
      arrange(GMU, ClusterID) %>%
      # rename("NewLocationID" = "NwLctID") %>%
      #'  Define NLCD landcover classifications using https://developers.google.com/earth-engine/datasets/catalog/USGS_NLCD_RELEASES_2021_REL_NLCD
      mutate(habitat_type = ifelse(NLCD_class == "11", "Other", NLCD_class),               #11: Open water
             habitat_type = ifelse(NLCD_class == "12", "Other", habitat_type),             #12: Perennial ice/snow
             habitat_type = ifelse(NLCD_class == "21", "Developed", habitat_type),         #21: Developed, open space: mixture of constructed materials, mostly vegetation in the form of lawn grasses
             habitat_type = ifelse(NLCD_class == "22", "Developed", habitat_type),         #22: Developed, low intensity: areas with a mixture of constructed materials and vegetation, most commonly single-family housing units
             habitat_type = ifelse(NLCD_class == "23", "Developed", habitat_type),         #21: Developed, medium intensity: areas with a mixture of constructed materials and vegetation, most commonly single-family housing units
             habitat_type = ifelse(NLCD_class == "24", "Developed", habitat_type),         #22: Developed, high intensity: highly developed areas where people reside or work in high numbers
             habitat_type = ifelse(NLCD_class == "31", "Other", habitat_type),             #31: Barren Land
             habitat_type = ifelse(NLCD_class == "41", "Forested", habitat_type),          #41: Deciduous Forest
             habitat_type = ifelse(NLCD_class == "42", "Forested", habitat_type),          #42: Evergreen Forest
             habitat_type = ifelse(NLCD_class == "43", "Forested", habitat_type),          #43: Mixed Forest
             habitat_type = ifelse(NLCD_class == "52", "Shrubland", habitat_type),         #52: Shrub/Scrub
             habitat_type = ifelse(NLCD_class == "71", "Grassland", habitat_type),         #71: Grassland/Herbaceous
             habitat_type = ifelse(NLCD_class == "81", "Agriculture", habitat_type),       #81: Pasture/Hay
             habitat_type = ifelse(NLCD_class == "82", "Agriculture", habitat_type),       #82: Cultivated Crops
             habitat_type = ifelse(NLCD_class == "90", "Riparian_woodland", habitat_type),      #90: Woody Wetlands
             habitat_type = ifelse(NLCD_class == "95", "Riparian_wetland", habitat_type)) %>%   #95: Emergent Herbaceous Wetlands)
      relocate(habitat_type, .after = "ClusterID")
    
    #'  Filter data to the dominant habitat class within each cluster
    #'  (dominant is defined as the landcover class comprising the largest proportion
    #'  of pixels within cluster)
    dominant_habitat <- reformatted %>%
      group_by(GMU, ClusterID) %>%
      slice_max(order_by = percentPix, n = 1) %>%
      ungroup() %>%
      arrange(GMU, ClusterID)
    
    #'  Review dominant habitat types
    print(table(dominant_habitat$habitat_type))
    
    #'  List both datasets together
    landcover_list <- list(dominant_habitat, reformatted)
    names(landcover_list) <- c("dominant_habitat", "percent_landcover")
    
    return(landcover_list)
    return()
  }
  landcover19 <- habitat_frequencies(nlcd19)
  landcover21 <- habitat_frequencies(nlcd21)
  
  #'  Load MTBS burn perimeter data 
  mtbs <- st_read(dsn = "./Shapefiles/MTBS_perimeter_data", layer = "mtbs_perims_DD") %>%
    st_transform(forest_proj)
  
  #'  Extract burn year of each pixel per cluster
  burn_perimeters <- function(polygon) {
    #'  Vectorize polygons
    polygon <- st_transform(polygon, forest_proj) #%>% vect(.)
    #'  Intersect overlapping clusters & burn perimeters
    mtbs_simple <- mtbs %>% 
      dplyr::select("Ig_Date", "Incid_Name", "Incid_Type") %>%
      mutate(Burn_year = format(as.Date(Ig_Date, format="%Y-%m-%d"),"%Y")) 
    mtbs_cluster_inter <- st_intersection(polygon, mtbs_simple) 
    burned_area <- as.numeric(st_area(mtbs_cluster_inter)/1000000) %>%
      as.data.frame()
    names(burned_area) <- "burned_area_km2"
    mtbs_cluster_inter <- bind_cols(mtbs_cluster_inter, burned_area)
    mapview::mapview(mtbs_cluster_inter, zcol = "Clusters")
    #'  Save burn year and area
    mtbs_burn_yr <- as.data.frame(mtbs_cluster_inter) %>%
      dplyr::select(c(GMU, Clusters, Burn_year, area_km2, burned_area_km2)) %>%
      rename("ClusterID" = "Clusters") 
    return(mtbs_burn_yr)
  }
  burned_gmu1 <- burn_perimeters(gmu1_poly)
  burned_gmu6 <- burn_perimeters(gmu6_poly)
  burned_gmu10a <- burn_perimeters(gmu10a_poly)
 
  
  
  #'  Add GEE data to larger covariate df
  format_cam_covs <- function(dat, landcov, season, camYr) {
    covs <- full_join(dat, gfc, by = "ClusterID") %>%
      full_join(gfc, by = "ClusterID") %>%
      relocate(Burn_year, .after = "percentPix") %>%
      mutate(Season = season,
             #'  Calculate number of years since burn/canopy loss
             YrsSinceBurn = camYr - as.numeric(Burn_year),
             YrsSinceBurn = ifelse(YrsSinceBurn <0, NA, YrsSinceBurn),
             YrsSinceLoss = camYr - CanopyLossYear,
             YrsSinceLoss = ifelse(YrsSinceLoss <0, NA, YrsSinceLoss),
             #'  Categorize years since burn/canopy loss following Barker et al. (2018) and Ganz et al. (2024)
             DisturbanceLoss = ifelse(YrsSinceLoss <= 20, "Loss_1_20", habitat_type),
             DisturbanceBurn = ifelse(YrsSinceBurn <= 5, "Burn_1_5", habitat_type),
             DisturbanceBurn = ifelse(YrsSinceBurn > 5 & YrsSinceBurn <=10, "Burn_6_10", DisturbanceBurn),
             DisturbanceBurn = ifelse(YrsSinceBurn > 10 & YrsSinceBurn <=15, "Burn_10_15", DisturbanceBurn),
             DisturbanceBurn = ifelse(YrsSinceBurn > 15 & YrsSinceBurn <=20, "Burn_16_20", DisturbanceBurn),
             DisturbanceBurn = ifelse(YrsSinceBurn > 20, "Burn_over20", DisturbanceBurn),
             #'  Generate a single habitat class covariate representing dominant habitat type and years since burn/canopy loss (if forested)
             Habitat_class = habitat_type,
             Habitat_class = ifelse(!is.na(DisturbanceLoss) & Habitat_class == "Forested", DisturbanceLoss, Habitat_class),
             Habitat_class = ifelse(!is.na(DisturbanceBurn) & Habitat_class == "Loss_1_20", DisturbanceBurn, Habitat_class),
             #'  Grab percentPix for sites with disturbance to indicate percent disturbed forest within last 20 years
             #'  (0 % disturbed forest if the dominant habitat class is unburned/unlogged forest or any other landcover type)
             PercDisturbedForest = ifelse(Habitat_class == "Loss_1_20" | Habitat_class == "Burn_1_5" | 
                                            Habitat_class == "Burn_6_10" | Habitat_class == "Burn_16_20", percentPix, 0)) %>%
      left_join(wsi, by = c("NewLocationID", "Season")) %>%
      relocate(Season, .after = "GMU") %>%
      filter(!is.na(GMU))
    
    #'  Review new habitat classes
    print(table(covs$Habitat_class))
    
    return(covs)
  }
  #'  Add GEE data to larger covariate dataframe; NOTE different landcover data applied to 2020 vs 2021 & 2022 data
  cams_eoe20s <- format_cam_covs(covariate_list[[1]], landcov = landcover19[[1]], camYr = 2020, year = "2020") 
  cams_eoe21s <- format_cam_covs(covariate_list[[2]], landcov = landcover21[[1]], camYr = 2021, year = "2021") 
  cams_eoe22s <- format_cam_covs(covariate_list[[3]], landcov = landcover21[[1]], camYr = 2022, year = "2022") 
  
  #'  List annual covariate data
  cam_covs_list <- list(cams_eoe20s, cams_eoe21s, cams_eoe22s)
  
  #'  ---------------------------------------------------
  ####  Calculate index of relative density per species  ####
  #'  ---------------------------------------------------
  #'  Join local abundance estimates with spatial cluster data
  cluster_RAI <- function(rai, clusters) {
    clustered_rai <- rai %>%
      left_join(clusters, by = c("NewLocationID", "CellID", "GMU", "Setup"))
    return(clustered_rai)
  }
  RN_abundance_sf <- lapply(RN_abundance, cluster_RAI, clusters = clusters_all)
  
  #'  Calculate index of relative density per species per cluster per year
  density_per_cluster <- function(rai, yr) {
    relative_density <- rai %>%
      group_by(GMU, ClusterID, Species) %>%
      reframe(SppN = sum(RN.n),
              SppDensity.km2 = SppN/area_km2,
              SppDensity.100km2 = SppDensity.km2*100,
              SppN.r = sum(round(RN.n, 0)),
              SppDensity.km2.r = SppN.r/area_km2,
              SppDensity.100km2.r = SppDensity.km2.r*100) %>% #,
      # NewLocationID = NewLocationID,
      # RN.n = RN.n) %>%
      unique() %>%
      mutate(Year = yr) %>%
      #'  Join spatial polygon data to relative density estimates
      left_join(cluster_poly, by = c("GMU", "ClusterID")) %>%
      st_as_sf()
    return(relative_density)
  }
  yr <- list("2020", "2021", "2022")
  cluster_density <- mapply(density_per_cluster, rai = RN_abundance_sf, yr = yr, SIMPLIFY = FALSE)
  
  #'  ---------------------------------
  ####  Explore density relationships  ####
  #'  ---------------------------------
  #'  Visualize with ggplot
  dat <- bind_rows(cluster_density[[1]], cluster_density[[2]], cluster_density[[3]])
  (wolf_density <- dat %>% filter(Species == "wolf") %>%
      ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (wolves/100km^2)"))
  (bear_density <- dat %>% filter(Species == "bear_black") %>%
      ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (black bear/100km^2)"))
  (lion_density <- dat %>% filter(Species == "mountain_lion") %>%
      ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "green") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (mountain lion/100km^2)"))
  (coy_density <- dat %>% filter(Species == "coyote") %>%
      ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "orange") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (coyote/100km^2)"))
  (elk_density <- dat %>% filter(Species == "elk") %>%
      ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "purple") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (elk/100km^2)"))
  (moose_density <- dat %>% filter(Species == "moose") %>%
      ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "gold") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (moose/100km^2)"))
  (wtd_density <- dat %>% filter(Species == "whitetailed_deer") %>%
      ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "black") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (white-tailed deer/100km^2)"))
  
  # ggsave("./Outputs/Relative_Abundance/RN_model/Figures/cluster_wolf_density.tiff", wolf_density,
  #        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  # ggsave("./Outputs/Relative_Abundance/RN_model/Figures/cluster_bear_density.tiff", bear_density,
  #        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  # ggsave("./Outputs/Relative_Abundance/RN_model/Figures/cluster_lion_density.tiff", lion_density,
  #        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  # ggsave("./Outputs/Relative_Abundance/RN_model/Figures/cluster_coy_density.tiff", coy_density,
  #        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  # ggsave("./Outputs/Relative_Abundance/RN_model/Figures/cluster_elk_density.tiff", elk_density,
  #        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  # ggsave("./Outputs/Relative_Abundance/RN_model/Figures/cluster_moose_density.tiff", moose_density,
  #        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  # ggsave("./Outputs/Relative_Abundance/RN_model/Figures/cluster_wtd_density.tiff", wtd_density,
  #        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  
  
  #'  ---------------------------------
  ####  Google Earth Engine data sets  ####
  #'  ---------------------------------
  #'  PRISM weather data
  #'  Source script to load & format average monthly precip & temp data
  #'  Produces standardized average monthly total precipitation and standardized 
  #'  average monthly minimum temp for each winter and GMU for past 50 years
  source("./Scripts/Structural_Equation_Models/Format_weather_data.R")   
  
  #'  Generate Winter Severity Index (WSI) per year and Cluster
  wsi <- full_join(wtr_totalPrecip, wtr_minTemp, by = c("GMU", "Clusters", "Season")) %>%
    #'  Multiply standardized precip by standardized temp data
    mutate(DecFeb_WSI = DecFeb_meanPPT_z * DecFeb_meanMinTemp_z) %>%
    dplyr::select(c("GMU", "Clusters", "Season", "DecFeb_WSI")) %>%
    #'  Filter to the winters preceding Summer 2020, 2021, and 2022 surveys
    filter(Season == "Wtr1920" | Season == "Wtr2021" | Season == "Wtr2122") %>%
    #'  Change Season column to reference annual summer sampling
    rename("PreviousWinter" = "Season") %>%
    mutate(year = ifelse(PreviousWinter == "Wtr1920", "yr1", "yr3"),
           year = ifelse(PreviousWinter == "Wtr2021", "yr2", year)) %>%
    relocate(year, .after = "Clusters") %>%
    rename("ClusterID" = "Clusters")
  print(wsi)
  
  
  #'  -----------------------------------------
  ####  Format local abundance estimates & SD  ####
  #'  -----------------------------------------
  #'  Join covariates together
  covs <- full_join(perc_forest, wsi, by = c("GMU", "ClusterID")) 
  
  #'  Long to wide data structure
  wide_data <- function(dat) {
    pivot_data_wide <- as.data.frame(dat) %>%
      #'  Create categorical year variable
      mutate(year = ifelse(Year == 2020, "yr1", Year),
             year = ifelse(Year == 2021, "yr2", year),
             year = ifelse(Year == 2022, "yr3", year)) %>% #,
             #' #'  Add small value to density estimates that are 0 (needed for boxcox trnasformatio below)
             #' SppDensity.100km2.r = ifelse(SppDensity.100km2.r == 0, 0.00001, SppDensity.100km2.r)) %>%
      dplyr::select(c(GMU, ClusterID, year, Species, SppDensity.100km2.r)) %>%
      #'  Create column per species with their site-specific local abundance 
      pivot_wider(names_from = "Species",
                  values_from = c("SppDensity.100km2.r")) %>%
      left_join(covs, by = c("GMU", "year", "ClusterID")) %>%
      dplyr::select(-PreviousWinter) %>%
      relocate(percent_forest, .after = year) %>%
      relocate(DecFeb_WSI, .after = percent_forest)
    return(pivot_data_wide)
  }
  density_wide <- lapply(cluster_density, wide_data)
  
  #'  Group data across years so t-1 data affects t data, regardless of year
  #'  Stack year 1 & year 2 data
  dat_t_minus_1 <- bind_rows(density_wide[[1]], density_wide[[2]])
  #'  Stack year 2 & year 3 data
  dat_t <- bind_rows(density_wide[[2]], density_wide[[3]])
  #'  List stacked data
  dat_stack_list <- list(dat_t_minus_1, dat_t)
  
  #'  Group all years together (stacking Year1, Year2, and Year3 together)
  density_wide_allyrs <- rbind(density_wide[[1]], density_wide[[2]], density_wide[[3]]) 
  
  #'  Visualize spread of relative density estimates across clusters per year
  plot_histogram <- function(dat, yr) {
    hist(dat$bear_black, main = paste("Relative bear density across clusters,", yr))
    hist(dat$coyote, main = paste("Relative coyote density across clusters,", yr))
    hist(dat$mountain_lion, main = paste("Relative lion density across clusters,", yr))
    hist(dat$wolf, main = paste("Relative wolf density across clusters,", yr))
    hist(dat$elk, main = paste("Relative elk density across clusters,", yr))
    hist(dat$moose, main = paste("Relative moose density across clusters,", yr))
    hist(dat$whitetailed_deer, main = paste("Relative wtd density across clusters,", yr))
    hist(dat$percent_forest, main = paste("Percent forested land cover,", yr))
    hist(dat$DecFeb_WSI, main = paste("Dec-Feb Winter Severity Index,", yr))
  }
  plot_histogram(density_wide[[1]], yr = "2020")
  plot_histogram(density_wide[[2]], yr = "2021")
  plot_histogram(density_wide[[3]], yr = "2022")
  
  #'  ----------------------------------------------------------
  ####  Format & transform data for different model structures  ####
  #'  ----------------------------------------------------------
  #'  ----------------------------------------------
  #####  Annual format: independent annual columns  #####
  #'  ----------------------------------------------
  #'  Retaining each year as independent variable allows SEMs to assess indirect effects
  #'  and let relationships vary annually
  #'  ----------------------------------------------
  #'  Create wide data structure but this time one column per year for each species & covariate
  wide_data_by_year <- function(dat, yr) {
    data_by_yr <- dat %>%
      #'  Add year identifier to each column name
      rename_with(.cols = percent_forest:wolf, function(x){paste0(x, ".", yr)}) %>% 
      dplyr::select(-year)
    return(data_by_yr)
  }
  density_wide_annual <- mapply(wide_data_by_year, dat = density_wide, yr = list("yr1", "yr2", "yr3"), SIMPLIFY = FALSE) #dat = full_bx_dat
  #'  Sneak peak of each year
  head(density_wide_annual[[1]])
  head(density_wide_annual[[2]])
  head(density_wide_annual[[3]])
  
  #'  Unlist as one single data frame (annual columns per species and covariates)
  density_wide_annual_20s_22s <- full_join(density_wide_annual[[1]], density_wide_annual[[2]], by = c("GMU", "ClusterID")) %>%
    full_join(density_wide_annual[[3]], by = c("GMU", "ClusterID")) %>%
    arrange(GMU, ClusterID)
  head(density_wide_annual_20s_22s)
  
  #'  Z-transform local abundance estimates (per year b/c annual estimates are stand alone variables)
  localN_z <- density_wide_annual_20s_22s %>%
    mutate(ClusterID = as.character(ClusterID)) %>%
    mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE))) 
    
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(contains(c(".yr1", ".yr2", ".yr3")))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z) # why is percent_forested and WSI so correlated???
  
  #'  Drop sites with NAs (missing 1+ years of data)
  localN_z <- drop_na(localN_z)
  
  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dplyr::select(dat, contains(c(".yr1", ".yr2", ".yr3"))) %>%
      as.matrix(.)
    for(i in 1:ncol(ndat)) {
      hist(ndat[,i])
    }
  }
  plot_histograms(localN_z)

  #'  Check out a few basic relationships
  summary(lm(whitetailed_deer.yr2 ~ mountain_lion.yr1, data = localN_z)) 
  summary(lm(whitetailed_deer.yr3 ~ mountain_lion.yr2, data = localN_z)) 
  summary(lm(elk.yr2 ~ mountain_lion.yr1, data = localN_z))  
  summary(lm(elk.yr3 ~ mountain_lion.yr2, data = localN_z))  
  summary(lm(whitetailed_deer.yr2 ~ wolf.yr1, data = localN_z))
  summary(lm(whitetailed_deer.yr3 ~ wolf.yr2, data = localN_z))
  summary(lm(elk.yr2 ~ wolf.yr1, data = localN_z))
  summary(lm(elk.yr3 ~ wolf.yr2, data = localN_z))
  summary(lm(moose.yr2 ~ wolf.yr1, data = localN_z))
  summary(lm(moose.yr3 ~ wolf.yr2, data = localN_z))
  summary(lm(whitetailed_deer.yr2 ~ coyote.yr1, data = localN_z))
  summary(lm(whitetailed_deer.yr3 ~ coyote.yr2, data = localN_z))
  summary(lm(mountain_lion.yr2 ~ wolf.yr1, data = localN_z))
  summary(lm(mountain_lion.yr3 ~ wolf.yr2, data = localN_z))
  summary(lm(bear_black.yr2 ~ wolf.yr1, data = localN_z))
  summary(lm(bear_black.yr3 ~ wolf.yr2, data = localN_z))
  summary(lm(coyote.yr2 ~ wolf.yr1, data = localN_z))
  summary(lm(coyote.yr3 ~ wolf.yr2, data = localN_z))
  summary(lm(coyote.yr2 ~ mountain_lion.yr1, data = localN_z))  
  summary(lm(coyote.yr3 ~ mountain_lion.yr2, data = localN_z))  
  summary(lm(moose.yr2 ~ DecFeb_WSI.yr1, data = localN_z))
  summary(lm(moose.yr3 ~ DecFeb_WSI.yr2, data = localN_z))
  summary(lm(bear_black.yr2 ~ percent_forest.yr1, data = localN_z))  
  summary(lm(bear_black.yr3 ~ percent_forest.yr2, data = localN_z))  
  
  plot(whitetailed_deer.yr3 ~ mountain_lion.yr2, data = localN_z)
  plot(moose.yr3 ~ wolf.yr2, data = localN_z)
  plot(whitetailed_deer.yr2 ~ coyote.yr1, data = localN_z)
  plot(whitetailed_deer.yr3 ~ coyote.yr2, data = localN_z)
  
  #'  ----------------------------------------------
  #####  Time period format: t-1 vs t across years  #####
  #'  ----------------------------------------------
  #'  Grouping data across years so t-1 data affects t data, regardless of year,
  #'  increases sample size and eliminates annual variation in estimated relationships
  #'  but prohibits estimation of indirect effects (yr 1 & 2 affect yr 2 & 3... yr2 on both sides)
  #'  Stack year 1 & year 2 data
  dat_t_minus_1 <- dat_stack_list[[1]] %>% #full_bx_dat_stack[[1]] %>%
    #'  Indicate whether observation was from first or second year in grouped data
    mutate(GroupYear = ifelse(year == "yr1", "first", "second")) %>%
    #'  Add time period identifier to each column name
    rename_with(.cols = percent_forest:wolf, function(x){paste0(x, ".Tminus1")}) %>%
    dplyr::select(-year)
  #'  Stack year 2 & year 3 data
  dat_t <- dat_stack_list[[2]] %>% #full_bx_dat_stack[[2]] %>%
    #'  Indicate whether observation was from first or second year in grouped data
    mutate(GroupYear = ifelse(year == "yr2", "first", "second")) %>%
    #'  Add time period identifier to each column name
    rename_with(.cols = percent_forest:wolf, function(x){paste0(x, ".T")}) %>%
    dplyr::select(-year)
  
  #'  Join t-1 and t data based on camera location
  density_wide_1YrLag_20s_22s <- full_join(dat_t_minus_1, dat_t, by = c("GMU", "ClusterID", "GroupYear")) %>%
    relocate(GroupYear, .after = ClusterID) %>%
    #'  Drop sites with NAs (missing 1+ years of data)
    na.omit(.)
  
  #'  Z-transform local abundance estimates (per year b/c annual estimates are stand alone variables)
  localN_z_1YrLag <- density_wide_1YrLag_20s_22s %>%
    mutate(ClusterID = as.character(ClusterID)) %>%
    mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    #mutate(across(starts_with(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(contains(c(".Tminus1", ".T")))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z_1YrLag) #'  Annual prey N correlated across years
  
  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dplyr::select(dat, contains(c(".Tminus1", ".T"))) %>% as.matrix(.)
    for(i in 1:ncol(ndat)) {
      hist(ndat[,i])
    }
  }
  plot_histograms(localN_z_1YrLag)
  
  #'  Check out a few basic relationships
  summary(lm(whitetailed_deer.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + coyote.Tminus1, data = localN_z_1YrLag))  
  summary(lm(elk.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag))  
  summary(lm(moose.T ~ wolf.Tminus1, data = localN_z_1YrLag))
  summary(lm(bear_black.T ~ wolf.Tminus1 + mountain_lion.T, data = localN_z_1YrLag))
  summary(lm(mountain_lion.T ~ wolf.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag))
  summary(lm(coyote.T ~ wolf.Tminus1 + mountain_lion.Tminus1, data = localN_z_1YrLag))
  
  plot(whitetailed_deer.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag)
  plot(elk.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag)
  plot(whitetailed_deer.T ~ wolf.Tminus1, data = localN_z_1YrLag)
  plot(moose.T ~ wolf.Tminus1, data = localN_z_1YrLag)
  plot(coyote.T ~ wolf.Tminus1, data = localN_z_1YrLag)
  
  #'  -------------------------------------------------------
  #####  Fully combined format: All years stacked together  #####
  #'  -------------------------------------------------------
  #'  Z-transform local abundance estimates (across years)
  localN_z_all <- density_wide_allyrs %>%
    mutate(ClusterID = as.character(ClusterID)) %>%
    mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    #mutate(across(starts_with(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat[,4:12]
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z_all) 

  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dat[,4:12] %>% as.matrix(.)
    for(i in 1:ncol(ndat)) {
      hist(ndat[,i])
    }
  }
  plot_histograms(localN_z_all)
  
  #'  Check out a few basic relationships
  summary(lm(whitetailed_deer ~ mountain_lion + wolf + coyote, data = localN_z_all))  
  summary(lm(elk ~ mountain_lion + wolf + bear_black, data = localN_z_all))  
  summary(lm(moose ~ wolf, data = localN_z_all))
  summary(lm(bear_black ~ wolf + mountain_lion, data = localN_z_all))
  summary(lm(mountain_lion ~ wolf + bear_black, data = localN_z_all))
  summary(lm(coyote ~ wolf + mountain_lion + bear_black, data = localN_z_all))
  
  plot(whitetailed_deer ~ mountain_lion, data = localN_z_all)
  plot(elk ~ mountain_lion, data = localN_z_all)
  plot(whitetailed_deer ~ wolf, data = localN_z_all)
  plot(moose ~ wolf, data = localN_z_all)
  plot(coyote ~ wolf, data = localN_z_all)
  
  #' #'  Save outputs
  #' save(localN_z, file = "./Data/Relative abundance data/RAI Phase 2/data_for_SEM_annual_n.RData")
  #' save(localN_z_1YrLag, file = "./Data/Relative abundance data/RAI Phase 2/data_for_SEM_1YrLag_n.RData")
  #' save(localN_z_all, file = "./Data/Relative abundance data/RAI Phase 2/data_for_SEM_allYrs_n.RData")
  
  
  print("IGNORE WARNINGS!")
  
  
  
  
