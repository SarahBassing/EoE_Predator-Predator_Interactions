  #'  -------------------------
  #'  Figures for publication
  #'  ICFWRU
  #'  Bassing
  #'  November 2024
  #'  -------------------------
  #'  Script to generate figures used in publication.
  #'  -------------------------
  
  #'  Load libraries
  library(sf)
  library(terra)
  library(ggplot2)
  library(tidyverse)
  
  #'  Load spatial data
  #'  Cluster polygons
  gmu1_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu1.shp") %>% mutate(GMU = "GMU1")
  gmu6_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu6.shp") %>% mutate(GMU = "GMU6")
  gmu10a_poly <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_cluster_polygons_gmu10a.shp") %>% mutate(GMU = "GMU10A")
  
  #'  Camera clusters
  clusters_gmu1 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu1.shp") 
  clusters_gmu6 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu6.shp") 
  clusters_gmu10a <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu10a.shp") 
  
  #'  Relative Density Index per cluster (loads as cluster_density)
  load("./Outputs/Relative_Abundance/RN_model/RelativeDensityIndex_per_SppCluster.RData")
  head(cluster_density[[1]]); head(cluster_density[[2]]); head(cluster_density[[3]])
  
  #'  Camera locations
  cams_20s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_20s_wgs84.shp")
  cams_21s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_21s_wgs84.shp")
  cams_22s_wgs84 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/cams_22s_wgs84.shp")
  cam_list <- list(cams_20s_wgs84, cams_21s_wgs84, cams_22s_wgs84)
  
  #'  Idaho GMUs and water
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
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
  

  #'  -----------------------
  ####  Format spatial data  ####
  #'  -----------------------
  #'  Merge cluster polygons across GMUs
  cluster_poly <- bind_rows(gmu1_poly, gmu6_poly, gmu10a_poly) %>%
    rename(ClusterID = Clusters)
  
  #'  Merge spatial camera cluster data together and reduce to one observation per camera
  clusters_all <- bind_rows(clusters_gmu1, clusters_gmu6, clusters_gmu10a) %>%
    dplyr::select(c(NwLctID, CellID, GMU, Setup, Clustrs, are_km2, geometry)) %>%
    rename(NewLocationID = NwLctID) %>%
    rename(ClusterID = Clustrs) %>%
    rename(area_km2 = are_km2) %>%
    group_by(NewLocationID) %>%
    slice(1L) %>%
    ungroup()
    
  #'  Merge RDI into single document
  rdi <- bind_rows(cluster_density) %>%
    filter(!is.na(ClusterID))
    
  #'  -------------------------------------------
  ####  Plot Relative Density Index per Cluster  ####
  #'  -------------------------------------------
  map_rdi <- function(sf_rdi, q) {
    
    #'  Define 5%, 25%, 50%, 75%, and 95% quantiles of RDI
    q0 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0), 0))
    q5 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.05), 0))
    q25  <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.25), 0))
    q50 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.5), 0))
    q75 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.75), 0))
    q95 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.95), 0))
    q100 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 1), 0))
    
    #'  Define breaks in continuous fill scale (if using this version)
    size_breaks <- c(q0, q5, q25, q50, q75, q95, q100)
    print(size_breaks)
    
    #'  Define spatial extent of study areas
    bbox <- st_bbox(eoe_gmu_wgs84)
    
    #'  Create figure
    spp_rdi <- ggplot() +
      # geom_sf(data = gmu_wgs84, fill = NA, col = "gray70") +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(data = sf_rdi, aes(fill = SppDensity.100km2.r, group = Year)) +  #, alpha = 8/10
      scale_fill_gradient2(low = "#f5f5f5", mid = "#54a1a1", high = "#1f6f6f", midpoint = q) +
      # scale_fill_gradient(low = "#f5f5f5", high = "#1f6f6f", breaks = size_breaks) +
      geom_sf(data = rivers_clip, color = "lightskyblue3") +
      geom_sf(data = bigwater_wgs84, fill = "lightskyblue2") +
      coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
      labs(fill = "Relative \nDensity \nIndex (RDI)", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 18)) +
      facet_wrap(~Year)
    
    return(spp_rdi)
  }
  (rdi_map_bear_2020 <- map_rdi(rdi[rdi$Species == "bear_black" & rdi$Year == "2020",], q = 9) + theme(axis.title.x = element_blank())) 
  (rdi_map_bear_2021 <- map_rdi(rdi[rdi$Species == "bear_black" & rdi$Year == "2021",], q = 8) + 
      guides(fill="none") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_bear_2022 <- map_rdi(rdi[rdi$Species == "bear_black" & rdi$Year == "2022",], q = 9) + 
      guides(fill="none") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_bear <- rdi_map_bear_2020 + rdi_map_bear_2021 + rdi_map_bear_2022 + plot_layout(ncol = 3, guides = 'collect') + 
    plot_annotation(title = "Relative density indices for black bears", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_coy_2020 <- map_rdi(rdi[rdi$Species == "coyote" & rdi$Year == "2020",], q = 11) + guides(fill="none") + theme(axis.title.x = element_blank()))
  (rdi_map_coy_2021 <- map_rdi(rdi[rdi$Species == "coyote" & rdi$Year == "2021",], q = 11) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_coy_2022 <- map_rdi(rdi[rdi$Species == "coyote" & rdi$Year == "2022",], q = 11) + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_coy <- rdi_map_coy_2020 + rdi_map_coy_2021 + rdi_map_coy_2022 + plot_layout(ncol = 3, guides = 'collect') + 
    plot_annotation(title = "Relative density indices for coyotes", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_lion_2020 <- map_rdi(rdi[rdi$Species == "mountain_lion" & rdi$Year == "2020",], q = 7) + theme(axis.title.x = element_blank()))
  (rdi_map_lion_2021 <- map_rdi(rdi[rdi$Species == "mountain_lion" & rdi$Year == "2021",], q = 7) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_lion_2022 <- map_rdi(rdi[rdi$Species == "mountain_lion" & rdi$Year == "2022",], q = 7) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_lion <- rdi_map_lion_2020 + rdi_map_lion_2021 + rdi_map_lion_2022 + plot_layout(ncol = 3, guides = 'collect') + 
    plot_annotation(title = "Relative density indices for mountain lions", theme = theme(plot.title = element_text(size = 18))))
  
  # (rdi_map_wolf_2020 <- map_rdi(rdi[rdi$Species == "wolf" & rdi$Year == "2020",], q = 2) + theme(legend.position="bottom", axis.title.x = element_blank()))
  # (rdi_map_wolf_2021 <- map_rdi(rdi[rdi$Species == "wolf" & rdi$Year == "2021",], q = 4) + theme(legend.position="bottom") + 
  #     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  # (rdi_map_wolf_2022 <- map_rdi(rdi[rdi$Species == "wolf" & rdi$Year == "2022",], q = 1) + theme(legend.position="bottom") + 
  #     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  # (rdi_map_wolf <- rdi_map_wolf_2020 + rdi_map_wolf_2021 + rdi_map_wolf_2022 + plot_layout(ncol = 3) +  #, guides = 'collect'
  #     plot_annotation(title = "Relative density indices for wolves", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_elk_2020 <- map_rdi(rdi[rdi$Species == "elk" & rdi$Year == "2020",], q = 16) + theme(axis.title.x = element_blank()))
  (rdi_map_elk_2021 <- map_rdi(rdi[rdi$Species == "elk" & rdi$Year == "2021",], q = 16) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_elk_2022 <- map_rdi(rdi[rdi$Species == "elk" & rdi$Year == "2022",], q = 16) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_elk <- rdi_map_elk_2020 + rdi_map_elk_2021 + rdi_map_elk_2022 + plot_layout(ncol = 3, guides = 'collect') + 
      plot_annotation(title = "Relative density indices for elk", theme = theme(plot.title = element_text(size = 18))))
  
  
  (rdi_map_moose_2020 <- map_rdi(rdi[rdi$Species == "moose" & rdi$Year == "2020",], q = 7) + guides(fill="none") + 
      theme(axis.title.x = element_blank()))
  (rdi_map_moose_2021 <- map_rdi(rdi[rdi$Species == "moose" & rdi$Year == "2021",], q = 7) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_moose_2022 <- map_rdi(rdi[rdi$Species == "moose" & rdi$Year == "2022",], q = 7) + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_moose <- rdi_map_moose_2020 + rdi_map_moose_2021 + rdi_map_moose_2022 + plot_layout(ncol = 3, guides = 'collect') + 
      plot_annotation(title = "Relative density indices for moose", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_wtd_2020 <- map_rdi(rdi[rdi$Species == "whitetailed_deer" & rdi$Year == "2020",], q = 20) + guides(fill="none") +
      theme(axis.title.x = element_blank()))
  (rdi_map_wtd_2021 <- map_rdi(rdi[rdi$Species == "whitetailed_deer" & rdi$Year == "2021",], q = 20) + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_wtd_2022 <- map_rdi(rdi[rdi$Species == "whitetailed_deer" & rdi$Year == "2022",], q = 20) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_wtd <- rdi_map_wtd_2020 + rdi_map_wtd_2021 + rdi_map_wtd_2022 + plot_layout(ncol = 3, guides = 'collect') + 
      plot_annotation(title = "Relative density indices for white-tailed deer", theme = theme(plot.title = element_text(size = 18))))
  
  #'  Plot wolf RDI across years to help with color gradient
  #'  Scale of RDI is so different in 2021 that it makes things hard to compare if
  #'  each year is plotted separately
  wolf_rdi <- function(sf_rdi, q) {
    
    #'  Define 5%, 25%, 50%, 75%, and 95% quantiles of RDI
    q0 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0), 0))
    q5 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.05), 0))
    q25  <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.25), 0))
    q50 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.5), 0))
    q75 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.75), 0))
    q95 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.95), 0))
    q100 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 1), 0))
    
    #'  Define breaks in continuous fill scale (if using this version)
    size_breaks <- c(q0, q5, q25, q50, q75, q95, q100)
    print(size_breaks)
    
    #'  Define spatial extent of study areas
    bbox <- st_bbox(eoe_gmu_wgs84)
    
    #'  Create figure
    spp_rdi <- ggplot() +
      geom_sf(data = eoe_gmu_wgs84, fill = NA) +
      geom_sf(data = sf_rdi, aes(fill = SppDensity.100km2.r, group = Year)) + #, alpha = 3/10
      scale_fill_gradient2(low = "#f5f5f5", mid = "#54a1a1", high = "#1f6f6f", midpoint = q) +
      # scale_fill_gradient(low = "#f5f5f5", high = "#1f6f6f", breaks = size_breaks) +
      geom_sf(data = rivers_clip, color = "lightskyblue3") +
      geom_sf(data = bigwater_wgs84, fill = "lightskyblue2") +
      # coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
      labs(fill = "Relative \nDensity \nIndex (RDI)", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 18)) +
      ggtitle("Relative density indices for wolves") + 
      facet_wrap(~Year)
    
    return(spp_rdi)
  }
  (rdi_map_wolf <- wolf_rdi(rdi[rdi$Species == "wolf",], q = 4) + 
      plot_annotation(theme = theme(plot.title = element_text(size = 18))))
  
  
  #'  Save
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_bear.tiff", rdi_map_bear,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_coyote.tiff", rdi_map_coy,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_lion.tiff", rdi_map_lion,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_wolf.tiff", rdi_map_wolf,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_elk.tiff", rdi_map_elk,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_moose.tiff", rdi_map_moose,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/Relative_Abundance/RN_model/Figures/RDI_wtd.tiff", rdi_map_wtd,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  
  
  