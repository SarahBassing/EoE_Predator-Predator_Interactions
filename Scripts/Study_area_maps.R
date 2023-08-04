  #'  ------------------------------
  #'  Study area maps
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  August 2023
  #'  -------------------------------
  #'  Script to create figures for study area maps
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(sf)
  library(terra)
  library(ggplot2)
  library(ggspatial)
  library(grid)
  library(cowplot)
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
  gmu_wgs84 <- st_transform(gmu, "+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Define projections and units
  gmu_proj <- crs(eoe_gmu) #crs("+init=EPSG:2243") #NAD83 / Idaho West (ftUS)
  id_proj <- crs(id)
  wgs84 <- crs(eoe_gmu_wgs84) 
  
  #'  Make camera location data spatial sf objects
  spatial_locs <- function(locs, proj) {
    locs <- arrange(locs, NewLocationID)
    sf_locs <- st_as_sf(locs, coords = c("Long", "Lat"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      dplyr::select(c("NewLocationID")) 
    return(sf_locs)
  }
  cams_wgs84 <- lapply(cams_list, spatial_locs, proj = wgs84)
  
  #'  Create text to overlap inset map
  grob <- grobTree(textGrob("Idaho, USA", x = 0.40,  y = 0.85, hjust = 0,
                            gp = gpar(col = "black", fontsize = 14, fontface = "italic")))
  
  #'  Plot state of Idaho for inset
  ID_state <- ggplot() + 
    geom_sf(data = gmu_wgs84, fill = "gray95", color = "gray50") +
    geom_sf(data = eoe_gmu_wgs84, color = "gray20", aes(fill = NAME), alpha = 0.5) +
    scale_fill_manual(values = c("#364B9A", "#FDB366", "#A50026")) +
    guides(fill = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
          #'  No margins around figure
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
    annotation_custom(grob)
  plot(ID_state)
  
  #'  Plot study areas against DEM with camera locations
  cams2021_wgs84 <- cams_wgs84[[2]]
  cam_SA_map <- ggplot() +
    # geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) + 
    #' #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    #' scale_alpha(range = c(0.3, 0.8)) +
    #' #'  Change colors of the raster
    #' scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area outlines and label with their names
    geom_sf(data = gmu_wgs84, fill = "gray95", color="gray50", size = 0.80) +
    geom_sf(data = eoe_gmu_wgs84, aes(fill = NAME), color="black", size = 0.80, alpha = 0.5) +
    scale_fill_manual(values = c("#364B9A", "#FDB366", "#A50026")) +
    geom_sf(data = cams2021_wgs84, color = "black", shape = 20, size = 2) +
    #'  Constrain plot to two study areas plus some room on the side & bottom
    coord_sf(xlim = c(-117.2, -113.5), ylim = c(46, 49), expand = TRUE) +
    labs(fill = "GMU") +
    #'  Get rid of lines and axis names
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
          axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 16), panel.border = element_blank(),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16)) +
      theme(legend.position = c(0.90, 0.3)) +
    labs(x = "Longitude", y = "Latitude") +
    #'  Add north arrow
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    #'  Add scale bar (be sure to double check the scale)
    annotation_scale(location = "br", width_hint = 0.5)
  plot(cam_SA_map)
  
  
  #'  Build plot with map of study areas and inset map of WA
  #'  https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  #'  Requires "cowplot" package
  #'  Don't use png or other calls to save while plotting- formatting gets messed up
  #'  Use export option in Plot window and formatting holds
  # tiff(file = "./Outputs/Figures/Maps/StudyAreas_Cameras1820_05.31.22.tiff",
  #     width = 1000, height = 691) 
  tiff(file = "./Outputs/Figures/StudyAreaMap_Cameras.tiff",
       units="in", width = 7, height = 9, res=800, compression = 'lzw') 
  StudyArea_Map <- ggdraw(cam_SA_map) +
    draw_plot(
      {
        ID_state #+
      },
      #'  Distance along a (0,1) x-axis to draw the left edge of the plot
      x = 0.6,
      #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
      y = 0.55,
      #'  Width & height of the plot expressed as proportion of the entire ggdraw object
      #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
      width = 0.35,
      height = 0.35)
  plot(StudyArea_Map)
  dev.off()
  