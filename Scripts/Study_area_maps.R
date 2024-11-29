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
  # library(elementalist)
  library(tidyterra)
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
  
  #' #'  Spatial data
  #' world.data <- map_data("world")
  #' world.usa <- filter(world.data, region == "USA")
  #' world.canada <- filter(world.data, region == "Canada")
  #' world.mexico <- filter(world.data, region == "Mexico")
  #' states.data <- map_data("state")
  #' states.idaho <- filter(states.data, region == "idaho")
  #' states.alaska <- filter(world.usa, subregion == "Alaska")
  
  gmu <- st_read("./Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp")
  eoe_gmu <- gmu[gmu$NAME == "1" | gmu$NAME == "6" | gmu$NAME == "10A",] 
  usa <- st_read("./Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
  id <- usa[usa$NAME == "Idaho",] 
  # na_wgs84 <- st_read("./Shapefiles/stanford-NorthAmerica-shapefile/ns372xw1938.shp")
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  id_wgs84 <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  gmu_wgs84 <- st_transform(gmu, "+proj=longlat +datum=WGS84 +no_defs")
  usa_wgs84 <- st_transform(usa, "+proj=longlat +datum=WGS84 +no_defs")
  
  dem_low <- rast("./Shapefiles/National Elevation Dataset (NED) NAD83/DEM_100m_res.tiff")
  dem_low_wgs84 <- project(dem_low, "+proj=longlat +datum=WGS84 +no_defs")
  #'  Crop low rez dem to extent of ID shapefile
  v <- vect(id_wgs84)
  dem_low_crop <- crop(dem_low_wgs84, v, mask = TRUE)
  
  #'  Define projections and units
  gmu_proj <- st_crs(eoe_gmu) #crs("+init=EPSG:2243") #NAD83 / Idaho West (ftUS)
  id_proj <- st_crs(id)
  wgs84 <- st_crs(eoe_gmu_wgs84) 
  
  #'  Identify bounding box of EoE study areas
  st_bbox(eoe_gmu_wgs84)
  
  #'  Make camera location data spatial sf objects
  spatial_locs <- function(locs, proj) {
    locs <- arrange(locs, NewLocationID)
    sf_locs <- st_as_sf(locs, coords = c("Long", "Lat"), crs = wgs84) %>%
      mutate(obs = seq(1:nrow(.))) %>%
      dplyr::select(c("NewLocationID")) 
    return(sf_locs)
  }
  cams_wgs84 <- lapply(cams_list, spatial_locs, proj = wgs84)
  
  #'  Order the GMUs
  eoe_gmu_wgs84 <- mutate(eoe_gmu_wgs84, NAME = factor(NAME, levels = c("1", "6", "10A")))
  
  #'  Study area colors
  eoe_color <- c("#364B9A", "#FDB366", "#A50026")
  
  #'  2021 camera locations
  cams2021_wgs84 <- cams_wgs84[[2]]
  
  #'  --------------------------------
  ####  Maps within maps within maps  ####
  #'  --------------------------------
  #####  Map #1  #####
  #'  Plot USA with Idaho highlighted
  usa_map <- ggplot() +
    geom_sf(data = usa_wgs84, fill = "white", color = "gray50", linewidth = 0.2) +
    geom_sf(data = id_wgs84, fill = "gray20", color = "gray20") +
    #'  Constrain plot to the lower 48 
    coord_sf(xlim = c(-124.2, -68.2), ylim = c(25, 49)) +
    theme_minimal() +  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
          #'  No margins around figure
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  #####  Map #2  #####
  #'  Text to label state of Idaho
  grob <- grobTree(textGrob("Idaho, USA", x = 0.45,  y = 0.95, hjust = 0,
                            gp = gpar(col = "black", fontsize = 10, fontface = "italic"))) #18
  
  #' #'  Order GMUs geographically (N to S)
  #' sa_gmu_wgs84$NAME <- factor(sa_gmu_wgs84$NAME, levels = c("4", "28", "33", "34", "35"))
  
  #'  Plot study areas within Idaho
  ID_study_areas <- ggplot() +
    geom_sf(data = gmu_wgs84, fill = "gray95", color="gray50", linewidth = 0.2) + #0.5
    geom_sf(data = eoe_gmu_wgs84, aes(fill = NAME, color = NAME), linewidth = 0.2, alpha = 0.7) + #0.75
    scale_fill_manual(values = eoe_color) +
    scale_color_manual(values = eoe_color) +
    labs(fill = "GMU", color = "GMU") +
    geom_rect(aes(xmin = -117.3, xmax = -114.15, ymin = 45.75, ymax = 49.25), color = "black", fill = NA, size = 0.2) + # 0.5
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(size = 6, family = "serif"), #18
          legend.text = element_text(size = 6, family = "serif"), #16
          legend.key.size = unit(0.5, "lines")) +
          # legend.margin=margin(0,0,-10,-10),
          # legend.box.margin=margin(0,0,0,0),
          # legend.margin=margin(t = -5, unit='cm'),
          # legend.spacing.x = unit(-1, "cm"))+
          # plot.margin = unit(c(0,-1,0,-1), "cm")) + 
    theme(legend.justification = c(1, 0)) 
  
  #####  Map #3  #####
  #'  Plot close up of study areas
  zoomed_in_gmus <- ggplot() +
    geom_sf(data = id_wgs84, fill = "gray85", color = "gray50", linewidth = 0.5) +
    geom_spatraster(data = dem_low_crop) + 
    scale_fill_continuous(low = "gray90", high = "gray25", na.value = "transparent") + 
    geom_sf(data = eoe_gmu_wgs84, color = "black", fill = c("#364B9A", "#FDB366", "#A50026"), alpha = 0.15, linewidth = 0.25, show.legend = FALSE) +  #0.75
    geom_sf(data = cams2021_wgs84, color = "black", shape = 20, size = 0.2) + #2
    #'  Constrain plot to study areas plus some room on the side & bottom
    coord_sf(xlim = c(-117.2, -114.25), ylim = c(45.75, 49), expand = TRUE) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.title.x = element_text(size = 8, family = "serif"), #18 
          axis.title.y = element_text(size = 8, family = "serif"), #18
          axis.text.x = element_text(size = 6, angle = 45, hjust = 1, colour = "black", family = "serif"),  #16
          axis.text.y = element_text(size = 6, colour = "black", family = "serif"),  #16
          axis.line = element_line(color = 'black', linewidth = 0.2),
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          legend.title = element_text(size = 8, family = "serif"), #18
          legend.text = element_text(size = 6, family = "serif"),
          legend.key.size = unit(0.75, "lines")) + #16
    labs(x = "Longitude", y = "Latitude", fill = "Elevation (m)") + 
    theme(legend.justification = c(1, 0),
          plot.margin = unit(c(0.1,0,0,0.1), "cm")) + 
    #'  Add north arrow
    annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.05, "cm"), pad_y = unit(0.4, "cm"), #0.2, "in"
                           height = unit(0.75, "cm"), width = unit(0.75, "cm"),
                           style = north_arrow_fancy_orienteering(fill = c("gray20", "white",
                                                                           line_col = "gray10"), 
                                                                  line_width = 0.5, text_size = 6, text_family = "serif")) +
    #'  Add scale bar (be sure to double check the scale)
    annotation_scale(location = "bl", width_hint = 0.5, bar_cols = c("gray10", "white"),
                     height = unit(0.1, "cm"), line_width = 0.5, text_cex = 0.6, text_family = "serif") #0.5
  
  #'  Build plot with map of study areas and inset map of WA
  #'  https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  #'  Requires "cowplot" package
  #'  Don't use png or other calls to save while plotting- formatting gets messed up
  #'  Use export option in Plot window and formatting holds
  # tiff(file = "./Outputs/Figures/StudyAreaMap_v1.tiff",
  #      units = "in", width = 7, height = 9, res = 800, compress = 'lzw')
  StudyArea_Map <- ggdraw(ID_study_areas + 
                            theme(legend.box.margin = margin(0,-8,0,-10))) + 
    draw_plot(
      {
        usa_map 
      },
      #'  Distance along a (0,1) x-axis to draw the left edge of the plot
      x = 0.40, 
      #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
      y = 0.50, 
      #'  Width & height of the plot expressed as proportion of the entire ggdraw object
      #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
      width = 0.55, 
      height = 0.55) + 
    theme(panel.background = element_rect(fill = "white", color = "black", linewidth = 0.15),
          plot.margin = unit(c(0.35, 0.5, 0.65, -0.75), "cm")) #c(0.5, 0.5, 1, -0.75), "cm")
  plot(StudyArea_Map)
  # dev.off()
  
  # tiff(file = "./Outputs/Figures/Map_StudyAreas_CameraSites.tiff",
  #      units = "in", width = 10, height = 11, res = 800, compress = 'lzw')
  tiff(file = "./Outputs/Figures/Map_StudyAreas_CameraSites.tiff",
       units = "cm", width = 8.5, height = 10, res = 600, compress = 'lzw')
  Full_Map <- ggdraw(zoomed_in_gmus + 
                       theme(legend.key.size = unit(0.55, "lines"),
                             legend.box.margin=margin(0,0,0,-5))) +
    draw_plot(
      {
        StudyArea_Map #+
      },
      #'  Distance along a (0,1) x-axis to draw the left edge of the plot
      x = 0.62,  #0.60
      #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
      y = 0.50, # 0.52
      #'  Width & height of the plot expressed as proportion of the entire ggdraw object
      #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
      width = 0.40, #0.37
      height = 0.52)  #0.45
  plot(Full_Map)
  dev.off()
  
  
  #' #'  Plot state of Idaho for inset
  #' ID_state <- ggplot() + 
  #'   geom_sf(data = na_wgs84, fill = "gray95") +
  #'   geom_sf(data = id_wgs84, fill = "gray40", color = "gray20") + #gmu_wgs84
  #'   geom_sf(data = eoe_gmu_wgs84, color = c("#364B9A", "#FDB366", "#A50026"), aes(fill = NAME)) + #, alpha = 0.5
  #'   scale_fill_manual(values = c("#364B9A", "#FDB366", "#A50026")) +
  #'   #'  Constrain plot to lower 48 plut a smidge of Canada 
  #'   coord_sf(xlim = c(-155, -107), ylim = c(24, 52)) +
  #'   #coord_sf(xlim = c(-123.2, -65.5), ylim = c(25, 52.5)) +
  #'   #'  Add rectangle to highlight study areas
  #'   geom_rect(aes(xmin = -117.5, xmax = -114, ymin = 45.8, ymax = 49.5), 
  #'             color = "black", fill = NA, linetype = "dashed") + #linetype = "dotted", 
  #'   labs(fill = "GMU") +
  #'   # guides(fill = "none") +
  #'   theme_bw() +
  #'   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #'         panel.background = element_rect(fill = "slategray2"),
  #'         axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
  #'         axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
  #'         legend.title = element_text(size = 16),
  #'         legend.text = element_text(size = 14),
  #'         #'  No margins around figure
  #'         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
  #'   theme(legend.position = c(0.925, 0.85), 
  #'         legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) + #c(1.1, 0.3)
  #'   #left line
  #'   geom_segment(aes(x = -127.2, xend = -117.5, y = 49.52, yend = 49.5),
  #'                linetype = "dashed", color = "black", size = 0.5) +
  #'   geom_segment(aes(x = -127.2, xend = -117.5, y = 24.51, yend = 45.8),
  #'                linetype = "dashed", color = "black", size = 0.5) 
  #' plot(ID_state)
  #' 
  #' #'  Plot study areas against DEM with camera locations
  #' cams2021_wgs84 <- cams_wgs84[[2]]
  #' cam_SA_map <- ggplot() +
  #'   #'  Add study area outlines and label with their names
  #'   geom_sf(data = na_wgs84, fill = "gray95", color="gray50") +
  #'   #'  Elevation gradient as background in Idaho
  #'   geom_spatraster(data = dem_low_crop) +
  #'   scale_fill_continuous(low = "gray90", high = "gray15", na.value = "transparent") +
  #'   geom_sf(data = id_wgs84, fill = "transparent", color="gray50") +
  #'   geom_sf(data = eoe_gmu_wgs84, fill = c("#364B9A", "#FDB366", "#A50026"), # note the flip in color order
  #'           size = 0.80, alpha = 0.15) +
  #'   # geom_sf(data = eoe_gmu_wgs84, aes(fill = NAME), size = 0.80, alpha = 0.5) +
  #'   # scale_fill_manual(values = c("#364B9A", "#FDB366", "#A50026")) +
  #'   geom_sf(data = cams2021_wgs84, color = "black", shape = 20, size = 2) +
  #'   #'  Constrain plot to two study areas plus some room on the side & bottom
  #'   coord_sf(xlim = c(-117.2, -113.5), ylim = c(46, 49), expand = TRUE) +
  #'   # labs(fill = "GMU") +
  #'   labs(fill = "Elevation (m)") +
  #'   # guides(fill = "none") +
  #'   #'  Get rid of lines and axis names
  #'   theme_bw() +
  #'   theme(axis.line = element_line(colour = "black"),
  #'         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  #'         axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
  #'         axis.text.x = element_text(size = 16),  #, angle = 45, hjust = 1
  #'         axis.text.y = element_text(size = 16), panel.border = element_blank(),
  #'         legend.title = element_text(size = 16),
  #'         legend.text = element_text(size = 14),
  #'         plot.background = element_rect(colour = "black", fill = "white", size = 0.8, 
  #'                                        linetype = "dashed")) +
  #'     theme(legend.position = c(0.80, 0.60),
  #'           legend.background = element_rect(fill = "gray95")) + #c(0.90, 0.3)
  #'   labs(x = "Longitude", y = "Latitude") +
  #'   # Sets lat/long breaks, expands to fill (e.g.)
  #'   scale_x_continuous(breaks = seq(-118, -113, by = 1), expand = c(0.02, 0.02)) +
  #'   scale_y_continuous(breaks = seq(46, 49, by = 1), expand = c(0.02, 0.02)) +
  #'   #'  Add north arrow
  #'   annotation_north_arrow(location = "tr", which_north = "true", 
  #'                          pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
  #'                          style = north_arrow_fancy_orienteering) +
  #'   #'  Add scale bar (be sure to double check the scale)
  #'   annotation_scale(location = "tr", width_hint = 0.5, height = unit(0.1, "in"), 
  #'                    pad_y = unit(0.2, "in"))
  #' plot(cam_SA_map)
  #' 
  #' 
  #' #'  Build plot with map of study areas and inset map of WA
  #' #'  https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  #' #'  Requires "cowplot" package
  #' #'  Don't use png or other calls to save while plotting- formatting gets messed up
  #' #'  Use export option in Plot window and formatting holds
  #' tiff(file = "./Outputs/Figures/StudyAreaMap_Cameras_v3.tiff",
  #'      units = "in", width = 10, height = 8, res = 800, compression = 'lzw') 
  #' StudyArea_Map <- ggdraw(ID_state) +
  #'   
  #'   draw_plot(
  #'     {
  #'       cam_SA_map
  #'     },
  #'     #'  Distance along a (0,1) x-axis to draw the left edge of the plot
  #'     x = -0.07,
  #'     #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
  #'     y = 0.10,
  #'     #'  Width & height of the plot expressed as proportion of the entire ggdraw object
  #'     #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
  #'     width = 0.75,
  #'     height = 0.75)
  #' plot(StudyArea_Map)
  #' dev.off()
  
  