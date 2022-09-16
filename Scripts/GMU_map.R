library(sf)
library(sp)
library(raster)
library(ggplot2)
library(tidyverse)

gmu <- st_read("C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp")
eoe_gmu <- gmu[gmu$NAME == "1" | gmu$NAME == "6" | gmu$NAME == "10A",]
usa <- st_read("C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
id <- usa[usa$NAME == "Idaho",]
gmu_proj <- projection(gmu)
projection(id)
id_reproj <- st_transform(id, gmu_proj) #"epsg:2243"
projection(id_reproj)
# st_write(eoe_gmu, dsn = "C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp", layer = "EoE_GMUs.shp", driver = "ESRI Shapefile")
# st_write(id_reproj, dsn = "C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/tl_2012_us_state/IdahoState.shp", layer = "IdahoState_2243.shp", driver = "ESRI Shapefile")

eoe_gmus <- st_read("C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp")
id <- st_read("C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/tl_2012_us_state/IdahoState.shp")
bbox <- st_bbox(eoe_gmus)

gmu_map <- ggplot() +
  geom_sf(data = gmu) +
  geom_sf(data = eoe_gmus, aes(fill = NAME)) +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  guides(fill=guide_legend(title="GMU")) +
  # coord_sf(xlim = c(-13050000, -12700000), ylim = c(5700000, 6274865), expand = TRUE) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.85))
ggsave(filename = "C:/Users/sarah/Desktop/Coop Unit Work/gmu_map.png", plot = gmu_map, 
       width = 9, height = 5, dpi = 300)    


            