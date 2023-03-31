#  ID CRU - Predator Interactions

library(sf)
library(sp)
library(raster)
library(ggplot2)
library(tidyverse)

gmu <- st_read("./Shapefiles/IDFG_Game_Management_Units/Game_Management_Units.shp")
eoe_gmu <- gmu[gmu$NAME == "1" | gmu$NAME == "6" | gmu$NAME == "10A",]
usa <- st_read("./Shapefiles/tl_2012_us_state/tl_2012_us_state.shp")
id <- usa[usa$NAME == "Idaho",]
gmu_proj <- projection(gmu)
wgs84 <- crs("+proj=longlat +datum=WGS84 +no_defs")
projection(id)
id_reproj <- st_transform(id, gmu_proj) #"epsg:2243"
projection(id_reproj)
id_wgs84 <- st_transform(id, wgs84)
projection(id_wgs84)
extent(id_wgs84)
# st_write(eoe_gmu, dsn = "C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp", layer = "EoE_GMUs.shp", driver = "ESRI Shapefile")
# st_write(id_reproj, dsn = "C:/Users/sarah/Desktop/Coop Unit Work/Shapefiles/tl_2012_us_state/IdahoState.shp", layer = "IdahoState_2243.shp", driver = "ESRI Shapefile")

eoe_gmus <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp")
id <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp")
bbox <- st_bbox(eoe_gmus)

units::set_units(st_area(eoe_gmus), km^2) #' area of each GMU in sq-km (order = GMU1, GMU6, GMU10A)
#'  Centroid of each gmu (lat/long)
gmu_center_latlong <- eoe_gmus %>%
  st_centroid() %>%
  st_transform(wgs84) %>%
  st_geometry() %>%
  unlist()
gmu_center_latlong


gmu_map <- ggplot() +
  geom_sf(data = gmu) +
  geom_sf(data = eoe_gmus, aes(fill = NAME)) +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99")) +
  guides(fill=guide_legend(title="GMU")) +
  # coord_sf(xlim = c(-13050000, -12700000), ylim = c(5700000, 6274865), expand = TRUE) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.85))
# ggsave(filename = "./gmu_map.png", plot = gmu_map, 
#        width = 9, height = 5, dpi = 300)    


            