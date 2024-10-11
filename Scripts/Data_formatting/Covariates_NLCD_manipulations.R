  #'  ---------------------------------------------------
  #'  Calculate Percent Landcover Type & Distance to Edge
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  ---------------------------------------------------
  #'  Script to reclassify NLCD 2019 landcover types to a smaller set of categories, 
  #'  and create rasters representing the percent of lancover type within 500m of
  #'  each pixel using a moving window analysis.
  #'  Reclassification and moving window code provided by L. Satterfield. Script
  #'  also reclassifies rasters to 1 and NA for habitat of interest and everything
  #'  else, which is used to calculate distance to edge in ArcGIS.
  #'
  #'  Initial NLCD 2019 Classifications:
  #'  11 Open water; 12 Perennial ice/snow; 21 Developed - open space; 
  #'  22 Developed - low intensity; 23 Developed - med intensity; 24 Developed - 
  #'  high intensity; 31 Barren land (rock/sand/clay); 41 Deciduous forest;
  #'  42 Evergreen forest; 43 Mixed forest; 52 Shrub/scrub; 71 Grassland/herbaceous;
  #'  81 Pasture/hay; 82 Cultivated crops; 90 Woody wetland; 95 Emergent herbaceous wetland
  #'  ---------------------------------------------------
  
  #'  Load libraries
  library(sf)
  # library(stars)
  # library(rgeos)
  # library(raster)
  library(terra)
  library(tidyverse)
  
  #' #'  Read in original landcover raster and Idaho state shapefile
  #' nlcd19 <- rast("./Shapefiles/National Land Cover Database (NCLD)/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img") 
  #' crs(nlcd19, describe = TRUE, proj = TRUE)
  #' res(nlcd19)
  #' nlcd_proj <- crs(nlcd19)
  #' 
  #' #'  Crop NLCD raster to Idaho State extent, plot, and save
  #' nlcd_id <- crop(nlcd19, id)
  #' plot(nlcd_id)
  #' plot(id[1], add = T)
  #' writeRaster(nlcd_id, filename = "./Shapefiles/National Land Cover Database (NCLD)/NLCD19_Idaho.tif", overwrite = TRUE)
  
  #'  Read in cropped landcover raster
  nlcd <- rast("./Shapefiles/National Land Cover Database (NCLD)/NLCD19_Idaho.tif")
  nlcd_proj <- crs(nlcd)
  
  id <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp") %>%
    st_transform(nlcd_proj)
  projection(id)
  
  #'  Generate terrain ruggedness index data from DEM
  elev <- rast("./Shapefiles/IDFG spatial data/Elevation__10m2.tif")
  tri <- terrain(elev, v = "TRI", filename = "./Shapefiles/IDFG spatial data/TRI__10m2.tif")
  
  
  #'  -----------------------------------  
  ####  Reclassify landcover categories  ####
  #'  -----------------------------------
  #'  E.g., Forest landcover types are labeled 41-43, so everything from = 40  
  #'  (non-inclusive), to = 43, (including 43) gets assigned value of 43. 
  #'  Anything from 43 (non-inclusive) through 95 gets a value of 0.
  
  #'  "Forest": 41 Deciduous forest; 42 Evergreen forest; 43 Mixed forest
  forest <- matrix(c(0,40,0,
                     40,43,1,
                     43,95,0), ncol = 3, byrow = TRUE)
  #'  Not forest: everything else
  notforest <- matrix(c(0,40,1,
                        40,43,0,
                        43,95,1), ncol = 3, byrow = TRUE)
  
  #'  Check to make sure it looks right
  forest
  notforest
  #'  Reclassify the raster based on a matrix
  forest19 <- classify(nlcd, forest, include.lowest = FALSE)
  notforest19 <- classify(nlcd, notforest, include.lowest = FALSE)
  
  #'  Plot to see how it looks - only forest areas (binary format)
  plot(forest19, main = "Forested landcover")
  plot(notforest19, main = "Non-forested landcover")
  writeRaster(forest19, filename = "./Shapefiles/National Land Cover Database (NCLD)/Forested_landcover.tif")
  writeRaster(notforest19, filename = "./Shapefiles/National Land Cover Database (NCLD)/NonForested_landcover.tif")
  
  #'  --------------------------
  ####  Moving window analysis  ####
  #'  --------------------------
  #'  Buffer is based on the resolution/projection of the input raster
  #'  If in UTMs then 500m, if in lat/long then 0.005 degrees is approx. 500m
  # buffer <- raster::focalWeight(forest19, 0.005, "circle") 
  buffer_500m <- raster::focalWeight(forest19, 500, "circle") # when projected
  buffer_100m <- raster::focalWeight(forest19, 100, "circle")
  
  #'  Create proportional cover rasters: 
  #'  Multiples the binary raster by the focal weight and then sums within the buffer
  perc_forest_500m <- focal(forest19, buffer_500m)
  perc_forest_100m <- focal(forest19, buffer_100m)
  
  #'  Check it out
  plot(perc_forest_500m)
  plot(perc_forest_100m)
  
  #'  Save
  writeRaster(perc_forest_500m, filename = "./Shapefiles/National Land Cover Database (NCLD)/PercentForest_500m.tif", overwrite = TRUE)
  writeRaster(perc_forest_100m, filename = "./Shapefiles/National Land Cover Database (NCLD)/PercentForest_100m.tif", overwrite = TRUE)
  
  
  
  
  
  
  