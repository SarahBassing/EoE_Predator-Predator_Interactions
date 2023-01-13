  #'  ------------------------------
  #'  Global Human Settlement Layer
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2023
  #'  ------------------------------
  #'  Script to convert Global Human Settlement Layer categories to binary 
  #'  rasters based on level of human settlement: high density clusters [3], 
  #'  low density clusters [2], rural grid cells [1], and inhabited areas [0]. 
  #'  This can be used to generate a proximity raster in ArcGIS.
  #'  
  #'  GHSL cropped to Idaho spatial extent and downloaded from Google Earth Engine
  #'  Data details: https://developers.google.com/earth-engine/datasets/catalog/JRC_GHSL_P2016_SMOD_POP_GLOBE_V1#description
  #'  ------------------------------
  
  #'  Load libraries
  library(sf)
  library(terra)
  
  #'  Load data
  id <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp")
  ihsl <- rast("./Shapefiles/GEE/HumanSettlement/Idaho_HSL_epsg8826.tif")
  
  #'  Reproject Idaho shapefile  
  (nad83 <- crs(elev, proj = TRUE)) # Idaho Transverse Mercator NAD83
  id_nad83 <- st_transform(id, nad83)
  
  plot(ihsl)
  plot(id_nad83, add = T)
  
  #'  Reformat raster so high density clusters are 1 and other values become NA
  hdc <- ihsl
  hdc[hdc == 0] <- NA
  hdc[hdc == 1] <- NA
  hdc[hdc == 2] <- NA
  hdc[hdc == 3] <- 1
  plot(hdc, col = "red")
  
  #'  Reformat raster so high/low density clusters are 1, rural/inhabited become NA
  urban <- ihsl
  urban[urban == 0] <- NA
  urban[urban == 1] <- NA
  urban[urban == 2] <- 1
  urban[urban == 3] <- 1
  plot(urban, col = "red")
  
  #'  Reformat raster so high/low density clusters & rural areas are 1, inhabited become NA
  developed <- ihsl
  developed[developed == 0] <- NA
  developed[developed == 1] <- 1
  developed[developed == 2] <- 1
  developed[developed == 3] <- 1
  plot(developed, col = "red")
  
  
  #'  Save!
  writeRaster(hdc, "./Shapefiles/GEE/HumanSettlement/IHSL_Urban.tif", overwrite = TRUE)
  writeRaster(urban, "./Shapefiles/GEE/HumanSettlement/IHSL_UrbanSuburban.tif", overwrite = TRUE)
  writeRaster(developed, "./Shapefiles/GEE/HumanSettlement/IHSL_UrbanRural.tif", overwrite = TRUE)
  
  