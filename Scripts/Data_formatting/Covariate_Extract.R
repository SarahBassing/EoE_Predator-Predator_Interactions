  #'  --------------------------------------
  #'  Extract covariates at camera locations
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  --------------------------------------
  #'  Read in spatial data and extract covariate values at camera locations.
  #'  
  #'  Requires:
  #'  Problem camera data frames from Detection_data_cleaning.R
  #'  Wolf detection and minimum group count data from Detection_histories_for_occmod.R
  #'  Reported mortality data from Covariates_NLCD_manipulations.R
  #'  And spatial layers
  #'  --------------------------------------
  
  #'  Load libraries
  library(sf)
  library(rgeos)
  library(raster)
  library(terra)
  library(tidyverse)
  library(ggplot2)
  library(grid)
  library(MASS)
  
  #'  ----------------
  ####  Spatial data  ####
  #'  ----------------
  #'  Load and review
  pforest <- rast("./Shapefiles/National Land Cover Database (NCLD)/PercentForeest_500m.tif")
  nlcd <- rast("./Shapefiles/National Land Cover Database (NCLD)/NLCD19_Idaho.tif")
  id <- st_read("./Shapefiles/tl_2012_us_state/IdahoState.shp")
  elev <- rast("./Shapefiles/IDFG spatial data/Elevation__10m2.tif")
  habclass <- rast("./Shapefiles/IDFG spatial data/HabLayer_30m2.tif")
  dist2suburbs <- rast("./Shapefiles/GEE/HumanSettlement/Dist2Suburbs.tif")
  dist2rural <- rast("./Shapefiles/GEE/HumanSettlement/Dist2Rural.tif")
  #'  IDFG Geodatabase with roads
  idfg_gdb <- "./Shapefiles/IDFG spatial data/IDFG Geodatabase.gdb"
  st_layers(idfg_gdb)
  rds <- sf::st_read(dsn = idfg_gdb, layer = "Road_OpenStreetMap")

  #'  Take a closer look  
  crs(pforest, describe = TRUE, proj = TRUE)
  crs(nlcd, describe = TRUE, proj = TRUE)
  crs(elev, describe = TRUE, proj = TRUE)
  crs(habclass, describe = TRUE, proj = TRUE)
  crs(dist2suburbs, describe = TRUE, proj = TRUE)
  # projection(id)
  projection(rds)
  
  res(pforest)
  res(nlcd)
  res(elev)
  res(habclass)
  res(dist2suburbs)
  
  #'  Define projections to use when reprojecting camera locations
  wgs84 <- crs("+proj=longlat +datum=WGS84 +no_defs")
  (aea <- crs(nlcd, proj = TRUE))
  (nad83 <- crs(elev, proj = TRUE)) # Idaho Transverse Mercator NAD83
  (hab_crs <- crs(habclass, proj = TRUE))
  
  #'  Reproject shapefiles
  id_aea <- st_transform(id, aea)
  id_nad83 <- st_transform(id, nad83)

  #'  Convert NAs to 0 - happens when pixel was suburban or rural and surrounded 
  #'  by other suburban or rural pixels so no distance to calculated (should be 0)
  dist2suburbs[is.na(dist2suburbs)] <- 0
  dist2rural[is.na(dist2rural)] <- 0
  
  #'  ---------------
  ####  Camera data  ####
  #'  ---------------
  #'  Load camera location data and format
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  #'  List camera locations
  cams_list <- list(eoe_probcams_20s, eoe_probcams_20w, eoe_probcams_21s)
  
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
  cams_aea <- lapply(cams_list, spatial_locs, proj = aea)
  cams_nad83 <- lapply(cams_list, spatial_locs, proj = nad83)
  cams_hab_crs <- lapply(cams_list, spatial_locs, proj = hab_crs)
  
  #'  Double check these are plotting correctly
  plot(pforest, main = "Camera locations over percent forested habitat, 500m radius")
  plot(id_aea[1], add = TRUE, col = NA)
  plot(cams_aea[[3]], add = TRUE, col = "black", cex = 0.75)

  
  #'  ------------------------
  ####  Additional data sets  ####
  #'  ------------------------
  #'  Wolf data generated from camera detections
  load("./Data/Wolf count data/count_eoe20s_wolf.RData")
  load("./Data/Wolf count data/min_group_size_eoe20s.RData") 
  
  load("./Data/Wolf count data/count_eoe20w_wolf.RData")
  load("./Data/Wolf count data/min_group_size_eoe20w.RData") 
  
  load("./Data/Wolf count data/count_eoe21s_wolf.RData")
  load("./Data/Wolf count data/min_group_size_eoe21s.RData")
  
  #'  Make sure everything is ordered correctly
  arrange_dat <- function(dat) {
    dat <- arrange(dat, NewLocationID)
    return(dat)
  }
  min_group_size_eoe20s <- arrange_dat(min_group_size_eoe20s)
  min_group_size_eoe20w <- arrange_dat(min_group_size_eoe20w)
  min_group_size_eoe21s <- arrange_dat(min_group_size_eoe21s)
  
  #'  Mortality data provided by IDFG
  load("./Data/IDFG BGMR data/mort_preSmr20.RData")
  load("./Data/IDFG BGMR data/mort_preWtr20.RData")
  load("./Data/IDFG BGMR data/mort_preSmr21.RData")
  
  reformat_mort_dat <- function(mort) {
    #'  Split predator mortality data by species
    mort <- dplyr::select(mort, -gmu_area)
    mort_bear <- filter(mort, Species == "Black Bear") %>%
      rename(Bear_mort_n = total_mortalities,
             Bear_mort_km2 = mortality_km2) %>%
      dplyr::select(-Species)
    mort_bob <- filter(mort, Species == "Bobcat") %>%
      rename(Bob_mort_n = total_mortalities,
             Bob_mort_km2 = mortality_km2) %>%
      dplyr::select(-c(Species, GMU))
    mort_lion <- filter(mort, Species == "Mountain Lion") %>%
      rename(Lion_mort_n = total_mortalities,
             Lion_mort_km2 = mortality_km2) %>%
      dplyr::select(-c(Species, GMU))
    mort_wolf <- filter(mort, Species == "Wolf") %>%
      rename(Wolf_mort_n = total_mortalities,
             Wolf_mort_km2 = mortality_km2) %>%
      dplyr::select(-c(Species, GMU))
    mort_df <- as.data.frame(cbind(mort_bear, mort_bob, mort_lion, mort_wolf))
    return(mort_df)
  }
  mort_Smr20_df <- reformat_mort_dat(mort_preSmr20) %>%
    filter(GMU != "GMU1")
  mort_Wtr20_df <- reformat_mort_dat(mort_preWtr20) %>%
    filter(GMU != "GMU1")
  mort_Smr21_df <- reformat_mort_dat(mort_preSmr21)
  
  #'  Relative abundance index (prey & human activity)
  #'  Using Hour of Detection as RA index b/c highly correlated with other 
  #'  definitions of independent detection events and consistent with Ausband et al. in review
  load("./Data/Relative abundance data/EoE_RelativeN_HrOfDetection.RData")
  reformat_relativeN_data <- function(dat) {
    RelativeN <- dat %>%
      #'  Create one column per species
      spread(Species, n_dets) %>%
      rowwise() %>%
      #'  Aggrigate relative abundance data into funcitonal groups
      mutate(human_plus = sum(human, cat_domestic, dog_domestic, horse, na.rm = TRUE),
             ungulate = sum(elk, moose, muledeer, whitetaileddeer, na.rm = TRUE),
             big_deer = sum(elk, moose, na.rm = TRUE),
             small_deer = sum(muledeer, whitetaileddeer, na.rm = TRUE)) %>%
      #'  Rename and drop columns
      rename(livestock = cattle_cow) %>%
      rename(lagomorphs = rabbit_hare) %>%
      dplyr::select(-c(cat_domestic, dog_domestic, horse)) %>%
      #'  Change all NAs introduced during spread to 0's (i.e., non-detection)
      replace(is.na(.), 0) %>%
      #'  Relocate columns so more intuitive order
      relocate(human_plus, .after = human) %>%
      relocate(lagomorphs, .before = moose) %>%
      relocate(livestock, .before = moose) %>%
      relocate(ungulate, .after = whitetaileddeer) %>%
      relocate(big_deer, .after = ungulate) %>%
      relocate(small_deer, .after = big_deer)
    return(RelativeN)
  }
  RA_Smr20_df <- reformat_relativeN_data(eoe_dethr_list[[1]])
  RA_Wtr20_df <- reformat_relativeN_data(eoe_dethr_list[[2]])
  RA_Smr21_df <- reformat_relativeN_data(eoe_dethr_list[[3]])
  
  
  #'  ----------------------------------
  ####  COVARIATE EXTRACTION & MERGING  ####
  #'  ----------------------------------
  cov_extract <- function(locs_aea, locs_nad83, locs_hab_crs, min_group_size, mort, relativeN) {
    
    #'  Extract covariate data for each camera site from spatial layers
    perc_forest <- terra::extract(pforest, vect(locs_aea)) %>%
      transmute(ID = ID, perc_forest = focal_sum) %>%
      mutate(perc_forest = round(perc_forest, 3))
    landcover <- terra::extract(nlcd, vect(locs_aea))
    elev <- terra::extract(elev, vect(locs_nad83))
    habitat <- terra::extract(habclass, vect(locs_hab_crs))
    dist2suburbs <- terra::extract(dist2suburbs, vect(locs_nad83))
    dist2rural <- terra::extract(dist2rural, vect(locs_nad83))
    #'  Find nearest road to each camera
    nearestrd <- st_nearest_feature(locs_nad83, rds)
    #'  Calculate distance to nearest road (in meters)
    dist2rd <- as.numeric(st_distance(locs_nad83, rds[nearestrd,], by_element = T))
    
    #'  Join each extracted covariate to the unique camera location data
    covs <- as.data.frame(locs_aea) %>%
      mutate(ID = seq(1:nrow(.))) %>%
      full_join(perc_forest, by = "ID") %>%
      full_join(landcover, by = "ID") %>%
      full_join(elev, by = "ID") %>%
      full_join(habitat, by = "ID") %>%
      full_join(dist2suburbs, by = "ID") %>%
      full_join(dist2rural, by = "ID") %>%
      cbind(dist2rd) %>%
      #'  Join with additional data sets already formatted for each camera site
      full_join(min_group_size, by = "NewLocationID") %>%
      mutate(GMU = sub("_.*", "", NewLocationID), 
             dist2rd = round(dist2rd, digits = 2)) %>%
      relocate(GMU, .after = NewLocationID) %>%
      full_join(relativeN, by = "NewLocationID") %>%
      #'  Change all NAs introduced during joining to 0's 
      #'  (MAKE SURE THIS ONLY AFFECTS RELATIVE ABUNDANCE DATA)
      replace(is.na(.), 0) %>%
      full_join(mort, by = "GMU") %>%
      dplyr::select(-c(geometry, ID, Lat, Long)) %>%
      arrange(NewLocationID)
    
     return(covs)
  }
  eoe_covs_20s <- cov_extract(locs_aea = cams_aea[[1]], locs_nad83 = cams_nad83[[1]], locs_hab_crs = cams_hab_crs[[1]], 
                              min_group_size = min_group_size_eoe20s, mort = mort_Smr20_df, relativeN = RA_Smr20_df)
  eoe_covs_20w <- cov_extract(locs_aea = cams_aea[[2]], locs_nad83 = cams_nad83[[2]], locs_hab_crs = cams_hab_crs[[2]], 
                              min_group_size = min_group_size_eoe20w, mort = mort_Wtr20_df, relativeN = RA_Wtr20_df)
  eoe_covs_21s <- cov_extract(locs_aea = cams_aea[[3]], locs_nad83 = cams_nad83[[3]], locs_hab_crs = cams_hab_crs[[3]], 
                              min_group_size = min_group_size_eoe21s, mort = mort_Smr21_df, relativeN = RA_Smr21_df)

  
  #'  -------------------------
  ####  Covariate exploration  ####
  #'  -------------------------
  #'  Histogram of covariate data
  spread_of_covariate_data <- function(covs, season) {
    hist(covs$Elevation__10m2, breaks =  20, main = paste("Frequency of elevation at cameras\n", season))
    hist(covs$perc_forest, breaks =  20, main = paste("Frequency of percent forest at cameras\n", season))
    hist(covs$Dist2Suburbs, breaks =  20, main = paste("Frequency of distance of camera to suburbs\n", season))
    hist(log(covs$dist2rd), breaks =  20, main = paste("Frequency of log distance of camera to nearest road\n", season))
    hist(covs$elk, breaks =  20, main = paste("Frequency of elk activity at cameras\n", season))
    hist(covs$human, breaks =  20, main = paste("Frequency of human activity at cameras\n", season))
    hist(covs$human_motorized, breaks =  30, main = paste("Frequency of motorized vehicles at cameras\n", season))
    hist(covs$lagomorphs, breaks =  20, main = paste("Frequency of lagomorph activity at cameras\n", season))
    hist(covs$livestock, breaks =  30, main = paste("Frequency of livestock activity at cameras\n", season))
    hist(covs$moose, breaks =  20, main = paste("Frequency of moose activity at cameras\n", season))
    hist(covs$muledeer, breaks =  20, main = paste("Frequency of mule deer activity at cameras\n", season))
    hist(covs$whitetaileddeer, breaks =  20, main = paste("Frequency of white-tailed deer activity at cameras\n", season))
    hist(covs$ungulate, breaks =  20, main = paste("Frequency of ungulate activity at cameras\n", season))
    hist(covs$big_deer, breaks = 20, main = paste("Frequency of elk & moose activity at cameras\n", season))
    hist(covs$small_deer, breaks = 20, main = paste("Frequency of mule deer & white-tail activity at cameras\n", season))
  }
  spread_of_covariate_data(eoe_covs_20s, season = "Summer 2020")
  spread_of_covariate_data(eoe_covs_20w, season = "Winter 2020-2021")
  spread_of_covariate_data(eoe_covs_21s, season = "Summer 2021")
  
  
  # ####  Is there a difference in relative abund by camera setup?  ####  
  # bunnies <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #   dplyr::select(c(NewLocationID, lagomorphs)) %>%
  #   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #          setup = factor(setup, levels = c("U", "P"))) %>%
  #   dplyr::select(-NewLocationID)
  # table(bunnies)
  # maxbunnies <- grobTree(textGrob("max count = 71 (U)", x=0.60,  y=0.95, hjust=0,
  #                             gp=gpar(col="red", fontsize=13, fontface="italic")))
  # ggplot(bunnies, aes(x = lagomorphs, fill = setup)) +
  #   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
  #   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #   labs(fill="") +
  #   ggtitle("Relative abundance of lagomorphs at Ungulate vs Predator cameras") +
  #   annotation_custom(maxbunnies)
  # 
  # ppl <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #   dplyr::select(c(NewLocationID, human)) %>%
  #   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #          setup = factor(setup, levels = c("U", "P"))) %>%
  #   dplyr::select(-NewLocationID)
  # table(ppl)
  # maxppl <- grobTree(textGrob("max count = 136 (P)", x=0.60,  y=0.95, hjust=0,
  #                           gp=gpar(col="red", fontsize=13, fontface="italic")))
  # ggplot(ppl, aes(x = human, fill = setup)) +
  #   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 3) +
  #   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #   labs(fill="") +
  #   ggtitle("Relative abundance of humans at Ungulate vs Predator cameras") +
  #   annotation_custom(maxppl)
  # 
  # elk <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #   dplyr::select(c(NewLocationID, elk)) %>%
  #   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #          setup = factor(setup, levels = c("U", "P"))) %>%
  #   dplyr::select(-NewLocationID)
  # table(elk)
  # maxelk <- grobTree(textGrob("max count = 145 (P)", x=0.60,  y=0.95, hjust=0,
  #                             gp=gpar(col="red", fontsize=13, fontface="italic")))
  # ggplot(elk, aes(x = elk, fill = setup)) +
  #   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 3) +
  #   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #   labs(fill="") +
  #   ggtitle("Relative abundance of elk at Ungulate vs Predator cameras") +
  #   annotation_custom(maxelk)
  # 
  # md <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #   dplyr::select(c(NewLocationID, muledeer)) %>%
  #   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #          setup = factor(setup, levels = c("U", "P"))) %>%
  #   dplyr::select(-NewLocationID)
  # table(md)
  # maxmd <- grobTree(textGrob("max count = 64 (U)", x=0.60,  y=0.95, hjust=0,
  #                             gp=gpar(col="red", fontsize=13, fontface="italic")))
  # ggplot(md, aes(x = muledeer, fill = setup)) +
  #   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
  #   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #   labs(fill="") +
  #   ggtitle("Relative abundance of mule deer at Ungulate vs Predator cameras") +
  #   annotation_custom(maxmd)
  # 
  # wtd <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #   dplyr::select(c(NewLocationID, whitetaileddeer)) %>%
  #   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #          setup = factor(setup, levels = c("U", "P"))) %>%
  #   dplyr::select(-NewLocationID)
  # table(wtd)
  # maxwtd <- grobTree(textGrob("max count = 225 (P)", x=0.60,  y=0.95, hjust=0,
  #                            gp=gpar(col="red", fontsize=13, fontface="italic")))
  # ggplot(wtd, aes(x = whitetaileddeer, fill = setup)) +
  #   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 3) +
  #   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #   labs(fill="") +
  #   ggtitle("Relative abundance of white-tails at Ungulate vs Predator cameras") +
  #   annotation_custom(maxwtd)
  # 
  # moose <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #   dplyr::select(c(NewLocationID, moose)) %>%
  #   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #          setup = factor(setup, levels = c("U", "P"))) %>%
  #   dplyr::select(-NewLocationID)
  # table(moose)
  # maxmoose <- grobTree(textGrob("max count = 44 (U)", x=0.60,  y=0.95, hjust=0,
  #                             gp=gpar(col="red", fontsize=13, fontface="italic")))
  # ggplot(moose, aes(x = moose, fill = setup)) +
  #   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
  #   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #   labs(fill="") +
  #   ggtitle("Relative abundance of moose at Ungulate vs Predator cameras") +
  #   annotation_custom(maxmoose)
  # 
  # livestock <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #   dplyr::select(c(NewLocationID, livestock)) %>%
  #   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #          setup = factor(setup, levels = c("U", "P"))) %>%
  #   dplyr::select(-NewLocationID)
  # table(livestock)
  # maxlivestock <- grobTree(textGrob("max count = 284 (U/P)", x=0.60,  y=0.95, hjust=0,
  #                                   gp=gpar(col="red", fontsize=13, fontface="italic")))
  # ggplot(livestock, aes(x = livestock, fill = setup)) +
  #   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 10) +
  #   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #   labs(fill="") +
  #   ggtitle("Relative abundance of livestock at Ungulate vs Predator cameras") +
  #   annotation_custom(maxlivestock)
  
  #'  -------------------------------------
  ####  Transform relative abundance data  ####
  #'  -------------------------------------
  #'  Function to test the Box-Cox transformation, square-root transformation,
  #'  and log transformation on each relative abundance data set
  transform_dat <- function(dat) {
    #'  Add 1 to all values so no zeros in data set
    dat1 <- dat + 1
    hist(dat + 1)
    
    #'  Find optimal lambda (tuning parameter) for Box-Cox transformation
    bc <- boxcox(dat1 ~ 1)
    print(lambda <- bc$x[which.max(bc$y)])

    #'  Transform data using Box-Cox transformation and optimal lambda
    bxcxdat <- round((dat1^lambda-1)/lambda, 2)

    #'  Visualize
    hist(bxcxdat)
    mod_bctrans <- lm(((dat1^lambda-1)/lambda) ~ 1)
    #'  QQplot of residuals (are they normal or at least more normal?)
    qqnorm(mod_bctrans$residuals, main = "boxcox data")
    qqline(mod_bctrans$residuals)
    
    #'  Inverse hyperbolic sine transformation
    ihs_dat <- lm(log(dat + sqrt(dat^2 + 1)) ~ 1)
    hist(log(dat + sqrt(dat^2 + 1)))
    qqnorm(ihs_dat$residuals, main = "inverse hyperbolic sine data")
    qqline(ihs_dat$residuals)
    ihsdat <- round(log(dat + sqrt(dat^2 + 1)), 2)

    #'  Is square rooting the data any better?
    sqrt_dat <- lm(sqrt(dat) ~ 1)
    hist(sqrt(dat))
    qqnorm(sqrt_dat$residuals, main = "sqare root data")
    qqline(sqrt_dat$residuals)
    sqrtdat <- round(sqrt(dat), 2)

    #'  What about logging it?
    log_dat <- lm(log(dat1) ~ 1)
    hist(log(dat + 1))
    qqnorm(log_dat$residuals, main = "log data + 1")
    qqline(log_dat$residuals)
    logdat1 <- round(log(dat1), 2)

    # return(bxcxdat)
    # return(ihsdat)
    return(sqrtdat)   #  Currently square-root transforming them all
    # return(logdat1)
  }
  covs <- eoe_covs_20s
  # covs <- eoe_covs_21s
  covs$lagomorphs_adj <- transform_dat(covs$lagomorphs) # 2020 all are very terrible, 2021 IHS or log better but definitely not normal
  covs$human_adj <- transform_dat(covs$human) # 20/21 IHS or log seem better but definitely not normal
  covs$human_motorized_adj <- transform_dat(covs$human_motorized) # 20/21 all are pretty terrible
  covs$elk_adj <- transform_dat(covs$elk) # 20/21 sqrt seems best but major tails
  covs$moose_adj <- transform_dat(covs$moose) # 20/21 sqrt seems best but definitely not normal
  covs$muledeer_adj <- transform_dat(covs$muledeer) # 2020 all are very terrible, 2021 all are slightly better but definitely  not normal
  covs$whitetailedeer_adj <- transform_dat(covs$whitetaileddeer) # 20/21 sqrt seems best but major tails
  covs$big_deer_adj <- transform_dat(covs$big_deer) # 20/21 sqrt seems best but major tails
  covs$small_deer_adj <- transform_dat(covs$small_deer) # 2020 sqrt seems best but major tails, Box-Cox seems best but still tails
  
  #'  List heavily skewed covaraites (relative abundance data)
  listnames <- c("lagomorphsT", "humanT", "human_motorizedT", "livestockT", "elkT", "mooseT", "muledeerT", "whitetaileddeerT", "ungulateT", "big_deerT", "small_deerT")
  cov_20s_list <- list(eoe_covs_20s$lagomorphs, eoe_covs_20s$human, eoe_covs_20s$human_motorized, eoe_covs_20s$livestock, eoe_covs_20s$elk, eoe_covs_20s$moose, eoe_covs_20s$muledeer, eoe_covs_20s$whitetaileddeer, eoe_covs_20s$ungulate, eoe_covs_20s$big_deer, eoe_covs_20s$small_deer)
  cov_20w_list <- list(eoe_covs_20w$lagomorphs, eoe_covs_20w$human, eoe_covs_20w$human_motorized, eoe_covs_20w$livestock, eoe_covs_20w$elk, eoe_covs_20w$moose, eoe_covs_20w$muledeer, eoe_covs_20w$whitetaileddeer, eoe_covs_20w$ungulate, eoe_covs_20w$big_deer, eoe_covs_20w$small_deer)
  cov_21s_list <- list(eoe_covs_21s$lagomorphs, eoe_covs_21s$human, eoe_covs_21s$human_motorized, eoe_covs_21s$livestock, eoe_covs_21s$elk, eoe_covs_21s$moose, eoe_covs_21s$muledeer, eoe_covs_21s$whitetaileddeer, eoe_covs_21s$ungulate, eoe_covs_21s$big_deer, eoe_covs_21s$small_deer)
  
  #'  Transform annual relative abundance data
  transf_covs_20s <- lapply(cov_20s_list, transform_dat)
  names(transf_covs_20s) <- listnames
  transf_covs_20s <- as.data.frame(transf_covs_20s)
  
  transf_covs_20w <- lapply(cov_20w_list, transform_dat)
  names(transf_covs_20w) <- listnames
  transf_covs_20w <- as.data.frame(transf_covs_20w)
  
  transf_covs_21s <- lapply(cov_21s_list, transform_dat)
  names(transf_covs_21s) <- listnames
  transf_covs_21s <- as.data.frame(transf_covs_21s)
  
  #'  Append adjusted relative abundance covariates to larger covariate data frame
  eoe_covs_20s <- cbind(eoe_covs_20s, transf_covs_20s)
  eoe_covs_20w <- cbind(eoe_covs_20w, transf_covs_20w)
  eoe_covs_21s <- cbind(eoe_covs_21s, transf_covs_21s)
  
  
  #'  -----------------------------------------------------------------
  ####  Identify outliers and cap covariates at 99th percentile value  ####
  #'  -----------------------------------------------------------------
  outliers <- function(cov, title) {
    #'  Summarize predicted values
    hist(cov, breaks = 100, main = title)
    boxplot(cov, main = title)
    #'  What value represents the 99th percentile in the covariate values
    quant <- quantile(cov, c(0.99), na.rm = TRUE)
    #'  Print that value and maximum value
    print(quant); print(max(cov, na.rm = TRUE))
    #'  Identify the 1% most extreme values and set to 99th percentile value
    extremevalues <- as.data.frame(cov) %>%
      mutate(outlier = ifelse(cov > quant, "outlier", "not_outlier"),
             adjusted_cov = ifelse(outlier == "outlier", quant, cov),
             adjusted_cov = round(adjusted_cov, 0))
    #'  How many covariate values are above the 99th percentile?
    outlier <- extremevalues[extremevalues$outlier == "outlier",]
    outlier <- filter(outlier, !is.na(outlier))
    print(nrow(outlier))
    #'  What percentage of observations are being forced to 99th percentile value
    print(nrow(outlier)/nrow(extremevalues))

    #'  Return these adjusted values
    adjusted_cov <- extremevalues$adjusted_cov
    return(adjusted_cov)
  }
  #' #'  Identify covariate outliers
  #' bunnies$lagomorphs_adj <- outliers(bunnies$lagomorphs, "Lagomorph outliers")
  #' ppl$humans_adj <- outliers(ppl$human, "Human outliers")
  #' elk$elk_adj <- outliers(elk$elk, "Elk outliers")
  #' moose$moose_adj <- outliers(moose$moose, "Moose outliers")
  #' md$muledeer_adj <- outliers(md$muledeer, "Mule deer outliers")
  #' wtd$whitetaileddeer_adj <- outliers(wtd$whitetaileddeer, "White-tailed deer outliers")

  #'  List heavily skewed covaraites (relative abundance data)
  listnames <- c("lagomorphs99", "human99", "human_motorized99", "livestock99", "elk99", "moose99", "muledeer99", "whitetaileddeer99", "ungulate99", "big_deer99", "small_deer99")
  cov_20s_list <- list(eoe_covs_20s$lagomorphs, eoe_covs_20s$human, eoe_covs_20s$human_motorized, eoe_covs_20s$livestock, eoe_covs_20s$elk, eoe_covs_20s$moose, eoe_covs_20s$muledeer, eoe_covs_20s$whitetaileddeer, eoe_covs_20s$ungulate, eoe_covs_20s$big_deer, eoe_covs_20s$small_deer)
  cov_20w_list <- list(eoe_covs_20w$lagomorphs, eoe_covs_20w$human, eoe_covs_20w$human_motorized, eoe_covs_20w$livestock, eoe_covs_20w$elk, eoe_covs_20w$moose, eoe_covs_20w$muledeer, eoe_covs_20w$whitetaileddeer, eoe_covs_20w$ungulate, eoe_covs_20w$big_deer, eoe_covs_20w$small_deer)
  cov_21s_list <- list(eoe_covs_21s$lagomorphs, eoe_covs_21s$human, eoe_covs_21s$human_motorized, eoe_covs_21s$livestock, eoe_covs_21s$elk, eoe_covs_21s$moose, eoe_covs_21s$muledeer, eoe_covs_21s$whitetaileddeer, eoe_covs_21s$ungulate, eoe_covs_21s$big_deer, eoe_covs_21s$small_deer)

  #'  Adjust annual relative abundance based on season-specific covariate ranges
  capped_covs_20s <- lapply(cov_20s_list, outliers, title = "Outliers")
  names(capped_covs_20s) <- listnames
  capped_covs_20s <- as.data.frame(capped_covs_20s)

  capped_covs_20w <- lapply(cov_20w_list, outliers, title = "Outliers")
  names(capped_covs_20w) <- listnames
  capped_covs_20w <- as.data.frame(capped_covs_20w)

  capped_covs_21s <- lapply(cov_21s_list, outliers, title = "Outliers")
  names(capped_covs_21s) <- listnames
  capped_covs_21s <- as.data.frame(capped_covs_21s)

  #'  Append adjusted relative abundance covariates to larger covariate data frame
  eoe_covs_20s <- cbind(eoe_covs_20s, capped_covs_20s)
  eoe_covs_20w <- cbind(eoe_covs_20w, capped_covs_20w)
  eoe_covs_21s <- cbind(eoe_covs_21s, capped_covs_21s)


  #' #'  Visualize capped relative abundance data based on camera setup
  #' bunnies <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #'   dplyr::select(c(NewLocationID, lagomorphs99)) %>%
  #'   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #'          setup = factor(setup, levels = c("U", "P"))) %>%
  #'   dplyr::select(-NewLocationID)
  #' table(bunnies)
  #' maxbunnies <- grobTree(textGrob("max count = 41 (U/P)", x=0.60,  y=0.95, hjust=0,
  #'                                 gp=gpar(col="red", fontsize=13, fontface="italic")))
  #' ggplot(bunnies, aes(x = lagomorphs99, fill = setup)) +
  #'   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
  #'   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #'   labs(fill="") +
  #'   ggtitle("Relative abundance of lagomorphs at Ungulate vs Predator cameras") +
  #'   annotation_custom(maxbunnies)
  #' 
  #' ppl <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #'   dplyr::select(c(NewLocationID, human99)) %>%
  #'   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #'          setup = factor(setup, levels = c("U", "P"))) %>%
  #'   dplyr::select(-NewLocationID)
  #' table(ppl)
  #' maxppl <- grobTree(textGrob("max count = 35 (P)", x=0.60,  y=0.95, hjust=0,
  #'                             gp=gpar(col="red", fontsize=13, fontface="italic")))
  #' ggplot(ppl, aes(x = human99, fill = setup)) +
  #'   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
  #'   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #'   labs(fill="") +
  #'   ggtitle("Relative abundance of humans at Ungulate vs Predator cameras") +
  #'   annotation_custom(maxppl)
  #' 
  #' elk <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #'   dplyr::select(c(NewLocationID, elk99)) %>%
  #'   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #'          setup = factor(setup, levels = c("U", "P"))) %>%
  #'   dplyr::select(-NewLocationID)
  #' table(elk)
  #' maxelk <- grobTree(textGrob("max count = 55 (U/P)", x=0.60,  y=0.95, hjust=0,
  #'                             gp=gpar(col="red", fontsize=13, fontface="italic")))
  #' ggplot(elk, aes(x = elk99, fill = setup)) +
  #'   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 1) +
  #'   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #'   labs(fill="") +
  #'   ggtitle("Relative abundance of elk at Ungulate vs Predator cameras") +
  #'   annotation_custom(maxelk)
  #' 
  #' livestock <- rbind(eoe_covs_20s, eoe_covs_21s) %>%
  #'   dplyr::select(c(NewLocationID, livestock99)) %>%
  #'   mutate(setup = ifelse(grepl("P", NewLocationID), "P", "U"),
  #'          setup = factor(setup, levels = c("U", "P"))) %>%
  #'   dplyr::select(-NewLocationID)
  #' table(livestock)
  #' maxlivestock <- grobTree(textGrob("max count = 91 (U/P)", x=0.60,  y=0.95, hjust=0,
  #'                             gp=gpar(col="red", fontsize=13, fontface="italic")))
  #' ggplot(livestock, aes(x = livestock99, fill = setup)) +
  #'   geom_histogram(color = "#e9ecef", alpha = 0.6, position = "identity", binwidth = 2) +
  #'   scale_fill_manual(values=c("#69b3a2", "#404080")) +
  #'   labs(fill="") +
  #'   ggtitle("Relative abundance of livestock at Ungulate vs Predator cameras") +
  #'   annotation_custom(maxlivestock)
  
 
  ####  Save  ###
  write.csv(eoe_covs_20s, file = "./Data/Covariates_extracted/Covariates_EoE_Smr20.csv")
  save(eoe_covs_20s, file = "./Data/Covariates_extracted/Covariates_EoE_Smr20.RData")
  
  write.csv(eoe_covs_20w, file = "./Data/Covariates_extracted/Covariates_EoE_Wtr20.csv")
  save(eoe_covs_20w, file = "./Data/Covariates_extracted/Covariates_EoE_Wtr20.RData")
  
  write.csv(eoe_covs_21s, file = "./Data/Covariates_extracted/Covariates_EoE_Smr21.csv")
  save(eoe_covs_21s, file = "./Data/Covariates_extracted/Covariates_EoE_Smr21.RData")
  