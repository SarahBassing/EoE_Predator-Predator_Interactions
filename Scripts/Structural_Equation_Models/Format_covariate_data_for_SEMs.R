  #'  --------------------------------
  #'  Prep Covariate Data for SEMs
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  October 2024
  #'  --------------------------------
  #'  Format covaraite data for SEMs. This includes:
  #'    1) Importing National Land Cover Dataset
  #'  and generating percent forest cover covariates for each cluster based on 
  #'  averaging all pixels within a cluster vs averaging pixel values within
  #'  a buffered area around each camera location.
  #'    2) Importing PRISM data (extracted from Google Earth Engine)
  #'  and combinging monthly minimum temp and monthly total precip to generate a
  #'  Winter Severity Index across each cluster
  #'    3) Importing Global Forest Change data (extracted from Google Earth Engine)
  #'  and calculating the accumulated canopy loss over time, the percentage of
  #'  forested habitat that's been disturbed over the last 20 years, and the 
  #'  percentage of forested habitat in each cluster
  #'    4) Importing IDFG BGMR harvest data (prepared by K. Peterson)
  #'  and calculating harvest density per cluster (number of wolves reported
  #'  harvested/km^2 from the previous year)
  #'  
  #'  Note the varying time-lags with covariate data: 
  #'    1) NLCD data generated 2019 so represents land cover classifications 1 - 3
  #'  years prior to annual RAI estimates.
  #'    2) Winter Severity Index represents winter conditions 6-months prior to 
  #'  annual RAI estimates, based on the hypothesis that winter severity of the 
  #'  previous winter affects probability of survival and thus abundance in summer.
  #'    3) Percentage of disturbed habitat and percentage of forested habitat
  #'  represent forest habitat conditions 1-year prior to annual RAI estimates,
  #'  based on the hypothesis that forage quality/availability the previous summer
  #'  affects body condition entering winter, which affects probability of 
  #'  survival and birth rates and thus abundance in summer.
  #'    4) Harvest intensity represents wolf harvest for the 12 months prior to
  #'  annual RAI estimates (e.g., summer RAI = June 1 - Sept. 15 so harvest
  #'  intensity includes any harvest that occurred from previous June 1 to 
  #'  current June 1), based on the hypothesis that harvest during the previous
  #'  season's rendezvousing (previous summer), dispersing (fall), breeding (winter), 
  #'  or denning (spring) periods affect abundance in summer.
  #'  --------------------------------
  
  library(sf)
  library(terra)
  library(tidyverse)
  
  
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
  
  #'  ---------------------------------------
  ####  NLCD Percent Forest Cover Covariate  ####
  #'  ---------------------------------------
  #'  Read in forested landcover
  perc_forest <- rast("./Shapefiles/National Land Cover Database (NCLD)/PercentForest_100m.tif")
  forest <- rast("./Shapefiles/National Land Cover Database (NCLD)/Forested_landcover.tif")
  notforest <- rast("./Shapefiles/National Land Cover Database (NCLD)/NonForested_landcover.tif")
  forest_proj <- crs(forest)
  
  #'  Calculate average percent forested habitat across cameras in each cluster
  extract_forest <- function(cams, gmu) {
    #'  Transform and reduce to a single observation per camera (don't need to
    #'  repeat information across years)
    cams <- st_transform(cams, forest_proj) %>%
      group_by(NewLocationID) %>%
      slice(1L) %>%
      ungroup()
    #'  Extract % forested habitat within 100m radius of each camera site and
    #'  average across cameras per cluster
    perc_forest_100m <- terra::extract(perc_forest, vect(cams)) %>%
      dplyr::select(-ID) %>%
      bind_cols(cams) %>%
      rename("percent_forest" = "focal_sum") %>%
      group_by(ClusterID) %>%
      summarise(percent_forest = mean(percent_forest)) %>%
      ungroup() %>%
      mutate(GMU = gmu)
    return(perc_forest_100m)
  }
  forest_gmu1 <- extract_forest(clusters_gmu1, gmu = "GMU1")
  forest_gmu6 <- extract_forest(clusters_gmu6, gmu = "GMU6")
  forest_gmu10a <- extract_forest(clusters_gmu10a, gmu = "GMU10A")

  percforest_cluster <- bind_rows(forest_gmu1, forest_gmu6, forest_gmu10a) %>%
    dplyr::select(c(GMU, ClusterID, percent_forest))
  
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
  
  
  #'  ---------------------------------
  ####  Google Earth Engine data sets  ####
  #'  ---------------------------------
  #'  -----------------------
  #####  PRISM weather data  #####
  #'  -----------------------
  #'  Source script to load & format average monthly precip & temp data
  #'  Produces standardized average monthly total precipitation and standardized 
  #'  average monthly minimum temp for each winter and GMU for past 50 years
  # source("./Scripts/Structural_Equation_Models/Format_weather_data.R")   
  
  #'  Load 50 years of monthly PRISM weather data, reformat and generate 
  #'  standardized monthly averages per winter from Dec. 1972 - Feb. 2023.
  #'  Standardizing across 50 years of data means each monthly average is 
  #'  represented relative to the 50-year average and 50-year variability.
  #'  Winter months = December, January, February
  
  #'  PRISM total monthly precipitation (averaged per GMU) per cluster
  precip <- read_csv("./Data/GEE outputs/PRISM_ClusterAvg_monthly_total_precip_1972_2023.csv") %>%  
    dplyr::select(c(featureID, Clusters, date, meanMonthlyValue)) %>%
    mutate(GMU = ifelse(str_detect(featureID, "1_1_"), "GMU1", featureID),
           GMU = ifelse(str_detect(GMU, "1_2_"), "GMU10A", GMU),
           GMU = ifelse(str_detect(GMU, "2_"), "GMU6", GMU))
  #'  PRISM minimum monthly temperature (averaged per GMU) per cluster
  temp <- read_csv("./Data/GEE outputs/PRISM_ClusterAvg_monthly_min_temp_1972_2023.csv") %>%  
    dplyr::select(c(featureID, Clusters, date, meanMonthlyValue)) %>%
    mutate(GMU = ifelse(str_detect(featureID, "1_1_"), "GMU1", featureID),
           GMU = ifelse(str_detect(GMU, "1_2_"), "GMU10A", GMU),
           GMU = ifelse(str_detect(GMU, "2_"), "GMU6", GMU))
  
  #'  Function to add season column to weather data
  add_season <- function(prism) {
    annual_winter_season <- prism %>%
      #'  Reformat
      mutate(Year = as.numeric(format(as.Date(date, format="%Y-%m-%d"),"%Y")),
             Month = as.numeric(format(as.Date(date, format="%Y-%m-%d"),"%m"))) %>%
      #'  Filter to winter months (Dec, Jan, Feb) of each year
      filter(Month <= 2 | Month == 12) %>%
      #'  Label monthly data with a yearly seasonal classification (e.g., Dec. '72 - Feb. '73 = Wtr7273)
      #'  Is there a more efficient way to do this? Probably.
      #'  Am I going to spend the time to figure it out? Nope.
      mutate(Season = ifelse(date <= "1973-02-01", "Wtr7273", "Wtr2223"),
             Season = ifelse(date >= "1973-12-01" & date <= "1974-02-01", "Wtr7374", Season),
             Season = ifelse(date >= "1974-12-01" & date <= "1975-02-01", "Wtr7475", Season),
             Season = ifelse(date >= "1975-12-01" & date <= "1976-02-01", "Wtr7576", Season),
             Season = ifelse(date >= "1976-12-01" & date <= "1977-02-01", "Wtr7677", Season),
             Season = ifelse(date >= "1977-12-01" & date <= "1978-02-01", "Wtr7778", Season),
             Season = ifelse(date >= "1978-12-01" & date <= "1979-02-01", "Wtr7879", Season),
             Season = ifelse(date >= "1979-12-01" & date <= "1980-02-01", "Wtr7980", Season),
             Season = ifelse(date >= "1980-12-01" & date <= "1981-02-01", "Wtr8081", Season),
             Season = ifelse(date >= "1981-12-01" & date <= "1982-02-01", "Wtr8182", Season),
             Season = ifelse(date >= "1982-12-01" & date <= "1983-02-01", "Wtr8283", Season),
             Season = ifelse(date >= "1983-12-01" & date <= "1984-02-01", "Wtr8384", Season),
             Season = ifelse(date >= "1984-12-01" & date <= "1985-02-01", "Wtr8485", Season),
             Season = ifelse(date >= "1985-12-01" & date <= "1986-02-01", "Wtr8586", Season),
             Season = ifelse(date >= "1986-12-01" & date <= "1987-02-01", "Wtr8687", Season),
             Season = ifelse(date >= "1987-12-01" & date <= "1988-02-01", "Wtr8788", Season),
             Season = ifelse(date >= "1988-12-01" & date <= "1989-02-01", "Wtr8889", Season),
             Season = ifelse(date >= "1989-12-01" & date <= "1990-02-01", "Wtr8990", Season),
             Season = ifelse(date >= "1990-12-01" & date <= "1991-02-01", "Wtr9091", Season),
             Season = ifelse(date >= "1991-12-01" & date <= "1992-02-01", "Wtr9192", Season),
             Season = ifelse(date >= "1992-12-01" & date <= "1993-02-01", "Wtr9293", Season),
             Season = ifelse(date >= "1993-12-01" & date <= "1994-02-01", "Wtr9394", Season),
             Season = ifelse(date >= "1994-12-01" & date <= "1995-02-01", "Wtr9495", Season),
             Season = ifelse(date >= "1995-12-01" & date <= "1996-02-01", "Wtr9596", Season),
             Season = ifelse(date >= "1996-12-01" & date <= "1997-02-01", "Wtr9697", Season),
             Season = ifelse(date >= "1997-12-01" & date <= "1998-02-01", "Wtr9798", Season),
             Season = ifelse(date >= "1998-12-01" & date <= "1999-02-01", "Wtr9899", Season),
             Season = ifelse(date >= "1999-12-01" & date <= "2000-02-01", "Wtr9900", Season),
             Season = ifelse(date >= "2000-12-01" & date <= "2001-02-01", "Wtr0001", Season),
             Season = ifelse(date >= "2001-12-01" & date <= "2002-02-01", "Wtr0102", Season),
             Season = ifelse(date >= "2002-12-01" & date <= "2003-02-01", "Wtr0203", Season),
             Season = ifelse(date >= "2003-12-01" & date <= "2004-02-01", "Wtr0304", Season),
             Season = ifelse(date >= "2004-12-01" & date <= "2005-02-01", "Wtr0405", Season),
             Season = ifelse(date >= "2005-12-01" & date <= "2006-02-01", "Wtr0506", Season),
             Season = ifelse(date >= "2006-12-01" & date <= "2007-02-01", "Wtr0607", Season),
             Season = ifelse(date >= "2007-12-01" & date <= "2008-02-01", "Wtr0708", Season),
             Season = ifelse(date >= "2008-12-01" & date <= "2009-02-01", "Wtr0809", Season),
             Season = ifelse(date >= "2009-12-01" & date <= "2010-02-01", "Wtr0910", Season),
             Season = ifelse(date >= "2010-12-01" & date <= "2011-02-01", "Wtr1011", Season),
             Season = ifelse(date >= "2011-12-01" & date <= "2012-02-01", "Wtr1112", Season),
             Season = ifelse(date >= "2012-12-01" & date <= "2013-02-01", "Wtr1213", Season),
             Season = ifelse(date >= "2013-12-01" & date <= "2014-02-01", "Wtr1314", Season),
             Season = ifelse(date >= "2014-12-01" & date <= "2015-02-01", "Wtr1415", Season),
             Season = ifelse(date >= "2015-12-01" & date <= "2016-02-01", "Wtr1516", Season),
             Season = ifelse(date >= "2016-12-01" & date <= "2017-02-01", "Wtr1617", Season),
             Season = ifelse(date >= "2017-12-01" & date <= "2018-02-01", "Wtr1718", Season),
             Season = ifelse(date >= "2018-12-01" & date <= "2019-02-01", "Wtr1819", Season),
             Season = ifelse(date >= "2019-12-01" & date <= "2020-02-01", "Wtr1920", Season),
             Season = ifelse(date >= "2020-12-01" & date <= "2021-02-01", "Wtr2021", Season),
             Season = ifelse(date >= "2021-12-01" & date <= "2022-02-01", "Wtr2122", Season),
             Season = ifelse(date >= "2022-12-01" & date <= "2023-02-01", "Wtr2223", Season))
    
    return(annual_winter_season)
  }
  precip_season <- add_season(precip)
  temp_season <- add_season(temp)
  
  #'  Function to reformat and summarize weather data
  format_weather <- function(prism) {
    avg_winter_weather <- prism %>%
      #'  Reformat columns
      transmute(Clusters = Clusters, #NewLocationID = NewLocationID,
                GMU = GMU,
                Season = Season,
                Date = date,
                meanWeather = meanMonthlyValue) %>%
      #'  Average monthly weather by year and site
      # group_by(NewLocationID, Season) %>%
      group_by(GMU, Clusters, Season) %>%
      summarise(DecFeb_meanWeather = mean(meanWeather),
                DecFeb_meanWeather_se = sd(meanWeather)/sqrt(nrow(.))) %>%
      ungroup() %>%
      #'  Standardize average monthly weather (mean = 0, SD = 1)
      #'  This represents the monthly mean relative to the 50 year average and SD
      mutate(DecFeb_meanWeather_z = as.numeric(scale(DecFeb_meanWeather)))
    
    #'  Double check standardized data
    print(mean(avg_winter_weather$DecFeb_meanWeather_z))
    print(sd(avg_winter_weather$DecFeb_meanWeather_z))
    
    return(avg_winter_weather)
  }
  wtr_totalPrecip <- format_weather(precip_season) %>%
    rename("DecFeb_meanPPT_mm" = "DecFeb_meanWeather") %>%
    rename("DecFeb_meanPPT_se" = "DecFeb_meanWeather_se") %>%
    rename("DecFeb_meanPPT_z" = "DecFeb_meanWeather_z")
  wtr_minTemp <- format_weather(temp_season) %>%
    rename("DecFeb_meanMinTemp_C" = "DecFeb_meanWeather") %>%
    rename("DecFeb_meanMinTemp_se" = "DecFeb_meanWeather_se") %>%
    rename("DecFeb_meanMinTemp_z" = "DecFeb_meanWeather_z") 
  
  #'  Generate Winter Severity Index (WSI) per year and Cluster
  wsi <- full_join(wtr_totalPrecip, wtr_minTemp, by = c("GMU", "Clusters", "Season")) %>%
    #'  Multiply standardized precip by standardized temp data
    mutate(DecFeb_WSI = DecFeb_meanPPT_z * DecFeb_meanMinTemp_z) %>%
    dplyr::select(c("GMU", "Clusters", "Season", "DecFeb_WSI")) %>%
    #'  Filter to the winters preceding Summer 2020, 2021, 2022, and 2023 surveys
    filter(Season == "Wtr1920" | Season == "Wtr2021" | Season == "Wtr2122" | Season == "Wtr2223") %>%
    rename("Winter" = "Season") %>%
    #'  Add column that references year of the previous summer (i.e., Wtr1920 follows summer 2019)
    #'  This ensures that data from most recent winter affects current summer RAI
    #'  (6-mo lag between WSI and summer RAI) based how data are stacked and lagged 
    #'  in SEM analyses (e.g., WSI Wtr2021 is hypothesized to affect Summer 2021 RAI)
    mutate(Year = ifelse(Winter == "Wtr1920", 2019, 2022),
           Year = ifelse(Winter == "Wtr2021", 2020, Year),
           Year = ifelse(Winter == "Wtr2122", 2021, Year)) %>%
    relocate(Year, .after = "Clusters") %>%
    rename("ClusterID" = "Clusters") %>%
    filter(Year != 2019)
  print(wsi)
  
  #'  -----------------------------
  #####  Percent disturbed forest  #####
  #'  -----------------------------
  #'  Hansen's Global Forest Change dataset (area of 2020 canopy cover per cluster)
  gfc <- read_csv("./Data/GEE outputs/GFC_2020_canopy_cover_area_clusters.csv") %>% 
    dplyr::select(-`.geo`) %>%
    #'  Relabel GEE indexing to corresponding GMUs
    mutate(GMU = ifelse(str_detect(`system:index`, "1_1_"), "GMU1", `system:index`),  # note the order is different with PRISM data
           GMU = ifelse(str_detect(GMU, "1_2_"), "GMU6", GMU),
           GMU = ifelse(str_detect(GMU, "2_"), "GMU10A", GMU)) %>%
    dplyr::select(-`system:index`) %>%
    relocate(CanopyCoverAreaSqKm, .after = area_km2) %>%
    rename("ClusterID" = "Clusters") %>%
    #'  Column to double check canopy cover area is less than or equal to area of polygon (if FALSE, that's a problem)
    mutate(areas_ok = ifelse(CanopyCoverAreaSqKm <= area_km2, TRUE, FALSE))
  
  #'  Hansen's Global Forest Change dataset (area and year of canopy loss per cluster)
  gfc_loss <- read_csv("./Data/GEE outputs/GFC_annual_canopy_loss_area_clusters.csv") %>% 
    dplyr::select(-`.geo`) %>%
    #'  Relabel GEE indexing to corresponding GMUs
    mutate(GMU = ifelse(str_detect(`system:index`, "_1_1_"), "GMU1", `system:index`),  # note the order is different with PRISM data
           GMU = ifelse(str_detect(GMU, "_1_2_"), "GMU6", GMU),
           GMU = ifelse(str_detect(GMU, "_2_"), "GMU10A", GMU)) %>%
    dplyr::select(-`system:index`) %>%
    full_join(gfc, by = c("GMU", "ClusterID")) %>%
    #'  Make canopy loss year easier to interpret
    mutate(CanopyLossYear = Year_add1 + 2001,
           #'  Calculate area canopy loss within each cluster
           CanopyLossArea_sq_km = CanopyLossArea_sq_m/1000000) %>%
    dplyr::select(c(GMU, ClusterID, area_km2, CanopyCoverAreaSqKm, CanopyLossYear, CanopyLossArea_sq_km)) 
  
  #'  Function to calculate remaining canopy cover after multiple years of loss
  gfc_loss_accumulation <- function(Hanson, season, yr) {
    gfc_loss <- Hanson %>%
      filter(CanopyLossYear <= yr) %>%
      #'  Accumulated canopy loss up to defined year
      group_by(GMU, ClusterID) %>%
      reframe(area_km2 = area_km2,
              CanopyCoverAreaSqKm = CanopyCoverAreaSqKm,
              TotalLoss_sq_km = sum(CanopyLossArea_sq_km),
              CanopyCover_sq_km = CanopyCoverAreaSqKm - TotalLoss_sq_km) %>%
      #'  Reduce to single observation per cluster
      group_by(GMU, ClusterID) %>%
      slice(1L) %>%
      arrange(GMU, ClusterID) 
    return(gfc_loss)
  }
  #'  Accumulated loss through 2001
  gfc_loss_thru_2001 <- gfc_loss_accumulation(gfc_loss, yr = 2001) %>%
    rename("TotalLoss2001_sq_km" = "TotalLoss_sq_km") %>%
    rename("CanopyCover2001_sq_km" = "CanopyCover_sq_km")
  #  Accumulated loss through 2002
  gfc_loss_thru_2002 <- gfc_loss_accumulation(gfc_loss, yr = 2002) %>%
    rename("TotalLoss2002_sq_km" = "TotalLoss_sq_km") %>%
    rename("CanopyCover2002_sq_km" = "CanopyCover_sq_km")
  #'  Accumulated loss over last 19 years (2000 - 2019)
  gfc_loss_thru_2019 <- gfc_loss_accumulation(gfc_loss, yr = 2019) %>%
    #'  Calculate proportion of forest cover that has been disturbed in last 20-years
    #'  Actually 19 years here b/c remotely sensed data only goes back to 2000
    mutate(DisturbedForest_last20Yrs = TotalLoss_sq_km/CanopyCoverAreaSqKm,
           #'  Calculate proportion of cluster that has canopy cover
           PercentCover = CanopyCover_sq_km/area_km2,
           #'  Year accumulated loss goes through (i.e., remaining canopy cover in this year)
           AccumulatedLoss_thru = 2019,
           #'  Year of current forest condition
           Year = 2019)
  #'  Accumulated loss over last 20 years (2000 - 2020)
  gfc_loss_thru_2020 <- gfc_loss_accumulation(gfc_loss, yr = 2020) %>%
    #'  Calculate proportion of forest cover that has been disturbed in last 20-years
    mutate(DisturbedForest_last20Yrs = TotalLoss_sq_km/CanopyCoverAreaSqKm,
           #'  Calculate proportion of cluster that has canopy cover
           PercentCover = CanopyCover_sq_km/area_km2,
           AccumulatedLoss_thru = 2020,
           Year = 2020) 
  #'  Accumulated loss over last 20 years (2001 - 2021)
  gfc_loss_thru_2021 <- gfc_loss_accumulation(gfc_loss, yr = 2021) %>%
    full_join(gfc_loss_thru_2001) %>%
    #'  Adjust accumulated loss through 2021 so TotalLoss is only for past 20 years
    #'  by updating CanopyCoverArea to that in 2001 and subtracting 2001 accumulated
    #'  loss from 2021 accumulated loss
    mutate(CanopyCoverAreaSqKm = CanopyCover2001_sq_km,
           TotalLoss_sq_km = TotalLoss_sq_km - TotalLoss2001_sq_km,
           CanopyCover_sq_km = CanopyCoverAreaSqKm - TotalLoss_sq_km,
           #'  Calculate proportion of forest cover that has been disturbed in last 20-years
           DisturbedForest_last20Yrs = TotalLoss_sq_km/CanopyCoverAreaSqKm,
           #'  Calculate proportion of cluster that has canopy cover
           PercentCover = CanopyCover_sq_km/area_km2,
           AccumulatedLoss_thru = 2021,
           Year = 2021) %>%
    dplyr::select(-c(TotalLoss2001_sq_km, CanopyCover2001_sq_km)) 
  #'  Accumulated loss over last 20 years (2002 - 2022)  
  gfc_loss_thru_2022 <- gfc_loss_accumulation(gfc_loss, yr = 2022) %>%
    full_join(gfc_loss_thru_2002) %>%
    #'  Adjust accumulated loss through 2022 so TotalLoss is only for past 20 years
    mutate(CanopyCoverAreaSqKm = CanopyCover2002_sq_km,
           TotalLoss_sq_km = TotalLoss_sq_km - TotalLoss2002_sq_km,
           CanopyCover_sq_km = CanopyCoverAreaSqKm - TotalLoss_sq_km,
           #'  Calculate proportion of forest cover that has been disturbed in last 20-years
           DisturbedForest_last20Yrs = TotalLoss_sq_km/CanopyCoverAreaSqKm,
           #'  Calculate proportion of cluster that has canopy cover
           PercentCover = CanopyCover_sq_km/area_km2,
           AccumulatedLoss_thru = 2022,
           Year = 2022) %>%
    dplyr::select(-c(TotalLoss2002_sq_km, CanopyCover2002_sq_km))
  
  #'  Join proportion of forest cover that has been disturbed in last 20-years for 2020-2022
  percDisturbedForest <- full_join(gfc_loss_thru_2020, gfc_loss_thru_2021) %>%
    full_join(gfc_loss_thru_2022) %>% 
    dplyr::select(c(GMU, ClusterID, Year, AccumulatedLoss_thru, DisturbedForest_last20Yrs, PercentCover)) 
  
  
  #'  -------------------------------
  ####  BGMR wolf harvest locations  ####
  #'  -------------------------------
  #'  Read in harvest data and reformat
  #'  Data prepped by Kat Peters
  KP_dat <- list.files(path = "./Data/IDFG BGMR data/Sibgroup_Harvest_Locations/", pattern = "\\.csv$", full.names = T) %>% 
    map_df(~read_csv(., col_types = cols(.default = "c"))) %>%
    mutate(Harvest_date = as.Date(MortDate, "%m/%d/%Y"),
           Latitude = as.numeric(Latitude),
           Longitude = as.numeric(Longitude),
           s_cell = NA) %>%
    dplyr::select(c(GMU, BGMR, Latitude, Longitude, Location_Description, Harvest_date, Source, A.Sex, s_cell)) 
  #'  Data from IDFG
  bgmr_15to21 <- read_csv("./Data/IDFG BGMR data/IDFG Wolf Harvest Data/BGMR_wolf_harv_June1520toDec2021-KO.csv") %>%
    dplyr::select(c(`Game Management Unit`, `Form Number`, Latitude, Longitude, `Kill Location`, `Kill Date`, `Mortality Agent`, Sex, s_cell)) %>%
    mutate(`Kill Date` = as.Date(`Kill Date`, "%Y-%m-%d"),
           Latitude = as.numeric(Latitude),
           Longitude = as.numeric(Longitude),
           Sex = ifelse(Sex == "M", "Male", "Female"),
           s_cell = ifelse(s_cell == 533, 53, s_cell))
  names(bgmr_15to21) <- c("GMU", "BGMR", "Latitude", "Longitude", "Location_Description", "Harvest_date", "Source", "A.Sex", "s_cell")
  bgmr_21to23 <- read_csv("./Data/IDFG BGMR data/IDFG Wolf Harvest Data/21-23 Occ SCell_regionsinput.csv") %>%
    dplyr::select(c(`Game Management Unit`, `Form Number`, Latitude, Longitude, `Kill Location`, `Kill Date`, `Mortality Agent`, Sex, `S-cell`)) %>%
    mutate(`Kill Date` = as.Date(`Kill Date`, "%m/%d/%Y"),
           Latitude = as.numeric(Latitude),
           Longitude = as.numeric(Longitude))
  names(bgmr_21to23) <- c("GMU", "BGMR", "Latitude", "Longitude", "Location_Description", "Harvest_date", "Source", "A.Sex", "s_cell")
  #'  Merge and format
  bgmr_dat <- bind_rows(bgmr_15to21, bgmr_21to23) %>%
    bind_rows(KP_dat) %>%
    #'  Filter to only observations from GMU1, 6, & 10A
    filter(GMU == 1 | GMU == 6 | GMU == "10A") %>%
    arrange(GMU, BGMR) %>%
    #'  Remove duplicates, keeping the observation with lat/long coordinates
    group_by(BGMR) %>%
    slice_max(!is.na(Latitude), with_ties = T) %>%
    slice_max(!is.na(s_cell), with_ties = T) %>%
    ungroup() %>%
    arrange(GMU, BGMR) %>%
    unique(.) %>%
    filter(Source != "Wildlife Services")
  
  #' #'  Save BGMR data and use wolf grid cells to help identify general harvest locations, then read in updated dataset
  #' write_csv(bgmr_dat, "./Data/IDFG BGMR data/IDFG Wolf Harvest Data/GMU1_6_10A_filtered.csv")
  #' s_cells <- st_read("./Shapefiles/IDFG spatial data/SCells_Strata.shp")
  #' library(mapview)
  #' mapview::mapview(s_cells, zcol = "SCell")
  
  #'  Read in updated BGMR data with approximate harvest coordinates
  bgmr_locs <- read_csv("././Data/IDFG BGMR data/IDFG Wolf Harvest Data/GMU1_6_10A_assigned_harvest_locations.csv") %>%
    mutate(Harvest_date = as.Date(Harvest_date, "%m/%d/%Y"),
           Latitude = as.numeric(Latitude),
           Longitude = as.numeric(Longitude)) %>%
    dplyr::select(c(-s_cell, BGMR)) %>%
    #'  Remove observations that were not legal harvest
    filter(Source != "Illegal Kill" & Source != "Depredation Kill") %>%
    #'  Exclude harvest occurring 1 year before study period
    filter(Harvest_date >= "2019-06-01") %>%
    #'  Exclude observations missing location coordinates
    filter(!is.na(Latitude))
  
  #'  Filter observations to the year prior to each study period and make spatial
  #'  e.g., June 1 of the previous yr to June 1 of the study year
  #'  Hypothesize harvest from within past year will influence relative wolf density that summer
  bgmr_Jun19_Jun20 <- filter(bgmr_locs, Harvest_date <= "2020-06-01") %>% mutate(ID = seq(1:nrow(.))) %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs") #%>% st_transform(crs = crs(cluster_poly))
  bgmr_Jun20_Jun21 <- filter(bgmr_locs, Harvest_date > "2020-06-01" & Harvest_date <= "2021-06-01") %>% mutate(ID = seq(1:nrow(.))) %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")
  bgmr_Jun21_Jun22 <- filter(bgmr_locs, Harvest_date > "2021-06-01" & Harvest_date <= "2022-06-01") %>% mutate(ID = seq(1:nrow(.))) %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")
  bgmr_Jun22_Jun23 <- filter(bgmr_locs, Harvest_date > "2022-06-01" & Harvest_date <= "2023-06-01") %>% mutate(ID = seq(1:nrow(.))) %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Intersect harvest locations with clusters, sum number of wolves harvested per cluster and divide by cluster area (harvest/km2)
  sf_use_s2(FALSE)
  cluster_harvest <- function(bgmr, season, yr) {
    clusters <- as.data.frame(cluster_poly) %>% dplyr::select(-geometry) 
    harvest_polygon <- st_intersection(bgmr, cluster_poly) %>% as.data.frame(.) %>%
      dplyr::select(c(GMU, BGMR, Location_Description, Harvest_date, Source, ClusterID, area_km2)) %>%
      mutate(GMU = paste0("GMU", GMU)) %>%
      group_by(GMU, ClusterID) %>%
      reframe(area_km2 = area_km2,
              annual_harvest = n(),
              harvest_sqKm = annual_harvest/area_km2) %>%
      group_by(GMU, ClusterID) %>%
      slice(1L) %>%
      ungroup() %>%
      full_join(clusters, harvest_polygon, by = c("GMU", "ClusterID", "area_km2")) %>%
      mutate(annual_harvest = ifelse(is.na(annual_harvest), 0, annual_harvest),
             harvest_sqKm = ifelse(is.na(harvest_sqKm), 0, harvest_sqKm), 
             harvest_sqKm = round(harvest_sqKm, 4),
             #'  Add column representing time period of harvest
             harvest_season = season,
             #'  Column representing start of harvest season (e.g., 2019 = summer 
             #'  that started harvest from June 2019 - June 2020). This ensures 
             #'  that data from past year (past 12 months of harvest) affects current 
             #'  summer RAI based how data are stacked and lagged in SEM analyses 
             #'  (e.g., harvest Jun19-Jun20 is hypothesized to affect Summer 2020 RAI)
             Year = yr) %>%
      arrange(GMU, ClusterID) %>%
      dplyr::select(-area_km2)
    return(harvest_polygon)
  }
  harvest_Jun19_Jun20 <- cluster_harvest(bgmr_Jun19_Jun20, season = "Jun19toJun20", yr = 2019)
  harvest_Jun20_Jun21 <- cluster_harvest(bgmr_Jun20_Jun21, season = "Jun20toJun21", yr = 2020)
  harvest_Jun21_Jun22 <- cluster_harvest(bgmr_Jun21_Jun22, season = "Jun21toJun22", yr = 2021)
  harvest_Jun22_Jun23 <- cluster_harvest(bgmr_Jun22_Jun23, season = "Jun22toJun23", yr = 2022)
  sf_use_s2(TRUE)
  harvest_intensity <- bind_rows(harvest_Jun19_Jun20, harvest_Jun20_Jun21, harvest_Jun21_Jun22, harvest_Jun22_Jun23) %>%
    relocate(harvest_season, .before = annual_harvest)
  
  #'  Merge forest and WSI data together
  covs <- full_join(percDisturbedForest, wsi, by = c("GMU", "ClusterID", "Year")) %>%
    full_join(harvest_intensity, by = c("GMU", "ClusterID", "Year")) %>%
    filter(GMU != "GMU1" | Year != 2020)
  
  
  
  ####  SCRAPPED  ####
  #' #'  National Landcover Database (frequency of landcover class per cluster)
  #' nlcd <- read_csv("./Data/GEE outputs/NLCD_frequencies_2019_2021_clusters.csv") %>%
  #'   #'  Relabel GEE indexing to corresponding GMUs
  #'   mutate(GMU = ifelse(str_detect(`system:index`, "1_1_"), "GMU1", `system:index`),  # note the order is different with PRISM data
  #'          GMU = ifelse(str_detect(GMU, "1_2_"), "GMU6", GMU),
  #'          GMU = ifelse(str_detect(GMU, "2_"), "GMU10A", GMU)) %>%
  #'   dplyr::select(c(GMU, Clusters, landcover_2019, landcover_2021)) %>%
  #'   rename("ClusterID" = "Clusters") %>%
  #'   mutate(landcover_2019 = str_replace_all(landcover_2019, "[[{}]]", ""),
  #'          landcover_2021 = str_replace_all(landcover_2021, "[[{}]]", ""))
  #' nlcd19 <- nlcd %>% dplyr::select(c(GMU, ClusterID, landcover_2019)) %>%
  #'   #'  Remove GMU1 cameras from this dataset
  #'   filter(GMU != "GMU1") %>%
  #'   #'  Separate wonky GEE character string grouping all landcover outputs together 
  #'   #'  (FYI separate does something weird with extra columns so be sure to provide 
  #'   #'  double the max number of possible landcover types). Ignore warning about missing pieces and NAs.
  #'   separate(landcover_2019, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17", "col18", "col19", "col20", "col21", "col22"), sep = "[, ]")
  #' nlcd21 <- nlcd %>% dplyr::select(c(GMU, ClusterID, landcover_2021)) %>%
  #'   separate(landcover_2021, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17", "col18", "col19", "col20", "col21", "col22"), sep = "[, ]")
  #' 
  #' #'  Function to reformat NLCD data extracted from GEE
  #' habitat_frequencies <- function(dat) {
  #'   reformatted <- dat %>%
  #'     #'  Convert to long format to organize columns better
  #'     pivot_longer(!c(GMU, ClusterID), names_to = "col_name", values_to = "landcover_class") %>%
  #'     #'  Split landcover class from number of pixels
  #'     separate(landcover_class, into = c("NLCD_class", "npixels"), sep = "[=]") %>%
  #'     #'  Drop unneeded columns and rows
  #'     dplyr::select(-col_name) %>%
  #'     filter(!is.na(npixels)) %>%
  #'     #'  Reformat npixels 
  #'     mutate(npixels = as.numeric(npixels)) %>%
  #'            # npixels = round(npixels, 2)) %>%
  #'     #'  Calculate proportion of pixels made up by each cover type near camera
  #'     group_by(GMU, ClusterID) %>%
  #'     reframe(NLCD_class = NLCD_class,
  #'             npixels = npixels,
  #'             totalPix = sum(npixels),
  #'             pixel_area_m2 = totalPix * (30*30),
  #'             pixel_area_km2 = pixel_area_m2/1000000,
  #'             percentPix = npixels/totalPix) %>%
  #'     ungroup() %>%
  #'     mutate(percentPix = round(percentPix, 3)) %>%
  #'     arrange(GMU, ClusterID) %>%
  #'     # rename("NewLocationID" = "NwLctID") %>%
  #'     #'  Define NLCD landcover classifications using https://developers.google.com/earth-engine/datasets/catalog/USGS_NLCD_RELEASES_2021_REL_NLCD
  #'     mutate(habitat_type = ifelse(NLCD_class == "11", "Other", NLCD_class),               #11: Open water
  #'            habitat_type = ifelse(NLCD_class == "12", "Other", habitat_type),             #12: Perennial ice/snow
  #'            habitat_type = ifelse(NLCD_class == "21", "Developed", habitat_type),         #21: Developed, open space: mixture of constructed materials, mostly vegetation in the form of lawn grasses
  #'            habitat_type = ifelse(NLCD_class == "22", "Developed", habitat_type),         #22: Developed, low intensity: areas with a mixture of constructed materials and vegetation, most commonly single-family housing units
  #'            habitat_type = ifelse(NLCD_class == "23", "Developed", habitat_type),         #21: Developed, medium intensity: areas with a mixture of constructed materials and vegetation, most commonly single-family housing units
  #'            habitat_type = ifelse(NLCD_class == "24", "Developed", habitat_type),         #22: Developed, high intensity: highly developed areas where people reside or work in high numbers
  #'            habitat_type = ifelse(NLCD_class == "31", "Other", habitat_type),             #31: Barren Land
  #'            habitat_type = ifelse(NLCD_class == "41", "Forested", habitat_type),          #41: Deciduous Forest
  #'            habitat_type = ifelse(NLCD_class == "42", "Forested", habitat_type),          #42: Evergreen Forest
  #'            habitat_type = ifelse(NLCD_class == "43", "Forested", habitat_type),          #43: Mixed Forest
  #'            habitat_type = ifelse(NLCD_class == "52", "Shrubland", habitat_type),         #52: Shrub/Scrub
  #'            habitat_type = ifelse(NLCD_class == "71", "Grassland", habitat_type),         #71: Grassland/Herbaceous
  #'            habitat_type = ifelse(NLCD_class == "81", "Agriculture", habitat_type),       #81: Pasture/Hay
  #'            habitat_type = ifelse(NLCD_class == "82", "Agriculture", habitat_type),       #82: Cultivated Crops
  #'            habitat_type = ifelse(NLCD_class == "90", "Riparian_woodland", habitat_type),      #90: Woody Wetlands
  #'            habitat_type = ifelse(NLCD_class == "95", "Riparian_wetland", habitat_type)) %>%   #95: Emergent Herbaceous Wetlands)
  #'     relocate(habitat_type, .after = "ClusterID")
  #'   
  #'   #'  Filter data to the dominant habitat class within each cluster
  #'   #'  (dominant is defined as the landcover class comprising the largest proportion
  #'   #'  of pixels within cluster)
  #'   dominant_habitat <- reformatted %>%
  #'     group_by(GMU, ClusterID) %>%
  #'     slice_max(order_by = percentPix, n = 1) %>%
  #'     ungroup() %>%
  #'     arrange(GMU, ClusterID)
  #'   
  #'   #'  Review dominant habitat types
  #'   print(table(dominant_habitat$habitat_type))
  #'   
  #'   #'  List both datasets together
  #'   landcover_list <- list(dominant_habitat, reformatted)
  #'   names(landcover_list) <- c("dominant_habitat", "percent_landcover")
  #'   
  #'   return(landcover_list)
  #'   return()
  #' }
  #' landcover19 <- habitat_frequencies(nlcd19)
  #' landcover21 <- habitat_frequencies(nlcd21)
  #' 
  #' #'  Load MTBS burn perimeter data 
  #' mtbs <- st_read(dsn = "./Shapefiles/MTBS_perimeter_data", layer = "mtbs_perims_DD") %>%
  #'   st_transform(forest_proj)
  #' 
  #' #'  Extract burn year of each pixel per cluster
  #' burn_perimeters <- function(polygon) {
  #'   #'  Vectorize polygons
  #'   polygon <- st_transform(polygon, forest_proj) #%>% vect(.)
  #'   #'  Intersect overlapping clusters & burn perimeters
  #'   mtbs_simple <- mtbs %>% 
  #'     dplyr::select("Ig_Date", "Incid_Name", "Incid_Type") %>%
  #'     mutate(Burn_year = format(as.Date(Ig_Date, format="%Y-%m-%d"),"%Y")) 
  #'   mtbs_cluster_inter <- st_intersection(polygon, mtbs_simple) 
  #'   burned_area <- as.numeric(st_area(mtbs_cluster_inter)/1000000) %>%
  #'     as.data.frame()
  #'   names(burned_area) <- "burned_area_km2"
  #'   mtbs_cluster_inter <- bind_cols(mtbs_cluster_inter, burned_area)
  #'   mapview::mapview(mtbs_cluster_inter, zcol = "Clusters")
  #'   #'  Save burn year and area
  #'   mtbs_burn_yr <- as.data.frame(mtbs_cluster_inter) %>%
  #'     dplyr::select(c(GMU, Clusters, Burn_year, area_km2, burned_area_km2)) %>%
  #'     rename("ClusterID" = "Clusters") 
  #'   return(mtbs_burn_yr)
  #' }
  #' burned_gmu1 <- burn_perimeters(gmu1_poly)
  #' burned_gmu6 <- burn_perimeters(gmu6_poly)
  #' burned_gmu10a <- burn_perimeters(gmu10a_poly)
  #' 
  #'  
  #' #'  Add GEE data to larger covariate df
  #' format_cam_covs <- function(dat, landcov, season, camYr) { 
  #'   covs <- full_join(dat, gfc, by = "ClusterID") %>%
  #'     full_join(gfc, by = "ClusterID") %>%
  #'     relocate(Burn_year, .after = "percentPix") %>%
  #'     mutate(Season = season,
  #'            #'  Calculate number of years since burn/canopy loss
  #'            YrsSinceBurn = camYr - as.numeric(Burn_year),
  #'            YrsSinceBurn = ifelse(YrsSinceBurn <0, NA, YrsSinceBurn),
  #'            YrsSinceLoss = camYr - CanopyLossYear,
  #'            YrsSinceLoss = ifelse(YrsSinceLoss <0, NA, YrsSinceLoss),
  #'            #'  Categorize years since burn/canopy loss following Barker et al. (2018) and Ganz et al. (2024)
  #'            DisturbanceLoss = ifelse(YrsSinceLoss <= 20, "Loss_1_20", habitat_type),
  #'            DisturbanceBurn = ifelse(YrsSinceBurn <= 5, "Burn_1_5", habitat_type),
  #'            DisturbanceBurn = ifelse(YrsSinceBurn > 5 & YrsSinceBurn <=10, "Burn_6_10", DisturbanceBurn),
  #'            DisturbanceBurn = ifelse(YrsSinceBurn > 10 & YrsSinceBurn <=15, "Burn_10_15", DisturbanceBurn),
  #'            DisturbanceBurn = ifelse(YrsSinceBurn > 15 & YrsSinceBurn <=20, "Burn_16_20", DisturbanceBurn),
  #'            DisturbanceBurn = ifelse(YrsSinceBurn > 20, "Burn_over20", DisturbanceBurn),
  #'            #'  Generate a single habitat class covariate representing dominant habitat type and years since burn/canopy loss (if forested)
  #'            Habitat_class = habitat_type,
  #'            Habitat_class = ifelse(!is.na(DisturbanceLoss) & Habitat_class == "Forested", DisturbanceLoss, Habitat_class),
  #'            Habitat_class = ifelse(!is.na(DisturbanceBurn) & Habitat_class == "Loss_1_20", DisturbanceBurn, Habitat_class),
  #'            #'  Grab percentPix for sites with disturbance to indicate percent disturbed forest within last 20 years
  #'            #'  (0 % disturbed forest if the dominant habitat class is unburned/unlogged forest or any other landcover type)
  #'            PercDisturbedForest = ifelse(Habitat_class == "Loss_1_20" | Habitat_class == "Burn_1_5" |
  #'                                           Habitat_class == "Burn_6_10" | Habitat_class == "Burn_16_20", percentPix, 0)) %>%
  #'     left_join(wsi, by = c("NewLocationID", "Season")) %>%
  #'     relocate(Season, .after = "GMU") %>%
  #'     filter(!is.na(GMU))
  #' 
  #'   #'  Review new habitat classes
  #'   print(table(covs$Habitat_class))
  #'   
  #'   return(covs)
  #' }
  #' #'  Add GEE data to larger covariate dataframe; NOTE different landcover data applied to 2020 vs 2021 & 2022 data
  #' cams_eoe20s <- format_cam_covs(covariate_list[[1]], landcov = landcover19[[1]], camYr = 2020, year = "2020")
  #' cams_eoe21s <- format_cam_covs(covariate_list[[2]], landcov = landcover21[[1]], camYr = 2021, year = "2021")
  #' cams_eoe22s <- format_cam_covs(covariate_list[[3]], landcov = landcover21[[1]], camYr = 2022, year = "2022")