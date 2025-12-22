  #'  --------------------------------
  #'  Prep Data for Structural Equation Models
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2024
  #'  --------------------------------
  #'  Format and explore data before building SEMs
  
  #' #'  Clean workspace
  #' rm(list = ls())
  
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
  clusters_gmu1 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu1_06.30.25.shp") %>%
    reformat_clusters(.)
  clusters_gmu6 <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu6_06.30.25.shp") %>%
    reformat_clusters(.)
  clusters_gmu10a <- st_read("./Shapefiles/IDFG spatial data/Camera_locations/Camera_clusters/cam_clusters_gmu10a_06.30.25.shp") %>%
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
  
  #' #'  Run script that formats covariate data
  #' source("./Scripts/Structural_Equation_Models/Format_covariate_data_for_SEMs.R")
  
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
      reframe(nCam = length(RN.n),
              SppN = sum(RN.n),
              meanSppN = SppN / nCam,
              SppDensity.km2 = meanSppN/area_km2,
              SppDensity.100km2 = SppDensity.km2*100,
              SppN.r = sum(round(RN.n, 0)),
              meanSppN.r = SppN.r / nCam,
              SppDensity.km2.r = meanSppN.r/area_km2,
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
  
  # save(cluster_density, file = "./Outputs/Relative_Abundance/RN_model/RelativeDensityIndex_per_SppCluster_06.30.25.RData")
  
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
  
 
  #'  -----------------------------------------
  ####  Format local abundance estimates & SD  ####
  #'  -----------------------------------------
  #'  Long to wide data structure
  wide_data <- function(dat) {
    pivot_data_wide <- as.data.frame(dat) %>%
      filter(!is.na(ClusterID)) %>%
      #'  Create categorical year variable
      mutate(Year = as.numeric(Year),
             year = ifelse(Year == 2020, "yr1", Year),
             year = ifelse(Year == 2021, "yr2", year),
             year = ifelse(Year == 2022, "yr3", year)) %>% #,
             #' #'  Add small value to density estimates that are 0 (needed for boxcox trnasformatio below)
             #' SppDensity.100km2.r = ifelse(SppDensity.100km2.r == 0, 0.00001, SppDensity.100km2.r)) %>%
      dplyr::select(c(GMU, ClusterID, Year, year, Species, SppDensity.100km2.r)) %>%
      #'  Create column per species with their site-specific local abundance 
      pivot_wider(names_from = "Species",
                  values_from = c("SppDensity.100km2.r")) %>%
      left_join(covs, by = c("GMU", "Year", "ClusterID"))
    #'  Review to make sure it time lags look right given these observations are 
    #'  hypothesized to affect the next set of observations
    print(pivot_data_wide)
    #'  Drop unneeded columns used to help double check time lags
    pivot_data_wide <- pivot_data_wide %>%
      dplyr::select(-c(AccumulatedLoss_thru, Winter, harvest_season))
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
    hist(dat$bear_black, main = paste("Relative bear density across clusters, \nsummer", yr))
    hist(dat$coyote, main = paste("Relative coyote density across clusters, \nsummer", yr))
    hist(dat$mountain_lion, main = paste("Relative lion density across clusters, \nsummer", yr))
    hist(dat$wolf, main = paste("Relative wolf density across clusters, \nsummer", yr))
    hist(dat$elk, main = paste("Relative elk density across clusters, \nsummer", yr))
    hist(dat$moose, main = paste("Relative moose density across clusters, \nsummer", yr))
    hist(dat$whitetailed_deer, main = paste("Relative wtd density across clusters, \nsummer", yr))
    hist(dat$DisturbedForest_last20Yrs, main = paste("Percent disturbed forest, \nsummer", yr))
    hist(dat$DecFeb_WSI, main = paste("Dec-Feb Winter Severity Index, \nstarting Dec", yr))
    hist(dat$harvest_sqKm, main = paste("Wolves harvested/km2 over past 12 months, \nstarting June", yr))
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
      rename_with(.cols = bear_black:harvest_sqKm, function(x){paste0(x, ".", yr)}) %>% 
      dplyr::select(-c(Year, year))
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
    arrange(GMU, ClusterID) %>%
    mutate(ClusterID = as.character(ClusterID),
           harvest_sqKm_quad.yr1 = harvest_sqKm.yr1^2,
           harvest_sqKm_quad.yr2 = harvest_sqKm.yr2^2)
  head(density_wide_annual_20s_22s); tail(density_wide_annual_20s_22s)
  
  #'  Z-transform local abundance estimates (per year b/c annual estimates are stand alone variables)
  localN_z <- density_wide_annual_20s_22s %>%
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
  summary(lm(bear_black.yr2 ~ DisturbedForest_last20Yrs.yr1, data = localN_z))  
  summary(lm(bear_black.yr3 ~ DisturbedForest_last20Yrs.yr2, data = localN_z))  
  summary(lm(wolf.yr2 ~ harvest_sqKm.yr1, data = localN_z))
  summary(lm(wolf.yr3 ~ harvest_sqKm.yr2, data = localN_z))
  summary(lm(wolf.yr2 ~ annual_harvest.yr1, data = localN_z))
  summary(lm(wolf.yr3 ~ annual_harvest.yr2, data = localN_z))
  
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
    rename_with(.cols = bear_black:harvest_sqKm, function(x){paste0(x, ".Tminus1")}) %>%
    dplyr::select(-year)
  #'  Stack year 2 & year 3 data
  dat_t <- dat_stack_list[[2]] %>% #full_bx_dat_stack[[2]] %>%
    #'  Indicate whether observation was from first or second year in grouped data
    mutate(GroupYear = ifelse(year == "yr2", "first", "second")) %>%
    #'  Add time period identifier to each column name
    rename_with(.cols = bear_black:harvest_sqKm, function(x){paste0(x, ".T")}) %>%
    dplyr::select(-year)
  
  #'  Join t-1 and t data based on camera location
  density_wide_1YrLag_20s_22s <- full_join(dat_t_minus_1, dat_t, by = c("GMU", "ClusterID", "GroupYear")) %>%
    relocate(GroupYear, .after = ClusterID) %>%
    #'  Drop sites with NAs (missing 1+ years of data)
    na.omit(.) %>%
    mutate(ClusterID = as.character(ClusterID),
           harvest_sqKm_quad.Tminus1 = harvest_sqKm.Tminus1^2) 
  
  #'  Z-transform local abundance estimates (per year b/c annual estimates are stand alone variables)
  localN_z_1YrLag <- density_wide_1YrLag_20s_22s %>%
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
  summary(lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = localN_z_1YrLag))
  summary(lm(wolf.T ~ wolf.Tminus1 + harvest_sqKm.Tminus1, data = localN_z_1YrLag))
  
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
    mutate(ClusterID = as.character(ClusterID),
           Year = as.character(Year)) %>%
    mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    #mutate(across(starts_with(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat[,5:18]
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z_all) 

  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dat[,5:16] %>% as.matrix(.)
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
  
  
  
  
