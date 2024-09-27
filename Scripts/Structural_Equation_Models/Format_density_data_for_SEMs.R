  #'  --------------------------------
  #'  Prep Data for Structural Equation Models
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2024
  #'  --------------------------------
  #'  Run SEMs
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(piecewiseSEM)
  library(labelled)
  library(lme4)
  library(MASS)
  #install.packages("multcompView", repos="http://R-Forge.R-project.org")
  # library(multcompView)
  library(tidyverse)
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
  
  #'  Read in EoE GMUs shapefile
  eoe_gmu_wgs84 <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  
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
      reframe(SppN = sum(RN.n),
              SppDensity.km2 = SppN/area_km2,
              SppDensity.100km2 = SppDensity.km2*100,
              SppN.r = sum(round(RN.n, 0)),
              SppDensity.km2.r = SppN.r/area_km2,
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
  
  #'  -----------------------------------------
  ####  Format local abundance estimates & SD  ####
  #'  -----------------------------------------
  #'  Long to wide data structure
  wide_data <- function(dat) {
    pivot_data_wide <- as.data.frame(dat) %>%
      #'  Create categorical year variable
      mutate(year = ifelse(Year == 2020, "yr1", Year),
             year = ifelse(Year == 2021, "yr2", year),
             year = ifelse(Year == 2022, "yr3", year)) %>%
      dplyr::select(c(GMU, ClusterID, year, Species, SppDensity.100km2.r)) %>%
      #'  Create column per species with their site-specific local abundance 
      pivot_wider(names_from = "Species",
                  values_from = c("SppDensity.100km2.r")) 
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
    hist(dat$bear_black, main = paste("Relative bear density across clusters,", yr))
    hist(dat$coyote, main = paste("Relative coyote density across clusters,", yr))
    hist(dat$mountain_lion, main = paste("Relative lion density across clusters,", yr))
    hist(dat$wolf, main = paste("Relative wolf density across clusters,", yr))
    hist(dat$elk, main = paste("Relative elk density across clusters,", yr))
    hist(dat$moose, main = paste("Relative moose density across clusters,", yr))
    hist(dat$whitetailed_deer, main = paste("Relative wtd density across clusters,", yr))
  }
  plot_histogram(density_wide[[1]], yr = "2020")
  plot_histogram(density_wide[[2]], yr = "2021")
  plot_histogram(density_wide[[3]], yr = "2022")
  
  #'  ----------------------------------------------------------
  ####  Format & transform data for different model structures  ####
  #'  ----------------------------------------------------------
  #'  Grab just the local abundance estimates per species
  dat_n <- function(dat) {
    newdat <- dat[,4:12]
    return(newdat)
  }
  #'  For annual data
  density_wide_n <- lapply(density_wide, dat_n)
  #'  For stacked data
  density_stack_n <- lapply(dat_stack_list, dat_n)
  #'  For all data combined
  density_all_n <- dat_n(density_wide_allyrs)
  
  #'  Grab corresponding column names
  dat_colnames <- function(dat) {
    density_cols <- names(dat)
    return(density_cols)
  }
  density_wide_colnames <- lapply(density_wide_n, dat_colnames)
  density_stack_colnames <- lapply(density_stack_n, dat_colnames)
  density_all_colnames <- dat_colnames(density_all_n)
  
  #'  Create empty matrix to hold transformed variables
  empty_matrix <- function(dat) {
    bx_dat <- matrix(NA, nrow = nrow(dat), ncol = length(dat))
    return(bx_dat)
  }
  bx_dat <- lapply(density_wide_n, empty_matrix)
  bx_dat_stack <- lapply(density_stack_n, empty_matrix)
  bx_dat_all <- empty_matrix(density_all_n)
  
  #'  ----------------------------------------------
  #####  Annual format: independent annual columns  #####
  #'  ----------------------------------------------
  #'  Retaining each year as independent variable allows SEMs to assess indirect effects
  #'  and let relationships vary annually
  #'  ----------------------------------------------
  ######  Box-Cox Transformation  ######
  #'  ----------------------------
  #'  Loop through each column and pass boxcox function to linear model of data
  #'  Why this won't work in a normal function... I don't know
  #'  YEAR 1
  for(i in 1:length(density_wide_n[[1]])) {
    x <- pull(density_wide_n[[1]][,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat[[1]][,i] <- new_x_exact
  }
  colnames(bx_dat[[1]]) <- density_wide_colnames[[1]]
  
  #'  YEAR 2
  for(i in 1:length(density_wide_n[[2]])) {   #### THROWING AN ERROR OVER 0s
    x <- pull(density_wide_n[[2]][,i])
    #'  Add tiny value to x when x = 0
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat[[2]][,i] <- new_x_exact
  }
  colnames(bx_dat[[2]]) <- density_wide_colnames[[2]]
  
  #'  YEAR 3
  for(i in 1:length(density_wide_n[[3]])) {
    x <- pull(density_wide_n[[3]][,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat[[3]][,i] <- new_x_exact
  }
  colnames(bx_dat[[3]]) <- density_wide_colnames[[3]]
  
  #'  ---------------------------------------------
  ######  Merge local abundance, SD, & covariates  ######
  #'  ---------------------------------------------
  #'  Add annual covariates and site info back to dataset   
  add_covs_back <- function(olddat, newdat, datsd) {
    covs <- dplyr::select(olddat, c("NewLocationID", "CellID", "GMU", "Setup", 
                                    "Season", "Year", "habitat_class", 
                                    "PercDisturbedForest", "DecFeb_WSI"))
    dat <- bind_cols(covs, newdat) %>%
      bind_cols(datsd) %>%
      relocate(habitat_class, .after = last_col()) %>%
      relocate(PercDisturbedForest, .after = last_col()) %>%
      relocate(DecFeb_WSI, .after = last_col())
    dat <- as.data.frame(dat)
    return(dat)
  }
  full_bx_dat <- mapply(add_covs_back, olddat = RN_wide, newdat = bx_dat, datsd = RN_wide_sd, SIMPLIFY = FALSE)
  
  #'  Create wide data structure but this time one column per year for each species & covariate
  wide_data_by_year <- function(dat, yr) {
    data_by_yr <- dat %>%
      #'  Add year identifier to each column name
      rename_with(.cols = RN.n_bear_black:DecFeb_WSI, function(x){paste0(x, ".", yr)}) %>% 
      dplyr::select(-c(Season, Year))
    return(data_by_yr)
  }
  # RN_wide_annual <- mapply(wide_data_by_year, dat = RN_wide, yr = list("yr1", "yr2", "yr3"), SIMPLIFY = FALSE)
  RN_wide_annual <- mapply(wide_data_by_year, dat = full_bx_dat, yr = list("yr1", "yr2", "yr3"), SIMPLIFY = FALSE)
  #'  Sneak peak of each year
  head(RN_wide_annual[[1]])
  head(RN_wide_annual[[2]])
  head(RN_wide_annual[[3]])
  
  #'  Unlist as one single data frame (annual columns per species and covariates)
  RN_wide_annual_20s_22s <- full_join(RN_wide_annual[[1]], RN_wide_annual[[2]], by = c("NewLocationID", "CellID", "GMU", "Setup")) %>%
    full_join(RN_wide_annual[[3]], by = c("NewLocationID", "CellID", "GMU", "Setup")) %>%
    arrange(GMU, NewLocationID)
  head(RN_wide_annual_20s_22s)
  tail(RN_wide_annual_20s_22s)
  
  #'  Z-transform local abundance estimates (per year b/c annual estimates are stand alone variables)
  localN_z <- RN_wide_annual_20s_22s %>%
    #mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    mutate(across(starts_with(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(contains(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z) #'  Annual prey N correlated across years
  
  #'  Drop sites with NAs (missing 1+ years of data)
  localN_z <- drop_na(localN_z)
  
  #'  Make sure habitat classes are categorical factors
  localN_z <- localN_z %>%
    mutate(habitat_class.yr1 = factor(habitat_class.yr1, levels = c("Forested", "Loss_1_20", "Shrubland", "Grassland")),
           habitat_class.yr2 = factor(habitat_class.yr2, levels = c("Forested", "Loss_1_20", "Shrubland", "Grassland")),
           habitat_class.yr3 = factor(habitat_class.yr3, levels = c("Forested", "Loss_1_20", "Shrubland", "Grassland")))
  
  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dplyr::select(dat, contains("RN.n_"))
    for(i in 1:length(ndat)) {
      hist(ndat[,i])
    }
  }
  plot_histograms(localN_z)
  
  #'  Remove RN.n_ from local abundance column names
  names(localN_z) <- gsub(pattern = "RN.n*_", replacement = "", x = names(localN_z))
  #'  Remove RN. from sd column names
  names(localN_z) <- gsub(pattern = "RN._*", replacement = "", x = names(localN_z))
  
  #'  Check out a few basic relationships
  summary(lm(whitetailed_deer.yr2 ~ mountain_lion.yr1, data = localN_z, weights = precision_whitetailed_deer.yr2)) 
  summary(lm(whitetailed_deer.yr3 ~ mountain_lion.yr2, data = localN_z, weights = precision_whitetailed_deer.yr3)) 
  summary(lm(elk.yr2 ~ mountain_lion.yr1, data = localN_z, weights = precision_elk.yr2))  
  summary(lm(elk.yr3 ~ mountain_lion.yr2, data = localN_z, weights = precision_elk.yr3))  
  summary(lm(whitetailed_deer.yr2 ~ wolf.yr1, data = localN_z, weights = precision_whitetailed_deer.yr2))
  summary(lm(whitetailed_deer.yr3 ~ wolf.yr2, data = localN_z, weights = precision_whitetailed_deer.yr3))
  summary(lm(elk.yr2 ~ wolf.yr1, data = localN_z, weights = precision_elk.yr2))
  summary(lm(elk.yr3 ~ wolf.yr2, data = localN_z, weights = precision_elk.yr3))
  summary(lm(moose.yr2 ~ wolf.yr1, data = localN_z, weights = precision_moose.yr2))
  summary(lm(moose.yr3 ~ wolf.yr2, data = localN_z, weights = precision_moose.yr3))
  summary(lm(whitetailed_deer.yr2 ~ coyote.yr1, data = localN_z, weights = precision_whitetailed_deer.yr2))
  summary(lm(whitetailed_deer.yr3 ~ coyote.yr2, data = localN_z, weights = precision_whitetailed_deer.yr3))
  summary(lm(mountain_lion.yr2 ~ wolf.yr1, data = localN_z, weights = precision_mountain_lion.yr2))
  summary(lm(mountain_lion.yr3 ~ wolf.yr2, data = localN_z, weights = precision_mountain_lion.yr3))
  summary(lm(bear_black.yr2 ~ wolf.yr1, data = localN_z, weights = precision_bear_black.yr2))
  summary(lm(bear_black.yr3 ~ wolf.yr2, data = localN_z, weights = precision_bear_black.yr3))
  summary(lm(coyote.yr2 ~ wolf.yr1, data = localN_z, weights = precision_coyote.yr2))
  summary(lm(coyote.yr3 ~ wolf.yr2, data = localN_z, weights = precision_coyote.yr3))
  summary(lm(coyote.yr2 ~ mountain_lion.yr1, data = localN_z, weights = precision_coyote.yr2))  
  summary(lm(coyote.yr3 ~ mountain_lion.yr2, data = localN_z, weights = precision_coyote.yr3))  
  summary(lm(bobcat.yr2 ~ coyote.yr1, data = localN_z, weights = precision_bobcat.yr2))
  summary(lm(bobcat.yr3 ~ coyote.yr2, data = localN_z, weights = precision_bobcat.yr3))
  
  plot(whitetailed_deer.yr2 ~ mountain_lion.yr1, data = localN_z)
  plot(elk.yr3 ~ mountain_lion.yr2, data = localN_z)
  plot(whitetailed_deer.yr2 ~ wolf.yr1, data = localN_z)
  plot(whitetailed_deer.yr3 ~ wolf.yr2, data = localN_z)
  plot(moose.yr2 ~ wolf.yr1, data = localN_z)
  plot(moose.yr3 ~ wolf.yr2, data = localN_z)
  
  #'  ----------------------------------------------
  #####  Time period format: t-1 vs t across years  #####
  #'  ----------------------------------------------
  ######  Box-Cox Transformation  ######
  #'  ----------------------------
  #'  Box-Cox Transformation on grouped data to make data more normal
  #'  T minus 1 group (per species)
  for(i in 1:length(RN_stack_n[[1]])) {
    x <- pull(RN_stack_n[[1]][,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat_stack[[1]][,i] <- new_x_exact
  }
  colnames(bx_dat_stack[[1]]) <- RN_stack_colnames[[1]]
  
  #'  T group (per species)
  for(i in 1:length(RN_stack_n[[2]])) {
    x <- pull(RN_stack_n[[2]][,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat_stack[[2]][,i] <- new_x_exact
  }
  colnames(bx_dat_stack[[2]]) <- RN_stack_colnames[[2]]
  
  #'  --------------------------------------------------
  ######  Merge t-1 & t local abundance and covariates  ######
  #'  --------------------------------------------------
  #'  Grouping data across years so t-1 data affects t data, regardless of year,
  #'  increases sample size and eliminates annual variation in estimated relationships
  #'  but prohibits estimation of indirect effects (yr 1 & 2 affect yr 2 & 3... yr2 on both sides)
  full_bx_dat_stack <- mapply(add_covs_back, olddat = dat_stack_list, newdat = bx_dat_stack, datsd = RN_stack_sd, SIMPLIFY = FALSE)
  #'  Stack year 1 & year 2 data
  dat_t_minus_1 <- full_bx_dat_stack[[1]] %>%
    #'  Indicate whether observation was from first or second year in grouped data
    mutate(GroupYear = ifelse(Year == "yr1", "first", "second")) %>%
    #'  Add time period identifier to each column name
    rename_with(.cols = RN.n_bear_black:DecFeb_WSI, function(x){paste0(x, ".Tminus1")})
  #'  Stack year 2 & year 3 data
  dat_t <- full_bx_dat_stack[[2]] %>%
    #'  Indicate whether observation was from first or second year in grouped data
    mutate(GroupYear = ifelse(Year == "yr2", "first", "second")) %>%
    #'  Add time period identifier to each column name
    rename_with(.cols = RN.n_bear_black:DecFeb_WSI, function(x){paste0(x, ".T")})
  
  #'  Join t-1 and t data based on camera location
  RN_wide_1YrLag_20s_22s <- full_join(dat_t_minus_1, dat_t, by = c("NewLocationID", "CellID", "GMU", "Setup", "GroupYear"))
  
  #'  Z-transform local abundance estimates (per year b/c annual estimates are stand alone variables)
  localN_z_1YrLag <- RN_wide_1YrLag_20s_22s %>%
    #mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    mutate(across(starts_with(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(contains(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z_1YrLag) #'  Annual prey N correlated across years
  
  #'  Drop sites with NAs (missing 1+ years of data)
  localN_z_1YrLag <- drop_na(localN_z_1YrLag)
  
  #'  Make sure habitat classes are categorical factors
  localN_z_1YrLag <- localN_z_1YrLag %>%
    mutate(habitat_class.Tminus1 = factor(habitat_class.Tminus1, levels = c("Forested", "Loss_1_20", "Shrubland", "Grassland")),
           habitat_class.T = factor(habitat_class.T, levels = c("Forested", "Loss_1_20", "Shrubland", "Grassland")))
  
  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dplyr::select(dat, contains("RN.n_"))
    for(i in 1:length(ndat)) {
      hist(ndat[,i])
    }
  }
  plot_histograms(localN_z_1YrLag)
  
  #'  Remove RN.n_ from local abundance column names
  names(localN_z_1YrLag) <- gsub(pattern = "RN.n*_", replacement = "", x = names(localN_z_1YrLag))
  #'  Remove RN. from sd column names
  names(localN_z_1YrLag) <- gsub(pattern = "RN._*", replacement = "", x = names(localN_z_1YrLag))
  
  #'  Check out a few basic relationships
  summary(lm(whitetailed_deer.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag, weights = precision_whitetailed_deer.Tminus1))  
  summary(lm(elk.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag, weights = precision_elk.Tminus1))  
  summary(lm(whitetailed_deer.T ~ wolf.Tminus1, data = localN_z_1YrLag, weights = precision_whitetailed_deer.Tminus1))
  summary(lm(elk.T ~ wolf.Tminus1, data = localN_z_1YrLag, weights = precision_elk.Tminus1))
  summary(lm(moose.T ~ wolf.Tminus1, data = localN_z_1YrLag, weights = precision_moose.Tminus1))
  summary(lm(whitetailed_deer.T ~ coyote.Tminus1, data = localN_z_1YrLag, weights = precision_whitetailed_deer.Tminus1))
  summary(lm(bear_black.T ~ wolf.Tminus1, data = localN_z_1YrLag, weights = precision_bear_black.Tminus1))
  summary(lm(mountain_lion.T ~ wolf.Tminus1, data = localN_z_1YrLag, weights = precision_mountain_lion.Tminus1))
  summary(lm(coyote.T ~ wolf.Tminus1, data = localN_z_1YrLag, weights = precision_coyote.Tminus1))
  summary(lm(coyote.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag, weights = precision_coyote.Tminus1))
  summary(lm(bobcat.T ~ coyote.Tminus1, data = localN_z_1YrLag, weights = precision_bobcat.Tminus1))
  
  plot(whitetailed_deer.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag)
  plot(elk.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag)
  plot(whitetailed_deer.T ~ wolf.Tminus1, data = localN_z_1YrLag)
  plot(moose.T ~ wolf.Tminus1, data = localN_z_1YrLag)
  plot(coyote.T ~ wolf.Tminus1, data = localN_z_1YrLag)
  
  
  #'  -------------------------------------------------------
  #####  Fully combined format: All years stacked together  #####
  #'  -------------------------------------------------------
  ######  Box-Cox Transformation  ######
  #'  ----------------------------
  #'  Loop through each column and pass boxcox function to linear model of data
  #'  Why this won't work in a normal function... I don't know
  #'  YEAR 1
  for(i in 1:length(RN_all_n)) {
    x <- pull(RN_all_n[,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat_all[,i] <- new_x_exact
  }
  colnames(bx_dat_all) <- RN_all_colnames
  
  #'  ----------------------------------------------
  ######  Merge all local abundance and covariates  ######
  #'  ----------------------------------------------
  #'  Merging tranformed fully combined data back with covaraites and precision estimates
  full_bx_dat_all <- add_covs_back(RN_wide_allyrs, newdat = bx_dat_all, datsd = RN_all_sd)
  
  #'  Z-transform local abundance estimates (across years)
  localN_z_all <- full_bx_dat_all %>%
    #mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    mutate(across(starts_with(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(contains(c("RN.n_", "PercDisturbedForest", "DecFeb_WSI")))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z_all) #' Looks good
  
  #'  Make sure habitat classes are categorical factors
  localN_z_all <- localN_z_all %>%
    mutate(habitat_class = factor(habitat_class, levels = c("Forested", "Loss_1_20", "Shrubland", "Grassland")))
  
  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dplyr::select(dat, contains("RN.n_"))
    for(i in 1:length(ndat)) {
      hist(ndat[,i])
    }
  }
  plot_histograms(localN_z_all)
  
  #'  Remove RN.n_ from local abundance column names
  names(localN_z_all) <- gsub(pattern = "RN.n*_", replacement = "", x = names(localN_z_all))
  #'  Remove RN. from sd column names
  names(localN_z_all) <- gsub(pattern = "RN._*", replacement = "", x = names(localN_z_all))
  
  #'  Check out a few basic relationships
  summary(lm(whitetailed_deer ~ mountain_lion, data = localN_z_all, weights = precision_whitetailed_deer))  
  summary(lm(elk ~ mountain_lion, data = localN_z_all, weights = precision_elk))  
  summary(lm(whitetailed_deer ~ wolf, data = localN_z_all, weights = precision_whitetailed_deer))
  summary(lm(elk ~ wolf, data = localN_z_all, weights = precision_elk))
  summary(lm(moose ~ wolf, data = localN_z_all, weights = precision_moose))
  summary(lm(whitetailed_deer ~ coyote, data = localN_z_all, weights = precision_whitetailed_deer))
  summary(lm(bear_black ~ wolf, data = localN_z_all, weights = precision_bear_black))
  summary(lm(mountain_lion ~ wolf, data = localN_z_all, weights = precision_mountain_lion))
  summary(lm(coyote ~ wolf, data = localN_z_all, weights = precision_coyote))
  summary(lm(coyote ~ mountain_lion, data = localN_z_all, weights = precision_mountain_lion))
  summary(lm(bobcat ~ coyote, data = localN_z_all, weights = precision_bobcat))
  
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
  
  
  
  
