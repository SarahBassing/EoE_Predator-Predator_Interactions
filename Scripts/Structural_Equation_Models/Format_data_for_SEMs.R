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
  
  #'  Load RN model local abundance estimates
  load("./Outputs/Relative_Abundance/RN_model/RN_abundance.RData")
  
  #'  ----------------------------------
  ####  Load and format covariate data  ####
  #'  ----------------------------------
  #'  Camera site covariate data
  load("./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022.RData") # make sure it's updated version with mtbs data
  
  #'  Google Earth Engine data sets
  #'  -----------------------------
  #'  PRISM weather data
  #'  Source script to load & format average monthly precip & temp data
  #'  Produces standardized average monthly total precipitation and standardized 
  #'  average monthly minimum temp for each winter and GMU for past 50 years
  source("./Scripts/Structural_Equation_Models/Format_weather_data.R")   
    
  #'  Generate Winter Severity Index (WSI) per year and NewLocationID
  wsi <- full_join(wtr_totalPrecip, wtr_minTemp, by = c("NewLocationID", "Season")) %>%
    #'  Multiply standardized precip by standardized temp data
    mutate(DecFeb_WSI = DecFeb_meanPPT_z * DecFeb_meanMinTemp_z) %>%
    dplyr::select(c("NewLocationID", "Season", "DecFeb_WSI")) %>%
    #'  Filter to the winters preceding Summer 2020, 2021, and 2022 surveys
    filter(Season == "Wtr1920" | Season == "Wtr2021" | Season == "Wtr2122") %>%
    #'  Change Season column to reference annual summer sampling
    rename("PreviousWinter" = "Season") %>%
    mutate(Season = ifelse(PreviousWinter == "Wtr1920", "Smr20", "Smr22"),
           Season = ifelse(PreviousWinter == "Wtr2021", "Smr21", Season)) %>%
    relocate(Season, .after = "NewLocationID")
  print(wsi)
  
  #'  Hansen's Global Forest Change dataset (area and year of canopy loss surrounding camera)
  gfc <- read_csv("./Data/GEE outputs/GFC_annual_canopy_loss_area.csv") %>% 
    #'  Make canopy loss year easier to interpret
    mutate(CanopyLossYear = Year_add1 + 2001,
           #'  Calculate proportion of canopy loss within area surrounding camera
           CanopyLossProp = CanopyLossArea_sq_m/BufferArea_sq_m,
           CanopyLossArea_sq_m = round(CanopyLossArea_sq_m, 3),
           CanopyLossProp = round(CanopyLossProp, 5)) %>%
    #'  Remove unnecessary columns and rows
    dplyr::select(c(NewLocationID, CanopyLossYear, CanopyLossArea_sq_m, CanopyLossProp)) %>%
    filter(CanopyLossArea_sq_m > 0) %>% 
    #'  Retain loss year where largest proportion of area was lost (not necessarily the most recent - want to represent the "dominant" habitat type so even if smaller more recent loss occurred area is best represented by larger, older loss)
    group_by(NewLocationID) %>%
    slice_max(order_by = CanopyLossProp, n = 1) %>% 
    ungroup() %>%
    arrange(NewLocationID)
    #' #'  Reformat so each loss year has a column identifying proportion of canopy lost that year (drop loss area) (only do if NOT slicing by year of max canopyloss)
    #' pivot_wider(!CanopyLossArea_sq_m, names_from = CanopyLossYear, values_from = CanopyLossProp) %>%
    #' transmute(NewLocationID = NewLocationID,
    #'           across(.cols = "2001":"2022",
    #'               .fns = ~ ., # Any data processing should go here
    #'               .names = "canopyloss_{.col}")) %>%
    #' arrange(NewLocationID)
  
  #'  National Landcover Database (frequency of landcover class surrounding camera)
  nlcd <- read_csv("./Data/GEE outputs/NLCD_frequencies_2019_2021.csv") %>%
    dplyr::select(c(NwLctID, landcover_2019, landcover_2021)) %>%
    mutate(landcover_2019 = str_replace_all(landcover_2019, "[[{}]]", ""),
           landcover_2021 = str_replace_all(landcover_2021, "[[{}]]", ""))
  nlcd19 <- nlcd %>% dplyr::select(c(NwLctID, landcover_2019)) %>%
    #'  Remove GMU1 cameras from this dataset
    mutate(GMU = sub("_.*", "", NwLctID)) %>%
    filter(GMU != "GMU1") %>%
    dplyr::select(-GMU) %>%
    #'  Separate wonky GEE character string grouping all landcover outputs together 
    #'  (FYI separate does something weird with extra columns so be sure to provide 
    #'  double the max number of possible landcover types). Ignore warning about missing pieces.
    separate(landcover_2019, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17", "col18", "col19", "col20", "col21", "col22"), sep = "[, ]")
  nlcd21 <- nlcd %>% dplyr::select(c(NwLctID, landcover_2021)) %>%
    separate(landcover_2021, into = c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17", "col18", "col19", "col20", "col21", "col22"), sep = "[, ]")
  
  #'  Function to reformat NLCD data extracted from GEE
  habitat_frequencies <- function(dat) {
    reformatted <- dat %>%
      #'  Convert to long format to organize columns better
      pivot_longer(!NwLctID, names_to = "col_name", values_to = "landcover_class") %>%
      #'  Split landcover class from number of pixels
      separate(landcover_class, into = c("NLCD_class", "npixels"), sep = "[=]") %>%
      #'  Drop unneeded columns and rows
      dplyr::select(-col_name) %>%
      filter(!is.na(npixels)) %>%
      #'  Reformat npixels 
      mutate(npixels = as.numeric(npixels),
             npixels = round(npixels, 2)) %>%
      #'  Calculate proportion of pixels made up by each cover type near camera
      group_by(NwLctID) %>%
      reframe(NLCD_class = NLCD_class,
              npixels = npixels,
              totalPix = sum(npixels),
              percentPix = npixels/totalPix) %>%
      ungroup() %>%
      mutate(percentPix = round(percentPix, 3)) %>%
      arrange(NwLctID) %>%
      rename("NewLocationID" = "NwLctID") %>%
      #'  Define NLCD landcover classifications using https://developers.google.com/earth-engine/datasets/catalog/USGS_NLCD_RELEASES_2021_REL_NLCD
      mutate(habitat_type = ifelse(NLCD_class == "11", "Other", NLCD_class),               #11: Open water
             habitat_type = ifelse(NLCD_class == "21", "Developed", habitat_type),         #21: Developed, open space: mixture of constructed materials, mostly vegetation in the form of lawn grasses
             habitat_type = ifelse(NLCD_class == "22", "Developed", habitat_type),         #22: Developed, low intensity: areas with a mixture of constructed materials and vegetation, most commonly single-family housing units
             habitat_type = ifelse(NLCD_class == "23", "Developed", habitat_type),         #21: Developed, medium intensity: areas with a mixture of constructed materials and vegetation, most commonly single-family housing units
             habitat_type = ifelse(NLCD_class == "24", "Developed", habitat_type),         #22: Developed, high intensity: highly developed areas where people reside or work in high numbers
             habitat_type = ifelse(NLCD_class == "31", "Other", habitat_type),             #31: Barren Land
             habitat_type = ifelse(NLCD_class == "41", "Forested", habitat_type),          #41: Deciduous Forest
             habitat_type = ifelse(NLCD_class == "42", "Forested", habitat_type),          #42: Evergreen Forest
             habitat_type = ifelse(NLCD_class == "43", "Forested", habitat_type),          #43: Mixed Forest
             habitat_type = ifelse(NLCD_class == "52", "Shrubland", habitat_type),         #52: Shrub/Scrub
             habitat_type = ifelse(NLCD_class == "71", "Grassland", habitat_type),         #71: Grassland/Herbaceous
             habitat_type = ifelse(NLCD_class == "81", "Agriculture", habitat_type),       #81: Pasture/Hay
             habitat_type = ifelse(NLCD_class == "82", "Agriculture", habitat_type),       #82: Cultivated Crops
             habitat_type = ifelse(NLCD_class == "90", "Riparian_woodland", habitat_type),      #90: Woody Wetlands
             habitat_type = ifelse(NLCD_class == "95", "Riparian_wetland", habitat_type)) %>%   #95: Emergent Herbaceous Wetlands)
      relocate(habitat_type, .after = "NewLocationID")
    
    #'  Filter data to the dominant habitat class within defined radius of camera
    #'  (dominant is defined as the landcover class comprising the largest proportion
    #'  of pixels within defined radius of camera)
    dominant_habitat <- reformatted %>%
      group_by(NewLocationID) %>%
      slice_max(order_by = percentPix, n = 1) %>%
      ungroup() %>%
      arrange(NewLocationID)

    #'  Review dominant habitat types
    print(table(dominant_habitat$habitat_type))

    #'  List both datasets together
    landcover_list <- list(dominant_habitat, reformatted)
    names(landcover_list) <- c("dominant_habitat", "percent_landcover")

    return(landcover_list)
    return()
  }
  landcover19 <- habitat_frequencies(nlcd19)
  landcover21 <- habitat_frequencies(nlcd21)
  
  #'  Add GEE data to larger covariate df
  format_cam_covs <- function(dat, landcov, season, camYr) {
    covs <- full_join(dat, landcov, by = "NewLocationID") %>%
      full_join(gfc, by = "NewLocationID") %>%
      relocate(Burn_year, .after = "percentPix") %>%
      mutate(Season = season,
             #'  Calculate number of years since burn/canopy loss
             YrsSinceBurn = camYr - as.numeric(Burn_year),
             YrsSinceBurn = ifelse(YrsSinceBurn <0, NA, YrsSinceBurn),
             YrsSinceLoss = camYr - CanopyLossYear,
             YrsSinceLoss = ifelse(YrsSinceLoss <0, NA, YrsSinceLoss),
             #'  Categorize years since burn/canopy loss following Barker et al. (2018) and Ganz et al. (2024)
             DisturbanceLoss = ifelse(YrsSinceLoss <= 20, "Loss_1_20", habitat_type),
             DisturbanceBurn = ifelse(YrsSinceBurn <= 5, "Burn_1_5", habitat_type),
             DisturbanceBurn = ifelse(YrsSinceBurn > 5 & YrsSinceBurn <=10, "Burn_6_10", DisturbanceBurn),
             DisturbanceBurn = ifelse(YrsSinceBurn > 10 & YrsSinceBurn <=15, "Burn_10_15", DisturbanceBurn),
             DisturbanceBurn = ifelse(YrsSinceBurn > 15 & YrsSinceBurn <=20, "Burn_16_20", DisturbanceBurn),
             DisturbanceBurn = ifelse(YrsSinceBurn > 20, "Burn_over20", DisturbanceBurn),
             #'  Generate a single habitat class covariate representing dominant habitat type and years since burn/canopy loss (if forested)
             Habitat_class = habitat_type,
             Habitat_class = ifelse(!is.na(DisturbanceLoss) & Habitat_class == "Forested", DisturbanceLoss, Habitat_class),
             Habitat_class = ifelse(!is.na(DisturbanceBurn) & Habitat_class == "Loss_1_20", DisturbanceBurn, Habitat_class),
             #'  Grab percentPix for sites with disturbance to indicate percent disturbed forest within last 20 years
             #'  (0 % disturbed forest if the dominant habitat class is unburned/unlogged forest or any other landcover type)
             PercDisturbedForest = ifelse(Habitat_class == "Loss_1_20" | Habitat_class == "Burn_1_5" | 
                                            Habitat_class == "Burn_6_10" | Habitat_class == "Burn_16_20", percentPix, 0)) %>%
      left_join(wsi, by = c("NewLocationID", "Season")) %>%
      relocate(Season, .after = "GMU") %>%
      filter(!is.na(GMU))
    
    #'  Review new habitat classes
    print(table(covs$Habitat_class))
    
    return(covs)
  }
  #'  Add GEE data to larger covariate dataframe; NOTE different landcover data applied to 2020 vs 2021 & 2022 data
  cams_eoe20s <- format_cam_covs(covariate_list[[1]], landcov = landcover19[[1]], camYr = 2020, season = "Smr20") 
  cams_eoe21s <- format_cam_covs(covariate_list[[2]], landcov = landcover21[[1]], camYr = 2021, season = "Smr21") 
  cams_eoe22s <- format_cam_covs(covariate_list[[3]], landcov = landcover21[[1]], camYr = 2022, season = "Smr22") 
  
  #'  List annual covariate data
  cam_covs_list <- list(cams_eoe20s, cams_eoe21s, cams_eoe22s)
    
  #'  Group habitat classes into fewer categories owing to few observations for some types of habitat
  #'  Grouping all burn categories with loss category, agriculture & riparian_wetland with grassland,
  #'  and riparian_woodland with forested
  lumped_habitat_class <- function(cov) {
    grouped_habitat <- cov %>%
      mutate(habitat_class = ifelse(Habitat_class == "Burn_1_5" | Habitat_class == "Burn_6_10" | 
                                    Habitat_class == "Burn_11_15" | Habitat_class == "Burn_16_20" | 
                                    Habitat_class == "Burn_over20", "Loss_1_20", Habitat_class),
             habitat_class = ifelse(Habitat_class == "Agriculture" | Habitat_class == "Riparian_wetland", "Grassland", habitat_class),
             habitat_class = ifelse(Habitat_class == "Riparian_woodland", "Forested", habitat_class))
    print(table(grouped_habitat$habitat_class))
    return(grouped_habitat)
  }
  cam_covs_list <- lapply(cam_covs_list, lumped_habitat_class)
  # save(cam_covs_list, file = "./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022_updated.RData")
  
  #'  Focus on only habitat class and weather variables
  sem_covs <- function(covs) {
    skinny_covs <- covs %>%
      dplyr::select(c(NewLocationID, GMU, Season, habitat_class, PercDisturbedForest, DecFeb_WSI))
    return(skinny_covs)
  }
  sem_covs_list <- lapply(cam_covs_list, sem_covs)
  
  #'  -----------------------------------------
  ####  Format local abundance estimates & SD  ####
  #'  -----------------------------------------
  #'  Long to wide data structure
  wide_data <- function(dat, covs) {
    pivot_data_wide <- dat %>%
      #'  Create categorical year variable
      mutate(Season = season, 
             Year = season,
             Year = ifelse(Year == "Smr20", "yr1", season),
             Year = ifelse(Year == "Smr21", "yr2", Year),
             Year = ifelse(Year == "Smr22", "yr3", Year),
             #'  Calculate the precision (1/var) of each estimate
             RN.precision = 1 / (RN.sd^2)) %>%
      #'  Drop extra columns
      dplyr::select(-c(season, RN.sd)) %>%
      #'  Create column per species with their site-specific local abundance & precision estimate
      pivot_wider(names_from = "Species",
                  values_from = c("RN.n", "RN.precision")) %>%
      left_join(covs, by = c("NewLocationID", "GMU", "Season"))
    return(pivot_data_wide)
  }
  RN_wide <- mapply(wide_data, RN_abundance, sem_covs_list, SIMPLIFY = FALSE)
  
  #'  Group data across years so t-1 data affects t data, regardless of year
  #'  Stack year 1 & year 2 data
  dat_t_minus_1 <- bind_rows(RN_wide[[1]], RN_wide[[2]])
  #'  Stack year 2 & year 3 data
  dat_t <- bind_rows(RN_wide[[2]], RN_wide[[3]])
  #'  List stacked data
  dat_stack_list <- list(dat_t_minus_1, dat_t)
  
  #'  Group all years together (stacking Year1, Year2, and Year3 together)
  RN_wide_allyrs <- rbind(RN_wide[[1]], RN_wide[[2]], RN_wide[[3]]) %>%
    filter(!is.na(DecFeb_WSI))
  
  #'  ----------------------------------------------------------
  ####  Format & transform data for different model structures  ####
  #'  ----------------------------------------------------------
  #'  Grab just the local abundance estimates per species
  dat_n <- function(dat) {
    newdat <- dat %>%
      dplyr::select(contains("RN.n_"))
    return(newdat)
  }
  #'  For annual data
  RN_wide_n <- lapply(RN_wide, dat_n)
  #'  For stacked data
  RN_stack_n <- lapply(dat_stack_list, dat_n)
  #'  For all data combined
  RN_all_n <- dat_n(RN_wide_allyrs)
  
  #'  Grab the precision estimates per species
  dat_sd <- function(dat) {
    newdat <- dat %>%
      dplyr::select(contains("RN.precision_")) 
    return(newdat)
  }
  RN_wide_sd <- lapply(RN_wide, dat_sd)
  RN_stack_sd <- lapply(dat_stack_list, dat_sd)
  RN_all_sd <- dat_sd(RN_wide_allyrs)
  
  #'  Grab corresponding column names
  dat_colnames <- function(dat) {
    RN_cols <- names(dat)
    return(RN_cols)
  }
  RN_wide_colnames <- lapply(RN_wide_n, dat_colnames)
  RN_stack_colnames <- lapply(RN_stack_n, dat_colnames)
  RN_all_colnames <- dat_colnames(RN_all_n)
  
  #'  Create empty matrix to hold transformed variables
  empty_matrix <- function(dat) {
    bx_dat <- matrix(NA, nrow = nrow(dat), ncol = length(dat))
    return(bx_dat)
  }
  bx_dat <- lapply(RN_wide_n, empty_matrix)
  bx_dat_stack <- lapply(RN_stack_n, empty_matrix)
  bx_dat_all <- empty_matrix(RN_all_n)
  
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
  for(i in 1:length(RN_wide_n[[1]])) {
    x <- pull(RN_wide_n[[1]][,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat[[1]][,i] <- new_x_exact
  }
  colnames(bx_dat[[1]]) <- RN_wide_colnames[[1]]
  
  #'  YEAR 2
  for(i in 1:length(RN_wide_n[[2]])) {
    x <- pull(RN_wide_n[[2]][,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat[[2]][,i] <- new_x_exact
  }
  colnames(bx_dat[[2]]) <- RN_wide_colnames[[2]]
  
  #'  YEAR 3
  for(i in 1:length(RN_wide_n[[3]])) {
    x <- pull(RN_wide_n[[3]][,i])
    b <- boxcox(lm(x ~ 1))
    #'  Extract lambda for each species
    lambda <- b$x[which.max(b$y)]
    print(lambda)
    #'  Transform variable based on exact lambda parameter
    new_x_exact <- (x^lambda-1)/lambda
    #'  Add transformed variable to the new matrix
    bx_dat[[3]][,i] <- new_x_exact
  }
  colnames(bx_dat[[3]]) <- RN_wide_colnames[[3]]
  
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
  
  
  
  
  