  #'  --------------------------------
  #'  Structural Equation Models
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
  source("./Scripts/Structural_Equation_Models/Format_weather_data.R")   ########## at some point redo this so each cam location has unique avg monthly data across years (at scale of PRISM ~4600m*4600m res)
    
  #'  Generate Winter Severity Index (WSI) per year and GMU ()
  wsi <- full_join(wtr_totalPrecip, wtr_minTemp, by = c("GMU", "Season")) %>%
    #'  Multiply standardized precip by standardized temp data
    mutate(DecFeb_WSI = DecFeb_meanPPT_z * DecFeb_meanMinTemp_z) %>%
    dplyr::select(c("GMU", "Season", "DecFeb_WSI")) %>%
    #'  Filter to the winters preceding Summer 2020, 2021, and 2022 surveys
    filter(Season == "Wtr1920" | Season == "Wtr2021" | Season == "Wtr2122") %>%
    #'  Change Season column to reference annual summer sampling
    rename("PreviousWinter" = "Season") %>%
    mutate(Season = ifelse(PreviousWinter == "Wtr1920", "Smr20", "Smr22"),
           Season = ifelse(PreviousWinter == "Wtr2021", "Smr21", Season)) %>%
    relocate(Season, .after = "GMU")
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
             Habitat_class = ifelse(!is.na(DisturbanceBurn) & Habitat_class == "Loss_1_20", DisturbanceBurn, Habitat_class)) %>%
      left_join(wsi, by = c("GMU", "Season")) %>%
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
  save(cam_covs_list, file = "./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022_updated.RData")
  
  #'  Focus on only habitat class and weather variables
  sem_covs <- function(covs) {
    skinny_covs <- covs %>%
      dplyr::select(c(NewLocationID, GMU, Season, habitat_class, DecFeb_WSI))
    return(skinny_covs)
  }
  sem_covs_list <- lapply(cam_covs_list, sem_covs)
  
  #'  -----------------------------------------------------
  ####  Format local abundance estimates for SEM analyses  ####
  #'  -----------------------------------------------------
  #'  Long to wide data structure
  wide_data <- function(dat, covs) {
    pivot_data_wide <- dat %>%
      #'  Create categorical year variable
      mutate(Season = season, 
             Year = season,
             Year = ifelse(Year == "Smr20", "yr1", season),
             Year = ifelse(Year == "Smr21", "yr2", Year),
             Year = ifelse(Year == "Smr22", "yr3", Year)) %>%
      #'  Drop extra columns
      dplyr::select(-c(season)) %>%
      #'  Create column per species with their site-specific local abundance estimate
      pivot_wider(names_from = "Species",
                  values_from = c("RN.n", "RN.sd")) %>%
      left_join(covs, by = c("NewLocationID", "GMU", "Season"))
    return(pivot_data_wide)
  }
  RN_wide <- mapply(wide_data, RN_abundance, sem_covs_list, SIMPLIFY = FALSE)
  
  #'  ---------------------------
  #####  Box Cox Transform data  #####
  #'  ---------------------------
  #'  Grab just the local abundance estimates per species
  dat_n <- function(dat) {
    newdat <- dat %>%
      dplyr::select(contains("RN.n_"))
    return(newdat)
  }
  RN_wide_n <- lapply(RN_wide, dat_n)
  
  dat_sd <- function(dat) {
    newdat <- dat %>%
      dplyr::select(contains("RN.sd_"))
    return(newdat)
  }
  RN_wide_sd <- lapply(RN_wide, dat_sd)
  
  #'  Grab corresponding column names
  dat_colnames <- function(dat) {
    RN_cols <- names(dat)
    return(RN_cols)
  }
  RN_wide_colnames <- lapply(RN_wide_n, dat_colnames)

  #'  Create empty matrix to hold transformed variables
  empty_matrix <- function(dat) {
    bx_dat <- matrix(NA, nrow = nrow(dat), ncol = length(dat))
    return(bx_dat)
  }
  bx_dat <- lapply(RN_wide_n, empty_matrix)
  
  #'  Loop through each column and pass boxcox function to linear model of data
  #'  Why this won't work in a normal function... I don't know
  #'  YEAR 1
  for(i in 1:length(RN_wide_skinny[[1]])) {
    x <- pull(RN_wide_skinny[[1]][,i])
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
  for(i in 1:length(RN_wide_skinny[[2]])) {
    x <- pull(RN_wide_skinny[[2]][,i])
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
  for(i in 1:length(RN_wide_skinny[[3]])) {
    x <- pull(RN_wide_skinny[[3]][,i])
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
  
  #'  Add annual covariates and site info back to dataset   
  add_covs_back <- function(olddat, newdat, datsd) {
    covs <- dplyr::select(olddat, c("NewLocationID", "CellID", "GMU", "Setup", 
                                    "Season", "Year", "habitat_class", "DecFeb_WSI"))
    dat <- bind_cols(covs, newdat) %>%
      bind_cols(datsd) %>%
      relocate(habitat_class, .after = last_col()) %>%
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
    mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(contains(c("RN.n_", "DecFeb_WSI")))
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
  
  
  #' #'  Reformat from wide to long while retaining species.yr category per observation
  #' RN_long_annual_20s_22s <- localN_z %>% 
  #'   #'  Drop covariates
  #'   dplyr::select(-c(perc_forest.yr1, perc_forest.yr2, perc_forest.yr3,
  #'                    Elevation__10m2.yr1, Elevation__10m2.yr2, Elevation__10m2.yr3,
  #'                    Dist2Suburbs.yr1, Dist2Suburbs.yr2, Dist2Suburbs.yr3,
  #'                    Dist2Rural.yr1, Dist2Rural.yr2, Dist2Rural.yr3, 
  #'                    dist2rd.yr1, dist2rd.yr2, dist2rd.yr3, NLCD_30m.yr1, NLCD_30m.yr2, NLCD_30m.yr3)) %>%
  #'   #'  Create columns for species-year and local N estimates
  #'   pivot_longer(!c(NewLocationID, CellID, GMU, Setup), names_to = "Species.yr", values_to = "localN") %>%
  #'   #'  Strip .yr data from Species.yr column and retain just species name
  #'   mutate(Species = gsub("\\..*", "", Species.yr)) %>%
  #'   relocate(Species, .after = "Setup")
  #' 
  #' #'  Rename b/c I'm lazy
  #' localN_z <- RN_long_annual_20s_22s
  
  #'  ---------------------------------------------
  ####  SEM based on hypothesized causal networks  ####
  #'  ---------------------------------------------
  #'  Models build around several broad ecological hypotheses
  #'  1) Top-down vs bottom-up system
  #'  - food web is structured through top-down factors (predator driven)
  #'  - food web is structured through bottom-up factors (habitat driven)
  #'  - food web is structured through both top-down (predator) and bottom-up (habitat) factors
  #'  2) Primary prey vs prey availability shape predator-prey relationships
  #'  - predator-prey interactions center on predators and their primary prey species
  #'  - predator-prey interactions center on predators and the most abundant prey species
  #'  3) Interference vs exploitative competition shape predator-predator relationships
  #'  - predators affect competitors through interference competition
  #'  - predators affect competitors through exploitative competition via their
  #'      primary prey or the most abundant prey species
  #'  ----------------------------------------------
  #####  Top-down system  #####
  #'  Models center around hypothesis that predators structure the food web. Models 
  #'  include a combination of possible species interactions within the top-down framework.
  #'  ----------------------------------------------
  ######  H.td1  ######
  #'  Predators affect their primary prey species; predators affect their competitors
  #'  through interference competition
  H.td1_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ wolf.yr1 + bear_black.yr1 + mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ wolf.yr2 + bear_black.yr2 + mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ wolf.yr1 + mountain_lion.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ wolf.yr2 + mountain_lion.yr2 + coyote.yr2, data = localN_z),
    lm(bobcat.yr2 ~ coyote.yr1 + bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ coyote.yr2 + bobcat.yr2, data = localN_z),
    lm(moose.yr2 ~ wolf.yr1 + moose.yr1, data = localN_z),
    lm(moose.yr3 ~ wolf.yr2 + moose.yr2, data = localN_z),
    lm(elk.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + elk.yr1, data = localN_z),
    lm(elk.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + elk.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  # basisSet(H.td1_psem)
  # dSep(H.td1_psem)
  summary(H.td1_psem)
  
  ######  H.td2  ######
  #'  Predators affect their primary prey species; predators affect their competitors
  #'  through exploitative competition via their primary prey
  H.td2_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ elk.yr1 + mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ elk.yr2 + mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ elk.yr1 + bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ elk.yr2 + bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ lagomorphs.yr2 + coyote.yr2, data = localN_z),
    lm(moose.yr2 ~ wolf.yr1 + moose.yr1, data = localN_z),
    lm(moose.yr3 ~ wolf.yr2 + moose.yr2, data = localN_z),
    lm(elk.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + elk.yr1, data = localN_z),
    lm(elk.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + elk.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.td2_psem)
    
  ######  H.td3  ######
  #'  Predators affect the most abundant prey species; predators affect their competitors
  #'  through interference competition
  H.td3_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ wolf.yr1 + bear_black.yr1 + mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ wolf.yr2 + bear_black.yr2 + mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ wolf.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ wolf.yr2 + coyote.yr2, data = localN_z),
    lm(bobcat.yr2 ~ coyote.yr1 + bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ coyote.yr2 + bobcat.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    data = localN_z)
  summary(H.td3_psem)
  
  ######  H.td4  ######
  #'  Predators affect the most abundant prey species; predators affect their competitors
  #'  through exploitative competition via the most abundant prey (wolves dominant)
  H.td4_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ bear_black.yr1 + mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ bear_black.yr2 + mountain_lion.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ whitetailed_deer.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.td4_psem)
  
  ######  H.td5  ######
  #'  Predators affect the most abundant prey species; predators affect their competitors
  #'  through exploitative competition via the most abundant prey (lions dominant)
  H.td5_psem <- psem(
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(wolf.yr2 ~ whitetailed_deer.yr1 + wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ whitetailed_deer.yr2 + wolf.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(coyote.yr3 ~ whitetailed_deer.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.td5_psem)
  
  #'  ----------------------------------------------
  #####  Bottom-up system  #####
  #'  Models center around hypothesis that prey structure the food web. Models include
  #'  a combination of possible species interactions within the bottom-up framework.
  #'  ----------------------------------------------
  ######  H.bu1  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by their primary prey species; predators affect their 
  #'  competitors through interference competition
  H.bu1_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1 + whitetailed_deer.yr1 + bear_black.yr1 + wolf.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + whitetailed_deer.yr2 + bear_black.yr2 + wolf.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + elk.yr1 + wolf.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + elk.yr2 + wolf.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + wolf.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + wolf.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu1_psem)
  
  ######  H.bu2  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by their primary prey species; predators affect their 
  #'  competitors through exploitative competition
  H.bu2_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + elk.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + elk.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu2_psem)
  
  ######  H.bu3  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through interference competition
  H.bu3_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu3_psem)
  
  ######  H.bu4  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (wolves dominant)
  H.bu4_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu4_psem)
  
  ######  H.bu5  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (lions dominant)
  H.bu5_psem <- psem(
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + mountain_lion.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + mountain_lion.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + lagomorphs.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + lagomorphs.yr2 + coyote.yr2, data = localN_z),
    data = localN_z)
  summary(H.bu5_psem)
  
  #'  ----------------------------------------------
  #####  Top-down and bottom-up system  #####
  #'  Models center around hypothesis that prey structure the food web. Models include
  #'  a combination of possible species interactions within the bottom-up framework.
  #'  ----------------------------------------------
  ######  H.tdbu1  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect their primary prey species; predators affect their 
  #'  competitors through interference competition
  H.tdbu1_psem <- psem(        # + DecFeb_WSI.yr2
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + whitetailed_deer.yr1, data = localN_z), # + bear_black.yr1 + mountain_lion.yr1
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + mountain_lion.yr2, data = localN_z), # + bear_black.yr2 + whitetailed_deer.yr2
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr2 + moose.yr1 + whitetailed_deer.yr1, data = localN_z), # + bear_black.yr1 + wolf.yr1 + moose.yr2 + whitetailed_deer.yr2
    lm(mountain_lion.yr3 ~ wolf.yr3, data = localN_z), # + bear_black.yr2 + wolf.yr2 + mountain_lion.yr2 + 
    lm(bear_black.yr2 ~ bear_black.yr1 + whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z), # + wolf.yr1
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z), # + whitetailed_deer.yr2
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + elk.yr1 + bear_black.yr2 + bear_black.yr1, data = localN_z), # + wolf.yr1 + mountain_lion.yr1 + wolf.yr2
    lm(coyote.yr3 ~ coyote.yr2 + bear_black.yr3 + wolf.yr3, data = localN_z), # + wolf.yr2 + mountain_lion.yr2 + whitetailed_deer.yr2 + elk.yr2 + bear_black.yr2
    lm(bobcat.yr2 ~ bobcat.yr1 + coyote.yr2 + mountain_lion.yr2 + whitetailed_deer.yr1, data = localN_z), # + coyote.yr1 + mountain_lion.yr1 + lagomorphs.yr1
    lm(bobcat.yr3 ~ bobcat.yr2 + coyote.yr3 + mountain_lion.yr2 + mountain_lion.yr3 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z), # + coyote.yr2
    lm(moose.yr2 ~ moose.yr1 + wolf.yr2, data = localN_z), # + wolf.yr1 + habitat_class.yr1
    lm(moose.yr3 ~ moose.yr2 + wolf.yr3, data = localN_z), # + wolf.yr2 + habitat_class.yr2
    lm(elk.yr2 ~ elk.yr1 + bear_black.yr2 + moose.yr1 + moose.yr2 + habitat_class.yr1, data = localN_z), # + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1
    lm(elk.yr3 ~ elk.yr2 + bear_black.yr3 + moose.yr2 + moose.yr3 + habitat_class.yr2, data = localN_z), # + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + mountain_lion.yr1+ coyote.yr2, data = localN_z), # + habitat_class.yr1 + bear_black.yr1 + wolf.yr1 + moose.yr2 
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2 + moose.yr3 + wolf.yr2 + coyote.yr3, data = localN_z), # + mountain_lion.yr2 + bear_black.yr2
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + habitat_class.yr1, data = localN_z), # + coyote.yr1 + bobcat.yr1
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + habitat_class.yr2, data = localN_z), #  + coyote.yr2 + bobcat.yr2
    data = localN_z)
  
  #'  Add correlated errors 
  #'  Accounting for autoregression at different time lags than t to t+1
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% wolf.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% mountain_lion.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, bear_black.yr3 %~~% bear_black.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% bobcat.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, moose.yr3 %~~% moose.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, whitetailed_deer.yr3 %~~% whitetailed_deer.yr1)
  
  #'  Among species already hypothesized to interact but at a different time lag
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% mountain_lion.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% wolf.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, coyote.yr3 %~~% wolf.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, moose.yr3 %~~% wolf.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bear_black.yr2 %~~% habitat_class.yr2)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bear_black.yr3 %~~% habitat_class.yr3)
  # H.tdbu1_psem <- update(H.tdbu1_psem, elk.yr3 %~~% habitat_class.yr1)
  
  #'  Among species/habitats not previously hypothesize to interact but interactions are biologically plausible
  #'  Predator abundances are likely correlated with certain habitat types related
  #'  to their own preferences or prey distributions
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr2 %~~% habitat_class.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr2 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% habitat_class.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr2 %~~% habitat_class.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr2 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% habitat_class.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% habitat_class.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% elk.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr2 %~~% habitat_class.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% habitat_class.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr2 %~~% habitat_class.yr2)
  # H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% habitat_class.yr3)
  
  #'  Among species not previously hypothesized to interact and not sure what's driving this correlation structure
  H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr2 %~~% bobcat.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, mountain_lion.yr3 %~~% bobcat.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr2 %~~% wolf.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% wolf.yr2)
  # H.tdbu1_psem <- update(H.tdbu1_psem, bobcat.yr3 %~~% wolf.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, lagomorphs.yr2 %~~% moose.yr1)
  # H.tdbu1_psem <- update(H.tdbu1_psem, lagomorphs.yr3 %~~% moose.yr2)
  H.tdbu1_psem <- update(H.tdbu1_psem, whitetailed_deer.yr3 %~~% elk.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, elk.yr2 %~~% whitetailed_deer.yr1)
  H.tdbu1_psem <- update(H.tdbu1_psem, coyote.yr3 %~~% mountain_lion.yr3)
  H.tdbu1_psem <- update(H.tdbu1_psem, wolf.yr3 %~~% bear_black.yr1)
  
  #'  Add causal relationships identified by Dsep that were missing in original model 
  #'  Among species already hypothesized to interact but at a different time lag
  #'  Bear t ~ habitat t --> nixed due to convergence issues (probably too many categories) so added correlated error
  #'  Moose t ~ wolf t
  #'  Bobcat t ~ coyote t + mountain_lion t
  #'  Elk t ~ bear_black t + moose t
  #'  Coyote t ~ bear_balck t + wolf t
  #'  Mountain lion t ~ wolf t
  #'  Whitetailed_deer t ~ moose t + coyote t
  
  #'  Among species not previously hypothesized to interact but interactions are biologically plausible
  #'  Bobcat t ~ Mountain lion (t-1) + whitetailed_deer (t-1) + lagomorphs (t-1)
  #'  Bear t ~ whitetailed_deer (t-1)
  #'  Coyote t ~ whitetailed_deer (t-1) + elk (t-1) + bear (t-1)
  #'  Mountain lion t ~ moose (t-1) + whitetailed_deer (t-1)
  #'  wolf t ~ moose (t-1) + whitetailed_deer (t-1) + mountain_lion (t-1) + bear (t-1) 
  #'  Elk t ~ moose (t-1)
  #'  Whitetailed_deer t ~ wolf (t-1)
  
  out <- summary(H.tdbu1_psem, getOption("max.print")); print(out)
  out_coeffs <- out$coefficients
  out_ds <- dSep(H.tdbu1_psem)
  
  ######  H.tdbu2  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect their primary prey species; predators affect their 
  #'  competitors through exploitative competition
  H.tdbu2_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1 + elk.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2 + elk.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu2_psem)
  
  ######  H.tdbu3  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect most abundant prey species; predators affect their 
  #'  competitors through interference competition
  H.tdbu3_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + wolf.yr1 + mountain_lion.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + wolf.yr2 + mountain_lion.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1 + coyote.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2 + coyote.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu3_psem)
  
  ######  H.tdbu4  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affect most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (wolves dominant)
  H.tdbu4_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu4_psem)
  
  ######  H.tdbu5  ######
  #'  Prey affected by previous summer's forage and recent winter severity;
  #'  predators affected by most abundant prey species; predators affect their 
  #'  competitors through exploitative competition (lions dominant)
  H.tdbu5_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2, data = localN_z),
    lm(bobcat.yr2 ~ bobcat.yr1, data = localN_z),
    lm(bobcat.yr3 ~ bobcat.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + habitat_class.yr1, data = localN_z), # + DecFeb_WSI.yr2
    lm(moose.yr3 ~ moose.yr2 + habitat_class.yr2, data = localN_z), # + DecFeb_WSI.yr3
    lm(elk.yr2 ~ elk.yr1 + habitat_class.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + habitat_class.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + habitat_class.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + habitat_class.yr2, data = localN_z),
    lm(lagomorphs.yr2 ~ lagomorphs.yr1 + coyote.yr1 + bobcat.yr1 + habitat_class.yr1, data = localN_z),
    lm(lagomorphs.yr3 ~ lagomorphs.yr2 + coyote.yr2 + bobcat.yr2 + habitat_class.yr2, data = localN_z),
    data = localN_z)
  summary(H.tdbu5_psem)
  
  
  
  
  
  
  
  
  
  
  
  