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
  #install.packages("multcompView", repos="http://R-Forge.R-project.org")
  library(multcompView)
  library(tidyverse)
  
  #'  Load RN model local abundance estimates
  load("./Outputs/Relative_Abundance/RN_model/RN_abundance.RData")
  
  #'  Load camera site covariate data
  load("./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022.RData")
  
  #'  Load & fromat covariate data extracted from Google Earth Engine
  gfc <- read_csv("./Data/GEE outputs/GFC_annual_canopy_loss_area.csv") %>%
    mutate(LossYear = Year_add1 + 2001,
           ProportionLost_sq_m = CanopyLossArea_sq_m/BufferArea_sq_m) %>%
    dplyr::select(c(NewLocationID, LossYear, CanopyLossArea_sq_m, ProportionLost_sq_m))
    
  nlcd <- read_csv("./Data/GEE outputs/NLCD_frequencies_2019_2021.csv")
  
  #'  Reformat camera covariates
  format_cam_covs <- function(dat, season) {
    covs <- dat %>%
      #'  Define NLCD landcover classifications using https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
      mutate(NLCD_30m = ifelse(Landcover_30m2 == "Herbaceous", "Grassland", Landcover_30m2),        #14: Herbaceous
             NLCD_30m = ifelse(Landcover_30m2 == "Hay/Pasture", "Agriculture", NLCD_30m),             #15: Hay/Pasture
             NLCD_30m = ifelse(Landcover_30m2 == "Developed, Open Space", "Agriculture", NLCD_30m),           #21: Developed, Open Space
             NLCD_30m = ifelse(Landcover_30m2 == "Barren Land", "Other", NLCD_30m),             #31: Barren Land
             NLCD_30m = ifelse(Landcover_30m2 == "Deciduous Forest", "Forested", NLCD_30m),                 #41: Deciduous Forest
             NLCD_30m = ifelse(Landcover_30m2 == "Evergreen Forest", "Forested", NLCD_30m),                 #42: Evergreen Forest
             NLCD_30m = ifelse(Landcover_30m2 == "Mixed Forest", "Forested", NLCD_30m),                     #43: Mixed Forest
             NLCD_30m = ifelse(Landcover_30m2 == "Shrub/Scrub", "Shrubland", NLCD_30m),                     #52: Shrub/Scrub
             NLCD_30m = ifelse(Landcover_30m2 == "Grassland/Herbaceous", "Grassland", NLCD_30m),    #71: Grassland/Herbaceous
             NLCD_30m = ifelse(Landcover_30m2 == "Pasture/Hay", "Agriculture", NLCD_30m),             #81: Pasture/Hay
             NLCD_30m = ifelse(Landcover_30m2 == "Cultivated Crops", "Agriculture", NLCD_30m),        #82: Cultivated Crops
             NLCD_30m = ifelse(Landcover_30m2 == "Woody Wetlands", "Riparian", NLCD_30m),                   #90: Woody Wetlands
             NLCD_30m = ifelse(Landcover_30m2 == "Emergent Herbaceous Wetlands", "Riparian", NLCD_30m), #95: Emergent Herbaceous Wetlands
             Season = season) %>%
      relocate(Season, .after = "GMU")
    return(covs)
  }
  cams_eoe20s <- format_cam_covs(covariate_list[[1]], season = "Smr20")
  cams_eoe21s <- format_cam_covs(covariate_list[[2]], season = "Smr21")
  cams_eoe22s <- format_cam_covs(covariate_list[[3]], season = "Smr22")
  
  cam_covs_list <- list(cams_eoe20s, cams_eoe21s, cams_eoe22s)
  save(cam_covs_list, file = "./Data/Relative abundance data/RAI Phase 2/site_covariates_2020-2022_updated.RData")
  
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
      dplyr::select(-c(RN.sd, season)) %>%
      #'  Create column per species with their site-specific local abundance estimate
      pivot_wider(names_from = "Species",
                  values_from = "RN.n") %>%
      left_join(covs, by = c("NewLocationID", "GMU", "Season")) %>% 
      mutate(NLCD_30m = factor(NLCD_30m, levels = c("Forested", "Shrub mix", "Grassland/wetland"))) %>%
      dplyr::select(-c(Landcover_30m2, HabLayer_30m2))
    return(pivot_data_wide)
  }
  # RN_wide <- lapply(RN_abundance, wide_data)
  RN_wide <- mapply(wide_data, RN_abundance, cam_covs_list, SIMPLIFY = FALSE)
  
  #' #'  Version 1: Unlist as one single data frame (one column per species)
  #' RN_wide_20s_22s <- do.call(rbind, RN_wide)
  #' head(RN_wide_20s_22s)
  #' 
  #' #'  Version 2: Keep each species's data in list form
  #' #'  List species
  #' spp_list <- list("bear_black", "bobcat", "coyote", "elk", "lagomorphs", "moose", "mountain_lion", "whitetailed_deer", "wolf")
  #' spp_names <- as.vector(unlist(spp_list))
  #' 
  #' #'  Append local abundance estimates across all years for each individual species
  #' spp_specific_n <- function(dat, spp) {
  #'   local_n_est <- do.call(rbind, dat) %>%
  #'     dplyr::select(c("GMU", "NewLocationID", "CellID", "Setup","Season", "Year", 
  #'                     "dist2rd", "NLCD_30m", "Dist2Rural", all_of(spp)))
  #'   return(local_n_est)
  #' }
  #' spp_specific_n_list <- lapply(spp_list, spp_specific_n, dat = RN_wide)
  #' names(spp_specific_n_list) <- spp_names
  
  #'  Version 3: Wide data structure but this time one column per year for each species & covariate
  wide_data_by_year <- function(dat, yr) {
    data_by_yr <- dat %>%
      #'  Add year identifier to each column name
      rename_with(.cols = bear_black:NLCD_30m, function(x){paste0(x, ".", yr)}) %>%
      dplyr::select(-c(Season, Year))
    return(data_by_yr)
  }
  RN_wide_annual <- mapply(wide_data_by_year, dat = RN_wide, yr = list("yr1", "yr2", "yr3"), SIMPLIFY = FALSE)
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
  # localN_z <- RN_wide_20s_22s %>%
  #   mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  localN_z <- RN_wide_annual_20s_22s %>%
    mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      # dplyr::select(c("bear_black", "bobcat", "coyote", "elk", "lagomorphs", 
      #                 "moose", "mountain_lion", "whitetailed_deer", "wolf",
      #                 "dist2rd", "Dist2Rural", "perc_forest"))
      dplyr::select(c("bear_black.yr1", "bobcat.yr1", "coyote.yr1", "elk.yr1", "lagomorphs.yr1", 
                    "moose.yr1", "mountain_lion.yr1", "whitetailed_deer.yr1", "wolf.yr1",
                    "bear_black.yr2", "bobcat.yr2", "coyote.yr2", "elk.yr2", "lagomorphs.yr2", 
                    "moose.yr2", "mountain_lion.yr2", "whitetailed_deer.yr2", "wolf.yr2",
                    "bear_black.yr3", "bobcat.yr3", "coyote.yr3", "elk.yr3", "lagomorphs.yr3", 
                    "moose.yr3", "mountain_lion.yr3", "whitetailed_deer.yr3", "wolf.yr3",
                    "dist2rd.yr3", "Dist2Rural.yr3", "perc_forest.yr3"))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z) #'  Annual prey N correlated across years
  
  #'  Drop sites with NAs (missing 1+ years of data)
  localN_z <- drop_na(localN_z)
  
  #'  Visualize data
  hist(localN_z$bear_black.yr1)
  hist(localN_z$bear_black.yr2)
  hist(localN_z$bear_black.yr3)
  hist(localN_z$bobcat.yr1)
  hist(localN_z$bobcat.yr2)
  hist(localN_z$bobcat.yr3)
  hist(localN_z$coyote.yr1)
  hist(localN_z$coyote.yr2)
  hist(localN_z$coyote.yr3)
  hist(localN_z$elk.yr1)
  hist(localN_z$elk.yr2)
  hist(localN_z$elk.yr3)
  hist(localN_z$lagomorphs.yr1)
  hist(localN_z$lagomorphs.yr2)
  hist(localN_z$lagomorphs.yr3)
  hist(localN_z$moose.yr1)
  hist(localN_z$moose.yr2)
  hist(localN_z$moose.yr3)
  hist(localN_z$mountain_lion.yr1)
  hist(localN_z$mountain_lion.yr2)
  hist(localN_z$mountain_lion.yr3)
  hist(localN_z$whitetailed_deer.yr1)
  hist(localN_z$whitetailed_deer.yr2)
  hist(localN_z$whitetailed_deer.yr3)
  hist(localN_z$wolf.yr1)
  hist(localN_z$wolf.yr2)
  hist(localN_z$wolf.yr3)
  hist(RN_wide_annual_20s_22s$wolf.yr2)
  hist(log(RN_wide_annual_20s_22s$wolf.yr2))
  #'  Most histograms have major right tail... 
  #'  Should I be logging and then z-transforming the data??? 
  #'  Kinda weird b/c response variable becomes explanatory variable
  
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
  
  #'  ---------------------------------------------
  #####  Causal relationships with NO time lag  #####
  #'  ------------------------------------------
  #' ######  DAG 1a  ###### 
  #' #'  Wolves directly negatively affect competitors (and moose), which indirectly 
  #' #'  affects subordinate predators and other prey
  #' dag1a_psem <- psem(
  #'   lm(bear_black ~ wolf, data = localN_z),
  #'   lm(coyote ~ wolf, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(mountain_lion ~ wolf, data = localN_z),
  #'   lm(bobcat ~ mountain_lion + coyote, data = localN_z),
  #'   lm(elk ~ NLCD_30m + mountain_lion + bear_black, data = localN_z),
  #'   lm(whitetailed_deer ~ NLCD_30m + mountain_lion + bear_black + bobcat + coyote, data = localN_z),
  #'   lm(lagomorphs ~ NLCD_30m + bobcat + coyote, data = localN_z),
  #'   data = localN_z
  #' )
  #' basisSet(dag1a_psem)
  #' dSep(dag1a_psem)
  #' summary(dag1a_psem)
  #' 
  #' #'  Incorporate spatial autocorrelation for paired random and trail cameras with 
  #' #'  a site-level random effect
  #' dag1a_auto_psem <- psem(
  #'   lmer(bear_black ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(coyote ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
  #'   lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = localN_z),
  #'   lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag1a_auto_psem)
  #' AIC(dag1a_psem, dag1a_auto_psem)
  #' 
  #' ######  DAG 1b  ######
  #' #'  Wolves directly negatively affect competitors and ungulate prey, which 
  #' #'  indirectly affects subordinate predators and prey
  #' dag1b_psem <- psem(
  #'   lm(bear_black ~ wolf, data = localN_z),
  #'   lm(coyote ~ wolf, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(mountain_lion ~ wolf, data = localN_z),
  #'   lm(bobcat ~ mountain_lion + coyote, data = localN_z),
  #'   lm(elk ~ mountain_lion + bear_black + wolf, data = localN_z),
  #'   lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + wolf, data = localN_z),
  #'   lm(lagomorphs ~ bobcat + coyote, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag1b_psem)
  #' 
  #' dag1b_auto_psem <- psem(
  #'   lmer(bear_black ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(coyote ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
  #'   lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = localN_z),
  #'   lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag1b_auto_psem)
  #' AIC(dag1b_psem, dag1b_auto_psem)
  #' 
  #' ######  DAG 2  ######
  #' #'  Wolves directly positively affect generalist predators via scavenging 
  #' #'  but negatively affect “specialist” competitor and ungulate prey, which 
  #' #'  indirectly affects subordinate predators and prey
  #' dag2_psem <- psem(
  #'   lm(bear_black ~ wolf + mountain_lion, data = localN_z),
  #'   lm(coyote ~ wolf + mountain_lion, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(mountain_lion ~ wolf, data = localN_z),
  #'   lm(bobcat ~ mountain_lion + coyote, data = localN_z),
  #'   lm(elk ~ mountain_lion + bear_black + wolf, data = localN_z),
  #'   lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + wolf, data = localN_z),
  #'   lm(lagomorphs ~ bobcat + coyote, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag2_psem)
  #' 
  #' dag2_auto_psem <- psem(
  #'   lmer(bear_black ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(coyote ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
  #'   lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = localN_z),
  #'   lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag2_auto_psem)
  #' AIC(dag2_psem, dag2_auto_psem)
  #' 
  #' ######  DAG 3  ######
  #' #'  Wolves directly affect prey which indirectly affects predators
  #' dag3_psem <- psem(
  #'   lm(whitetailed_deer ~ wolf, data = localN_z),
  #'   lm(elk ~ wolf, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(bear_black ~ whitetailed_deer + elk, data = localN_z),
  #'   lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = localN_z),
  #'   lm(coyote ~  whitetailed_deer, data = localN_z),
  #'   lm(mountain_lion ~ whitetailed_deer, data = localN_z),
  #'   lm(lagomorphs ~ coyote, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag3_psem, .progressBar = FALSE)
  #' 
  #' dag3_auto_psem <- psem(
  #'   lmer(whitetailed_deer ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(elk ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(bear_black ~ whitetailed_deer + elk + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = localN_z),
  #'   lmer(coyote ~  whitetailed_deer + (1 | CellID), data = localN_z),
  #'   lmer(mountain_lion ~ whitetailed_deer + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag3_auto_psem, .progressBar = FALSE)
  #' AIC(dag3_psem, dag3_auto_psem)
  #' 
  #' ######  DAG 4  ######
  #' #'  Wolves directly affect closest competitor and primary prey, which indirectly 
  #' #'  affects other predators and secondary prey
  #' dag4_psem <- psem(
  #'   lm(mountain_lion ~ wolf, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(elk ~ wolf, data = localN_z),
  #'   lm(whitetailed_deer ~ mountain_lion, data = localN_z),
  #'   lm(bear_black ~ elk + whitetailed_deer, data = localN_z),
  #'   lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = localN_z),
  #'   lm(coyote ~ whitetailed_deer, data = localN_z),
  #'   lm(lagomorphs ~ coyote, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag4_psem)
  #' 
  #' dag4_auto_psem <- psem(
  #'   lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(elk ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(whitetailed_deer ~ mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(bear_black ~ elk + whitetailed_deer + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = localN_z),
  #'   lmer(coyote ~ whitetailed_deer + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag4_auto_psem)
  #' AIC(dag4_psem, dag4_auto_psem)
  #' 
  #' ######  DAG 5 & 6  ######
  #' #'  Lions have direct (positive or negative) effect on competitors (positively 
  #' #'  via scavenging, negatively via interference) which affects ungulate prey 
  #' #'  directly and indirectly through other predators
  #' dag5_psem <- psem(
  #'   lm(wolf ~ mountain_lion, data = localN_z),
  #'   lm(bear_black ~ mountain_lion, data = localN_z),
  #'   lm(bobcat ~ mountain_lion + coyote, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(coyote ~ mountain_lion + wolf, data = localN_z),
  #'   lm(elk ~ mountain_lion + wolf + bear_black, data = localN_z),
  #'   lm(whitetailed_deer ~ mountain_lion + wolf + bear_black + coyote + bobcat, data = localN_z),
  #'   lm(lagomorphs ~ coyote + bobcat, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag5_psem)
  #' 
  #' dag5_auto_psem <- psem(
  #'   lmer(wolf ~ mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(bear_black ~ mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(coyote ~ mountain_lion + wolf + (1 | CellID), data = localN_z),
  #'   lmer(elk ~ mountain_lion + wolf + bear_black + (1 | CellID), data = localN_z),
  #'   lmer(whitetailed_deer ~ mountain_lion + wolf + bear_black + coyote + bobcat + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ coyote + bobcat + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag5_auto_psem)
  #' AIC(dag5_psem, dag5_auto_psem)
  #' 
  #' ######  DAG 7  ######
  #' #'  Lions directly affect prey which indirectly affects predators
  #' dag7_psem <- psem(
  #'   lm(elk ~ mountain_lion, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(whitetailed_deer ~ mountain_lion, data = localN_z),
  #'   lm(lagomorphs ~ coyote, data = localN_z),
  #'   lm(bear_black ~ elk + whitetailed_deer, data = localN_z),
  #'   lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = localN_z),
  #'   lm(wolf ~ elk + whitetailed_deer, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag7_psem)
  #' 
  #' dag7_auto_psem <- psem(
  #'   lmer(elk ~ mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(whitetailed_deer ~ mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
  #'   lmer(bear_black ~ elk + whitetailed_deer + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = localN_z),
  #'   lmer(wolf ~ elk + whitetailed_deer + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag7_auto_psem)
  #' AIC(dag7_psem, dag7_auto_psem)
  #' 
  #' ######  DAG 8  ######
  #' #'  Apex predators directly affect prey which indirectly affects predators
  #' dag8_psem <- psem(
  #'   lm(elk ~ wolf + mountain_lion, data = localN_z),
  #'   lm(whitetailed_deer ~ wolf + mountain_lion, data = localN_z),
  #'   lm(moose ~ wolf, data = localN_z),
  #'   lm(bear_black ~ elk, data = localN_z),
  #'   lm(bobcat ~ whitetailed_deer + lagomorphs + coyote, data = localN_z),
  #'   lm(coyote ~ whitetailed_deer, data = localN_z),
  #'   lm(lagomorphs ~ coyote, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag8_psem)
  #' 
  #' dag8_auto_psem <- psem(
  #'   lmer(elk ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(whitetailed_deer ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
  #'   lmer(moose ~ wolf + (1 | CellID), data = localN_z),
  #'   lmer(bear_black ~ elk + (1 | CellID), data = localN_z),
  #'   lmer(bobcat ~ whitetailed_deer + lagomorphs + coyote + (1 | CellID), data = localN_z),
  #'   lmer(coyote ~ whitetailed_deer + (1 | CellID), data = localN_z),
  #'   lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(dag8_auto_psem)
  #' AIC(dag8_psem, dag8_auto_psem)
  #' 
  #' 
  #' #'  ---------------------------------------
  #' #####  Causal relationships with time lag  #####
  #' #'  ---------------------------------------
  #' ######  Reformat data to include 1-year time lag  ######
  #' #'  ----------------------------------------------
  #' #'  Lag local abundance of each species based on when they start in the causal relationships
  #' #'  Local abundance of wolves and lions start yr 1 when they are the dominant predator
  #' wolf <- spp_specific_n_list[[9]] %>%
  #'   filter(Year != "yr3") %>%
  #'   #'  Column to connect causal year to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  #' lion <- spp_specific_n_list[[7]] %>%
  #'   filter(Year != "yr3") %>%
  #'   #'  Column to connect causal year to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  #' #'  Local abundance of species in year 1
  #' yr1 <- RN_wide_20s_22s %>%
  #'   filter(Year == "yr1") %>%
  #'   #'  Column to connect causal year to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  #' #'  Local abundance of species being affected in year 2 by year 1 species
  #' lag_yr2 <- RN_wide_20s_22s %>%
  #'   filter(Year == "yr2") %>%
  #'   #'  Column to connect causal year to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = ifelse(Season == "Smr21", "yr1_effects_yr2", "yr2_effects_yr3"))
  #' #'  Local abundance of species being affected in year 3 by year 1 and year 2 species
  #' lag_yr3 <- RN_wide_20s_22s %>% 
  #'   filter(Year == "yr3") %>%
  #'   #'  Column to connect causal year to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = "yr2_effects_yr3")
  #' 
  #' #'  Format data frames one of two ways:
  #' #'  -----------------------------------
  #' #'  1) Merge local N for species across a 3-year time-lagged relationship
  #' #'  Species A in Year 1 affects species B in Year 2 which affect species C in Year 3
  #' #'  Rows in df follow one of two patterns:
  #' #'  Smr20 --> Smr21 --> Smr22 
  #' #'  sppA  --> sppB  --> sppC
  #' #'            sppA  --> sppB
  #' localN_lag_3yrs <- function(center_of_universe, Smr21, Smr22, spp_drop1, spp_drop2, keep1, keep2) {
  #'   #'  Merge species A yr1 (Smr20) with species B yr2 (Smr21)
  #'   sppA_smr20 <- center_of_universe[center_of_universe$Season == "Smr20",]
  #'   sppB_smr21 <- Smr21 %>% dplyr::select(-all_of(starts_with(spp_drop1)))
  #'   sppA.B <- left_join(sppA_smr20, sppB_smr21, by = c("GMU", "NewLocationID", "CellID", "Setup")) %>%
  #'     filter(!is.na(Season.y))
  #'   #'  Merge species A & B yr1 & yr2 (Smr20 & Smr21) above with species C yr3 (Smr22)
  #'   sppC_smr22 <- Smr22 %>% dplyr::select(-all_of(starts_with(spp_drop2))) 
  #'   sppA.B.C <- left_join(sppA.B, sppC_smr22, by = c("GMU", "NewLocationID", "CellID", "Setup")) %>%
  #'     dplyr::select(all_of(starts_with(keep1))) %>%
  #'     mutate(start_yr = "Smr20") %>%
  #'     relocate(start_yr, .after = "Setup")
  #'   #'  Merge species A yr2 (Smr21) with species C yr3 (Smr22)
  #'   sppA_smr21 <- center_of_universe[center_of_universe$Season == "Smr21",]
  #'   sppB_smr22 <- Smr22 %>% dplyr::select(-all_of(spp_drop1))
  #'   sppA.B.nextyr <- left_join(sppA_smr21, sppB_smr22, by = c("GMU", "NewLocationID", "CellID", "Setup")) %>%
  #'     dplyr::select(all_of(starts_with(keep2))) %>%
  #'     mutate(start_yr = "Smr21") %>%
  #'     relocate(start_yr, .after = "Setup")
  #'   #'  Merge both lagged data sets and z-transform local abundance estimates
  #'   lagged_data <- bind_rows(sppA.B.C, sppA.B.nextyr) %>%
  #'     arrange(GMU, NewLocationID, CellID) %>%
  #'     #'  Center and scale all numeric variables
  #'     mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  #'   return(lagged_data)
  #' }
  #' #'  df1) yr 1: wolf --> yr 2: bear, coyote, lion, moose --> yr 3: bobcat, elk, whitetail 
  #' df1_lagA.B.C <- localN_lag_3yrs(center_of_universe = wolf, Smr21 = lag_yr2, Smr22 = lag_yr3, 
  #'                        spp_drop1 = c("wolf", "bobcat", "elk", "whitetailed_deer", "lagomorphs"), 
  #'                        spp_drop2 = c("wolf", "bear_black", "coyote", "mountain_lion", "moose", "lagomorphs"),
  #'                        keep1 = c("GMU", "NewLocationID", "CellID", "Setup", "wolf", "bear_black", "coyote", "moose", "mountain_lion", "bobcat", "elk", "whitetailed_deer"),
  #'                        keep2 = c("GMU", "NewLocationID", "CellID", "Setup", "wolf", "bear_black", "coyote", "moose", "mountain_lion"))
  #' #' #'  Drop any row with missing data.... OH MY GOD WE LOSE A LOT OF DATA
  #' #' df1_lagA.B.C <- df1_lagA.B.C[complete.cases(df1_lagA.B.C),]
  #' 
  #' 
  #' #'  2) Merge local N for species that affect another species the following year
  #' #'  Year 1 species affects Year 2 species which affect Year 3 species but only 
  #' #'  focus on pairwise partial sequences (i.e., yr1 to yr2 or yr2 to yr3) 
  #' #'  Smr20 --> Smr21 --> Smr22 
  #' #'  sppA  --> sppB
  #' #'            sppB  -->  sppC   # Note: data for sppB here is same as sppB data above!
  #' #'            sppA  -->  sppB   # Start sequence over in Smr21, despite lacking Smr23 data for sppC 
  #' #'  sppB  --> sppC              # Start mid-sequence where sppB in Smr20 affects sppC in Smr21, despite lacking Smr19 data for sppA 
  #' localN_lag_1yr <- function(center_of_universe, Smr20, Smr21, Smr22, spp_drop1, spp_drop2, keep1, keep2, drop_cols) {
  #'   #'  Merge yr1 species that affect yr2 species
  #'   sppA_cause_Smr20 <- center_of_universe[center_of_universe$Season == "Smr20",]
  #'   sppA_cause_Smr21 <- center_of_universe[center_of_universe$Season == "Smr21",]
  #'   sppB_effect_Smr21 <- Smr21 %>% dplyr::select(-all_of(starts_with(spp_drop1)))
  #'   sppB_effect_Smr22 <- Smr22 %>% dplyr::select(-all_of(starts_with(spp_drop1)))
  #'   sppA.B_Smr20.21 <- left_join(sppA_cause_Smr20, sppB_effect_Smr21, by = c("GMU", "NewLocationID", "CellID", "Setup"))
  #'   sppA.B_Smr21.22 <- left_join(sppA_cause_Smr21, sppB_effect_Smr22, by = c("GMU", "NewLocationID", "CellID", "Setup"))
  #'   sppA.B_lag <- bind_rows(sppA.B_Smr20.21, sppA.B_Smr21.22) %>% filter(!is.na(Season.y)) %>%
  #'     mutate(start_yr = Season.x) %>%
  #'     dplyr::select(-starts_with(drop_cols)) %>%
  #'     relocate(start_yr, .after ="Setup") %>%
  #'     arrange(GMU, NewLocationID, CellID) %>%
  #'     #'  Center and scale all continuous variables
  #'     mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  #'   #'  Merge yr2 species that affect yr3 species
  #'   sppB_cause_Smr20 <- Smr20 %>% dplyr::select(-all_of(starts_with(spp_drop1))) 
  #'   sppB_cause_Smr21 <- Smr21 %>% dplyr::select(-all_of(starts_with(spp_drop1)))
  #'   sppC_effect_Smr21 <- Smr21 %>% dplyr::select(-all_of(starts_with(spp_drop2)))
  #'   sppC_effect_Smr22 <- Smr22 %>% dplyr::select(-all_of(starts_with(spp_drop2)))
  #'   sppB.C_Smr20.21 <- left_join(sppB_cause_Smr20, sppC_effect_Smr21, by = c("GMU", "NewLocationID", "CellID", "Setup"))
  #'   sppB.C_Smr21.22 <- left_join(sppB_cause_Smr21, sppC_effect_Smr22, by = c("GMU", "NewLocationID", "CellID", "Setup"))
  #'   sppB.C_lag <- bind_rows(sppB.C_Smr20.21, sppB.C_Smr21.22) %>% filter(!is.na(Season.y)) %>%
  #'     mutate(start_yr = Season.x) %>%
  #'     dplyr::select(-all_of(starts_with(drop_cols))) %>%
  #'     relocate(start_yr, .after ="Setup") %>%
  #'     mutate(bear_black_b = bear_black, coyote_b = coyote, moose_b = moose, mountain_lion_b = mountain_lion) %>%
  #'     arrange(GMU, NewLocationID, CellID) %>%
  #'     #'  Center and scale all numeric variables
  #'     mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
  #'   #'  List two lagged data sets
  #'   lagged_data <- list(sppA.B_lag, sppB.C_lag) 
  #'   return(lagged_data)
  #' }
  #' #'  df1: wolf yr 1 --> bear, coyote, lion, moose yr 2 --> bobcat, elk, whitetail yr 3
  #' df1_lagA.B_lagB.C <- localN_lag_1yr(center_of_universe = wolf, Smr20 = yr1, Smr21 = lag_yr2, Smr22 = lag_yr3, 
  #'                        spp_drop1 = c("wolf", "bobcat", "elk", "whitetailed_deer", "lagomorphs"), 
  #'                        spp_drop2 = c("wolf", "bear_black", "coyote", "mountain_lion", "moose", "lagomorphs"),
  #'                        drop_cols = c("Season.x", "Season.y", "Year.x", "Year.y", "Cause_effect_yrs.x", "Cause_effect_yrs.y"))
  #' head(df1_lagA.B_lagB.C)
  #' 
  #' 
  #' #'  Correlation matrix for all continuous covariates
  #' cov_correlation <- function(dat) {
  #'   covs <- dat %>%
  #'     dplyr::select(-c("GMU", "NewLocationID", "CellID", "Setup", "start_yr"))
  #'   cor_matrix <- cor(covs, use = "complete.obs")
  #'   return(cor_matrix)
  #' }
  #' cov_correlation(df1_lagA.B.C)
  #' cov_correlation(df1_lagA.B_lagB.C[[1]]); cov_correlation(df1_lagA.B_lagB.C[[2]])
  #' 
  #' ######  DAG 1a time lag  ######
  #' #'  DAG 1a: wolf yr 1 --> bear, coyote, lion, moose yr 2 --> bobcat, elk, whitetail yr 3
  #' #'  lag across 3 years
  #' dag1a_lag <- psem(
  #'   lmer(moose ~ wolf + (1 | CellID), data = df1_lagA.B.C),
  #'   lmer(coyote ~ wolf + (1 | CellID), data = df1_lagA.B.C),
  #'   lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = df1_lagA.B.C),
  #'   lmer(elk ~ bear_black + mountain_lion + (1 | CellID), data = df1_lagA.B.C),
  #'   lmer(whitetailed_deer ~ bear_black + mountain_lion + (1 | CellID), data = df1_lagA.B.C),
  #'   data = df1_lagA.B.C
  #' )
  #' summary(dag1a_lag)
  #' 
  #' 
  #' 
  #' 
  #' 
  #' 
  #' 
  #' #'  pairwise lag...................NOT WORKING WITH LIST DATA
  #' dag1a_lag <- psem(
  #'   lmer(elk ~ bear_black + mountain_lion + (1 | CellID), data = df1_lagA.B_lagB.C[[2]]),
  #'   lmer(whitetailed_deer ~ bear_black + mountain_lion + (1 | CellID), data = df1_lagA.B_lagB.C[[2]]),
  #'   lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = df1_lagA.B_lagB.C[[2]]),
  #'   lmer(coyote ~ wolf + (1 | CellID), data = df1_lagA.B_lagB.C[[1]]),
  #'   lmer(moose ~ wolf + (1 | CellID), data = df1_lagA.B_lagB.C[[1]]),
  #'   data = list(df1_lagA.B_lagB.C[[1]], df1_lagA.B_lagB.C[[2]])
  #' )
  #' summary(dag1a_lag)
  
  
  
  
  
  
  
  
  
  