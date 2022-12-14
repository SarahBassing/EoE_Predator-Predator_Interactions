  #'  ------------------------------
  #'  Human-caused mortality metrics
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2022
  #'  ------------------------------
  #'  Calculate human-caused mortality (harvest/removals) per GMU using BGMR data.
  #'  Metrics to consider include annual count of animals removed/GMU and the
  #'  annual average count of animals removed per sq-km/GMU. 
  #'  
  #'  Data provided by IDFG (12/15/2022): reported black bear, bobcat, mountain 
  #'  lion, and wolf mortality in GMU1, 6, & 10A 2018 - 2022. No coyote data
  #'  available.
  #'  ------------------------------
  
  #'  Load libraries
  library(sf) 
  library(units)
  library(tidyverse)
  
  #'  Load camera locations
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe20w_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  
  #'  Load mortality data
  mort_bear <- read.csv("./Data/IDFG BGMR data/blackbear_GMU1_6_10A.csv")
  
  mort_bob <- read.csv("./Data/IDFG BGMR data/bobcat_GMU1_6_10A.csv") %>%
    transmute(Species = "Bobcat",
              Mortality.Agent = NA,
              Kill.Date = DATE_TAKEN,
              Game.Management.Unit = UNIT,
              Method = METHOD_TAKE_DESC, 
              Comments = COMMENT) %>%
    mutate(Mortality.Agent = ifelse(Method == "TRAPPING", "Trapping", Method),
           Mortality.Agent = ifelse(Mortality.Agent == "FOUND" | Mortality.Agent == "", "Other", Mortality.Agent),
           Mortality.Agent = ifelse(Mortality.Agent != "Trapping" & Mortality.Agent != "Other", "Hunting", Mortality.Agent))
  mort_lion <- read.csv("./Data/IDFG BGMR data/mountainlion_GMU1_6_10A.csv")
  mort_wolf <- read.csv("./Data/IDFG BGMR data/wolf_GMU1_6_10A.csv")
  
  #'  List all mortality data together
  mort_predators <- list(mort_bear, mort_bob, mort_lion, mort_wolf)
  
  #'  -----------------------------------------------------
  ####  Filter mortality data to relevant sources & times  ####
  #'  -----------------------------------------------------
  #'  Function to remove illegal take, road kill, other, & unknown sources of 
  #'  mortality owing to incomplete reporting.
  legal_mortality <- function(mort) {
    mort <- mort %>%
      filter(Mortality.Agent != "Road Kill") %>%
      filter(Mortality.Agent != "Unknown") %>%
      filter(Mortality.Agent != "Illegal Kill") %>%
      filter(Mortality.Agent != "Other") %>%
      filter(Mortality.Agent != "Natural Mortality") 
    print(unique(mort$Mortality.Agent))
    return(mort)
  }
  mort_predators_legal <- lapply(mort_predators, legal_mortality)
  
  #'  Filter data to specific date range - focus on mortality that occurred one
  #'  year prior to each sampling period (e.g., fall 2019 - spring 2020 for Smr20 
  #'  data; winter 2020 - fall 2020 for Wrt20; fall 2020 - spring 2021 for Smr21 data)
  annual_mortality <- function(mort, start_date, end_date) {
    mort <- mort %>%
      mutate(Kill.Date = gsub(" ", "", Kill.Date), # get rid of white spaces in date
             Kill.Date = as.Date(Kill.Date, format = "%m/%d/%Y")) %>%
      filter(Kill.Date >= start_date) %>%
      filter(Kill.Date <= end_date)
    return(mort)
  }
  mort_pre_Smr20 <- lapply(mort_predators_legal, annual_mortality, start_date = "2019-09-01", end_date = "2020-06-30")
  mort_pre_Wtr20 <- lapply(mort_predators_legal, annual_mortality, start_date = "2020-02-01", end_date = "2020-11-30")
  mort_pre_Smr21 <- lapply(mort_predators_legal, annual_mortality, start_date = "2020-09-01", end_date = "2021-06-30")
  
  
  #'  ---------------------
  ####  Mortality metrics  ####
  #'  ---------------------
  #'  1) Total mortality per year and GMU
  total_morts <- function(mort) {
    count_dead <- mort %>%
      #'  Sum number of reported mortalities for each GMU
      group_by(Game.Management.Unit) %>%
      summarize(total_mortalities = n()) %>%
      ungroup() %>% 
      #'  Extra formatting so the resulting table looks nice
      rename(GMU = Game.Management.Unit) %>%
      mutate(Species = unique(mort$Species),
             GMU = paste0("GMU", GMU)) %>%
      relocate(Species, .before = (GMU))
      
    return(count_dead)
  }
  mort_preSmr20_n <- lapply(mort_pre_Smr20, total_morts) %>%
    do.call(rbind, .)
  mort_preWtr20_n <- lapply(mort_pre_Wtr20, total_morts) %>%
    do.call(rbind, .)
  mort_preSmr21_n <- lapply(mort_pre_Smr21, total_morts) %>%
    do.call(rbind, .)
  

  #'  2) Area-weighted mortality per year and GMU
  #'  Read in shapefile of EoE GMUs
  eoe_gmus <- st_read("./Shapefiles/IDFG_Game_Management_Units/EoE_GMUs.shp")
  
  #'  Calculate area of each GMU in sq-km
  (gmu_area <- set_units(st_area(eoe_gmus), "km2"))
  eoe_area <- as.data.frame(cbind(gmu_area, eoe_gmus$NAME)) %>%
    mutate(gmu_area = round(as.numeric(gmu_area), 2),
           V2 = paste0("GMU", V2)) %>%
    rename(GMU = V2) %>%
    relocate(GMU, .before = gmu_area)
  
  #'  Calculate mortality per square-km in each GMU
  morts_per_area <- function(mort, gmu_area) {
    #'  Join total mortality per GMU with GMU area
    km2_morts <- full_join(mort, gmu_area, by = "GMU") %>%
      mutate(mortality_km2 = total_mortalities/gmu_area,
             mortality_km2 = round(mortality_km2, 3))
    
    return(km2_morts)
  }
  mort_preSmr20 <- morts_per_area(mort_preSmr20_n, gmu_area = eoe_area) 
  mort_preWtr20 <- morts_per_area(mort_preWtr20_n, gmu_area = eoe_area)
  mort_preSmr21 <- morts_per_area(mort_preSmr21_n, gmu_area = eoe_area) 
  
  #'  Save
  # save(mort_preSmr20, file = "./Data/IDFG BGMR data/mort_preSmr20.RData")
  # write.csv(mort_preSmr20, "./Data/IDFG BGMR data/mort_preSmr20.csv")
  # save(mort_preWtr20, file = "./Data/IDFG BGMR data/mort_preWtr20.RData")
  # write.csv(mort_preWtr20, "./Data/IDFG BGMR data/mort_preWtr20.csv")
  # save(mort_preSmr21, file = "./Data/IDFG BGMR data/mort_preSmr21.RData")
  # write.csv(mort_preSmr21, "./Data/IDFG BGMR data/mort_preSmr21.csv")
  
  
  #'  Next stop, read these into the Covaraite_Extract.R script
  