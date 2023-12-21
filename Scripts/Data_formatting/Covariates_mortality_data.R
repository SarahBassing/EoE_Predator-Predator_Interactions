  #'  ------------------------------
  #'  Human-caused mortality metrics
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  December 2023
  #'  ------------------------------
  #'  Calculate human-caused mortality (harvest/removals) per GMU using BGMR data.
  #'  Metrics to consider include annual count of animals removed/GMU and the
  #'  annual average count of animals removed per sq-km/GMU. 
  #'  
  #'  Data provided by IDFG (12/15/2022): reported black bear, bobcat, mountain 
  #'  lion, and wolf mortality in GMUs 1, 6, & 10A 2018 - 2022. Coyote harvest is 
  #'  not reported in a systematic way so coyote mortality data are massive 
  #'  underestimate and likely spatially biased representation of actual level of
  #'  human-caused mortality for coyote.
  #'  ------------------------------
  
  #'  Load libraries
  library(sf) 
  library(units)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(tidyverse)
  
  #'  Load camera locations
  load("./Data/IDFG camera data/Problem cams/eoe20s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe21s_problem_cams.RData")
  load("./Data/IDFG camera data/Problem cams/eoe22s_problem_cams.RData")
  
  #'  Load mortality data
  mort_bear <- read.csv("./Data/IDFG BGMR data/Black Bear_BGMR_Data (2).csv") %>% 
    transmute(Species = Species,
              Kill.Date = as_date(Kill.Date, format = "%Y-%m-%d"),
              Game.Management.Unit = Game.Management.Unit,
              Method = Mortality.Agent, 
              Season = Harvest.Year,
              LICYEAR = NA) %>% 
    mutate(Method = str_to_title(Method)) %>%
    filter(Game.Management.Unit == "1" | Game.Management.Unit == "6" | Game.Management.Unit == "10A") %>%
    filter(Season >= "2017" & Season < "2023")
  
  mort_bob <- read.csv("./Data/IDFG BGMR data/Bobcats All Harvest.csv") %>% 
    transmute(Species = "Bobcat",
              Kill.Date = DATE_TAKEN,
              Game.Management.Unit = UNIT,
              Method = METHOD_TAKE_DESC, 
              Season = SEASON,
              LICYEAR = TLIC_Year) %>%
    filter(Game.Management.Unit == "1" | Game.Management.Unit == "6" | Game.Management.Unit == "10A") %>%
    filter(Season >= "2017" & Season < "2023") %>%
    mutate(Kill.Date = gsub(" ", "", Kill.Date),
           Kill.Date = as_date(Kill.Date, format = "%m/%d/%Y"),
           Method = str_to_title(Method))
  
  mort_coy <- read.csv("./Data/IDFG BGMR data/Coyotes Harvested by County.csv") %>%
    transmute(Species = "Coyote",
              Kill.Date = paste0(LICYEAR,"-01-01"),
              Game.Management.Unit = NA,
              Method = "Trapping",
              Season = NA,
              LICYEAR = LICYEAR,
              County = COUNYTRA) %>%
    filter(LICYEAR >= "2017" & LICYEAR < "2023") %>%
    mutate(Game.Management.Unit = ifelse(County == "BOUNDARY" | County == "BONNER", "1", Game.Management.Unit),
           Game.Management.Unit = ifelse(County == "KOOTENAI" | County == "BENEWAH"  | County == "LATAH", "6", Game.Management.Unit),
           Game.Management.Unit = ifelse(County == "IDAHO" | County == "LEWIS", "10A", Game.Management.Unit),
           Game.Management.Unit = ifelse(County == "SHOSHONE" | County == "CLEARWATER", "6, 10A", Game.Management.Unit)) %>%
    filter(County != "NEZ PERCE")
  
  mort_lion <- read.csv("./Data/IDFG BGMR data/Mountain Lion_BGMR_Data (2).csv") %>% 
    transmute(Species = Species,
              Kill.Date = Kill.Date,
              Game.Management.Unit = Game.Management.Unit,
              Method = Mortality.Agent, 
              Season = Harvest.Year,
              LICYEAR = NA) %>%
    filter(Game.Management.Unit == "1" | Game.Management.Unit == "6" | Game.Management.Unit == "10A") %>%
    filter(Season >= "2017" & Season < "2023") %>%
    mutate(Kill.Date = gsub(" ", "", Kill.Date),
           Kill.Date = as_date(Kill.Date, format = "%m/%d/%Y"),
           Method = str_to_title(Method))
  
  mort_wolf <- read.csv("./Data/IDFG BGMR data/wolf_GMU1_6_10A.csv") %>% # "Wolves Harvested by County.csv" is less detailed
    transmute(Species = Species,
              Kill.Date = Kill.Date,
              Game.Management.Unit = Game.Management.Unit,
              Method = Mortality.Agent, 
              Season = NA,
              LICYEAR = NA) %>%
    filter(Game.Management.Unit == "1" | Game.Management.Unit == "6" | Game.Management.Unit == "10A") %>%
    mutate(Kill.Date = gsub(" ", "", Kill.Date),
           Kill.Date = as_date(Kill.Date, format = "%m/%d/%Y"),
           Season = year(Kill.Date),
           Method = str_to_title(Method)) %>%
    filter(Season >= "2017" & Season < "2023")
  
  #'  List all mortality data together
  mort_predators <- list(mort_bear, mort_bob, mort_coy, mort_lion, mort_wolf) 
  
  #'  -----------------------------------------------------
  ####  Filter mortality data to relevant sources & times  ####
  #'  -----------------------------------------------------
  #'  Focus on sources of mortality that are reported in a relatively reliabl way
  legal_mortality <- function(mort) {
    print(unique(mort$Species))
    print(unique(mort$Method))
    mort <- mort %>%
      filter(Method == "General Hunt" | Method == "Controlled Hunt" | Method == "With Hounds" | 
               Method == "Incidental Hunting" | Method == "Trapping" | Method =="Incidental Trapping" | 
               Method == "Calling" | Method == "1" | Method == "3") 
    #' NOTE: 2022 bobcat data are just numbers but I assume Methods 1 & 3 represent 
    #' trapping/hunting b/c most common values in data
    return(mort)
  }
  mort_predators_legal <- lapply(mort_predators, legal_mortality)
  
  #'  Filter data to specific date range - focus on mortality that occurred one
  #'  year prior to each sampling period (e.g., fall 2020 - spring 2021 for Smr21 data)
  annual_mortality <- function(mort, start_date, end_date) {
    mort <- mort %>%
      filter(Kill.Date >= start_date) %>%
      filter(Kill.Date <= end_date)
    return(mort)
  }
  mort_pre_Smr20 <- lapply(mort_predators_legal, annual_mortality, start_date = "2019-09-16", end_date = "2020-06-30")
  mort_pre_Smr21 <- lapply(mort_predators_legal, annual_mortality, start_date = "2020-09-16", end_date = "2021-06-30")
  mort_pre_Smr22 <- lapply(mort_predators_legal, annual_mortality, start_date = "2021-09-16", end_date = "2022-06-30")
  
  
  #'  ---------------------
  ####  Mortality metrics  ####
  #'  ---------------------
  #'  1) Total mortality per year and GMU
  total_morts <- function(mort) {
    count_dead <- mort %>%
      #'  Sum number of reported harvested animals per GMU
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
  mort_preSmr21_n <- lapply(mort_pre_Smr21, total_morts) %>%
    do.call(rbind, .)
  mort_preSmr22_n <- lapply(mort_pre_Smr22, total_morts) %>%
    do.call(rbind, .)
  
  #'  Adjust coyote trapping data in GMU 6 & 10A
  #'  Trapping data recorded at COUNTY level but several counties overlap both GMUs
  #'  Splitting county-level trapping data based on proportional size of GMU 
  #'  (not sure how much of each county overlaps each GMU so this is the best we got)
  #'  and adding a portion of that county's trapping count to each GMU's total mortality count
  #'  GMU6: 5905.44 km2; GMU10A: 8527.31 km2
  proportion_gmu6 <- 5905.44/(5905.44 + 8527.31)
  proportion_gmu10a <- 8527.31/(5905.44 + 8527.31)
  coy20s_gmu_split <- mort_preSmr20_n[mort_preSmr20_n$GMU == "GMU6, 10A",] %>%
    mutate(gmu6_split = round(total_mortalities * proportion_gmu6, 0),
           gmu10a_split = round(total_mortalities * proportion_gmu10a, 0),
           gmu6_total = gmu6_split + mort_preSmr20_n$total_mortalities[mort_preSmr20_n$GMU == "GMU6" & mort_preSmr20_n$Species == "Coyote"], 
           gmu10a_total = gmu10a_split + mort_preSmr20_n$total_mortalities[mort_preSmr20_n$GMU == "GMU10A" & mort_preSmr20_n$Species == "Coyote"])
  coy21s_gmu_split <- mort_preSmr21_n[mort_preSmr21_n$GMU == "GMU6, 10A",] %>%
    mutate(gmu6_split = round(total_mortalities * proportion_gmu6, 0),
           gmu10a_split = round(total_mortalities * proportion_gmu10a, 0), 
           gmu6_total = gmu6_split + mort_preSmr21_n$total_mortalities[mort_preSmr21_n$GMU == "GMU6" & mort_preSmr21_n$Species == "Coyote"], 
           gmu10a_total = gmu10a_split + mort_preSmr21_n$total_mortalities[mort_preSmr21_n$GMU == "GMU10A" & mort_preSmr21_n$Species == "Coyote"])
  coy22s_gmu_split <- mort_preSmr22_n[mort_preSmr22_n$GMU == "GMU6, 10A",] %>%
    mutate(gmu6_split = round(total_mortalities * proportion_gmu6, 0),
           gmu10a_split = round(total_mortalities * proportion_gmu10a, 0),
           gmu6_total = gmu6_split + mort_preSmr22_n$total_mortalities[mort_preSmr22_n$GMU == "GMU6" & mort_preSmr22_n$Species == "Coyote"], 
           gmu10a_total = gmu10a_split + mort_preSmr22_n$total_mortalities[mort_preSmr22_n$GMU == "GMU10A" & mort_preSmr22_n$Species == "Coyote"])
  
  mort_preSmr20_n$total_mortalities[mort_preSmr20_n$GMU == "GMU6" & mort_preSmr20_n$Species == "Coyote"] <- coy20s_gmu_split$gmu6_total
  mort_preSmr20_n$total_mortalities[mort_preSmr20_n$GMU == "GMU10A" & mort_preSmr20_n$Species == "Coyote"] <- coy20s_gmu_split$gmu10a_total
  mort_preSmr20_n <- filter(mort_preSmr20_n, GMU != "GMU6, 10A")
  mort_preSmr21_n$total_mortalities[mort_preSmr21_n$GMU == "GMU6" & mort_preSmr21_n$Species == "Coyote"] <- coy21s_gmu_split$gmu6_total
  mort_preSmr21_n$total_mortalities[mort_preSmr21_n$GMU == "GMU10A" & mort_preSmr21_n$Species == "Coyote"] <- coy21s_gmu_split$gmu10a_total
  mort_preSmr21_n <- filter(mort_preSmr21_n, GMU != "GMU6, 10A")
  mort_preSmr22_n$total_mortalities[mort_preSmr22_n$GMU == "GMU6" & mort_preSmr22_n$Species == "Coyote"] <- coy22s_gmu_split$gmu6_total
  mort_preSmr22_n$total_mortalities[mort_preSmr22_n$GMU == "GMU10A" & mort_preSmr22_n$Species == "Coyote"] <- coy22s_gmu_split$gmu10a_total
  mort_preSmr22_n <- filter(mort_preSmr22_n, GMU != "GMU6, 10A")
  

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
  morts_per_area <- function(mort, gmu_area, season) {
    #'  Join total mortality per GMU with GMU area
    km2_morts <- full_join(mort, gmu_area, by = "GMU") %>%
      mutate(mortality_km2 = total_mortalities/gmu_area,
             mortality_km2 = round(mortality_km2, 3),
             #'  Which camera trapping dataset does this apply to (harvest occurs before camera trapping)
             harvest_period = season) %>%
      relocate(harvest_period, .after = "GMU") %>%
      arrange(GMU, Species)
    
    return(km2_morts)
  }
  (mort_preSmr20 <- morts_per_area(mort_preSmr20_n, gmu_area = eoe_area, season = "pre Smr20") %>%
      filter(GMU != "GMU1")) 
  (mort_preSmr21 <- morts_per_area(mort_preSmr21_n, gmu_area = eoe_area, season = "pre Smr21")) 
  (mort_preSmr22 <- morts_per_area(mort_preSmr22_n, gmu_area = eoe_area, season = "pre Smr22"))
  mort_preSmr <- rbind(mort_preSmr20, mort_preSmr21, mort_preSmr22)
  
  #'  Save
  save(mort_preSmr, file = "./Data/IDFG BGMR data/mort_preSmr.RData")
  write.csv(mort_preSmr, "./Data/IDFG BGMR data/mort_preSmr.csv")
  

  #'  -------------------------
  ####  Visualize and compare  ####
  #'  -------------------------
  mort_preSmr <- mutate(mort_preSmr,
                        year = ifelse(harvest_period == "pre Smr20", "2019-2020", harvest_period),
                        year = ifelse(year == "pre Smr21", "2020-2021", year),
                        year = ifelse(year == "pre Smr22", "2021-2022", year))
  harvest_summary <- ggplot(mort_preSmr, aes(x = year, y = mortality_km2, color = Species, group = Species)) +
    geom_point(position=position_dodge(.9)) + geom_line(position=position_dodge(.9)) +
    facet_wrap(~GMU) +
    theme_light() +
    theme(text = element_text(size = 13)) +
    labs(title = "Reported annual harvest prior to camera trapping season", 
         x = "Harvest season (Sept 16 - June 30)", y = "Harvest (animals/km2)")
  
  ggsave("./Outputs/Relative_Abundance/Figures/Harvest_summary_data.tiff", harvest_summary,
         units = "in", width = 10, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  #'  Issues: 
  #'  2022 bobcat data based on Methods "1" abd "3" so need further clarification on what to actually include
  #'  Coyote is based on reported trapping but many sources of mortality unaccounted for
  #'  No information about hound hunting for lions, only general hunt 
  #'  Currently not adjusted for hunter effort at all
    
  #'  Compare to mean density per GMU and species estimated from the RN model
  mean_RN_gmu1 <- read_csv("./Outputs/Relative_Abundance/RN_model/RN_density_gmu1.csv") %>% mutate(GMU = "GMU1")
  mean_RN_gmu6 <- read_csv("./Outputs/Relative_Abundance/RN_model/RN_density_gmu6.csv") %>% mutate(GMU = "GMU6")
  mean_RN_gmu10A <- read_csv("./Outputs/Relative_Abundance/RN_model/RN_density_gmu10A.csv") %>% mutate(GMU = "GMU10A")
  mean_RNden <- rbind(mean_RN_gmu1, mean_RN_gmu6, mean_RN_gmu10A) %>% arrange(Year, GMU, Species) %>%
    mutate(Species = str_to_title(Species))
  names(mean_RNden) <- c("Species", "Year", "Mean", "SD", "GMU")
  
  #'  Visualize RN density estimates over time
  RN_density_plot <- ggplot(mean_RNden, aes(x = Year, y = Mean, color = Species, group = Species)) +
    geom_point(position=position_dodge(.9)) + geom_line(position=position_dodge(.9)) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position=position_dodge(.9)) +
    facet_wrap(~GMU) +
    theme_light() +
    theme(text = element_text(size = 13)) +
    labs(title = "Estimated annual predator density from camera data", 
         x = "Sampling season (July 1 - Sept 15)", y = "Density (animals/km2)")
  
  ggsave("./Outputs/Relative_Abundance/Figures/RN_density_summary_plot.tiff", RN_density_plot,
         units = "in", width = 10, height = 6, dpi = 600, device = "tiff", compression = "lzw")
  
  #'  Patchworked together
  harvest_density_comparison_plot <- (harvest_summary + theme(legend.position = "none")) / RN_density_plot +
    plot_layout(guides = 'collect')
  
  ggsave("./Outputs/Relative_Abundance/Figures/Comparison_harvest_to_RN_density.tiff", harvest_density_comparison_plot,
         units = "in", width = 10, height = 12, dpi = 600, device = "tiff", compression = "lzw")
  
  
  #'  Next stop, read these into the Covaraite_Extract.R script
  