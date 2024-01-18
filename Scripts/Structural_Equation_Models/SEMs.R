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
  library(lme4)
  library(tidyverse)
  
  #'  Load RN model local abundance estimates
  load("./Outputs/Relative_Abundance/RN_model/RN_abundance.RData")
  
  #'  -----------------------------------------------------
  ####  Format local abundance estimates for SEM analyses  ####
  #'  -----------------------------------------------------
  #'  Long to wide data structure
  wide_data <- function(dat) {
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
                  values_from = "RN.n")
    return(pivot_data_wide)
  }
  RN_wide <- lapply(RN_abundance, wide_data)
  
  #'  Unlist as one single dataframe
  RN_wide_20s_22s <- do.call(rbind, RN_wide)
  #'  Sneak peak of each year
  head(RN_wide[[1]])
  head(RN_wide[[2]])
  head(RN_wide[[3]])
  
  #'  List species
  spp_list <- list("bear_black", "bobcat", "coyote", "elk", "lagomorphs", "moose", "mountain_lion", "wolf", "whitetailed_deer")
  spp_names <- as.vector(unlist(spp_list))
  
  #'  Append local abundance estimates across all years for each individual species
  spp_specific_n <- function(dat, spp) {
    local_n_est <- do.call(rbind, dat) %>%
      dplyr::select(c("GMU", "NewLocationID", "CellID", "Setup","Season", "Year", all_of(spp)))
    return(local_n_est)
  }
  spp_specific_n_list <- lapply(spp_list, spp_specific_n, dat = RN_wide)
  names(spp_specific_n_list) <- spp_names
  
  #'  Lag local abundance of subordinate predators by one year
  #'  Keep local abundance of dominant predators for first two years
  wolf_dom <- spp_specific_n_list[[8]] %>%
    filter(Year != "yr3") %>%
    #'  Column to connect causal year of dominant predator to effect year of subordinate species
    mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  lion_dom <- spp_specific_n_list[[7]] %>%
    filter(Year != "yr3") %>%
    #'  Column to connect causal year of dominant predator to effect year of subordinate species
    mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  #'  Keep local abundance of subordinate predators for second two years
  lag_community <- RN_wide_20s_22s %>%
    filter(Year != "yr1") %>%
    #'  Column to connect causal year of dominant predator to effect year of subordinate species
    mutate(Cause_effect_yrs = ifelse(Season == "Smr21", "yr1_effects_yr2", "yr2_effects_yr3"))
  
  #'  Combine dominant predator and lagged subordinate predator data
  #'  Wolves are dominant over all other species
  wolf_centric <- left_join(wolf_dom, lag_community, by = c("GMU", "NewLocationID", "CellID", "Setup", "Cause_effect_yrs")) %>%
    rename(wolf = wolf.x) %>%
    rename(Cause_season = Season.x) %>%
    rename(Cause_year = Year.x) %>%
    rename(Effect_season = Season.y) %>%
    rename(Effect_year = Year.y) %>%
    dplyr::select(-c(wolf.y, Cause_year, Effect_year)) %>%
    relocate(wolf, .before = "bear_black") %>%
    relocate(Cause_effect_yrs, .after = "Setup") %>%
    filter(!is.na(Effect_season))

  #'  Mountain lions are dominant over all other species
  lion_centric <- left_join(lion_dom, lag_community, by = c("GMU", "NewLocationID", "CellID", "Setup", "Cause_effect_yrs")) %>%
    rename(mountain_lion = mountain_lion.x) %>%
    rename(Cause_season = Season.x) %>%
    rename(Cause_year = Year.x) %>%
    rename(Effect_season = Season.y) %>%
    rename(Effect_year = Year.y) %>%
    dplyr::select(-c(mountain_lion.y, Cause_year, Effect_year)) %>%
    relocate(mountain_lion, .before = "bear_black") %>%
    relocate(Cause_effect_yrs, .after = "Setup") %>%
    filter(!is.na(Effect_season))
  
  #'  Wolves and mountain lions are both dominant over other species
  wolf_lion_centric <- left_join(wolf_dom, lion_dom, by = c("GMU", "NewLocationID", "CellID", "Setup", "Year", "Season", "Cause_effect_yrs")) %>%
    left_join(lag_community, by = c("GMU", "NewLocationID", "CellID", "Setup", "Cause_effect_yrs")) %>%
    rename(mountain_lion = mountain_lion.x) %>%
    rename(wolf = wolf.x) %>%
    rename(Cause_season = Season.x) %>%
    rename(Cause_year = Year.x) %>%
    rename(Effect_season = Season.y) %>%
    rename(Effect_year = Year.y) %>%
    dplyr::select(-c(wolf.y, mountain_lion.y, Cause_year, Effect_year)) %>%
    relocate(mountain_lion, .before = "bear_black") %>%
    relocate(wolf, .before = "mountain_lion") %>%
    relocate(Cause_effect_yrs, .after = "Setup") %>%
    filter(!is.na(Effect_season))
  
  #'  -------------------------------------------------
  ####  Run SEM based on hypothesized causal networks  ####
  #'  -------------------------------------------------
  #'  DAG 3: Wolf-centric, exploitation competition
  #'  wolf --> wtd 
  #'           wtd --> bear
  #'           wtd --> lion
  #'           wtd --> coyote
  #'                   coyote --> bobcat
  #'                   coyote --> lagomorphs
  #'                              lagomorphs --> bobcat
  #'           wtd --> bobcat 
  #'  wolf --> elk 
  #'           elk --> bear
  #'  wolf --> moose
  
  #'  DAG 1a: Wolves directly negatively affect competitors (and moose), 
  #'  which indirectly affects subordinate predators and other prey
  #'  Assume normality with linear format
  dag1a_psem <- psem(
    lm(bear_black ~ wolf, data = wolf_centric),
    lm(coyote ~ wolf, data = wolf_centric),
    lm(moose ~ wolf, data = wolf_centric),
    lm(mountain_lion ~ wolf, data = wolf_centric),
    lm(bobcat ~ mountain_lion + coyote, data = wolf_centric),
    lm(elk ~ mountain_lion + bear_black, data = wolf_centric),
    lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote, data = wolf_centric),
    lm(lagomorphs ~ bobcat + coyote, data = wolf_centric),
    data = wolf_centric
  )
  basisSet(dag1a_psem)
  dSep(dag1a_psem)
  summary(dag1a_psem)

  #'  Incorporate spatial autocorrelation for paired random and trail cameras with 
  #'  a site-level random effect
  dag1a_auto_psem <- psem(
    lmer(bear_black ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(coyote ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(moose ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(mountain_lion ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = wolf_centric),
    lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = wolf_centric),
    lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = wolf_centric),
    lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = wolf_centric),
    data = wolf_centric
  )
  summary(dag1a_auto_psem)
  AIC(dag1a_psem, dag1a_auto_psem)
  
  #'  DAG 1b: Wolves directly negatively affect competitors and ungulate prey, 
  #'  which indirectly affects subordinate predators and prey
  dag1b_psem <- psem(
    lm(bear_black ~ wolf, data = wolf_centric),
    lm(coyote ~ wolf, data = wolf_centric),
    lm(moose ~ wolf, data = wolf_centric),
    lm(mountain_lion ~ wolf, data = wolf_centric),
    lm(bobcat ~ mountain_lion + coyote, data = wolf_centric),
    lm(elk ~ mountain_lion + bear_black + wolf, data = wolf_centric),
    lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + wolf, data = wolf_centric),
    lm(lagomorphs ~ bobcat + coyote, data = wolf_centric),
    data = wolf_centric
  )
  summary(dag1b_psem)
  
  dag1b_auto_psem <- psem(
    lmer(bear_black ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(coyote ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(moose ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(mountain_lion ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = wolf_centric),
    lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = wolf_centric),
    lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = wolf_centric),
    lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = wolf_centric),
    data = wolf_centric
  )
  summary(dag1b_auto_psem)
  AIC(dag1b_psem, dag1b_auto_psem)
  
  #'  DAG 2: Wolves directly positively affect generalist predators via scavenging 
  #'  but negatively affect “specialist” competitor and ungulate prey, which 
  #'  indirectly affects subordinate predators and prey
  dag2_psem <- psem(
    lm(bear_black ~ wolf + mountain_lion, data = wolf_centric),
    lm(coyote ~ wolf + mountain_lion, data = wolf_centric),
    lm(moose ~ wolf, data = wolf_centric),
    lm(mountain_lion ~ wolf, data = wolf_centric),
    lm(bobcat ~ mountain_lion + coyote, data = wolf_centric),
    lm(elk ~ mountain_lion + bear_black + wolf, data = wolf_centric),
    lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + wolf, data = wolf_centric),
    lm(lagomorphs ~ bobcat + coyote, data = wolf_centric),
    data = wolf_centric
  )
  summary(dag2_psem)
  
  dag2_auto_psem <- psem(
    lmer(bear_black ~ wolf + mountain_lion + (1 | CellID), data = wolf_centric),
    lmer(coyote ~ wolf + mountain_lion + (1 | CellID), data = wolf_centric),
    lmer(moose ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(mountain_lion ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = wolf_centric),
    lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = wolf_centric),
    lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = wolf_centric),
    lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = wolf_centric),
    data = wolf_centric
  )
  summary(dag2_auto_psem)
  AIC(dag2_psem, dag2_auto_psem)
  
  #'  DAG 3: Wolves directly affect prey which indirectly affects predators
  dag3_psem <- psem(
    lm(whitetailed_deer ~ wolf, data = wolf_centric),
    lm(elk ~ wolf, data = wolf_centric),
    lm(moose ~ wolf, data = wolf_centric),
    lm(bear_black ~ whitetailed_deer + elk, data = wolf_centric),
    lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = wolf_centric),
    lm(coyote ~  whitetailed_deer, data = wolf_centric),
    lm(mountain_lion ~ whitetailed_deer, data = wolf_centric),
    lm(lagomorphs ~ coyote, data = wolf_centric),
    data = wolf_centric
  )
  summary(dag3_psem, .progressBar = FALSE)
  
  dag3_auto_psem <- psem(
    lmer(whitetailed_deer ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(elk ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(moose ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(bear_black ~ whitetailed_deer + elk + (1 | CellID), data = wolf_centric),
    lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = wolf_centric),
    lmer(coyote ~  whitetailed_deer + (1 | CellID), data = wolf_centric),
    lmer(mountain_lion ~ whitetailed_deer + (1 | CellID), data = wolf_centric),
    lmer(lagomorphs ~ coyote + (1 | CellID), data = wolf_centric),
    data = wolf_centric
  )
  summary(dag3_auto_psem, .progressBar = FALSE)
  AIC(dag3_psem, dag3_auto_psem)
  
  
  #'  DAG 4: Wolves directly affect closest competitor and primary prey, which 
  #'  indirectly affects other predators and secondary prey
  dag4_psem <- psem(
    lm(mountain_lion ~ wolf, data = wolf_centric),
    lm(moose ~ wolf, data = wolf_centric),
    lm(elk ~ wolf, data = wolf_centric),
    lm(whitetailed_deer ~ mountain_lion, data = wolf_centric),
    lm(bear_black ~ elk + whitetailed_deer, data = wolf_centric),
    lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = wolf_centric),
    lm(coyote ~ whitetailed_deer, data = wolf_centric),
    lm(lagomorphs ~ coyote, data = wolf_centric),
    data = wolf_centric
  )
  summary(dag4_psem)
  
  dag4_auto_psem <- psem(
    lmer(mountain_lion ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(moose ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(elk ~ wolf + (1 | CellID), data = wolf_centric),
    lmer(whitetailed_deer ~ mountain_lion + (1 | CellID), data = wolf_centric),
    lmer(bear_black ~ elk + whitetailed_deer + (1 | CellID), data = wolf_centric),
    lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = wolf_centric),
    lmer(coyote ~ whitetailed_deer + (1 | CellID), data = wolf_centric),
    lmer(lagomorphs ~ coyote + (1 | CellID), data = wolf_centric),
    data = wolf_centric
  )
  summary(dag4_auto_psem)
  AIC(dag4_psem, dag4_auto_psem)
  
  
  #'  DAG 5: Lions directly negatively affect competitors and ungulate prey, 
  #'  which indirectly affects subordinate predators and prey
  dag5_psem <- psem(
    lm(wolf ~ mountain_lion, data = wolf_centric),
    lm(bear_black ~ mountain_lion, data = wolf_centric),
    lm(bobcat ~ mountain_lion + coyote, data = wolf_centric),
    lm(moose ~ wolf, data = wolf_centric),
    lm(coyote ~ mountain_lion + wolf, data = wolf_centric),
    lm(elk ~ mountain_lion + wolf + bear_black, data = wolf_centric),
    lm(whitetailed_deer ~ mountain_lion + wolf + bear_black + coyote + bobcat, data = wolf_centric),
    lm(lagomorphs ~ coyote + bobcat, data = wolf_centric),
    data = wolf_centric
  )
  summary(dag5_psem)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #' #'  Adjust RN estimates to count data
  #' wolf_centric_round <- wolf_centric %>%
  #'   mutate(wolf = round(wolf, digits = 0),
  #'          bear_black = round(bear_black, 0),
  #'          bobcat = round(bobcat, 0),
  #'          coyote = round(coyote, 0),
  #'          elk = round(elk, 0),
  #'          lagomorphs = round(lagomorphs, 0),
  #'          moose = round(moose, 0),
  #'          mountain_lion = round(mountain_lion, 0),
  #'          whitetailed_deer = round(whitetailed_deer, 0))
  #' 
  #' dag3_rnd_count_psem <- psem(
  #'   glmer(whitetailed_deer ~ wolf + (1 | CellID), data = wolf_centric_round),
  #'   glmer(elk ~ wolf + (1 | CellID), family = poisson(link = "log"), data = wolf_centric_round),
  #'   glmer(moose ~ wolf + (1 | CellID), family = poisson(link = "log"), data = wolf_centric_round),
  #'   glmer(bear_black ~ whitetailed_deer + elk + (1 | CellID), family = poisson(link = "log"), data = wolf_centric_round),
  #'   glmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), family = poisson(link = "log"), data = wolf_centric_round),
  #'   glmer(coyote ~  whitetailed_deer + (1 | CellID), family = poisson(link = "log"), data = wolf_centric_round),
  #'   glmer(mountain_lion ~ whitetailed_deer + (1 | CellID), family = poisson(link = "log"), data = wolf_centric_round)
  #' )
  #' summary(dag1_rnd_count_psem, .progressBar = FALSE)
  ########### NOTE CONVERGING  ###################
  
  