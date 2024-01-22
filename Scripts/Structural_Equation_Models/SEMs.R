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
  spp_list <- list("bear_black", "bobcat", "coyote", "elk", "lagomorphs", "moose", "mountain_lion", "whitetailed_deer", "wolf")
  spp_names <- as.vector(unlist(spp_list))
  
  #'  Append local abundance estimates across all years for each individual species
  spp_specific_n <- function(dat, spp) {
    local_n_est <- do.call(rbind, dat) %>%
      dplyr::select(c("GMU", "NewLocationID", "CellID", "Setup","Season", "Year", all_of(spp)))
    return(local_n_est)
  }
  spp_specific_n_list <- lapply(spp_list, spp_specific_n, dat = RN_wide)
  names(spp_specific_n_list) <- spp_names
  
  #'  Z-transform local abundance estimates
  localN_z <- RN_wide_20s_22s %>%
    mutate(zbear_black = scale(bear_black),
           zbobcat = scale(bobcat),
           zcoyote = scale(coyote),
           zelk = scale(elk),
           zlagomorphs = scale(lagomorphs),
           zmoose = scale(moose),
           zmountain_lion = scale(mountain_lion),
           zwhitetailed_deer = scale(whitetailed_deer),
           zwolf = scale(wolf))
  
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(c("zbear_black", "zbobcat", "zcoyote", "zelk", "zlagomorphs", 
                      "zmoose", "zmountain_lion", "zwhitetailed_deer", "zwolf"))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  (localN_corr <- cov_correlation(localN_z)) # Absolutely NOTHING's correlated
  
  
  
  #' #'  Lag local abundance of subordinate predators by one year
  #' #'  Keep local abundance of dominant predators for first two years
  #' wolf_dom <- spp_specific_n_list[[9]] %>%
  #'   filter(Year != "yr3") %>%
  #'   #'  Column to connect causal year of dominant predator to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  #' lion_dom <- spp_specific_n_list[[7]] %>%
  #'   filter(Year != "yr3") %>%
  #'   #'  Column to connect causal year of dominant predator to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  #' #'  Keep local abundance of subordinate predators for second two years
  #' lag_community <- RN_wide_20s_22s %>%
  #'   filter(Year != "yr1") %>%
  #'   #'  Column to connect causal year of dominant predator to effect year of subordinate species
  #'   mutate(Cause_effect_yrs = ifelse(Season == "Smr21", "yr1_effects_yr2", "yr2_effects_yr3"))
  #' 
  #' #'  Combine dominant predator and lagged subordinate predator data
  #' #'  Wolves are dominant over all other species
  #' wolf_centric <- left_join(wolf_dom, lag_community, by = c("GMU", "NewLocationID", "CellID", "Setup", "Cause_effect_yrs")) %>%
  #'   rename(wolf = wolf.x) %>%
  #'   rename(Cause_season = Season.x) %>%
  #'   rename(Cause_year = Year.x) %>%
  #'   rename(Effect_season = Season.y) %>%
  #'   rename(Effect_year = Year.y) %>%
  #'   dplyr::select(-c(wolf.y, Cause_year, Effect_year)) %>%
  #'   relocate(wolf, .before = "bear_black") %>%
  #'   relocate(Cause_effect_yrs, .after = "Setup") %>%
  #'   filter(!is.na(Effect_season))
  #' 
  #' #'  Mountain lions are dominant over all other species
  #' lion_centric <- left_join(lion_dom, lag_community, by = c("GMU", "NewLocationID", "CellID", "Setup", "Cause_effect_yrs")) %>%
  #'   rename(mountain_lion = mountain_lion.x) %>%
  #'   rename(Cause_season = Season.x) %>%
  #'   rename(Cause_year = Year.x) %>%
  #'   rename(Effect_season = Season.y) %>%
  #'   rename(Effect_year = Year.y) %>%
  #'   dplyr::select(-c(mountain_lion.y, Cause_year, Effect_year)) %>%
  #'   relocate(mountain_lion, .before = "bear_black") %>%
  #'   relocate(Cause_effect_yrs, .after = "Setup") %>%
  #'   filter(!is.na(Effect_season))
  #' 
  #' #'  Wolves and mountain lions are both dominant over other species
  #' wolf_lion_centric <- left_join(wolf_dom, lion_dom, by = c("GMU", "NewLocationID", "CellID", "Setup", "Year", "Season", "Cause_effect_yrs")) %>%
  #'   left_join(lag_community, by = c("GMU", "NewLocationID", "CellID", "Setup", "Cause_effect_yrs")) %>%
  #'   rename(mountain_lion = mountain_lion.x) %>%
  #'   rename(wolf = wolf.x) %>%
  #'   rename(Cause_season = Season.x) %>%
  #'   rename(Cause_year = Year.x) %>%
  #'   rename(Effect_season = Season.y) %>%
  #'   rename(Effect_year = Year.y) %>%
  #'   dplyr::select(-c(wolf.y, mountain_lion.y, Cause_year, Effect_year)) %>%
  #'   relocate(mountain_lion, .before = "bear_black") %>%
  #'   relocate(wolf, .before = "mountain_lion") %>%
  #'   relocate(Cause_effect_yrs, .after = "Setup") %>%
  #'   filter(!is.na(Effect_season))
  
  #'  ---------------------------------------------
  ####  SEM based on hypothesized causal networks  ####
  #'  ---------------------------------------------
  #####  Causal relationships with NO time lag  #####
  #'  ------------------------------------------
  ######  DAG 1a  ###### 
  #'  Wolves directly negatively affect competitors (and moose), which indirectly 
  #'  affects subordinate predators and other prey
  dag1a_psem <- psem(
    lm(bear_black ~ wolf, data = localN_z),
    lm(coyote ~ wolf, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(mountain_lion ~ wolf, data = localN_z),
    lm(bobcat ~ mountain_lion + coyote, data = localN_z),
    lm(elk ~ mountain_lion + bear_black, data = localN_z),
    lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote, data = localN_z),
    lm(lagomorphs ~ bobcat + coyote, data = localN_z),
    data = localN_z
  )
  basisSet(dag1a_psem)
  dSep(dag1a_psem)
  summary(dag1a_psem)

  #'  Incorporate spatial autocorrelation for paired random and trail cameras with 
  #'  a site-level random effect
  dag1a_auto_psem <- psem(
    lmer(bear_black ~ wolf + (1 | CellID), data = localN_z),
    lmer(coyote ~ wolf + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
    lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
    lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = localN_z),
    lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag1a_auto_psem)
  AIC(dag1a_psem, dag1a_auto_psem)
  
  ######  DAG 1b  ######
  #'  Wolves directly negatively affect competitors and ungulate prey, which 
  #'  indirectly affects subordinate predators and prey
  dag1b_psem <- psem(
    lm(bear_black ~ wolf, data = localN_z),
    lm(coyote ~ wolf, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(mountain_lion ~ wolf, data = localN_z),
    lm(bobcat ~ mountain_lion + coyote, data = localN_z),
    lm(elk ~ mountain_lion + bear_black + wolf, data = localN_z),
    lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + wolf, data = localN_z),
    lm(lagomorphs ~ bobcat + coyote, data = localN_z),
    data = localN_z
  )
  summary(dag1b_psem)
  
  dag1b_auto_psem <- psem(
    lmer(bear_black ~ wolf + (1 | CellID), data = localN_z),
    lmer(coyote ~ wolf + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
    lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
    lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = localN_z),
    lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag1b_auto_psem)
  AIC(dag1b_psem, dag1b_auto_psem)
  
  ######  DAG 2  ######
  #'  Wolves directly positively affect generalist predators via scavenging 
  #'  but negatively affect “specialist” competitor and ungulate prey, which 
  #'  indirectly affects subordinate predators and prey
  dag2_psem <- psem(
    lm(bear_black ~ wolf + mountain_lion, data = localN_z),
    lm(coyote ~ wolf + mountain_lion, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(mountain_lion ~ wolf, data = localN_z),
    lm(bobcat ~ mountain_lion + coyote, data = localN_z),
    lm(elk ~ mountain_lion + bear_black + wolf, data = localN_z),
    lm(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + wolf, data = localN_z),
    lm(lagomorphs ~ bobcat + coyote, data = localN_z),
    data = localN_z
  )
  summary(dag2_psem)
  
  dag2_auto_psem <- psem(
    lmer(bear_black ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
    lmer(coyote ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
    lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
    lmer(elk ~ mountain_lion + bear_black + (1 | CellID), data = localN_z),
    lmer(whitetailed_deer ~ mountain_lion + bear_black + bobcat + coyote + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ bobcat + coyote + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag2_auto_psem)
  AIC(dag2_psem, dag2_auto_psem)
  
  ######  DAG 3  ######
  #'  Wolves directly affect prey which indirectly affects predators
  dag3_psem <- psem(
    lm(whitetailed_deer ~ wolf, data = localN_z),
    lm(elk ~ wolf, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(bear_black ~ whitetailed_deer + elk, data = localN_z),
    lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = localN_z),
    lm(coyote ~  whitetailed_deer, data = localN_z),
    lm(mountain_lion ~ whitetailed_deer, data = localN_z),
    lm(lagomorphs ~ coyote, data = localN_z),
    data = localN_z
  )
  summary(dag3_psem, .progressBar = FALSE)
  
  dag3_auto_psem <- psem(
    lmer(whitetailed_deer ~ wolf + (1 | CellID), data = localN_z),
    lmer(elk ~ wolf + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(bear_black ~ whitetailed_deer + elk + (1 | CellID), data = localN_z),
    lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = localN_z),
    lmer(coyote ~  whitetailed_deer + (1 | CellID), data = localN_z),
    lmer(mountain_lion ~ whitetailed_deer + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag3_auto_psem, .progressBar = FALSE)
  AIC(dag3_psem, dag3_auto_psem)
  
  ######  DAG 4  ######
  #'  Wolves directly affect closest competitor and primary prey, which indirectly 
  #'  affects other predators and secondary prey
  dag4_psem <- psem(
    lm(mountain_lion ~ wolf, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(elk ~ wolf, data = localN_z),
    lm(whitetailed_deer ~ mountain_lion, data = localN_z),
    lm(bear_black ~ elk + whitetailed_deer, data = localN_z),
    lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = localN_z),
    lm(coyote ~ whitetailed_deer, data = localN_z),
    lm(lagomorphs ~ coyote, data = localN_z),
    data = localN_z
  )
  summary(dag4_psem)
  
  dag4_auto_psem <- psem(
    lmer(mountain_lion ~ wolf + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(elk ~ wolf + (1 | CellID), data = localN_z),
    lmer(whitetailed_deer ~ mountain_lion + (1 | CellID), data = localN_z),
    lmer(bear_black ~ elk + whitetailed_deer + (1 | CellID), data = localN_z),
    lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = localN_z),
    lmer(coyote ~ whitetailed_deer + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag4_auto_psem)
  AIC(dag4_psem, dag4_auto_psem)
  
  ######  DAG 5 & 6  ######
  #'  Lions have direct (positive or negative) effect on competitors (positively 
  #'  via scavenging, negatively via interference) which affects ungulate prey 
  #'  directly and indirectly through other predators
  dag5_psem <- psem(
    lm(wolf ~ mountain_lion, data = localN_z),
    lm(bear_black ~ mountain_lion, data = localN_z),
    lm(bobcat ~ mountain_lion + coyote, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(coyote ~ mountain_lion + wolf, data = localN_z),
    lm(elk ~ mountain_lion + wolf + bear_black, data = localN_z),
    lm(whitetailed_deer ~ mountain_lion + wolf + bear_black + coyote + bobcat, data = localN_z),
    lm(lagomorphs ~ coyote + bobcat, data = localN_z),
    data = localN_z
  )
  summary(dag5_psem)
  
  dag5_auto_psem <- psem(
    lmer(wolf ~ mountain_lion + (1 | CellID), data = localN_z),
    lmer(bear_black ~ mountain_lion + (1 | CellID), data = localN_z),
    lmer(bobcat ~ mountain_lion + coyote + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(coyote ~ mountain_lion + wolf + (1 | CellID), data = localN_z),
    lmer(elk ~ mountain_lion + wolf + bear_black + (1 | CellID), data = localN_z),
    lmer(whitetailed_deer ~ mountain_lion + wolf + bear_black + coyote + bobcat + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ coyote + bobcat + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag5_auto_psem)
  AIC(dag5_psem, dag5_auto_psem)
  
  ######  DAG 7  ######
  #'  Lions directly affect prey which indirectly affects predators
  dag7_psem <- psem(
    lm(elk ~ mountain_lion, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(whitetailed_deer ~ mountain_lion, data = localN_z),
    lm(lagomorphs ~ coyote, data = localN_z),
    lm(bear_black ~ elk + whitetailed_deer, data = localN_z),
    lm(bobcat ~ whitetailed_deer + coyote + lagomorphs, data = localN_z),
    lm(wolf ~ elk + whitetailed_deer, data = localN_z),
    data = localN_z
  )
  summary(dag7_psem)
  
  dag7_auto_psem <- psem(
    lmer(elk ~ mountain_lion + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(whitetailed_deer ~ mountain_lion + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
    lmer(bear_black ~ elk + whitetailed_deer + (1 | CellID), data = localN_z),
    lmer(bobcat ~ whitetailed_deer + coyote + lagomorphs + (1 | CellID), data = localN_z),
    lmer(wolf ~ elk + whitetailed_deer + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag7_auto_psem)
  AIC(dag7_psem, dag7_auto_psem)
  
  ######  DAG 8  ######
  #'  Apex predators directly affect prey which indirectly affects predators
  dag8_psem <- psem(
    lm(elk ~ wolf + mountain_lion, data = localN_z),
    lm(whitetailed_deer ~ wolf + mountain_lion, data = localN_z),
    lm(moose ~ wolf, data = localN_z),
    lm(bear_black ~ elk, data = localN_z),
    lm(bobcat ~ whitetailed_deer + lagomorphs + coyote, data = localN_z),
    lm(coyote ~ whitetailed_deer, data = localN_z),
    lm(lagomorphs ~ coyote, data = localN_z),
    data = localN_z
  )
  summary(dag8_psem)
  
  dag8_auto_psem <- psem(
    lmer(elk ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
    lmer(whitetailed_deer ~ wolf + mountain_lion + (1 | CellID), data = localN_z),
    lmer(moose ~ wolf + (1 | CellID), data = localN_z),
    lmer(bear_black ~ elk + (1 | CellID), data = localN_z),
    lmer(bobcat ~ whitetailed_deer + lagomorphs + coyote + (1 | CellID), data = localN_z),
    lmer(coyote ~ whitetailed_deer + (1 | CellID), data = localN_z),
    lmer(lagomorphs ~ coyote + (1 | CellID), data = localN_z),
    data = localN_z
  )
  summary(dag8_auto_psem)
  AIC(dag8_psem, dag8_auto_psem)
  
  
  #'  ---------------------------------------
  #####  Causal relationships with time lag  #####
  #'  ---------------------------------------
  #'  Lag local abundance of each species based on where they appear in the 
  #'  causal relationships
  wolf_yr1 <- spp_specific_n_list[[9]] %>%
    filter(Year != "yr3") %>%
    #'  Column to connect causal year of dominant predator to effect year of subordinate species
    mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  lion_yr1 <- spp_specific_n_list[[7]] %>%
    filter(Year != "yr3") %>%
    #'  Column to connect causal year of dominant predator to effect year of subordinate species
    mutate(Cause_effect_yrs = ifelse(Season == "Smr20", "yr1_effects_yr2", "yr2_effects_yr3"))
  #'  Keep local abundance of subordinate species for second two years
  lag_yr2 <- RN_wide_20s_22s %>%
    filter(Year != "yr1") %>%
    #'  Column to connect causal year of dominant predator to effect year of subordinate species
    mutate(Cause_effect_yrs = ifelse(Season == "Smr21", "yr1_effects_yr2", "yr2_effects_yr3"))
  #'  Keep local abundance of subordinate species for third year
  lag_yr3 <- localN_yr2 %>%
    filter(Year == "yr3") %>%
    #'  Column to connect causal year of dominant predator to effect year of subordinate species
    mutate(Cause_effect_yrs = "yr2_effects_yr3")
  
  ######  DAG 1a time lag  ######
  #'  wolf yr 1 --> bear, coyote, lion, moose yr 2 --> bobcat, elk, whitetail yr 3
  yr1 <- wolf_yr1
  yr2 <- lag_yr2 %>% dplyr::select(-c(wolf, bobcat, elk, whitetailed_deer))
  yr3 <- lag_yr3 %>% dplyr::select(-c(wolf, bear_black, coyote, mountain_lion, moose))
  timelag1 <- left_join(wolf_yr1, yr2, by = c("GMU", "NewLocationID", "CellID", "Setup", "Cause_effect_yrs"))
  
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
  
  dag1a_lag <- psem(
    lm(elk ~ bear_black + wolf, data = timelag1),
    lm(whitetailed_deer ~ bear_black + mountain_lion + wolf, data = timelag1),
    lm(bobcat ~ mountain_lion + coyote + wolf, data = timelag1),
    lm(coyote ~ wolf, data = timelag1),
    lm(moose ~ wolf, data = timelag1),
    data = timelag1
  )
  
  
  
  
  
  
  
  
  
  
  