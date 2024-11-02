  #'  --------------------------------
  #'  Structural Equation Models
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2024
  #'  --------------------------------
  #'  Source data formatting script and run structural equation models (SEM) to
  #'  test hypotheses about how predator-prey and predator-predator interactions
  #'  influence wildlife populations in northern Idaho.
  #'  --------------------------------
  
  #'  Clean workspace
  rm(list = ls())

  library(piecewiseSEM)
  library(semEff)
  library(labelled)
  library(DiagrammeR)
  library(lme4)
  library(tidyverse)
  
  #'  Run script that formats covariate data
  source("./Scripts/Structural_Equation_Models/Format_covariate_data_for_SEMs.R") 
  
  #'  Run script that formats density data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_density_data_for_SEMs.R")
  
  #'  Take a quick look
  head(localN_z); head(localN_z_1YrLag); head(localN_z_all)
  localN_z_all <- localN_z_all %>%
    mutate(year = as.factor(year),
           year = as.numeric(year)) 
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  
  #'  -------------------------------------------------
  ####  SEM with 1-year time lag, no annual variation  ####
  #'  -------------------------------------------------
  #'  ---------------------------------------
  #####  Top down, interference competition  #####
  #'  ---------------------------------------
  top_down_inter <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_inter)
  AIC_psem(top_down_inter, AIC.type = "loglik")
  dSep(top_down_inter)

  #'  Reduced model
  #'  For all linear models at once-
  #'    1. Remove non-significant relationships that were initially hypothesized, rerun
  #'    2. Add newly significant relationships identified by d-sep test that were initially hypothesized, rerun one at a time
  #'    3. Add newly significant relationships identified by d-sep test that were not initially hypothesized but biologically plausible, rerun one at a time
  #'    4. Update unspecified correlation for newly significant relationships identified by d-sep test that are not expected to have a causal, unidirectional relationship, rerun
  top_down_inter.a <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    #'  Update with unspecified correlated errors - focus on exogenous variables 
    #'  (year T variables) that are not presumed to have a causal & unidirectional
    #'  relationship but are related by an underlying factor making them appear 
    #'  correlated (e.g., similar habitat use)
    elk.T %~~% whitetailed_deer.T,
    bear_black.T %~~% whitetailed_deer.T,
    coyote.T %~~% mountain_lion.T,
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_inter.a)
  
  
  #'  -------------------------------------------
  #####  Top down, interference, simpler system  #####
  #'  -------------------------------------------
  top_down_inter_simple <- psem(
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = localN_z_1YrLag
  )
  summary(top_down_inter_simple)
  
  #'  Reduced model
  top_down_inter_simple.a <- psem(
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_inter_simple.a)
  
  
  #'  ---------------------------------------
  #####  Top-down, exploitation competition  #####
  #'  ---------------------------------------
  top_down_exploit <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  ) 
  summary(top_down_exploit)
  
  top_down_exploit.a <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    elk.T %~~% whitetailed_deer.T,
    bear_black.T %~~% whitetailed_deer.T,
    coyote.T %~~% mountain_lion.T,
    wolf.T %~~% elk.T,
    data = density_wide_1YrLag_20s_22s
  ) 
  summary(top_down_exploit.a)
  
  
  #'  -------------------------------------------
  #####  Top down, exploitation, simpler system  #####
  #'  -------------------------------------------
  top_down_exploit_simple <- psem(
    lm(elk.T ~ elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_exploit_simple)
  
  top_down_exploit_simple.a <- psem(
    lm(elk.T ~ elk.Tminus1 + bear_black.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_exploit_simple.a) # Converges on the exact same model as the simplified top-down interference model
  
  
  #'  ----------------------------------------
  #####  Bottom up, interference competition  #####
  #'  ----------------------------------------
  bottom_up_inter <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s), 
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_inter)
  AIC_psem(bottom_up_inter, AIC.type = "loglik")
  
  #'  Reduced model
  bottom_up_inter.a <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s), 
    elk.T %~~% whitetailed_deer.T,
    coyote.T %~~% mountain_lion.T,
    whitetailed_deer.T %~~% bear_black.T,
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_inter.a)
  AIC_psem(bottom_up_inter.a, AIC.type = "loglik")
  
  
  #'  --------------------------------------------
  #####  Bottom up, interference, simpler system  #####
  #'  --------------------------------------------
  bottom_up_inter_simple <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(elk.T ~ elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s), 
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_inter_simple)
  
  bottom_up_inter_simple.a <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(elk.T ~ elk.Tminus1 + bear_black.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s), 
    elk.T %~~% wolf.T,
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_inter_simple.a)
  
  
  #'  --------------
  #####  Bottom up  #####
  #'  --------------
  bottom_up <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up)
  AIC_psem(bottom_up, AIC.type = "loglik")
  
  #'  Reduced model
  bottom_up.a <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    coyote.T %~~% mountain_lion.T, 
    elk.T %~~% whitetailed_deer.T,
    whitetailed_deer.T %~~% bear_black.T,
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up.a)
  AIC_psem(bottom_up.a, AIC.type = "loglik")
  
  
  #'  ------------------------------
  #####  Bottom up, simpler system  #####
  #'  ------------------------------
  bottom_up_simple <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(elk.T ~ elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + DisturbedForest_last20Yrs.Tminus1 + DecFeb_WSI.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_simple)
  
  bottom_up_simple.a <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(elk.T ~ elk.Tminus1 + bear_black.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_simple.a)
  
  
  #'  -------------------------------------------
  #####  Top-down AND bottom-up w/o competition  #####
  #'  -------------------------------------------
  top_down_bottom_up <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1 + DecFeb_WSI.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + DecFeb_WSI.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1 + DecFeb_WSI.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + harvest_sqKm.Tminus1 + harvest_sqKm_quad.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_bottom_up)
  
  top_down_bottom_up.a <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    elk.T %~~% whitetailed_deer.T,
    bear_black.T %~~% whitetailed_deer.T,
    coyote.T %~~% mountain_lion.T,
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_bottom_up.a)
  
  
  #'  -----------------------------
  ####  Model selection using AIC  ####
  #'  -----------------------------
  #'  Include all permutations of these models (original & reduced via dSep)
  modSelect <- AIC(top_down_inter, top_down_inter.a, top_down_inter_simple, top_down_inter_simple.a, 
                   top_down_exploit, top_down_exploit.a, top_down_exploit_simple, top_down_exploit_simple.a, 
                   bottom_up_inter, bottom_up_inter.a, bottom_up, bottom_up.a, 
                   bottom_up_simple, bottom_up_simple.a, bottom_up_inter_simple, bottom_up_inter_simple.a, 
                   top_down_bottom_up, top_down_bottom_up.a) 
  mod_names <- c("top_down_inter", "top_down_inter_reduced", "top_down_inter_simple", "top_down_inter_simple_reduced", 
                 "top_down_exploit", "top_down_exploit_reduced", "top_down_exploit_simple", "top_down_exploit_simple_reduced", 
                 "bottom_up_inter", "bottom_up_inter_reduced", "bottom_up", "bottom_up_reduced", 
                 "bottom_up_simple", "bottom_up_simple_reduced", "bottom_up_inter_simple", "bottom_up_inter_simple_reduced", 
                 "top_down_bottom_up", "top_down_bottom_up_reduced") 
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("model", "AIC", "K", "n")
  arrange(modSelect, AIC, decreasing = TRUE)
  
  #'  Just the reduced models 
  modSelect <- AIC(top_down_inter.a, top_down_inter_simple.a, top_down_exploit.a, 
                   top_down_exploit_simple.a, bottom_up_inter.a, bottom_up.a, 
                   bottom_up_simple.a, bottom_up_inter_simple.a, top_down_bottom_up.a)  
  mod_names <- c("top_down_inter_reduced", "top_down_inter_simple_reduced", "top_down_exploit_reduced", 
                 "top_down_exploit_simple_reduced", "bottom_up_inter_reduced", "bottom_up_reduced", 
                 "bottom_up_simple_reduced", "bottom_up_inter_simple_reduced", "top_down_bottom_up_reduced") 
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("model", "AIC", "K", "n")
  arrange(modSelect, AIC, decreasing = TRUE)
  
  #'  Review Fisher's C and Chi-square test statistics for each model
  LLchisq(top_down_inter_simple.a); fisherC(top_down_inter_simple.a)
  LLchisq(top_down_exploit_simple.a); fisherC(top_down_exploit_simple.a)
  LLchisq(bottom_up_inter.a); fisherC(bottom_up_inter.a)
  LLchisq(top_down_inter.a); fisherC(top_down_inter.a)
  LLchisq(top_down_exploit.a ); fisherC(top_down_exploit.a )
  
  #'  ----------------------
  ####  Dig into top model  ####
  #'  ----------------------
  #'  Dig into top model based on AIC, Fisher's C
  #'  Check for multicollinearity
  RVIF(top_down_inter_simple.a[[1]]) 
  RVIF(top_down_inter_simple.a[[2]]) 
  RVIF(top_down_inter_simple.a[[3]]) 
  RVIF(top_down_inter_simple.a[[4]]) 
  RVIF(top_down_inter_simple.a[[5]]) 
  
  #'  Run some basic model diagnostics on individuals models
  #'  This is handy: http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/
  # wtd_mod <- lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s)
  # plot(wtd_mod)
  elk_mod <- lm(elk.T ~ elk.Tminus1 + bear_black.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(elk_mod)
  moose_mod <- lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(moose_mod)
  lion_mod <- lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(lion_mod)
  bear_mod <- lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(bear_mod)
  wolf_mod <- lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(wolf_mod)
  # coy_mod <- lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s)
  # plot(coy_mod)
  
  #'  Visualize SEM (top but also next top models that are all within 10 deltaAIC)
  piecewiseSEM:::plot.psem(top_down_inter_simple.a, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"),
                           layout = "circle",
                           output = "visNetwork") 
  piecewiseSEM:::plot.psem(top_down_inter_simple.a, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  
  piecewiseSEM:::plot.psem(bottom_up_inter_simple.a, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(bottom_up_inter_simple.a, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  
  
  
  #'  Calculate direct, indirect, total, and mediator effects (SE & 95% CI) for 
  #'  all endogenous (response) variables using semEff package
  #'  https://murphymv.github.io/semEff/articles/semEff.html
  #'  First bootstrap standardized model coefficients (necessary for calculating SEs)
  #'  THIS TAKES AWHILE! 
  top_down_inter_simple.a_bootEff <- bootEff(top_down_inter_simple.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(top_down_inter_simple.a_bootEff, file = paste0("./Outputs/SEM/top_down_inter_simple.a_bootEff_", Sys.Date(), ".RData"))
  
  bottom_up_inter_simple.a_bootEff <- bootEff(bottom_up_inter_simple.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(bottom_up_inter_simple.a_bootEff, file = paste0("./Outputs/SEM/bottom_up_inter_simple.a_bootEff_", Sys.Date(), ".RData"))
  
  
  #'  Second calculate standardized effects for all casual pathways
  #'  Note, these are standardized unique effects (i.e., adjusted for multicollinearity;
  #'  i.e., semipartial correlations), allowing us to fully partition effects in the system
  #'  These tend to be smaller than the unadjusted standardized coefficients
  #'  If there are large differences btwn the two then consideration should be given
  #'  to the impact and relevance of multicollinearity in the system (check with RVIF())
  (top_down_inter_simple.a_semEff <- semEff(top_down_inter_simple.a_bootEff))
  summary(top_down_inter_simple.a_semEff)
  save(top_down_inter_simple.a_semEff, file = paste0("./Outputs/SEM/top_down_inter_simple.a_semEff_", Sys.Date(), ".RData"))
  
  (bottom_up_inter_simple.a_semEff <- semEff(bottom_up_inter_simple.a_bootEff))
  summary(bottom_up_inter_simple.a_semEff)
  save(bottom_up_inter_simple.a_semEff, file = paste0("./Outputs/SEM/bottom_up_inter_simple.a_semEff_", Sys.Date(), ".RData"))
  
  
  #'  Notes:
  #'  Many of the top-down and bottom-up hypothesized networks reduced down to 
  #'  the same refined model (e.g., top-down interference simple, top-down exploitation 
  #'  simple, and bottom-up simple converged on the same network).
  #'  A second model (bottom-up interference simple) was < 3 deltaAIC from this top set.
  #'  All non-simplified models were ~350 deltaAIC away from the top models, likely
  #'  because the WTD and coyote models added no additional inference about the 
  #'  causal relationship among species.
  #'  Habitat, WSI, and wolf harvest variables were not influential in any model.
  
  
  
  
  
  
  
  
  

  #'  ---------------------------------------------------
  ####  SEM with 1-year time lag, with annual variation  ####
  #'  ---------------------------------------------------
  #'  ---------------------------------------
  #####  Top down, interference competition  #####
  #'  ---------------------------------------
  top_down_inter_yr <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + coyote.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + coyote.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + harvest_sqKm.yr1 + harvest_sqKm_quad.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + harvest_sqKm.yr2 + harvest_sqKm_quad.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + wolf.yr1 + mountain_lion.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + wolf.yr2 + mountain_lion.yr2, data = localN_z),
    data = localN_z)
  summary(top_down_inter_yr)
  AIC_psem(top_down_inter_yr, AIC.type = "loglik")
  
  top_down_inter_yr.r <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    # lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    # lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
    moose.yr3 %~~% moose.yr1,
    moose.yr3 %~~% elk.yr1,
    elk.yr3 %~~% moose.yr1,
    elk.yr3 %~~% mountain_lion.yr1,
    coyote.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr3,
    bear_black.yr3 %~~% bear_black.yr1,
    elk.yr2 %~~% whitetailed_deer.yr2,
    mountain_lion.yr3 %~~% elk.yr1,
    data = localN_z
  )
  summary(top_down_inter_yr.r)     
  AIC_psem(top_down_inter_yr.r, AIC.type = "loglik")
  
  top_down_inter_yr <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr1, data = localN_z),
    lm(wolf.yr2 ~ moose.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr1 + mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + bear_black.yr1, data = localN_z),
    data = localN_z)
  summary(top_down_inter_yr)
  AIC_psem(top_down_inter_yr, AIC.type = "loglik")
  
  
  #'  -------------------------------------------
  #####  Top down, interference, simpler system  #####
  #'  -------------------------------------------
  top_down_inter_simple_yr <- psem(
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + harvest_sqKm.yr1 + harvest_sqKm_quad.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + harvest_sqKm.yr2 + harvest_sqKm_quad.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2, data = localN_z),
    data = localN_z
  )
  summary(top_down_inter_simple_yr)
  
  top_down_inter_simple_yr.r <- psem(
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    lm(wolf.yr2 ~ 1, data = localN_z),
    lm(wolf.yr3 ~ 1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    mountain_lion.yr3 %~~% elk.yr1,
    moose.yr3 %~~% elk.yr1,
    moose.yr3 %~~% moose.yr1,
    data = localN_z
  )
  summary(top_down_inter_simple_yr.r)
  
  
  #'  -----------------------
  #####  Top-down prey only  #####
  #'  -----------------------
  top_down_prey_yr <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + coyote.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + coyote.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2, data = localN_z),
    data = localN_z
  )
  summary(top_down_prey_yr)
  AIC_psem(top_down_prey_yr, AIC.type = "loglik")
  
  top_down_prey_yr.r <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    elk.yr2 %~~% whitetailed_deer.yr2,
    elk.yr3 %~~% moose.yr1,
    elk.yr1 %~~% moose.yr3,
    moose.yr3 %~~% moose.yr1,
    data = localN_z
  )
  summary(top_down_prey_yr.r)
  AIC_psem(top_down_prey_yr.r, AIC.type = "loglik")
  
  #'  ---------------------------------------
  #####  Top-down, exploitation competition  #####
  #'  ---------------------------------------
  top_down_exploit_yr <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + coyote.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + coyote.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + harvest_sqKm.yr1 + harvest_sqKm_quad.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + harvest_sqKm.yr2 + harvest_sqKm_quad.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + elk.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + elk.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2, data = localN_z),
    data = localN_z
  )
  summary(top_down_exploit_yr)
  AIC_psem(top_down_exploit_yr, AIC.type = "loglik")
  
  top_down_exploit_yr.r <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    lm(wolf.yr2 ~ 1, data = localN_z),
    lm(wolf.yr3 ~ 1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + wolf.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
    elk.yr2 %~~% whitetailed_deer.yr2,
    elk.yr3 %~~% moose.yr1,
    elk.yr1 %~~% moose.yr3,
    moose.yr3 %~~% moose.yr1,
    bear_black.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr3,
    mountain_lion.yr3 %~~% whitetailed_deer.yr1,
    data = localN_z
  )
  summary(top_down_exploit_yr.r)
  AIC_psem(top_down_exploit_yr.r, AIC.type = "loglik")
  
  #'  -------------------------------------------
  #####  Top-down, exploitation, simpler system  #####
  #'  -------------------------------------------
  top_down_exploit_simple_yr <- psem(
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1 + harvest_sqKm.yr1 + harvest_sqKm_quad.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + harvest_sqKm.yr2 + harvest_sqKm_quad.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + elk.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + elk.yr2, data = localN_z),
    data = localN_z
  )
  summary(top_down_exploit_simple_yr)
  
  top_down_exploit_simple_yr.r <- psem(
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    lm(wolf.yr2 ~ 1, data = localN_z),
    lm(wolf.yr3 ~ 1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    moose.yr3 %~~% elk.yr1,
    moose.yr3 %~~% moose.yr1,
    data = localN_z
  )
  summary(top_down_exploit_simple_yr.r)
  
  #'  ----------------------------------------
  #####  Bottom up, interference competition  #####
  #'  ----------------------------------------
  bottom_up_inter_yr <- psem(
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1 + elk.yr1 + whitetailed_deer.yr1 + DisturbedForest_last20Yrs.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + elk.yr2 + whitetailed_deer.yr2 + DisturbedForest_last20Yrs.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + wolf.yr1 + mountain_lion.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + wolf.yr2 + mountain_lion.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + DisturbedForest_last20Yrs.yr1 + DecFeb_WSI.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + DisturbedForest_last20Yrs.yr2 + DecFeb_WSI.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + DisturbedForest_last20Yrs.yr1 + DecFeb_WSI.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + DisturbedForest_last20Yrs.yr2 + DecFeb_WSI.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + DisturbedForest_last20Yrs.yr1 + DecFeb_WSI.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + DisturbedForest_last20Yrs.yr2 + DecFeb_WSI.yr2, data = localN_z),
    data = localN_z
  )
  summary(bottom_up_inter_yr)
  AIC_psem(bottom_up_inter_yr, AIC.type = "loglik")
  
  bottom_up_inter_yr.r <- psem(
    lm(wolf.yr2 ~ moose.yr1, data = localN_z),
    lm(wolf.yr3 ~ 1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + DisturbedForest_last20Yrs.yr2, data = localN_z),
    moose.yr3 %~~% moose.yr1,
    # moose.yr3 %~~% elk.yr1,
    # bear_black.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr3,
    # mountain_lion.yr3 %~~% whitetailed_deer.yr1,
    mountain_lion.yr2 %~~% wolf.yr2,
    whitetailed_deer.yr2 %~~% elk.yr2,
    data = localN_z
  )
  summary(bottom_up_inter_yr.r)
  AIC_psem(bottom_up_inter_yr.r, AIC.type = "loglik")
  
  
  #'  --------------------------------------------
  #####  Bottom up, interference, simpler system  #####
  #'  --------------------------------------------
  bottom_up_inter_simple_yr <- psem(
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1 + elk.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1 + elk.yr1 + DisturbedForest_last20Yrs.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + elk.yr2 + DisturbedForest_last20Yrs.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + DisturbedForest_last20Yrs.yr1 + DecFeb_WSI.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + DisturbedForest_last20Yrs.yr2 + DecFeb_WSI.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + DisturbedForest_last20Yrs.yr1 + DecFeb_WSI.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + DisturbedForest_last20Yrs.yr2 + DecFeb_WSI.yr2, data = localN_z),
    data = localN_z
  )
  summary(bottom_up_inter_simple_yr)
  
  bottom_up_inter_simple_yr.r <- psem(
    lm(wolf.yr2 ~ moose.yr1, data = localN_z),
    lm(wolf.yr3 ~ 1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + DisturbedForest_last20Yrs.yr2, data = localN_z),
    moose.yr3 %~~% moose.yr1,
    mountain_lion.yr3 %~~% elk.yr1,
    mountain_lion.yr2 %~~% wolf.yr2,
    data = localN_z
  )
  summary(bottom_up_inter_simple_yr.r)

  
  #'  ----------------------------------------
  #####  Bottom up, exploitative competition  #####
  #'  ----------------------------------------
  bottom_up_exploit_yr <- psem(
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1 + elk.yr1 + whitetailed_deer.yr1 + DisturbedForest_last20Yrs.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + elk.yr2 + whitetailed_deer.yr2 + DisturbedForest_last20Yrs.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + wolf.yr1 + mountain_lion.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + wolf.yr2 + mountain_lion.yr2 + whitetailed_deer.yr2, data = localN_z),
    # lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    # lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    # lm(elk.yr2 ~ elk.yr1 + wolf.yr1, data = localN_z),
    # lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    # lm(moose.yr2 ~ moose.yr1, data = localN_z),
    # lm(moose.yr3 ~ moose.yr2, data = localN_z),
    data = localN_z
  )
  summary(bottom_up_exploit_yr)
  AIC_psem(bottom_up_exploit_yr, AIC.type = "loglik")
  
  bottom_up_exploit_yr.r <- psem(
    lm(wolf.yr2 ~ moose.yr1, data = localN_z),
    lm(wolf.yr3 ~ 1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
    # lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    # lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    # lm(elk.yr2 ~ elk.yr1 + wolf.yr1, data = localN_z),
    # lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    # lm(moose.yr2 ~ moose.yr1, data = localN_z),
    # lm(moose.yr3 ~ moose.yr2, data = localN_z),
    # bear_black.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr1,
    mountain_lion.yr2 %~~% wolf.yr2,
    mountain_lion.yr3 %~~% coyote.yr2,
    coyote.yr3 %~~% bear_black.yr3,
    data = localN_z
  )
  summary(bottom_up_exploit_yr.r)
  AIC_psem(bottom_up_exploit_yr.r, AIC.type = "loglik")
  
  
  piecewiseSEM:::plot.psem(top_down_inter_yr.r, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(top_down_yr.r, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(top_down_exploit_yr.r, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(bottom_up_inter_yr.r, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(bottom_up_exploit_yr.r, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  
  
  
  
  top_down_inter_yr.r_bootEff <- bootEff(top_down_inter_yr.r, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(top_down_inter_yr.r_bootEff, file = paste0("./Outputs/SEM/top_down_inter_yr.r_bootEff_", Sys.Date(), ".RData"))
  (top_down_inter_yr.r_semEff <- semEff(top_down_inter_yr.r_bootEff))
  summary(top_down_inter_yr.r_semEff)
  save(top_down_inter_yr.r_semEff, file = paste0("./Outputs/SEM/top_down_inter_yr.r_semEff_", Sys.Date(), ".RData"))
  
  top_down_yr.r_bootEff <- bootEff(top_down_yr.r, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(top_down_yr.r_bootEff, file = paste0("./Outputs/SEM/top_down_yr.r_bootEff_", Sys.Date(), ".RData"))
  (top_down_yr.r_semEff <- semEff(top_down_yr.r_bootEff))
  summary(top_down_yr.r_semEff)
  save(top_down_yr.r_semEff, file = paste0("./Outputs/SEM/top_down_yr.r_semEff_", Sys.Date(), ".RData"))
  
  top_down_exploit_yr.r_bootEff <- bootEff(top_down_exploit_yr.r, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(top_down_exploit_yr.r_bootEff, file = paste0("./Outputs/SEM/top_down_exploit_yr.r_bootEff_", Sys.Date(), ".RData"))
  (top_down_exploit_yr.r_semEff <- semEff(top_down_exploit_yr.r_bootEff))
  summary(top_down_exploit_yr.r_semEff)
  save(top_down_exploit_yr.r_semEff, file = paste0("./Outputs/SEM/top_down_exploit_yr.r_semEff_", Sys.Date(), ".RData"))
  
  bottom_up_inter_yr.r_bootEff <- bootEff(bottom_up_inter_yr.r, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(bottom_up_inter_yr.r_bootEff, file = paste0("./Outputs/SEM/bottom_up_inter_yr.r_bootEff_", Sys.Date(), ".RData"))
  (bottom_up_inter_yr.r_semEff <- semEff(bottom_up_inter_yr.r_bootEff))
  summary(bottom_up_inter_yr.r_semEff)
  save(bottom_up_inter_yr.r_semEff, file = paste0("./Outputs/SEM/bottom_up_inter_yr.r_semEff_", Sys.Date(), ".RData"))
  
  bottom_up_exploit_yr.r_bootEff <- bootEff(bottom_up_exploit_yr.r, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(bottom_up_exploit_yr.r_bootEff, file = paste0("./Outputs/SEM/bottom_up_exploit_yr.r_bootEff_", Sys.Date(), ".RData"))
  (bottom_up_exploit_yr.r_semEff <- semEff(bottom_up_exploit_yr.r_bootEff))
  summary(bottom_up_exploit_yr.r_semEff)
  save(bottom_up_exploit_yr.r_semEff, file = paste0("./Outputs/SEM/bottom_up_exploit_yr.r_semEff_", Sys.Date(), ".RData"))
  
  
  