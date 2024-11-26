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
  head(density_wide_1YrLag_20s_22s)
  # head(localN_z); head(localN_z_1YrLag); head(localN_z_all)
  # localN_z_all <- localN_z_all %>%
  #   mutate(year = as.factor(year),
  #          year = as.numeric(year)) 
  
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
  
  top_down_inter.b <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.T, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + moose.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_inter.b)
  
  
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
  
  
  top_down_exploit.b <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + mountain_lion.Tminus1 + elk.T, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + moose.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  ) 
  summary(top_down_exploit.b)
  
  
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
  
  
  bottom_up_inter.b <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s), 
    coyote.T %~~% mountain_lion.T,
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_inter.b)
  
  
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
  
  
  bottom_up.b <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up.b)
  
  
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
  
  top_down_bottom_up.b <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.T, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + moose.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_bottom_up.b)
  
  
  #'  -----------------------------
  ####  Model selection using AIC  ####
  #'  -----------------------------
  #'  Include all permutations of these models (original & reduced via dSep)
  modSelect <- AIC(top_down_inter, top_down_inter.a, #top_down_inter_simple, top_down_inter_simple.a, 
                   top_down_exploit, top_down_exploit.a, #top_down_exploit_simple, top_down_exploit_simple.a, 
                   bottom_up_inter, bottom_up_inter.a, #bottom_up_inter_simple, bottom_up_inter_simple.a, 
                   bottom_up, bottom_up.a, #bottom_up_simple, bottom_up_simple.a, 
                   top_down_bottom_up, top_down_bottom_up.a, 
                   top_down_inter.b, top_down_exploit.b, bottom_up_inter.b, bottom_up.b, top_down_bottom_up.b) 
  mod_names <- c("top_down_inter", "top_down_inter_reduced", #"top_down_inter_simple", "top_down_inter_simple_reduced", 
                 "top_down_exploit", "top_down_exploit_reduced", #"top_down_exploit_simple", "top_down_exploit_simple_reduced", 
                 "bottom_up_inter", "bottom_up_inter_reduced", #"bottom_up_inter_simple", "bottom_up_inter_simple_reduced", 
                 "bottom_up", "bottom_up_reduced", #"bottom_up_simple", "bottom_up_simple_reduced", 
                 "top_down_bottom_up", "top_down_bottom_up_reduced",
                 "top_down_inter.b", "top_down_exploit.b", "bottom_up_inter.b", "bottom_up.b", "top_down_bottom_up.b") 
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("model", "AIC", "K", "n")
  arrange(modSelect, AIC, decreasing = TRUE)
  
  #'  Just the reduced models 
  modSelect <- AIC(#top_down_inter.a, top_down_exploit.a, #top_down_inter_simple.a, top_down_exploit_simple.a, 
                   #bottom_up_inter.a, bottom_up.a, #bottom_up_inter_simple.a, bottom_up_simple.a, 
                   #top_down_bottom_up.a,
                   top_down_inter.b, top_down_exploit.b, bottom_up_inter.b, bottom_up.b, top_down_bottom_up.b)  
  mod_names <- c(#"top_down_inter_reduced", "top_down_exploit_reduced", #"top_down_inter_simple_reduced", "top_down_exploit_simple_reduced", 
                 #"bottom_up_inter_reduced", "bottom_up_reduced", #"bottom_up_inter_simple_reduced", "bottom_up_simple_reduced", 
                 #"top_down_bottom_up_reduced",
                 "top_down_inter.b", "top_down_exploit.b", "bottom_up_inter.b", "bottom_up.b", "top_down_bottom_up.b") 
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("Model", "AIC", "K", "n")
  (sem_aic <- arrange(modSelect, AIC, decreasing = TRUE))
  
  #'  Review Fisher's C and Chi-square test statistics for each model
  sem_gof <- function(mod, gof, modname) {
    model <- modname
    chisq <- LLchisq(mod)
    chisq_stat <- round(chisq[[1]], 2)
    chisq_df <- chisq[[2]]
    chisq_pval <- round(chisq[[3]], 2)
    fishers <- fisherC(mod)
    fishers_stat <- round(fishers[[1]], 2)
    fishers_df <- fishers[[2]]
    fishers_pval <- round(fishers[[3]], 2)
    GoF <- data.frame(model, chisq_stat, chisq_df, chisq_pval, fishers_stat, fishers_df, fishers_pval)
    names(GoF) <- c("Model", "Chi2", "Chi2 df", "Chi2 p-value", "Fisher's C", "Fisher's C df", "Fisher's C p-value")
    return(GoF)
  }
  mod_list <- list(top_down_inter.b, top_down_exploit.b, bottom_up_inter.b, bottom_up.b, top_down_bottom_up.b) 
  mod_names <- c("top_down_inter.b", "top_down_exploit.b", "bottom_up_inter.b", "bottom_up.b", "top_down_bottom_up.b")
  sem_GoF <- mapply(sem_gof, mod_list, modname = mod_names, SIMPLIFY = FALSE) %>% bind_rows()
  
  modSel_table <- full_join(sem_aic, sem_GoF, by = "Model") %>%
    mutate(AIC = round(AIC, 2),
           deltaAIC = AIC - first(x = AIC),
           Model = ifelse(Model == "top_down_inter.b", "Top-down, interference competition", Model),
           Model = ifelse(Model == "top_down_exploit.b", "Top-down, exploitation competition", Model),
           Model = ifelse(Model == "bottom_up_inter.b", "Bottom-up, interference competition", Model),
           Model = ifelse(Model == "bottom_up.b", "Bottom-up, exploitation competition", Model),
           Model = ifelse(Model == "top_down_bottom_up.b", "Top-down, bottom-up", Model)) %>%
    relocate(deltaAIC, .after = AIC) %>%
    dplyr::select(-n)
  write_csv(modSel_table, "./Outputs/SEM/ModSel_GoF_table.csv")
  
  #'  ----------------------
  ####  Dig into top model  ####
  #'  ----------------------
  #'  Dig into top model based on AIC, Fisher's C
  #'  Check for multicollinearity
  RVIF(bottom_up.b[[1]]) 
  RVIF(bottom_up.b[[2]]) 
  RVIF(bottom_up.b[[3]]) 
  RVIF(bottom_up.b[[4]]) 
  RVIF(bottom_up.b[[5]]) 
  RVIF(bottom_up.b[[6]]) 
  RVIF(bottom_up.b[[7]]) 
  
  RVIF(bottom_up_inter.b[[1]]) 
  RVIF(bottom_up_inter.b[[2]]) 
  RVIF(bottom_up_inter.b[[3]]) 
  RVIF(bottom_up_inter.b[[4]]) 
  RVIF(bottom_up_inter.b[[5]]) 
  RVIF(bottom_up_inter.b[[6]]) 
  RVIF(bottom_up_inter.b[[7]]) 
  
  #'  Run some basic model diagnostics on individuals models
  #'  This is handy: http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/
  wtd_mod <- lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s)
  plot(wtd_mod)
  elk_mod <- lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s)
  plot(elk_mod)
  moose_mod <- lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(moose_mod)
  lion_mod <- lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(lion_mod)
  bear_mod <- lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(bear_mod)
  wolf_mod <- lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(wolf_mod)
  coy_mod <- lm(coyote.T ~ coyote.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s)
  plot(coy_mod)
  
  #'  ------------------
  #####  Visualize SEM  #####
  #'  ------------------
  #'  (top but also next top models that are all within 10 deltaAIC)
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
  
  
  
  piecewiseSEM:::plot.psem(bottom_up.b, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  piecewiseSEM:::plot.psem(bottom_up_inter.b, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  piecewiseSEM:::plot.psem(top_down_inter.b, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  piecewiseSEM:::plot.psem(top_down_bottom_up.b, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  piecewiseSEM:::plot.psem(top_down_exploit.b, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  
  #'  -------------------------------
  #####  DIRECT vs INDIRECT EFFECTS  #####
  #'  -------------------------------
  #'  Calculate direct, indirect, total, and mediator effects (SE & 95% CI) for 
  #'  all endogenous (response) variables using semEff package
  #'  https://murphymv.github.io/semEff/articles/semEff.html
  #'  
  #'  First: bootstrap standardized model coefficients (necessary for calculating SEs)
  #'  THIS TAKES AWHILE! 
  #'  
  #'  Second: calculate standardized effects for all casual pathways
  #'  Note, these are standardized unique effects (i.e., adjusted for multicollinearity;
  #'  i.e., semipartial correlations), allowing us to fully partition effects in the system
  #'  These tend to be smaller than the unadjusted standardized coefficients
  #'  If there are large differences btwn the two then consideration should be given
  #'  to the impact and relevance of multicollinearity in the system (check with RVIF())
  bottom_up.b_bootEff <- bootEff(bottom_up.b, R = 10000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  (bottom_up.b_semEff <- semEff(bottom_up.b_bootEff))
  summary(bottom_up.b_semEff)
  save(bottom_up.b_bootEff, file = paste0("./Outputs/SEM/bottom_up.b_bootEff_", Sys.Date(), ".RData"))
  save(bottom_up.b_semEff, file = paste0("./Outputs/SEM/bottom_up.b_semEff_", Sys.Date(), ".RData"))
  
  bottom_up_inter.b_bootEff <- bootEff(bottom_up_inter.b, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  (bottom_up_inter.b_semEff <- semEff(bottom_up_inter.b_bootEff))
  summary(bottom_up_inter.b_semEff)
  save(bottom_up_inter.b_bootEff, file = paste0("./Outputs/SEM/bottom_up_inter.b_bootEff_", Sys.Date(), ".RData"))
  save(bottom_up_inter.b_semEff, file = paste0("./Outputs/SEM/bottom_up_inter.b_semEff_", Sys.Date(), ".RData"))
  
  top_down_inter.b_bootEff <- bootEff(top_down_inter.b, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  (top_down_inter.b_semEff <- semEff(top_down_inter.b_bootEff))
  summary(top_down_inter.b_semEff)
  save(top_down_inter.b_bootEff, file = paste0("./Outputs/SEM/top_down_inter.b_bootEff_", Sys.Date(), ".RData"))
  save(top_down_inter.b_semEff, file = paste0("./Outputs/SEM/top_down_inter.b_semEff_", Sys.Date(), ".RData"))
  
  top_down_bottom_up.b_bootEff <- bootEff(top_down_bottom_up.b, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  (top_down_bottom_up.b_semEff <- semEff(top_down_bottom_up.b_bootEff))
  summary(top_down_bottom_up.b_semEff)
  save(top_down_bottom_up.b_bootEff, file = paste0("./Outputs/SEM/top_down_bottom_up.b_bootEff_", Sys.Date(), ".RData"))
  save(top_down_bottom_up.b_semEff, file = paste0("./Outputs/SEM/top_down_bottom_up.b_semEff_", Sys.Date(), ".RData"))
  
  top_down_exploit.b_bootEff <- bootEff(top_down_exploit.b, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  (top_down_exploit.b_semEff <- semEff(top_down_exploit.b_bootEff))
  summary(top_down_exploit.b_semEff)
  save(top_down_exploit.b_bootEff, file = paste0("./Outputs/SEM/top_down_exploit.b_bootEff_", Sys.Date(), ".RData"))
  save(top_down_exploit.b_semEff, file = paste0("./Outputs/SEM/top_down_exploit.b_semEff_", Sys.Date(), ".RData"))
  
  top_down_inter_simple.a_bootEff <- bootEff(top_down_inter_simple.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(top_down_inter_simple.a_bootEff, file = paste0("./Outputs/SEM/top_down_inter_simple.a_bootEff_", Sys.Date(), ".RData"))
  
  bottom_up_inter_simple.a_bootEff <- bootEff(bottom_up_inter_simple.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(bottom_up_inter_simple.a_bootEff, file = paste0("./Outputs/SEM/bottom_up_inter_simple.a_bootEff_", Sys.Date(), ".RData"))
  
  #'  ------------------
  #####  Result tables  #####
  #'  ------------------
  ######  Direct & indirect effects table  ######
  #'  -------------------------------------
  #'  Extract standardized direct, indirect, etc. effects from SEMs after bootstrapping
  library(stringi)
  #'  Create results tables
  tbl_std_est <- function(mod, mod_name) {
    #'  Extract standardized results
    predictor_wtd <- mod$Summary$whitetailed.deer.T
    predictor_elk <- mod$Summary$elk.T
    predictor_moose <- mod$Summary$moose.T
    predictor_wolf <- mod$Summary$wolf.T
    predictor_lion <- mod$Summary$mountain.lion.T
    predictor_bear <- mod$Summary$bear.black.T
    predictor_coy <- mod$Summary$coyote.T
    
    #'  Bind into a single data frame and add column for endogenous variable
    tbl_wtd <- as.data.frame(predictor_wtd) %>% bind_cols("White-tailed deer t")
    tbl_elk <- as.data.frame(predictor_elk) %>% bind_cols("Elk t")
    tbl_moose <- as.data.frame(predictor_moose) %>% bind_cols("Moose t")
    tbl_wolf <- as.data.frame(predictor_wolf) %>% bind_cols("Wolf t")
    tbl_lion <- as.data.frame(predictor_lion) %>% bind_cols("Mountain lion t")
    tbl_bear <- as.data.frame(predictor_bear) %>% bind_cols("Black bear t")
    tbl_coy <- as.data.frame(predictor_coy) %>% bind_cols("Coyote t")
    
    #'  Rename columns
    col_names <- c("Effect type", "Exogenous_variable", "space", "Standardized_effect", "space", "Bias", "space", "Std.Error", "space", "lower_CI", "upper_CI", "space", "Signif", "Endogenous_variable")
    names(tbl_wtd) <- col_names
    names(tbl_elk) <- col_names
    names(tbl_moose) <- col_names
    names(tbl_wolf) <- col_names
    names(tbl_lion) <- col_names
    names(tbl_bear) <- col_names
    names(tbl_coy) <- col_names
    
    #'  Bind results together and clean up final table
    full_tbl <- bind_rows(tbl_wtd, tbl_elk, tbl_moose, tbl_wolf, tbl_lion, tbl_bear, tbl_coy) %>%
      relocate(Endogenous_variable, .before = Exogenous_variable) %>%
      dplyr::select(c("Endogenous_variable", "Exogenous_variable", "Standardized_effect", 
                      "Std.Error", "lower_CI", "upper_CI", "Signif")) %>%
      #'  Remove extra spaces before or end of words
      mutate(across(everything(), ~stri_trim(.))) %>%
      #'  Replace empty cells with NA
      mutate(across(everything(), ~na_if(., ""))) %>%
      #'  Remove any rows with `NA`
      # filter(!if_any(everything(), is.na)) %>%
      filter(!is.na(Exogenous_variable)) %>%
      filter(Exogenous_variable != "n/a") %>%
      #'  Remove duplicate information
      distinct(.) %>%
      #'  Create single column for confidence intervals
      mutate("95% CI" = paste(" ", lower_CI, "-", upper_CI)) %>%
      dplyr::select(-c(lower_CI, upper_CI)) %>%
      arrange(Endogenous_variable, Exogenous_variable) %>%
      #'  Add column reporting model name
      mutate(Model = mod_name) %>%
      relocate(Model, .before = Endogenous_variable) %>%
      relocate(Signif, .after = "95% CI") %>%
      #'  Remove "." and adjust time period indicator
      mutate(Exogenous_variable = ifelse(Exogenous_variable == "bear.black.Tminus1", "bear black.Tminus1", Exogenous_variable),
             Exogenous_variable = ifelse(Exogenous_variable == "mountain.lion.Tminus1", "mountain lion.Tminus1", Exogenous_variable),
             Exogenous_variable = ifelse(Exogenous_variable == "whitetailed.deer.Tminus1", "white-tailed deer.Tminus1", Exogenous_variable),
             Exogenous_variable = gsub(".Tminus1", " t-1", Exogenous_variable),
             Exogenous_variable = str_to_sentence(Exogenous_variable))
    
    return(full_tbl)
  }
  result_tbl_bottom_up.b <- tbl_std_est(bottom_up.b_semEff, mod_name = "Bottom-up, Exploit")
  write_csv(result_tbl_bottom_up.b, "./Outputs/SEM/result_table_top_model_std_effects.csv")
  
  #'  Tables of less supported models
  result_tbl_bottom_up_inter.b <- tbl_std_est(bottom_up_inter.b_semEff, mod_name = "Bottom-up, Inter")
  result_tbl_top_down_inter.b <- tbl_std_est(top_down_inter.b_semEff, mod_name = "Top-down, Exploit")
  result_tbl_top_down_exploit.b <- tbl_std_est(top_down_exploit.b_semEff, mod_name = "Top-down, Inter")
  result_tbl_top_down_bottom_up.b <- tbl_std_est(top_downbottom_up.b_semEff, mod_name = "Top-down, bottom-up")
  result_tbl_top_model <- bind_rows(result_tbl_bottom_up.b, result_tbl_bottom_up_inter.b, result_tbl_top_down_inter.b, result_tbl_top_down_bottom_up.b, result_tbl_top_down_exploit.b)
  write_csv(result_tbl_top_model, "./Outputs/SEM/result_table_top_model_std_effects.csv")
  
  #'  --------------------------------------------
  ######  Unstandardized regression coefficients  ######
  #'  --------------------------------------------
  #'  Run linear regression models from best supported SEM
  wtd_mod <- lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s)
  elk_mod <- lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s)
  moose_mod <- lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  wolf_mod <- lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  lion_mod <- lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s)
  bear_mod <- lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s)
  coy_mod <- lm(coyote.T ~ coyote.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s)
  
  #'  List models
  regression_list <- list(wtd_mod, elk_mod, moose_mod, wolf_mod, lion_mod, bear_mod, coy_mod)
  response_var <- list("White-tailed deer t", "Elk t", "Moose t", "Wolf t", "Mountain lion t", "Black bear t", "Coyote t")
  
  #'  Create results table of unstandardized regression coefficients and R2 values
  SEM_unstandardized_tbl <- function(mod, response_var) {
    #'  Grab coefficient estimates, SE, and p-value from linear model
    Params <- rownames(as.data.frame(mod$coefficients))
    mod_summary <- summary(mod)
    coefEst <- round(mod_summary$coefficients[,1], 2)
    stdErr <- round(mod_summary$coefficients[,2], 2)
    pVal <- round(mod_summary$coefficients[,4], 5)
    r2 <- round(mod_summary$r.squared, 2)
    
    #'  Merge into a single data frame
    lm_out <- bind_cols(response_var, Params, coefEst, stdErr, pVal, r2) 
    names(lm_out) <- c("Response", "Parameter", "Estimate", "SE", "p-value", "R2")
    #'  Clean up parameter names - remove "." and adjust time period indicator
    lm_out <- lm_out %>%
      mutate(Parameter = ifelse(Parameter == "bear_black.Tminus1", "bear black.Tminus1", Parameter),
             Parameter = ifelse(Parameter == "bear_black.T", "bear black.T", Parameter),
             Parameter = ifelse(Parameter == "mountain_lion.Tminus1", "mountain lion.Tminus1", Parameter),
             Parameter = ifelse(Parameter == "mountain_lion.T", "mountain lion.T", Parameter),
             Parameter = ifelse(Parameter == "whitetailed_deer.Tminus1", "white-tailed deer.Tminus1", Parameter),
             Parameter = ifelse(Parameter == "whitetailed_deer.T", "white-tailed deer.T", Parameter),
             Parameter = gsub(".Tminus1", " t-1", Parameter),
             Parameter = gsub(".T", " t", Parameter),
             Parameter = str_to_sentence(Parameter))
   
    return(lm_out) 
  }
  sem_regression_tbl <- mapply(SEM_unstandardized_tbl, regression_list, response_var = response_var, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  #'  Save!
  write_csv(sem_regression_tbl, "./Outputs/SEM/SEM_Regression_coefs_tbl.csv")
  
  
  
  
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
    lm(elk.yr2 ~ elk.yr1 + whitetailed_deer.yr1 + whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + mountain_lion.yr1, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr1, data = localN_z),
    lm(wolf.yr2 ~ moose.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr1 + wolf.yr2 + mountain_lion.yr1 + elk.yr1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr2, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr1 + mountain_lion.yr2 + wolf.yr2 + elk.yr2 + whitetailed_deer.yr1, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + mountain_lion.yr3, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + bear_black.yr3 + mountain_lion.yr3 + elk.yr3, data = localN_z),
    data = localN_z)
  summary(top_down_inter_yr.r)
  
  top_down_inter_yr.r_bootEff <- bootEff(top_down_inter_yr.r, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  (top_down_inter_yr.r_semEff <- semEff(top_down_inter_yr.r_bootEff))
  summary(top_down_inter_yr.r_semEff)
  save(top_down_inter_yr.r_bootEff, file = paste0("./Outputs/SEM/top_down_inter_yr.r_bootEff_", Sys.Date(), ".RData"))
  save(top_down_inter_yr.r_semEff, file = paste0("./Outputs/SEM/top_down_inter_yr.r_semEff_", Sys.Date(), ".RData"))
  
  piecewiseSEM:::plot.psem(top_down_inter_yr.r, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(top_down_inter_yr.r, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  
  
  
  # top_down_inter_yr.r <- psem(
  #   lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
  #   lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
  #   lm(elk.yr2 ~ elk.yr1, data = localN_z),
  #   lm(elk.yr3 ~ elk.yr2, data = localN_z),
  #   lm(moose.yr2 ~ moose.yr1, data = localN_z),
  #   lm(moose.yr3 ~ moose.yr2, data = localN_z),
  #   # lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
  #   # lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
  #   lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
  #   lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
  #   lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
  #   lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
  #   lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
  #   lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
  #   moose.yr3 %~~% moose.yr1,
  #   moose.yr3 %~~% elk.yr1,
  #   elk.yr3 %~~% moose.yr1,
  #   elk.yr3 %~~% mountain_lion.yr1,
  #   coyote.yr3 %~~% bear_black.yr1,
  #   coyote.yr3 %~~% bear_black.yr3,
  #   bear_black.yr3 %~~% bear_black.yr1,
  #   elk.yr2 %~~% whitetailed_deer.yr2,
  #   mountain_lion.yr3 %~~% elk.yr1,
  #   data = localN_z
  # )
  # summary(top_down_inter_yr.r)     
  # AIC_psem(top_down_inter_yr.r, AIC.type = "loglik")
  # 
  # top_down_inter_yr <- psem(
  #   lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
  #   lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
  #   lm(elk.yr2 ~ elk.yr1, data = localN_z),
  #   lm(elk.yr3 ~ elk.yr2, data = localN_z),
  #   lm(moose.yr2 ~ moose.yr1, data = localN_z),
  #   lm(moose.yr3 ~ moose.yr1, data = localN_z),
  #   lm(wolf.yr2 ~ moose.yr1, data = localN_z),
  #   lm(wolf.yr3 ~ wolf.yr1, data = localN_z),
  #   lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
  #   lm(mountain_lion.yr3 ~ mountain_lion.yr1 + mountain_lion.yr2 + elk.yr2, data = localN_z),
  #   lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
  #   lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
  #   lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1, data = localN_z),
  #   lm(coyote.yr3 ~ coyote.yr2 + bear_black.yr1, data = localN_z),
  #   data = localN_z)
  # summary(top_down_inter_yr)
  # AIC_psem(top_down_inter_yr, AIC.type = "loglik")
  
  
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
  
  top_down_exploit_yr.r <- psem(
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
  summary(top_down_exploit_yr.r)
  
  
 
  # top_down_exploit_yr <- psem(
  #   lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + coyote.yr1, data = localN_z),
  #   lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + coyote.yr2, data = localN_z),
  #   lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1, data = localN_z),
  #   lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2, data = localN_z),
  #   lm(moose.yr2 ~ moose.yr1 + wolf.yr1, data = localN_z),
  #   lm(moose.yr3 ~ moose.yr2 + wolf.yr2, data = localN_z),
  #   lm(wolf.yr2 ~ wolf.yr1 + harvest_sqKm.yr1 + harvest_sqKm_quad.yr1, data = localN_z),
  #   lm(wolf.yr3 ~ wolf.yr2 + harvest_sqKm.yr2 + harvest_sqKm_quad.yr2, data = localN_z),
  #   lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
  #   lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
  #   lm(bear_black.yr2 ~ bear_black.yr1 + elk.yr1, data = localN_z),
  #   lm(bear_black.yr3 ~ bear_black.yr2 + elk.yr2, data = localN_z),
  #   lm(coyote.yr2 ~ coyote.yr1 + whitetailed_deer.yr1, data = localN_z),
  #   lm(coyote.yr3 ~ coyote.yr2 + whitetailed_deer.yr2, data = localN_z),
  #   data = localN_z
  # )
  # summary(top_down_exploit_yr)
  # AIC_psem(top_down_exploit_yr, AIC.type = "loglik")
  # 
  # top_down_exploit_yr.r <- psem(
  #   lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
  #   lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
  #   lm(elk.yr2 ~ elk.yr1, data = localN_z),
  #   lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
  #   lm(moose.yr2 ~ moose.yr1, data = localN_z),
  #   lm(moose.yr3 ~ moose.yr2, data = localN_z),
  #   lm(wolf.yr2 ~ 1, data = localN_z),
  #   lm(wolf.yr3 ~ 1, data = localN_z),
  #   lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
  #   lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + wolf.yr2, data = localN_z),
  #   lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
  #   lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
  #   lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
  #   lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
  #   elk.yr2 %~~% whitetailed_deer.yr2,
  #   elk.yr3 %~~% moose.yr1,
  #   elk.yr1 %~~% moose.yr3,
  #   moose.yr3 %~~% moose.yr1,
  #   bear_black.yr3 %~~% bear_black.yr1,
  #   coyote.yr3 %~~% bear_black.yr1,
  #   coyote.yr3 %~~% bear_black.yr3,
  #   mountain_lion.yr3 %~~% whitetailed_deer.yr1,
  #   data = localN_z
  # )
  # summary(top_down_exploit_yr.r)
  # AIC_psem(top_down_exploit_yr.r, AIC.type = "loglik")
  
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
  
  
  