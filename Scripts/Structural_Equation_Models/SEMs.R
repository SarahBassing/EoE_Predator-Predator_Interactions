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
  
  library(piecewiseSEM)
  library(labelled)
  library(lme4)
  library(tidyverse)
  
  #'  Run script to prepare data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_data_for_SEMs.R")
  
  #'  Take a quick look
  head(localN_z); head(localN_z_1YrLag)
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 999999)
  
  
  #'  -------------------------------------------------
  ####  SEM with 1-year time lag, no annual variation  ####
  #'  -------------------------------------------------
  #'  Starting with saturated model (all possible relationships, except where 
  #'  feedback loops arise) and reducing model to only significant reltionships. 
  #'  This version assumes relationships between species are the consistent 
  #'  across years (no annual variation in causal effects) but allows for a
  #'  1-year time lag where last year's population affects this year's population
  #'  of various species.
  
  #####  Saturated model  #####
  saturated_mod <- psem(
    lm(elk.T ~ elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.T + bobcat.Tminus1 + coyote.T + coyote.Tminus1 + mountain_lion.T + mountain_lion.Tminus1 + wolf.T + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + elk.T + elk.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.T + bobcat.Tminus1 + coyote.T + coyote.Tminus1 + mountain_lion.T + mountain_lion.Tminus1 + wolf.T + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + elk.T + elk.Tminus1 + moose.T + moose.Tminus1 + lagomorphs.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.T + bobcat.Tminus1 + coyote.T + coyote.Tminus1 + mountain_lion.T + mountain_lion.Tminus1 + wolf.T + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + elk.T + elk.Tminus1 + moose.T + moose.Tminus1 + whitetailed_deer.T + whitetailed_deer.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.T + bobcat.Tminus1 + coyote.T + coyote.Tminus1 + mountain_lion.T + mountain_lion.Tminus1 + wolf.T + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bobcat.Tminus1 + coyote.Tminus1 + mountain_lion.Tminus1 + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + bear_black.Tminus1 + coyote.Tminus1 + mountain_lion.Tminus1 + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.T + bobcat.Tminus1 + mountain_lion.Tminus1 + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_coyote.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.T + bobcat.Tminus1 + coyote.T + coyote.Tminus1 + wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.T + bobcat.Tminus1 + coyote.T + coyote.Tminus1 + mountain_lion.T + mountain_lion.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(saturated_mod)
  dSep(saturated_mod)
  # Global goodness-of-fit:
  # Chi-Squared = 0 with P-value = 1 and on 0 degrees of freedom
  # Fisher's C = NA with P-value = NA and on 0 degrees of freedom
  
  ####  Reduced model  #####
  #'  Method to this madness:
  #'  For all linear models at once-
  #'    1. Remove non-significant relationships that are biologically not meaningful anyway, rerun
  #'    2. Remove non-significant relationships that were not initially hypothesized, rerun
  #'    3. Remove non-significant relationships that were initially hypothesized, rerun
  #'    4. Add newly significant relationships identified by d-sep test that were initially hypothesized, one at a time & rerun
  #'    5. Add newly significant relationships identified by d-sep test that were not initially hypothesized, one at a time & rerun
  #'    6. Update unspecified correlation for newly significant relationships identified by d-sep test that are not initially biologically meaningful or can't be included elsewhere owing to feedback loops, rerun
  reduced_mod <- psem(
    lm(elk.T ~ elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + bear_black.T + bobcat.Tminus1 + coyote.T + wolf.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + elk.T + elk.Tminus1 + bobcat.Tminus1 + coyote.T + mountain_lion.T + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + elk.T + elk.Tminus1 + moose.T + moose.Tminus1 + lagomorphs.Tminus1 + bear_black.Tminus1 + coyote.T + coyote.Tminus1+ wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + bear_black.T + bobcat.Tminus1 + mountain_lion.T + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + lagomorphs.Tminus1 + coyote.Tminus1 + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + elk.Tminus1 +lagomorphs.Tminus1 + bear_black.T + bobcat.T + bobcat.Tminus1 + PercDisturbedForest.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bobcat.T + bobcat.Tminus1 + wolf.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + bobcat.T + coyote.T + coyote.Tminus1 + mountain_lion.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(reduced_mod)
  reduced_mod <- update(reduced_mod, bobcat.T %~~% lagomorphs.T)
  reduced_mod <- update(reduced_mod, bobcat.T %~~% elk.T)
  reduced_mod <- update(reduced_mod, lagomorphs.T %~~% wolf.T)
  summary(reduced_mod)
  
  ####  Extra reduced model  #####
  #'  Method to this madness:
  #'  Using the reduced model above
  #'    1.  
  reduced_mod <- psem(
    lm(elk.T ~ elk.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + bear_black.T + bobcat.Tminus1 + coyote.T + wolf.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + elk.T + elk.Tminus1 + bobcat.Tminus1 + coyote.T + mountain_lion.T + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + elk.T + elk.Tminus1 + moose.T + moose.Tminus1 + lagomorphs.Tminus1 + bear_black.Tminus1 + coyote.T + coyote.Tminus1+ wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + bear_black.T + bobcat.Tminus1 + mountain_lion.T + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + lagomorphs.Tminus1 + coyote.Tminus1 + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + elk.Tminus1 +lagomorphs.Tminus1 + bear_black.T + bobcat.T + bobcat.Tminus1 + PercDisturbedForest.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bobcat.T + bobcat.Tminus1 + wolf.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + bobcat.T + coyote.T + coyote.Tminus1 + mountain_lion.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(reduced_mod)
  reduced_mod <- update(reduced_mod, bobcat.T %~~% lagomorphs.T)
  reduced_mod <- update(reduced_mod, bobcat.T %~~% elk.T)
  reduced_mod <- update(reduced_mod, lagomorphs.T %~~% wolf.T)
  summary(reduced_mod)

  #'  ---------------------------------------------------
  ####  SEM with 1-year time lag, with annual variation  ####
  #'  ---------------------------------------------------
  #'  Starting with saturated model (all possible relationships, except where 
  #'  feedback loops arise) and reducing model to only significant reltionships. 
  #'  This version allows annual local abundance estimates to be stand alone 
  #'  variables and does not assume relationships between species are the same
  #'  from one year to the next (allows more dynamic relationships to occur).
  
  #####  Saturated model  #####
  
  
  ######  H.td1  ######
  #'  Predators affect their primary prey species; predators affect their competitors
  #'  through interference competition
  H.td1_psem <- psem(
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z, weights = precision_wolf.yr2),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z, weights = precision_wolf.yr3),
    lm(mountain_lion.yr2 ~ wolf.yr1 + bear_black.yr1 + mountain_lion.yr1, data = localN_z, weights = precision_mountain_lion.yr2),
    lm(mountain_lion.yr3 ~ wolf.yr2 + bear_black.yr2 + mountain_lion.yr2, data = localN_z, weights = precision_mountain_lion.yr3),
    lm(bear_black.yr2 ~ wolf.yr1 + bear_black.yr1, data = localN_z, weights = precision_bear_black.yr2),
    lm(bear_black.yr3 ~ wolf.yr2 + bear_black.yr2, data = localN_z, weights = precision_bear_black.yr3),
    lm(coyote.yr2 ~ wolf.yr1 + mountain_lion.yr1 + coyote.yr1, data = localN_z, weights = precision_coyote.yr2),
    lm(coyote.yr3 ~ wolf.yr2 + mountain_lion.yr2 + coyote.yr2, data = localN_z, weights = precision_coyote.yr3),
    lm(bobcat.yr2 ~ coyote.yr1 + bobcat.yr1, data = localN_z, weights = precision_bobcat.yr2),
    lm(bobcat.yr3 ~ coyote.yr2 + bobcat.yr2, data = localN_z, weights = precision_bobcat.yr3),
    lm(moose.yr2 ~ wolf.yr1 + moose.yr1, data = localN_z, weights = precision_moose.yr2),
    lm(moose.yr3 ~ wolf.yr2 + moose.yr2, data = localN_z, weights = precision_moose.yr3),
    lm(elk.yr2 ~ wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + elk.yr1, data = localN_z, weights = precision_elk.yr2),
    lm(elk.yr3 ~ wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + elk.yr2, data = localN_z, weights = precision_elk.yr3),
    lm(whitetailed_deer.yr2 ~ mountain_lion.yr1 + bear_black.yr1 + whitetailed_deer.yr1, data = localN_z, weights = precision_whitetailed_deer.yr2),
    lm(whitetailed_deer.yr3 ~ mountain_lion.yr2 + bear_black.yr2 + whitetailed_deer.yr2, data = localN_z, weights = precision_whitetailed_deer.yr3),
    lm(lagomorphs.yr2 ~ coyote.yr1 + bobcat.yr1 + lagomorphs.yr1, data = localN_z, weights = precision_lagomorphs.yr2),
    lm(lagomorphs.yr3 ~ coyote.yr2 + bobcat.yr2 + lagomorphs.yr2, data = localN_z, weights = precision_lagomorphs.yr3),
    data = localN_z)
  # basisSet(H.td1_psem)
  # dSep(H.td1_psem)
  summary(H.td1_psem)
  
  
  
  
  