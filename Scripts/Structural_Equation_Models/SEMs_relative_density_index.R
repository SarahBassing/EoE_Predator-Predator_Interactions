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
  library(semEff)
  library(labelled)
  library(DiagrammeR)
  library(lme4)
  library(tidyverse)
  
  #'  Run script to prepare data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_density_data_for_SEMs.R")
  
  #'  Take a quick look
  head(localN_z); head(localN_z_1YrLag); head(localN_z_all)
  localN_z_all2 <- localN_z_all %>%
    mutate(year = as.factor(year),
           year = as.numeric(year)) 
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  
  #'  -------------------------------------------------
  ####  SEM with 1-year time lag, no annual variation  ####
  #'  -------------------------------------------------
  #'  Top down hypothesis, interference competition btwn predators
  top_down1 <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(top_down1)
  dSep(top_down1)
  
  #'  Reduced model
  #'  For all linear models at once-
  #'    1. Remove non-significant relationships that were initially hypothesized, rerun
  #'    2. Add newly significant relationships identified by d-sep test that were initially hypothesized, rerun one at a time
  #'    3. Add newly significant relationships identified by d-sep test that were not initially hypothesized but biologically plausible, rerun one at a time
  #'    4. Update unspecified correlation for newly significant relationships identified by d-sep test that are not initially biologically meaningful or can't be included elsewhere owing to feedback loops, rerun
  top_down1a <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(top_down1a)
  
  
  
  