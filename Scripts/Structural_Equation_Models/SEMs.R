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
  source("./Scripts/Structural_Equation_Models/Format_data_for_SEMs.R")
  
  #'  Take a quick look
  head(localN_z); head(localN_z_1YrLag)
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  
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
  #'    4. Add newly significant relationships identified by d-sep test that were initially hypothesized, rerun one at a time
  #'    5. Add newly significant relationships identified by d-sep test that were not initially hypothesized, rerun one at a time
  #'    6. Update unspecified correlation for newly significant relationships identified by d-sep test that are not initially biologically meaningful or can't be included elsewhere owing to feedback loops, rerun
  reduced_mod <- psem(
    lm(elk.T ~ elk.Tminus1 + moose.Tminus1 + bear_black.T + bear_black.Tminus1 + wolf.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + elk.T + mountain_lion.T + mountain_lion.Tminus1 + wolf.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + elk.T + elk.Tminus1 + moose.T + moose.Tminus1 + bear_black.Tminus1 + coyote.T + coyote.Tminus1+ wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + bobcat.Tminus1 + mountain_lion.T + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + elk.Tminus1 + lagomorphs.T + lagomorphs.Tminus1 + bobcat.T + bobcat.Tminus1 + PercDisturbedForest.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + wolf.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + mountain_lion.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(reduced_mod)
  #'  Update with unspecified correlations
  reduced_mod <- update(reduced_mod, bobcat.T %~~% lagomorphs.T) #' Unspecified only b/c a feedback loop is created if included in model
  reduced_mod <- update(reduced_mod, lagomorphs.T %~~% wolf.T) #'  Unspecified b/c not a biological hypothesis (from here down)
  reduced_mod <- update(reduced_mod, elk.T %~~% bobcat.Tminus1) 
  reduced_mod <- update(reduced_mod, elk.T %~~% coyote.T)
  reduced_mod <- update(reduced_mod, moose.T %~~% bobcat.Tminus1)
  reduced_mod <- update(reduced_mod, moose.T %~~% coyote.T)
  reduced_mod <- update(reduced_mod, moose.T %~~% coyote.Tminus1)
  reduced_mod <- update(reduced_mod, whitetailed_deer.T %~~% lagomorphs.Tminus1)
  reduced_mod <- update(reduced_mod, lagomorphs.T %~~% bear_black.T)
  reduced_mod <- update(reduced_mod, lagomorphs.T %~~% bear_black.Tminus1)
  reduced_mod <- update(reduced_mod, bear_black.T %~~% lagomorphs.Tminus1)
  reduced_mod <- update(reduced_mod, bear_black.T %~~% coyote.Tminus1)
  reduced_mod <- update(reduced_mod, mountain_lion.T %~~% bobcat.T)
  reduced_mod <- update(reduced_mod, mountain_lion.T %~~% bobcat.Tminus1)
  reduced_mod <- update(reduced_mod, wolf.T %~~% bobcat.T)
  reduced_mod <- update(reduced_mod, wolf.T %~~% coyote.T)
  reduced_mod <- update(reduced_mod, wolf.T %~~% coyote.Tminus1)
  summary(reduced_mod)
  
  #'  Check for multicollinearity
  RVIF(reduced_mod[[6]]) # double check what [[6]] indexes
  
  #'  Visualize SEM
  piecewiseSEM:::plot.psem(reduced_mod, 
                           node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "orange", fontcolor = "black"), 
                           layout = "tree")
  piecewiseSEM:::plot.psem(reduced_mod, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  
  #'  Run some basic model diagnostics on individuals models
  #'  This is handy: http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/
  elk_mod <- lm(elk.T ~ elk.Tminus1 + moose.Tminus1 + bear_black.T + bear_black.Tminus1 + wolf.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_elk.T, data = localN_z_1YrLag)
  plot(elk_mod)
  moose_mod <- lm(moose.T ~ moose.Tminus1 + elk.T + mountain_lion.T + mountain_lion.Tminus1 + wolf.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag)
  plot(moose_mod)
  wtd_mod <- lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + elk.T + elk.Tminus1 + moose.T + moose.Tminus1 + bear_black.Tminus1 + coyote.T + coyote.Tminus1+ wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag)
  plot(wtd_mod)
  lago_mod <- lm(lagomorphs.T ~ lagomorphs.Tminus1 + bobcat.Tminus1 + mountain_lion.T + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag)
  plot(lago_mod)
  bear_mod <- lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag)
  plot(bear_mod)
  bob_mod <- lm(bobcat.T ~ bobcat.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bear_black.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_bobcat.T, data = localN_z_1YrLag)
  plot(bob_mod)
  coy_mod <- lm(coyote.T ~ coyote.Tminus1 + elk.Tminus1 + lagomorphs.T + lagomorphs.Tminus1 + bobcat.T + bobcat.Tminus1 + PercDisturbedForest.T, weights = precision_coyote.T, data = localN_z_1YrLag)
  plot(coy_mod)
  lion_mod <- lm(mountain_lion.T ~ mountain_lion.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + wolf.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag)
  plot(lion_mod)
  wolf_mod <- lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + mountain_lion.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag)
  plot(wolf_mod)
  #'  Well this is a hot mess. FML.
  
  
  #'  Taking reduced model and including all those correlated errors back 
  #'  into the main model even if they aren't a biological hypothesis to get rid
  #'  of these model updates. Seeing if these unhypothesized relationships are
  #'  showing up as indirect effects.
  reduced_modv2 <- psem(
    lm(elk.T ~ elk.Tminus1 + moose.Tminus1 + bear_black.T + bear_black.Tminus1 + bobcat.Tminus1 + wolf.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_elk.T, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + elk.T + elk.Tminus1 + bobcat.Tminus1 + coyote.Tminus1 + mountain_lion.T + DecFeb_WSI.Tminus1, weights = precision_moose.T, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + elk.T + elk.Tminus1 + moose.T + moose.Tminus1 + lagomorphs.Tminus1 + bear_black.Tminus1 + coyote.T + coyote.Tminus1+ wolf.Tminus1 + PercDisturbedForest.T + PercDisturbedForest.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_whitetailed_deer.T, data = localN_z_1YrLag),
    lm(lagomorphs.T ~ lagomorphs.Tminus1 + bear_black.T + bobcat.Tminus1 + mountain_lion.T + PercDisturbedForest.T + DecFeb_WSI.Tminus1, weights = precision_lagomorphs.T, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + lagomorphs.Tminus1 + coyote.Tminus1 + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.T + lagomorphs.Tminus1 + bear_black.T + wolf.T + wolf.Tminus1 + PercDisturbedForest.T, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + elk.T + moose.T + moose.Tminus1 + lagomorphs.T + lagomorphs.Tminus1 + bobcat.T + bobcat.Tminus1 + wolf.T, weights = precision_coyote.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + lagomorphs.Tminus1 + bobcat.Tminus1 + wolf.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + whitetailed_deer.Tminus1 + bobcat.Tminus1 + coyote.Tminus1 + mountain_lion.Tminus1 + DecFeb_WSI.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(reduced_modv2)
  reduced_mod <- reduced_modv2
  #'  Yeah that didn't work. Some of these relationships are popping as direct effects
  #'  even though not biologically meaningful (e.g., bobcat --> elk)
  
  
  #'  Ultra simple model with only predators affecting predators
  reduced_modv3 <- psem(
    lm(bear_black.T ~ bear_black.Tminus1 + coyote.Tminus1 + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag),
    lm(bobcat.T ~ bobcat.Tminus1 + bear_black.T + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bobcat.T, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + bobcat.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_coyote.T, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + bear_black.Tminus1 + wolf.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1 + bear_black.Tminus1 + mountain_lion.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(reduced_modv3)
  reduced_modv3 <- update(reduced_modv3, wolf.T %~~% coyote.Tminus1)
  reduced_modv3 <- update(reduced_modv3, mountain_lion.T %~~% bobcat.Tminus1)
  reduced_modv3 <- update(reduced_modv3, wolf.T %~~% bobcat.Tminus1)
  reduced_modv3 <- update(reduced_modv3, mountain_lion.T %~~% bear_black.T)
  reduced_modv3 <- update(reduced_modv3, wolf.T %~~% bear_black.T)
  reduced_modv3 <- update(reduced_modv3, coyote.T %~~% bobcat.T)
  reduced_modv3 <- update(reduced_modv3, mountain_lion.T %~~% bobcat.T)
  reduced_modv3 <- update(reduced_modv3, wolf.T %~~% bobcat.T)
  reduced_modv3 <- update(reduced_modv3, wolf.T %~~% coyote.T)
  reduced_modv3 <- update(reduced_modv3, wolf.T %~~% mountain_lion.T)
  summary(reduced_modv3)
  
  piecewiseSEM:::plot.psem(reduced_modv3, 
                           node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "orange", fontcolor = "black"), 
                           layout = "tree")
  
  bear_mod <- lm(bear_black.T ~ bear_black.Tminus1 + coyote.Tminus1 + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bear_black.T, data = localN_z_1YrLag)
  summary(bear_mod); plot(bear_mod)
  bob_mod <- lm(bobcat.T ~ bobcat.Tminus1 + bear_black.T + wolf.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_bobcat.T, data = localN_z_1YrLag)
  summary(bob_mod); plot(bob_mod)
  coy_mod <- lm(coyote.T ~ coyote.Tminus1 + bobcat.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_coyote.T, data = localN_z_1YrLag)
  summary(coy_mod); plot(coy_mod)
  lion_mod <- lm(mountain_lion.T ~ mountain_lion.Tminus1 + bear_black.Tminus1 + wolf.Tminus1, weights = precision_mountain_lion.T, data = localN_z_1YrLag)
  summary(lion_mod); plot(lion_mod)
  wolf_mod <- lm(wolf.T ~ wolf.Tminus1 + bear_black.Tminus1 + mountain_lion.Tminus1 + PercDisturbedForest.Tminus1, weights = precision_wolf.T, data = localN_z_1YrLag)
  summary(wolf_mod); plot(wolf_mod)
  
  
  
  #'  Calculate direct, indirect, total, and mediator effects (SE & 95% CI) for 
  #'  all endogenous (response) variables using semEff package
  #'  https://murphymv.github.io/semEff/articles/semEff.html
  #'  First bootstrap standardized model coefficients (necessary for calculating SEs)
  #'  THIS TAKES AWHILE! 
  reduced_mod_bootEff <- bootEff(reduced_mod, R = 5000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(reduced_mod_bootEff, file = paste0("./Outputs/SEM/reduced_mod_bootEff_", Sys.Date(), ".RData"))
  
  #'  Second calculate standardized effects for all casual pathways
  #'  Note, these are standardized unique effects (i.e., adjusted for multicollinearity;
  #'  i.e., semipartial correlations), allowing us to fully partition effects in the system
  #'  These tend to be smaller than the unadjusted standardized coefficients
  #'  If there are large differences btwn the two then consideration should be given
  #'  to the impact and relevance of multicollinearity in the system (check with RVIF())
  (reduced_mod_semEff <- semEff(reduced_mod_bootEff))
  summary(reduced_mod_semEff)
  save(reduced_mod_semEff, file = paste0("./Outputs/SEM/reduced_mod_semEff_", Sys.Date(), ".RData"))
  
  #'  Pull out individual direct, indirect, total, and mediator effects
  directEff <- getDirEff(reduced_mod_semEff)
  indirectEff <- getIndEff(reduced_mod_semEff)
  totalEff <- getTotEff(reduced_mod_semEff)
  mediatorEff <- getMedEff(reduced_mod_semEff)
  
  load("./Outputs/SEM/reduced_mod_semEff_2024-06-04.RData") #'  Original reduced model
  load("./Outputs/SEM/reduced_mod_semEff_2024-06-05.RData") #'  reduced model v2
  
  plot(reduced_mod, return = TRUE)
  plot(reduced_mod, node_attrs = data.frame(shape = "rectangle", fixedsize = FALSE, color = "black", fillcolor = "lightblue", fontcolor = "black", fontsize = 14),
       edge_attrs = data.frame(style = "solid", color = "gray15", arrowsize = 0.75, fontsize = 14),
       ns_dashed = TRUE, 
       alpha = 0.05,
       show = "std",
       digits = 2, 
       add_edge_label_spaces = TRUE)
  
  
  
  
  #' #'  ---------------------------------------------------
  #' ####  SEM with 1-year time lag, with annual variation  ####
  #' #'  ---------------------------------------------------
  #' #'  Starting with saturated model (all possible relationships, except where 
  #' #'  feedback loops arise) and reducing model to only significant relationships. 
  #' #'  This version allows annual local abundance estimates to be stand alone 
  #' #'  variables and does not assume relationships between species are the same
  #' #'  from one year to the next (allows more dynamic relationships to occur).
  #' 
  #' #####  Saturated model  #####
  #' saturated_mod_annual <- psem(
  #'   lm(elk.yr2 ~ elk.yr1 + moose.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + bear_black.yr1 + bobcat.yr1 + coyote.yr1 + mountain_lion.yr1 + wolf.yr1, weights = precision_elk.yr2, data = localN_z),
  #'   lm(elk.yr3 ~ elk.yr2 + moose.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + bear_black.yr2 + bobcat.yr2 + coyote.yr2 + mountain_lion.yr2 + wolf.yr2, weights = precision_elk.yr3, data = localN_z),
  #'   lm(moose.yr2 ~ moose.yr1 + lagomorphs.yr1 + bobcat.yr1 + coyote.yr1 + wolf.yr1, weights = precision_moose.yr2, data = localN_z),
  #'   lm(moose.yr3 ~ moose.yr2 + lagomorphs.yr2 + bobcat.yr2 + coyote.yr2 + mountain_lion.yr2 + wolf.yr2, weights = precision_moose.yr3, data = localN_z),
  #'   lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + elk.yr1 + lagomorphs.yr1 + bear_black.yr1 + bobcat.yr1 + coyote.yr1 + mountain_lion.yr1 + wolf.yr1, weights = precision_whitetailed_deer.yr2, data = localN_z),
  #'   lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + moose.yr2 + lagomorphs.yr2 + bear_black.yr2 + bobcat.yr2 + coyote.yr2 + mountain_lion.yr2 + wolf.yr2, weights = precision_whitetailed_deer.yr3, data = localN_z),
  #'   lm(lagomorphs.yr2 ~ lagomorphs.yr1 + elk.yr1 + moose.yr1 + whitetailed_deer.yr1 + bear_black.yr1 + bobcat.yr1 + coyote.yr1 + mountain_lion.yr1, weights = precision_lagomorphs.yr2, data = localN_z),
  #'   lm(lagomorphs.yr3 ~ lagomorphs.yr2 + elk.yr2 + moose.yr2 + whitetailed_deer.yr2 + bear_black.yr2 + bobcat.yr2 + coyote.yr2, weights = precision_lagomorphs.yr3, data = localN_z),
  #'   lm(bear_black.yr2 ~ bear_black.yr1 + elk.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + bobcat.yr1 + mountain_lion.yr1 + wolf.yr1, weights = precision_bear_black.yr2, data = localN_z),
  #'   lm(bear_black.yr3 ~ bear_black.yr2 + elk.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + bobcat.yr2 + wolf.yr2, weights = precision_bear_black.yr3, data = localN_z),
  #'   lm(bobcat.yr2 ~ bobcat.yr1 + elk.yr1 + moose.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr1 + mountain_lion.yr1 + wolf.yr1, weights = precision_bobcat.yr2, data = localN_z),
  #'   lm(bobcat.yr3 ~ bobcat.yr2 + elk.yr2 + moose.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + bear_black.yr2 + coyote.yr2 + mountain_lion.yr2 + wolf.yr2, weights = precision_bobcat.yr3, data = localN_z),
  #'   lm(coyote.yr2 ~ coyote.yr1 + moose.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + bobcat.yr1 + mountain_lion.yr1 + wolf.yr1, weights = precision_coyote.yr2, data = localN_z),
  #'   lm(coyote.yr3 ~ coyote.yr2 + moose.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + bobcat.yr2 + mountain_lion.yr2 + wolf.yr2, weights = precision_coyote.yr3, data = localN_z),
  #'   lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1 + moose.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + bear_black.yr1 + bobcat.yr1 + wolf.yr1, weights = precision_mountain_lion.yr2, data = localN_z),
  #'   lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + bear_black.yr2 + bobcat.yr2 + wolf.yr2, weights = precision_mountain_lion.yr3, data = localN_z),
  #'   lm(wolf.yr2 ~ wolf.yr1 + elk.yr1 + moose.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + bobcat.yr1 + mountain_lion.yr1, weights = precision_wolf.yr2, data = localN_z),
  #'   lm(wolf.yr3 ~ wolf.yr2 + elk.yr2 + moose.yr2 + whitetailed_deer.yr2 + lagomorphs.yr2 + bear_black.yr2 + bobcat.yr2 + coyote.yr2 + mountain_lion.yr2, weights = precision_wolf.yr3, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(saturated_mod_annual)
  #' 
  #' #####  Reduced model  #####
  #' #'  Method to this madness:
  #' #'  For all linear models at once-
  #' #'    1. Remove non-significant relationships that are biologically not meaningful anyway, rerun
  #' #'    2. Remove non-significant relationships that were not initially hypothesized, rerun
  #' #'    3. Remove non-significant relationships that were initially hypothesized, rerun
  #' #'    4. Add 2-yr lag autocorrelation for individual species (i.e. sppA yr1 --> sppB yr3), rerun & remove newly nonsig relationships
  #' #'    5. Add significant relationships identified by d-sep test that were initially hypothesized (including 0-yr/2-yr lag relationships), rerun one at a time
  #' #'    6. Add significant relationships identified by d-sep test that were not initially hypothesized, rerun one at a time
  #' #'    7. Remove significant biologically implausible relationships and update remaining significant relationships as unspecified correlation structure
  #' reduced_mod_annual <- psem(
  #'   lm(elk.yr2 ~ elk.yr1 + moose.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + bobcat.yr1 + mountain_lion.yr1, weights = precision_elk.yr2, data = localN_z),
  #'   lm(elk.yr3 ~ elk.yr1 + elk.yr2 + moose.yr2 + whitetailed_deer.yr3 + bear_black.yr1 + bear_black.yr3 + coyote.yr2 + wolf.yr2, weights = precision_elk.yr3, data = localN_z),
  #'   lm(moose.yr2 ~ moose.yr1 + lagomorphs.yr1 + bobcat.yr1 + coyote.yr1, weights = precision_moose.yr2, data = localN_z),
  #'   lm(moose.yr3 ~ moose.yr1 + moose.yr2 + bobcat.yr2 + coyote.yr2 + mountain_lion.yr2, weights = precision_moose.yr3, data = localN_z),
  #'   lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + elk.yr1 + elk.yr2 + moose.yr1 + moose.yr2 + lagomorphs.yr1 + bobcat.yr1 + mountain_lion.yr1 + wolf.yr1, weights = precision_whitetailed_deer.yr2, data = localN_z),
  #'   lm(whitetailed_deer.yr3 ~ moose.yr2 + moose.yr3 + lagomorphs.yr2 + bobcat.yr2 + coyote.yr1 + mountain_lion.yr2 + wolf.yr2, weights = precision_whitetailed_deer.yr3, data = localN_z),
  #'   lm(lagomorphs.yr2 ~ lagomorphs.yr1 + bear_black.yr1 + bobcat.yr1 + mountain_lion.yr1, weights = precision_lagomorphs.yr2, data = localN_z),
  #'   lm(lagomorphs.yr3 ~ lagomorphs.yr2 + bear_black.yr2 + bobcat.yr2 + mountain_lion.yr1 + mountain_lion.yr2 + mountain_lion.yr3, weights = precision_lagomorphs.yr3, data = localN_z),
  #'   lm(bear_black.yr2 ~ bear_black.yr1 + whitetailed_deer.yr1 + mountain_lion.yr1 + wolf.yr1, weights = precision_bear_black.yr2, data = localN_z),
  #'   lm(bear_black.yr3 ~ bear_black.yr1 + bear_black.yr2 + whitetailed_deer.yr3 + lagomorphs.yr2 + wolf.yr2, weights = precision_bear_black.yr3, data = localN_z),
  #'   lm(bobcat.yr2 ~ bobcat.yr1 + whitetailed_deer.yr1 + lagomorphs.yr1 + coyote.yr2 + mountain_lion.yr1 + wolf.yr1 + wolf.yr2, weights = precision_bobcat.yr2, data = localN_z),
  #'   lm(bobcat.yr3 ~ bobcat.yr1 + lagomorphs.yr1 + lagomorphs.yr3 + bear_black.yr2 + coyote.yr2 + mountain_lion.yr2 + wolf.yr2, weights = precision_bobcat.yr3, data = localN_z),
  #'   lm(coyote.yr2 ~ coyote.yr1 + moose.yr1 + whitetailed_deer.yr1 + whitetailed_deer.yr2 + lagomorphs.yr1 + bobcat.yr1 + mountain_lion.yr1, weights = precision_coyote.yr2, data = localN_z),
  #'   lm(coyote.yr3 ~ whitetailed_deer.yr1 + whitetailed_deer.yr3 + lagomorphs.yr1 + lagomorphs.yr2 + lagomorphs.yr3 + bobcat.yr1 + bobcat.yr3 + mountain_lion.yr3 + wolf.yr2, weights = precision_coyote.yr3, data = localN_z),
  #'   lm(mountain_lion.yr2 ~ mountain_lion.yr1 + moose.yr1 + whitetailed_deer.yr1 + bobcat.yr1 + wolf.yr1, weights = precision_mountain_lion.yr2, data = localN_z),
  #'   lm(mountain_lion.yr3 ~ mountain_lion.yr1 + elk.yr2 + bobcat.yr2 + wolf.yr1, weights = precision_mountain_lion.yr3, data = localN_z),
  #'   lm(wolf.yr2 ~ wolf.yr1 + moose.yr2 + whitetailed_deer.yr1 + lagomorphs.yr1, weights = precision_wolf.yr2, data = localN_z),
  #'   lm(wolf.yr3 ~ wolf.yr1 + wolf.yr2 + elk.yr1 + moose.yr1 + whitetailed_deer.yr1 + bear_black.yr1 + mountain_lion.yr3, weights = precision_wolf.yr3, data = localN_z),
  #'   data = localN_z
  #' )
  #' summary(reduced_mod_annual)
  #' reduced_mod_annual <- update(reduced_mod_annual, elk.yr3 %~~% moose.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, elk.yr3 %~~% whitetailed_deer.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, elk.yr3 %~~% lagomorphs.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, moose.yr3 %~~% elk.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, moose.yr3 %~~% elk.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual, moose.yr3 %~~% elk.yr3)
  #' reduced_mod_annual <- update(reduced_mod_annual, moose.yr3 %~~% bear_black.yr3)
  #' reduced_mod_annual <- update(reduced_mod_annual, whitetailed_deer.yr3 %~~% elk.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, lagomorphs.yr3 %~~% lagomorphs.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, lagomorphs.yr3 %~~% bear_black.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, lagomorphs.yr3 %~~% bear_black.yr3)
  #' reduced_mod_annual <- update(reduced_mod_annual, bear_black.yr2 %~~% elk.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual, bear_black.yr2 %~~% lagomorphs.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual, bear_black.yr2 %~~% coyote.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, bear_black.yr3 %~~% bobcat.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, bobcat.yr3 %~~% elk.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, bobcat.yr3 %~~% elk.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual, coyote.yr3 %~~% moose.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, coyote.yr3 %~~% bear_black.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, coyote.yr3 %~~% bear_black.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual, coyote.yr3 %~~% mountain_lion.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual, mountain_lion.yr2 %~~% bobcat.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual, mountain_lion.yr3 %~~% coyote.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual, mountain_lion.yr3 %~~% moose.yr3)
  #' reduced_mod_annual <- update(reduced_mod_annual,  wolf.yr2 %~~% coyote.yr1)
  #' reduced_mod_annual <- update(reduced_mod_annual,  wolf.yr2 %~~% coyote.yr2)
  #' reduced_mod_annual <- update(reduced_mod_annual,  wolf.yr3 %~~% mountain_lion.yr2)
  #' #'  Weird correlation structures that are biologically impossible but need to be accounted for (yr3 affects yr2)
  #' # reduced_mod_annual <- update(reduced_mod_annual, whitetailed_deer.yr2 %~~% elk.yr3)
  #' # reduced_mod_annual <- update(reduced_mod_annual, whitetailed_deer.yr2 %~~% whitetailed_deer.yr3)
  #' # reduced_mod_annual <- update(reduced_mod_annual, whitetailed_deer.yr2 %~~% bear_black.yr3)
  #' # reduced_mod_annual <- update(reduced_mod_annual, lagomorphs.yr2 %~~% elk.yr3)
  #' summary(reduced_mod_annual)
  #' 
  #' #'  I don't trust this version. Some weird correlations (yr3 pop affects yr2 pop) 
  #' #'  and SO MUCH unspecified correlation structure. Plus d-sep test indicates
  #' #'  many of these correlations exist but they are not significant when estimated 
  #' #'  in the model.
  
  
  
  
  
  
  