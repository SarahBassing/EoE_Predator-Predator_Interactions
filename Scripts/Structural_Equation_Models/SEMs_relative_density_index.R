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
  #'  ---------------------------------------
  #####  Top down, interference competition  #####
  #'  ---------------------------------------
  top_down_inter <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag),
    lm(wolf.T ~ wolf.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = localN_z_1YrLag),
    data = localN_z_1YrLag
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
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1, data = localN_z_1YrLag),
    lm(coyote.T ~ coyote.Tminus1, data = localN_z_1YrLag),
    #'  Update with unspecified correlated errors - focus on exogenous variables 
    #'  (year T variables) that are not presumed to have a causal & unidirectional
    #'  relationship but are related by an underlying factor making them appear 
    #'  correlated (e.g., similar habitat use)
    elk.T %~~% whitetailed_deer.T,
    coyote.T %~~% mountain_lion.T,
    bear_black.T %~~% whitetailed_deer.T,
    bear_black.T %~~% elk.T,
    data = localN_z_1YrLag
  )
  summary(top_down_inter.a)
  AIC_psem(top_down_inter.a, AIC.type = "loglik")
  
  #'  ---------------------------------------
  #####  Top-down, exploitative competition  #####
  #'  ---------------------------------------
  top_down_exploit <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  ) 
  summary(top_down_exploit)
  AIC_psem(top_down_exploit, AIC.type = "loglik")
  
  #'  Reduced model
  #'  For all linear models at once-
  #'    1. Remove non-significant relationships that were initially hypothesized, rerun
  #'    2. Add newly significant relationships identified by d-sep test that were initially hypothesized, rerun one at a time
  #'    3. Add newly significant relationships identified by d-sep test that were not initially hypothesized but biologically plausible, rerun one at a time
  #'    4. Update unspecified correlation for newly significant relationships identified by d-sep test that are not expected to have a causal, unidirectional relationship, rerun
  top_down_exploit.a <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1, data = localN_z_1YrLag),
    #'  Update with unspecified correlated errors - focus on exogenous variables 
    #'  (year T variables) that are not presumed to have a causal & unidirectional
    #'  relationship but are related by an underlying factor making them appear 
    #'  correlated (e.g., similar habitat use)
    elk.T %~~% whitetailed_deer.T,
    # coyote.T %~~% mountain_lion.T,
    data = localN_z_1YrLag
  ) 
  summary(top_down_exploit.a)
  AIC_psem(top_down_exploit.a, AIC.type = "loglik")
  
  #' #'  Need to account for potentially missing relationships with mountain_lion.T, bear_black.T, and coyote.T
  #' #'  before top_down_exploit.a can be compared to top_down_inter.a (basically nested models)
  #' top_down_exploit.a.new <- update(top_down_exploit.a, 
  #'                                  mountain_lion.T ~ 1,
  #'                                  bear_black.T ~ 1,
  #'                                  coyote.T ~ 1)
  #' summary(top_down_exploit.a.new)
  #' AIC(top_down_inter.a, top_down_exploit.a.new)

  #'  ----------------------------------------
  #####  Bottom up, interference competition  #####
  #'  ----------------------------------------
  bottom_up_inter <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag), 
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1, data = localN_z_1YrLag), 
    data = localN_z_1YrLag
  )
  summary(bottom_up_inter)
  AIC_psem(bottom_up_inter, AIC.type = "loglik")
  
  #'  Reduced model
  #'  For all linear models at once-
  #'    1. Remove non-significant relationships that were initially hypothesized, rerun
  #'    2. Add newly significant relationships identified by d-sep test that were initially hypothesized, rerun one at a time
  #'    3. Add newly significant relationships identified by d-sep test that were not initially hypothesized but biologically plausible, rerun one at a time
  #'    4. Update unspecified correlation for newly significant relationships identified by d-sep test that are not expected to have a causal, unidirectional relationship, rerun
  bottom_up_inter.a <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1, data = localN_z_1YrLag), 
    lm(coyote.T ~ coyote.Tminus1, data = localN_z_1YrLag),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + bear_black.Tminus1 + mountain_lion.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1, data = localN_z_1YrLag), 
    #'  Update with unspecified correlated errors - focus on exogenous variables 
    #'  (year T variables) that are not presumed to have a causal & unidirectional
    #'  relationship but are related by an underlying factor making them appear 
    #'  correlated (e.g., similar habitat use)
    elk.T %~~% whitetailed_deer.T,
    coyote.T %~~% mountain_lion.T,
    whitetailed_deer.T %~~% bear_black.T,
    elk.T %~~% bear_black.T,
    data = localN_z_1YrLag
  )
  summary(bottom_up_inter.a)
  AIC_psem(bottom_up_inter.a, AIC.type = "loglik")
  
  #'  ----------------------------------------
  #####  Bottom up, exploitative competition  #####
  #'  ----------------------------------------
  bottom_up_exploit <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag), 
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(bottom_up_exploit)
  AIC_psem(bottom_up_exploit, AIC.type = "loglik")
  
  #'  Reduced model
  #'  For all linear models at once-
  #'    1. Remove non-significant relationships that were initially hypothesized, rerun
  #'    2. Add newly significant relationships identified by d-sep test that were initially hypothesized, rerun one at a time
  #'    3. Add newly significant relationships identified by d-sep test that were not initially hypothesized but biologically plausible, rerun one at a time
  #'    4. Update unspecified correlation for newly significant relationships identified by d-sep test that are not expected to have a causal, unidirectional relationship, rerun
  bottom_up_exploit.a <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1, data = localN_z_1YrLag), 
    lm(coyote.T ~ coyote.Tminus1, data = localN_z_1YrLag),
    #'  Update with unspecified correlated errors - focus on exogenous variables 
    #'  (year T variables) that are not presumed to have a causal & unidirectional
    #'  relationship but are related by an underlying factor making them appear 
    #'  correlated (e.g., similar habitat use)
    coyote.T %~~% mountain_lion.T, 
    data = localN_z_1YrLag
  )
  summary(bottom_up_exploit.a)
  AIC_psem(bottom_up_exploit.a, AIC.type = "loglik")
  
  #' #'  Need to account for potentially missing relationships with whitetailed_deer.T, elk.T, and moose.T
  #' #'  before bottom_up_exploit.a can be compared to bottom_up_inter.a (basically nested models)
  #' bottom_up_exploit.a.new <- update(bottom_up_exploit.a,
  #'                                   whitetailed_deer.T ~ 1,
  #'                                   elk.T ~ 1,
  #'                                   moose.T ~ 1)
  #' summary(top_down_exploit.a.new)
  #' AIC(bottom_up_inter.a, bottom_up_exploit.a.new)
  
  
  #'  --------------------------------------------------------------
  #####  Bottom up, exploitative competition via white-tailed deer  #####
  #'  --------------------------------------------------------------
  #'  White-tailed deer economy
  bottom_up_exploit_wtd <- psem(
    lm(wolf.T ~ wolf.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag), 
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    data = localN_z_1YrLag
  )
  summary(bottom_up_exploit_wtd)
  AIC_psem(bottom_up_exploit_wtd, AIC.type = "loglik")
  
  bottom_up_exploit_wtd.a <- psem(
    lm(mountain_lion.T ~ mountain_lion.Tminus1, data = localN_z_1YrLag),
    lm(bear_black.T ~ bear_black.Tminus1, data = localN_z_1YrLag), 
    lm(coyote.T ~ coyote.Tminus1, data = localN_z_1YrLag),
    coyote.T %~~% mountain_lion.T,
    data = localN_z_1YrLag
  )
  summary(bottom_up_exploit_wtd.a) #.... not very informative
  
  #'  -----------------------------
  ####  Model selection using AIC  ####
  #'  -----------------------------
  #'  Include all permutations of these models (original & reduced via dSep)
  modSelect <- AIC(top_down_inter, top_down_inter.a, top_down_exploit, top_down_exploit.a, 
                   bottom_up_inter, bottom_up_inter.a, bottom_up_exploit, bottom_up_exploit.a, 
                   bottom_up_exploit_wtd, bottom_up_exploit_wtd.a)
  mod_names <- c("top_down_inter", "top_down_inter_reduced", "top_down_exploit", "top_down_exploit_reduced", 
                 "bottom_up_inter", "bottom_up_inter_reduced", "bottom_up_exploit", "bottom_up_exploit_reduced", 
                 "bottom_up_exploit_wtd", "bottom_up_exploit_wtd_reduced")
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("model", "AIC", "K", "n")
  arrange(modSelect, AIC, decreasing = TRUE)
  
  #'  Just the reduced models 
  modSelect <- AIC(top_down_inter.a, top_down_exploit.a, bottom_up_inter.a, 
                   bottom_up_exploit.a, bottom_up_exploit_wtd.a)
  mod_names <- c("top_down_inter_reduced", "top_down_exploit_reduced", "bottom_up_inter_reduced", 
                 "bottom_up_exploit_reduced", "bottom_up_exploit_wtd_reduced")
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("model", "AIC", "K", "n")
  arrange(modSelect, AIC, decreasing = TRUE)
  
  #'  Review Fisher's C and Chi-square test statistics for each model
  LLchisq(top_down_inter.a); fisherC(top_down_inter.a)
  LLchisq(top_down_exploit.a); fisherC(top_down_exploit.a)
  LLchisq(bottom_up_inter.a); fisherC(bottom_up_inter.a)
  LLchisq(bottom_up_exploit.a); fisherC(bottom_up_exploit.a)
  LLchisq(bottom_up_exploit_wtd.a); fisherC(bottom_up_exploit_wtd.a)
  
  #'  Dig into top model based on AIC, Fisher's C
  getDAG(top_down_exploit.a)
  residuals(top_down_exploit.a)
  