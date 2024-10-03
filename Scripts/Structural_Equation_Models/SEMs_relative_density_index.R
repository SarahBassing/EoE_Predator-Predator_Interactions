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
  
  #'  -------------
  #####  Top-down  #####
  #'  -------------
  top_down_exploit <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = localN_z_1YrLag),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = localN_z_1YrLag),
    # lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    # lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    # lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = localN_z_1YrLag),
    # lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1 + mountain_lion.Tminus1, data = localN_z_1YrLag),
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
    # lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = localN_z_1YrLag),
    # lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = localN_z_1YrLag),
    # lm(bear_black.T ~ bear_black.Tminus1, data = localN_z_1YrLag),
    # lm(coyote.T ~ coyote.Tminus1, data = localN_z_1YrLag),
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
  
  #'  -----------------------------
  ####  Model selection using AIC  ####
  #'  -----------------------------
  #'  Include all permutations of these models (original & reduced via dSep)
  modSelect <- AIC(top_down_inter, top_down_inter.a, top_down_exploit, top_down_exploit.a, 
                   bottom_up_inter, bottom_up_inter.a, bottom_up_exploit, bottom_up_exploit.a)
  mod_names <- c("top_down_inter", "top_down_inter_reduced", "top_down_exploit", "top_down_exploit_reduced", 
                 "bottom_up_inter", "bottom_up_inter_reduced", "bottom_up_exploit", "bottom_up_exploit_reduced")
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("model", "AIC", "K", "n")
  arrange(modSelect, AIC, decreasing = TRUE)
  
  #'  Just the reduced models 
  modSelect <- AIC(top_down_inter.a, top_down_exploit.a, bottom_up_inter.a, 
                   bottom_up_exploit.a)
  mod_names <- c("top_down_inter_reduced", "top_down_exploit_reduced", 
                 "bottom_up_inter_reduced", "bottom_up_exploit_reduced")
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("model", "AIC", "K", "n")
  arrange(modSelect, AIC, decreasing = TRUE)
  
  #'  Review Fisher's C and Chi-square test statistics for each model
  LLchisq(top_down_inter.a); fisherC(top_down_inter.a)
  LLchisq(top_down_exploit.a); fisherC(top_down_exploit.a)
  LLchisq(bottom_up_inter.a); fisherC(bottom_up_inter.a)
  LLchisq(bottom_up_exploit.a); fisherC(bottom_up_exploit.a)
  
  #'  ----------------------
  ####  Dig into top model  ####
  #'  ----------------------
  #'  Dig into top model based on AIC, Fisher's C
  #'  Check for multicollinearity
  RVIF(top_down_exploit.a[[1]]) # white-tailed deer model
  RVIF(top_down_exploit.a[[2]]) # elk model
  RVIF(top_down_exploit.a[[3]]) # moose model
  
  #'  Run some basic model diagnostics on individuals models
  #'  This is handy: http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/
  wtd_mod <- lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = localN_z_1YrLag)
  plot(wtd_mod)
  elk_mod <- lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = localN_z_1YrLag)
  plot(elk_mod)
  moose_mod <- lm(moose.T ~ moose.Tminus1, data = localN_z_1YrLag)
  plot(moose_mod)
  
  #'  Visualize SEM
  piecewiseSEM:::plot.psem(top_down_exploit.a, 
                           node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "orange", fontcolor = "black"), 
                           layout = "tree")
  piecewiseSEM:::plot.psem(top_down_exploit.a, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  
  piecewiseSEM:::plot.psem(top_down_inter.a, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(bottom_up_inter.a, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(bottom_up_exploit.a, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  
  #'  Calculate direct, indirect, total, and mediator effects (SE & 95% CI) for 
  #'  all endogenous (response) variables using semEff package
  #'  https://murphymv.github.io/semEff/articles/semEff.html
  #'  First bootstrap standardized model coefficients (necessary for calculating SEs)
  #'  THIS TAKES AWHILE! 
  top_down_exploit.a_bootEff <- bootEff(top_down_exploit.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  save(top_down_exploit.a_bootEff, file = paste0("./Outputs/SEM/top_down_exploit.a_bootEff_", Sys.Date(), ".RData"))
  
  top_down_inter.a_bootEff <- bootEff(top_down_inter.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  bottom_up_inter.a_bootEff <- bootEff(bottom_up_inter.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  bottom_up_exploit.a_bootEff <- bootEff(bottom_up_exploit.a, R = 1000, seed = 13, type = "nonparametric", parallel = "multicore", ncpus = 5) 
  
  #'  Second calculate standardized effects for all casual pathways
  #'  Note, these are standardized unique effects (i.e., adjusted for multicollinearity;
  #'  i.e., semipartial correlations), allowing us to fully partition effects in the system
  #'  These tend to be smaller than the unadjusted standardized coefficients
  #'  If there are large differences btwn the two then consideration should be given
  #'  to the impact and relevance of multicollinearity in the system (check with RVIF())
  (top_down_exploit.a_semEff <- semEff(top_down_exploit.a_bootEff))
  summary(top_down_exploit.a_semEff)
  save(top_down_exploit.a_semEff, file = paste0("./Outputs/SEM/top_down_exploit.a_semEff_", Sys.Date(), ".RData"))
  
  (top_down_inter.a_semEff <- semEff(top_down_inter.a_bootEff))
  summary(top_down_inter.a_semEff)
  (bottom_up_inter.a_semEff <- semEff(bottom_up_inter.a_bootEff))
  summary(bottom_up_inter.a_semEff)
  (bottom_up_exploit.a_semEff <- semEff(bottom_up_exploit.a_bootEff))
  summary(bottom_up_exploit.a_semEff)
  
  
  
  
  
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
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + wolf.yr1 + mountain_lion.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + wolf.yr2 + mountain_lion.yr2, data = localN_z),
    data = localN_z
  )
  summary(top_down_inter_yr)
  AIC_psem(top_down_inter_yr, AIC.type = "loglik")
  
  top_down_inter_yr.r <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
    moose.yr3 %~~% elk.yr1,
    elk.yr3 %~~% moose.yr1,
    coyote.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr3,
    elk.yr2 %~~% whitetailed_deer.yr2,
    data = localN_z
  )
  summary(top_down_inter_yr.r)     
  AIC_psem(top_down_inter_yr.r, AIC.type = "loglik")
  
  
  top_down_exploit_yr <- psem(
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1 + coyote.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2 + coyote.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1 + mountain_lion.yr1 + bear_black.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2 + mountain_lion.yr2 + bear_black.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1 + wolf.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2 + wolf.yr2, data = localN_z),
    lm(wolf.yr2 ~ wolf.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2, data = localN_z),
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
    lm(elk.yr3 ~ elk.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + elk.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + elk.yr2, data = localN_z),
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
    mountain_lion.yr3 %~~% elk.yr1,
    data = localN_z
  )
  summary(top_down_exploit_yr.r)
  AIC_psem(top_down_exploit_yr.r, AIC.type = "loglik")
  
  
  
  bottom_up_inter_yr <- psem(
    lm(wolf.yr2 ~ wolf.yr1 + moose.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(wolf.yr3 ~ wolf.yr2 + moose.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1 + wolf.yr1 + bear_black.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + bear_black.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1 + wolf.yr1 + elk.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2 + wolf.yr2 + elk.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(coyote.yr2 ~ coyote.yr1 + wolf.yr1 + mountain_lion.yr1 + whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2 + wolf.yr2 + mountain_lion.yr2 + whitetailed_deer.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ elk.yr1 + wolf.yr1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    data = localN_z
  )
  summary(bottom_up_inter_yr)
  AIC_psem(bottom_up_inter_yr, AIC.type = "loglik")
  
  bottom_up_inter_yr <- psem(
    lm(wolf.yr2 ~ moose.yr1, data = localN_z),
    lm(wolf.yr3 ~ 1, data = localN_z),
    lm(mountain_lion.yr2 ~ mountain_lion.yr1, data = localN_z),
    lm(mountain_lion.yr3 ~ mountain_lion.yr2 + wolf.yr2 + elk.yr2, data = localN_z),
    lm(bear_black.yr2 ~ bear_black.yr1, data = localN_z),
    lm(bear_black.yr3 ~ bear_black.yr2, data = localN_z),
    lm(coyote.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(coyote.yr3 ~ coyote.yr2, data = localN_z),
    lm(whitetailed_deer.yr2 ~ whitetailed_deer.yr1, data = localN_z),
    lm(whitetailed_deer.yr3 ~ whitetailed_deer.yr2, data = localN_z),
    lm(elk.yr2 ~ 1, data = localN_z),
    lm(elk.yr3 ~ elk.yr2 + wolf.yr2, data = localN_z),
    lm(moose.yr2 ~ moose.yr1, data = localN_z),
    lm(moose.yr3 ~ moose.yr2, data = localN_z),
    # moose.yr3 %~~% moose.yr1,
    # bear_black.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr1,
    coyote.yr3 %~~% bear_black.yr3,
    mountain_lion.yr3 %~~% whitetailed_deer.yr1,
    mountain_lion.yr2 %~~% wolf.yr2,
    data = localN_z
  )
  summary(bottom_up_inter_yr)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  