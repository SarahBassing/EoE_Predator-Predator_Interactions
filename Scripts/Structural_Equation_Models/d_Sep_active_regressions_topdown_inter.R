  #'  ----------------------------------------------------------------
  #'  Active regressions for d-Sep iterations: top-down interference
  #'  Aug 2026
  #'  ----------------------------------------------------------------
  #'  List active regressions to iterate over for each d-separation test. 
  #'  Necessary objects for each active regression include:
  #'  
  #'  -- dSep_test: integer that indicates which regression (1:7) within the SEM 
  #'     gets the custom regression where: 
  #'            1 = cougar 
  #'            2 = wolf 
  #'            3 = black bear 
  #'            4 = coyote 
  #'            5 = elk 
  #'            6 = moose
  #'            7 = wtd 
  #'  -- covariates: vector that include the covariates included in the SEM 
  #'     a priori, plus the one needed for the conditional independence test.
  #'     NOTE: the last covariate in this vector should always be the explanatory
  #'     variable you are testing the independence claim on. This will make life
  #'     easier when using the p.rope() and bayes_pvalue() functions in 
  #'     d_Sep_test_and_FishersC_for_Bayesian_SEM.R
  #'     Note: need to include the AR1 term in covariates list for each regression -
  #'     basic_set seems to drop it b/c it can't distinguish from the response 
  #'     variable with the [nSite, nYear] indexing. Exceptions noted for specific
  #'     regressions (lion regressions).
  #'  -- spp: vector containing the suffix of each parameter name to be appended
  #'     to "beta" in the model template. Must include a period before the term's
  #'     name to work.
  #'  -- indices: vector containing the index value for each parameter. For most
  #'     regressions, each parameter's name is unique so the index is 1. But for some,
  #'     the parameter name is used twice in the model representing different coefficients
  #'     (e.g., beta.harvest[1] and beta.harvest[2] correspond to the effects of 
  #'     wolfHarv and lionHarv in the same regression). Double check indexing for 
  #'     each regression to ensure appropriate coefficient estimates are saved and 
  #'     used in ROPE / Bayesian p-value d-Sep tests.
  #'  ---------------------------------------------------------------- 

  dSep_iterations_topdown_int <- list(
    # Regression #1
    # bs_topdown_inter[[1]] # "wtd.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #2 
    # bs_topdown_inter[[2]] # "wtd.tmin1"  "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wtd.latent"), spp = c(".elk", ".lion", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #3 
    # bs_topdown_inter[[3]] # "wtd.tmin1"  "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".lion", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #4 
    # bs_topdown_inter[[4]] # "wtd.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #5 
    # bs_topdown_inter[[5]] # "wtd.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #6 
    # bs_topdown_inter[[6]] # "wtd.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #7 
    # bs_topdown_inter[[7]] # "moose.tmin1" "wtd.t"       "wtd.tmin1"   "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "moose.latent"), spp = c(".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #8 
    # bs_topdown_inter[[8]] # "moose.tmin1" "elk.t"       "elk.tmin1"   "lion.tmin1"  "wolf.tmin1"
    list(dSep_test = 5, covaraites = c("elk.latent", "lion.latent", "wolf.latent", "moose.latent"), spp = c(".elkt", ".lion", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #9 
    # bs_topdown_inter[[9]] # "moose.tmin1" "coy.t"       "lion.tmin1"  "coy.tmin1"   "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "moose.latent"), spp = c(".lion", ".coy", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #10 
    # bs_topdown_inter[[10]] # "moose.tmin1"    "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "moose.latent"), spp = c(".harvest", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #11 
    # bs_topdown_inter[[11]] # "moose.tmin1"    "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #12 
    # bs_topdown_inter[[12]] # "moose.tmin1"    "lion.t"         "lionHarv.tmin1" "wolf.tmin1" 
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #13 
    # bs_topdown_inter[[13]] # "elk.tmin1"  "wtd.t"      "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.latent", "lion.latent", "elk.latent"), spp = c(".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #14 
    # bs_topdown_inter[[14]] # "elk.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elk.latent"), spp = c(".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #15 
    # bs_topdown_inter[[15]] "elk.tmin1"  "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "elk.latent"), spp = c(".lion", ".coy", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #16
    # bs_topdown_inter[[16]] # "elk.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1" 
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "elk.latent"), spp = c(".harvest", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #17 
    # bs_topdown_inter[[17]] # "elk.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #18 
    # bs_topdown_inter[[18]] # "elk.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #19 
    # bs_topdown_inter[[19]] # "lion.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lion.latent"), spp = c(".moose", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #20 
    # bs_topdown_inter[[20]] # "lion.tmin1"     "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1" 
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lion.latent"), spp = c(".harvest", ".bear", ".wolf", ".lion"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #21 
    # bs_topdown_inter[[21]] # "lion.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"    
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #22 
    # bs_topdown_inter[[22]] # "lion.tmin1"     "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # # Regression #NA         #  NOT POSSIBLE - t cannot affect t-1 
    # # bs_topdown_inter[[23]] # "wtd.t"      "coy.tmin1"  "wtd.tmin1"  "lion.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # # Regression #NA         #  NOT POSSIBLE - t cannot affect t-1  
    # # bs_topdown_inter[[24]] # "wtd.t"          "bearHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # # Regression #NA         #  NOT POSSIBLE - t cannot affect t-1  
    # # bs_topdown_inter[[25]] # "wtd.t"      "bear.tmin1" "wtd.tmin1"  "lion.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # # Regression #NA         #  NOT POSSIBLE - t cannot affect t-1   
    # # bs_topdown_inter[[26]] # "wtd.t"          "wolfHarv.tmin1" "wtd.tmin1"      "lion.tmin1" 
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()), 
    # # Regression #NA         #  NOT POSSIBLE - t cannot affect t-1   
    # # bs_topdown_inter[[27]] # "wtd.t"          "lionHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # # Regression #NA         #  NOT POSSIBLE - t cannot affect t-1   
    # # bs_topdown_inter[[28]] # "wtd.t"      "wolf.tmin1" "wtd.tmin1"  "lion.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #23 
    # bs_topdown_inter[[29]] # "wtd.t"       "moose.t"     "wtd.tmin1"   "lion.tmin1"  "moose.tmin1" "wolf.tmin1" 
    list(dSep_test = 6, covariates = c("wtd.latent", "lion.latent", "moose.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #24 
    # bs_topdown_inter[[30]] # "wtd.t"      "elk.t"      "wtd.tmin1"  "lion.tmin1" "elk.tmin1"  "wolf.tmin1"
    list(dSep_test = 5, covariates = c("wtd.latent", "lion.latent", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regerssion #25 
    # bs_topdown_inter[[31]] # "wtd.t"      "coy.t"      "wtd.tmin1"  "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("wtd.latent", "lion.latent", "coy.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #26 
    # bs_topdown_inter[[32]] # "wtd.t"          "bear.t"         "wtd.tmin1"      "lion.tmin1"     "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("wtd.latent", "lion.latent", "bearHarv", "bear.latent", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".harvest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #27
    # bs_topdown_inter[[33]] # "wtd.t"          "wolf.t"         "wtd.tmin1"      "lion.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wtd.latent", "lion.latent", "wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #28 
    # bs_topdown_inter[[34]] # "wtd.t"          "lion.t"         "wtd.tmin1"      "lion.tmin1"     "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("wtd.latent", "lion.latent", "lionHarv", "wolf.latent", "wtd.latent"), spp = c(".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #29
    # bs_topdown_inter[[35]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #30
    # bs_topdown_inter[[36]] # "coy.tmin1"  "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "coy.latent"), spp = c(".elk", ".lion", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #31 
    # bs_topdown_inter[[37]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv",  "bear.latent", "wolf.latent", "coy.latent"), spp = c(".harvest", ".bear", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #32 
    # bs_topdown_inter[[38]] # "coy.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #33 
    # bs_topdown_inter[[39]] # "coy.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1" 
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #34
    # bs_topdown_inter[[40]] # "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1" 
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bearHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #35
    # bs_topdown_inter[[41]] # "bearHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"   
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bearHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #36 
    # bs_topdown_inter[[42]] # "bearHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #37
    # bs_topdown_inter[[43]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"  
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #38 
    # bs_topdown_inter[[44]] # "bearHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"    
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #39
    # bs_topdown_inter[[45]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #40 
    # bs_topdown_inter[[46]] # "bear.tmin1" "elk.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bear.latent"), spp = c(".elk", ".lion", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #41 
    # bs_topdown_inter[[47]] # "bear.tmin1" "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bear.latent"), spp = c(".lion", ".coy", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #42 
    # bs_topdown_inter[[48]] # "bear.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #43 
    # bs_topdown_inter[[49]] # "bear.tmin1"     "lion.t"         "lionHarv.tmin1" "wolf.tmin1"   
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #44 
    # bs_topdown_inter[[50]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wolfHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #45 
    # bs_topdown_inter[[51]] # "wolfHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wolfHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #46 
    # bs_topdown_inter[[52]] # "wolfHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"    
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #47 
    # bs_topdown_inter[[53]] # "wolfHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #48 
    # bs_topdown_inter[[54]] # "wolfHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #49 
    # bs_topdown_inter[[55]] # "lionHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lionHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #50 
    # bs_topdown_inter[[56]] # "lionHarv.tmin1" "elk.t"          "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elk.latent", "lion.latent", "wolf.latent", "lionHarv"), spp = c(".elk", ".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #51 
    # bs_topdown_inter[[57]] # "lionHarv.tmin1" "coy.t"          "lion.tmin1"     "coy.tmin1"      "wolf.tmin1" 
    list(dSep_test = 4, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv"), spp = c(".lion", ".coy", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #52 
    # bs_topdown_inter[[58]] # "lionHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #53 
    # bs_topdown_inter[[59]] # "lionHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"    
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lionHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #54 
    # bs_topdown_inter[[60]] # "moose.t"     "elk.t"       "moose.tmin1" "wolf.tmin1"  "elk.tmin1"   "lion.tmin1" 
    list(dSep_test = 5, covariates = c("moose.latent", "wolf.latent", "elk.latent", "lion.latent", "moose.latent"), spp = c(".moose", ".wolf", ".elk", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #55 
    # bs_topdown_inter[[61]] # "moose.t"     "coy.t"       "moose.tmin1" "wolf.tmin1"  "lion.tmin1"  "coy.tmin1"
    list(dSep_test = 4, covariates = c("moose.latent", "wolf.latent", "lion.latent", "coy.latent", "moose.latent"), spp = c(".moose", ".wolf", ".lion", ".coy", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #56 
    # bs_topdown_inter[[62]] # "moose.t"        "bear.t"         "moose.tmin1"    "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("moose.latent", "wolf.latent", "bearHarv", "bear.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".bear", ".moose"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #57 
    # bs_topdown_inter[[63]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolf.latent", "wolfHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #58 
    # bs_topdown_inter[[64]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "wolf.latent", "lionHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #59 
    # bs_topdown_inter[[65]] # "elk.t"      "coy.t"      "elk.tmin1"  "lion.tmin1" "wolf.tmin1" "coy.tmin1" 
    list(dSep_test = 4, covariates = c("elk.latent", "lion.latent", "wolf.latent", "coy.latent", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".coy", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #60 
    # bs_topdown_inter[[66]] # "elk.t"          "bear.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"    
    list(dSep_test = 3, covariates = c("elk.latent", "lion.latent", "wolf.latent", "bearHarv", "bear.latent", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".bear", ".elk"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #61 
    # bs_topdown_inter[[67]] # "elk.t"          "wolf.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("elk.latent", "lion.latent", "wolf.latent", "wolfHarv", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #62 
    # bs_topdown_inter[[68]] # "elk.t"          "lion.t"         "elk.tmin1"      "lion.tmin1"     "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("elk.latent", "lion.latent", "wolf.latent", "lionHarv", "elk.latent"), spp = c(".elk", ".lion", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #63 
    # bs_topdown_inter[[69]] # "coy.t"          "bear.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("lion.latent", "coy.latent", "wolf.latent", "bearHarv", "bear.latent", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #64 
    # bs_topdown_inter[[70]] # "coy.t"          "wolf.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("lion.latent", "coy.latent", "wolf.latent", "wolfHarv", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".coy"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #65 
    # bs_topdown_inter[[71]] # "coy.t"          "lion.t"         "lion.tmin1"     "coy.tmin1"      "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lion.latent", "coy.latent", "wolf.latent", "lionHarv", "coy.latent"), spp = c(".lion", ".coy", ".wolf", ".harvest", ".coy"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #66 
    # bs_topdown_inter[[72]] # "bear.t"         "wolf.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("bearHarv", "bear.latent", "wolf.latent", "wolfHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #67
    # bs_topdown_inter[[73]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("bearHarv", "bear.latent", "wolf.latent", "lionHarv", "bear.latent"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y"))
  )
  