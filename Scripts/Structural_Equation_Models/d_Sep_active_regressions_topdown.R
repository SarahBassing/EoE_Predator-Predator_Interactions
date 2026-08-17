  #'  ----------------------------------------------------------------
  #'  Active regressions for d-Sep iterations: top-down
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
  
  dSep_iterations_topdown <- list(
    #'  Conditional independence claims
    # Regression #1
    # bs_topdown[[1]] # "deerHarv.tmin1" "coy.t"          "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "deerHarv"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #2
    # bs_topdown[[2]] # "deerHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "deerHarv"), spp = c(".harvest", ".bear", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #3
    # bs_topdown[[3]] # "deerHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "deerHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #4
    # bs_topdown[[4]] # "deerHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "deerHarv"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #5
    # bs_topdown[[5]] # "deerHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "deerHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #6
    # bs_topdown[[6]] # "deerHarv.tmin1" "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "deerHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # Regression #7
    # bs_topdown[[7]] # "wtd.tmin1" "coy.t"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent"), spp = c(".coy", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #8
    # bs_topdown[[8]] # "wtd.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "wtd.latent"), spp = c(".harvest", ".bear", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #9
    # bs_topdown[[9]] # "wtd.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wtd.latent"), spp = c(".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #10
    # bs_topdown[[10]] # "wtd.tmin1"     "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #11
    # bs_topdown[[11]] # "wtd.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #12
    # bs_topdown[[12]] # "wtd.tmin1"      "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wtd.latent"), spp = c(".harvest", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #13
    # bs_topdown[[13]] # "moose.tmin1"    "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "moose.latent"), spp = c(".harvest", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #14
    # bs_topdown[[14]] # "moose.tmin1" "coy.t"       "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "moose.latent"), spp = c(".coy", ".moose"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #15
    # bs_topdown[[15]] # "moose.tmin1"    "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "moose.latent"), spp = c(".harvest", ".bear", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #16
    # bs_topdown[[16]] # "moose.tmin1"   "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "moose.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #17
    # bs_topdown[[17]] # "moose.tmin1"    "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "moose.latent"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #18
    # bs_topdown[[18]] # "moose.tmin1"    "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "moose.latent"), spp = c(".harvest",  ".moose"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #19
    # bs_topdown[[19]] # "elkHarv.tmin1"  "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elkHarv"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #20
    # bs_topdown[[20]] # "elkHarv.tmin1" "coy.t"         "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "elkHarv"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #21
    # bs_topdown[[21]] # "elkHarv.tmin1" "moose.t"       "moose.tmin1"   "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elkHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #22
    # bs_topdown[[22]] # "elkHarv.tmin1" "moose.t"       "moose.tmin1"   "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elkHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #23
    # bs_topdown[[23]] # "elkHarv.tmin1"  "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elkHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #24
    # bs_topdown[[24]] # "elkHarv.tmin1"  "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv",  "elkHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # Regression #25
    # bs_topdown[[25]] # "lion.tmin1" "coy.t"      "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "lion.latent"), spp = c(".coy", ".lion"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #26
    # bs_topdown[[26]] # "lion.tmin1"     "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "lion.latent"), spp = c(".harvest", ".bear", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #27
    # bs_topdown[[27]] # "lion.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "lion.latent"), spp = c(".moose", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #28
    # bs_topdown[[28]] # "lion.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "lion.latent"), spp = c(".harvest", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #29
    # bs_topdown[[29]] # "lion.tmin1"     "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "lion.latent"), spp = c(".harvest", ".lion"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    
    ##########
    # Regression #30   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[30]] # "wtd.t"          "elk.tmin1"      "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.latent", "wtd.latent", "lion.latent", "elk.latent"), spp = c(".harvest", ".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #31   #  NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[31]] # "wtd.t"          "coy.tmin1"      "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.latent", "wtd.latent", "lion.latent", "coy.latent"), spp = c(".harvest", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    ##########
    
    # Regression #32
    # bs_topdown[[32]] # "wtd.t"          "coy.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("deerHarv", "wtd.latent", "lion.latent", "coy.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".coy", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    
    ##########
    # Regression #33   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[33]] # "wtd.t"          "bearHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.latent", "wtd.latent", "lion.latent", "bearHarv.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y-1")),
    # Regression #NA   # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS TEMPORALLY VALID
    # bs_topdown[[34]] # "wtd.t"          "bear.tmin1"     "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.latent", "wtd.latent", "lion.latent", "bear.latent"), spp = c(".harvest", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    ##########
    
    # Regression #31
    # bs_topdown[[35]] # "wtd.t"          "bear.t"         "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("deerHarv", "wtd.latent", "lion.latent", "bearHarv", "bear.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".bear", ".wtd"), indices = as.integer(c(1,1,1,2,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[36]] # "wtd.t"          "wolfHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    # # Regression #NA   NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[37]] # "wtd.t"          "wolf.tmin1"     "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #32
    # bs_topdown[[38]] # "wtd.t"          "moose.t"        "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("deerHarv", "wtd.latent", "lion.latent", "moose.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #33
    # bs_topdown[[39]] # "wtd.t"          "elk.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("deerHarv", "wtd.latent", "lion.latent", "elkHarv", "elk.latent", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #34
    # bs_topdown[[40]] # "wtd.t"          "wolf.t"         "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("deerHarv", "wtd.latent", "lion.latent", "wolfHarv", "wolf.latent", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[41]] #  "wtd.t"          "lionHarv.tmin1" "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    # list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #35
    # bs_topdown[[42]] # "wtd.t"          "lion.t"         "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("deerHarv", "wtd.latent", "lion.latent", "lionHarv", "wtd.latent"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".wtd"), indices = as.integer(c(1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #36
    # bs_topdown[[43]] # "elk.tmin1" "coy.t"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.latent", "elk.latent"), spp = c(".coy", ".elk"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #37
    # bs_topdown[[44]] # "elk.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "elk.latent"), spp = c(".harvest", ".bear", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #38
    # bs_topdown[[45]] # "elk.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "elk.latent"), spp = c(".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #39
    # bs_topdown[[46]] # "elk.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "elk.latent"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #40
    # bs_topdown[[47]] # "elk.tmin1"      "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "elk.latent"), spp = c(".harvest", ".elk"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #41
    # bs_topdown[[48]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv", "bear.latent", "coy.latent"), spp = c(".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #42
    # bs_topdown[[49]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "coy.latent"), spp = c(".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #43
    # bs_topdown[[50]] # "coy.tmin1"     "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "coy.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #44
    # bs_topdown[[51]] # "coy.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "coy.latent"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #45
    # bs_topdown[[52]] # "coy.tmin1"      "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "coy.latent"), spp = c(".harvest", ".coy"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[53]] # "coy.t"          "bearHarv.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[54]] # "coy.t"      "bear.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #46
    # bs_topdown[[55]] # "coy.t"          "bear.t"         "coy.tmin1"      "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("coy.latent", "bearHarv", "bear.latent", "coy.latent"), spp = c(".coy", ".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[56]] # "coy.t"          "wolfHarv.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[57]] # "coy.t"      "wolf.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #47
    # bs_topdown[[58]] # "coy.t"       "moose.t"     "coy.tmin1"   "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("coy.latent", "moose.latent", "wolf.latent", "coy.latent"), spp = c(".coy", ".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #48
    # bs_topdown[[59]] #  "coy.t"         "elk.t"         "coy.tmin1"     "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("coy.latent", "elkHarv", "lion.latent", "elk.latent", "wolf.latent", "coy.latent"), spp = c(".coy", ".harvest", ".lion", ".elk", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #49
    # bs_topdown[[60]] # "coy.t"          "wolf.t"         "coy.tmin1"      "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("coy.latent", "wolfHarv", "wolf.latent", "coy.latent"), spp = c(".coy", ".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[61]] # "coy.t"          "lionHarv.tmin1" "coy.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #50
    # bs_topdown[[62]] # "coy.t"          "lion.t"         "coy.tmin1"      "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("coy.latent", "lionHarv", "coy.latent"), spp = c(".coy", ".harvest", ".coy"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y")),
    # Regression #51
    # bs_topdown[[63]] #  "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bearHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #52
    # bs_topdown[[64]] #  "bearHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "bearHarv"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #53
    # bs_topdown[[65]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bearHarv"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y-1")),
    # Regression #54
    # bs_topdown[[66]] # "bearHarv.tmin1" "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "bearHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # Regression #55
    # bs_topdown[[67]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "bear.latent"), spp = c(".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #56
    # bs_topdown[[68]] # "bear.tmin1"    "elk.t"         "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"     "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "bear.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #57
    # bs_topdown[[69]] # "bear.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #58
    # bs_topdown[[70]] # "bear.tmin1"     "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "bear.latent"), spp = c(".harvest", ".bear"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[71]] # "bear.t"         "wolfHarv.tmin1" "bearHarv.tmin1" "bear.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[72]] # "bear.t"         "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #59
    # bs_topdown[[73]] # "bear.t"         "moose.t"        "bearHarv.tmin1" "bear.tmin1"     "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("bearHarv", "bear.latent", "moose.latent", "wolf.latent", "bear.latent"), spp = c(".harvest", ".bear", ".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    # Regression #60
    # bs_topdown[[74]] # "bear.t"         "elk.t"          "bearHarv.tmin1" "bear.tmin1"     "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("bearHarv", "bear.latent", "elkHarv", "lion.latent", "elk.latent", "wolf.latent", "bear.latent"), spp = c(".harvest", ".bear", ".harvest", ".lion", ".elk", ".wolf", ".bear"), indices = as.integer(c(1,1,2,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #61
    # bs_topdown[[75]] # "bear.t"         "wolf.t"         "bearHarv.tmin1" "bear.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("bearHarv", "bear.latent", "wolfHarv", "wolf.latent", "bear.latent"), spp = c(".harvest", ".bear", ".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,2,1,2)), lags = c("y-1","y-1","y-1","y-1","y")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[76]] # "bear.t"         "lionHarv.tmin1" "bearHarv.tmin1" "bear.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #62
    # bs_topdown[[77]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("bearHarv", "bear.latent", "lionHarv", "bear.latent"), spp = c(".harvest", ".bear", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #63
    # bs_topdown[[78]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.latent", "wolf.latent", "wolfHarv"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    # Regression #64
    # bs_topdown[[79]] # "wolfHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "wolfHarv"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #65
    # bs_topdown[[80]] #  "wolfHarv.tmin1" "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolfHarv"), spp = c(".harvest", ".harvest"), indices = as.integer(c(1,2)), lags = c("y-1","y-1")),
    # Regression #66
    # bs_topdown[[81]] #  "wolf.tmin1"     "lion.t"         "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv", "wolf.latent"), spp = c(".harvest", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    # Regression #67
    # bs_topdown[[82]] # "moose.t"       "elk.t"         "moose.tmin1"   "wolf.tmin1"    "elkHarv.tmin1" "lion.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("moose.latent", "wolf.latent", "elkHarv", "lion.latent", "elk.latent", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".lion", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #68
    # bs_topdown[[83]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("moose.latent", "wolf.latent", "wolfHarv", "moose.latent"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[84]] # "moose.t"        "lionHarv.tmin1" "moose.tmin1"    "wolf.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #69
    # bs_topdown[[85]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "wolf.latent", "lionHarv", "moose.latent"), spp = c(".moose", ".wtd", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #70
    # bs_topdown[[86]] # "elk.t"          "wolf.t"         "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "wolfHarv", "elk.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    
    ##########
    # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[87]] # "elk.t"          "lionHarv.tmin1" "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c())),
    ##########
    
    # Regression #71
    # bs_topdown[[88]] # "elk.t"          "lion.t"         "elkHarv.tmin1"  "lion.tmin1"     "elk.tmin1"      "wolf.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("elkHarv", "lion.latent", "elk.latent", "wolf.latent", "lionHarv", "elk.latent"), spp = c(".harvest", ".lion", ".elk", ".wolf", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y"))
    
    ##########
    # # # Regression #NA   #  NOT POSSIBLE - t cannot affect t-1
    # # bs_topdown[[89]] # "wolf.t"         "lionHarv.tmin1" "wolfHarv.tmin1" "wolf.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()))
    ##########
    
    #'  Unconditional independence claims
    
  )