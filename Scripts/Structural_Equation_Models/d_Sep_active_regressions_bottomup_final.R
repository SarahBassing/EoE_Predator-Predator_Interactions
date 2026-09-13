  #'  ----------------------------------------------------------------
  #'  Active regression for d-Sep iterations: bottom-up exploitative FINAL
  #'  Sept 2026
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
  #'     NOTE: in cases where independence claim is testing whether variable A at 
  #'     time t affects variable B at time t-1 (which is not possible given our
  #'     current understanding of space and time), I switched which variable was
  #'     the response vs explanatory variable. This effectively tests whether 
  #'     there is a reasonable correlation between the variables. Although the 
  #'     exact independence claim cannot be included in the SEM (due to how time
  #'     works), it allows us to include the independence claim when calculating
  #'     GoF. This is important for accurately assessing whether the hypothesized 
  #'     causal model fits the data.
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
  #'  -- lags: vector containing character strings that indicate whether the year
  #'     index should be lagged by one year. "y-1" applies a lag to the variable;
  #'     "y" does not apply a lag.
  #'  
  #'  Some independence claims cannot be assess with this script given how the
  #'  SEM regressions are written. Those are included here for consistency but 
  #'  are commented out and assessed in other scripts: 
  #'  d_Sep_active_regression_bottomup_tmin1_only_final.R  and 
  #'  d_Sep_active_regression_bottomup_exog_only_final.R
  #'  ---------------------------------------------------------------- 
  
  dSep_iterations_bottomup_final <- list(
    # #########
    # # Regression #1
    # # bs_bottomup_final[[1]] # "wsi.tmin1" "coy.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #2
    # # bs_bottomup_final[[2]] # "wsi.tmin1"    "forest.tmin1"
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),      # exog on exog
    # #########
    # # Regression #3
    # # bs_bottomup_final[[3]] # "wsi.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #4
    # # bs_bottomup_final[[4]] # "wsi.tmin1"   "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #5
    # # bs_bottomup_final[[5]] # "wsi.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #6
    # # bs_bottomup_final[[6]] # "wsi.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #7
    # bs_bottomup_final[[7]] # "wsi.tmin1" "coy.t"     "coy.tmin1" "wtd.tmin1" "wtd.t"    
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "wsi"), spp = c(".coy", ".wtd", ".wtd", ".wsi"), indices = as.integer(c(1,1,2,1)), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #8
    # # bs_bottomup_final[[8]] # "wsi.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #9
    # bs_bottomup_final[[9]] # "wsi.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"   
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "wsi"), spp = c(".forest", ".bear", ".wtd", ".elk", ".wsi"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y","y-1","y-1")),
    # Regression #10
    # bs_bottomup_final[[10]] # "wsi.tmin1"   "wolf.t"      "moose.tmin1" "moose.t"     "wolf.tmin1"  "wtd.t"       "elk.tmin1"   "elk.t"  
    list(dSep_test = 2, covariates = c("moose.latent", "moose.latent", "wolf.latent", "wtd.latent", "elk.latent", "elk.latent", "wsi"), spp = c(".moose", ".moose", ".wolf", ".wtd", ".elk", ".elk", ".wsi"), indices = as.integer(c(1,2,1,1,1,2,1)), lags = c("y-1","y","y-1","y","y-1","y","y-1")),
    # Regression #11
    # bs_bottomup_final[[11]] # "wsi.tmin1" "lion.t"    "wtd.tmin1" "wtd.t"     "elk.tmin1" "elk.t"
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "elk.latent", "wsi"), spp = c(".wtd", ".wtd", ".elk", ".elk", ".wsi"), indices = as.integer(c(1,2,1,2,1)), lags = c("y-1","y","y-1","y","y-1")),
    # #########
    # # Regression #12          # NOTE: I'VE FLIPPED THE RESPONSE & EXPLANATORY VARIABLES HERE SO THE INDEP. CLAIM IS FOCUSED ON SPP NOT FOREST
    # # bs_bottomup_final[[12]] # "coy.tmin1"    "forest.tmin1"
    # list(dSep_test = 4, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #13
    # # bs_bottomup_final[[13]] # "coy.tmin1"  "bear.tmin1"
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #14
    # # bs_bottomup_final[[14]] # "coy.tmin1"   "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #15
    # bs_bottomup_final[[15]] # "coy.tmin1"    "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1" 
    list(dSep_test = 6, covariates = c("wsi", "forest", "moose.latent", "coy.latent"), spp = c(".wsi", ".forest", ".moose", ".coy"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #16
    # # bs_bottomup_final[[16]] # "coy.tmin1"  "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #17
    # # bs_bottomup_final[[17]] # "coy.tmin1" "wtd.tmin1"
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #18
    # # bs_bottomup_final[[18]] # "coy.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #19
    # bs_bottomup_final[[19]] # "coy.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"   
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "coy.latent"), spp = c(".forest", ".bear", ".wtd", ".elk", ".coy"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y","y-1","y-1")),
    # Regression #20
    # bs_bottomup_final[[20]] # "coy.tmin1"    "elk.t"        "wsi.tmin1"    "forest.tmin1" "bear.tmin1"   "elk.tmin1"   
    list(dSep_test = 5, covariates = c("wsi", "forest", "bear.latent", "elk.latent", "coy.latent"), spp = c(".wsi", ".forest", ".bear", ".elk", ".coy"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #21
    # bs_bottomup_final[[21]] # "coy.tmin1"   "wolf.t"      "moose.tmin1" "moose.t"     "wolf.tmin1"  "wtd.t"       "elk.tmin1"   "elk.t"
    list(dSep_test = 2, covariates = c("moose.latent", "moose.latent", "wolf.latent", "wtd.latent", "elk.latent", "elk.latent", "coy.latent"), spp = c(".moose", ".moose", ".wolf", ".wtd", ".elk", ".elk", ".coy"), indices = as.integer(c(1,2,1,1,1,2,1)), lags = c("y-1","y","y-1","y","y-1","y","y-1")),
    # Regression #22
    # bs_bottomup_final[[22]] # "coy.tmin1" "lion.t"    "wtd.tmin1" "wtd.t"     "elk.tmin1" "elk.t"    
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "elk.latent", "coy.latent"), spp = c(".wtd", ".wtd", ".elk", ".elk", ".coy"), indices = as.integer(c(1,2,1,2,1)), lags = c("y-1","y","y-1","y","y-1")),
    # #########
    # # Regression #23
    # # bs_bottomup_final[[23]] # "forest.tmin1" "bear.tmin1"  
    # list(dSep_test = 3, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #24
    # # bs_bottomup_final[[34]] # "forest.tmin1" "moose.tmin1" 
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #25
    # # bs_bottomup_final[[25]] # "forest.tmin1" "wolf.tmin1"  
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #26
    # # bs_bottomup_final[[26]] # "forest.tmin1" "wtd.tmin1"   
    # list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),      # exog on exog
    # Regression #27
    # bs_bottomup_final[[27]] # "forest.tmin1" "coy.t"        "coy.tmin1"    "wtd.tmin1"    "wtd.t"     
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "forest"), spp = c(".coy", ".wtd", ".wtd", ".forest"), indices = as.integer(c(1,1,2,1)), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #28
    # # bs_bottomup_final[[28]] # "forest.tmin1" "elk.tmin1"   
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #29
    # bs_bottomup_final[[29]] # "forest.tmin1" "wolf.t"       "moose.tmin1"  "moose.t"      "wolf.tmin1"   "wtd.t"        "elk.tmin1"    "elk.t"       
    list(dSep_test = 2, covariates = c("moose.latent", "moose.latent", "wolf.latent", "wtd.latent", "elk.latent", "elk.latent", "forest"), spp = c(".moose", ".moose", ".wolf", ".wtd", ".elk", ".elk", ".forest"), indices = as.integer(c(1,2,1,1,1,2,1)), lags = c("y-1","y","y-1","y","y-1","y","y-1")),
    # Regression #30
    # bs_bottomup_final[[30]] # "forest.tmin1" "lion.t"       "wtd.tmin1"    "wtd.t"        "elk.tmin1"    "elk.t"       
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "elk.latent", "forest"), spp = c(".wtd", ".wtd", ".elk", ".elk", ".forest"), indices = as.integer(c(1,2,1,2,1)), lags = c("y-1","y","y-1","y","y-1")),
    # #########
    # # Regression #31
    # # bs_bottomup_final[[31]] # "bear.tmin1"  "moose.tmin1"
    # list(dSep_test = 6, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #32
    # bs_bottomup_final[[32]] # "bear.tmin1"   "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1" 
    list(dSep_test = 6, covariates = c("wsi", "forest", "moose.latent", "bear.latent"), spp = c(".wsi", ".forest", ".moose", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    # #########
    # # Regression #33
    # # bs_bottomup_final[[33]] # "bear.tmin1" "wolf.tmin1"
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #34
    # # bs_bottomup_final[[34]] # "bear.tmin1" "wtd.tmin1" 
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #35
    # bs_bottomup_final[[35]] # "bear.tmin1" "coy.t"      "coy.tmin1"  "wtd.tmin1"  "wtd.t"     
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "bear.latent"), spp = c(".coy", ".wtd", ".wtd", '.bear'), indices = as.integer(c(1,1,2,1)), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #36
    # # bs_bottomup_final[[36]] # "bear.tmin1" "elk.tmin1" 
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #37
    # bs_bottomup_final[[37]] # "bear.tmin1"  "wolf.t"      "moose.tmin1" "moose.t"     "wolf.tmin1"  "wtd.t"       "elk.tmin1"   "elk.t"   
    list(dSep_test = 2, covariates = c("moose.latent", "moose.latent", "wolf.latent", "wtd.latent", "elk.latent", "elk.latent", "bear.latent"), spp = c(".moose", ".moose", ".wolf", ".wtd", ".elk", ".elk", ".bear"), indices = as.integer(c(1,2,1,1,1,2,1)), lags = c("y-1","y","y-1","y","y-1","y","y-1")),
    # Regression #38
    # bs_bottomup_final[[38]] # "bear.tmin1" "lion.t"     "wtd.tmin1"  "wtd.t"      "elk.tmin1"  "elk.t"     
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "elk.latent", "bear.latent"), spp = c(".wtd", ".wtd", ".elk", ".elk", ".bear"), indices = as.integer(c(1,2,1,2,1)), lags = c("y-1","y","y-1","y","y-1")),
    # #########
    # # Regression #39
    # # bs_bottomup_final[[39]] # "moose.tmin1" "wolf.tmin1" 
    # list(dSep_test = 2, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # #########
    # # Regression #40
    # # bs_bottomup_final[[40]] # "moose.tmin1" "wtd.tmin1"  
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #41
    # bs_bottomup_final[[41]] # "moose.tmin1"  "wtd.t"        "wsi.tmin1"    "coy.tmin1"    "forest.tmin1" "bear.tmin1"   "wtd.tmin1"   
    list(dSep_test = 7, covariates = c("wsi", "coy.latent", "forest", "bear.latent", "wtd.latent", "moose.latent"), spp = c(".wsi", ".coy", ".forest", ".bear", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #42
    # bs_bottomup_final[[42]] # "moose.tmin1" "coy.t"       "coy.tmin1"   "wtd.tmin1"   "wtd.t"    
    list(dSep_test = 4, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "moose.latent"), spp = c(".coy", ".wtd", ".wtd", ".moose"), indices = as.integer(c(1,1,2,1)), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #43
    # # bs_bottomup_final[[43]] # "moose.tmin1" "elk.tmin1"  
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #44
    # bs_bottomup_final[[44]] # "moose.tmin1"  "bear.t"       "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"   
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "moose.latent"), spp = c(".forest", ".bear", ".wtd", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y","y-1","y-1")),
    # Regression #45
    # bs_bottomup_final[[45]] # "moose.tmin1"  "elk.t"        "wsi.tmin1"    "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi", "forest", "bear.latent", "elk.latent", "moose.latent"), spp = c(".wsi", ".forest", ".bear", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #46
    # bs_bottomup_final[[46]] # "moose.tmin1" "lion.t"      "wtd.tmin1"   "wtd.t"       "elk.tmin1"   "elk.t"    
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "elk.latent", "moose.latent"), spp = c(".wtd", ".wtd", ".elk", ".elk", ".moose"), indices = as.integer(c(1,2,1,2,1)), lags = c("y-1","y","y-1","y","y-1")),
    # Regression #47
    # bs_bottomup_final[[47]] # "moose.t"      "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "moose.tmin1" 
    list(dSep_test = 2, covariates = c("wsi", "forest", "moose.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #48
    # bs_bottomup_final[[48]] # "moose.t"      "wtd.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1" 
    list(dSep_test = 7, covariates = c("wsi", "forest", "moose.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #49
    # bs_bottomup_final[[49]] # "moose.t"      "wtd.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "coy.tmin1"    "bear.tmin1"   "wtd.tmin1"   
    list(dSep_test = 7, covariates = c("wsi", "forest", "moose.latent", "coy.latent", "bear.latent", "wtd.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".coy", ".bear", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #50
    # bs_bottomup_final[[50]] # "moose.t"      "coy.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "coy.tmin1"    "wtd.tmin1"    "wtd.t"       
    list(dSep_test = 4, covariates = c("wsi", "forest", "moose.latent", "coy.latent", "wtd.latent", "wtd.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".coy", ".wtd", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y","y")),
    # Regression #51
    # bs_bottomup_final[[51]] # "moose.t"      "elk.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1" 
    list(dSep_test = 5, covariates = c("wsi", "forest", "moose.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".moose"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    # Regression #52
    # bs_bottomup_final[[52]] # "moose.t"      "bear.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "bear.tmin1"   "wtd.t"        "elk.tmin1" 
    list(dSep_test = 3, covariates = c("wsi", "forest", "moose.latent", "bear.latent", "wtd.latent", "elk.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".bear", ".wtd", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y","y-1","y")),
    # Regression #53
    # bs_bottomup_final[[53]] # "moose.t"      "elk.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "bear.tmin1"   "elk.tmin1"   
    list(dSep_test = 5, covariates = c("wsi", "forest", "moose.latent", "bear.latent", "elk.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".bear", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #54
    # bs_bottomup_final[[54]] # "moose.t"      "lion.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wtd.tmin1"    "wtd.t"        "elk.tmin1"    "elk.t"       
    list(dSep_test = 1, covariates = c("wsi", "forest", "moose.latent", "wtd.latent", "wtd.latent", "elk.latent", "elk.latent", "moose.latent"), spp = c(".wsi", ".forest", ".moose", ".wtd", ".wtd", ".elk", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,2,1,2,2)), lags = c("y-1","y-1","y-1","y-1","y","y-1","y","y")),
    # #########
    # # Regression #55
    # # bs_bottomup_final[[55]] # "wolf.tmin1" "wtd.tmin1" 
    # list(dSep_test = 7, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #56
    # bs_bottomup_final[[56]] # "wolf.tmin1"   "wtd.t"        "wsi.tmin1"    "coy.tmin1"    "forest.tmin1" "bear.tmin1"   "wtd.tmin1"   
    list(dSep_test = 7, covariates = c("wsi", "coy.latent", "forest", "bear.latent", "wtd.latent", "wolf.latent"), spp = c(".wsi", ".coy", ".forest", ".bear", ".wtd", ".wolf"), indices = as.integer(c(1,1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1")),
    # Regression #57
    # bs_bottomup_final[[57]] # "wolf.tmin1" "coy.t"      "coy.tmin1"  "wtd.tmin1"  "wtd.t"     
    list(dSep_test = 4, covariates = c("coy.latent", 'wtd.latent', "wtd.latent", "wolf.latent"), spp = c(".coy", ".wtd", ".wtd", ".wolf"), indices = as.integer(c(1,1,2,1)), lags = c("y-1","y-1","y","y-1")),
    # #########
    # # Regression #58
    # # bs_bottomup_final[[58]] # "wolf.tmin1" "elk.tmin1" 
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #59
    # bs_bottomup_final[[59]] # "wolf.tmin1"   "bear.t"       "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"   
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "wolf.latent"), spp = c(".forest", ".bear", ".wtd", ".elk", ".wolf"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y","y-1","y-1")),
    # Regression #60
    # bs_bottomup_final[[60]] # "wolf.tmin1"   "elk.t"        "wsi.tmin1"    "forest.tmin1" "bear.tmin1"   "elk.tmin1"   
    list(dSep_test = 5, covariates = c("wsi", "forest", "bear.latent", "elk.latent", "wolf.latent"), spp = c(".wsi", ".forest", ".bear", ".elk", ".wolf"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #61
    # bs_bottomup_final[[61]] # "wolf.tmin1" "lion.t"     "wtd.tmin1"  "wtd.t"      "elk.tmin1"  "elk.t"   
    list(dSep_test = 1, covariates = c("wtd.latent", "wtd.latent", "elk.latent", "elk.latent", "wolf.latent"), spp = c(".wtd", ".wtd", ".elk", ".elk", ".wolf"), indices = as.integer(c(1,2,1,2,1)), lags = c("y-1","y","y-1","y","y-1")),
    # #########
    # # Regression #62
    # # bs_bottomup_final[[62]] # "wtd.tmin1" "elk.tmin1"
    # list(dSep_test = 5, covariates = c(), spp = c(), indices = as.integer(c()), lags = c()),
    # Regression #63
    # bs_bottomup_final[[63]] # "wtd.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"   
    list(dSep_test = 3, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "wtd.latent"), spp = c(".forest", ".bear", ".wtd", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1,2)), lags = c("y-1","y-1","y","y-1","y-1")),
    # Regression #64
    # bs_bottomup_final[[64]] # "wtd.tmin1"    "elk.t"        "wsi.tmin1"    "forest.tmin1" "bear.tmin1"   "elk.tmin1"   
    list(dSep_test = 5, covariates = c("wsi", "forest", "bear.latent", "elk.latent", "wtd.latent"), spp = c(".wsi", ".forest", ".bear", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    # Regression #65
    # bs_bottomup_final[[65]] # "wtd.tmin1"   "wolf.t"      "moose.tmin1" "moose.t"     "wolf.tmin1"  "wtd.t"       "elk.tmin1"   "elk.t"      
    list(dSep_test = 2, covariates = c("moose.latent", "moose.latent", "wolf.latent", "wtd.latent", "elk.latent", "elk.latent", "wtd.latent"), spp = c(".moose", ".moose", ".wolf", ".wtd", ".elk", ".elk", ".wtd"), indices = as.integer(c(1,2,1,1,1,2,2)), lags = c("y-1","y","y-1","y","y-1","y","y-1")),
    # Regression #66
    # bs_bottomup_final[[66]] # "wtd.t"        "elk.tmin1"    "wsi.tmin1"    "coy.tmin1"    "forest.tmin1" "bear.tmin1"   "wtd.tmin1"   
    list(dSep_test = 5, covariates = c("wsi", "coy.latent", "forest", "bear.latent", "wtd.latent", "wtd.latent"), spp = c(".wsi", ".coy", ".forest", ".bear", ".wtd", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #67
    # bs_bottomup_final[[67]] # "wtd.t"        "elk.t"        "wsi.tmin1"    "coy.tmin1"    "forest.tmin1" "bear.tmin1"   "wtd.tmin1"    "elk.tmin1"   
    list(dSep_test = 5, covariates = c("wsi", "coy.latent", "forest", "bear.latent", "wtd.latent", "elk.latent", "wtd.latent"), spp = c(".wsi", ".coy", ".forest", ".bear", ".wtd", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y-1","y")),
    # Regression #68
    # bs_bottomup_final[[68]] # "coy.t"     "elk.tmin1" "coy.tmin1" "wtd.tmin1" "wtd.t"    
    list(dSep_test = 5, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "coy.latent"), spp = c(".coy", ".wtd", ".wtd", ".coy"), indices = as.integer(c(1,1,2,2)), lags = c("y-1","y-1","y","y")),
    # Regression #69
    # bs_bottomup_final[[69]] # "coy.t"        "bear.t"       "coy.tmin1"    "wtd.tmin1"    "wtd.t"        "forest.tmin1" "bear.tmin1"   "elk.tmin1"   
    list(dSep_test = 3, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "forest", "bear.latent", "elk.latent", "coy.latent"), spp = c(".coy", ".wtd", ".wtd", ".forest", ".bear", ".elk", ".coy"), indices = as.integer(c(1,1,2,1,1,1,2)), lags = c("y-1","y-1","y","y-1","y-1","y-1","y")),
    # Regression #70
    # bs_bottomup_final[[70]] # "coy.t"        "elk.t"        "coy.tmin1"    "wtd.tmin1"    "wtd.t"        "wsi.tmin1"    "forest.tmin1" "bear.tmin1"   "elk.tmin1"   
    list(dSep_test = 5, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "wsi", "forest", "bear.latent", "elk.latent", "coy.latent"), spp = c(".coy", ".wtd", ".wtd", ".wsi", ".forest", ".bear", ".elk", ".coy"), indices = as.integer(c(1,1,2,1,1,1,1,2)), lags = c("y-1","y-1","y","y-1","y-1","y-1","y-1","y")),
    # Regression #71
    # bs_bottomup_final[[71]] # "coy.t"       "wolf.t"      "coy.tmin1"   "wtd.tmin1"   "wtd.t"       "moose.tmin1" "moose.t"     "wolf.tmin1"  "elk.tmin1"   "elk.t"      
    list(dSep_test = 2, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "moose.latent", "moose.latent", "wolf.latent", "elk.latent", "elk.latent", "coy.latent"), spp = c(".coy", ".wtd", ".wtd", ".moose", ".moose", ".wolf", ".elk", ".elk", ".coy"), indices = as.integer(c(1,1,2,1,2,1,1,2,2)), lags = c("y-1","y-1","y","y-1","y","y-1","y-1","y","y")),
    # Regression #72
    # bs_bottomup_final[[72]] # "coy.t"     "lion.t"    "coy.tmin1" "wtd.tmin1" "wtd.t"     "elk.tmin1" "elk.t"    
    list(dSep_test = 1, covariates = c("coy.latent", "wtd.latent", "wtd.latent", "elk.latent", "elk.latent" ,"coy.latent"), spp = c(".coy", ".wtd", ".wtd", ".elk", ".elk", ".coy"), indices = as.integer(c(1,1,2,1,2,2)), lags = c("y-1","y-1","y","y-1","y","y")),
    # Regression #73
    # bs_bottomup_final[[73]] # "bear.t"       "elk.t"        "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"    "wsi.tmin1"   
    list(dSep_test = 5, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "wsi", "bear.latent"), spp = c(".forest", ".bear", ".wtd", ".elk", ".wsi", ".bear"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y","y-1","y-1","y")),
    # Regression #74
    # bs_bottomup_final[[74]] # "bear.t"       "wolf.t"       "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"    "moose.tmin1"  "moose.t"      "wolf.tmin1"   "elk.t"    
    list(dSep_test = 2, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "moose.latent", "moose.latent", "wolf.latent", "elk.latent", "bear.latent"), spp = c(".forest", ".bear", ".wtd", ".elk", ".moose", ".moose", ".wolf", ".elk", ".bear"), indices = as.integer(c(1,1,1,1,1,2,1,2,2)), lags = c("y-1","y-1","y","y-1","y-1","y","y-1","y","y")),
    # Regression #75
    # bs_bottomup_final[[75]] # "bear.t"       "lion.t"       "forest.tmin1" "bear.tmin1"   "wtd.t"        "elk.tmin1"    "wtd.tmin1"    "elk.t"      
    list(dSep_test = 1, covariates = c("forest", "bear.latent", "wtd.latent", "elk.latent", "wtd.latent", "elk.latent", "bear.latent"), spp = c(".forest", ".bear", ".wtd", ".elk", ".wtd", ".elk", ".bear"), indices = as.integer(c(1,1,1,1,2,2,2)), lags = c("y-1","y-1","y","y-1","y-1","y","y")),
    # Regression #76
    # bs_bottomup_final[[76]] # "wolf.t"      "lion.t"      "moose.tmin1" "moose.t"     "wolf.tmin1"  "wtd.t"       "elk.tmin1"   "elk.t"       "wtd.tmin1"
    list(dSep_test = 1, covariates = c("moose.latent", "moose.latent", "wolf.latent", "wtd.latent", "elk.latent", "elk.latent", "wtd.latent", "wolf.latent"), spp = c(".moose", ".moose", ".wolf", ".wtd", ".elk", ".elk", ".wtd", ".wolf"), indices = as.integer(c(1,2,1,1,1,2,2,2)), lags = c("y-1","y","y-1","y","y-1","y","y-1","y"))
  )
  
  
  
  