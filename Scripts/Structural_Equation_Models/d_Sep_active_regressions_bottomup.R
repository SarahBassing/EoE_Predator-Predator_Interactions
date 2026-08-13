  #'  ----------------------------------------------------------------
  #'  Active regression for d-Sep iterations: bottom-up
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
  #'     easier when using the p.rope function in d_Sep_test_ROPE_method.R
  #'  -- spp: vector containing the suffex of each parameter name to be appended
  #'     to "beta" in the model template. Must include a period before the term's
  #'     name to work.
  #'  -- indices: vector containing the index value for each parameter. For most
  #'     regressions, each parameter's name is unique so the index is 1. But for some,
  #'     the parameter name is used twice in the model representing different coefficients
  #'     (e.g., beta.harvest[1] and beta.harvest[2] correspond to the effects of 
  #'     wolfHarvets.tmin1 and lionHarvest.tmin1 in the same regression). Double check
  #'     indexing for each regression to ensure appropriate coefficient estimates are
  #'     saved and used in ROPE method.
  #'  ---------------------------------------------------------------- 
  
  dSep_iterations_bottomup <- list(
    # Regression #1
    # bs_bottomup[[1]] # "wsi.tmin1" "coy.t"     "coy.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wtd.tmin1", "wsi.tmin1"), spp = c(".coy", ".wtd", ".wsi"), indices = as.integer(c(1,1,1))),
    # Regression #2
    # bs_bottomup[[2]] # "wsi.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "elk.tmin1", "wsi.tmin1"), spp = c(".forest", ".bear", ".elk", ".wsi"), indices = as.integer(c(1,1,1,1))),
    # bs_bottomup[[3]] # "wsi.tmin1"   "wolf.t"      "moose.tmin1" "wolf.tmin1"  "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1", "wsi.tmin1"), spp = c(".moose", ".wolf", ".elk", ".wsi"), indices = as.integer(c(1,1,1,1))),
    # Regression #4
    # bs_bottomup[[4]] # "wsi.tmin1"  "lion.t"     "wtd.tmin1"  "elk.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "elk.tmin1", "lion.tmin1", "wsi.tmin1"), spp = c(".wtd", ".elk", ".lion", ".wsi"), indices = as.integer(c(1,1,1,1))),
    # Regression #5
    # bs_bottomup[[5]] # "coy.tmin1"    "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "coy.tmin1"), spp = c(".wsi", ".forest", ".moose", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #6
    # bs_bottomup[[6]] # "coy.tmin1"    "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "coy.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #7
    # bs_bottomup[[7]] # "coy.tmin1"    "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "coy.tmin1"), spp = c(".wsi", ".forest", ".elk", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #8
    # bs_bottomup[[8]] # "coy.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "elk.tmin1", "coy.tmin1"), spp = c(".forest", ".bear", ".elk", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #9
    # bs_bottomup[[9]] # "coy.tmin1"   "wolf.t"      "moose.tmin1" "wolf.tmin1"  "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1", "coy.tmin1"), spp = c(".moose", ".wolf", ".elk", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #10
    # bs_bottomup[[10]] # "coy.tmin1"  "lion.t"     "wtd.tmin1"  "elk.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "elk.tmin1", "lion.tmin1", "coy.tmin1"), spp = c(".wtd", ".elk", ".lion", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #11
    # bs_bottomup[[11]] # "forest.tmin1" "coy.t"        "coy.tmin1"    "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wtd.tmin1", "forest.tmin1"), spp = c(".coy", ".wtd", ".forest"), indices = as.integer(c(1,1,1))),
    # Regression #12
    # bs_bottomup[[12]] # "forest.tmin1" "wolf.t"       "moose.tmin1"  "wolf.tmin1"   "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1", "forest.tmin1"), spp = c(".moose", ".wolf", ".elk", ".forest"), indices = as.integer(c(1,1,1,1))),
    # Regression #13
    # bs_bottomup[[13]] # "forest.tmin1" "lion.t"       "wtd.tmin1"    "elk.tmin1"    "lion.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "elk.tmin1", "lion.tmin1", "forest.tmin1"), spp = c(".wtd", ".elk", ".lion", ".forest"), indices = as.integer(c(1,1,1,1))),
    # Regression #14
    # bs_bottomup[[14]] # "bear.tmin1"   "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "bear.tmin1"), spp = c(".wsi", ".forest", ".moose", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #15
    # bs_bottomup[[15]] # "bear.tmin1"   "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "bear.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #16
    # bs_bottomup[[16]] # "bear.tmin1" "coy.t"      "coy.tmin1"  "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wtd.tmin1", "bear.tmin1"), spp = c(".coy", ".wtd", ".bear"), indices = as.integer(c(1,1,1))),
    # Regression #17
    # bs_bottomup[[17]] # "bear.tmin1"   "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "bear.tmin1"), spp = c(".wsi", ".forest", ".elk", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #18
    # bs_bottomup[[18]] # "bear.tmin1"  "wolf.t"      "moose.tmin1" "wolf.tmin1"  "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1", "bear.tmin1"), spp = c(".moose", ".wolf", ".elk", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #19
    # bs_bottomup[[19]] # "bear.tmin1" "lion.t"     "wtd.tmin1"  "elk.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "elk.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".wtd", ".elk", ".lion", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #20
    # bs_bottomup[[20]] # "moose.tmin1"  "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "moose.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1))),
    # Regression #21
    # bs_bottomup[[21]] # "moose.tmin1" "coy.t"       "coy.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wtd.tmin1", "moose.tmin1"), spp = c(".coy", ".wtd", ".moose"), indices = as.integer(c(1,1,1))),
    # Regression #22
    # bs_bottomup[[22]] # "moose.tmin1"  "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "moose.tmin1"), spp = c(".wsi", ".forest", ".elk", ".moose"), indices = as.integer(c(1,1,1,1))),
    # Regression #23
    # bs_bottomup[[23]] # "moose.tmin1"  "bear.t"       "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "elk.tmin1", "moose.tmin1"), spp = c(".forest", ".bear", ".elk", ".moose"), indices = as.integer(c(1,1,1,1))),
    # Regression #24
    # bs_bottomup[[24]] # "moose.tmin1" "lion.t"      "wtd.tmin1"   "elk.tmin1"   "lion.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "elk.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".wtd", ".elk", ".lion", ".moose"), indices = as.integer(c(1,1,1,1))),
    #### bs_bottomup[[25]] # NOT POSSIBLE- t cannot affect t-1 # "moose.t"      "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    #### bs_bottomup[[26]] # NOT POSSIBLE- t cannot affect t-1 # "moose.t"      "wtd.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #25
    # bs_bottomup[[27]] # "moose.t"      "wtd.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "wtd.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # Regression #26
    # bs_bottomup[[28]] # "moose.t"      "coy.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "coy.tmin1"    "wtd.tmin1"
    list(dSep_test = 4, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "coy.tmin1", "wtd.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".coy", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup[[29]] # NOT POSSIBLE- t cannot affect t-1 # "moose.t"      "elk.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #27
    # bs_bottomup[[30]] # "moose.t"      "elk.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "elk.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # Regression #28
    # bs_bottomup[[31]] # "moose.t"      "bear.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 3, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "bear.tmin1", "elk.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".bear", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #29
    # bs_bottomup[[32]] # "moose.t"      "wolf.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wolf.tmin1"   "elk.tmin1"
    list(dSep_test = 2, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "wolf.tmin1", "elk.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".wolf", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup[[33]] # NOT POSSIBLE- t cannot affect t-1 # "moose.t"      "lion.tmin1"   "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #30
    # bs_bottomup[[34]] # "moose.t"      "lion.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wtd.tmin1"    "elk.tmin1"    "lion.tmin1"
    list(dSep_test = 1, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "wtd.tmin1", "elk.tmin1", "lion.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".wtd", ".elk", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,1,1,2))),
    # Regression #31
    # bs_bottomup[[35]] # "wolf.tmin1"   "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "wolf.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".wolf"), indices = as.integer(c(1,1,1,1))),
    # Regression #32
    # bs_bottomup[[36]] # "wolf.tmin1" "coy.t"      "coy.tmin1"  "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wtd.tmin1", "wolf.tmin1"), spp = c(".coy", ".wtd", ".wolf"), indices = as.integer(c(1,1,1))),
    # Regression # 33
    # bs_bottomup[[37]] # "wolf.tmin1"   "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "wolf.tmin1"), spp = c(".wsi", ".forest", ".elk", ".wolf"), indices = as.integer(c(1,1,1,1))),
    # Regression #34
    # bs_bottomup[[38]] # "wolf.tmin1"   "bear.t"       "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "elk.tmin1", "wolf.tmin1"), spp = c(".forest", ".bear", ".elk", ".wolf"), indices = as.integer(c(1,1,1,1))),
    # Regression #35
    # bs_bottomup[[39]] # "wolf.tmin1" "lion.t"     "wtd.tmin1"  "elk.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "elk.tmin1", "lion.tmin1", "wolf.tmin1"), spp = c(".wtd", ".elk", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1))),
    # Regression #36
    # bs_bottomup[[40]] # "wtd.tmin1"    "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "wtd.tmin1"), spp = c(".wsi", ".forest", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1))),
    # Regression #37
    # bs_bottomup[[41]] # "wtd.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "elk.tmin1", "wtd.tmin1"), spp = c(".forest", ".bear", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1))),
    # Regression #38
    # bs_bottomup[[42]] # "wtd.tmin1"   "wolf.t"      "moose.tmin1" "wolf.tmin1"  "elk.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1", "wtd.tmin1"), spp = c(".moose", ".wolf", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1))),
    # Regression #39
    # bs_bottomup[[43]] # "wtd.t"        "coy.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "coy.tmin1"
    list(dSep_test = 4, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "coy.tmin1", "wtd.t"), spp = c(".wsi", ".forest", ".wtd", ".coy", ".wtd"), indices = as.integer(c(1,1,1,1,2))),
    #### bs_bottomup[[44]] # NOT POSSIBLE- t cannot affect t-1 # "wtd.t"        "elk.tmin1"    "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #40
    # bs_bottomup[[45]] # "wtd.t"        "elk.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "elk.tmin1", "wtd.t"), spp = c(".wsi", ".forest", ".wtd", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1,2))),
    # Regression #41
    # bs_bottomup[[46]] # "wtd.t"        "bear.t"       "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "bear.tmin1"   "elk.tmin1" 
    list(dSep_test = 3, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "bear.tmin1", "elk.tmin1", "wtd.t"), spp = c(".wsi", ".forest", ".wtd", ".bear", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #42 
    # bs_bottomup[[47]] # "wtd.t"        "wolf.t"       "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "moose.tmin1"  "wolf.tmin1"   "elk.tmin1"
    list(dSep_test = 2, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "moose.tmin1", "wolf.tmin1", "elk.tmin1", "wtd.t"), spp = c(".wsi", ".forest", ".wtd", ".moose", ".wolf", ".elk", ".wtd"), indices = as.integer(c(1,1,1,1,1,1,2))),
    #### bs_bottomup[[48]] # NOT POSSIBLE- t cannot affect t-1 # "wtd.t"        "lion.tmin1"   "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    ####list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #43 
    # bs_bottomup[[49]] # "wtd.t"        "lion.t"       "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "elk.tmin1"    "lion.tmin1"
    list(dSep_test = 2, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "elk.tmin1", "lion.tmin1", "wtd.t"), spp = c(".wsi", ".forest", ".wtd", ".elk", ".lion", ".wtd"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup[[50]] # NOT POSSIBLE- t cannot affect t-1  # "coy.t"     "elk.tmin1" "coy.tmin1" "wtd.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #44 
    # bs_bottomup[[51]] # "coy.t"        "elk.t"        "coy.tmin1"    "wtd.tmin1"    "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("coy.tmin1", "wtd.tmin1", "wsi.tmin1", "forest.tmin1", "elk.tmin1", "coy.t"), spp = c(".coy", ".wtd", ".wsi", ".forest", ".elk", ".coy"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #45 
    # bs_bottomup[[52]] # "coy.t"        "bear.t"       "coy.tmin1"    "wtd.tmin1"    "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    list(dSep_test = 3, covariates = c("coy.tmin1", "wtd.tmin1", "forest.tmin1", "bear.tmin1", "elk.tmin1", "coy.t"), spp = c(".coy", ".wtd", ".forest", ".bear", ".elk", ".coy"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #46 
    # bs_bottomup[[53]] # "coy.t"       "wolf.t"      "coy.tmin1"   "wtd.tmin1"   "moose.tmin1" "wolf.tmin1"  "elk.tmin1"
    list(dSep_test = 2, covariates = c("coy.tmin1", "wtd.tmin1", "moose.tmin1", "wolf.tmin1", "elk.tmin1", "coy.t"), spp = c(".coy", ".wtd", ".moose", ".wolf", ".elk", ".coy"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup[[54]] # NOT POSSIBLE- t cannot affect t-1  # "coy.t"      "lion.tmin1" "coy.tmin1"  "wtd.tmin1" 
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #47 
    # bs_bottomup[[55]] # "coy.t"      "lion.t"     "coy.tmin1"  "wtd.tmin1"  "elk.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("coy.tmin1", "wtd.tmin1", "elk.tmin1", "lion.tmin1", "coy.t"), spp = c(".coy", ".wtd", ".elk", ".lion", ".coy"), indices = as.integer(c(1,1,1,1,2))),
    # Regression #48 
    # bs_bottomup[[56]] # "elk.t"        "bear.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "bear.tmin1" 
    list(dSep_test = 3, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "bear.tmin1", "elk.t"), spp = c(".wsi", ".forest", ".elk", ".bear", ".elk"), indices = as.integer(c(1,1,1,1,2))),
    # Regression #49 
    # bs_bottomup[[57]] # "elk.t"        "wolf.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "moose.tmin1"  "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "moose.tmin1", "wolf.tmin1", "elk.t"), spp = c('.wsi', ".forest", ".elk", ".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup[[58]] # NOT POSSIBLE- t cannot affect t-1  # "elk.t"        "lion.tmin1"   "wsi.tmin1"    "forest.tmin1" "elk.tmin1" 
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #50 
    # bs_bottomup[[59]] # "elk.t"        "lion.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "wtd.tmin1"    "lion.tmin1"
    list(dSep_test = 1, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "wtd.tmin1", "lion.tmin1", "elk.t"), spp = c(".wsi", ".forest", ".elk", ".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #51 
    # bs_bottomup[[60]] # "bear.t"       "wolf.t"       "forest.tmin1" "bear.tmin1"   "elk.tmin1"    "moose.tmin1"  "wolf.tmin1" 
    list(dSep_test = 2, covariates = c("forest.tmin1", "bear.tmin1", "elk.tmin1", "moose.tmin1", "wolf.tmin1", "bear.t"), spp = c(".forest", ".bear", ".elk", ".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup[[61]] # NOT POSSIBLE- t cannot affect t-1  # "bear.t"       "lion.tmin1"   "forest.tmin1" "bear.tmin1"   "elk.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #52 
    # bs_bottomup[[62]] # "bear.t"       "lion.t"       "forest.tmin1" "bear.tmin1"   "elk.tmin1"    "wtd.tmin1"    "lion.tmin1"
    list(dSep_test = 1, covariates = c("forest.tmin1", "bear.tmin1", "elk.tmin1", "wtd.tmin1", "lion.tmin1", "bear.t"), spp = c(".forest", ".bear", ".elk", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup[[63]] # NOT POSSIBLE- t cannot affect t-1  # "wolf.t"      "lion.tmin1"  "moose.tmin1" "wolf.tmin1"  "elk.tmin1" 
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c())),
    # Regression #53 
    # bs_bottomup[[64]] # "wolf.t"      "lion.t"      "moose.tmin1" "wolf.tmin1"  "elk.tmin1"   "wtd.tmin1"   "lion.tmin1"
    list(dSep_test = 1, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1", "wtd.tmin1", "lion.tmin1", "wolf.t"), spp = c(".moose", ".wolf", ".elk", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1,1,2)))
  )
  
  