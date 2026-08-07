  #'  ----------------------------------------------------------------
  #'  Active regression for d-Sep iterations: bottom-up interference
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
  
  dSep_iterations_bottomup_int <- list(
    # Regression #1
    # bs_bottomup_inter[[1]] # "wsi.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "wolf.tmin1", "wsi.tmin1"), spp = c(".forest", ".bear", ".wolf", ".wsi"), indices = as.integer(c(1,1,1,1))),
    # Regression #2
    # bs_bottomup_inter[[2]] # "wsi.tmin1"   "wolf.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "elk.tmin1", "wolf.tmin1", "wsi.tmin1"), spp = c(".moose", ".elk", ".wolf", ".wsi"), indices = as.integer(c(1,1,1,1))),
    # Regression #3
    # bs_bottomup_inter[[3]] # "wsi.tmin1"  "coy.t"      "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "wtd.tmin1", "wsi.tmin1"), spp = c(".coy", ".wolf", ".wtd", ".wsi"), indices = as.integer(c(1,1,1,1))),
    # Regression #4
    # bs_bottomup_inter[[4]] # "wsi.tmin1"  "lion.t"     "wolf.tmin1" "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wolf.tmin1", "wtd.tmin1", "lion.tmin1", "wsi.tmin1"), spp = c(".wolf", ".wtd", ".lion", ".wsi"), indices = as.integer(c(1,1,1,1))),
    # Regression #5
    # bs_bottomup_inter[[5]] # "coy.tmin1"    "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "coy.tmin1"), spp = c(".wsi", ".forest", ".moose", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #6
    # bs_bottomup_inter[[6]] # "coy.tmin1"    "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "coy.tmin1"), spp = c(".wsi", ".forest", ".elk", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #7
    # bs_bottomup_inter[[7]] # "coy.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "wolf.tmin1", "coy.tmin1"), spp = c(".forest", ".bear", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #8
    # bs_bottomup_inter[[8]] # "coy.tmin1"   "wolf.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "elk.tmin1", "wolf.tmin1", "coy.tmin1"), spp = c(".moose", ".elk", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #9
    # bs_bottomup_inter[[9]] # "coy.tmin1"    "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "coy.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #10
    # bs_bottomup_inter[[10]] # "coy.tmin1"  "lion.t"     "wolf.tmin1" "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wolf.tmin1", "wtd.tmin1", "lion.tmin1", "coy.tmin1"), spp = c(".wolf", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1))),
    # Regression #11
    # bs_bottomup_inter[[11]] # "forest.tmin1" "wolf.t"       "moose.tmin1"  "elk.tmin1"    "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "elk.tmin1", "wolf.tmin1", "forest.tmin1"), spp = c(".moose", ".elk", ".wolf", ".forest"), indices = as.integer(c(1,1,1,1))),
    # Regression #12
    # bs_bottomup_inter[[12]] # "forest.tmin1" "coy.t"        "coy.tmin1"    "wolf.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "wtd.tmin1", "forest.tmin1"), spp = c(".coy", ".wolf", ".wtd", ".forest"), indices = as.integer(c(1,1,1,1))),
    # Regression #13
    # bs_bottomup_inter[[13]] # "forest.tmin1" "lion.t"       "wolf.tmin1"   "wtd.tmin1"    "lion.tmin1"
    list(dSep_test = 1, covariates = c("wolf.tmin1", "wtd.tmin1", "lion.tmin1", "forest.tmin1"), spp = c(".wolf", ".wtd", ".lion", ".forest"), indices = as.integer(c(1,1,1,1))),
    # Regression #14
    # bs_bottomup_inter[[14]] # "bear.tmin1"   "moose.t"      "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    list(dSep_test = 6, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "bear.tmin1"), spp = c(".wsi", ".forest", ".moose", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #15
    # bs_bottomup_inter[[15]] # "bear.tmin1"   "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "bear.tmin1"), spp = c(".wsi", ".forest", ".elk", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #16
    # bs_bottomup_inter[[16]] # "bear.tmin1"  "wolf.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "elk.tmin1", "wolf.tmin1", "bear.tmin1"), spp = c(".moose", ".elk", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #17
    # bs_bottomup_inter[[17]] # "bear.tmin1"   "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "bear.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #18
    # bs_bottomup_inter[[18]] # "bear.tmin1" "coy.t"      "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "wtd.tmin1", "bear.tmin1"), spp = c(".coy", ".wolf", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #19
    # bs_bottomup_inter[[19]] # "bear.tmin1" "lion.t"     "wolf.tmin1" "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wolf.tmin1", "wtd.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".wolf", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1))),
    # Regression #20
    # bs_bottomup_inter[[20]] # "moose.tmin1"  "elk.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "moose.tmin1"), spp = c(".wsi", ".forest", ".elk", ".moose"), indices = as.integer(c(1,1,1,1))),
    # Regression #21
    # bs_bottomup_inter[[21]] # "moose.tmin1"  "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "wolf.tmin1", "moose.tmin1"), spp = c(".forest", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1))),
    # Regression #22
    # bs_bottomup_inter[[22]] # "moose.tmin1"  "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "moose.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1))),
    # Regression #23
    # bs_bottomup_inter[[23]] # "moose.tmin1" "coy.t"       "coy.tmin1"   "wolf.tmin1"  "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1",  "wtd.tmin1", "moose.tmin1"), spp = c(".coy", ".wolf", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1))),
    # Regression #24
    # bs_bottomup_inter[[24]] # "moose.tmin1" "lion.t"      "wolf.tmin1"  "wtd.tmin1"   "lion.tmin1"
    list(dSep_test = 1, covariates = c("wolf.tmin1", "wtd.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".wolf", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1))),
    #### bs_bottomup_inter[[25]] # NOT POSSIBLE - t cannot affect t-1 # "moose.t"      "elk.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #25
    # bs_bottomup_inter[[26]] # "moose.t"      "elk.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "elk.tmin1"
    list(dSep_test = 5, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "elk.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".elk", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    #### bs_bottomup_inter[[27]] # NOT POSSIBLE - t cannot affect t-1  # "moose.t"      "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #26
    # bs_bottomup_inter[[28]] # "moose.t"      "bear.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "bear.tmin1", "wolf.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #27
    # bs_bottomup_inter[[29]] # "moose.t"      "wolf.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "elk.tmin1"    "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "elk.tmin1", "wolf.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".elk", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup_inter[[30]] # NOT POSSIBLE - t cannot affect t-1  # "moose.t"      "wtd.tmin1"    "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #28
    # bs_bottomup_inter[[31]] # "moose.t"      "wtd.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "wtd.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # Regression #29
    # bs_bottomup_inter[[32]] # "moose.t"      "coy.t"        "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "coy.tmin1"    "wolf.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "coy.tmin1", "wolf.tmin1", "wtd.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".coy", ".wolf", ".wtd", ".moose"), indices = as.integer(c(1,1,1,1,1,1,2))),
    #### bs_bottomup_inter[[33]] # NOT POSSIBLE - t cannot affect t-1   # "moose.t"      "lion.tmin1"   "wsi.tmin1"    "forest.tmin1" "moose.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #30
    # bs_bottomup_inter[[34]] # "moose.t"      "lion.t"       "wsi.tmin1"    "forest.tmin1" "moose.tmin1"  "wolf.tmin1"   "wtd.tmin1"    "lion.tmin1"
    list(dSep_test = 1, covariates = c("wsi.tmin1", "forest.tmin1", "moose.tmin1", "wolf.tmin1", "wtd.tmin1", "lion.tmin1", "moose.t"), spp = c(".wsi", ".forest", ".moose", ".wolf", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,1,1,2))),
    # Regression #31
    # bs_bottomup_inter[[35]] # "elk.tmin1"    "bear.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("forest.tmin1", "bear.tmin1", "wolf.tmin1", "elk.tmin1"), spp = c(".forest", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1))),
    # Regression #32
    # bs_bottomup_inter[[36]] # "elk.tmin1"    "wtd.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "elk.tmin1"), spp = c(".wsi", ".forest", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1))),
    # Regression #33
    # bs_bottomup_inter[[37]] # "elk.tmin1"  "coy.t"      "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "wtd.tmin1", "elk.tmin1"), spp = c(".coy", ".wolf", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1))),
    # Regression #34
    # bs_bottomup_inter[[38]] # "elk.tmin1"  "lion.t"     "wolf.tmin1" "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wolf.tmin1", "wtd.tmin1", "lion.tmin1", "elk.tmin1"), spp = c(".wolf", ".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1,1))),
    #### bs_bottomup_inter[[39]] # NOT POSSIBLE - t cannot affect t-1  # "elk.t"        "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #35
    # bs_bottomup_inter[[40]] # "elk.t"        "bear.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "bear.tmin1"   "wolf.tmin1"
    list(dSep_test = 3, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "bear.tmin1", "wolf.tmin1", "elk.t"), spp = c(".wsi", ".forest", ".elk", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #36 
    # bs_bottomup_inter[[41]] # "elk.t"        "wolf.t"       "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "moose.tmin1"  "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "moose.tmin1", "wolf.tmin1", "elk.t"), spp = c(".wsi", ".forest", ".elk", ".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup_inter[[42]] # NOT POSSIBLE - t cannot affect t-1  # "elk.t"        "wtd.tmin1"    "wsi.tmin1"    "forest.tmin1" "elk.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #37 
    # bs_bottomup_inter[[43]] # "elk.t"        "wtd.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "wtd.tmin1"
    list(dSep_test = 7, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "wtd.tmin1", "elk.t"), spp = c(".wsi", ".forest", ".elk", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1,2))),
    # Regression #38 
    # bs_bottomup_inter[[44]] # "elk.t"        "coy.t"        "wsi.tmin1"    "forest.tmin1" "elk.tmin1"    "coy.tmin1"    "wolf.tmin1"   "wtd.tmin1"  
    list(dSep_test = 4, covariates = c("wsi.tmin1", "forest.tmin1", "elk.tmin1", "coy.tmin1", "wolf.tmin1", "wtd.tmin1", "elk.t"), spp = c(".wsi", ".forest", ".elk", ".coy", ".wolf", ".wtd", ".elk"), indices = as.integer(c(1,1,1,1,1,1,2))),
    #### bs_bottomup_inter[[45]] # NOT POSSIBLE - t cannot affect t-1 # "elk.t"        "lion.tmin1"   "wsi.tmin1"    "forest.tmin1" "elk.tmin1" 
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #39 
    # bs_bottomup_inter[[50]] # "bear.t"       "wtd.t"        "forest.tmin1" "bear.tmin1"   "wolf.tmin1"   "wsi.tmin1"    "wtd.tmin1"
    list(dSep_test = 7, covariates = c("forest.tmin1", "bear.tmin1", "wolf.tmin1", "wsi.tmin1", "wtd.tmin1", "bear.t"), spp = c(".forest", ".bear", ".wolf", ".wsi", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #40 
    # bs_bottomup_inter[[51]] # "bear.t"       "coy.t"        "forest.tmin1" "bear.tmin1"   "wolf.tmin1"   "coy.tmin1"    "wtd.tmin1"
    list(dSep_test = 4, covariates = c("forest.tmin1", "bear.tmin1", "wolf.tmin1", "coy.tmin1", "wtd.tmin1", "bear.t"), spp = c(".forest", ".bear", ".wolf", ".coy", ".wtd", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup_inter[[52]] # NOT POSSIBLE - t cannot affect t-1 # "bear.t"       "lion.tmin1"   "forest.tmin1" "bear.tmin1"   "wolf.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #41 
    # bs_bottomup_inter[[53]] # "bear.t"       "lion.t"       "forest.tmin1" "bear.tmin1"   "wolf.tmin1"   "wtd.tmin1"    "lion.tmin1"
    list(dSep_test = 1, covariates = c("forest.tmin1", "bear.tmin1", "wolf.tmin1", "wtd.tmin1", "lion.tmin1", "bear.t"), spp = c(".forest", ".bear", ".wolf", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup_inter[[54]] # NOT POSSIBLE - t cannot affect t-1 # "wolf.t"      "wtd.tmin1"   "moose.tmin1" "elk.tmin1"   "wolf.tmin1" 
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #42 
    # bs_bottomup_inter[[55]] # "wolf.t"       "wtd.t"        "moose.tmin1"  "elk.tmin1"    "wolf.tmin1"   "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    list(dSep_test = 7, covariates = c("moose.tmin1", "elk.tmin1", "wolf.tmin1", "wsi.tmin1", "forest.tmin1", "wtd.tmin1", "wolf.t"), spp = c(".moose", ".elk", ".wolf", ".wsi", ".forest", ".wtd", ".wolf"), indices = as.integer(c(1,1,1,1,1,1,2))),
    # Regression #43 
    # bs_bottomup_inter[[56]] # "wolf.t"      "coy.t"       "moose.tmin1" "elk.tmin1"   "wolf.tmin1"  "coy.tmin1"   "wtd.tmin1"
    list(dSep_test = 4, covariates = c("moose.tmin1", "elk.tmin1", "wolf.tmin1", "coy.tmin1", "wtd.tmin1", "wolf.t"), spp = c(".moose", ".elk", ".wolf", ".coy", ".wtd", ".wolf"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup_inter[[57]] # NOT POSSIBLE - t cannot affect t-1 # "wolf.t"      "lion.tmin1"  "moose.tmin1" "elk.tmin1"   "wolf.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #44 
    # bs_bottomup_inter[[58]] # "wolf.t"      "lion.t"      "moose.tmin1" "elk.tmin1"   "wolf.tmin1"  "wtd.tmin1"   "lion.tmin1" 
    list(dSep_test = 1, covariates = c("moose.tmin1", "elk.tmin1", "wolf.tmin1", "wtd.tmin1", "lion.tmin1", "wolf.t"), spp = c(".moose", ".elk", ".wolf", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1,1,2))),
    # Regression #45 
    # bs_bottomup_inter[[59]] # "wtd.t"        "coy.t"        "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "coy.tmin1"    "wolf.tmin1" 
    list(dSep_test = 4, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "coy.tmin1", "wolf.tmin1", "wtd.t"), spp = c(".wsi", ".forest", ".wtd", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup_inter[[60]] # NOT POSSIBLE - t cannot affect t-1 # "wtd.t"        "lion.tmin1"   "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #46 
    # bs_bottomup_inter[[61]] # "wtd.t"        "lion.t"       "wsi.tmin1"    "forest.tmin1" "wtd.tmin1"    "wolf.tmin1"   "lion.tmin1"
    list(dSep_test = 1, covariates = c("wsi.tmin1", "forest.tmin1", "wtd.tmin1", "wolf.tmin1", "lion.tmin1", "wtd.t"), spp = c(".wsi", ".forest", ".wtd", ".wolf", ".lion", ".wtd"), indices = as.integer(c(1,1,1,1,1,2))),
    #### bs_bottomup_inter[[62]] # NOT POSSIBLE - t cannot affect t-1 # "coy.t"      "lion.tmin1" "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # Regression #47 
    # bs_bottomup_inter[[63]] # "coy.t"      "lion.t"     "coy.tmin1"  "wolf.tmin1" "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("coy.tmin1", "wolf.tmin1", "wtd.tmin1", "lion.tmin1", "coy.t"), spp = c(".coy", ".wolf", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1,2)))
  )
  
  