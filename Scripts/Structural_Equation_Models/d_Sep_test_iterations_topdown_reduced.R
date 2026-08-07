  #'  ----------------------------------------------------------------
  #'  Active regression for d-Sep iterations: top-down
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
  
  dSep_iterations_topdown <- list(
    # bs_topdown[[1]] # "deerHarv.tmin1" "coy.t"          "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "deerHarv.tmin1"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1))),
    # bs_topdown[[2]] # "deerHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "deerHarv.tmin1"), spp = c(".harvest", ".bear", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[3]] # "deerHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "deerHarv.tmin1"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[4]] # "deerHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "deerHarv.tmin1"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[5]] # "deerHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "deerHarv.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown[[6]] # "deerHarv.tmin1" "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "deerHarv.tmin1"), spp = c(".harvest", ".lion", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[7]] # "wtd.tmin1" "coy.t"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wtd.tmin1"), spp = c(".coy", ".wtd"), indices = as.integer(c(1,1))),
    # bs_topdown[[8]] # "wtd.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "wtd.tmin1"), spp = c(".harvest", ".bear", ".wtd"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[9]] # "wtd.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "wtd.tmin1"), spp = c(".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[10]] # "wtd.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "wtd.tmin1"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[11]] # "wtd.tmin1"     "elk.t"         "elkHarv.tmin1" "elk.tmin1"     "wolf.tmin1"    "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "wtd.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".wtd"), indices = as.integer(c(1,1,1,1,1))),
    # bs_topdown[[12]] # "wtd.tmin1"      "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "wtd.tmin1"), spp = c(".harvest", ".lion", ".wtd"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[13]] # "moose.tmin1" "coy.t"       "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "moose.tmin1"), spp = c(".coy", ".moose"), indices = as.integer(c(1,1))),
    # bs_topdown[[14]] # "moose.tmin1"    "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "moose.tmin1"), spp = c(".harvest", ".bear", ".moose"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[15]] # "moose.tmin1"    "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "moose.tmin1"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[16]] # "moose.tmin1"    "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown[[17]] # "moose.tmin1"   "elk.t"         "elkHarv.tmin1" "elk.tmin1"     "wolf.tmin1"    "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,1))),
    # bs_topdown[[18]] # "moose.tmin1"    "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".harvest", ".lion", ".moose"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[19]] # "elkHarv.tmin1" "coy.t"         "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "elkHarv.tmin1"), spp = c(".coy", ".harvest"), indices = as.integer(c(1,1))),
    # bs_topdown[[20]] # "elkHarv.tmin1"  "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "elkHarv.tmin1"), spp = c(".harvest", ".bear", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[21]] # "elkHarv.tmin1" "moose.t"       "moose.tmin1"   "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "elkHarv.tmin1"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[22]] # "elkHarv.tmin1"  "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "elkHarv.tmin1"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[23]] # "elkHarv.tmin1"  "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "elkHarv.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown[[24]] # "elkHarv.tmin1"  "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "elkHarv.tmin1"), spp = c(".harvest", ".lion", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[25]] # "elk.tmin1" "coy.t"     "coy.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "elk.tmin1"), spp = c(".coy", ".elk"), indices = as.integer(c(1,1))),
    # bs_topdown[[26]] # "elk.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "elk.tmin1"), spp = c(".harvest", ".bear", ".elk"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[27]] # "elk.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1"), spp = c(".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[28]] # "elk.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "elk.tmin1"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[29]] # "elk.tmin1"      "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "elk.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown[[30]] # "elk.tmin1"      "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "elk.tmin1"), spp = c(".harvest", ".lion", ".elk"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[31]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "coy.tmin1"), spp = c(".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[32]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "coy.tmin1"), spp = c(".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[33]] # "coy.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "coy.tmin1"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[34]] # "coy.tmin1"      "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "coy.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown[[35]] # "coy.tmin1"     "elk.t"         "elkHarv.tmin1" "elk.tmin1"     "wolf.tmin1"    "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "coy.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".coy"), indices = as.integer(c(1,1,1,1,1))),
    # bs_topdown[[36]] # "coy.tmin1"      "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "coy.tmin1"), spp = c(".harvest", ".lion", ".coy"), indices = as.integer(c(1,1,1))),
    #### bs_topdown[[37]] #  NOT POSSIBLE EXCLUDE # "coy.t"          "bearHarv.tmin1" "coy.tmin1"
    ####list(dSep_test = NA, covariates = c("coy.tmin1", "coy.t"), spp = c(".coy", ".coy"), indices = as.integer(c(1,2))),
    #### bs_topdown[[38]] #  NOT POSSIBLE - t cannot affect t-1 # "coy.t"      "bear.tmin1" "coy.tmin1"
    ####list(dSep_test = 3, covariates = c("coy.tmin1", "coy.t"), spp = c(".coy", ".coy"), indices = as.integer(c(1,2))),
    # bs_topdown[[39]] # "coy.t"          "bear.t"         "coy.tmin1"      "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("coy.tmin1", "bearHarv.tmin1", "bear.tmin1", "coy.t"), spp = c(".coy", ".harvest", ".bear", ".coy"), indices = as.integer(c(1,1,1,2))),
    #### bs_topdown[[40]] #  NOT POSSIBLE EXCLUDE # "coy.t"          "wolfHarv.tmin1" "coy.tmin1"
    ####list(dSep_test = NA, covariates = c("coy.tmin1", ".coy"), spp = c(".coy", ".coy"), indices = as.integer(c(1,2))),
    #### bs_topdown[[41]] #  NOT POSSIBLE - t cannot affect t-1 # "coy.t"      "wolf.tmin1" "coy.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # bs_topdown[[42]] # "coy.t"       "moose.t"     "coy.tmin1"   "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("coy.tmin1", "moose.tmin1", "wolf.tmin1", "coy.t"), spp = c(".coy", ".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown[[43]] # "coy.t"          "wolf.t"         "coy.tmin1"      "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("coy.tmin1", "wolfHarv.tmin1", "wolf.tmin1", "coy.t"), spp = c(".coy", ".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1,2))),
    #### bs_topdown[[44]] #  NOT POSSIBLE EXCLUDE  # "coy.t"          "lionHarv.tmin1" "coy.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    #### bs_topdown[[45]] #  NOT POSSIBLE - t cannot affect t-1  # "coy.t"      "lion.tmin1" "coy.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # bs_topdown[[46]] # "coy.t"          "wtd.t"          "coy.tmin1"      "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("coy.tmin1", "deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "coy.t"), spp = c(".coy", ".harvest", ".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown[[47]] # "coy.t"         "elk.t"         "coy.tmin1"     "elkHarv.tmin1" "elk.tmin1"     "wolf.tmin1"    "lion.tmin1"
    list(dSep_test = 5, covariates = c("coy.tmin1", "elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "coy.t"), spp = c(".coy", ".harvest", ".elk", ".wolf", ".lion", ".coy"), indices = as.integer(c(1,1,1,1,1,2))),
    # bs_topdown[[48]] # "coy.t"          "lion.t"         "coy.tmin1"      "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("coy.tmin1", "lionHarv.tmin1", "lion.tmin1", "coy.t"), spp = c(".coy", ".harvest", ".lion", ".coy"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown[[49]] # "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "bearHarv.tmin1"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[50]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "bearHarv.tmin1"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[51]] # "bearHarv.tmin1" "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "bearHarv.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown[[52]] # "bearHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "bearHarv.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown[[53]] # "bearHarv.tmin1" "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "bearHarv.tmin1"), spp = c(".harvest", ".lion", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[54]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "bear.tmin1"), spp = c(".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[55]] # "bear.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "bear.tmin1"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[56]] # "bear.tmin1"     "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown[[57]] # "bear.tmin1"    "elk.t"         "elkHarv.tmin1" "elk.tmin1"     "wolf.tmin1"    "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".bear"), indices = as.integer(c(1,1,1,1,1))),
    # bs_topdown[[58]] # "bear.tmin1"     "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".harvest", ".lion", ".bear"), indices = as.integer(c(1,1,1))),
    #### bs_topdown[[59]] #  NOT POSSIBLE EXCLUDE  # "bear.t"         "wolfHarv.tmin1" "bearHarv.tmin1" "bear.tmin1"
    ####list(dSep_test = NA, covariates = c("wolfHarv.tmin1", "bearHarv.tmin1", "bear.tmin1", "bear.t"), spp = c(), indices = as.integer(c()))
    #### bs_topdown[[60]] #  NOT POSSIBLE - t cannot affect t-1  # "bear.t"         "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    #### list(dSep_test = , covariates = c(), spp = c(), indices = as.integer(c()))
    # bs_topdown[[61]] # "bear.t"         "moose.t"        "bearHarv.tmin1" "bear.tmin1"     "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("bearHarv.tmin1", "bear.tmin1", "moose.tmin1", "wolf.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown[[62]] # "bear.t"         "wolf.t"         "bearHarv.tmin1" "bear.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolfHarv.tmin1", "wolf.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,2,1,2))),
    #### bs_topdown[[63]] #  NOT POSSIBLE EXCLUDE  # "bear.t"         "lionHarv.tmin1" "bearHarv.tmin1" "bear.tmin1"
    #### list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()))
    #### bs_topdown[[64]] #  NOT POSSIBLE - t cannot affect t-1  # "bear.t"         "lion.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    #### list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()))
    # bs_topdown[[65]] # "bear.t"         "wtd.t"          "bearHarv.tmin1" "bear.tmin1"     "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("bearHarv.tmin1", "bear.tmin1", "deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".harvest", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,2,1,1,2))),
    # bs_topdown[[66]] # "bear.t"         "elk.t"          "bearHarv.tmin1" "bear.tmin1"     "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("bearHarv.tmin1", "bear.tmin1", "elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".harvest", ".elk", ".wolf", ".lion", ".bear"), indices = as.integer(c(1,1,2,1,1,1,2))),
    # bs_topdown[[67]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("bearHarv.tmin1", "bear.tmin1", "lionHarv.tmin1", "lion.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".harvest", ".lion", ".bear"), indices = as.integer(c(1,1,2,1,2))),
    # bs_topdown[[68]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "wolfHarv.tmin1"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[69]] # "wolfHarv.tmin1" "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "wolfHarv.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown[[70]] # "wolfHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "wolfHarv.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown[[71]] # "wolfHarv.tmin1" "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "wolfHarv.tmin1"), spp = c(".harvest", ".lion", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown[[72]] # "wolf.tmin1"     "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "wolf.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown[[73]] # "wolf.tmin1"     "lion.t"         "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "lion.tmin1", "wolf.tmin1"), spp = c(".harvest", ".lion", ".wolf"), indices = as.integer(c(1,1,1))),
    # bs_topdown[[74]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "wolfHarv.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2))),
    #### bs_topdown[[75]] #  NOT POSSIBLE EXCLUDE  # "moose.t"        "lionHarv.tmin1" "moose.tmin1"    "wolf.tmin1"
    #### list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()))
    #### bs_topdown[[76]] #  NOT POSSIBLE EXCLUDE  # "moose.t"     "lion.tmin1"  "moose.tmin1" "wolf.tmin1"
    #### list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()))
    # bs_topdown[[77]] # "moose.t"        "wtd.t"          "moose.tmin1"    "wolf.tmin1"     "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("moose.tmin1", "wolf.tmin1", "deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".harvest", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,1,2))),
    # bs_topdown[[78]] # "moose.t"       "elk.t"         "moose.tmin1"   "wolf.tmin1"    "elkHarv.tmin1" "elk.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("moose.tmin1", "wolf.tmin1", "elkHarv.tmin1", "elk.tmin1", "lion.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".harvest", ".elk", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,1,2))),
    # bs_topdown[[79]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("moose.tmin1", "wolf.tmin1", "lionHarv.tmin1", "lion.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".harvest", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    #### bs_topdown[[80]] #  NOT POSSIBLE EXCLUDE  # "wolf.t"         "lionHarv.tmin1" "wolfHarv.tmin1" "wolf.tmin1"
    #### list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()))
    #### bs_topdown[[81]] #  NOT POSSIBLE EXCLUDE  # "wolf.t"         "lion.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    #### list(dSep_test = NA, covariates = c(), spp = c(), indices = as.integer(c()))
    # bs_topdown[[82]] # "wolf.t"         "wtd.t"          "wolfHarv.tmin1" "wolf.tmin1"     "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".harvest", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,2,1,1,2))),
    # bs_topdown[[83]] # "wolf.t"         "elk.t"          "wolfHarv.tmin1" "wolf.tmin1"     "elkHarv.tmin1"  "elk.tmin1"      "lion.tmin1"
    list(dSep_test = 5, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "elkHarv.tmin1", "elk.tmin1", "lion.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".harvest", ".elk", ".lion", ".wolf"), indices = as.integer(c(1,1,2,1,1,2))),
    # bs_topdown[[84]] # "wolf.t"         "lion.t"         "wolfHarv.tmin1" "wolf.tmin1"     "lionHarv.tmin1" "lion.tmin1" 
    list(dSep_test = 1, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "lionHarv.tmin1", "lion.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".harvest", ".lion", ".wolf"), indices = as.integer(c(1,1,2,1,2))),
    # bs_topdown[[85]] # "lionHarv.tmin1" "wtd.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "lionHarv.tmin1"), spp = c(".harvest", ".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown[[86]] # "lionHarv.tmin1" "elk.t"          "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "lionHarv.tmin1"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown[[87]] # "wtd.t"          "elk.t"          "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"
    list(dSep_test = 5, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "wtd.t"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,2,1,1,2))),
    # bs_topdown[[88]] # "wtd.t"          "lion.t"         "deerHarv.tmin1" "wtd.tmin1"      "lion.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("deerHarv.tmin1", "wtd.tmin1", "lion.tmin1", "lionHarv.tmin1", "wtd.t"), spp = c(".harvest", ".wtd", ".lion", ".harvest", ".wtd"), indices = as.integer(c(1,1,1,2,2))),
    # bs_topdown[[89]] # "elk.t"          "lion.t"         "elkHarv.tmin1"  "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("elkHarv.tmin1", "elk.tmin1", "wolf.tmin1", "lion.tmin1", "lionHarv.tmin1", "elk.t"), spp = c(".harvest", ".elk", ".wolf", ".lion", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2,2)))
  )