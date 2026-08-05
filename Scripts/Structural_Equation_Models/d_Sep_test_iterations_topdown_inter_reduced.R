  #'  ----------------------------------------------------------------
  #'  Active regression for d-Sep iterations: top-down interference
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

  dSep_iterations_topdown_int <- list(
    # bs_topdown_inter[[1]] # wtd.tmin1 and moose.t are independent given moose.tmin1 and wolf.tmin1
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "wtd.tmin1"), spp = c(".moose", ".wolf", ".wtd"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[2]] # wtd.tmin1 and bear.t are independent given bearHarv.tmin1, bear.tmin1 and wolf.tmin1
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "wtd.tmin1"), spp = c(".harvest", ".bear", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[3]] # wtd.tmin1 and wolf.t are independent given wolfHarv.tmin1 and wolf.tmin1
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "wtd.tmin1"), spp = c(".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[4]] # wtd.tmin1 and elk.t are independent given elk.tmin1, wolf.tmin1 and lion.tmin1
    list(dSep_test = 5, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "wtd.tmin1"), spp = c(".elk", ".wolf", ".lion", ".wtd"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[5]] # "wtd.tmin1" and "coy.t" are independent given "coy.tmin1"  "wolf.tmin1" and "lion.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "wtd.tmin1"), spp = c(".coy", ".wolf", ".lion", ".wtd"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[6]] # "wtd.tmin1" and "lion.t" are independent given "lionHarv.tmin1" "wolf.tmin1" and "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "wolf.tmin1", "lion.tmin1", "wtd.tmin1"), spp = c(".harvest", ".wolf", ".lion", ".wtd"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[7]] # "moose.tmin1" and "bear.t" are independent given "bearHarv.tmin1" "bear.tmin1" and "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "moose.tmin1"), spp = c(".harvest", ".bear", ".wolf", ".moose"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[8]] # "moose.tmin1" and "wolf.t" are independent given "wolfHarv.tmin1" and "wolf.tmin1"
    list(dSep_test = 2, covaraites = c("wolfHarv.tmin1", "wolf.tmin1", "moose.tmin1"), spp = c(".harvest", ".wolf", ".moose"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[9]] # "moose.tmin1" and "wtd.t" are independent given "wtd.tmin1" and "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[10]] # "moose.tmin1" and "elk.t" are independent given "elk.tmin1" "wolf.tmin1" and "lion.tmin1"
    list(dSep_test = 5, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".elk", ".wolf", ".lion", ".moose"), indices = as.integer(c(1,1,1,1))),
    # # bs_topdown_inter[[11]] # "moose.tmin1" and "coy.t" are independent given "coy.tmin1"   "wolf.tmin1" and "lion.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".coy", ".wolf", ".lion", ".moose"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[12]] # "moose.tmin1" and "lion.t" are independent given "lionHarv.tmin1" "wolf.tmin1" and "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "wolf.tmin1", "lion.tmin1", "moose.tmin1"), spp = c(".harvest", ".wolf", ".lion", ".moose"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[13]] # "elk.tmin1" and "moose.t" are independent given "moose.tmin1" and "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1"), spp = c(".moose", ".wolf", ".elk"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[14]] # "elk.tmin1" and "bear.t" are independent given "bearHarv.tmin1" "bear.tmin1" and "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "elk.tmin1"), spp = c(".harvest", ".bear", ".wolf", ".elk"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[15]] "elk.tmin1" and "wolf.t" are independent given "wolfHarv.tmin1" and "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "elk.tmin1"), spp = c(".harvest", ".wolf", ".elk"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[16]] # "elk.tmin1" and "wtd.t" are independent given "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "lion.tmin1", "elk.tmin1"), spp = c(".wtd", ".lion", ".elk"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[17]] # "elk.tmin1"  "coy.t"      "coy.tmin1"  "wolf.tmin1" "lion.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "elk.tmin1"), spp = c(".coy", ".wolf", ".lion", ".elk"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[18]] # "elk.tmin1"      "lion.t"         "lionHarv.tmin1" "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "wolf.tmin1", "lion.tmin1", "elk.tmin1"), spp = c(".harvest", ".wolf", ".lion", ".elk"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[19]] # "coy.tmin1"   "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "coy.tmin1"), spp = c(".moose", ".wolf", ".coy"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[20]] # "coy.tmin1"      "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "coy.tmin1"), spp = c(".harvest", ".bear", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[21]] # "coy.tmin1"      "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "coy.tmin1"), spp = c(".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[22]] # "coy.tmin1"  "wtd.t"      "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.tmin1", "lion.tmin1", "coy.tmin1"), spp = c(".wtd", ".lion", ".coy"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[23]] # "coy.tmin1"  "elk.t"      "elk.tmin1"  "wolf.tmin1" "lion.tmin1"
    list(dSep_test = 5, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "coy.tmin1"), spp = c(".elk", ".wolf", ".lion", ".coy"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[24]] # "bearHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "bearHarv.tmin1"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[25]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "bearHarv.tmin1"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[26]] # "bearHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "bearHarv.tmin1"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2))), # notice different values for harvest indices
    # bs_topdown_inter[[27]] # "bearHarv.tmin1" "wtd.t"          "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.tmin1", "lion.tmin1", "bearHarv.tmin1"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[28]] # "bearHarv.tmin1" "elk.t"          "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "bearHarv.tmin1"), spp = c(".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[29]] # "bearHarv.tmin1" "coy.t"          "coy.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "bearHarv.tmin1"), spp = c(".coy", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[30]] # "bearHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "wolf.tmin1", "lion.tmin1", "bearHarv.tmin1"), spp = c(".harvest", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown_inter[[31]] # "bear.tmin1"  "moose.t"     "moose.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "bear.tmin1"), spp = c(".moose", ".wolf", ".bear"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[32]] # "bear.tmin1"     "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "bear.tmin1"), spp = c(".harvest", ".wolf", ".bear"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[33]] # "bear.tmin1" "wtd.t"      "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[34]] # "bear.tmin1" "elk.t"      "elk.tmin1"  "wolf.tmin1" "lion.tmin1"
    list(dSep_test = 5, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".elk", ".wolf", ".lion", ".bear"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[35]] # "bear.tmin1" "coy.t"      "coy.tmin1"  "wolf.tmin1" "lion.tmin1"
    list(dSep_test = 2, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".coy", ".wolf", ".lion", ".bear"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[36]] # "bear.tmin1"     "lion.t"         "lionHarv.tmin1" "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "wolf.tmin1", "lion.tmin1", "bear.tmin1"), spp = c(".harvest", ".wolf", ".lion", ".bear"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[37]] # "wolfHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1",  "wolf.tmin1", "wolfHarv.tmin1"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[38]] # "wolfHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "wolfHarv.tmin1"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown_inter[[39]] # "wolfHarv.tmin1" "wtd.t"          "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.tmin1", "lion.tmin1", "wolfHarv.tmin1"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[40]] # "wolfHarv.tmin1" "elk.t"          "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "wolfHarv.tmin1"), spp = c(".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[41]] # "wolfHarv.tmin1" "coy.t"          "coy.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "wolfHarv.tmin1"), spp = c(".coy", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[42]] # "wolfHarv.tmin1" "lion.t"         "lionHarv.tmin1" "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 1, covariates = c("lionHarv.tmin1", "wolf.tmin1", "lion.tmin1", "wolfHarv.tmin1"), spp = c(".harvest", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown_inter[[43]] # "lionHarv.tmin1" "moose.t"        "moose.tmin1"    "wolf.tmin1"
    list(dSep_test = 6, covariates = c("moose.tmin1", "wolf.tmin1", "lionHarv.tmin1"), spp = c(".moose", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[44]] # "lionHarv.tmin1" "bear.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    list(dSep_test = 3, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "lionHarv.tmin1"), spp = c(".harvest", ".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1,2))),
    # bs_topdown_inter[[45]] # "lionHarv.tmin1" "wolf.t"         "wolfHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 2, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "lionHarv.tmin1"), spp = c(".harvest", ".wolf", ".harvest"), indices = as.integer(c(1,1,2))),
    # bs_topdown_inter[[46]] # "lionHarv.tmin1" "wtd.t"          "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.tmin1", "lion.tmin1", "lionHarv.tmin1"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[47]] # "lionHarv.tmin1" "elk.t"          "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 5, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "lionHarv.tmin1"), spp = c(".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[48]] # "lionHarv.tmin1" "coy.t"          "coy.tmin1"      "wolf.tmin1"     "lion.tmin1"
    list(dSep_test = 4, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "lionHarv.tmin1"), spp = c(".coy", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[49]] # "wolf.tmin1" "wtd.t"      "wtd.tmin1"  "lion.tmin1"
    list(dSep_test = 7, covariates = c("wtd.tmin1", "lion.tmin1", "wolf.tmin1"), spp = c(".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[50]] # "moose.t"        "bear.t"         "moose.tmin1"    "wolf.tmin1"     "bearHarv.tmin1" "bear.tmin1"
    list(dSep_test = 3, covariates = c("moose.tmin1", "wolf.tmin1", "bearHarv.tmin1", "bear.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".harvest", ".bear", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[51]] # "moose.t"        "wolf.t"         "moose.tmin1"    "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("moose.tmin1", "wolf.tmin1", "wolfHarv.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".harvest", ".moose"), indices = as.integer(c(1,1,1,2))),
    #### bs_topdown_inter[[52]] # NOT POSSIBLE: t cannot affect t-1 # "moose.t"     "lion.tmin1"  "moose.tmin1" "wolf.tmin1"
    ####list(dSep_test = 1, covariates = c("moose.tmin1", "wolf.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".moose"), indices = as.integer(c(1,1,1))),
    # bs_topdown_inter[[53]] # "moose.t"     "wtd.t"       "moose.tmin1" "wolf.tmin1"  "wtd.tmin1"   "lion.tmin1"
    list(dSep_test = 7, covariates = c("moose.tmin1", "wolf.tmin1", "wtd.tmin1", "lion.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".wtd", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[54]] # "moose.t"     "elk.t"       "moose.tmin1" "wolf.tmin1"  "elk.tmin1"   "lion.tmin1"
    list(dSep_test = 5, covariates = c("moose.tmin1", "wolf.tmin1", "elk.tmin1", "lion.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".elk", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[55]] # "moose.t"     "coy.t"       "moose.tmin1" "wolf.tmin1"  "coy.tmin1"   "lion.tmin1"
    list(dSep_test = 4, covariates = c("moose.tmin1", "wolf.tmin1", "coy.tmin1", "lion.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".coy", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[56]] # "moose.t"        "lion.t"         "moose.tmin1"    "wolf.tmin1"     "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("moose.tmin1", "wolf.tmin1", "lionHarv.tmin1", "lion.tmin1", "moose.t"), spp = c(".moose", ".wolf", ".harvest", ".lion", ".moose"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[57]] # "bear.t"         "wolf.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "wolfHarv.tmin1"
    list(dSep_test = 2, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "wolfHarv.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1,2,2))),
    #### bs_topdown_inter[[58]] # NOT POSSIBLE: t cannot affect t-1 # "bear.t"         "lion.tmin1"     "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"
    ####list(dSep_test = 1, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1))),
    # bs_topdown_inter[[59]] # "bear.t"         "wtd.t"          "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "wtd.tmin1", "lion.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".wolf", ".wtd", ".lion", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    # bs_topdown_inter[[60]] # "bear.t"         "elk.t"          "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "elk.tmin1"      "lion.tmin1"
    list(dSep_test = 5, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "elk.tmin1", "lion.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".wolf", ".elk", ".lion", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    # bs_topdown_inter[[61]] # "bear.t"         "coy.t"          "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "coy.tmin1"      "lion.tmin1"
    list(dSep_test = 4, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "coy.tmin1", "lion.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".wolf", ".coy", ".lion", ".bear"), indices = as.integer(c(1,1,1,1,1,2))),
    # bs_topdown_inter[[62]] # "bear.t"         "lion.t"         "bearHarv.tmin1" "bear.tmin1"     "wolf.tmin1"     "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("bearHarv.tmin1", "bear.tmin1", "wolf.tmin1", "lionHarv.tmin1", "lion.tmin1", "bear.t"), spp = c(".harvest", ".bear", ".wolf", ".harvest", ".lion", ".bear"), indices = as.integer(c(1,1,1,2,1,2))),
    #### bs_topdown_inter[[63]] # NOT POSSIBLE: t cannot affect t-1 # "wolf.t"         "lion.tmin1"     "wolfHarv.tmin1" "wolf.tmin1"
    ###list(dSep_test = 1, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".wolf"), indices = as.integer(c(1,1,1)))
    # bs_topdown_inter[[64]] # "wolf.t"         "wtd.t"          "wolfHarv.tmin1" "wolf.tmin1"     "wtd.tmin1"      "lion.tmin1"
    list(dSep_test = 7, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "wtd.tmin1", "lion.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".wtd", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[65]] # "wolf.t"         "elk.t"          "wolfHarv.tmin1" "wolf.tmin1"     "elk.tmin1"      "lion.tmin1"
    list(dSep_test = 5, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "elk.tmin1", "lion.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".elk", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[66]] # "wolf.t"         "coy.t"          "wolfHarv.tmin1" "wolf.tmin1"     "coy.tmin1"      "lion.tmin1"
    list(dSep_test = 4, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "coy.tmin1", "lion.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".coy", ".lion", ".wolf"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[67]] # "wolf.t"         "lion.t"         "wolfHarv.tmin1" "wolf.tmin1"     "lionHarv.tmin1" "lion.tmin1"
    list(dSep_test = 1, covariates = c("wolfHarv.tmin1", "wolf.tmin1", "lionHarv.tmin1", "lion.tmin1", "wolf.t"), spp = c(".harvest", ".wolf", ".harvest", ".lion", ".wolf"), indices = as.integer(c(1,1,2,1,2))),
    # bs_topdown_inter[[68]] # "wtd.t"      "elk.t"      "wtd.tmin1"  "lion.tmin1" "elk.tmin1"  "wolf.tmin1"
    list(dSep_test = 5, covariates = c("wtd.tmin1", "lion.tmin1", "elk.tmin1", "wolf.tmin1", "wtd.t"), spp = c(".wtd", ".lion", ".elk", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[69]] # "wtd.t"      "coy.t"      "wtd.tmin1"  "lion.tmin1" "coy.tmin1"  "wolf.tmin1"
    list(dSep_test = 4, covariates = c("wtd.tmin1", "lion.tmin1", "coy.tmin1", "wolf.tmin1", "wtd.t"), spp = c(".wtd", ".lion", ".coy", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[70]] # "wtd.t"          "lion.t"         "wtd.tmin1"      "lion.tmin1"     "lionHarv.tmin1" "wolf.tmin1"
    list(dSep_test = 1, covariates = c("wtd.tmin1", "lion.tmin1", "lionHarv.tmin1", "wolf.tmin1", "wtd.t"), spp = c(".wtd", ".lion", ".harvest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[71]] # "elk.t"      "coy.t"      "elk.tmin1"  "wolf.tmin1" "lion.tmin1" "coy.tmin1"
    list(dSep_test = 4, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "coy.tmin1", "elk.t"), spp = c(".elk", ".wolf", ".lion", ".coy", ".elk"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[72]] # "elk.t"          "lion.t"         "elk.tmin1"      "wolf.tmin1"     "lion.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "lionHarv.tmin1", "elk.t"), spp = c(".elk", ".wolf", ".lion", ".harvest", ".elk"), indices = as.integer(c(1,1,1,1,2))),
    # bs_topdown_inter[[73]] # "coy.t"          "lion.t"         "coy.tmin1"      "wolf.tmin1"     "lion.tmin1"     "lionHarv.tmin1"
    list(dSep_test = 1, covariates = c("coy.tmin1", "wolf.tmin1", "lion.tmin1", "lionHarv.tmin1", "coy.t"), spp = c(".coy", ".wolf", ".lion", ".harvest", ".coy"), indices = as.integer(c(1,1,1,1,2)))
  )
  