  #'  ----------------------------------------------------------------
  #'  Inputs for for d-Sep iterations for cor(exog, exog): bottom-up
  #'  Aug 2026
  #'  ----------------------------------------------------------------
  #'  List inputs for active regressions to iterate over for each d-separation 
  #'  test when testing correlation between two exogenous variables flagged in 
  #'  the basic set. Because these are predictor variables, we do not have to 
  #'  worry about dealing with the latent process or propagating uncertainty.
  #'  It is also worth noting that there is not causal ordering implied here. 
  #'  This is simply to test for correlations between these variables to be 
  #'  consistent with the dSep process and to incorporate into the Fisher's C 
  #'  statistic if desired.
  #'  
  #'  Input for each iteration are: 
  #'  y_array = observed predictor variable 
  #'  x_array = observed predictor variable 
  #'  
  #'  This list of inputs differs from the input structures of the other two active 
  #'  regression scripts because we are not building an entirely new model with 
  #'  each iteration here. Instead, we are simply testing the correlation between 
  #'  two exogenous variables included in the larger SEM. 
  #'  ----------------------------------------------------------------
  
  dSep_iterations_bottomup_exog_only <- list(
    # Regression #2
    # bs_topdown[[2]] # "wsi.tmin1"    "forest.tmin1"
    list(y_array = data_JAGS_bundle_topdown$forest, x_array = data_JAGS_bundle_topdown$wsi)
  )
  

  