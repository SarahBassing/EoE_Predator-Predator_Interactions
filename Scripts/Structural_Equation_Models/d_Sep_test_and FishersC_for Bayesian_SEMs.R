  #'  -----------------------------------
  #'  D-separation tests with ROPE method
  #'  Sarah Bassing & Matt Falcy
  #'  July 2026
  #'  -----------------------------------
  #'  Using Region of Practical Equivalence (ROPE) to test d-sep claims. This is 
  #'  a bit of a hack to generate pseudo p-values for Fisher's C. Method adapted
  #'  from Kruschke and Liddell (2018), who used "0.1" of 1 SD to define a ROPE:
  #'  https://link.springer.com/article/10.3758/s13423-016-1221-4
  #'  
  #'  Script first formats all data by sourcing scripts to format covariates and
  #'  posteriors from the Royle-Nichols model analyses, then formats and bundles
  #'  input data to fit SEMs in JAGS. Data bundling, initial values, and parameters
  #'  to monitor are general to allow for my dynamism with the d-separation tests.
  #'  This means some data, inits, and parameters will not be used in a given model.
  #'  This produces warnings that can be ignored.
  #'  
  #'  Script then conducts d-separation tests iteratively It first defines the 
  #'  ROPE function and identifies the basic set for each SEM. It then sources a
  #'  template JAGS model and updates the template for each individual d-sep test
  #'  by creating a set of custom regressions in the JAGS template for each test*.
  #'  It then fits each iteration* of the model and calculates a p-value using the
  #'  ROPE method.
  #'  
  #'  *Coding assistance and trouble shooting conducted with the help of Claude.ai
  #'  -----------------------------------
  
  #'  Load libraries
  # install.packages("BiocManager") # Needed for 'graph' which is needed by ggm
  # BiocManager::install("graph")
  library(ggm) # for DAG() and basiSet(), which retrieves the basis set
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  library(future.apply)
  
  #'  -----------------------------------
  ####  Format data and bundle for JAGS  ####
  #'  -----------------------------------
  #'  Formats covariate data
  source("./Scripts/Structural_Equation_Models/Format_spatial_covariates_for_SEMs.R")
  #'  Formats density data for SEMs
  source("./Scripts/Structural_Equation_Models/Format_RNmodel_Posteriors_for_SEM.R")
  #'  Basic set and d-Separation test functions
  source("./Scripts/Structural_Equation_Models/d_Sep_functions.R")
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  #' #'  Source functions to setup data for JAGS
  #' source("./Scripts/Structural_Equation_Models/JAGS_setup_functions.R")
  
  #'  Parameters to be monitored
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk", 
              "beta.moose", "beta.wtd", "beta.harvest", "beta.wsi","beta.forest", "sigma.spp")  
  
  #'  MCMC settings
  nc <- 3
  ni <- 100000
  nb <- 50000
  nt <- 10
  na <- 5000
  
  #'  Run all iterations in parallel
  plan(multisession, workers = parallel::detectCores() - 1)
  
  #'  -------------------------------
  ####  Basic set for original SEMs  ####
  #'  -------------------------------
  #'  -------------------
  #####  Top-down model  #####
  #'  -------------------
  #'  Generate DAG
  dag_topdown <- DAG(lion.t ~ lionHarv.tmin1, # lion.tmin1 + 
                     wolf.t ~ wolf.tmin1 + wolfHarv.tmin1,
                     bear.t ~ bear.tmin1 + bearHarv.tmin1,
                     coy.t ~ coy.tmin1,
                     elk.t ~ elk.tmin1 + wolf.tmin1 + lion.tmin1 + elkHarv.tmin1,
                     moose.t ~ moose.tmin1 + wolf.tmin1,
                     wtd.t ~ wtd.tmin1 + lion.tmin1 + deerHarv.tmin1)
  
  #'  Generate basic set
  bs_topdown <- basic_set(dag_topdown)
  
  #'  --------------------------------
  #####  Top-down interference model  #####
  #'  --------------------------------
  dag_topdown_inter <- DAG(lion.t ~ wolf.tmin1 + lionHarv.tmin1, # lion.tmin1 + 
                           wolf.t ~ wolf.tmin1 + wolfHarv.tmin1,
                           bear.t ~ bear.tmin1 + wolf.tmin1 + bearHarv.tmin1,
                           coy.t ~ coy.tmin1 + wolf.tmin1 + lion.tmin1,
                           elk.t ~ elk.tmin1 + wolf.tmin1 + lion.tmin1,
                           moose.t ~ moose.tmin1 + wolf.tmin1,
                           wtd.t ~ wtd.tmin1 + lion.tmin1)
  
  bs_topdown_inter <- basic_set(dag_topdown_inter)
  
  #'  --------------------
  #####  Bottom-up model  ####
  #'  --------------------
  dag_bottomup <- DAG(lion.t ~ elk.tmin1 + wtd.tmin1, #lion.tmin1 + 
                      wolf.t ~ wolf.tmin1 + elk.tmin1 + moose.tmin1,
                      bear.t ~ bear.tmin1 + elk.tmin1 + forest.tmin1,
                      coy.t ~ coy.tmin1 + wtd.tmin1,
                      elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1,
                      moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                      wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1)
  bs_bottomup <- basic_set(dag_bottomup)
  
  #'  ---------------------------------
  #####  Bottom-up interference model  ####
  #'  ---------------------------------
  dag_bottomup_inter <- DAG(lion.t ~ elk.tmin1 + wtd.tmin1, #+ wolf.tmin1 lion.tmin1 + 
                            wolf.t ~ wolf.tmin1 + elk.tmin1 + moose.tmin1,
                            bear.t ~ bear.tmin1 + elk.tmin1 + forest.tmin1 + wolf.tmin1,
                            coy.t ~ coy.tmin1 + wtd.tmin1 + wolf.tmin1,
                            elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1,
                            moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                            wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1)
  
  bs_bottomup_inter <- basic_set(dag_bottomup_inter)
  
  #'  -------------------------------
  ####  Iterate through d-Sep tests  ####
  #'  -------------------------------
  #'  ------------------------------
  #####  Top-down model iterations  #####
  #'  ------------------------------
  #'  Fit independence claims for variables where t-1 --> t or t --> t
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.latent
    list(covs = c("lionHarv"), spp = c(".harvest"), indices = as.integer(c(1)), lags = c("y-1")),
    #'  Regression 2: wolf.latent
    list(covs = c("wolf.latent", "wolfHarv"), spp = c(".wolf", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 3: bear.latent
    list(covs = c("bear.latent", "bearHarv"), spp = c(".bear", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 4: coy.latent
    list(covs = c("coy.latent"), spp = c(".coy"), indices = as.integer(c(1)), lags = c("y-1")),
    #'  Regression 5: elk.latent
    list(covs = c("elk.latent", "wolf.latent", "lion.latent", "elkHarv"), spp = c(".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    #'  Regression 6: moose.latent
    list(covs = c("moose.latent", "wolf.latent"), spp = c(".moose", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 7: wtd.latent
    list(covs = c("wtd.latent", "lion.latent", "deerHarv"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1"))
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown.R")
  
  #'  Bundle data and draw inits using functions in in Format_RNmodel_Posteriors_for_SEM.R
  data_JAGS_bundle_topdown <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                         dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                         covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                         covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                         nwolf = 4, nlion = 3, nbear = 2, ncoy = 2, nelk = 2, 
                                         nmoose = 2, nwtd = 2, nharv = 6, nfor = 0, nwsi = 0)
  num.chains <- 3
  initsList_topdown <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_topdown[[i]] <- generate_inits(nwolf = 4, nlion = 3, nbear = 2, ncoy = 2, nelk = 2, nmoose = 2, 
                                             nwtd = 2, nharv = 6, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  #'  Fit and save model iterations
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_topdown, listInits = initsList_topdown, model_name = "TopDown"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_tmin1_only.R")
  
  ### MAKE SURE SEM_TOPDOWN IS IN WORKING DIRECTORY  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_tmin1_only),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_topdown_tmin1_only, 
                              og_fit = SEM_topdown, nSites = 23, nYear = 4, model_name = "TopDown",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_exog_only.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_exog_only),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_topdown_exog_only, 
                                    model_name = "TopDown_exog", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  -------------------------------------------
  #####  Top-down interference model iterations  #####
  #'  -------------------------------------------
  #'  Fit independence claims for variables where t-1 --> t or t --> t
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.latent
    list(covs = c("lionHarv", "wolf.latent"), spp = c(".harvest", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")), 
    #'  Regression 2: wolf.latent
    list(covs = c("wolf.latent", "wolfHarv"), spp = c(".wolf", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 3: bear.latent
    list(covs = c("bear.latent", "bearHarv", "wolf.latent"), spp = c(".bear", ".harvest", ".wolf"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 4: coy.latent
    list(covs = c("coy.latent", "wolf.latent", "lion.latent"), spp = c(".coy", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 5: elk.latent
    list(covs = c("elk.latent", "wolf.latent", "lion.latent"), spp = c(".elk", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 6: moose.latent
    list(covs = c("moose.latent", "wolf.latent"), spp = c(".moose", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 7: wtd.latent
    list(covs = c("wtd.latent", "lion.latent"), spp = c(".wtd", ".lion"), indices = as.integer(c(1,1)), lags = c("y-1","y-1"))
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_inter.R")
  
  data_JAGS_bundle_topdown_inter <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                               dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                               covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                               covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                               nwolf = 7, nlion = 4, nbear = 2, ncoy = 2, nelk = 2, 
                                               nmoose = 2, nwtd = 2, nharv = 4, nfor = 0, nwsi = 0)
  num.chains <- 3
  initsList_topdown_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_topdown_inter[[i]] <- generate_inits(nwolf = 7, nlion = 4, nbear = 2, ncoy = 2, nelk = 2, nmoose = 2, 
                                                   nwtd = 2, nharv = 4, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  #'  Fit and save model iterations
  start.time = Sys.time()
  saved_paths <- future_lapply(
    #'  Apply across every element in list of active regressions
    seq_along(dSep_iterations_topdown_int),
    #'  Call run_dSep_iterations function using specified active regression list, model template, and data/inits prepared for JAGS
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown_int, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_topdown_inter, listInits = initsList_topdown_inter, model_name = "TopDown_Interference"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_inter_tmin1_only.R")
  
  ### MAKE SURE SEM_TOPDOWN_INTER IS IN WORKING DIRECTORY  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_inter_tmin1_only),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_topdown_inter_tmin1_only, 
                              og_fit = SEM_topdown_inter, nSites = 23, nYear = 4, model_name = "TopDown_Interference",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_inter_exog_only.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_inter_exog_only),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_topdown_inter_exog_only, 
                                    model_name = "TopDown_Interference_exog", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  -------------------------------
  #####  Bottom-up model iterations  #####
  #'  -------------------------------
  #'  Fit independence claims for variables where t-1 --> t or t --> t
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.t
    list(covs = c("elk.latent", "wtd.latent"), spp = c(".elk", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 2: wolf.t
    list(covs = c("wolf.latent", "elk.latent", "moose.latent"), spp = c(".wolf", ".elk", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 3: bear.t
    list(covs = c("bear.latent", "elk.latent", "forest"), spp = c(".bear", ".elk", ".forest"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 4: coy.t
    list(covs = c("coy.latent", "wtd.latent"), spp = c(".coy", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 5: elk.t
    list(covs = c("elk.latent", "forest", "wsi"), spp = c(".elk", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 6: moose.t
    list(covs = c("moose.latent", "forest", "wsi"), spp = c(".moose", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 7: wtd.t
    list(covs = c("wtd.latent", "forest", "wsi"), spp = c(".wtd", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1"))
    )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup.R")
  
  #'  Bundle data and draw inits using functions in in Format_RNmodel_Posteriors_for_SEM.R
  data_JAGS_bundle_bottomup <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                          dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                          covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                          covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                          nwolf = 2, nlion = 1, nbear = 2, ncoy = 2, nelk = 5, 
                                          nmoose = 3, nwtd = 4, nharv = 0, nfor = 5, nwsi = 5)
                                         
  num.chains <- 3
  initsList_bottomup <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup[[i]] <- generate_inits(nwolf = 2, nlion = 1, nbear = 2, ncoy = 2, nelk = 5, nmoose = 3, 
                                             nwtd = 4, nharv = 0, nfor = 5, nwsi = 5, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_bottomup, listInits = initsList_bottomup, model_name = "BottomUp"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_tmin1_only.R")
  
  ### MAKE SURE SEM_BOTTOMUP IS IN WORKING DIRECTORY  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_tmin1_only),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_bottomup_tmin1_only, 
                              og_fit = SEM_bottomup, nSites = 23, nYear = 4, model_name = "BottomUp",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_exog_only.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_exog_only),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_bottomup_exog_only, 
                                    model_name = "BottomUp_exog", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.latent
    list(covs = c("elk.latent", "wtd.latent"), spp = c(".elk", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 2: wolf.latent
    list(covs = c("wolf.latent", "elk.latent", "moose.latent"), spp = c(".wolf", ".elk", ".moose"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 3: bear.latent
    list(covs = c("bear.latent", "elk.latent", "forest", "wolf.latent"), spp = c(".bear", ".elk", ".forest", ".wolf"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 4: coy.latent
    list(covs = c("coy.latent", "wtd.latent", "wolf.latent"), spp = c(".coy", ".wtd", ".wolf"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 5: elk.latent
    list(covs = c("elk.latent", "forest", "wsi"), spp = c(".elk", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 6: moose.latent
    list(covs = c("moose.latent", "forest", "wsi"), spp = c(".moose", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 7: wtd.latent
    list(covs = c("wtd.latent", "forest", "wsi"), spp = c(".wtd", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1"))
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_inter_updated.R")
  
  #'  Bundle data and draw inits using functions in in Format_RNmodel_Posteriors_for_SEM.R
  data_JAGS_bundle_bottomup_inter <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s,
                                                covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                nwolf = 4, nlion = 1, nbear = 2, ncoy = 2, nelk = 5, 
                                                nmoose = 3, nwtd = 4, nharv = 0, nfor = 5, nwsi = 4)
  num.chains <- 3
  initsList_bottomup_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup_inter[[i]] <- generate_inits(nwolf = 4, nlion = 1, nbear = 2, ncoy = 2, nelk = 5, nmoose = 3, 
                                                    nwtd = 4, nharv = 0, nfor = 5, nwsi = 4, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_inter), 
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup_inter, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_bottomup_inter, listInits = initsList_bottomup_inter, model_name = "BottomUp_Interference"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_inter_tmin1_only_updated.R")
  
  ### MAKE SURE SEM_BOTTOMUP_INTER IS IN WORKING DIRECTORY  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_inter_tmin1_only),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_bottomup_inter_tmin1_only, 
                              og_fit = SEM_bottomup_inter, nSites = 23, nYear = 4, model_name = "BottomUp_Interference",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_inter_exog_only.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_inter_exog_only),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_bottomup_inter_exog_only, 
                                    model_name = "BottomUp_Interference_exog", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  -------------------------------------
  ####  Calculate p-values and Fisher's C  ####
  #'  -------------------------------------
  #'  ------------------------------
  #####  Top-down model iterations  #####
  #'  ------------------------------
  #'  Load iterations of the JAGS model
  #'  Note: this list of outputs is based on the number of independence claims 
  #'  assessed using the d_Sep_active_regression_topdown.R active regression list. 
  #'  This is not a complete list of all independence claims being tested for this SEM.
  all_results_topdown <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/TopDown", full.names = TRUE), readRDS)
  
  #'  Rename data bundle
  data_JAGS_bundle <- data_JAGS_bundle_topdown
  #'  Create list of "observed" values of focal response variable, one per d-Sep test               # instances where x was used as y in d-Sep test noted below
  y_list <- list(data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,   #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,    #
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,   #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,    #
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat,     #     
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,    #
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat,     #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,  #
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,   #
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$wtd.hat,    # wtd flipped
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$wtd.hat,      # wtd flipped both times      
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wtd.hat,     # wtd fillped both times
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,    # wtd flipped
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$lion.hat,    # wtd flipped
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,   #
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat,   #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$coy.hat,     # coy flipped both times
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$coy.hat,     # coy flipped both times
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat,   # coy flipped
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,    #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$bear.hat,   # bear flipped both times
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat,  # bear flipped
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$lion.hat,    #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$moose.hat,   # moose flipped
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$elk.hat,    # elk flipped
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat)   # wolf flipped
  #'  Results in 90 d-Sep tests that were t-1 --> t or t --> t
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta.harvest indices. 
  #'  Most will be [,1] but some will be [,2] where the same beta name was used 
  #'  twice in the same regression.
  mod_out <- list()
  for(i in 1:length(all_results_topdown)) {
    mod_out[[i]] <- all_results_topdown[[i]]$fit$sims.list
  }
  #'  Note: indexing is based on the number of independence claims assessed using the
  #'  d_Sep_active_regression_topdown.R active regression list. This is not a complete
  #'  list of all independence claims being tested.
  post_list_topdown <- list(mod_out[[1]]$beta.harvest[,1], mod_out[[2]]$beta.harvest[,2], mod_out[[3]]$beta.harvest[,1],     # note indexing
                            mod_out[[4]]$beta.harvest[,2], mod_out[[5]]$beta.harvest[,2], mod_out[[6]]$beta.harvest[,2],     # note indexing
                            mod_out[[7]]$beta.wtd[,1], mod_out[[8]]$beta.wtd[,1], mod_out[[9]]$beta.wtd[,1],
                            mod_out[[10]]$beta.wtd[,1], mod_out[[11]]$beta.wtd[,1], mod_out[[12]]$beta.wtd[,1],
                            mod_out[[13]]$beta.moose[,1], mod_out[[14]]$beta.moose[,1], mod_out[[15]]$beta.moose[,1],
                            mod_out[[16]]$beta.moose[,1], mod_out[[17]]$beta.moose[,1], mod_out[[18]]$beta.moose[,1],
                            mod_out[[19]]$beta.harvest[,2], mod_out[[20]]$beta.harvest[,1], mod_out[[21]]$beta.harvest[,1],  # note indexing
                            mod_out[[22]]$beta.harvest[,1], mod_out[[23]]$beta.harvest[,2], mod_out[[24]]$beta.harvest[,2],  # note indexing
                            mod_out[[25]]$beta.lion[,1], mod_out[[26]]$beta.lion[,1], mod_out[[27]]$beta.lion[,1], 
                            mod_out[[28]]$beta.lion[,1], mod_out[[29]]$beta.lion[,1], mod_out[[30]]$beta.elk[,1],            # note elk flipped
                            mod_out[[31]]$beta.coy[,1], mod_out[[32]]$beta.wtd[,2], mod_out[[33]]$beta.harvest[,2],          # note indexing and coy & harvest both flipped
                            mod_out[[34]]$beta.bear[,1], mod_out[[35]]$beta.wtd[,2], mod_out[[36]]$beta.harvest[,2],         # note indexing and bear & harvest both flipped
                            mod_out[[37]]$beta.harvest[,1], mod_out[[38]]$beta.wtd[,2], mod_out[[39]]$beta.wtd[,2],          # note indexing and harvest flipped              
                            mod_out[[40]]$beta.wtd[,2], mod_out[[41]]$beta.harvest[,2], mod_out[[42]]$beta.wtd[,2],          # note indexing and harvest flipped
                            mod_out[[43]]$beta.elk[,1], mod_out[[44]]$beta.elk[,1], mod_out[[45]]$beta.elk[,1],  
                            mod_out[[46]]$beta.elk[,1], mod_out[[47]]$beta.elk[,1], mod_out[[48]]$beta.coy[,1],              
                            mod_out[[49]]$beta.coy[,1], mod_out[[50]]$beta.coy[,1], mod_out[[51]]$beta.coy[,1],              
                            mod_out[[52]]$beta.coy[,1], mod_out[[53]]$beta.harvest[,1], mod_out[[54]]$beta.coy[,1],          # note harvest flipped
                            mod_out[[55]]$beta.coy[,2], mod_out[[56]]$beta.harvest[,1], mod_out[[57]]$beta.wolf[,1],         # note indexing and harvest & wolf both flipped  
                            mod_out[[58]]$beta.coy[,2], mod_out[[59]]$beta.coy[,2], mod_out[[60]]$beta.coy[,2],              # note indexing
                            mod_out[[61]]$beta.harvest[,1], mod_out[[62]]$beta.coy[,2], mod_out[[63]]$beta.harvest[,1],      # note indexing and 1st harvest flipped
                            mod_out[[64]]$beta.harvest[,2], mod_out[[65]]$beta.harvest[,2], mod_out[[66]]$beta.harvest[,2],  # note indexing
                            mod_out[[67]]$beta.bear[,1], mod_out[[68]]$beta.bear[,1], mod_out[[69]]$beta.bear[,1],           
                            mod_out[[70]]$beta.bear[,1], mod_out[[71]]$beta.harvest[,2], mod_out[[72]]$beta.wolf[,1],        # note indexing and harvest & wolf both flipped
                            mod_out[[73]]$beta.bear[,2], mod_out[[74]]$beta.bear[,2], mod_out[[75]]$beta.bear[,2],           # note indexing
                            mod_out[[76]]$beta.harvest[,2], mod_out[[77]]$beta.bear[,2], mod_out[[78]]$beta.harvest[,1],     # note indexing and 1st harvest flipped
                            mod_out[[79]]$beta.harvest[,2], mod_out[[80]]$beta.harvest[,2], mod_out[[81]]$beta.wolf[,1],     # note indexing
                            mod_out[[82]]$beta.moose[,2], mod_out[[83]]$beta.moose[,2], mod_out[[84]]$beta.harvest[,1],      # note indexing and harvest flipped 
                            mod_out[[85]]$beta.moose[,2], mod_out[[86]]$beta.elk[,2], mod_out[[87]]$beta.harvest[,2],        # note indexing and harvest flipped
                            mod_out[[88]]$beta.elk[,2], mod_out[[89]]$beta.harvest[,2], mod_out[[90]]$beta.wolf[,2])         # note indexing and harvest flipped
  
  #'  Load more iterations of the JAGS model
  #'  Note: this list of outputs is based on the number of independence claims 
  #'  assessed using the d_Sep_active_regression_topdown_tmin1_only.R list. 
  #'  This is not a complete list of all independence claims being tested.
  all_results_topdown_tmin1 <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/TopDown", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test                 # instances where x was used as y in d-Sep test noted below
  y_list2 <- list(data_JAGS_bundle$wtd.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$lion.hat,
                  data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat,
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wolf.hat,
                  data_JAGS_bundle$lion.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,
                  data_JAGS_bundle$wtd.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wtd.hat,
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$moose.hat,
                  data_JAGS_bundle$lion.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat, 
                  data_JAGS_bundle$moose.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$lion.hat,
                  data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat,
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,
                  data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat,
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat,
                  data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$elk.hat,
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,
                  data_JAGS_bundle$bear.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$wolf.hat,
                  data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,
                  data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$bear.hat,
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wolf.hat)
  #'  Leaves you with 56 d-Sep tests that were t-1 --> t-1
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Generally don't need to worry about indexing here because most of these independence
  #'  claims are really marginal independence claims (not conditional ind. claims) 
  #'  so there is typically only 1 explanatory variable to consider and index will be [1]
  mod_out2 <- list()
  for(i in 1:length(all_results_topdown_tmin1)) {
    mod_out2[[i]] <- all_results_topdown_tmin1[[i]]$fit$sims.list
  }
  post_list_topdown_tmin1 <- list(mod_out2[[1]]$b_x, mod_out2[[2]]$b_x,  mod_out2[[3]]$b_x,
                                  mod_out2[[4]]$b_x,  mod_out2[[5]]$b_x,  mod_out2[[6]]$b_x,
                                  mod_out2[[7]]$b_x,  mod_out2[[8]]$b_x,  mod_out2[[9]]$b_x,
                                  mod_out2[[10]]$b_x,  mod_out2[[11]]$b_x,  mod_out2[[12]]$b_x, 
                                  mod_out2[[13]]$b_x,  mod_out2[[14]]$b_x,  mod_out2[[15]]$b_x, 
                                  mod_out2[[16]]$b_x,  mod_out2[[17]]$b_x,  mod_out2[[18]]$b_x,
                                  mod_out2[[19]]$b_x,  mod_out2[[20]]$b_x,  mod_out2[[21]]$b_x, 
                                  mod_out2[[22]]$b_x,  mod_out2[[23]]$b_x,  mod_out2[[24]]$b_x,
                                  mod_out2[[25]]$b_x,  mod_out2[[26]]$b_x,  mod_out2[[27]]$b_x,
                                  mod_out2[[28]]$b_x,  mod_out2[[29]]$b_x,  mod_out2[[30]]$b_x,
                                  mod_out2[[31]]$b_x,  mod_out2[[32]]$b_x,  mod_out2[[33]]$b_x,
                                  mod_out2[[34]]$b_x,  mod_out2[[35]]$b_x,  mod_out2[[36]]$b_x, 
                                  mod_out2[[37]]$b_x,  mod_out2[[38]]$b_x,  mod_out2[[39]]$b_x, 
                                  mod_out2[[40]]$b_x,  mod_out2[[41]]$b_x,  mod_out2[[42]]$b_x,
                                  mod_out2[[43]]$b_x,  mod_out2[[44]]$b_x,  mod_out2[[45]]$b_x, 
                                  mod_out2[[46]]$b_x,  mod_out2[[47]]$b_x,  mod_out2[[48]]$b_x, 
                                  mod_out2[[49]]$b_x, mod_out2[[50]]$b_x, mod_out2[[51]]$b_x, 
                                  mod_out2[[52]]$b_x, mod_out2[[53]]$b_x, mod_out2[[54]]$b_x,
                                  mod_out2[[55]]$b_x, mod_out2[[56]]$b_x)
  
  #'  Load more iterations of the JAGS model (this time assessing correlation between exogenous variables)
  all_results_topdown_exog <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/TopDown_exog", full.names = TRUE), readRDS)
  y_list3 <- list(data_JAGS_bundle$elkHarv, data_JAGS_bundle$bearHarv, data_JAGS_bundle$wolfHarv,
                  data_JAGS_bundle$lionHarv, data_JAGS_bundle$bearHarv, data_JAGS_bundle$wolfHarv,
                  data_JAGS_bundle$lionHarv, data_JAGS_bundle$wolfHarv, data_JAGS_bundle$lionHarv,
                  data_JAGS_bundle$wolfHarv) # y_array
  #'  Leaves you with 10 d-Sep tests that were exog --> exog
  
  mod_out3 <- list()
  for(i in 1:length(all_results_topdown_exog)) {
    mod_out3[[i]] <- all_results_topdown_exog[[i]]$fit$sims.list
  }
  post_list_topdown_exog <- list(mod_out3[[1]]$b_x, mod_out3[[2]]$b_x, mod_out3[[3]]$b_x,
                                 mod_out3[[4]]$b_x, mod_out3[[5]]$b_x, mod_out3[[6]]$b_x,
                                 mod_out3[[7]]$b_x, mod_out3[[8]]$b_x, mod_out3[[9]]$b_x,
                                 mod_out3[[10]]$b_x) # correspond to x_array
  
  
  #'  -------------------------
  ######  ROPE method p-value  ######
  #'  -------------------------
  #'  Calculate p.rope value for each iteration of the d-Sep test
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_topdown_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_topdown, SIMPLIFY = FALSE)
  p.rope_topdown_tmin1_list <- mapply(p_rope_iterations, y_dat = y_list2, post_beta = post_list_topdown_tmin1, SIMPLIFY = FALSE)
  p.rope_topdown_exog_list <- mapply(p_rope_iterations, y_dat = y_list3, post_beta = post_list_topdown_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_topdown_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_topdown_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_topdown_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_topdown_tmin1_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_topdown_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_topdown_exog_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_topdown_df <- stack(p.rope_topdown_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "normal")
  p.rope_topdown_tmin1_df <- stack(p.rope_topdown_tmin1_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "tmin1")
  p.rope_topdown_exog_df <- stack(p.rope_topdown_exog_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "exog")
  
  #'  ----------------------
  ######  Bayesian p-value  ######
  #'  ----------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_topdown_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown, SIMPLIFY = FALSE)
  bayes.p_topdown_tmin1_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown_tmin1, SIMPLIFY = FALSE)
  bayes.p_topdown_exog_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_topdown_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_topdown_list)[i] <- list_name
  }
  bayes.p_topdown_df <- stack(bayes.p_topdown_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "normal")
  
  for(i in 1:length(bayes.p_topdown_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_topdown_tmin1_list)[i] <- list_name
  }
  bayes.p_topdown_tmin1_df <- stack(bayes.p_topdown_tmin1_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "tmin1")
  
  for(i in 1:length(bayes.p_topdown_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_topdown_exog_list)[i] <- list_name
  }
  bayes.p_topdown_exog_df <- stack(bayes.p_topdown_exog_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "exog")
  
  #'  Join both d-Sep test p-values and save
  p.val_topdown_df <- full_join(p.rope_topdown_df, bayes.p_topdown_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_topdown_tmin1_df <- full_join(p.rope_topdown_tmin1_df, bayes.p_topdown_tmin1_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_topdown_exog_df <- full_join(p.rope_topdown_exog_df, bayes.p_topdown_exog_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_topdown_all_df <- bind_rows(p.val_topdown_df, p.val_topdown_tmin1_df, p.val_topdown_exog_df)
  
  write_csv(p.val_topdown_all_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_topdown_all_claims.csv")
  
  #'  ----------------
  ######  Fisher's C  ######
  #'  ----------------
  fishers.C_topdown <- fishers_C(pval = p.val_topdown_all_df$bayes.p, n_iter = nrow(p.val_topdown_all_df))
  print(fishers.C_topdown)
  
  
  
  #'  -------------------------------------------
  #####  Top-down interference model iterations  #####
  #'  -------------------------------------------
  #'  Load all iterations of the JAGS model
  #'  Note: this list of outputs is based on the number of independence claims 
  #'  assessed using the d_Sep_active_regression_topdown_inter.R active regression list. 
  #'  This is not a complete list of all independence claims being tested for this SEM.
  all_results_topdown_inter <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/TopDown_Interference", full.names = TRUE), readRDS)
  
  #'  Rename data bundle
  data_JAGS_bundle <- data_JAGS_bundle_topdown_inter
  #'  Create list of "observed" values of focal response variable, one per d-Sep test               # instances where x was used as y in d-Sep test noted below
  y_list <- list(data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,    #
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,   #
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,      #
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,   # 
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$coy.hat,    #      
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,   #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,  #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$wtd.hat,     # wtd flipped both times
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$wtd.hat,      # wtd flipped all three times
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,    # wtd flipped    
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,    #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,   #      
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,   #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,    #
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat,  #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$wolf.hat,     #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,   #     
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat,    #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,    #     
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$elk.hat,    #
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,    #    
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat,    #  
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat,   #
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$wolf.hat,    #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$lion.hat)                                                           
  #'  Leaves you with 74 d-Sep tests that were t-1 --> t or t --> t
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta.harvest indices. 
  #'  Most will be [,1] but some will be [,2] where the same beta name was used 
  #'  twice in the same regression.
  mod_out <- list()
  for(i in 1:length(all_results_topdown_inter)) {
    mod_out[[i]] <- all_results_topdown_inter[[i]]$fit$sims.list
  }
  #'  Note: indexing is based on the number of independence claims assessed using the
  #'  d_Sep_active_regression_topdown_inter.R active regression list. This is not a complete
  #'  list of all independence claims being tested.
  post_list_topdown_inter <- list(mod_out[[1]]$beta.wtd[,1], mod_out[[2]]$beta.wtd[,1], mod_out[[3]]$beta.wtd[,1],
                                  mod_out[[4]]$beta.wtd[,1], mod_out[[5]]$beta.wtd[,1], mod_out[[6]]$beta.wtd[,1],
                                  mod_out[[7]]$beta.moose[,1], mod_out[[8]]$beta.moose[,1], mod_out[[9]]$beta.moose[,1],
                                  mod_out[[10]]$beta.moose[,1], mod_out[[11]]$beta.moose[,1], mod_out[[12]]$beta.moose[,1],
                                  mod_out[[13]]$beta.elk[,1], mod_out[[14]]$beta.elk[,1], mod_out[[15]]$beta.elk[,1],
                                  mod_out[[16]]$beta.elk[,1], mod_out[[17]]$beta.elk[,1], mod_out[[18]]$beta.elk[,1],
                                  mod_out[[19]]$beta.lion[,1], mod_out[[20]]$beta.lion[,1], mod_out[[21]]$beta.lion[,1],
                                  mod_out[[22]]$beta.lion[,1], mod_out[[23]]$beta.coy[,1], mod_out[[24]]$beta.harvest[,1],        
                                  mod_out[[25]]$beta.bear[,1], mod_out[[26]]$beta.harvest[,1], mod_out[[27]]$beta.harvest[,1],             
                                  mod_out[[28]]$beta.wolf[,1], mod_out[[29]]$beta.wtd[,2], mod_out[[30]]$beta.wtd[,2],             # note the different indexing
                                  mod_out[[31]]$beta.wtd[,2], mod_out[[32]]$beta.wtd[,2], mod_out[[33]]$beta.wtd[,2],              # note the different indexing
                                  mod_out[[34]]$beta.wtd[,2], mod_out[[35]]$beta.coy[,1], mod_out[[36]]$beta.coy[,1],              # note the different indexing
                                  mod_out[[37]]$beta.coy[,1], mod_out[[38]]$beta.coy[,1], mod_out[[39]]$beta.coy[,1],    
                                  mod_out[[40]]$beta.harvest[,1], mod_out[[41]]$beta.harvest[,1], mod_out[[42]]$beta.harvest[,1], 
                                  mod_out[[43]]$beta.harvest[,2], mod_out[[44]]$beta.harvest[,2], mod_out[[45]]$beta.bear[,1],     # note the different indexing
                                  mod_out[[46]]$beta.bear[,1], mod_out[[47]]$beta.bear[,1], mod_out[[48]]$beta.bear[,1],  
                                  mod_out[[49]]$beta.bear[,1], mod_out[[50]]$beta.harvest[,1], mod_out[[51]]$beta.harvest[,1], 
                                  mod_out[[52]]$beta.harvest[,1], mod_out[[53]]$beta.harvest[,2], mod_out[[54]]$beta.harvest[,2], # note the different indexing
                                  mod_out[[55]]$beta.harvest[,1], mod_out[[56]]$beta.harvest[,1], mod_out[[57]]$beta.harvest[,1],       
                                  mod_out[[58]]$beta.harvest[,2], mod_out[[59]]$beta.harvest[,2], mod_out[[60]]$beta.moose[,2],   # note the different indexing
                                  mod_out[[61]]$beta.moose[,2], mod_out[[62]]$beta.moose[,2], mod_out[[63]]$beta.moose[,2],       # note the different indexing
                                  mod_out[[64]]$beta.moose[,2], mod_out[[65]]$beta.elk[,2], mod_out[[66]]$beta.elk[,2],           # note the different indexing
                                  mod_out[[67]]$beta.elk[,2], mod_out[[68]]$beta.elk[,2], mod_out[[69]]$beta.coy[,2],             # note the different indexing
                                  mod_out[[70]]$beta.coy[,2], mod_out[[71]]$beta.coy[,2], mod_out[[72]]$beta.bear[,2],            # note the different indexing
                                  mod_out[[73]]$beta.bear[,2], mod_out[[74]]$beta.wolf[,2])                                       # note the different indexing
  
  #'  Load more iterations of the JAGS model
  #'  Note: this list of outputs is based on the number of independence claims 
  #'  assessed using the d_Sep_active_regression_topdown_inter_tmin1_only.R list. 
  #'  This is not a complete list of all independence claims being tested.
  all_results_topdown_inter_tmin1 <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/TopDown_Interference", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test                 # instances where x was used as y in d-Sep test noted below
  y_list2 <- list(data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$lion.hat,
                  data_JAGS_bundle$coy.hat, data_JAGS_bundle$bearHarv, data_JAGS_bundle$bear.hat,     # harvest flipped
                  data_JAGS_bundle$wtd.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$wolf.hat,      # wtd both flipped
                  data_JAGS_bundle$elk.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat,
                  data_JAGS_bundle$moose.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,  # moose both flipped
                  data_JAGS_bundle$moose.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$coy.hat,    # moose flipped
                  data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$elk.hat,      # elk both flipped
                  data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$coy.hat,      # elk flipped
                  data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat,    # lion both flipped
                  data_JAGS_bundle$lion.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$coy.hat,     # lion and coy flipped
                  data_JAGS_bundle$bear.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$coy.hat,      # coy both flipped
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, 
                  data_JAGS_bundle$bear.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,    # bear both flipped
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wolf.hat)
  #'  Leaves you with 41 d-Sep tests that were t-1 --> t-1
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Generally don't need to worry about indexing here because most of these independence
  #'  claims are really marginal independence claims (not conditional ind. claims) 
  #'  so there is typically only 1 explanatory variable to consider and index will be [1]
  mod_out2 <- list()
  for(i in 1:length(all_results_topdown_inter_tmin1)) {
    mod_out2[[i]] <- all_results_topdown_inter_tmin1[[i]]$fit$sims.list
  }
  post_list_topdown_inter_tmin1 <- list(mod_out2[[1]]$b_x, mod_out2[[2]]$b_x,  mod_out2[[3]]$b_x,
                                        mod_out2[[4]]$b_x,  mod_out2[[5]]$b_x,  mod_out2[[6]]$b_x,
                                        mod_out2[[7]]$b_x,  mod_out2[[8]]$b_x,  mod_out2[[9]]$b_x,
                                        mod_out2[[10]]$b_x,  mod_out2[[11]]$b_x,  mod_out2[[12]]$b_x, 
                                        mod_out2[[13]]$b_x,  mod_out2[[14]]$b_x,  mod_out2[[15]]$b_x, 
                                        mod_out2[[16]]$b_x,  mod_out2[[17]]$b_x,  mod_out2[[18]]$b_x,
                                        mod_out2[[19]]$b_x,  mod_out2[[20]]$b_x,  mod_out2[[21]]$b_x, 
                                        mod_out2[[22]]$b_x,  mod_out2[[23]]$b_x,  mod_out2[[24]]$b_x,
                                        mod_out2[[25]]$b_x,  mod_out2[[26]]$b_x,  mod_out2[[27]]$b_x,
                                        mod_out2[[28]]$b_x,  mod_out2[[29]]$b_x,  mod_out2[[30]]$b_x,
                                        mod_out2[[31]]$b_x,  mod_out2[[32]]$b_x,  mod_out2[[33]]$b_x,
                                        mod_out2[[34]]$b_x,  mod_out2[[35]]$b_x,  mod_out2[[36]]$b_x, 
                                        mod_out2[[37]]$b_x,  mod_out2[[38]]$b_x,  mod_out2[[39]]$b_x, 
                                        mod_out2[[40]]$b_x,  mod_out2[[41]]$b_x)
                                  
  
  #'  Load more iterations of the JAGS model (this time assessing correlation between exogenous variables)
  all_results_topdown_inter_exog <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/TopDown_Interference_exog", full.names = TRUE), readRDS)
  y_list3 <- list(data_JAGS_bundle$wolfHarv, data_JAGS_bundle$lionHarv, data_JAGS_bundle$lionHarv) # y_array
  #'  Leaves you with 3 d-Sep tests that were exog --> exog
  
  mod_out3 <- list()
  for(i in 1:length(all_results_topdown_inter_exog)) {
    mod_out3[[i]] <- all_results_topdown_inter_exog[[i]]$fit$sims.list
  }
  post_list_topdown_inter_exog <- list(mod_out3[[1]]$b_x, mod_out3[[2]]$b_x, mod_out3[[3]]$b_x) # correspond to x_array
  
  
  #'  -------------------------
  ######  ROPE method p-value  ######
  #'  -------------------------
  #'  Calculate p.rope value for each iteration of the d-Sep test
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_topdown_inter_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_topdown_inter, SIMPLIFY = FALSE)
  p.rope_topdown_inter_tmin1_list <- mapply(p_rope_iterations, y_dat = y_list2, post_beta = post_list_topdown_inter_tmin1, SIMPLIFY = FALSE)
  p.rope_topdown_inter_exog_list <- mapply(p_rope_iterations, y_dat = y_list3, post_beta = post_list_topdown_inter_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_topdown_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_topdown_inter_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_topdown_inter_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_topdown_inter_tmin1_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_topdown_inter_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_topdown_inter_exog_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_topdown_inter_df <- stack(p.rope_topdown_inter_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "normal")
  p.rope_topdown_inter_tmin1_df <- stack(p.rope_topdown_inter_tmin1_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "tmin1")
  p.rope_topdown_inter_exog_df <- stack(p.rope_topdown_inter_exog_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "exog")
  
  #'  ----------------------
  ######  Bayesian p-value  ######
  #'  ----------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_topdown_inter_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown_inter, SIMPLIFY = FALSE)
  bayes.p_topdown_inter_tmin1_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown_inter_tmin1, SIMPLIFY = FALSE)
  bayes.p_topdown_inter_exog_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown_inter_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_topdown_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_topdown_inter_list)[i] <- list_name
  }
  for(i in 1:length(bayes.p_topdown_inter_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_topdown_inter_tmin1_list)[i] <- list_name
  }
  for(i in 1:length(bayes.p_topdown_inter_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_topdown_inter_exog_list)[i] <- list_name
  }
  bayes.p_topdown_inter_df <- stack(bayes.p_topdown_inter_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "normal")
  bayes.p_topdown_inter_tmin1_df <- stack(bayes.p_topdown_inter_tmin1_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "tmin1")
  bayes.p_topdown_inter_exog_df <- stack(bayes.p_topdown_inter_exog_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "exog")
  
  #'  Join both d-Sep test p-values and save
  p.val_topdown_inter_df <- full_join(p.rope_topdown_inter_df, bayes.p_topdown_inter_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_topdown_inter_tmin1_df <- full_join(p.rope_topdown_inter_tmin1_df, bayes.p_topdown_inter_tmin1_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_topdown_inter_exog_df <- full_join(p.rope_topdown_inter_exog_df, bayes.p_topdown_inter_exog_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_topdown_inter_all_df <- bind_rows(p.val_topdown_inter_df, p.val_topdown_inter_tmin1_df, p.val_topdown_inter_exog_df)
  
  write_csv(p.val_topdown_inter_all_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_topdown_inter_all_claims.csv")
  
  #'  ----------------
  ######  Fisher's C  ######
  #'  ----------------
  fishers.C_topdown_inter <- fishers_C(pval = p.val_topdown_inter_all_df$bayes.p, n_iter = nrow(p.val_topdown_inter_all_df))
  print(fishers.C_topdown_inter)
  
  
  
  #'  -------------------------------
  #####  Bottom-up model iterations  #####
  #'  -------------------------------
  #'  Load all iterations of the JAGS model
  all_results_bottomup <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/BottomUp", full.names = TRUE), readRDS)
  
  #'  Rename data bundle
  data_JAGS_bundle <- data_JAGS_bundle_bottomup
  #'  Create list of "observed" values of focal response variable, one per d-Sep test              # instances where x was used as y in d-Sep test noted below
  y_list <- list(data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wtd.hat,  #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wtd.hat,  #     
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat,    #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$coy.hat,    #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat,   #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wtd.hat, # moose both flipped        
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,   # moose flipped    
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,  #
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$elk.hat,     # 
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$elk.hat,   #      
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$coy.hat,   #      
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat,    # wtd flipped     
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat,   # coy flipped   
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,   #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,  #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,  #
                 data_JAGS_bundle$lion.hat)                                                        #                                  
  #'  Leaves you with 58 d-Sep tests that were t-1 --> t or t --> t
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta.spp indices. 
  #'  Most will be [,1] but some will be [,2] where the same beta name was used 
  #'  twice in the same regression.
  mod_out <- list()
  for(i in 1:length(all_results_bottomup)) {
    mod_out[[i]] <- all_results_bottomup[[i]]$fit$sims.list
  }
  post_list_bottomup <- list(mod_out[[1]]$beta.wsi[,1], mod_out[[2]]$beta.wsi[,1], mod_out[[3]]$beta.wsi[,1],
                             mod_out[[4]]$beta.wsi[,1], mod_out[[5]]$beta.coy[,1], mod_out[[6]]$beta.coy[,1],
                             mod_out[[7]]$beta.coy[,1], mod_out[[8]]$beta.coy[,1], mod_out[[9]]$beta.coy[,1],
                             mod_out[[10]]$beta.coy[,1], mod_out[[11]]$beta.forest[,1], mod_out[[12]]$beta.forest[,1],
                             mod_out[[13]]$beta.forest[,1], mod_out[[14]]$beta.bear[,1], mod_out[[15]]$beta.bear[,1],
                             mod_out[[16]]$beta.bear[,1], mod_out[[17]]$beta.bear[,1], mod_out[[18]]$beta.bear[,1],
                             mod_out[[19]]$beta.bear[,1], mod_out[[20]]$beta.moose[,1], mod_out[[21]]$beta.moose[,1], 
                             mod_out[[22]]$beta.moose[,1], mod_out[[23]]$beta.moose[,1], mod_out[[24]]$beta.moose[,1], 
                             mod_out[[25]]$beta.wolf[,1], mod_out[[26]]$beta.wtd[,1], mod_out[[27]]$beta.moose[,2],    # note the indexing
                             mod_out[[28]]$beta.moose[,2], mod_out[[29]]$beta.elk[,1], mod_out[[30]]$beta.moose[,2],   # note the indexing
                             mod_out[[31]]$beta.moose[,2], mod_out[[32]]$beta.moose[,2], mod_out[[33]]$beta.moose[,2], # note the indexing
                             mod_out[[34]]$beta.wolf[,1], mod_out[[35]]$beta.wolf[,1], mod_out[[36]]$beta.wolf[,1],       
                             mod_out[[37]]$beta.wolf[,1], mod_out[[38]]$beta.wolf[,1], mod_out[[39]]$beta.wtd[,1],       
                             mod_out[[40]]$beta.wtd[,1], mod_out[[41]]$beta.wtd[,1], mod_out[[42]]$beta.wtd[,2],       # note the indexing
                             mod_out[[43]]$beta.elk[,1], mod_out[[44]]$beta.wtd[,2], mod_out[[45]]$beta.wtd[,2],       # note the indexing
                             mod_out[[46]]$beta.wtd[,2], mod_out[[47]]$beta.wtd[,2], mod_out[[48]]$beta.elk[,1],       # note the indexing
                             mod_out[[49]]$beta.coy[,2], mod_out[[50]]$beta.coy[,2], mod_out[[51]]$beta.coy[,2],       # note the indexing
                             mod_out[[52]]$beta.coy[,2], mod_out[[53]]$beta.elk[,2], mod_out[[54]]$beta.elk[,2],       # note the indexing
                             mod_out[[55]]$beta.elk[,2], mod_out[[56]]$beta.bear[,2], mod_out[[56]]$beta.bear[,2],     # note the indexing
                             mod_out[[57]]$beta.wolf[,2])                                                              # note the indexing                                             
  
  
  #'  Load more iterations of the JAGS model
  #'  Note: this list of outputs is based on the number of independence claims 
  #'  assessed using the d_Sep_active_regression_bottomup_tmin1_only.R list. 
  #'  This is not a complete list of all independence claims being tested.
  all_results_bottomup_tmin1 <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/BottomUp", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test                 # instances where x was used as y in d-Sep test noted below
  y_list2 <- list(data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,    #
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      #
                  data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,    # coy flipped
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      #
                  data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wolf.hat,   #
                  data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$moose.hat,     #
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      # 
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      #
                  data_JAGS_bundle$wtd.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat)       #
  #'  Leaves you with 27 d-Sep tests that were t-1 --> t-1
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Generally don't need to worry about indexing here because most of these independence
  #'  claims are really marginal independence claims (not conditional ind. claims) 
  #'  so there is typically only 1 explanatory variable to consider and index will be [1]
  mod_out2 <- list()
  for(i in 1:length(all_results_bottomup_tmin1)) {
    mod_out2[[i]] <- all_results_bottomup_tmin1[[i]]$fit$sims.list
  }
  post_list_bottomup_tmin1 <- list(mod_out2[[1]]$b_x, mod_out2[[2]]$b_x,  mod_out2[[3]]$b_x,
                                   mod_out2[[4]]$b_x,  mod_out2[[5]]$b_x,  mod_out2[[6]]$b_x,
                                   mod_out2[[7]]$b_x,  mod_out2[[8]]$b_x,  mod_out2[[9]]$b_x,
                                   mod_out2[[10]]$b_x,  mod_out2[[11]]$b_x,  mod_out2[[12]]$b_x, 
                                   mod_out2[[13]]$b_x,  mod_out2[[14]]$b_x,  mod_out2[[15]]$b_x, 
                                   mod_out2[[16]]$b_x,  mod_out2[[17]]$b_x,  mod_out2[[18]]$b_x,
                                   mod_out2[[19]]$b_x,  mod_out2[[20]]$b_x,  mod_out2[[21]]$b_x, 
                                   mod_out2[[22]]$b_x,  mod_out2[[23]]$b_x,  mod_out2[[24]]$b_x,
                                   mod_out2[[25]]$b_x,  mod_out2[[26]]$b_x,  mod_out2[[27]]$b_x)
  
  #'  Load more iterations of the JAGS model (this time assessing correlation between exogenous variables)
  all_results_bottomup_exog <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/BottomUp_exog", full.names = TRUE), readRDS)
  y_list3 <- list(data_JAGS_bundle$forest) # y_array
  #'  Leaves you with 1 d-Sep tests that were exog --> exog
  
  mod_out3 <- list()
  for(i in 1:length(all_results_bottomup_exog)) {
    mod_out3[[i]] <- all_results_bottomup_exog[[i]]$fit$sims.list
  }
  post_list_bottomup_exog <- list(mod_out3[[1]]$b_x) # correspond to x_array
  
  
  #'  -------------------------
  ######  ROPE method p-value  ######
  #'  -------------------------
  #'  Calculate p.rope value for each iteration of the d-Sep test
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_bottomup_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_bottomup, SIMPLIFY = FALSE)
  p.rope_bottomup_tmin1_list <- mapply(p_rope_iterations, y_dat = y_list2, post_beta = post_list_bottomup_tmin1, SIMPLIFY = FALSE)
  p.rope_bottomup_exog_list <- mapply(p_rope_iterations, y_dat = y_list3, post_beta = post_list_bottomup_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_bottomup_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_bottomup_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_tmin1_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_bottomup_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_exog_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_bottomup_df <- stack(p.rope_bottomup_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "normal")
  p.rope_bottomup_tmin1_df <- stack(p.rope_bottomup_tmin1_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "tmin1")
  p.rope_bottomup_exog_df <- stack(p.rope_bottomup_exog_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "exog")
  
  #'  ----------------------
  ######  Bayesian p-value  ######
  #'  ----------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_bottomup_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup, SIMPLIFY = FALSE)
  bayes.p_bottomup_tmin1_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup_tmin1, SIMPLIFY = FALSE)
  bayes.p_bottomup_exog_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_bottomup_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_list)[i] <- list_name
  }
  for(i in 1:length(bayes.p_bottomup_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_tmin1_list)[i] <- list_name
  }
  for(i in 1:length(bayes.p_bottomup_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_exog_list)[i] <- list_name
  }
  bayes.p_bottomup_df <- stack(bayes.p_bottomup_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "normal")
  bayes.p_bottomup_tmin1_df <- stack(bayes.p_bottomup_tmin1_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "tmin1")
  bayes.p_bottomup_exog_df <- stack(bayes.p_bottomup_exog_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "exog")
  
  #'  Join both d-Sep test p-values and save
  p.val_bottomup_df <- full_join(p.rope_bottomup_df, bayes.p_bottomup_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_bottomup_tmin1_df <- full_join(p.rope_bottomup_tmin1_df, bayes.p_bottomup_tmin1_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_bottomup_exog_df <- full_join(p.rope_bottomup_exog_df, bayes.p_bottomup_exog_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_bottomup_all_df <- bind_rows(p.val_bottomup_df, p.val_bottomup_tmin1_df, p.val_bottomup_exog_df)
  
  write_csv(p.val_bottomup_all_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_bottomup_all_claims.csv")
  
  #'  ----------------
  ######  Fisher's C  ######
  #'  ----------------
  fishers.C_bottomup <- fishers_C(pval = p.val_bottomup_all_df$bayes.p, n_iter = nrow(p.val_bottomup_all_df))
  print(fishers.C_bottomup)
  
  
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Load all iterations of the JAGS model
  all_results_bottomup_inter <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/BottomUp_Interference", full.names = TRUE), readRDS)
  
  #'  Rename data bundle
  data_JAGS_bundle <- data_JAGS_bundle_bottomup_inter
  #'  Create list of "observed" values of focal response variable, one per d-Sep test                   # instances where x was used as y in d-Sep test noted below
  y_list <- list(data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,        #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wtd.hat,       #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,        #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$wolf.hat,        #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wtd.hat,       #    
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$wolf.hat,         #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$coy.hat,         #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat,        #
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wtd.hat,      # moose both flipped
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,        # moose flipped 
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,       #      
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$lion.hat,         #
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,        # 
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,          # wtd flipped 
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,       #        
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$bear.hat,         # coy flipped
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat,       #
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$wolf.hat,       #
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$lion.hat)                                  
  #'  Leaves you with 56 d-Sep tests that were possible given the constructs of space and time and our data
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta.spp indices. 
  #'  Most will be [,1] but some will be [,2] where the same beta name was used 
  #'  twice in the same regression.
  mod_out <- list()
  for(i in 1:length(all_results_bottomup_inter)) {
    mod_out[[i]] <- all_results_bottomup_inter[[i]]$fit$sims.list
  }
  post_list_bottomup_inter <- list(mod_out[[1]]$beta.wsi[1], mod_out[[2]]$beta.wsi[1], mod_out[[3]]$beta.wsi[1],
                                   mod_out[[4]]$beta.wsi[1], mod_out[[5]]$beta.coy[,1], mod_out[[6]]$beta.coy[,1],
                                   mod_out[[7]]$beta.coy[,1], mod_out[[8]]$beta.coy[,1], mod_out[[9]]$beta.coy[,1],
                                   mod_out[[10]]$beta.coy[,1], mod_out[[11]]$beta.forest[1], mod_out[[12]]$beta.forest[1],
                                   mod_out[[13]]$beta.forest[1], mod_out[[14]]$beta.bear[,1], mod_out[[15]]$beta.bear[,1],
                                   mod_out[[16]]$beta.bear[,1], mod_out[[17]]$beta.bear[,1], mod_out[[18]]$beta.bear[,1],
                                   mod_out[[19]]$beta.bear[,1], mod_out[[20]]$beta.moose[,1], mod_out[[21]]$beta.moose[,1],
                                   mod_out[[22]]$beta.moose[,1], mod_out[[23]]$beta.moose[,1], mod_out[[24]]$beta.moose[,1], 
                                   mod_out[[25]]$beta.wolf[,1], mod_out[[26]]$beta.wtd[,1], mod_out[[27]]$beta.moose[,2],    # note the different indexing
                                   mod_out[[28]]$beta.moose[,2], mod_out[[29]]$beta.elk[,1], mod_out[[30]]$beta.moose[,2],   # note the different indexing
                                   mod_out[[31]]$beta.moose[,2], mod_out[[32]]$beta.moose[,2], mod_out[[33]]$beta.moose[,2], # note the different indexing
                                   mod_out[[34]]$beta.wolf[,1], mod_out[[35]]$beta.wolf[,1], mod_out[[36]]$beta.wolf[,1],
                                   mod_out[[37]]$beta.wtd[,1], mod_out[[38]]$beta.wtd[,1], mod_out[[39]]$beta.wtd[,1],       
                                   mod_out[[40]]$beta.wtd[,2], mod_out[[41]]$beta.elk[,1], mod_out[[42]]$beta.wtd[,2],       # note the different indexing
                                   mod_out[[43]]$beta.wtd[,2], mod_out[[44]]$beta.wtd[,2], mod_out[[45]]$beta.wtd[,2],       # note the different indexing
                                   mod_out[[46]]$beta.elk[,1], mod_out[[47]]$beta.coy[,2], mod_out[[48]]$beta.coy[,2],       # note the different indexing
                                   mod_out[[49]]$beta.coy[,2], mod_out[[50]]$beta.coy[,2], mod_out[[51]]$beta.elk[,2],       # note the different indexing
                                   mod_out[[52]]$beta.elk[,2], mod_out[[53]]$beta.elk[,2], mod_out[[54]]$beta.bear[,2],      # note the different indexing
                                   mod_out[[55]]$beta.bear[,2], mod_out[[56]]$beta.wolf[,2])                                 # note the different indexing
  
  #'  Load more iterations of the JAGS model
  #'  Note: this list of outputs is based on the number of independence claims 
  #'  assessed using the d_Sep_active_regression_bottomup_inter_tmin1_only.R list. 
  #'  This is not a complete list of all independence claims being tested.
  all_results_bottomup_inter_tmin1 <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/BottomUp_Interference", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test                 # instances where x was used as y in d-Sep test noted below
  y_list2 <- list(data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,    #
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      #
                  data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat,    # coy flipped
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      #
                  data_JAGS_bundle$bear.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$wolf.hat,   #
                  data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$moose.hat,     #
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      # 
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat,      #
                  data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$elk.hat)       #
  #'  Leaves you with 27 d-Sep tests that were t-1 --> t-1
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Generally don't need to worry about indexing here because most of these independence
  #'  claims are really marginal independence claims (not conditional ind. claims) 
  #'  so there is typically only 1 explanatory variable to consider and index will be [1]
  mod_out2 <- list()
  for(i in 1:length(all_results_bottomup_inter_tmin1)) {
    mod_out2[[i]] <- all_results_bottomup_inter_tmin1[[i]]$fit$sims.list
  }
  post_list_bottomup_inter_tmin1 <- list(mod_out2[[1]]$b_x, mod_out2[[2]]$b_x,  mod_out2[[3]]$b_x,
                                         mod_out2[[4]]$b_x,  mod_out2[[5]]$b_x,  mod_out2[[6]]$b_x,
                                         mod_out2[[7]]$b_x,  mod_out2[[8]]$b_x,  mod_out2[[9]]$b_x,
                                         mod_out2[[10]]$b_x,  mod_out2[[11]]$b_x,  mod_out2[[12]]$b_x, 
                                         mod_out2[[13]]$b_x,  mod_out2[[14]]$b_x,  mod_out2[[15]]$b_x, 
                                         mod_out2[[16]]$b_x,  mod_out2[[17]]$b_x,  mod_out2[[18]]$b_x,
                                         mod_out2[[19]]$b_x,  mod_out2[[20]]$b_x,  mod_out2[[21]]$b_x, 
                                         mod_out2[[22]]$b_x,  mod_out2[[23]]$b_x,  mod_out2[[24]]$b_x,
                                         mod_out2[[25]]$b_x,  mod_out2[[26]]$b_x,  mod_out2[[27]]$b_x)
  
  #'  Load more iterations of the JAGS model (this time assessing correlation between exogenous variables)
  all_results_bottomup_inter_exog <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/BottomUp_Interference_exog", full.names = TRUE), readRDS)
  y_list3 <- list(data_JAGS_bundle$forest) # y_array
  #'  Leaves you with 1 d-Sep tests that were exog --> exog
  
  mod_out3 <- list()
  for(i in 1:length(all_results_bottomup_inter_exog)) {
    mod_out3[[i]] <- all_results_bottomup_inter_exog[[i]]$fit$sims.list
  }
  post_list_bottomup_inter_exog <- list(mod_out3[[1]]$b_x) # correspond to x_array
  
  #'  -------------------------
  ######  ROPE method p-value  ######
  #'  -------------------------
  #'  Calculate p.rope value for each iteration of the d-Sep test
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_bottomup_inter_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_bottomup_inter, SIMPLIFY = FALSE)
  p.rope_bottomup_inter_tmin1_list <- mapply(p_rope_iterations, y_dat = y_list2, post_beta = post_list_bottomup_inter_tmin1, SIMPLIFY = FALSE)
  p.rope_bottomup_inter_exog_list <- mapply(p_rope_iterations, y_dat = y_list3, post_beta = post_list_bottomup_inter_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_bottomup_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_inter_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_bottomup_inter_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_inter_tmin1_list)[i] <- list_name
  }
  for(i in 1:length(p.rope_bottomup_inter_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_inter_exog_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_bottomup_inter_df <- stack(p.rope_bottomup_inter_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "normal")
  p.rope_bottomup_inter_tmin1_df <- stack(p.rope_bottomup_inter_tmin1_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "tmin1")
  p.rope_bottomup_inter_exog_df <- stack(p.rope_bottomup_inter_exog_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4),
              basicset = "exog")
  
  #'  ----------------------
  ######  Bayesian p-value  ######
  #'  ----------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_bottomup_inter_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup_inter, SIMPLIFY = FALSE)
  bayes.p_bottomup_inter_tmin1_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup_inter_tmin1, SIMPLIFY = FALSE)
  bayes.p_bottomup_inter_exog_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup_inter_exog, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_bottomup_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_inter_list)[i] <- list_name
  }
  bayes.p_bottomup_inter_df <- stack(bayes.p_bottomup_inter_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "normal")
  
  for(i in 1:length(bayes.p_bottomup_inter_tmin1_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_inter_tmin1_list)[i] <- list_name
  }
  bayes.p_bottomup_inter_tmin1_df <- stack(bayes.p_bottomup_inter_tmin1_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "tmin1")
  
  for(i in 1:length(bayes.p_bottomup_inter_exog_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_inter_exog_list)[i] <- list_name
  }
  bayes.p_bottomup_inter_exog_df <- stack(bayes.p_bottomup_inter_exog_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4),
              basicset = "exog")
  
  #'  Join both d-Sep test p-values and save
  p.val_bottomup_inter_df <- full_join(p.rope_bottomup_inter_df, bayes.p_bottomup_inter_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_bottomup_inter_tmin1_df <- full_join(p.rope_bottomup_inter_tmin1_df, bayes.p_bottomup_inter_tmin1_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_bottomup_inter_exog_df <- full_join(p.rope_bottomup_inter_exog_df, bayes.p_bottomup_inter_exog_df, by = c("iteration", "basicset")) %>% relocate("basicset", .after = "bayes.p")
  p.val_bottomup_inter_all_df <- bind_rows(p.val_bottomup_inter_df, p.val_bottomup_inter_tmin1_df, p.val_bottomup_inter_exog_df)
  
  write_csv(p.val_bottomup_inter_all_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_bottomup_inter_all_claims_updated.csv")
  
  #'  ----------------
  ######  Fisher's C  ######
  #'  ----------------
  fishers.C_bottomup_inter <- fishers_C(pval = p.val_bottomup_inter_all_df$bayes.p, n_iter = nrow(p.val_bottomup_inter_all_df))
  print(fishers.C_bottomup_inter)
  
  
  
  #'  --------------------------------------
  ####  Final SEMs following d-Sep updates  ####
  #'  --------------------------------------
  #'  --------------------------------
  #####  Top-down exploitative model  #####
  #'  --------------------------------
  #'  Generate DAG
  dag_topdown_final <- DAG(lion.t ~ lionHarv.tmin1 + elk.tmin1 + wtd.tmin1 + wtd.t, 
                           wolf.t ~ wolf.tmin1 + wolfHarv.tmin1 + moose.tmin1 + elk.tmin1 + wtd.tmin1 + wtd.t,
                           bear.t ~ bear.tmin1 + bearHarv.tmin1 + wtd.t,
                           coy.t ~ coy.tmin1 + wtd.tmin1,
                           elk.t ~ elk.tmin1 + wolf.tmin1 + lion.tmin1 + elkHarv.tmin1 + bear.tmin1,
                           moose.t ~ moose.tmin1 + wolf.tmin1,
                           wtd.t ~ wtd.tmin1 + lion.tmin1 + deerHarv.tmin1 + wolf.tmin1 + coy.tmin1)
                     
  
  #'  Generate basic set
  bs_topdown_final <- basic_set(dag_topdown_final)
  
  #'  ----------------------
  ######  d-Sep iterations  ######
  #'  ----------------------
  #'  Fit independence claims for variables where t-1 --> t or t --> t
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.latent
    list(covs = c("lionHarv", "elk.latent", "wtd.latent", "wtd.latent"), spp = c(".harvest", ".elk", ".wtd", ".wtd"), indices = as.integer(c(1,1,1,2)), lags = c("y-1","y-1","y-1","y")),
    #'  Regression 2: wolf.latent
    list(covs = c("wolf.latent", "wolfHarv", "moose.latent", "elk.latent", "wtd.latent", "wtd.latent"), spp = c(".wolf", ".harvest", ".moose", ".elk", ".wtd", ".wtd"), indices = as.integer(c(1,1,1,1,1,2)), lags = c("y-1","y-1","y-1","y-1","y-1","y")),
    #'  Regression 3: bear.latent
    list(covs = c("bear.latent", "bearHarv", "wtd.latent"), spp = c(".bear", ".harvest", ".wtd"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y")),
    #'  Regression 4: coy.latent
    list(covs = c("coy.latent", "wtd.latent"), spp = c(".coy", ".wtd"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 5: elk.latent
    list(covs = c("elk.latent", "wolf.latent", "lion.latent", "elkHarv", "bear.latent"), spp = c(".elk", ".wolf", ".lion", ".harvest", ".bear"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1")),
    #'  Regression 6: moose.latent
    list(covs = c("moose.latent", "wolf.latent"), spp = c(".moose", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 7: wtd.latent
    list(covs = c("wtd.latent", "lion.latent", "deerHarv", "wolf.latent", "coy.latent"), spp = c(".wtd", ".lion", ".harvest", ".wolf", ".coy"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1"))
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_final.R")
  
  #'  Bundle data and draw inits using functions in in Format_RNmodel_Posteriors_for_SEM.R
  data_JAGS_bundle_topdown_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                               dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                               covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                               covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                               nwolf = 5, nlion = 3, nbear = 3, ncoy = 3, nelk = 4, 
                                               nmoose = 3, nwtd = 9, nharv = 6, nfor = 0, nwsi = 0)
                                         
  num.chains <- 3
  initsList_topdown_final <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_topdown_final[[i]] <- generate_inits(nwolf = 5, nlion = 3, nbear = 3, ncoy = 3, nelk = 4, nmoose = 3, 
                                             nwtd = 9, nharv = 6, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  #'  Fit and save model iterations
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_final),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown_final, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_topdown_final, listInits = initsList_topdown_final, model_name = "TopDown_Exploitative_Final"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_tmin1_only_final.R")
  
  ### MAKE SURE SEM_TOPDOWN_FINAL IS IN WORKING DIRECTORY  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_tmin1_only_final),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_topdown_tmin1_only_final, 
                              og_fit = SEM_topdown_final, nSites = 23, nYear = 4, model_name = "TopDown_Exploitative_Final",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_exog_only_final.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_exog_only_final),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_topdown_exog_only_final, 
                                    model_name = "TopDown_Exploitative_Final", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  
  
  
  
  
  
  #'  --------------------------------
  #####  Top-down interference model  #####
  #'  --------------------------------
  dag_topdown_inter_final <- DAG(lion.t ~ wolf.tmin1 + lionHarv.tmin1, 
                                 wolf.t ~ wolf.tmin1 + wolfHarv.tmin1 + bear.t,
                                 bear.t ~ bear.tmin1 + bearHarv.tmin1 + wolf.tmin1,
                                 coy.t ~ coy.tmin1 + wolf.tmin1 + lion.tmin1,
                                 elk.t ~ elk.tmin1 + wolf.tmin1 + lion.tmin1,
                                 moose.t ~ moose.tmin1 + wolf.tmin1,
                                 wtd.t ~ wtd.tmin1 + lion.tmin1 + coy.tmin1 + bear.tmin1)
                           
  
  bs_topdown_inter_final <- basic_set(dag_topdown_inter_final)
  
  #'  ----------------------
  ######  d-Sep iterations  ######
  #'  ----------------------
  #'  Fit independence claims for variables where t-1 --> t or t --> t
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.latent
    list(covs = c("wolf.latent", "lionHarv"), spp = c(".wolf", ".harvest"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")), 
    #'  Regression 2: wolf.latent
    list(covs = c("wolf.latent", "wolfHarv", "bear.latent"), spp = c(".wolf", ".harvest", ".bear"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y")),
    #'  Regression 3: bear.latent
    list(covs = c("bear.latent", "bearHarv", "wolf.latent"), spp = c(".bear", ".harvest", ".wolf"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 4: coy.latent
    list(covs = c("coy.latent", "wolf.latent", "lion.latent"), spp = c(".coy", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 5: elk.latent
    list(covs = c("elk.latent", "wolf.latent", "lion.latent"), spp = c(".elk", ".wolf", ".lion"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 6: moose.latent
    list(covs = c("moose.latent", "wolf.latent"), spp = c(".moose", ".wolf"), indices = as.integer(c(1,1)), lags = c("y-1","y-1")),
    #'  Regression 7: wtd.latent
    list(covs = c("wtd.latent", "lion.latent", "coy.latent", "bear.latent"), spp = c(".wtd", ".lion", ".coy", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1"))
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_inter_final.R")
  
  data_JAGS_bundle_topdown_inter_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                     dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                                     covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                     covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                     nwolf = 7, nlion = 4, nbear = 4, ncoy = 3, nelk = 2, 
                                                     nmoose = 2, nwtd = 2, nharv = 4, nfor = 0, nwsi = 0)
                                               
  num.chains <- 3
  initsList_topdown_inter_final <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_topdown_inter_final[[i]] <- generate_inits(nwolf = 7, nlion = 4, nbear = 4, ncoy = 3, nelk = 2, nmoose = 2, 
                                                   nwtd = 2, nharv = 4, nfor = 0, nwsi = 0, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  #'  Fit and save model iterations
  start.time = Sys.time()
  saved_paths <- future_lapply(
    #'  Apply across every element in list of active regressions
    seq_along(dSep_iterations_topdown_inter_final),
    #'  Call run_dSep_iterations function using specified active regression list, model template, and data/inits prepared for JAGS
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown_inter_final, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_topdown_inter_final, listInits = initsList_topdown_inter_final, model_name = "TopDown_Interference_Final"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_inter_tmin1_only_final.R")
  
  ### MAKE SURE SEM_TOPDOWN_INTER_FINAL IS IN GLOBAL ENVI and spp.latent params were monitored  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_inter_tmin1_only_final),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_topdown_inter_tmin1_only_final, 
                              og_fit = SEM_topdown_inter_final, nSites = 23, nYear = 4, model_name = "TopDown_Interference_Final",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_topdown_inter_exog_only_final.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown_inter_exog_only_final),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_topdown_inter_exog_only_final, 
                                    model_name = "TopDown_Interference_exog_final", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  ---------------------------------
  #####  Bottom-up exploitative model  #####
  #'  ---------------------------------
  dag_bottomup_final <- DAG(lion.t ~ elk.tmin1 + elk.t + wtd.tmin1 + wtd.t, 
                            wolf.t ~ wolf.tmin1 + elk.tmin1 + elk.t + moose.tmin1 + moose.t + wtd.t,
                            bear.t ~ bear.tmin1 + elk.tmin1 + forest.tmin1 + wtd.t,
                            coy.t ~ coy.tmin1 + wtd.tmin1 + wtd.t,
                            elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1 + bear.tmin1,
                            moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                            wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1 + coy.tmin1 + bear.tmin1)
                      
  bs_bottomup_final <- basic_set(dag_bottomup_final)
  
  #'  ----------------------
  ######  d-Sep iterations  ######
  #'  ----------------------
  #'  Fit independence claims for variables where t-1 --> t or t --> t
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.t
    list(covs = c("elk.latent", "wtd.latent", "wtd.latent"), spp = c(".elk", ".wtd", ".wtd"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y")),
    #'  Regression 2: wolf.t
    list(covs = c("wolf.latent", "elk.latent", "elk.latent", "moose.latent", "moose.latent", "wtd.latent"), spp = c(".wolf", ".elk", ".elk", ".moose", ".moose", ".wtd"), indices = as.integer(c(1,1,2,1,2,1)), lags = c("y-1","y-1","y","y-1","y","y")),
    #'  Regression 3: bear.t
    list(covs = c("bear.latent", "elk.latent", "forest", "wtd.latent"), spp = c(".bear", ".elk", ".forest", ".wtd"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y")),
    #'  Regression 4: coy.t
    list(covs = c("coy.latent", "wtd.latent", "wtd.latent"), spp = c(".coy", ".wtd", ".wtd"), indices = as.integer(c(1,1,2)), lags = c("y-1","y-1","y")),
    #'  Regression 5: elk.t
    list(covs = c("elk.latent", "forest", "wsi", "bear.latent"), spp = c(".elk", ".forest", ".wsi", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    #'  Regression 6: moose.t
    list(covs = c("moose.latent", "forest", "wsi"), spp = c(".moose", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 7: wtd.t
    list(covs = c("wtd.latent", "forest", "wsi", "coy.latent", "bear.latent"), spp = c(".wtd", ".forest", ".wsi", ".coy", ".bear"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y-1"))
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_final.R")
  
  #'  Bundle data and draw inits using functions in in Format_RNmodel_Posteriors_for_SEM.R
  data_JAGS_bundle_bottomup_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
                                                covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                nwolf = 2, nlion = 1, nbear = 2, ncoy = 2, nelk = 5, 
                                                nmoose = 3, nwtd = 4, nharv = 0, nfor = 5, nwsi = 5)
                            
  num.chains <- 3
  initsList_bottomup_final <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup_final[[i]] <- generate_inits(nwolf = 2, nlion = 1, nbear = 2, ncoy = 2, nelk = 5, nmoose = 3, 
                                              nwtd = 4, nharv = 0, nfor = 5, nwsi = 5, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_final),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup_final, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_bottomup_final, listInits = initsList_bottomup_final, model_name = "BottomUp_Exploitative_Final"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_tmin1_only_final.R")
  
  ### MAKE SURE SEM_BOTTOMUP_FINAL IS IN GLOBAL ENVIRO and spp.latent params were monitored  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_tmin1_only_final),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_bottomup_tmin1_only_final, 
                              og_fit = SEM_bottomup_final, nSites = 23, nYear = 4, model_name = "BottomUp_Exploitative_Final",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_exog_only_final.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_exog_only_final),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_bottomup_exog_only_final, 
                                    model_name = "BottomUp_Exploitative_exog_final", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  ---------------------------------
  #####  Bottom-up interference model  #####
  #'  ---------------------------------
  dag_bottomup_inter_final <- DAG(lion.t ~ elk.tmin1 + elk.t + wtd.tmin1 + wtd.t + bear.tmin1, 
                                  wolf.t ~ wolf.tmin1 + elk.tmin1 + elk.t + moose.tmin1 + moose.t + bear.t,
                                  bear.t ~ bear.tmin1 + elk.tmin1 + forest.tmin1 + wolf.tmin1 + wtd.t,
                                  coy.t ~ coy.tmin1 + wtd.tmin1 + wolf.tmin1 + bear.tmin1,
                                  elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1,
                                  moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                                  wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1)
                            
  
  bs_bottomup_inter_final <- basic_set(dag_bottomup_inter_final)
  
  #'  ----------------------
  ######  d-Sep iterations  ######
  #'  ----------------------
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.latent
    list(covs = c("elk.latent", "elk.latent", "wtd.latent", "wtd.latent", "bear.latent"), spp = c(".elk", ".elk", ".wtd", ".wtd", ".bear"), indices = as.integer(c(1,2,1,2,1)), lags = c("y-1","y","y-1","y","y-1")),
    #'  Regression 2: wolf.latent
    list(covs = c("wolf.latent", "elk.latent", "elk.latent", "moose.latent", "moose.latent", "bear.latent"), spp = c(".wolf", ".elk", ".elk", ".moose", ".moose", ".bear"), indices = as.integer(c(1,1,2,1,2,1)), lags = c("y-1","y-1","y","y-1","y","y")),
    #'  Regression 3: bear.latent
    list(covs = c("bear.latent", "elk.latent", "forest", "wolf.latent", "wtd.latent"), spp = c(".bear", ".elk", ".forest", ".wolf", ".wtd"), indices = as.integer(c(1,1,1,1,1)), lags = c("y-1","y-1","y-1","y-1","y")),
    #'  Regression 4: coy.latent
    list(covs = c("coy.latent", "wtd.latent", "wolf.latent", "bear.latent"), spp = c(".coy", ".wtd", ".wolf", ".bear"), indices = as.integer(c(1,1,1,1)), lags = c("y-1","y-1","y-1","y-1")),
    #'  Regression 5: elk.latent
    list(covs = c("elk.latent", "forest", "wsi"), spp = c(".elk", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 6: moose.latent
    list(covs = c("moose.latent", "forest", "wsi"), spp = c(".moose", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1")),
    #'  Regression 7: wtd.latent
    list(covs = c("wtd.latent", "forest", "wsi"), spp = c(".wtd", ".forest", ".wsi"), indices = as.integer(c(1,1,1)), lags = c("y-1","y-1","y-1"))
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_inter_updated_final.R")
  
  #'  Bundle data and draw inits using functions in in Format_RNmodel_Posteriors_for_SEM.R
  data_JAGS_bundle_bottomup_inter_final <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                      dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s,
                                                      covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                      covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                      nwolf = 4, nlion = 1, nbear = 5, ncoy = 2, nelk = 7, 
                                                      nmoose = 4, nwtd = 6, nharv = 0, nfor = 5, nwsi = 4)
                                                
  num.chains <- 3
  initsList_bottomup_inter_final <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup_inter_final[[i]] <- generate_inits(nwolf = 4, nlion = 1, nbear = 5, ncoy = 2, nelk = 7, nmoose = 4, 
                                                    nwtd = 6, nharv = 0, nfor = 5, nwsi = 4, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_inter_final), 
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup_inter_final, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_bottomup_inter_final, listInits = initsList_bottomup_inter_final, model_name = "BottomUp_Interference_final"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  Source second d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_inter_updated_tmin1_only_final.R")
  
  ### MAKE SURE SEM_BOTTOMUP_INTER_FINAL IS IN GLOBAL ENVIRO and spp.latent params were monitored  ###
  
  #'  Fit independence claims for variables where t-1 --> t-1 
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_inter_tmin1_only),
    function(i) fit_aux_claim(i, iterations = dSep_iterations_bottomup_inter_tmin1_only, 
                              og_fit = SEM_bottomup_inter_final, nSites = 23, nYear = 4, model_name = "BottomUp_Interference_Final",
                              n.chains = nc, n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  #'  Source third d-Sep custom regressions for iterative d-separation tests -
  #'  this time to simply test correlation between exogenous variables flagged 
  #'  in the basic set
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_inter_updated_exog_only.R")
  #'  Fit independence claims for pairs of exogenous variables
  start.time = Sys.time()
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_inter_exog_only_final),
    function(i) fit_covariate_claim(i, iterations = dSep_iterations_bottomup_inter_exog_only_final, 
                                    model_name = "BottomUp_Interference_exog_final", n.chains = nc, 
                                    n.adapt = na, n.burnin = nb, n.iter = ni, n.thin = nt),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
