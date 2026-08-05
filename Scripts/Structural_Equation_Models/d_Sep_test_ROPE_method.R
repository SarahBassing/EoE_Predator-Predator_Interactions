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
  
  #'  Load libries
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
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  #'  Source functions to setup data for JAGS
  source("./Scripts/Structural_Equation_Models/JAGS_setup_functions.R")
  
  #'  Bundle data for JAGS
  data_JAGS_bundle <- bundle_dat(post_summaries, covs = covs_ztransformed, 
                                 nwolf = 2, nlion = 2, nbear = 2, ncoy = 2, nelk = 2, 
                                 nmoose = 2, nwtd = 2, nharv = 2, nfor = 1, nwsi = 1)
  
  #'  Draw initial values for each chain in JAGS
  num.chains <- 3
  #'  Create empty list
  initsList <- vector('list', num.chains)
  #'  Setting seed for reproducibility
  set.seed(9714)
  #'  Loop through generate_inits function 3 times (1 for each chain) 
  for(i in 1:num.chains){
    initsList[[i]] <- generate_inits(nwolf = 2, nlion = 2, nbear = 2, ncoy = 2, nelk = 2, 
                                     nmoose = 2, nwtd = 2, nharv = 2, nfor = 1, nwsi = 1)
  }
  
  #'  Parameters to be monitored
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk", 
              "beta.moose", "beta.wtd", "beta.harvest", "beta.wsi","beta.forest", "beta.road", "beta.public",
              "sigma.spp", "sigma.spp.tmin1", "sigma.cluster", "cluster.randeff")  
  
  #'  MCMC settings
  nc <- 3
  ni <- 100000
  nb <- 50000
  nt <- 10
  na <- 5000
  
  #'  -------------------------------
  ####  ROPE method for d-Sep tests  ####
  #'  -------------------------------
  #'  Define function to establish a region of practical equivalence (ROPE) around
  #'  the null value. This expresses a range of parameter values considered equivalent
  #'  to the null value (0). This range is -pct*sd(y) to pct*sd(y)
  #'  If 95% CRI falls entirely within this range, 95% of the posterior is practically
  #'  equivalent to the null value. Taking the mean calculates the proportion
  #'  of the posterior that falls inside the ROPE. In other words, what is the 
  #'  probability that the effect is small enough to be ignored? Larger p.rope 
  #'  values indicate a large percentage of the posterior is "practically" 0; 
  #'  these variables are conditionally independent. Small p.rope values (i.e., < 0.05) 
  #'  indicate a small percentage of the posterior falls inside the ROPE and thus  
  #'  it is not practically equivalent to 0 (i.e., it is "practically" significant).
  #'  These variables are NOT conditionally independent and this relationship is 
  #'  missing in the model. Another way to think about it is, by defining the ROPE 
  #'  as a small percentage of the standard deviation of y, we have defined a threshold 
  #'  of practical significance. If the variable of interest only moves the needle 
  #'  by a 10% of a standard deviation, we can consider its effect to be noise.
  #'  
  #'  ---> p.rope is the probability that the parameter is actually trivial given 
  #'  ---> the data, where small values indicate the prob of being trivial is low. 
  #'  
  #'  y represents the SD of the observed response variable
  #'  post represents the effect of the coefficient of interest - if the parameter 
  #'  effect is "meaningful" it should be larger than the typical "noise" or trivial
  #'  variation associated with the response variable
  p.rope <- function(pct = 0.1, y = NULL, post = NULL){
    #'  Average number of posterior draws that fall above lower ROPE value AND below upper ROPE value
    mean(-pct*sd(y) < post & pct*sd(y) > post)
    
  }
  
  #'  Function to generate and simplify the basic set 
  basic_set <- function(dag) {
    #'  Generate the basic set 
    basicset <- basiSet(dag)
    print(length(basicset))
    #'  Return basic set with conditional independence claims involving 3+ variables
    basicset_3plus <- sapply(basicset, length) >= 3
    basicset_skinny <- basicset[basicset_3plus]
    print(length(basicset_skinny))
    View(basicset_skinny)
    return(basicset_skinny)
  }
  
  #'  -------------------------------------
  #####  Top-down model (reduced version)  #####
  #'  -------------------------------------
  #'  Generate DAG
  dag_topdown <- DAG(lion.t ~ lion.tmin1 + lionHarv.tmin1,  
                     wolf.t ~ wolf.tmin1 + wolfHarv.tmin1,  
                     bear.t ~ bear.tmin1 + bearHarv.tmin1,  
                     coy.t ~ coy.tmin1,  
                     elk.t ~ elk.tmin1 + wolf.tmin1 + lion.tmin1 + elkHarv.tmin1,  
                     moose.t ~ moose.tmin1 + wolf.tmin1,  
                     wtd.t ~ wtd.tmin1 + lion.tmin1 + deerHarv.tmin1)  
  
  #'  Generate basic set
  bs_topdown <- basic_set(dag_topdown)
  
  #'  --------------------------------------------------
  #####  Top-down interference model (reduced version)  #####
  #'  --------------------------------------------------
  dag_topdown_inter <- DAG(lion.t ~ lion.tmin1 + wolf.tmin1 + lionHarv.tmin1,  
                           wolf.t ~ wolf.tmin1 + wolfHarv.tmin1, 
                           bear.t ~ bear.tmin1 + wolf.tmin1 + bearHarv.tmin1,  
                           coy.t ~ coy.tmin1 + wolf.tmin1 + lion.tmin1,  
                           elk.t ~ elk.tmin1 + wolf.tmin1 + lion.tmin1,  
                           moose.t ~ moose.tmin1 + wolf.tmin1,  
                           wtd.t ~ wtd.tmin1 + lion.tmin1) 
  
  bs_topdown_inter <- basic_set(dag_topdown_inter)
  
  #'  --------------------------------------
  #####  Bottom-up model (reduced version)  ####
  #'  --------------------------------------
  dag_bottomup <- DAG(lion.t ~ lion.tmin1 + elk.tmin1 + wtd.tmin1,
                      wolf.t ~ wolf.tmin1 + elk.tmin1 + moose.tmin1,
                      bear.t ~ bear.tmin1 + elk.tmin1 + forest.tmin1,
                      coy.t ~ coy.tmin1 + wtd.tmin1,
                      elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1,
                      moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                      wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1)
  bs_bottomup <- basic_set(dag_bottomup)
  
  #'  ---------------------------------------------------
  #####  Bottom-up interference model (reduced version)  ####
  #'  ---------------------------------------------------
  dag_bottomup_inter <- DAG(lion.t ~ lion.tmin1 + wtd.tmin1 + wolf.tmin1, 
                            wolf.t ~ wolf.tmin1 + elk.tmin1 + moose.tmin1,
                            bear.t ~ bear.tmin1 + forest.tmin1 + wolf.tmin1, 
                            coy.t ~ coy.tmin1 + wtd.tmin1 + wolf.tmin1,
                            elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1,
                            moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                            wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1)
  
  bs_bottomup_inter <- basic_set(dag_bottomup_inter)
  
  #'  ---------------------------------------------
  ####  Iterate JAGS model for d-separation tests  ####
  #'  ---------------------------------------------
  #'  Source JAGS template
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_dsep_template.R")
  
  #'  Function to build custom regressions for each iteration
  build_individual_submodels <- function(reg_num, covariates = NULL, spp = NULL, indices = NULL, timestep = "") {  
    #'  Ensures only valid time step suffixes are provided
    if(!timestep %in% c("", ".tmin1")) {
      stop('timestep must be "" or ".tmin1')
    }
    
    #'  Create intercept and slope parameter names based on time step
    beta0_array <- paste0("beta.int", timestep) # beta.int or beta.int.tmin1
    beta_array <- paste0("beta", timestep)      # beta or beta.tmin1
    
    #'  If-Else statement for time t regressions
    #'  If covariate is null or 0 (i.e., intercept only regressions)
    if(is.null(covariates) || length(covariates) == 0) {
      #'  Build intercept only regressions for time t and t-1
      sprintf("%s[%d]", beta0_array, reg_num)
    } else {
      #'  Build one string per "beta * covariate" by filling placeholders with 
      #'  specified character strings or integers (sprintf is vectorized so does 
      #'  this in order that strings/integers are provided)
      #'  %s is placeholder for character strings; %d is placeholder for integers
      terms <- paste0(sprintf("%s%s[%d] * %s[i]", beta_array, spp, indices, covariates),
                      #'  And link regression terms into one single string
                      collapse = " + ")
      #'  Build full regression with intercept plus regression terms defined above
      sprintf("%s[%d] + %s", beta0_array, reg_num, terms)
    }
  }
  
  #'  Function to assemble full model string for a single iteration
  #'  Requires having sourced the JAGS template model
  #'  mu_lines[(r * 2) - 1] indexing maps each species regression (1:7) onto two 
  #'  consecutive slots in the mu_lines vector (for t and t-1 versions) - ensures 
  #'  regressions are in correct order and match template's 14 %s placeholders
  build_model_string <- function(iter_config, template) {
    mu_lines <- character(14)  # 7 regressions for time t and t-1
    for(r in 1:7) {
      if (r == iter_config$dSep_test) {
        #'  Grab covariate, species, and index numbers of terms included in d-sep test
        covs <- iter_config$covariates
        spp <- iter_config$spp
        indices <- iter_config$indices
        #'  time t: create custom regression for d-sep test
        mu_lines[(r * 2) - 1] <- build_individual_submodels(r, covs, spp, indices, timestep = "")
      } else {
        #'  time t: intercept-only for all other regressions
        mu_lines[(r * 2) - 1] <- build_individual_submodels(r, NULL, timestep = "") 
      }
      
      #'  time t-1: always intercept-only for all regressions and iterations
      mu_lines[r * 2] <- build_individual_submodels(r, NULL, timestep = ".tmin1")
    }
    
    do.call(sprintf, c(list(template), as.list(mu_lines)))
  }
  
  #'  Function to call JAGS and run a single iteration of the model
  run_dSep_iterations <- function(i, iterations, template, data_bundle, listInits, model_name) {  
    
    iter_config <- iterations[[i]]
    
    #'  Call function to build the full model string with custom regressions for d-Sep test
    model_string <- build_model_string(iter_config, template)
    
    #'  Create temporary directory to save all iterations of the template
    temp_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/temp_models", model_name)
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Create name new name for each model where %03d is a placeholder for the 
    #'  iteration number, padded by up to three 0's (e.g., temp_model_001, temp_model_012)
    temp_file <- file.path(temp_dir, sprintf("temp_model_%03d.txt", i))
    writeLines(model_string, temp_file)
    
    #'  Snag bundled data
    data_i <- data_bundle
    
    #'  Fit model in JAGS
    SEM_dSep <- jagsUI::jags(data = data_i, inits = listInits, params, model.file = temp_file, 
                             n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                             n.burnin = nb, parallel = FALSE, verbose = FALSE)
                             # chain-level parallelism is off - parallelizing across
                             # iterations instead
    
    #'  Flag convergence issues
    max_rhat <- suppressWarnings(max(unlist(SEM_dSep$Rhat), na.rm = TRUE))
    if(is.finite(max_rhat) && max_rhat > 1.1) {
      warning(sprintf("Iteration %d: max Rhat = %.3f -- check convergence", i, max_rhat))
    }
    
    #'  Temporary directory to save JAGS outputs
    out_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/Results", model_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Save JAGS output for each iteration
    jags_out <- file.path(out_dir, sprintf("iter_%03d.rds", i))
    saveRDS(list(fit = SEM_dSep, config = iter_config, max_rhat = max_rhat), jags_out)
    
    jags_out
  }
  
  #'  Run all iterations in parallel
  plan(multisession, workers = parallel::detectCores() - 1)
  
  #'  Source dSep_interactions: list of custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_topdown_inter_reduced.R")
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_topdown_reduced.R")
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_bottomup_reduced.R")
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_bottomup_inter_reduced.R")
  
  #'  -------------------------------------------
  #####  Top-down interference model iterations  #####
  #'  -------------------------------------------
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    #'  Apply across every element in list of active regressions
    seq_along(dSep_iterations_topdown_int),
    #'  Call run_dSep_iterations function using specified active regression list, model template, and data/inits prepared for JAGS
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown_int, template = model_template, 
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "TopDown_Interference"),
    future.seed = TRUE
  )
  
  #'  ------------------------------
  #####  Top-down model iterations  #####
  #'  ------------------------------
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown, template = model_template, 
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "TopDown"),
    future.seed = TRUE
  )
  
  #'  -------------------------------
  #####  Bottom-up model iterations  #####
  #'  -------------------------------
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup, template = model_template, 
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "BottomUp"),
    future.seed = TRUE
  )
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_inter),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup_inter, template = model_template, 
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "BottomUp_Interference"),
    future.seed = TRUE
  )
 
  #'  -----------------------------
  ####  Directed-Separation Tests  ####
  #'  -----------------------------
  #'  -------------------------------------------
  #####  Top-down interference model iterations  #####
  #'  -------------------------------------------
  #'  Load all iterations of the JAGS model
  all_results_topdown_inter <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/TopDown_Interference", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values
  y_list <- list(data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$elk.t_hat)
  
  #'  Create list of posterior distributions for coefficient of interest
  mod_out <- list()
  for(i in 1:length(all_results_topdown_inter)) {
    mod_out[[i]] <- all_results_topdown_inter[[i]]$fit$sims.list
  }
  post_list_topdown_inter <- list(mod_out[[1]]$beta.wtd[,1], mod_out[[2]]$beta.wtd[,1])
  
  #'  Calculate p.rope value for each iteration of the d-Sep test
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_topdown_inter_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_topdown_inter, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_topdown_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(p.rope_topdown_inter_list)[i] <- list_name
  }
  head(p.rope_topdown_inter_list)
  
  save(p.rope_topdown_inter_list, file = "./Outputs/SEM/JAGS_out/d_Sep/p.ROPE_topedown_inter.RData")
  
  
  
  
