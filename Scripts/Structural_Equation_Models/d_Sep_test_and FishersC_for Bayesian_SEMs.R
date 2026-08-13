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
  ni <- 150000
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
  
  #'  Bayesian p-value for d-Sep test and Fisher's C statistic
  bayes_pvalue <- function(post) {
    p_pos <- mean(post > 0)
    p_neg <- mean(post < 0)
    2 * min(p_pos, p_neg)
  }
  
  #'  Fisher's C using Bayesian p-values
  fishers_C <- function(pval, n_iter = NULL) {
    if(is.null(n_iter)) {
      eps <- 1e-10
    } else {
      eps <- 1 / n_iter # smallest resolvable p-value given draws
    }
    pval[pval <= 0] <- eps
    pval[pval >= 1] <- 1 - eps
    
    k <- length(pval)
    C <- -2 * sum(log(pval))
    df <- 2 * k
    pC <- 1 - pchisq(C, df)
    
    out <- list(claims.p = pval, k = k, C = C, df = df, p.value = pC)
    return(out)
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
  
  #'  ---------------------------
  ####  Basic set for each SEM  ####
  #'  ---------------------------
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
  
  #'  After 1st set of d-Sep tests
  dag_bottomup_inter_updated <- DAG(lion.t ~ lion.tmin1 + wtd.tmin1 + wolf.tmin1, 
                            wolf.t ~ wolf.tmin1 + elk.tmin1 + moose.tmin1,
                            bear.t ~ bear.tmin1 + forest.tmin1 + wolf.tmin1, 
                            coy.t ~ coy.tmin1 + wtd.tmin1 + wolf.tmin1,
                            elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1,
                            moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                            wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1 + bear.t + coy.tmin1)
  
  bs_bottomup_inter_updated <- basic_set(dag_bottomup_inter_updated)
  
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
  
  #'  ------------------------------------------------------
  ####  Build iterative JAGS models for d-separation tests  ####
  #'  ------------------------------------------------------
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
  build_model_string <- function(iter_config, template, registry) {
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
        #'  Non-focal time t: use original SEM covariates from the registry
        orig <- registry[[(r * 2) - 1]]
        #'  Grab covariate, species, and index numbers of terms NOT included in d-sep test
        covs <- orig$covs
        spp <- orig$spp
        indices <- orig$indices
        mu_lines[(r * 2) - 1] <- build_individual_submodels(r, covs, spp, indices, timestep = "")
      }
      
      #'  time t-1: always intercept-only for all regressions and iterations
      mu_lines[r * 2] <- build_individual_submodels(r, NULL, timestep = ".tmin1")
    }
    
    do.call(sprintf, c(list(template), as.list(mu_lines)))
  }
  
  #'  Function to call JAGS and run a single iteration of the model
  run_dSep_iterations <- function(i, iterations, template, registry, data_bundle, listInits, model_name) {  
    
    iter_config <- iterations[[i]]
    
    #'  Call function to build the full model string with custom regressions for d-Sep test
    model_string <- build_model_string(iter_config, template, registry)
    
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
  
  
  #'  -------------------------------
  ####  Iterate through d-Sep tests  ####
  #'  -------------------------------
  #'  -------------------------------------------
  #####  Top-down interference model iterations  #####
  #'  -------------------------------------------
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.t
    list(covs = c("lion.tmin1", "wolf.tmin1", "lionHarv.tmin1"), spp = c(".lion", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    #'  Regression 2: lion.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 3: wolf.t
    list(covs = c("wolf.tmin1", "wolfHarv.tmin1"), spp = c(".wolf", ".harvest"), indices = as.integer(c(1,1))),
    #'  Regression 4: wolf.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 5: bear.t
    list(covs = c("bear.tmin1", "wolf.tmin1", "bearHarv.tmin1"), spp = c(".bear", ".wolf", ".harvest"), indices = as.integer(c(1,1,1))),
    #'  Regression 6: bear.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 7: coy.t
    list(covs = c("coy.tmin1", "wolf.tmin1", "lion.tmin1"), spp = c(".coy", ".wolf", ".lion"), indices = as.integer(c(1,1,1))),
    #'  Regression 8: coy.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 9: elk.t
    list(covs = c("elk.tmin1", "wolf.tmin1", "lion.tmin1"), spp = c(".elk", ".wolf", ".lion"), indices = as.integer(c(1,1,1))),
    #'  Regression 10: elk.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 11: moose.t
    list(covs = c("moose.tmin1", "wolf.tmin1"), spp = c(".moose", ".wolf"), indices = as.integer(c(1,1))),
    #'  Regression 12: moose.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 13: wtd.t
    list(covs = c("wtd.tmin1", "lion.tmin1"), spp = c(".wtd", ".lion"), indices = as.integer(c(1,1))),
    #'  Regerssion 14: wtd.tmin1
    list(covs = NULL, spp = NULL, indices = NULL)
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_topdown_inter_reduced.R")
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    #'  Apply across every element in list of active regressions
    seq_along(dSep_iterations_topdown_int),
    #'  Call run_dSep_iterations function using specified active regression list, model template, and data/inits prepared for JAGS
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown_int, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "TopDown_Interference"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  ------------------------------
  #####  Top-down model iterations  #####
  #'  ------------------------------
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.t
    list(covs = c("lion.tmin1", "lionHarv.tmin1"), spp = c(".lion", ".harvest"), indices = as.integer(c(1,1))),
    #'  Regression 2: lion.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 3: wolf.t
    list(covs = c("wolf.tmin1", "wolfHarv.tmin1"), spp = c(".wolf", ".harvest"), indices = as.integer(c(1,1))),
    #'  Regression 4: wolf.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 5: bear.t
    list(covs = c("bear.tmin1", "bearHarv.tmin1"), spp = c(".bear", ".harvest"), indices = as.integer(c(1,1))),
    #'  Regression 6: bear.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 7: coy.t
    list(covs = c("coy.tmin1"), spp = c(".coy"), indices = as.integer(c(1))),
    #'  Regression 8: coy.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 9: elk.t
    list(covs = c("elk.tmin1", "wolf.tmin1", "lion.tmin1", "elkHarv.tmin1"), spp = c(".elk", ".wolf", ".lion", ".harvest"), indices = as.integer(c(1,1,1,1))),
    #'  Regression 10: elk.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 11: moose.t
    list(covs = c("moose.tmin1", "wolf.tmin1"), spp = c(".moose", ".wolf"), indices = as.integer(c(1,1))),
    #'  Regression 12: moose.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 13: wtd.t
    list(covs = c("wtd.tmin1", "lion.tmin1", "deerHarv.tmin1"), spp = c(".wtd", ".lion", ".harvest"), indices = as.integer(c(1,1,1))),
    #'  Regerssion 14: wtd.tmin1
    list(covs = NULL, spp = NULL, indices = NULL)
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_topdown_reduced.R")
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_topdown),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_topdown, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "TopDown"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  # source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_dsep_template_registry_bottomup_inter_reduced.R")
  sem_registry <- list(
    #'  Regression 1: lion.t
    list(covs = c("lion.tmin1", "wtd.tmin1", "wolf.tmin1"), spp = c(".lion", ".wtd", ".wolf"), indices = as.integer(c(1,1,1))),
    #'  Regression 2: lion.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 3: wolf.t
    list(covs = c("wolf.tmin1", "elk.tmin1", "moose.tmin1"), spp = c(".wolf", ".elk", ".moose"), indices = as.integer(c(1,1,1))),
    #'  Regression 4: wolf.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 5: bear.t
    list(covs = c("bear.tmin1", "forest.tmin1", "wolf.tmin1"), spp = c(".bear", ".forest", ".wolf"), indices = as.integer(c(1,1,1))),
    #'  Regression 6: bear.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 7: coy.t
    list(covs = c("coy.tmin1", "wtd.tmin1", "wolf.tmin1"), spp = c(".coy", ".wtd", ".wolf"), indices = as.integer(c(1,1,1))),
    #'  Regression 8: coy.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 9: elk.t
    list(covs = c("elk.tmin1", "forest.tmin1", "wsi.tmin1"), spp = c(".elk", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regression 10: elk.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 11: moose.t
    list(covs = c("moose.tmin1", "forest.tmin1", "wsi.tmin1"), spp = c(".moose", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regression 12: moose.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 13: wtd.t
    list(covs = c("wtd.tmin1", "forest.tmin1", "wsi.tmin1"), spp = c(".wtd", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regerssion 14: wtd.tmin1
    list(covs = NULL, spp = NULL, indices = NULL)
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_bottomup_inter_reduced.R")
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_int),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup_int, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "BottomUp_Interference"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  -------------------------------
  #####  Bottom-up model iterations  #####
  #'  -------------------------------
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.t
    list(covs = c("lion.tmin1", "elk.tmin1", "wtd.tmin1"), spp = c(".lion", ".elk", ".wtd"), indices = as.integer(c(1,1,1))),
    #'  Regression 2: lion.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 3: wolf.t
    list(covs = c("wolf.tmin1", "elk.tmin1", "moose.tmin1"), spp = c(".wolf", ".elk", ".moose"), indices = as.integer(c(1,1,1))),
    #'  Regression 4: wolf.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 5: bear.t
    list(covs = c("bear.tmin1", "elk.tmin1", "forest.tmin1"), spp = c(".bear", ".elk", ".forest"), indices = as.integer(c(1,1,1))),
    #'  Regression 6: bear.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 7: coy.t
    list(covs = c("coy.tmin1", "wtd.tmin1"), spp = c(".coy", ".wtd"), indices = as.integer(c(1,1))),
    #'  Regression 8: coy.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 9: elk.t
    list(covs = c("elk.tmin1", "forest.tmin1", "wsi.tmin1"), spp = c(".elk", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regression 10: elk.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 11: moose.t
    list(covs = c("moose.tmin1", "forest.tmin1", "wsi.tmin1"), spp = c(".moose", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regression 12: moose.tmin1
    list(covs = NULL, spp = NULL, indices = NULL),
    #'  Regression 13: wtd.t
    list(covs = c("wtd.tmin1", "forest.tmin1", "wsi.tmin1"), spp = c(".wtd", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regerssion 14: wtd.tmin1
    list(covs = NULL, spp = NULL, indices = NULL)
  )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_test_iterations_bottomup_reduced.R")
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle, listInits = initsList, model_name = "BottomUp"),
    future.seed = TRUE
  )
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  
  #'  ---------------------------------------------
  ####  Calculate p-values for d-Separation tests  ####
  #'  ---------------------------------------------
  #'  -------------------------------------------
  #####  Top-down interference model iterations  #####
  #'  -------------------------------------------
  #'  Load all iterations of the JAGS model
  all_results_topdown_inter <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/TopDown_Interference", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test
  y_list <- list(data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,    #15 (iteration/regression number)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$moose.t_hat,
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,       #30
                 data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wtd.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,    #45
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, 
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,      #51 (skipped #52)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, 
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat,      #59 (skipped #58)
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,       #62 (skipped #63)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$coy.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$lion.t_hat)                                                               #73
  #'  Leaves you with 70 d-Sep tests that were possible given the constructs of space and time and our data
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta indices. Most will be [,1]
  #'  but some will be [,2] where the same beta name was used twice in the same regression.
  mod_out <- list()
  for(i in 1:length(all_results_topdown_inter)) {
    mod_out[[i]] <- all_results_topdown_inter[[i]]$fit$sims.list
  }
  post_list_topdown_inter <- list(mod_out[[1]]$beta.wtd[,1], mod_out[[2]]$beta.wtd[,1], mod_out[[3]]$beta.wtd[,1],
                                  mod_out[[4]]$beta.wtd[,1], mod_out[[5]]$beta.wtd[,1], mod_out[[6]]$beta.wtd[,1],
                                  mod_out[[7]]$beta.moose[,1], mod_out[[8]]$beta.moose[,1], mod_out[[9]]$beta.moose[,1],
                                  mod_out[[10]]$beta.moose[,1], mod_out[[11]]$beta.moose[,1], mod_out[[12]]$beta.moose[,1],
                                  mod_out[[13]]$beta.elk[,1], mod_out[[14]]$beta.elk[,1], mod_out[[15]]$beta.elk[,1],
                                  mod_out[[16]]$beta.elk[,1], mod_out[[17]]$beta.elk[,1], mod_out[[18]]$beta.elk[,1],
                                  mod_out[[19]]$beta.coy[,1], mod_out[[20]]$beta.coy[,1], mod_out[[21]]$beta.coy[,1],
                                  mod_out[[22]]$beta.coy[,1], mod_out[[23]]$beta.coy[,1], mod_out[[24]]$beta.harvest[,1],
                                  mod_out[[25]]$beta.harvest[,1], mod_out[[26]]$beta.harvest[,2], mod_out[[27]]$beta.harvest[,1], # note the different indexing
                                  mod_out[[28]]$beta.harvest[,1], mod_out[[29]]$beta.harvest[,1], mod_out[[30]]$beta.harvest[,2], # note the different indexing
                                  mod_out[[31]]$beta.bear[,1], mod_out[[32]]$beta.bear[,1], mod_out[[33]]$beta.bear[,1],
                                  mod_out[[34]]$beta.bear[,1], mod_out[[35]]$beta.bear[,1], mod_out[[36]]$beta.bear[,1],
                                  mod_out[[37]]$beta.harvest[,1], mod_out[[38]]$beta.harvest[,2], mod_out[[39]]$beta.harvest[,1], # note the different indexing
                                  mod_out[[40]]$beta.harvest[,1], mod_out[[41]]$beta.harvest[,1], mod_out[[42]]$beta.harvest[,2], # note the different indexing
                                  mod_out[[43]]$beta.harvest[,1], mod_out[[44]]$beta.harvest[,2], mod_out[[45]]$beta.harvest[,2], # note the different indexing
                                  mod_out[[46]]$beta.harvest[,1], mod_out[[47]]$beta.harvest[,1], mod_out[[48]]$beta.harvest[,1], 
                                  mod_out[[49]]$beta.wolf[,1], mod_out[[50]]$beta.moose[,2], mod_out[[51]]$beta.moose[,2],        # note the different indexing
                                  mod_out[[52]]$beta.moose[,2], mod_out[[53]]$beta.moose[,2], mod_out[[54]]$beta.moose[,2],       # note the different indexing
                                  mod_out[[55]]$beta.moose[,2], mod_out[[56]]$beta.bear[,2], mod_out[[57]]$beta.bear[,2],         # note the different indexing
                                  mod_out[[58]]$beta.bear[,2], mod_out[[59]]$beta.bear[,2], mod_out[[60]]$beta.bear[,2],          # note the different indexing
                                  mod_out[[61]]$beta.wolf[,2], mod_out[[62]]$beta.wolf[,2], mod_out[[63]]$beta.wolf[,2],          # note the different indexing
                                  mod_out[[64]]$beta.wolf[,2], mod_out[[65]]$beta.wtd[,2], mod_out[[66]]$beta.wtd[,2],            # note the different indexing
                                  mod_out[[67]]$beta.wtd[,2], mod_out[[68]]$beta.elk[,2], mod_out[[69]]$beta.elk[,2],             # note the different indexing
                                  mod_out[[70]]$beta.coy[,2])
  
  #'  ----------------------------------------------------------------
  #####  Calculate p.rope value for each iteration of the d-Sep test  #####
  #'  ----------------------------------------------------------------
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_topdown_inter_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_topdown_inter, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_topdown_inter_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(p.rope_topdown_inter_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_topdown_inter_df <- stack(p.rope_topdown_inter_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4))
  
  #'  --------------------------------------------
  ####  Calculate Bayesian p-value as d-Sep test  ####
  #'  --------------------------------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_topdown_inter_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown_inter, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_topdown_inter_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(bayes.p_topdown_inter_list)[i] <- list_name
  }
  bayes.p_topdown_inter_df <- stack(bayes.p_topdown_inter_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4))
  
  #'  Join both d-Sep test p-values and save
  p.val_topdown_inter_df <- full_join(p.rope_topdown_inter_df, bayes.p_topdown_inter_df, by = "iteration")
  
  write_csv(p.val_topdown_inter_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_topdown_inter.csv")
  
  #'  ---------------
  #####  Fisher's C  #####
  #'  ---------------
  fishers.C_topdown_inter <- fishers_C(pval = bayes.p_topdown_inter_df$bayes.p, n_iter = nrow(bayes.p_topdown_inter_df))
  print(fishers.C_topdown_inter)
  
  
  #'  ------------------------------
  #####  Top-down model iterations  #####
  #'  ------------------------------
  #'  Load all iterations of the JAGS model
  all_results_topdown <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/TopDown", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test
  y_list <- list(data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$moose.t_hat,
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$moose.t_hat,
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,      #15 (iteration number)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$moose.t_hat,
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$moose.t_hat,
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$lion.t_hat,      #30
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wolf.t_hat,
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wolf.t_hat,    #39 (skipped #37, #38, #40, #41)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat,       #46 (skipped #44 & #45)
                 data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat,    
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$moose.t_hat, 
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat,      
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wolf.t_hat,    #62 (skipped #59 & #60)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat,       #67 (skipped #63 & #64)
                 data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat,       
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat,       #78 (skipped #75 & #76)
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat,       #83 (skipped #80 & #81)
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$elk.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$lion.t_hat)      #89 
  #'  Leaves you with 75 d-Sep tests that were possible given the constructs of space and time and our data
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta indices. Most will be [,1]
  #'  but some will be [,2] where the same beta name was used twice in the same regression.
  mod_out <- list()
  for(i in 1:length(all_results_topdown)) {
    mod_out[[i]] <- all_results_topdown[[i]]$fit$sims.list
  }
  post_list_topdown <- list(mod_out[[1]]$beta.harvest[,1], mod_out[[2]]$beta.harvest[,2], mod_out[[3]]$beta.harvest[,1],     # note the different indexing
                            mod_out[[4]]$beta.harvest[,2], mod_out[[5]]$beta.harvest[,2], mod_out[[6]]$beta.harvest[,2],     # note the different indexing
                            mod_out[[7]]$beta.wtd[,1], mod_out[[8]]$beta.wtd[,1], mod_out[[9]]$beta.wtd[,1],
                            mod_out[[10]]$beta.wtd[,1], mod_out[[11]]$beta.wtd[,1], mod_out[[12]]$beta.wtd[,1],
                            mod_out[[13]]$beta.moose[,1], mod_out[[14]]$beta.moose[,1], mod_out[[15]]$beta.moose[,1],
                            mod_out[[16]]$beta.moose[,1], mod_out[[17]]$beta.moose[,1], mod_out[[18]]$beta.moose[,1],
                            mod_out[[19]]$beta.harvest[,1], mod_out[[20]]$beta.harvest[,2], mod_out[[21]]$beta.harvest[,1],  # note the different indexing
                            mod_out[[22]]$beta.harvest[,2], mod_out[[23]]$beta.harvest[,2], mod_out[[24]]$beta.harvest[,1],  # note the different indexing
                            mod_out[[25]]$beta.elk[,1], mod_out[[26]]$beta.elk[,1], mod_out[[27]]$beta.elk[,1], 
                            mod_out[[28]]$beta.elk[,1], mod_out[[29]]$beta.elk[,1], mod_out[[30]]$beta.elk[,1], 
                            mod_out[[31]]$beta.coy[,1], mod_out[[32]]$beta.coy[,1], mod_out[[33]]$beta.coy[,1],
                            mod_out[[34]]$beta.coy[,1], mod_out[[35]]$beta.coy[,1], mod_out[[36]]$beta.coy[,1],
                            mod_out[[37]]$beta.coy[,2], mod_out[[38]]$beta.coy[,2], mod_out[[39]]$beta.coy[,2],              # note the different indexing
                            mod_out[[40]]$beta.coy[,1], mod_out[[41]]$beta.coy[,1], mod_out[[42]]$beta.coy[,1],              # note the different indexing
                            mod_out[[43]]$beta.harvest[,1], mod_out[[44]]$beta.harvest[,2], mod_out[[45]]$beta.harvest[,2],  # note the different indexing
                            mod_out[[46]]$beta.harvest[,2], mod_out[[47]]$beta.harvest[,2], mod_out[[48]]$beta.bear[,1],     # note the different indexing
                            mod_out[[49]]$beta.bear[,1], mod_out[[50]]$beta.bear[,1], mod_out[[51]]$beta.bear[,1],        
                            mod_out[[52]]$beta.bear[,1], mod_out[[53]]$beta.bear[,2], mod_out[[54]]$beta.bear[,2],           # note the different indexing
                            mod_out[[55]]$beta.bear[,2], mod_out[[56]]$beta.bear[,2], mod_out[[57]]$beta.bear[,2],           # note the different indexing
                            mod_out[[58]]$beta.harvest[,1], mod_out[[59]]$beta.harvest[,2], mod_out[[60]]$beta.harvest[,2],  # note the different indexing
                            mod_out[[61]]$beta.harvest[,2], mod_out[[62]]$beta.wolf[,1], mod_out[[63]]$beta.wolf[,1],        # note the different indexing
                            mod_out[[64]]$beta.moose[,2], mod_out[[65]]$beta.moose[,2], mod_out[[66]]$beta.moose[,2],        # note the different indexing
                            mod_out[[67]]$beta.moose[,2], mod_out[[68]]$beta.wolf[,2], mod_out[[69]]$beta.wolf[,2],          # note the different indexing
                            mod_out[[70]]$beta.wolf[,2], mod_out[[71]]$beta.harvest[,2], mod_out[[72]]$beta.harvest[,2],
                            mod_out[[73]]$beta.wtd[,2], mod_out[[74]]$beta.wtd[,2], mod_out[[75]]$beta.elk[,2])
  
  #'  Calculate p.rope value for each iteration of the d-Sep test 
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_topdown_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_topdown, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_topdown_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(p.rope_topdown_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_topdown_df <- stack(p.rope_topdown_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4))
  
  #'  --------------------------------------------
  ####  Calculate Bayesian p-value as d-Sep test  ####
  #'  --------------------------------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_topdown_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_topdown_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(bayes.p_topdown_list)[i] <- list_name
  }
  bayes.p_topdown_df <- stack(bayes.p_topdown_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4))
  
  #'  Join both d-Sep test p-values and save
  p.val_topdown_df <- full_join(p.rope_topdown_df, bayes.p_topdown_df, by = "iteration")
  
  write_csv(p.val_topdown_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_topdown.csv")
  
  #'  ---------------
  #####  Fisher's C  #####
  #'  ---------------
  fishers.C_topdown <- fishers_C(pval = bayes.p_topdown_df$bayes.p, n_iter = nrow(bayes.p_topdown_df))
  print(fishers.C_topdown)
  
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Load all iterations of the JAGS model
  all_results_bottomup_inter <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/BottomUp_Interference", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test
  y_list <- list(data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$coy.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$elk.t_hat,
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$coy.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$elk.t_hat,     #15 (iteration/regression number)
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$bear.t_hat,
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,      #29 (skipped #25)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat,       #34 (skipped #30)
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,     #41 (skipped #39)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$wtd.t_hat,        #50 (skipped #42)
                 data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wtd.t_hat,       #55 (skipped #52 & #54)
                 data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$coy.t_hat,       #59 (skipped #57)
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$lion.t_hat)                                  #63 (skipped #60 & #62)
  #'  Leaves you with 47 d-Sep tests that were possible given the constructs of space and time and our data
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta indices. Most will be [,1]
  #'  but some will be [,2] where the same beta name was used twice in the same regression.
  #'  Indexing for beta.wsi[1] and beta.forest[1] are slightly different b/c only one coeff
  #'  was generated per parameter for these.
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
                                   mod_out[[25]]$beta.moose[,2], mod_out[[26]]$beta.moose[,2], mod_out[[27]]$beta.moose[,2], # note the different indexing
                                   mod_out[[28]]$beta.moose[,2], mod_out[[29]]$beta.moose[,2], mod_out[[30]]$beta.moose[,2], # note the different indexing
                                   mod_out[[31]]$beta.elk[,1], mod_out[[32]]$beta.elk[,1], mod_out[[33]]$beta.elk[,1],
                                   mod_out[[34]]$beta.elk[,1], mod_out[[35]]$beta.elk[,2], mod_out[[36]]$beta.elk[,2],       # note the different indexing
                                   mod_out[[37]]$beta.elk[,2], mod_out[[38]]$beta.elk[,2], mod_out[[39]]$beta.bear[,2],      # note the different indexing
                                   mod_out[[40]]$beta.bear[,2], mod_out[[41]]$beta.bear[,2], mod_out[[42]]$beta.wolf[,2],    # note the different indexing
                                   mod_out[[43]]$beta.wolf[,2], mod_out[[44]]$beta.wolf[,2], mod_out[[45]]$beta.wtd[,2],     # note the different indexing
                                   mod_out[[46]]$beta.wtd[,2], mod_out[[47]]$beta.coy[,2])                                   # note the different indexing
  
  #'  Calculate p.rope value for each iteration of the d-Sep test
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_bottomup_inter_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_bottomup_inter, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_bottomup_inter_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(p.rope_bottomup_inter_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_bottomup_inter_df <- stack(p.rope_bottomup_inter_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4))
  
  #'  --------------------------------------------
  ####  Calculate Bayesian p-value as d-Sep test  ####
  #'  --------------------------------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_bottomup_inter_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup_inter, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_bottomup_inter_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(bayes.p_bottomup_inter_list)[i] <- list_name
  }
  bayes.p_bottomup_inter_df <- stack(bayes.p_bottomup_inter_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4))
  
  #'  Join both d-Sep test p-values and save
  p.val_bottomup_inter_df <- full_join(p.rope_bottomup_inter_df, bayes.p_bottomup_inter_df, by = "iteration")
  
  write_csv(p.val_bottomup_inter_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_bottomup_inter.csv")
  
  #'  ---------------
  #####  Fisher's C  #####
  #'  ---------------
  fishers.C_bottomup_inter <- fishers_C(pval = bayes.p_bottomup_inter_df$bayes.p, n_iter = nrow(bayes.p_bottomup_inter_df))
  print(fishers.C_bottomup_inter)
  
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Load all iterations of the JAGS model
  all_results_bottomup <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/BottomUp", full.names = TRUE), readRDS)
  
  #'  Create list of "observed" values of focal response variable, one per d-Sep test
  y_list <- list(data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wtd.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$wolf.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$moose.t_hat, data_JAGS_bundle$wtd.t_hat,     #15 (iteration/regression number)
                 data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$wolf.t_hat,
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat,
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$lion.t_hat,
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$elk.t_hat,        #30 (skipped #25, #26, #29)
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$lion.t_hat,     #34 (skipped #33)
                 data_JAGS_bundle$wtd.t_hat, data_JAGS_bundle$coy.t_hat, data_JAGS_bundle$elk.t_hat,
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$elk.t_hat,     
                 data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$coy.t_hat,        
                 data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$bear.t_hat, data_JAGS_bundle$wolf.t_hat,      #47 (skipped #48 & #50)
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$elk.t_hat, data_JAGS_bundle$bear.t_hat,      #52 (skipped #57)
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$bear.t_hat,     #56 (skipped #54)
                 data_JAGS_bundle$wolf.t_hat, data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$wolf.t_hat,     #60 (skipped #58)
                 data_JAGS_bundle$lion.t_hat, data_JAGS_bundle$lion.t_hat)                                  #64 (skipped #61 & #63)
  #'  Leaves you with 53 d-Sep tests that were possible given the constructs of space and time and our data
  
  #'  Create list of posterior distributions for coefficient of interest, one per d-Sep test
  #'  Pay close attention to the indexing, especially with the beta indices. Most will be [,1]
  #'  but some will be [,2] where the same beta name was used twice in the same regression.
  #'  Indexing for beta.wsi[1] and beta.forest[1] are slightly different b/c only one coeff
  #'  was generated per parameter for these.
  mod_out <- list()
  for(i in 1:length(all_results_bottomup)) {
    mod_out[[i]] <- all_results_bottomup[[i]]$fit$sims.list
  }
  post_list_bottomup <- list(mod_out[[1]]$beta.wsi[1], mod_out[[2]]$beta.wsi[1], mod_out[[3]]$beta.wsi[1],
                             mod_out[[4]]$beta.wsi[1], mod_out[[5]]$beta.coy[,1], mod_out[[6]]$beta.coy[,1],
                             mod_out[[7]]$beta.coy[,1], mod_out[[8]]$beta.coy[,1], mod_out[[9]]$beta.coy[,1],
                             mod_out[[10]]$beta.coy[,1], mod_out[[11]]$beta.forest[1], mod_out[[12]]$beta.forest[1],
                             mod_out[[13]]$beta.forest[1], mod_out[[14]]$beta.bear[,1], mod_out[[15]]$beta.bear[,1],
                             mod_out[[16]]$beta.bear[,1], mod_out[[17]]$beta.bear[,1], mod_out[[18]]$beta.bear[,1],
                             mod_out[[19]]$beta.bear[,1], mod_out[[20]]$beta.moose[,1], mod_out[[21]]$beta.moose[,1],
                             mod_out[[22]]$beta.moose[,1], mod_out[[23]]$beta.moose[,1], mod_out[[24]]$beta.moose[,1], 
                             mod_out[[25]]$beta.moose[,2], mod_out[[26]]$beta.moose[,2], mod_out[[27]]$beta.moose[,2], # note the different indexing
                             mod_out[[28]]$beta.moose[,2], mod_out[[29]]$beta.moose[,2], mod_out[[30]]$beta.moose[,2], # note the different indexing
                             mod_out[[31]]$beta.wolf[,1], mod_out[[32]]$beta.wolf[,1], mod_out[[33]]$beta.wolf[,1],
                             mod_out[[34]]$beta.wolf[,1], mod_out[[35]]$beta.wolf[,1], mod_out[[36]]$beta.wtd[,1],       
                             mod_out[[37]]$beta.wtd[,1], mod_out[[38]]$beta.wtd[,1], mod_out[[39]]$beta.wtd[,2],       # note the different indexing
                             mod_out[[40]]$beta.wtd[,2], mod_out[[41]]$beta.wtd[,2], mod_out[[42]]$beta.wtd[,2],       # note the different indexing
                             mod_out[[43]]$beta.wtd[,2], mod_out[[44]]$beta.coy[,2], mod_out[[45]]$beta.coy[,2],       # note the different indexing
                             mod_out[[46]]$beta.coy[,2], mod_out[[47]]$beta.coy[,2], mod_out[[47]]$beta.elk[,2],       # note the different indexing
                             mod_out[[43]]$beta.elk[,2], mod_out[[44]]$beta.elk[,2], mod_out[[45]]$beta.bear[,2],      # note the different indexing
                             mod_out[[46]]$beta.bear[,2], mod_out[[47]]$beta.wolf[,2])                                 # note the different indexing                            
  
  #'  Calculate p.rope value for each iteration of the d-Sep test
  p_rope_iterations <- function(y_dat, post_beta) { 
    p.rope_val <- p.rope(y = y_dat, post = post_beta)
    print(p.rope_val)
    return(p.rope_val)
  }
  p.rope_bottomup_list <- mapply(p_rope_iterations, y_dat = y_list, post_beta = post_list_bottomup, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_bottomup_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(p.rope_bottomup_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_bottomup_df <- stack(p.rope_bottomup_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4))
  
  #'  --------------------------------------------
  ####  Calculate Bayesian p-value as d-Sep test  ####
  #'  --------------------------------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_bottomup_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_bottomup_list)) {
    list_name <- sprintf("p.rope.%03d", i)
    names(bayes.p_bottomup_list)[i] <- list_name
  }
  bayes.p_bottomup_df <- stack(bayes.p_bottomup_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4))
  
  #'  Join both d-Sep test p-values and save
  p.val_bottomup_df <- full_join(p.rope_bottomup_df, bayes.p_bottomup_df, by = "iteration")
  
  write_csv(p.val_bottomup_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_bottomup.csv")
  
  #'  ---------------
  #####  Fisher's C  #####
  #'  ---------------
  fishers.C_bottomup <- fishers_C(pval = bayes.p_bottomup_df$bayes.p, n_iter = nrow(bayes.p_bottomup_df))
  print(fishers.C_bottomup)
  
