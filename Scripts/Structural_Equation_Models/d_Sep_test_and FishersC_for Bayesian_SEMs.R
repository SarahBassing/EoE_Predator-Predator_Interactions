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
  #' #'  Source functions to setup data for JAGS
  #' source("./Scripts/Structural_Equation_Models/JAGS_setup_functions.R")
  
  #' #'  Bundle data for JAGS
  #' data_JAGS_bundle <- bundle_dat(post_summaries, covs = covs_ztransformed, 
  #'                                nwolf = 2, nlion = 2, nbear = 2, ncoy = 2, nelk = 2, 
  #'                                nmoose = 2, nwtd = 2, nharv = 2, nfor = 1, nwsi = 1)
  #' 
  #' #'  Draw initial values for each chain in JAGS
  #' num.chains <- 3
  #' #'  Create empty list
  #' initsList <- vector('list', num.chains)
  #' #'  Setting seed for reproducibility
  #' set.seed(9714)
  #' #'  Loop through generate_inits function 3 times (1 for each chain) 
  #' for(i in 1:num.chains){
  #'   initsList[[i]] <- generate_inits(nwolf = 2, nlion = 2, nbear = 2, ncoy = 2, nelk = 2, 
  #'                                    nmoose = 2, nwtd = 2, nharv = 2, nfor = 1, nwsi = 1)
  #' }
  
  #'  Parameters to be monitored
  params <- c("beta.int", "beta.int.tmin1", "beta.wolf", "beta.lion", "beta.bear", "beta.coy", "beta.elk", 
              "beta.moose", "beta.wtd", "beta.harvest", "beta.wsi","beta.forest", "beta.road", "beta.public",
              "sigma.spp", "sigma.spp.tmin1")  
  
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
    #'  Ensure y is a vector
    y_vec <- as.vector(y)
    #'  Calculate SD, after removing any NAs in y
    sd_y <- sd(y_vec, na.rm = TRUE)
    #'  Average number of posterior draws that fall above lower ROPE value AND below upper ROPE value
    mean(-pct * sd_y < post & pct * sd_y > post)
    
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
    View(basicset)
    #' #'  Return basic set with conditional independence claims involving 3+ variables
    #' basicset_3plus <- sapply(basicset, length) >= 3
    #' basicset_skinny <- basicset[basicset_3plus]
    #' #'  Return basic set of unconditional independence tests (marginal independence)
    #' basicset_marginal_ind <- sapply(basicset, length) < 3
    #' basicset_marginal <- basicset[basicset_marginal_ind]
    #' print(length(basicset_skinny))
    #' View(basicset_skinny)
    #' basicset_to_test <- list(basicset_skinny, basicset_marginal)
    return(basicset)
  }
  
  #'  ---------------------------
  ####  Basic set for each SEM  ####
  #'  ---------------------------
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
  
  #'  ---------------------------------
  #####  Bottom-up interference model  ####
  #'  ---------------------------------
  dag_bottomup_inter <- DAG(lion.t ~ wtd.tmin1 + wolf.tmin1, # lion.tmin1 + 
                            wolf.t ~ wolf.tmin1 + elk.tmin1 + moose.tmin1,
                            bear.t ~ bear.tmin1 + forest.tmin1 + wolf.tmin1,
                            coy.t ~ coy.tmin1 + wtd.tmin1 + wolf.tmin1,
                            elk.t ~ elk.tmin1 + forest.tmin1 + wsi.tmin1,
                            moose.t ~ moose.tmin1 + forest.tmin1 + wsi.tmin1,
                            wtd.t ~ wtd.tmin1 + forest.tmin1 + wsi.tmin1)
  
  bs_bottomup_inter <- basic_set(dag_bottomup_inter)
  
  
  
  #'  --------------------
  #####  Bottom-up model  ####
  #'  --------------------
  dag_bottomup <- DAG(lion.t ~ elk.tmin1 + wtd.tmin1, lion.tmin1 + 
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
  #'  -------------------------------------------------------
  #####  Functions for when variable A t-1 --> variable B t  #####
  #'  These functions also work for when variable A t --> variable B t
  #'  -------------------------------------------------------
  #'  Source JAGS template
  source("./Scripts/Structural_Equation_Models/Bayesian_SEM/JAGS_SEM_dsep_template.R")
  
  #'  Function to build custom regressions for each iteration
  build_individual_submodels <- function(reg_num, covariates = NULL, spp = NULL, 
                                         indices = NULL, lags = NULL) {   #, beta_prefix = "beta"
    
    
    #' #'  Create intercept and slope parameter names based on time step
    #' #'  beta_prefix indicates whether to use main model terms (beta.int) or the 
    #' #'  auxilary beta terms (beta.aux). Only needed for d-Sep claims where t-1
    #' #'  affects t-1. This helps keep the rest of the main model beta arrays and 
    #' #'  number of regressions unchanged.
    #' beta0_array <- if(beta_prefix == "beta") "beta.int" else "beta0.aux"
    beta0_array <- "beta.int"    # intercept array 
    
    #' #'  If-Else statement for regressions
    #' if(is.null(covariates) || length(covariates) == 0) {
    #'   sprintf("%s[%d]", beta0_array, reg_num)
    #' } else {
    #'   if(is.null(lags)) lags <- rep("y-1", length(covariates))
    #'   terms <- paste0(sprintf("%s%s[%d] * %s[i,%s]", beta_prefix, spp, indices, covariates, lags),
    #'                   collapse = " + ")
    #'   if(beta_prefix == "beta") {
    #'     sprintf("%s[%d] + %s", beta0_array, reg_num, terms)
    #'   } else {
    #'     sprintf("%s + %s", beta0_array, terms)  # aux intercept has no [reg_num] index
    #'   }
    #' }
    
    #'  If covariate is null or 0 (i.e., intercept only regressions)
    if(is.null(covariates) || length(covariates) == 0) {
      sprintf("%s[%d]", beta0_array, reg_num)
    } else {
      #'  Build one string per "beta * covariate" by filling placeholders with
      #'  specified character strings or integers (sprintf is vectorized so does
      #'  this in order that strings/integers are provided)
      #'  %s is placeholder for character strings; %d is placeholder for integers
      if(is.null(lags)) lags <- rep("y-1", length(covariates))
      if (length(lags) != length(covariates)) {
        stop("lags must be NULL or the same length as covariates")
      }

    terms <- paste0(sprintf("beta%s[%d] * %s[i,%s]", spp, indices, covariates, lags),
                    collapse = " + ")
    sprintf("%s[%d] + %s", beta0_array, reg_num, terms)
    }
  }
  
  #' #'  Function to build auxiliary regressions to test same-year independence claims
  #' #'  Essentially generates additional priors and likelihood for cases where the
  #' #'  original model structure did not include a regression for the potential 
  #' #'  relationship flagged in the basic set.
  #' build_aux_submodel <- function(outcome_spp_name, covariates, spp, indices, lags) {
  #'   #'  outcome_spp_name must match those used in .hat / .sigma_hat array
  #'   #'  lags will typically be "y" indicating same year, as opposed to a 1-yr 
  #'   #'  time lag ("y-1").
  #'   linpred <- build_individual_submodels(reg_num = NA, covariates = covariates,
  #'                                         spp = spp, indices = indices, lags = lags,
  #'                                         beta_prefix = "beta.aux")
  #'   n_terms <- length(covariates)
  #'   #'  Generate priors and likelihood for aux data. beta.aux%s, only works if 
  #'   #'  there is a single aux covariate based on how it is currently formulated.
  #'   sprintf("
  #'     beta0.aux ~ dnorm(0, 0.01)
  #'     for (a in 1:%d) { beta.aux%s[a] ~ dnorm(0, 0.01) }
  #'     sigma.aux ~ dunif(0, 2)
  #'     tau.aux <- 1 / pow(sigma.aux, 2)
  #' 
  #'     for (i in 1:nSites) {
  #'       for (y in 1:nYear) {
  #'         %s.hat[i,y] ~ dnorm(%s.aux[i,y], %s.tau_hat[i,y])
  #'         %s.aux[i,y] ~ dnorm(mu.aux[i,y], tau.aux)
  #'         mu.aux[i,y] <- %s
  #'       }
  #'     }
  #'     ", n_terms, paste(unique(spp), collapse = ""),  # a little hacky but it works for this situation
  #'     outcome_spp_name, outcome_spp_name, outcome_spp_name, outcome_spp_name, linpred)
  #' }
  
  #'  Function to assemble full model string for a single iteration
  #'  Requires having sourced the JAGS template model
  #'  mu_lines[(r * 2) - 1] indexing maps each species regression (1:7) onto two 
  #'  consecutive slots in the mu_lines vector (for t and t-1 versions) - ensures 
  #'  regressions are in correct order and match template's 14 %s placeholders
  build_model_string <- function(iter_config, template, registry) {
    #' #'  Identify if regression is for the main independence claims or an auxiliary claim
    #' mode <- if(is.null(iter_config$mode)) "main" else iter_config$mode

    #'  Create empty strings to be filled with each regressions terms
    mu_lines <- character(7)  
    # mu_lines <- character(14)  # 7 regressions for time t and t-1

    #' for(r in 1:7) {
    #'   if(mode == "main" && r == iter_config$dSep_test) {
    #'     #'  Focal species for a MAIN-mode test
    #'     covs <- iter_config$covariates
    #'     spp <- iter_config$spp
    #'     indices <- iter_config$indices
    #'     lags <- iter_config$lags
    #'     mu_lines[r] <- build_individual_submodels(r, covs, spp, indices, lags, beta_prefix = "beta")
    #'   } else {
    #'     #'  Non-focal species, use original SEM covariates from the registry 
    #'     #'  (always, for main and aux modes)
    #'     orig <- registry[[r]]
    #'     #'  Grab covariate, species, and index numbers of terms NOT included in d-sep test
    #'     covs <- orig$covs
    #'     spp <- orig$spp
    #'     indices <- orig$indices
    #'     lags <- orig$lags
    #'     mu_lines[r] <- build_individual_submodels(r, covs, spp, indices, lags)
    #'   }
    #' }
    
    for(r in 1:7) {   # for each regression
      #'  If the regression index is the same as the dSep_test value then...
      if(r == iter_config$dSep_test) {
        covs <- iter_config$covariates
        spp <- iter_config$spp
        indices <- iter_config$indices
        lags <- iter_config$lags        # if NULL, defaults to "y-1
        mu_lines[r] <- build_individual_submodels(r, covs, spp, indices, lags)
      } else {
        #'  Non-focal time t: use original SEM covariates from the registry
        orig <- registry[[r]]
        #'  Grab covariate, species, and index numbers of terms NOT included in d-sep test
        covs <- orig$covs
        spp <- orig$spp
        indices <- orig$indices
        lags <- orig$lags
        mu_lines[r] <- build_individual_submodels(r, covs, spp, indices, lags)
      }
    }
    
    #' #'  8th placeholder: empty for main mode, auxiliary block for aux mode
    #' aux_block <- if (mode == "aux") {
    #'   outcome_spp_name <- iter_config$outcome_spp_name
    #'   covariates <- iter_config$covariates
    #'   spp <- iter_config$spp
    #'   indices <- iter_config$indices
    #'   lags <- iter_config$lags
    #'   build_aux_submodel(outcome_spp_name, covariates, spp, indices, lags)
    #' } else {
    #'   ""
    #' }
    #' 
    #' all_lines <- c(mu_lines, aux_block)
    
    do.call(sprintf, c(list(template), as.list(mu_lines))) #all_lines
  }
  
  #' #'  Function to indicate which parameters to monitor (it will depend on whether
  #' #'  main or aux independence claim is being tested)
  #' get_monitor_params <- function(iter_config, main_params) {
  #'   mode <- if (is.null(iter_config$mode)) "main" else iter_config$mode
  #'   if (mode == "main") {
  #'     main_params  
  #'   } else {
  #'     #'  Aux mode: monitor the intercept + this iteration's specific
  #'     #'  aux beta term(s), whose names depend on this claim's spp/indices
  #'     aux_names <- sprintf("beta.aux%s[%d]", iter_config$spp, iter_config$indices)
  #'     c("beta0.aux", aux_names)
  #'   }
  #' }
  
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
    #' #'  Indicate which paramters to monitor
    #' monitor_params <- get_monitor_params(iter_config, params)
    
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
  
  #'  ---------------------------------------------------------
  #####  Functions for when variable A t-1 --> variable B t-1  #####
  #'  ---------------------------------------------------------
  #'  Function to extract the latent posterior mean and sd from the main model 
  #'  that was originally fitted. This requires that each spp.latent was monitored 
  #'  when the original model was fit.
  extract_latent_post_summaries <- function(og_fit, spp, nSites, nYear) {
    #'  og_fit: fitted jagsUI output from original SEM (fit in Bayesian_SEMs_relative_density_index_1yLag.R)
    #'  spp: species name - must match naming convention use din spp.latent array name
    
    #'  Grab posteriors
    samples_matrix <- as.matrix(og_fit$samples)
    latent_mean <- matrix(NA_real_, nSites, nYear)
    latent_sd <- matrix(NA_real_, nSites, nYear)
    
    for(i in 1:nSites) {
      for(y in 1:nYear) {
        #'  Create column name for specific spp.latent variable indexed by [nSite,nYear]
        col_name <- sprintf("%s.latent[%d,%d]", spp, i, y)
        #'  If created column name is in the colnames extracted from the og_fit samples
        if(col_name %in% colnames(samples_matrix)) {
          #'  Snag those draws from the posterior and save the mean and sd from each iteration
          draws <- samples_matrix[, col_name]
          latent_mean[i,y] <- mean(draws)
          latent_sd[i,y] <- sd(draws)
        }
      }
    }
    list(mean = latent_mean, sd = latent_sd)
  }
  
  #'  Function to fit a standalone regression using the posterior summaries as
  #'  the noisy "observed" response
  #'  Note: z = NULL is a placeholder for if any of these independence claims had
  #'  a condition set of variables. But in this case, z is never needed.
  fit_one_dSep_claim <- function(y, x, z = NULL, n.chains = n.chains, n.adapt = n.adapt, 
                                 n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin, 
                                 model_name, iter_num) {
    ncondvars <- if(is.null(z) || ncol(as.matrix(z)) == 0) 0 else ncol(as.matrix(z))
    #'  Create and fill in input data for JAGS (in list format)
    jd <- list(y = y, x = x, N = length(y))
    if(ncondvars > 0) {
      jd$z <- as.matrix(z)
      jd$ncondvars <- ncondvars
    }
    
    #'  Create custom regression with added variable for independence claims
    if(ncondvars > 0) {
      #'  Build conditioning variable terms
      z_terms <- paste(sprintf("b_z[%d] * z[i,%d]", 1:ncondvars, 1:ncondvars), collapse = " + ")
      
      #'  Create model string for JAGS that can be updated dynamically for each ind. claim
      #'  IF the independence claim includes a condition set of predictors (the 
      #'  "given blah blah blah" variables) create the first model_string. ELSE 
      #'  create the second model_string with only the y (spp.latent) and x (focal predictor)
      model_string <- sprintf("
                              model {
                              #'  Likelihood to be appended with conditioning claim
                              for(i in 1:N) {
                                y[i] ~ dnorm(mu[i], tau)
                                mu[i] <- b0 + b_x * x[i] + %s }
                              
                              #'  Priors
                              b0 ~ dnorm(0, 1e-4)
                              b_x ~ dnorm(0, 1e-4)
                                for(j in 1:ncondvars) {
                                  b_z[j] ~ dnorm(0, 1e-4) }
                              tau ~ dgamma(0.01, 0.01)
                              }
                              ", z_terms)
    } else {
      model_string <- "
      model {
      #'  Likelihood for all other regressions
      for(i in 1:N) {
        y[i] ~ dnorm(mu[i], tau)
        mu[i] <- b0 + b_x * x[i] }
      #'  Priors
      b0 ~ dnorm(0, 1e-4)
      b_x ~ dnorm(0, 1e-4)
      tau ~ dgamma(0.01, 0.01) 
      }
      "
    }
    
    #' #'  Create a temporary file to hold the model_string for jagsUI
    #' temp_file <- tempfile(fileext = ".txt")
    #' writeLines(model_string, temp_file)
    #' on.exit(unlink(temp_file), add = TRUE)
    
    #'  Create temporary directory to save all iterations of the template
    temp_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/temp_models/tmin1", model_name)
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Create name new name for each model where %03d is a placeholder for the 
    #'  iteration number, padded by up to three 0's (e.g., temp_model_001, temp_model_012)
    temp_file <- file.path(temp_dir, sprintf("temp_model_%03d.txt", iter_num))
    writeLines(model_string, temp_file)
    
    #'  Refit model with added independence claim
    fit <- jagsUI::jags(data = jd, inits = NULL, parameters.to.save = "b_x",
                        model.file = temp_file, n.chains = n.chains, n.adapt = n.adapt, 
                        n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin, 
                        parallel = FALSE, verbose = FALSE)
    
    fit
  }
  
  #'  Function tying this together by extracting posteriors and refitting the model
  fit_aux_claim <- function(i, iterations, og_fit, nSites, nYear, model_name,
                            n.chains, n.adapt, n.burnin, n.iter, n.thin) { #spp, covariate_array, 
    
    #'  Grab details for focal iteration
    iter_deets <- iterations[[i]]
    spp <- iter_deets$spp
    covariate_array <- iter_deets$covariate_array
    
    #'  Grab the posterior samples from the specified variable in the ind. claim
    post <- extract_latent_post_summaries(og_fit, spp, nSites, nYear)
    #'  Grab the spp.latent posterior mean and the "observed" predictor 
    y_vec <- as.vector(post$mean)
    x_vec <- as.vector(covariate_array)
    keep <- !is.na(y_vec) & !is.na(x_vec)
    #print(sum(keep)) # Better not be 0 or close to 0 (means very little data going into model)
    #'  Refit model with added independence claim using spp.latent posterior mean 
    #'  and specified predictor as y and x
    mod_out <- fit_one_dSep_claim(y = y_vec[keep], x = x_vec[keep], z = NULL, 
                                  n.chains = n.chains, n.adapt = n.adapt, n.burnin = n.burnin, 
                                  n.iter = n.iter, n.thin = n.thin, model_name, iter_num = i)
    
    #'  Flag convergence issues
    max_rhat <- suppressWarnings(max(unlist(mod_out$Rhat), na.rm = TRUE))
    if(is.finite(max_rhat) && max_rhat > 1.1) {
      warning(sprintf("Iteration %d: max Rhat = %.3f -- check convergence", i, max_rhat))
    }
    
    #'  Temporary directory to save JAGS outputs
    out_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1", model_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Save JAGS output for each iteration
    jags_out <- file.path(out_dir, sprintf("iter_%03d.rds", i))
    saveRDS(list(fit = mod_out, b_x = as.numeric(mod_out$sims.list$b_x), 
                 config = iter_deets, max_rhat = max_rhat), jags_out)
  }
  
  
  #'  ------------------------------------------------------------
  #####  Functions for ind. claims with only exogenous variables  #####
  #'  ------------------------------------------------------------
  #'  Function to grab covariate data and call fit_one_dSep_claim()
  fit_covariate_claim <- function(i, iterations, model_name, #y_array, x_array, iter_num,
                                  n.chains, n.adapt, n.burnin, n.iter, n.thin) {
    iter_deets <- iterations[[i]]
    y_vec <- as.vector(iter_deets$y_array)
    x_vec <- as.vector(iter_deets$x_array)
    keep <- !is.na(y_vec) & !is.na(x_vec)
  
    #'  Use simplified JAGS code and regression in fit_one_dSep_claim() to test claim
    mod_out <- fit_one_dSep_claim(y = y_vec[keep], x = x_vec[keep], z = NULL,
                       n.chains = n.chains, n.adapt = n.adapt, n.burnin = n.burnin,
                       n.iter = n.iter, n.thin = n.thin, model_name = model_name,
                       iter_num = i)
    
    #'  Temporary directory to save JAGS outputs
    out_dir <- file.path("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1", model_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    #'  Save JAGS output for each iteration
    jags_out <- file.path(out_dir, sprintf("iter_%03d.rds", i))
    saveRDS(list(fit = mod_out, b_x = as.numeric(mod_out$sims.list$b_x), 
                 config = iter_deets), jags_out)
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
  
  ### MAKE SURE SEM_TOPDOWN IS IN WORKING DIRECTORY  ###
  
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
  
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Model registry that defines the original regressions in SEM to be updated
  #'  with each iteration of d-Sep testing
  sem_registry <- list(
    #'  Regression 1: lion.latent
    list(covs = c("wtd.latent"), spp = c(".wtd"), indices = as.integer(c(1))),
    #'  Regression 2: wolf.latent
    list(covs = c("wolf.latent", "elk.latent", "moose.latent"), spp = c(".wolf", ".elk", ".moose"), indices = as.integer(c(1,1,1))),
    #'  Regression 3: bear.latent
    list(covs = c("bear.latent", "forest", "wolf.latent"), spp = c(".bear", ".forest", ".wolf"), indices = as.integer(c(1,1,1))),
    #'  Regression 4: coy.latent
    list(covs = c("coy.latent", "wtd.latent", "wolf.latent"), spp = c(".coy", ".wtd", ".wolf"), indices = as.integer(c(1,1,1))),
    #'  Regression 5: elk.latent
    list(covs = c("elk.latent", "forest", "wsi"), spp = c(".elk", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regression 6: moose.latent
    list(covs = c("moose.latent", "forest", "wsi"), spp = c(".moose", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    #'  Regression 7: wtd.latent
    list(covs = c("wtd.latent", "forest", "wsi"), spp = c(".wtd", ".forest", ".wsi"), indices = as.integer(c(1,1,1))),
    )
  #'  Source d-Sep custom regressions for iterative d-separation tests 
  source("./Scripts/Structural_Equation_Models/d_Sep_active_regressions_bottomup_inter.R")
  
  #'  Bundle data and draw inits using functions in in Format_RNmodel_Posteriors_for_SEM.R
  data_JAGS_bundle_bottomup_inter <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
                                                dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s,
                                                covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
                                                covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
                                                nwolf = 4, nlion = 0, nbear = 2, ncoy = 2, nelk = 3, 
                                                nmoose = 3, nwtd = 4, nharv = 0, nfor = 4, nwsi = 3)
  num.chains <- 3
  initsList_bottomup_inter <- vector('list', num.chains) 
  for(i in 1:num.chains) {
    initsList_bottomup_inter[[i]] <- generate_inits(nwolf = 4, nlion = 0, nbear = 2, ncoy = 2, nelk = 3, nmoose = 3, 
                                             nwtd = 4, nharv = 0, nfor = 4, nwsi = 3, nSpp = 7, nSites = 23, nYear = 4)
  }
  
  start.time = Sys.time()
  #'  Fit and save model iterations
  saved_paths <- future_lapply(
    seq_along(dSep_iterations_bottomup_int),
    function(i) run_dSep_iterations(i, iterations = dSep_iterations_bottomup_int, template = model_template, registry = sem_registry,
                                    data_bundle = data_JAGS_bundle_bottomup_inter, listInits = initsList_bottomup_inter, model_name = "BottomUp_Interference"),
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
  
  #'  -------------------------------------
  ####  Calculate p-values and Fisher's C  ####
  #'  -------------------------------------
  #'  -------------------------------------------
  #####  Top-down interference model iterations  #####
  #'  -------------------------------------------
  #'  Load all iterations of the JAGS model
  all_results_topdown_inter <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/TopDown_Interference", full.names = TRUE), readRDS)
  
  #'  Rename data bundle
  data_JAGS_bundle <- data_JAGS_bundle_topdown_inter
  #'  Create list of "observed" values of focal response variable, one per d-Sep test               # instances where x was used as y in d-Sep test noted below
  y_list <- list(data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, 
                 data_JAGS_bundle$wtd.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$coy.hat,      
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,     
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,       
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat,
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat,
                 data_JAGS_bundle$elk.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$wolf.hat,
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat,     
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$lion.hat, 
                 data_JAGS_bundle$moose.hat, data_JAGS_bundle$elk.hat, data_JAGS_bundle$coy.hat,      
                 data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat, data_JAGS_bundle$elk.hat, 
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat, data_JAGS_bundle$wolf.hat,      
                 data_JAGS_bundle$lion.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat,      
                 data_JAGS_bundle$wolf.hat, data_JAGS_bundle$lion.hat, data_JAGS_bundle$bear.hat,
                 data_JAGS_bundle$coy.hat, data_JAGS_bundle$coy.hat, data_JAGS_bundle$bear.hat,
                 data_JAGS_bundle$bear.hat)                                                           
  #'  Leaves you with 67 d-Sep tests that were possible given the constructs of space and time and our data
  
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
                                  mod_out[[19]]$beta.lion[,1], mod_out[[20]]$beta.lion[,1], mod_out[[21]]$beta.lion[,1],
                                  mod_out[[22]]$beta.lion[,1], mod_out[[23]]$beta.wtd[,2], mod_out[[24]]$beta.wtd[,2],            # note the different indexing
                                  mod_out[[25]]$beta.wtd[,2], mod_out[[26]]$beta.wtd[,2], mod_out[[27]]$beta.wtd[,2],             # note the different indexing
                                  mod_out[[28]]$beta.wtd[,2], mod_out[[29]]$beta.coy[,1], mod_out[[30]]$beta.coy[,1],             # note the different indexing
                                  mod_out[[31]]$beta.coy[,1], mod_out[[32]]$beta.coy[,1], mod_out[[33]]$beta.coy[,1],
                                  mod_out[[34]]$beta.harvest[,1], mod_out[[35]]$beta.harvest[,1], mod_out[[36]]$beta.harvest[,1],
                                  mod_out[[37]]$beta.harvest[,2], mod_out[[38]]$beta.harvest[,2], mod_out[[39]]$beta.bear[,1],    # note the different indexing
                                  mod_out[[40]]$beta.bear[,1], mod_out[[41]]$beta.bear[,1], mod_out[[42]]$beta.bear[,1], 
                                  mod_out[[43]]$beta.bear[,1], mod_out[[44]]$beta.harvest[,1], mod_out[[45]]$beta.harvest[,1], 
                                  mod_out[[46]]$beta.harvest[,1], mod_out[[47]]$beta.harvest[,2], mod_out[[48]]$beta.harvest[,2], # note the different indexing
                                  mod_out[[49]]$beta.harvest[,1], mod_out[[50]]$beta.harvest[,1], mod_out[[51]]$beta.harvest[,1], 
                                  mod_out[[52]]$beta.harvest[,2], mod_out[[53]]$beta.harvest[,2], mod_out[[54]]$beta.moose[,2],   # note the different indexing
                                  mod_out[[55]]$beta.moose[,2], mod_out[[56]]$beta.moose[,2], mod_out[[57]]$beta.moose[,2],       # note the different indexing
                                  mod_out[[58]]$beta.moose[,2], mod_out[[59]]$beta.elk[,2], mod_out[[60]]$beta.elk[,2],           # note the different indexing
                                  mod_out[[61]]$beta.elk[,2], mod_out[[62]]$beta.elk[,2], mod_out[[63]]$beta.coy[,2],             # note the different indexing
                                  mod_out[[64]]$beta.coy[,2], mod_out[[65]]$beta.coy[,2], mod_out[[66]]$beta.bear[,2],            # note the different indexing
                                  mod_out[[67]]$beta.bear[,2])                                                                    # note the different indexing
  
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
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_topdown_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_topdown_inter_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_topdown_inter_df <- stack(p.rope_topdown_inter_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4))
  
  #'  ----------------------
  ######  Bayesian p-value  ######
  #'  ----------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_topdown_inter_list <- mapply(bayes_p_iterations, post_beta = post_list_topdown_inter, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_topdown_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_topdown_inter_list)[i] <- list_name
  }
  bayes.p_topdown_inter_df <- stack(bayes.p_topdown_inter_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4))
  
  #'  Join both d-Sep test p-values and save
  p.val_topdown_inter_df <- full_join(p.rope_topdown_inter_df, bayes.p_topdown_inter_df, by = "iteration")
  
  write_csv(p.val_topdown_inter_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_topdown_inter.csv")
  
  #'  ----------------
  ######  Fisher's C  ######
  #'  ----------------
  fishers.C_topdown_inter <- fishers_C(pval = bayes.p_topdown_inter_df$bayes.p, n_iter = nrow(bayes.p_topdown_inter_df))
  print(fishers.C_topdown_inter)
  
  
  #'  ------------------------------
  #####  Top-down model iterations  #####
  #'  ------------------------------
  #'  Load iterations of the JAGS model
  #'  Note: this list of outputs is based on the number of independence claims 
  #'  assessed using the d_Sep_active_regression_topdown.R active regression list. 
  #'  This is not a complete list of all independence claims being tested.
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
                  data_JAGS_bundle$wolf.hat, data_JAGS_bundle$nWolf)
  
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
  # post_list_topdown_tmin1 <- list(mod_out2[[1]]$beta.harvest[,1], mod_out2[[2]]$beta.harvest[,1],  mod_out2[[3]]$beta.harvest[,1],
  #                                 mod_out2[[4]]$beta.harvest[,1],  mod_out2[[5]]$beta.harvest[,1],  mod_out2[[6]]$beta.harvest[,1],
  #                                 mod_out2[[7]]$beta.harvest[,1],  mod_out2[[8]]$beta.wtd[,1],  mod_out2[[9]]$beta.harvest[,1],
  #                                 mod_out2[[10]]$beta.wtd[,1],  mod_out2[[11]]$beta.wtd[,1],  mod_out2[[12]]$beta.wtd[,1], 
  #                                 mod_out2[[13]]$beta.harvest[,1],  mod_out2[[14]]$beta.wtd[,1],  mod_out2[[15]]$beta.harvest[,1], 
  #                                 mod_out2[[16]]$beta.wtd[,1],  mod_out2[[17]]$beta.harvest[,1],  mod_out2[[18]]$beta.harvest[,1],
  #                                 mod_out2[[19]]$beta.moose[,1],  mod_out2[[20]]$beta.moose[,1],  mod_out2[[21]]$beta.moose[,1], 
  #                                 mod_out2[[22]]$beta.harvest[,1],  mod_out2[[23]]$beta.moose[,1],  mod_out2[[24]]$beta.harvest[,1],
  #                                 mod_out2[[25]]$beta.moose[,1],  mod_out2[[26]]$beta.harvest[,1],  mod_out2[[27]]$beta.harvest[,1],
  #                                 mod_out2[[28]]$beta.harvest[,1],  mod_out2[[29]]$beta.harvest[,1],  mod_out2[[30]]$beta.harvest[,1],
  #                                 mod_out2[[31]]$beta.harvest[,1],  mod_out2[[32]]$beta.lion[,1],  mod_out2[[33]]$beta.lion[,1],
  #                                 mod_out2[[34]]$beta.harvest[,1],  mod_out2[[35]]$beta.lion[,1],  mod_out2[[36]]$beta.harvest[,1], 
  #                                 mod_out2[[37]]$beta.lion[,1],  mod_out2[[38]]$beta.harvest[,1],  mod_out2[[39]]$beta.elk[,1], 
  #                                 mod_out2[[40]]$beta.harvest[,1],  mod_out2[[41]]$beta.elk[,1],  mod_out2[[42]]$beta.harvest[,1],
  #                                 mod_out2[[43]]$beta.elk[,1],  mod_out2[[44]]$beta.harvest[,1],  mod_out2[[45]]$beta.harvest[,1], 
  #                                 mod_out2[[46]]$beta.coy[,1],  mod_out2[[47]]$beta.harvest[,1],  mod_out2[[48]]$beta.coy[,1], 
  #                                 mod_out2[[49]]$beta.harvest[,1], mod_out2[[50]]$beta.harvest[,1], mod_out2[[51]]$beta.harvest[,1], 
  #                                 mod_out2[[52]]$beta.harvest[,1], mod_out2[[53]]$beta.bear[,1], mod_out2[[54]]$beta.harvest[,1],
  #                                 mod_out2[[55]]$beta.harvest[,1], mod_out2[[56]]$beta.harvest[,1])
  
  #'  Load more iterations of the JAGS model (this time assessing correlation between exogenous variables)
  all_results_topdown_exog <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/tmin1/TopDown_exog", full.names = TRUE), readRDS)
  y_list3 <- list(data_JAGS_bundle$elkHarv, data_JAGS_bundle$bearHarv, data_JAGS_bundle$wolfHarv,
                  data_JAGS_bundle$lionHarv, data_JAGS_bundle$bearHarv, data_JAGS_bundle$wolfHarv,
                  data_JAGS_bundle$lionHarv, data_JAGS_bundle$wolfHarv, data_JAGS_bundle$lionHarv,
                  data_JAGS_bundle$wolfHarv) ######### UPDATE THESE
  mod_out3 <- list()
  for(i in 1:length(all_results_topdown_exog)) {
    mod_out3[[i]] <- all_results_topdown_exog[[i]]$fit$sims.list
  }
  post_list_topdown_exog <- list(mod_out3[[1]]$b_x, mod_out3[[2]]$b_x, mod_out3[[3]]$b_x,
                                  mod_out3[[4]]$b_x, mod_out3[[5]]$b_x, mod_out3[[6]]$b_x,
                                  mod_out3[[7]]$b_x, mod_out3[[8]]$b_x, mod_out3[[9]]$b_x,
                                  mod_out3[[10]]$b_x)
  
  
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
  fishers.C_topdown <- fishers_C(pval = bayes.p_topdown_df$bayes.p, n_iter = nrow(bayes.p_topdown_df))
  print(fishers.C_topdown)
  
  
  #'  --------------------------------------------
  #####  Bottom-up interference model iterations  #####
  #'  --------------------------------------------
  #'  Load all iterations of the JAGS model
  all_results_bottomup_inter <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/BottomUp_Interference", full.names = TRUE), readRDS)
  
  #'  Rename data bundle
  data_JAGS_bundle <- data_JAGS_bundle_bottomup_inter
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
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_bottomup_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_inter_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_bottomup_inter_df <- stack(p.rope_bottomup_inter_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4))
  
  #'  ----------------------
  ######  Bayesian p-value  ######
  #'  ----------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_bottomup_inter_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup_inter, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_bottomup_inter_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_inter_list)[i] <- list_name
  }
  bayes.p_bottomup_inter_df <- stack(bayes.p_bottomup_inter_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4))
  
  #'  Join both d-Sep test p-values and save
  p.val_bottomup_inter_df <- full_join(p.rope_bottomup_inter_df, bayes.p_bottomup_inter_df, by = "iteration")
  
  write_csv(p.val_bottomup_inter_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_bottomup_inter.csv")
  
  #'  ----------------
  ######  Fisher's C  ######
  #'  ----------------
  fishers.C_bottomup_inter <- fishers_C(pval = bayes.p_bottomup_inter_df$bayes.p, n_iter = nrow(bayes.p_bottomup_inter_df))
  print(fishers.C_bottomup_inter)
  
  
  #'  -------------------------------
  #####  Bottom-up model iterations  #####
  #'  -------------------------------
  #'  Load all iterations of the JAGS model
  all_results_bottomup <- lapply(list.files("./Outputs/SEM/JAGS_out/d_Sep/Results/BottomUp", full.names = TRUE), readRDS)
  
  #'  Rename data bundle
  data_JAGS_bundle <- data_JAGS_bundle_bottomup
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
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(p.rope_bottomup_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(p.rope_bottomup_list)[i] <- list_name
  }
  
  #'  Convert list to a data frame
  p.rope_bottomup_df <- stack(p.rope_bottomup_list) %>%
    transmute(iteration = ind,
              p.rope = round(values, 4))
  
  #'  ----------------------
  ######  Bayesian p-value  ######
  #'  ----------------------
  bayes_p_iterations <- function(post_beta) {
    bayes.p_val <- bayes_pvalue(post = post_beta)
    print(bayes.p_val)
    return(bayes.p_val)
  }
  bayes.p_bottomup_list <- mapply(bayes_p_iterations, post_beta = post_list_bottomup, SIMPLIFY = FALSE)
  
  #'  Rename objects in the list based on iteration 
  for(i in 1:length(bayes.p_bottomup_list)) {
    list_name <- sprintf("regression.%03d", i)
    names(bayes.p_bottomup_list)[i] <- list_name
  }
  bayes.p_bottomup_df <- stack(bayes.p_bottomup_list) %>%
    transmute(iteration = ind,
              bayes.p = round(values, 4))
  
  #'  Join both d-Sep test p-values and save
  p.val_bottomup_df <- full_join(p.rope_bottomup_df, bayes.p_bottomup_df, by = "iteration")
  
  write_csv(p.val_bottomup_df, "./Outputs/SEM/JAGS_out/d_Sep/p_val_bottomup.csv")
  
  #'  ----------------
  ######  Fisher's C  ######
  #'  ----------------
  fishers.C_bottomup <- fishers_C(pval = bayes.p_bottomup_df$bayes.p, n_iter = nrow(bayes.p_bottomup_df))
  print(fishers.C_bottomup)
  
