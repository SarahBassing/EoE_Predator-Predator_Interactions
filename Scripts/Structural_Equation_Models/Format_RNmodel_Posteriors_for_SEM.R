  #'  ----------------------------------------------------
  #'  Format Royle-Nichols model results for Bayesian SEM
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  January 2026
  #'  ----------------------------------------------------
  #'  Load posterior distributions from RN models and format for as input for
  #'  Bayesian structural equation models
  #'  This script is sourced by Bayesian_SEMs_relative_density_index_1yLag.R
  #'  ----------------------------------------------------
  
  #' #'  Clean workspace
  #' rm(list = ls())
  
  #'  Load libraries
  library(tidyverse)
  library(ggplot2)
  
  #'  Load JAGS outputs for species and cluster-specific posterior distributions
  file_names <- list.files(path = "./Outputs/Relative_Abundance/RN_model/JAGS_out/", 
                           pattern = "\\.RData$", full.names = T)
  for(i in 1:length(file_names)) {
    load(file_names[i])
  }
  
  #'  List model outputs per year
  mods_yr1 <- list(RN_wolf_20s, RN_lion_20s, RN_bear_20s, RN_coy_20s, RN_elk_20s, RN_moose_20s, RN_wtd_20s)
  mods_yr2 <- list(RN_wolf_21s, RN_lion_21s, RN_bear_21s, RN_coy_21s, RN_elk_21s, RN_moose_21s, RN_wtd_21s)
  mods_yr3 <- list(RN_wolf_22s, RN_lion_22s, RN_bear_22s, RN_coy_22s, RN_elk_22s, RN_moose_22s, RN_wtd_22s)
  mods_yr4 <- list(RN_wolf_23s, RN_lion_23s, RN_bear_23s, RN_coy_23s, RN_elk_23s, RN_moose_23s, RN_wtd_23s)
  
  #'  ------------------------------------------
  ####  Format posterior means and SDs for SEM  ####
  #'  ------------------------------------------
  #'  Function to save posterior means and SD for each site per species, cluster, and year
  cluster_posterior_summary <- function(mod_post, cluster_out, cluster_integer) {
    #' #'  Empty lists to hold model outputs of interest
    #' posterior_mu <- posterior_sd <- c()
    #' #'  Loop through posteriors and save outputs of interest that match cluster names
    #' for(i in 1:length(cluster_out)) {
    #'   posterior_mu[i] <- unlist(mod_post$mean[names(mod_post$mean) %in% cluster_out[i]])
    #'   posterior_sd[i] <- unlist(mod_post$sd[names(mod_post$sd) %in% cluster_out[i]])
    #' }
    #' 
    #' #'  Bind mu and sd outputs into a matrix
    #' cluster_posteriors <- as.data.frame(cbind(posterior_mu, posterior_sd)) %>%
    #'   cbind(cluster_out) %>%
    #'   rename(cluster = cluster_out)
    
    cluster_posteriors <- data.frame(
      posterior_mu = unlist(mod_post$mean[cluster_out]),
      posterior_sd = unlist(mod_post$sd[cluster_out]),
      cluster_id = cluster_out,
      cluster = cluster_integer,
      row.names = NULL
    )
      
    return(cluster_posteriors)
  }
  #'  Names of estimated posteriors of interest per year
  rdi.cl_2020 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", #"rdi.cl9",
                   "rdi.cl19", "rdi.cl20","rdi.cl21", "rdi.cl22", "rdi.cl23") # , "rdi.cl24"
  rdi.cl_2021.2023 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9", "rdi.cl10",
                        "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", "rdi.cl15", "rdi.cl16", "rdi.cl17", "rdi.cl18", "rdi.cl19", "rdi.cl20",
                        "rdi.cl21", "rdi.cl22", "rdi.cl23") #, "rdi.cl24"
  rdi.cl_2020_int <- c(1, 2, 3, 4, 5, 6, 7, 8, 19, 20, 21, 22, 23) 
  rdi.cl_2021.2023_int <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23) 
  
  #'  Run function for each year and species
  posteriors_20s <- lapply(mods_yr1, cluster_posterior_summary, cluster_out = rdi.cl_2020, cluster_integer = rdi.cl_2020_int)
  posteriors_21s <- lapply(mods_yr2, cluster_posterior_summary, cluster_out = rdi.cl_2021.2023, cluster_integer = rdi.cl_2021.2023_int)
  posteriors_22s <- lapply(mods_yr3, cluster_posterior_summary, cluster_out = rdi.cl_2021.2023, cluster_integer = rdi.cl_2021.2023_int)
  posteriors_23s <- lapply(mods_yr4, cluster_posterior_summary, cluster_out = rdi.cl_2021.2023, cluster_integer = rdi.cl_2021.2023_int)
  
  #'  Name lists by species
  spp <- c("wolf", "mountain_lion", "bear_black", "coyote", "elk", "moose", "whitetailed_deer")
  names(posteriors_20s) <- spp
  names(posteriors_21s) <- spp
  names(posteriors_22s) <- spp
  names(posteriors_23s) <- spp
  
  #'  Double check these look right
  #'  Mean and SD matrix vs cluster specific summary outputs
  print(posteriors_20s[[1]]) # wolf (2020)
  RN_wolf_20s$mean[10:11]; RN_wolf_20s$sd[10:11]  # 1st two clusters
  print(posteriors_21s[[3]]) # bear (2021)
  RN_bear_21s$mean[10:11]; RN_bear_21s$sd[10:11]  # 1st two clusters
  print(posteriors_22s[[5]]) # elk (2022)
  RN_elk_22s$mean[31:32]; RN_elk_22s$sd[31:32]  # last two clusters
  print(posteriors_23s[[7]]) # white-tailed deer (2023)
  RN_wtd_23s$mean[31:32]; RN_wtd_23s$sd[31:32]  # last two clusters
  
  #' #'  Create new posteriors for 2021 that exclude GMU 1 RDI posteriors 
  #' #'  Needed below when stacking data for time t vs t-1
  #' drop_gmu1 <- function(post_rdi, gmu) {  
  #'   #'  Filter out GMU1 clusters
  #'   noGMU1 <- post_rdi %>%
  #'     filter(!cluster %in% gmu)
  #' 
  #'   return(noGMU1)
  #' }
  #' gmu1 <- c("rdi.cl9", "rdi.cl10", "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", 
  #'           "rdi.cl15", "rdi.cl16", "rdi.cl17", "rdi.cl18")
  #' posteriors_21s_noGMU1 <- lapply(posteriors_21s, drop_gmu1, gmu = gmu1) 
  #' 
  # posteriors <- list(posteriors_20s, posteriors_21s, posteriors_22s, posteriors_23s)
  # names(posteriors) <- c("post_2020", "post_2021", "post_2022", "post_2023")
  
  #'  -------------------------
  ####  Format covariate data  ####
  #'  -------------------------
  #'  All covariates generated in Format_spatial_covariates_for_SMEs.R script
  
  #'  Rename CLuster_unique to cluster
  cluster_poly_covs_df <- rename(cluster_poly_covs_df, cluster = Cluster_unique)
  
  #'  Split out covariate data that do not need to be time-lagged (e.g., wolf_2020 is applied to RDI 2020)
  wolf_harvest <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, wolfharvest_2020_per100km, wolfharvest_2021_per100km, 
                  wolfharvest_2022_per100km, wolfharvest_2023_per100km)
  wsi <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, WSI_20_gmu, WSI_21_gmu, WSI_22_gmu, WSI_23_gmu)
  forest <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, distForest20, distForest21, distForest22, distForest23)
  roads <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, roaddens)
  public <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, proppub)
  
  #'  Split out covariate data that do need to be time lagged (e.g., bear_2019 is applied to RDI 2020)
  bear <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, Bear_per100km_2019, Bear_per100km_2020, Bear_per100km_2021, Bear_per100km_2022)
  lion <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, Lion_per100km_2019, Lion_per100km_2020, Lion_per100km_2021, Lion_per100km_2022)
  elk <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, Elk_per100km_2019, Elk_per100km_2020, Elk_per100km_2021, Elk_per100km_2022)
  moose <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, Moose_per100km_2019, Moose_per100km_2020, Moose_per100km_2021, Moose_per100km_2022)
  deer <- cluster_poly_covs_df %>%
    dplyr::select(cluster, GMU, Deer_per100km_2019, Deer_per100km_2020, Deer_per100km_2021, Deer_per100km_2022)
  
 covs_2020 <- data.frame(cluster = wolf_harvest$cluster, wolfHarvest = wolf_harvest$wolfharvest_2020_per100km,
                         bearHarvest = bear$Bear_per100km_2019, lionHarvest = lion$Lion_per100km_2019, 
                         elkHarvest = elk$Elk_per100km_2019, mooseHarvest = moose$Moose_per100km_2019, 
                         deerHarvest = deer$Deer_per100km_2019, forest = forest$distForest20,
                         wsi = wsi$WSI_20_gmu)
 covs_2021 <- data.frame(cluster = wolf_harvest$cluster, wolfHarvest = wolf_harvest$wolfharvest_2021_per100km, 
                         bearHarvest = bear$Bear_per100km_2020, lionHarvest = lion$Lion_per100km_2020, 
                         elkHarvest = elk$Elk_per100km_2020, mooseHarvest = moose$Moose_per100km_2020, 
                         deerHarvest = deer$Deer_per100km_2020, forest = forest$distForest21,
                         wsi = wsi$WSI_21_gmu)
 covs_2022 <- data.frame(cluster = wolf_harvest$cluster, wolfHarvest = wolf_harvest$wolfharvest_2021_per100km, 
                         bearHarvest = bear$Bear_per100km_2021, lionHarvest = lion$Lion_per100km_2021, 
                         elkHarvest = elk$Elk_per100km_2021, mooseHarvest = moose$Moose_per100km_2021, 
                         deerHarvest = deer$Deer_per100km_2021, forest = forest$distForest22, 
                         wsi = wsi$WSI_22_gmu)
 covs_2023 <- data.frame(cluster = wolf_harvest$cluster, wolfHarvest = wolf_harvest$wolfharvest_2021_per100km, 
                         bearHarvest = bear$Bear_per100km_2022, lionHarvest = lion$Lion_per100km_2022, 
                         elkHarvest = elk$Elk_per100km_2022, mooseHarvest = moose$Moose_per100km_2022, 
                         deerHarvest = deer$Deer_per100km_2022, forest = forest$distForest23, 
                         wsi = wsi$WSI_23_gmu)
 
 # covs <- list(covs_2020, covs_2021, covs_2022, covs_2023)
  
  #'  Function to reformat data for JAGS so indexed by cluster and year
  reshape_to_site_year <- function(year_list, value_col, nSites, nYear) {
    mat <- matrix(NA_real_, nrow = nSites, ncol = nYear)
    for (y in seq_len(nYear)) {
      df <- year_list[[y]]
      
      #'  Fail with a clear, specific message rather than base R's
      #'  cryptic "subscript out of bounds" if the expected column is
      #'  missing -- tells you exactly which column/year is the problem
      #'  and what columns ARE actually available, so you don't have to
      #'  hunt for it separately.
      if (!value_col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found (year index %d). Available columns: %s",
          value_col, y, paste(names(df), collapse = ", ")
        ))
      }
      if (!"cluster" %in% names(df)) {
        stop(sprintf(
          "Column 'cluster' not found (year index %d). Available columns: %s",
          y, paste(names(df), collapse = ", ")
        ))
      }
      
      site_idx <- as.integer(as.character(df$cluster))
      mat[site_idx, y] <- df[[value_col]]
    }
    mat
  }
  
  bundle_dat <- function(dat_yr1, dat_yr2, dat_yr3, dat_yr4, covs_yr1, covs_yr2, 
                         covs_yr3, covs_yr4, nwolf, nlion, nbear, ncoy, nelk, 
                         nmoose, nwtd, nharv, nfor, nwsi, nSites = 23, nYear = 4) {
    
    #'  --------------------------------------------------------------
    ####  Reshape each species' posterior mu/sd into [nSites, nYear]  ####
    #'  --------------------------------------------------------------
    spp_names <- c("wolf", "mountain_lion", "bear_black", "coyote", "elk", "moose", "whitetailed_deer")
    
    hat_list <- list()
    sigma_hat_list <- list()
    
    for (sp in spp_names) {
      #'  Create list of species specific data across all years of sampling
      yr_list <- list(dat_yr1[[sp]], dat_yr2[[sp]], dat_yr3[[sp]], dat_yr4[[sp]])  
      
      hat_mat   <- reshape_to_site_year(yr_list, "posterior_mu", nSites, nYear)
      sigma_mat <- reshape_to_site_year(yr_list, "posterior_sd", nSites, nYear)
      
      #'  Standardize (z-score) the posterior mean, pooling across all clusters
      #'  and years per species. SD is rescaled by the same scaling factor but 
      #'  NOT re-centered, since a standard deviation describes spread, not location.
      #'  This preserves variability and magnitude of changes in density across time. 
      pooled_mean <- mean(hat_mat, na.rm = TRUE)
      pooled_sd   <- sd(hat_mat, na.rm = TRUE)
      
      hat_mat   <- (hat_mat - pooled_mean) / pooled_sd
      sigma_mat <- sigma_mat / pooled_sd
      
      #'  Fill missing data with median observed sigma
      #'  GMU1 not sampled in year1, resulting in NAs for those clusters. But
      #'  cannot have NAs in sigma_hat in JAGS. Median observed sigma is a
      #'  placeholder here; it is inconsequential, since it only sets the spread
      #'  of an imputed .hat draw that has no downstream data to inform.
      sigma_mat[is.na(sigma_mat)] <- median(sigma_mat, na.rm = TRUE)
      
      hat_list[[sp]] <- hat_mat
      sigma_hat_list[[sp]] <- sigma_mat
    }
  
  #'  -------------------------------------------
  ####  Reshape covariates into [nSites, nYear]  ####
  #'  -------------------------------------------
  cov_years    <- list(covs_yr1, covs_yr2, covs_yr3, covs_yr4)
  
  wsi_mat    <- reshape_to_site_year(cov_years, "wsi", nSites, nYear)
  forest_mat <- reshape_to_site_year(cov_years, "forest", nSites, nYear)
  wolf_mat    <- reshape_to_site_year(cov_years, "wolfHarvest", nSites, nYear)
  lion_mat <- reshape_to_site_year(cov_years, "lionHarvest", nSites, nYear)
  bear_mat    <- reshape_to_site_year(cov_years, "bearHarvest", nSites, nYear)
  elk_mat <- reshape_to_site_year(cov_years, "elkHarvest", nSites, nYear)
  moose_mat    <- reshape_to_site_year(cov_years, "mooseHarvest", nSites, nYear)
  deer_mat <- reshape_to_site_year(cov_years, "deerHarvest", nSites, nYear)
  
  #'  -------------------------
  ####  Assemble final bundle  ####
  #'  -------------------------
  bundled <- list(
    nWolf = nwolf, nLion = nlion, nBear = nbear, nCoy = ncoy,
    nElk = nelk, nMoose = nmoose, nWtd = nwtd,
    nharvest = nharv, nforest = nfor, nWSI = nwsi,
    nSites = nSites, nYear = nYear, nSpp = 7,
    
    #'  Posterior means/SDs, now [nSites, nYear]
    wolf.hat = hat_list$wolf,             wolf.sigma_hat = sigma_hat_list$wolf,
    lion.hat = hat_list$mountain_lion,    lion.sigma_hat = sigma_hat_list$mountain_lion,
    bear.hat = hat_list$bear_black,       bear.sigma_hat = sigma_hat_list$bear_black,
    coy.hat  = hat_list$coy,              coy.sigma_hat  = sigma_hat_list$coy,
    elk.hat  = hat_list$elk,              elk.sigma_hat  = sigma_hat_list$elk,
    moose.hat = hat_list$moose,           moose.sigma_hat = sigma_hat_list$moose,
    wtd.hat  = hat_list$whitetailed_deer, wtd.sigma_hat  = sigma_hat_list$whitetailed_deer,
    
    #'  Covariates, now [nSites, nYear]
    wsi = wsi_mat,
    forest = forest_mat,
    wolfHarv  = wolf_mat,
    lionHarv  = lion_mat,
    bearHarv  = bear_mat,
    elkHarv   = elk_mat,
    mooseHarv = moose_mat,
    deerHarv  = deer_mat
  )
  
  str(bundled)
  return(bundled)
  
  }
  
  # data_JAGS_bundle_tst <- bundle_dat(dat_yr1 = posteriors_20s, dat_yr2 = posteriors_21s, 
  #                                                    dat_yr3 = posteriors_22s, dat_yr4 = posteriors_23s, 
  #                                                    covs_yr1 = covs_2020, covs_yr2 = covs_2021, 
  #                                                    covs_yr3 = covs_2022, covs_yr4 = covs_2023, 
  #                                                    nwolf = 4, nlion = 1, nbear = 1, ncoy = 1, nelk = 2, 
  #                                                    nmoose = 2, nwtd = 3, nharv = 0, nfor = 4, nwsi = 3)
  
  #'  Function to generate inits for intercepts, slopes and latent-state array
  generate_inits <- function(nwolf, nlion, nbear, ncoy, nelk, nmoose, nwtd, nharv, 
                             nfor, nwsi, nSpp, nSites, nYear) {
    beta_inits <- list(beta.int = runif(nSpp, -0.5, 0.5),
                       beta.int.tmin1 = runif(nSpp, -0.5, 0.5),
                       beta.wolf = runif(nwolf, -0.5, 0.5),
                       beta.lion = runif(nlion, -0.5, 0.5),
                       beta.bear = runif(nbear, -0.5, 0.5),
                       beta.coy = runif(ncoy, -0.5, 0.5),
                       beta.elk = runif(nelk, -0.5, 0.5),
                       beta.moose = runif(nmoose, -0.5, 0.5),
                       beta.wtd = runif(nwtd, -0.5, 0.5), 
                       beta.harvest = runif(nharv, -0.5, 0.5),
                       beta.forest = runif(nfor, -0.5, 0.5),
                       beta.wsi = runif(nwsi, -0.5, 0.5))
    
    #'  Inits for latent-state array [nSites, nYear], one per species
    #'  Important because NAs in GMU1 year 1 data so JAGS needs something to start
    #'  with since its missing data to inform how to start
    latent_inits <- list(wolf.latent = matrix(runif(nSites * nYear, -0.5, 0.5), nrow = nSites, ncol = nYear), 
                         lion.latent = matrix(runif(nSites * nYear, -0.5, 0.5), nrow = nSites, ncol = nYear),
                         bear.latent = matrix(runif(nSites * nYear, -0.5, 0.5), nrow = nSites, ncol = nYear),
                         coy.latent = matrix(runif(nSites * nYear, -0.5, 0.5), nrow = nSites, ncol = nYear),
                         elk.latent = matrix(runif(nSites * nYear, -0.5, 0.5), nrow = nSites, ncol = nYear),
                         moose.latent = matrix(runif(nSites * nYear, -0.5, 0.5), nrow = nSites, ncol = nYear),
                         wtd.latent = matrix(runif(nSites * nYear, -0.5, 0.5), nrow = nSites, ncol = nYear))
    
    #'  Contain all inits togther
    c(beta_inits, latent_inits)
  }
  
  
  
  
  #' #'  Stack species-specific estimated posterior summaries across iterations 
  #' #'  (e.g., time t-1 = wolf 2020, 2021, 2022 affects time t = wolf 2021, 2022, 2023)
  #' #'  yr1 2020 --> yr2 2021
  #' #'  yr2 2021 --> yr3 2022
  #' #'  yr3 2022 --> yr4 2023
  #' stacked_posteriors <- function(yr1, yr2, yr3, yr4, yr2_noGMU1) {
  #'   #'  Posteriors from time t (yr2, yr3, yr4)
  #'   #'  NOTE: yr2 for time_t excludes GMU 1 posteriors because there were no
  #'   #'  posteriors for GMU 1 in the corresponding time_tmin1 (yr1, 2020)
  #'   time_t <- rbind(yr2_noGMU1, yr3, yr4)
  #'   rownames(time_t) <- NULL
  #'   
  #'   #'  Posteriors from time t-1 (yr1, yr2, yr3)
  #'   #'  NOTE: yr2 for time_tmin1 includes GMU 1 posteriors because the corresponding 
  #'   #'  yr3 posteriors in time_t also include GMU 1
  #'   time_tmin1 <- rbind(yr1, yr2, yr3)
  #'   rownames(time_tmin1) <- NULL
  #'   
  #'   #'  List posteriors from time t and t-1
  #'   timelag_list <- list(time_t, time_tmin1)
  #'   names(timelag_list) <- c("time_t", "time_tmin1") 
  #'   return(timelag_list)
  #' }
  #' #'  Run function for each species
  #' #'  NOTE: indexing is non-intuitive here (does not follow chronological time)
  #' #'  List order: time_t [[1]], time_tmin1 [[2]] even thought time_tmin1 --> time_t
  #' wolf_timelag <- stacked_posteriors(posteriors_20s[[1]], posteriors_21s[[1]], posteriors_22s[[1]], posteriors_23s[[1]], posteriors_21s_noGMU1[[1]])
  #' lion_timelag <- stacked_posteriors(posteriors_20s[[2]], posteriors_21s[[2]], posteriors_22s[[2]], posteriors_23s[[2]], posteriors_21s_noGMU1[[2]])
  #' bear_timelag <- stacked_posteriors(posteriors_20s[[3]], posteriors_21s[[3]], posteriors_22s[[3]], posteriors_23s[[3]], posteriors_21s_noGMU1[[3]])
  #' coy_timelag <- stacked_posteriors(posteriors_20s[[4]], posteriors_21s[[4]], posteriors_22s[[4]], posteriors_23s[[4]], posteriors_21s_noGMU1[[4]])
  #' elk_timelag <- stacked_posteriors(posteriors_20s[[5]], posteriors_21s[[5]], posteriors_22s[[5]], posteriors_23s[[5]], posteriors_21s_noGMU1[[5]])
  #' moose_timelag <- stacked_posteriors(posteriors_20s[[6]], posteriors_21s[[6]], posteriors_22s[[6]], posteriors_23s[[6]], posteriors_21s_noGMU1[[6]])
  #' wtd_timelag <- stacked_posteriors(posteriors_20s[[7]], posteriors_21s[[7]], posteriors_22s[[7]], posteriors_23s[[7]], posteriors_21s_noGMU1[[7]])
  #' 
  #' #'  Double check everything looks right (both time steps should have 37 rows)
  #' #'  [[1]] = time t (yr2, yr3, yr4); [[2]] = time t-1 (yr1, yr2, yr3)
  #' print(wolf_timelag[[1]]) # observations [c(1:8, 9:13),] in [[1]] should be observations [c(14:21, 32:36),] in [[2]]
  #' print(wolf_timelag[[2]]) 
  #' print(elk_timelag[[1]])
  #' print(elk_timelag[[2]])
  #' 
  #' #'  -------------------------
  #' ####  Format covariate data  ####
  #' #'  -------------------------
  #' #'  All covariates generated in Format_spatial_covariates_for_SMEs.R script
  #' 
  #' #'  Rename CLuster_unique to cluster
  #' cluster_poly_covs_df <- rename(cluster_poly_covs_df, cluster = Cluster_unique)
  #' 
  #' #'  Split out covariate data that do not need to be time-lagged (e.g., wolf_2020 is applied to RDI 2020)
  #' wolf_harvest <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, wolfharvest_2020_per100km, wolfharvest_2021_per100km, 
  #'                 wolfharvest_2022_per100km, wolfharvest_2023_per100km)
  #' wsi <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, WSI_20_gmu, WSI_21_gmu, WSI_22_gmu, WSI_23_gmu)
  #' forest <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, distForest20, distForest21, distForest22, distForest23)
  #' roads <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, roaddens)
  #' public <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, proppub)
  #' 
  #' #'  Split out covariate data that do need to be time lagged (e.g., bear_2019 is applied to RDI 2020)
  #' bear <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, Bear_per100km_2019, Bear_per100km_2020, Bear_per100km_2021, Bear_per100km_2022)
  #' lion <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, Lion_per100km_2019, Lion_per100km_2020, Lion_per100km_2021, Lion_per100km_2022)
  #' elk <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, Elk_per100km_2019, Elk_per100km_2020, Elk_per100km_2021, Elk_per100km_2022)
  #' moose <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, Moose_per100km_2019, Moose_per100km_2020, Moose_per100km_2021, Moose_per100km_2022)
  #' deer <- cluster_poly_covs_df %>%
  #'   dplyr::select(cluster, GMU, Deer_per100km_2019, Deer_per100km_2020, Deer_per100km_2021, Deer_per100km_2022)
  #' 
  #' cov_list <- list(wolf_harvest, wsi, forest, roads, public, bear, lion, elk, moose, deer)
  #' 
  #' #'  Call drop_gmu1() function to remove GMU1 observations 
  #' gmu1_numeric <- c(9, 10, 11, 12, 13, 14, 15, 16, 17, 18)
  #' names(cov_list) <- c("wolf_harvest", "wsi", "forest", "roads", "public", "bear", "lion", "elk", "moose", "deer")
  #' cov_list_noGMU1 <- lapply(cov_list, drop_gmu1, gmu = gmu1_numeric)
  #' #'  NOTE: this removed GMU1 observations from ALL years. OK to have them removed
  #' #'  from 2020 and 2021 data, but will need them for 2022 and 2023 data (and 2021
  #' #'  data depending on which annual time lag is beeing stacked)
  #' 
  #' #'  View an example of covaraite data w/o GMU1
  #' cov_list_noGMU1[[1]]
  #' 
  #' #'  Call function to format covariate data for stacked_posteriors() function
  #' #'  Function arguments include: yr1, yr2, yr3, yr4, yr2_noGMU1, cov_name
  #' #'    - yr1, yr2_noGMU1 need to come from cov_list_noGMU1 list
  #' #'    - yr2, yr3, and yr4 need to come from cov_list
  #' reformat_covs <- function(yr1, yr2, yr3, yr4, yr2_noGMU1, cov_name) {
  #'   #'  Convert vector to data frame
  #'   yr1 <- as.data.frame(yr1)
  #'   yr2 <- as.data.frame(yr2)
  #'   yr3 <- as.data.frame(yr3)
  #'   yr4 <- as.data.frame(yr4)
  #'   yr2_noGMU1 <- as.data.frame(yr2_noGMU1)
  #'   #'  Rename single column in each data frame to match covariate name
  #'   names(yr1) <- names(yr2) <- names(yr3) <- names(yr4) <- names(yr2_noGMU1) <- cov_name
  #'   #'  Relist annual covariate data
  #'   covs <- list(yr1, yr2, yr3, yr4, yr2_noGMU1)
  #'   #'  Rename list based on data set time period to make life easier below
  #'   names(covs) <- c("yr1", "yr2", "yr3", "yr4", "yr2_noGMU1")
  #'   return(covs)
  #' }
  #' #'  For reference: yr1, yr2, yr3, yr4, yr2_noGMU1, cov_name
  #' #'  cov_list indexing order: [[1]] wolf_harvest, [[2]] wsi, [[3]] forest, [[4]] roads, [[5]]public, 
  #' #'                           [[6]] bear, [[7]] lion, [[8]] elk, [[9]] moose, [[10]] deer
  #' wolf_harvest_cov <- reformat_covs(cov_list_noGMU1[[1]]$wolfharvest_2020_per100km, 
  #'                                   cov_list[[1]]$wolfharvest_2021_per100km, 
  #'                                   cov_list[[1]]$wolfharvest_2022_per100km, 
  #'                                   cov_list[[1]]$wolfharvest_2023_per100km, 
  #'                                   cov_list_noGMU1[[1]]$wolfharvest_2021_per100km, cov_name = "wolf_harvest")
  #' wsi_cov <- reformat_covs(cov_list_noGMU1[[2]]$WSI_20_gmu, cov_list[[2]]$WSI_21_gmu, cov_list[[2]]$WSI_22_gmu, cov_list[[2]]$WSI_23_gmu, 
  #'                          cov_list_noGMU1[[2]]$WSI_21_gmu, cov_name = "wsi")
  #' forest_cov <- reformat_covs(cov_list_noGMU1[[3]]$distForest20, cov_list[[3]]$distForest21, cov_list[[3]]$distForest22, cov_list[[3]]$distForest23, 
  #'                             cov_list_noGMU1[[3]]$distForest21, cov_name = "forest")
  #' roads_cov <- reformat_covs(cov_list_noGMU1[[4]]$roaddens, cov_list[[4]]$roaddens, cov_list[[4]]$roaddens, cov_list[[4]]$roaddens, 
  #'                            cov_list_noGMU1[[4]]$roaddens, cov_name = "roaddens")
  #' public_cov <- reformat_covs(cov_list_noGMU1[[5]]$proppub, cov_list[[5]]$proppub, cov_list[[5]]$proppub, cov_list[[5]]$proppub, 
  #'                             cov_list_noGMU1[[5]]$proppub, cov_name = "roaddens")
  #' bear_cov <- reformat_covs(cov_list_noGMU1[[6]]$Bear_per100km_2019, cov_list[[6]]$Bear_per100km_2020, cov_list[[6]]$Bear_per100km_2021, cov_list[[6]]$Bear_per100km_2022, 
  #'                           cov_list_noGMU1[[6]]$Bear_per100km_2020, cov_name = "bear_harvest")
  #' lion_cov <- reformat_covs(cov_list_noGMU1[[7]]$Lion_per100km_2019, cov_list[[7]]$Lion_per100km_2020, cov_list[[7]]$Lion_per100km_2021, cov_list[[7]]$Lion_per100km_2022,
  #'                           cov_list_noGMU1[[7]]$Lion_per100km_2020, cov_name = "lion_harvest")
  #' elk_cov <- reformat_covs(cov_list_noGMU1[[8]]$Elk_per100km_2019, cov_list[[8]]$Elk_per100km_2020, cov_list[[8]]$Elk_per100km_2021, cov_list[[8]]$Elk_per100km_2022,
  #'                          cov_list_noGMU1[[8]]$Elk_per100km_2020, cov_name = "elk_harvest")
  #' moose_cov <- reformat_covs(cov_list_noGMU1[[9]]$Moose_per100km_2019, cov_list[[9]]$Moose_per100km_2020, cov_list[[9]]$Moose_per100km_2021, cov_list[[9]]$Moose_per100km_2022,
  #'                           cov_list_noGMU1[[9]]$Moose_per100km_2020, cov_name = "moose_harvest")
  #' deer_cov <- reformat_covs(cov_list_noGMU1[[10]]$Deer_per100km_2019, cov_list[[10]]$Deer_per100km_2020, cov_list[[10]]$Deer_per100km_2021, cov_list[[10]]$Deer_per100km_2022,
  #'                           cov_list_noGMU1[[10]]$Deer_per100km_2020, cov_name = "deer_harvest")
  #' 
  #' #'  Call stacked_posteriors() function to stack covaraite data based on year and 
  #' #'  time lag for each covariate 
  #' #'  Function arguments include: yr1, yr2, yr3, yr4, yr2_noGMU1
  #' wolf_harvest_timelag <- stacked_posteriors(wolf_harvest_cov$yr1, wolf_harvest_cov$yr2, wolf_harvest_cov$yr3,
  #'                                            wolf_harvest_cov$yr4, wolf_harvest_cov$yr2_noGMU1)
  #' wsi_timelag <- stacked_posteriors(wsi_cov$yr1, wsi_cov$yr2, wsi_cov$yr3, wsi_cov$yr4, wsi_cov$yr2_noGMU1)
  #' forest_timelag <- stacked_posteriors(forest_cov$yr1, forest_cov$yr2, forest_cov$yr3, forest_cov$yr4, forest_cov$yr2_noGMU1)
  #' roads_timelag <- stacked_posteriors(roads_cov$yr1, roads_cov$yr2, roads_cov$yr3, roads_cov$yr4, roads_cov$yr2_noGMU1)
  #' public_timelag <- stacked_posteriors(public_cov$yr1, public_cov$yr2, public_cov$yr3, public_cov$yr4, public_cov$yr2_noGMU1)
  #' bear_harvest_timelag <- stacked_posteriors(bear_cov$yr1, bear_cov$yr2, bear_cov$yr3, bear_cov$yr4, bear_cov$yr2_noGMU1)
  #' lion_harvest_timelag <- stacked_posteriors(lion_cov$yr1, lion_cov$yr2, lion_cov$yr3, lion_cov$yr4, lion_cov$yr2_noGMU1)
  #' elk_harvest_timelag <- stacked_posteriors(elk_cov$yr1, elk_cov$yr2, elk_cov$yr3, elk_cov$yr4, elk_cov$yr2_noGMU1)
  #' moose_harvest_timelag <- stacked_posteriors(moose_cov$yr1, moose_cov$yr2, moose_cov$yr3, moose_cov$yr4, moose_cov$yr2_noGMU1)
  #' deer_harvest_timelag <- stacked_posteriors(deer_cov$yr1, deer_cov$yr2, deer_cov$yr3, deer_cov$yr4, deer_cov$yr2_noGMU1)
  #' 
  #' #'  Double check everything looks right (both time steps should have 37 rows)
  #' #'  [[1]] = time t (yr2, yr3, yr4); [[2]] = time t-1 (yr1, yr2, yr3)
  #' print(wolf_harvest_timelag[[1]]) # observations [c(1:8, 9:13),] in [[1]] should be observations [c(14:21, 32:36),] in [[2]]
  #' print(wolf_harvest_timelag[[2]]) 
  #' print(wsi_timelag[[1]][c(1:8, 9:13),]); print(wsi_timelag[[2]][c(14:21, 32:36),])                   # same value for entire GMU
  #' print(forest_timelag[[1]][c(1:8, 9:13),]); print(forest_timelag[[2]][c(14:21, 32:36),])
  #' print(roads_timelag[[1]][c(1:8, 9:13),]); print(roads_timelag[[2]][c(14:21, 32:36),])
  #' print(public_timelag[[1]][c(1:8, 9:13),]); print(public_timelag[[2]][c(14:21, 32:36),])
  #' print(bear_harvest_timelag[[1]][c(1:8, 9:13),]); print(bear_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  #' print(lion_harvest_timelag[[1]][c(1:8, 9:13),]); print(lion_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  #' print(elk_harvest_timelag[[1]][c(1:8, 9:13),]); print(elk_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  #' print(moose_harvest_timelag[[1]][c(1:8, 9:13),]); print(moose_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  #' print(deer_harvest_timelag[[1]][c(1:8, 9:13),]); print(deer_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  #' 
  #' #'  ------------------------
  #' ####  Standardize all data  ####
  #' #'  ------------------------
  #' #'  Standardize posterior means and SD
  #' z_transform_posteriors <- function(post_summary) {
  #'   #'  Find the empirical mean and SD of the posterior means
  #'   post_mean_mu <- mean(post_summary$posterior_mu)
  #'   post_mean_sd <- sd(post_summary$posterior_mu)
  #'   
  #'   #'  Standardize the posterior means AND standard deviations
  #'   #'  Standardizing posterior SD ensures variance scales with the variable
  #'   post_summary <- post_summary %>%
  #'     mutate(posterior_mu_z = (posterior_mu - post_mean_mu) / post_mean_sd,
  #'            posterior_sd_z = posterior_sd / post_mean_sd)
  #'   
  #'   return(post_summary)
  #' }
  #' wolf_timelag_z <- lapply(wolf_timelag, z_transform_posteriors)
  #' lion_timelag_z <- lapply(lion_timelag, z_transform_posteriors)
  #' bear_timelag_z <- lapply(bear_timelag, z_transform_posteriors)
  #' coy_timelag_z <- lapply(coy_timelag, z_transform_posteriors)
  #' elk_timelag_z <- lapply(elk_timelag, z_transform_posteriors)
  #' moose_timelag_z <- lapply(moose_timelag, z_transform_posteriors)
  #' wtd_timelag_z <- lapply(wtd_timelag, z_transform_posteriors)
  #' 
  #' #'  Take a look
  #' #'  Note: standardized values will be different between t and t-1 for same raw
  #' #'  posterior means because different years are stacked together, leading to 
  #' #'  different empirical means and SD for t vs t-1
  #' print(wolf_timelag_z$time_t)
  #' print(wolf_timelag_z$time_tmin1)
  #' print(wtd_timelag_z$time_t)
  #' print(wtd_timelag_z$time_tmin1)
  #' 
  #' #'  Standardize covariate data
  #' z_transform_covs <- function(covs, cov_name) {
  #'   covs_z <- covs %>%
  #'     mutate(cov_z = scale(.),
  #'            cov_z = as.numeric(cov_z))
  #'   names(covs_z) <- c(cov_name, paste0(cov_name, "_z"))
  #'   return(covs_z)
  #' }
  #' wolf_harv_timelag_z <- lapply(wolf_harvest_timelag, z_transform_covs, cov_name = "wolf_harv")
  #' wsi_timelag_z <- lapply(wsi_timelag, z_transform_covs, cov_name = "wsi")
  #' forest_timelag_z <- lapply(forest_timelag, z_transform_covs, cov_name = "prop_disturb")
  #' roads_timelag_z <- lapply(roads_timelag, z_transform_covs, cov_name = "road_density")
  #' public_timelag_z <- lapply(public_timelag, z_transform_covs, cov_name = "public_land")
  #' bear_harv_timelag_z <- lapply(bear_harvest_timelag, z_transform_covs, cov_name = "bear_harv")
  #' lion_harv_timelag_z <- lapply(lion_harvest_timelag, z_transform_covs, cov_name = "lion_harv")
  #' elk_harv_timelag_z <- lapply(elk_harvest_timelag, z_transform_covs, cov_name = "elk_harv")
  #' moose_harv_timelag_z <- lapply(moose_harvest_timelag, z_transform_covs, cov_name = "moose_harv")
  #' deer_harv_timelag_z <- lapply(deer_harvest_timelag, z_transform_covs, cov_name = "deer_harv")
  #' 
  #' #'  List species-specific lists 
  #' #'  First list indexes by species (wolf [[1]] .... wtd [[7]])
  #' #'  Second list (per species) indexes by time (time_t [[1]] & time_tmin1 [[2]])
  #' post_summaries <- list(wolf_timelag_z, lion_timelag_z, bear_timelag_z, coy_timelag_z, 
  #'                        elk_timelag_z, moose_timelag_z, wtd_timelag_z)
  #' names(post_summaries) <- c("wolf", "lion", "bear", "coy", "elk", "moose", "wtd")
  #' #'  List z-transformed covariates following similar indexing
  #' #'  Only really need time_tmin1 [[2]] for time lagged effect of covs in SEM
  #' covs_ztransformed <- list(wolf_harv_timelag_z, wsi_timelag_z, forest_timelag_z,
  #'                           roads_timelag_z, public_timelag_z, bear_harv_timelag_z,
  #'                           lion_harv_timelag_z, elk_harv_timelag_z, moose_harv_timelag_z, 
  #'                           deer_harv_timelag_z)
  #' names(covs_ztransformed) <- c("wolf_harv", "wsi", "prop_disturbed", "road_density", "public_land",
  #'                               "bear_harv", "lion_harv", "elk_harv", "moose_harv", "deer_harv")
  
  
  
  
  

  