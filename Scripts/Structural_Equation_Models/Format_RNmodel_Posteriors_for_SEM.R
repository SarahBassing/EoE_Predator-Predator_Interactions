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
  cluster_posterior_summary <- function(mod_post, cluster_out) {
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
      cluster = cluster_out,
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
  
  #'  Run function for each year and species
  posteriors_20s <- lapply(mods_yr1, cluster_posterior_summary, cluster_out = rdi.cl_2020)
  posteriors_21s <- lapply(mods_yr2, cluster_posterior_summary, cluster_out = rdi.cl_2021.2023)
  posteriors_22s <- lapply(mods_yr3, cluster_posterior_summary, cluster_out = rdi.cl_2021.2023)
  posteriors_23s <- lapply(mods_yr4, cluster_posterior_summary, cluster_out = rdi.cl_2021.2023)
  
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
  
  #'  Create new posteriors for 2021 that exclude GMU 1 RDI posteriors 
  #'  Needed below when stacking data for time t vs t-1
  drop_gmu1 <- function(post_rdi, gmu) {  
    #'  Filter out GMU1 clusters
    noGMU1 <- post_rdi %>%
      filter(!cluster %in% gmu)

    return(noGMU1)
  }
  gmu1 <- c("rdi.cl9", "rdi.cl10", "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", 
            "rdi.cl15", "rdi.cl16", "rdi.cl17", "rdi.cl18")
  posteriors_21s_noGMU1 <- lapply(posteriors_21s, drop_gmu1, gmu = gmu1) 
  
  #'  Stack species-specific estimated posterior summaries across iterations 
  #'  (e.g., time t-1 = wolf 2020, 2021, 2022 affects time t = wolf 2021, 2022, 2023)
  #'  yr1 2020 --> yr2 2021
  #'  yr2 2021 --> yr3 2022
  #'  yr3 2022 --> yr4 2023
  stacked_posteriors <- function(yr1, yr2, yr3, yr4, yr2_noGMU1) {
    #'  Posteriors from time t (yr2, yr3, yr4)
    #'  NOTE: yr2 for time_t excludes GMU 1 posteriors because there were no
    #'  posteriors for GMU 1 in the corresponding time_tmin1 (yr1, 2020)
    time_t <- rbind(yr2_noGMU1, yr3, yr4)
    rownames(time_t) <- NULL
    
    #'  Posteriors from time t-1 (yr1, yr2, yr3)
    #'  NOTE: yr2 for time_tmin1 includes GMU 1 posteriors because the corresponding 
    #'  yr3 posteriors in time_t also include GMU 1
    time_tmin1 <- rbind(yr1, yr2, yr3)
    rownames(time_tmin1) <- NULL
    
    #'  List posteriors from time t and t-1
    timelag_list <- list(time_t, time_tmin1)
    names(timelag_list) <- c("time_t", "time_tmin1") 
    return(timelag_list)
  }
  #'  Run function for each species
  #'  NOTE: indexing is non-intuitive here (does not follow chronological time)
  #'  List order: time_t [[1]], time_tmin1 [[2]] even thought time_tmin1 --> time_t
  wolf_timelag <- stacked_posteriors(posteriors_20s[[1]], posteriors_21s[[1]], posteriors_22s[[1]], posteriors_23s[[1]], posteriors_21s_noGMU1[[1]])
  lion_timelag <- stacked_posteriors(posteriors_20s[[2]], posteriors_21s[[2]], posteriors_22s[[2]], posteriors_23s[[2]], posteriors_21s_noGMU1[[2]])
  bear_timelag <- stacked_posteriors(posteriors_20s[[3]], posteriors_21s[[3]], posteriors_22s[[3]], posteriors_23s[[3]], posteriors_21s_noGMU1[[3]])
  coy_timelag <- stacked_posteriors(posteriors_20s[[4]], posteriors_21s[[4]], posteriors_22s[[4]], posteriors_23s[[4]], posteriors_21s_noGMU1[[4]])
  elk_timelag <- stacked_posteriors(posteriors_20s[[5]], posteriors_21s[[5]], posteriors_22s[[5]], posteriors_23s[[5]], posteriors_21s_noGMU1[[5]])
  moose_timelag <- stacked_posteriors(posteriors_20s[[6]], posteriors_21s[[6]], posteriors_22s[[6]], posteriors_23s[[6]], posteriors_21s_noGMU1[[6]])
  wtd_timelag <- stacked_posteriors(posteriors_20s[[7]], posteriors_21s[[7]], posteriors_22s[[7]], posteriors_23s[[7]], posteriors_21s_noGMU1[[7]])
  
  #'  Double check everything looks right (both time steps should have 37 rows)
  #'  [[1]] = time t (yr2, yr3, yr4); [[2]] = time t-1 (yr1, yr2, yr3)
  print(wolf_timelag[[1]]) # observations [c(1:8, 9:13),] in [[1]] should be observations [c(14:21, 32:36),] in [[2]]
  print(wolf_timelag[[2]]) 
  print(elk_timelag[[1]])
  print(elk_timelag[[2]])
  
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
  
  cov_list <- list(wolf_harvest, wsi, forest, roads, public, bear, lion, elk, moose, deer)
  
  #'  Call drop_gmu1() function to remove GMU1 observations 
  gmu1_numeric <- c(9, 10, 11, 12, 13, 14, 15, 16, 17, 18)
  names(cov_list) <- c("wolf_harvest", "wsi", "forest", "roads", "public", "bear", "lion", "elk", "moose", "deer")
  cov_list_noGMU1 <- lapply(cov_list, drop_gmu1, gmu = gmu1_numeric)
  #'  NOTE: this removed GMU1 observations from ALL years. OK to have them removed
  #'  from 2020 and 2021 data, but will need them for 2022 and 2023 data (and 2021
  #'  data depending on which annual time lag is beeing stacked)
  
  #'  View an example of covaraite data w/o GMU1
  cov_list_noGMU1[[1]]
  
  #'  Call function to format covariate data for stacked_posteriors() function
  #'  Function arguments include: yr1, yr2, yr3, yr4, yr2_noGMU1, cov_name
  #'    - yr1, yr2_noGMU1 need to come from cov_list_noGMU1 list
  #'    - yr2, yr3, and yr4 need to come from cov_list
  reformat_covs <- function(yr1, yr2, yr3, yr4, yr2_noGMU1, cov_name) {
    #'  Convert vector to data frame
    yr1 <- as.data.frame(yr1)
    yr2 <- as.data.frame(yr2)
    yr3 <- as.data.frame(yr3)
    yr4 <- as.data.frame(yr4)
    yr2_noGMU1 <- as.data.frame(yr2_noGMU1)
    #'  Rename single column in each data frame to match covariate name
    names(yr1) <- names(yr2) <- names(yr3) <- names(yr4) <- names(yr2_noGMU1) <- cov_name
    #'  Relist annual covariate data
    covs <- list(yr1, yr2, yr3, yr4, yr2_noGMU1)
    #'  Rename list based on data set time period to make life easier below
    names(covs) <- c("yr1", "yr2", "yr3", "yr4", "yr2_noGMU1")
    return(covs)
  }
  #'  For reference: yr1, yr2, yr3, yr4, yr2_noGMU1, cov_name
  #'  cov_list indexing order: [[1]] wolf_harvest, [[2]] wsi, [[3]] forest, [[4]] roads, [[5]]public, 
  #'                           [[6]] bear, [[7]] lion, [[8]] elk, [[9]] moose, [[10]] deer
  wolf_harvest_cov <- reformat_covs(cov_list_noGMU1[[1]]$wolfharvest_2020_per100km, 
                                    cov_list[[1]]$wolfharvest_2021_per100km, 
                                    cov_list[[1]]$wolfharvest_2022_per100km, 
                                    cov_list[[1]]$wolfharvest_2023_per100km, 
                                    cov_list_noGMU1[[1]]$wolfharvest_2021_per100km, cov_name = "wolf_harvest")
  wsi_cov <- reformat_covs(cov_list_noGMU1[[2]]$WSI_20_gmu, cov_list[[2]]$WSI_21_gmu, cov_list[[2]]$WSI_22_gmu, cov_list[[2]]$WSI_23_gmu, 
                           cov_list_noGMU1[[2]]$WSI_21_gmu, cov_name = "wsi")
  forest_cov <- reformat_covs(cov_list_noGMU1[[3]]$distForest20, cov_list[[3]]$distForest21, cov_list[[3]]$distForest22, cov_list[[3]]$distForest23, 
                              cov_list_noGMU1[[3]]$distForest21, cov_name = "forest")
  roads_cov <- reformat_covs(cov_list_noGMU1[[4]]$roaddens, cov_list[[4]]$roaddens, cov_list[[4]]$roaddens, cov_list[[4]]$roaddens, 
                             cov_list_noGMU1[[4]]$roaddens, cov_name = "roaddens")
  public_cov <- reformat_covs(cov_list_noGMU1[[5]]$proppub, cov_list[[5]]$proppub, cov_list[[5]]$proppub, cov_list[[5]]$proppub, 
                              cov_list_noGMU1[[5]]$proppub, cov_name = "roaddens")
  bear_cov <- reformat_covs(cov_list_noGMU1[[6]]$Bear_per100km_2019, cov_list[[6]]$Bear_per100km_2020, cov_list[[6]]$Bear_per100km_2021, cov_list[[6]]$Bear_per100km_2022, 
                            cov_list_noGMU1[[6]]$Bear_per100km_2020, cov_name = "bear_harvest")
  lion_cov <- reformat_covs(cov_list_noGMU1[[7]]$Lion_per100km_2019, cov_list[[7]]$Lion_per100km_2020, cov_list[[7]]$Lion_per100km_2021, cov_list[[7]]$Lion_per100km_2022,
                            cov_list_noGMU1[[7]]$Lion_per100km_2020, cov_name = "lion_harvest")
  elk_cov <- reformat_covs(cov_list_noGMU1[[8]]$Elk_per100km_2019, cov_list[[8]]$Elk_per100km_2020, cov_list[[8]]$Elk_per100km_2021, cov_list[[8]]$Elk_per100km_2022,
                           cov_list_noGMU1[[8]]$Elk_per100km_2020, cov_name = "elk_harvest")
  moose_cov <- reformat_covs(cov_list_noGMU1[[9]]$Moose_per100km_2019, cov_list[[9]]$Moose_per100km_2020, cov_list[[9]]$Moose_per100km_2021, cov_list[[9]]$Moose_per100km_2022,
                            cov_list_noGMU1[[9]]$Moose_per100km_2020, cov_name = "moose_harvest")
  deer_cov <- reformat_covs(cov_list_noGMU1[[10]]$Deer_per100km_2019, cov_list[[10]]$Deer_per100km_2020, cov_list[[10]]$Deer_per100km_2021, cov_list[[10]]$Deer_per100km_2022,
                            cov_list_noGMU1[[10]]$Deer_per100km_2020, cov_name = "deer_harvest")
  
  #'  Call stacked_posteriors() function to stack covaraite data based on year and 
  #'  time lag for each covariate 
  #'  Function arguments include: yr1, yr2, yr3, yr4, yr2_noGMU1
  wolf_harvest_timelag <- stacked_posteriors(wolf_harvest_cov$yr1, wolf_harvest_cov$yr2, wolf_harvest_cov$yr3,
                                             wolf_harvest_cov$yr4, wolf_harvest_cov$yr2_noGMU1)
  wsi_timelag <- stacked_posteriors(wsi_cov$yr1, wsi_cov$yr2, wsi_cov$yr3, wsi_cov$yr4, wsi_cov$yr2_noGMU1)
  forest_timelag <- stacked_posteriors(forest_cov$yr1, forest_cov$yr2, forest_cov$yr3, forest_cov$yr4, forest_cov$yr2_noGMU1)
  roads_timelag <- stacked_posteriors(roads_cov$yr1, roads_cov$yr2, roads_cov$yr3, roads_cov$yr4, roads_cov$yr2_noGMU1)
  public_timelag <- stacked_posteriors(public_cov$yr1, public_cov$yr2, public_cov$yr3, public_cov$yr4, public_cov$yr2_noGMU1)
  bear_harvest_timelag <- stacked_posteriors(bear_cov$yr1, bear_cov$yr2, bear_cov$yr3, bear_cov$yr4, bear_cov$yr2_noGMU1)
  lion_harvest_timelag <- stacked_posteriors(lion_cov$yr1, lion_cov$yr2, lion_cov$yr3, lion_cov$yr4, lion_cov$yr2_noGMU1)
  elk_harvest_timelag <- stacked_posteriors(elk_cov$yr1, elk_cov$yr2, elk_cov$yr3, elk_cov$yr4, elk_cov$yr2_noGMU1)
  moose_harvest_timelag <- stacked_posteriors(moose_cov$yr1, moose_cov$yr2, moose_cov$yr3, moose_cov$yr4, moose_cov$yr2_noGMU1)
  deer_harvest_timelag <- stacked_posteriors(deer_cov$yr1, deer_cov$yr2, deer_cov$yr3, deer_cov$yr4, deer_cov$yr2_noGMU1)
  
  #'  Double check everything looks right (both time steps should have 37 rows)
  #'  [[1]] = time t (yr2, yr3, yr4); [[2]] = time t-1 (yr1, yr2, yr3)
  print(wolf_harvest_timelag[[1]]) # observations [c(1:8, 9:13),] in [[1]] should be observations [c(14:21, 32:36),] in [[2]]
  print(wolf_harvest_timelag[[2]]) 
  print(wsi_timelag[[1]][c(1:8, 9:13),]); print(wsi_timelag[[2]][c(14:21, 32:36),])                   # same value for entire GMU
  print(forest_timelag[[1]][c(1:8, 9:13),]); print(forest_timelag[[2]][c(14:21, 32:36),])
  print(roads_timelag[[1]][c(1:8, 9:13),]); print(roads_timelag[[2]][c(14:21, 32:36),])
  print(public_timelag[[1]][c(1:8, 9:13),]); print(public_timelag[[2]][c(14:21, 32:36),])
  print(bear_harvest_timelag[[1]][c(1:8, 9:13),]); print(bear_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  print(lion_harvest_timelag[[1]][c(1:8, 9:13),]); print(lion_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  print(elk_harvest_timelag[[1]][c(1:8, 9:13),]); print(elk_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  print(moose_harvest_timelag[[1]][c(1:8, 9:13),]); print(moose_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  print(deer_harvest_timelag[[1]][c(1:8, 9:13),]); print(deer_harvest_timelag[[2]][c(14:21, 32:36),]) # same value for entire GMU
  
  #'  ------------------------
  ####  Standardize all data  ####
  #'  ------------------------
  #'  Standardize posterior means and SD
  z_transform_posteriors <- function(post_summary) {
    #'  Find the empirical mean and SD of the posterior means
    post_mean_mu <- mean(post_summary$posterior_mu)
    post_mean_sd <- sd(post_summary$posterior_mu)
    
    #'  Standardize the posterior means AND standard deviations
    #'  Standardizing posterior SD ensures variance scales with the variable
    post_summary <- post_summary %>%
      mutate(posterior_mu_z = (posterior_mu - post_mean_mu) / post_mean_sd,
             posterior_sd_z = posterior_sd / post_mean_sd)
    
    return(post_summary)
  }
  wolf_timelag_z <- lapply(wolf_timelag, z_transform_posteriors)
  lion_timelag_z <- lapply(lion_timelag, z_transform_posteriors)
  bear_timelag_z <- lapply(bear_timelag, z_transform_posteriors)
  coy_timelag_z <- lapply(coy_timelag, z_transform_posteriors)
  elk_timelag_z <- lapply(elk_timelag, z_transform_posteriors)
  moose_timelag_z <- lapply(moose_timelag, z_transform_posteriors)
  wtd_timelag_z <- lapply(wtd_timelag, z_transform_posteriors)
 
  #'  Take a look
  #'  Note: standardized values will be different between t and t-1 for same raw
  #'  posterior means because different years are stacked together, leading to 
  #'  different empirical means and SD for t vs t-1
  print(wolf_timelag_z$time_t)
  print(wolf_timelag_z$time_tmin1)
  print(wtd_timelag_z$time_t)
  print(wtd_timelag_z$time_tmin1)
  
  #'  Standardize covariate data
  z_transform_covs <- function(covs, cov_name) {
    covs_z <- covs %>%
      mutate(cov_z = scale(.),
             cov_z = as.numeric(cov_z))
    names(covs_z) <- c(cov_name, paste0(cov_name, "_z"))
    return(covs_z)
  }
  wolf_harv_timelag_z <- lapply(wolf_harvest_timelag, z_transform_covs, cov_name = "wolf_harv")
  wsi_timelag_z <- lapply(wsi_timelag, z_transform_covs, cov_name = "wsi")
  forest_timelag_z <- lapply(forest_timelag, z_transform_covs, cov_name = "prop_disturb")
  roads_timelag_z <- lapply(roads_timelag, z_transform_covs, cov_name = "road_density")
  public_timelag_z <- lapply(public_timelag, z_transform_covs, cov_name = "public_land")
  bear_harv_timelag_z <- lapply(bear_harvest_timelag, z_transform_covs, cov_name = "bear_harv")
  lion_harv_timelag_z <- lapply(lion_harvest_timelag, z_transform_covs, cov_name = "lion_harv")
  elk_harv_timelag_z <- lapply(elk_harvest_timelag, z_transform_covs, cov_name = "elk_harv")
  moose_harv_timelag_z <- lapply(moose_harvest_timelag, z_transform_covs, cov_name = "moose_harv")
  deer_harv_timelag_z <- lapply(deer_harvest_timelag, z_transform_covs, cov_name = "deer_harv")
  
  #'  List species-specific lists 
  #'  First list indexes by species (wolf [[1]] .... wtd [[7]])
  #'  Second list (per species) indexes by time (time_t [[1]] & time_tmin1 [[2]])
  post_summaries <- list(wolf_timelag_z, lion_timelag_z, bear_timelag_z, coy_timelag_z, 
                         elk_timelag_z, moose_timelag_z, wtd_timelag_z)
  names(post_summaries) <- c("wolf", "lion", "bear", "coy", "elk", "moose", "wtd")
  #'  List z-transformed covariates following similar indexing
  #'  Only really need time_tmin1 [[2]] for time lagged effect of covs in SEM
  covs_ztransformed <- list(wolf_harv_timelag_z, wsi_timelag_z, forest_timelag_z,
                            roads_timelag_z, public_timelag_z, bear_harv_timelag_z,
                            lion_harv_timelag_z, elk_harv_timelag_z, moose_harv_timelag_z, 
                            deer_harv_timelag_z)
  names(covs_ztransformed) <- c("wolf_harv", "wsi", "prop_disturbed", "road_density", "public_land",
                                "bear_harv", "lion_harv", "elk_harv", "moose_harv", "deer_harv")
  
  #' ####  Save full posteriors (all iterations)  ####
  #' ####  Not necessary for how uncertainty from RN models are being propagated in SEM  ####
  #' #'  Function to save estimated posterior distributions of interest for each species,
  #' #'  cluster, and year
  #' cluster_mcmc_matrix <- function(mod_post, cluster_out) {
  #'   #'  Empty list to hold model outputs of interest
  #'   mod_out <- NULL
  #'   #'  Loop through posteriors and save outputs of interest
  #'   for(i in 1:length(cluster_out)) {
  #'     mod_out[[i]] <- mod_post$sims.list[names(mod_post$sims.list) %in% cluster_out[i]]
  #'   }
  #'   
  #'   #'  Create an empty matrix with dimensions specific to clusters & iterations in model
  #'   cluster_posteriors <- matrix(data = NA, nrow = length(cluster_out), ncol = length(mod_post$sims.list[[1]]), byrow = TRUE)
  #'   #'  Loop through matrix and fill with estimated posterior distribution for each RDI cluster
  #'   for(i in 1:nrow(cluster_posteriors)) {
  #'     cluster_posteriors[i,] <- mod_out[[i]][[1]]
  #'   }
  #' 
  #'   #' #'  Turn matrix into a data frame with uniqueCluster as first column
  #'   #' uniqueCluster <- as.data.frame(seq(1:nrow(cluster_posteriors)))
  #'   #' names(uniqueCluster) <- "uniqueCluster"
  #'   #' cluster_posteriors <- cbind(uniqueCluster, cluster_posteriors) 
  #'   
  #'   ####  ADD IN COVARIATE DATA HERE!!!!
  #'   
  #'   return(cluster_posteriors)
  #' }
  #' #'  Names of estimated posteriors of interest per year
  #' rdi.cl_2020 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9",
  #'                  "rdi.cl21", "rdi.cl22", "rdi.cl23", "rdi.cl24")
  #' rdi.cl_2021.2022 <- c("rdi.cl1", "rdi.cl2", "rdi.cl3", "rdi.cl4", "rdi.cl5", "rdi.cl6", "rdi.cl7", "rdi.cl8", "rdi.cl9", "rdi.cl10", 
  #'                       "rdi.cl11", "rdi.cl12", "rdi.cl13", "rdi.cl14", "rdi.cl15", "rdi.cl16", "rdi.cl17", "rdi.cl18", "rdi.cl19", "rdi.cl20", 
  #'                       "rdi.cl21", "rdi.cl22", "rdi.cl23", "rdi.cl24")
  #' #'  Run function for each year and species
  #' posteriors_20s <- lapply(mods_yr1, cluster_mcmc_matrix, cluster_out = rdi.cl_2020)
  #' posteriors_21s <- lapply(mods_yr2, cluster_mcmc_matrix, cluster_out = rdi.cl_2021.2022)
  #' posteriors_22s <- lapply(mods_yr3, cluster_mcmc_matrix, cluster_out = rdi.cl_2021.2022)
  #' 
  #' #'  Name lists by species
  #' spp <- c("wolf", "mountain_lion", "bear_black", "coyote", "elk", "moose", "whitetailed_deer")
  #' names(posteriors_20s) <- spp
  #' names(posteriors_21s) <- spp
  #' names(posteriors_22s) <- spp
  #' 
  #' #'  Double check these look right
  #' #'  First 10 iterations for each cluster in matrix vs model output
  #' print(posteriors_20s[[1]][, 1:10]) # wolf (2020)
  #' RN_wolf_20s$sims.list$rdi.cl1[1:10]; RN_wolf_20s$sims.list$rdi.cl24[1:10] # note that rdi.cl10 - rdi.cl20 are not in 2020 posteriors
  #' print(posteriors_21s[[3]][, 1:10]) # bear (2021)
  #' RN_bear_21s$sims.list$rdi.cl1[1:10]; RN_bear_21s$sims.list$rdi.cl24[1:10]
  #' print(posteriors_22s[[7]][, 1:10]) # white-tailed deer (2022)
  #' RN_wtd_22s$sims.list$rdi.cl1[1:10]; RN_wtd_22s$sims.list$rdi.cl24[1:10]
  #' 
  #' #'  Create new posteriors for 2021 that exclude GMU 1 RDI posteriors 
  #' #'  Needed below when stacking data for time t vs t-1
  #' #'  posteriors_21s[[i]][10:20,] needs to go (corresponds to GMU 1 clusters)
  #' drop_gmu1 <- function(post_rdi) {
  #'   #'  Snag posteriors from GMU 10A and GMU 6
  #'   gmu10a_rdi <- post_rdi[1:9,]
  #'   gmu6_rdi <- post_rdi[21:24,]
  #'   
  #'   #'  Bind into a single matrix that no long includes posteriors from GMU 1
  #'   rdi_2021_no_gmu1 <- rbind(gmu10a_rdi, gmu6_rdi)
  #'   return(rdi_2021_no_gmu1)
  #' }
  #' posteriors_21s_noGMU1 <- lapply(posteriors_21s, drop_gmu1)
  #' 
  #' #'  Stack species-specific estimated posteriors across iterations 
  #' #'  (e.g., time t = wolf 2021 and 2022 vs. time t-1 = wolf 2020 and 2021)
  #' stacked_posteriors <- function(yr1, yr2, yr3, yr2_noGMU1) {
  #'   #'  Posteriors from time t 
  #'   #'  NOTE: yr2 for time_t excludes GMU 1 posteriors because there were no
  #'   #'  posteriors for GMU 1 in the corresponding time_tmin1 (yr1, 2020)
  #'   time_t <- rbind(yr2_noGMU1, yr3)
  #'   
  #'   #'  Posteriors from time t-1, given an annual lag from time 
  #'   #'  NOTE: yr2 for time_tmin1 includes GMU 1 posteriors because the corresponding 
  #'   #'  yr3 posteriors in time_t also include GMU 1
  #'   time_tmin1 <- rbind(yr1, yr2)
  #'   
  #'   #'  List posteriors from time t and t-1
  #'   timelag_list <- list(time_t, time_tmin1)
  #'   return(timelag_list)
  #' }
  #' #'  Run function for each species
  #' wolf_timelag <- stacked_posteriors(posteriors_20s[[1]], posteriors_21s[[1]], posteriors_22s[[1]], posteriors_21s_noGMU1[[1]])
  #' lion_timelag <- stacked_posteriors(posteriors_20s[[2]], posteriors_21s[[2]], posteriors_22s[[2]], posteriors_21s_noGMU1[[2]])
  #' bear_timelag <- stacked_posteriors(posteriors_20s[[3]], posteriors_21s[[3]], posteriors_22s[[3]], posteriors_21s_noGMU1[[3]])
  #' coy_timelag <- stacked_posteriors(posteriors_20s[[4]], posteriors_21s[[4]], posteriors_22s[[4]], posteriors_21s_noGMU1[[4]])
  #' elk_timelag <- stacked_posteriors(posteriors_20s[[5]], posteriors_21s[[5]], posteriors_22s[[5]], posteriors_21s_noGMU1[[5]])
  #' moose_timelag <- stacked_posteriors(posteriors_20s[[6]], posteriors_21s[[6]], posteriors_22s[[6]], posteriors_21s_noGMU1[[6]])
  #' wtd_timelag <- stacked_posteriors(posteriors_20s[[7]], posteriors_21s[[7]], posteriors_22s[[7]], posteriors_21s_noGMU1[[7]])
  #' 
  #' #'  Double check everything looks right (both time steps should have 37 rows)
  #' #'  [[1]] = time t (yr2 and yr3); [[2]] = time t-1 (yr1 and yr2)
  #' print(wolf_timelag[[1]][,1:10])
  #' print(wolf_timelag[[2]][,1:10])
  #' print(elk_timelag[[1]][,1:10])
  #' print(elk_timelag[[2]][,1:10])
  
  

  