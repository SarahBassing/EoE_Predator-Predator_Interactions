  #'  -----------------------------------
  #'  Co-Occurrence plots 
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  Februray 2023
  #'  -----------------------------------
  #'  Plot species interactions using multispecies occupancy model results.
  #'  ------------------------------------
  
  #'  Load libraries
  library(unmarked)
  library(ggplot2)
  library(tidyverse)
  library(khroma)
  library(patchwork)
  
  #'  Read in models
  # load(file = "./Outputs/MultiSpp_OccMod_Outputs/MultiSpp_CoOcc_Models_2spp.RData")
  load(file = "./Outputs/MultiSpp_OccMod_Outputs/MultiSpp_CoOcc_Models_2spp_PredatorCamsOnly.RData")
  
  #'  Read in covariate data
  load("./Data/Covariates_extracted/Covariate_skinny_EoE20s.RData")
  load("./Data/Covariates_extracted/Covariate_skinny_EoE20w.RData")
  load("./Data/Covariates_extracted/Covariate_skinny_EoE21s.RData")
  
  scale_cov <- function(cov) {
    #'  Identify range of covariate of interest
    r <- range(cov)
    #'  Create sequence of values starting and ending with range values
    x <- seq(r[1], r[2], length.out = 100)
    #'  Scale values based on covariate mean and standard deviation
    x_scaled <- (x - mean(cov)) / (sd(cov))
    #'  Make single data frame of scaled and unscaled covariate
    x_cov <- as.data.frame(cbind(x, x_scaled))
    
    return(x_cov)
  }
  scaled_elev <- lapply(list(stations_skinny_eoe20s$Elev, stations_skinny_eoe21s$Elev), scale_cov)
  scaled_forest <- lapply(list(stations_skinny_eoe20s$PercForest, stations_skinny_eoe21s$PercForest), scale_cov)
  scaled_dist2burbs <- lapply(list(stations_skinny_eoe20s$Dist2Burbs, stations_skinny_eoe21s$Dist2Burbs), scale_cov)
  scaled_lognearestrd <- lapply(list(stations_skinny_eoe20s$logNearestRd, stations_skinny_eoe21s$logNearestRd), scale_cov)
  
  scale_cov_round <- function(cov) {
    #'  Identify range of covariate of interest
    r <- range(cov)
    #'  Create sequence of values starting and ending with range values
    x <- seq(r[1], r[2], length.out = 100)
    #'  Make sure the values are integers (can't have partial individuals)
    x <- round(x, 0)
    #'  Remove duplicate integers if range is <100 individuals
    x <- unique(x)
    #'  Scale values based on covariate mean and standard deviation
    x_scaled <- (x - mean(cov)) / (sd(cov))
    #'  Make single data frame of scaled and unscaled covariate
    x_cov <- as.data.frame(cbind(x, x_scaled))
    
    return(x_cov)
  }
  scaled_elk <- lapply(list(stations_skinny_eoe20s$Nelk, stations_skinny_eoe21s$Nelk), scale_cov_round)
  scaled_human <- lapply(list(stations_skinny_eoe20s$Nhuman, stations_skinny_eoe21s$Nhuman), scale_cov_round)
  scaled_lagomorph <- lapply(list(stations_skinny_eoe20s$Nlagomorph, stations_skinny_eoe21s$Nlagomorph), scale_cov_round)
  scaled_livestock <- lapply(list(stations_skinny_eoe20s$Nlivestock, stations_skinny_eoe21s$Nlivestock), scale_cov_round)
  scaled_moose <- lapply(list(stations_skinny_eoe20s$Nmoose, stations_skinny_eoe21s$Nmoose), scale_cov_round)
  scaled_md <- lapply(list(stations_skinny_eoe20s$Nmd, stations_skinny_eoe21s$Nmd), scale_cov_round)
  scaled_wtd <- lapply(list(stations_skinny_eoe20s$Nwtd, stations_skinny_eoe21s$Nwtd), scale_cov_round)
  scaled_bigdeer <- lapply(list(stations_skinny_eoe20s$Nbig_deer, stations_skinny_eoe21s$Nbig_deer), scale_cov_round)
  scaled_smalldeer <- lapply(list(stations_skinny_eoe20s$Nsmall_deer, stations_skinny_eoe21s$Nsmall_deer), scale_cov_round)
  
  
  #'  ----------------------------------------------------------
  ####  Covariate effects on predator-predator interaction term  ####
  #'  ----------------------------------------------------------
  #'  Function to predict species interactions in response to covariate of interest
  spp_interactions_g <- function(mod, elev, forest, dist2sub, nearestrd, elk, human, 
                                 bunny, cow, moose, md, wtd, bigdeer, smdeer, spp1, spp2, cov) {
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, PercForest = forest, Dist2Burbs = dist2sub, 
                         logNearestRd = nearestrd, Nelk = elk, Nhuman = human,
                         Nlagomorph = bunny, Nlivestock = cow, Nmoose = moose,
                         Nmd = md, Nwtd = wtd, Nbig_deer = bigdeer, Nsmall_deer = smdeer)
    #'  Create characters for each species that include a "-", necessary for cond 
    #'  argument in predict when species is not present
    no_spp2 <- paste0("-",spp2); no_spp1 <- paste0("-",spp1)
    
    #'  Predict conditional occupancy when spp2 is absent
    spp2_absent <- predict(mod, type = "state", species = spp1, cond = no_spp2,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^4) %>%    #change to 10^5 when doing it for real
      mutate(Species = spp1,
             Interaction = paste0(spp2, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp2 is present
    spp2_present <- predict(mod, type = "state", species = spp1, cond = spp2,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^4) %>%   #change to 10^5 when doing it for real
      mutate(Species = spp1,
             Interaction = paste0(spp2, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Predict conditional occupancy when spp1 is absent
    spp1_absent <- predict(mod, type = "state", species = spp2, cond = no_spp1,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^4) %>%    #change to 10^5 when doing it for real
      mutate(Species = spp2,
             Interaction = paste0(spp1, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp1 is present
    spp1_present <- predict(mod, type = "state", species = spp2, cond = spp1,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^4) %>%   #change to 10^5 when doing it for real
      mutate(Species = spp2,
             Interaction = paste0(spp1, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Create one large data frame with all marginal probabilities
    sppX_df <- rbind(spp2_absent, spp2_present, spp1_absent, spp1_present)
    
    return(sppX_df)
  }
  #'  Predict cond. occupancy for each pairwise interaction over range of covariate values
  #'  NOTE: setting all other continuous variables to 0 (their mean value) and 
  #'  setting categorical variables (if any) as either 0 or 1
  #'  Indexing: first list [[1]] = 2020 data, [[2]] = 2021 data; 
  #'            second list [1] = raw covariate data, [2] = scaled covariate data
  
  ####  Predict co-occurrence across covariate values  ####
  #####  Wolf-Bobcat predictions  #####
  #' 2020 significant covariate predictions
  wb_nearestrd_20s <- spp_interactions_g(wb_20s_top, elev = 0, forest = 0, dist2sub = 0, 
                                       nearestrd = scaled_lognearestrd[[1]][,2], 
                                       elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                       md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                       spp1 = "wolf", spp2 = "bobcat", cov = scaled_lognearestrd[[1]][,1])
  wb_dist2burbs_20s <- spp_interactions_g(wb_20s_top, elev = 0, forest = 0,
                                          dist2sub = scaled_dist2burbs[[1]][,2], 
                                          nearestrd = 0, elk = 0, human = 0, bunny = 0, 
                                          cow = 0, moose = 0, md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                          spp1 = "wolf", spp2 = "bobcat", cov = scaled_dist2burbs[[1]][,1])
  wb_livestock_20s <- spp_interactions_g(wb_20s_top, elev = 0, forest = 0, dist2sub = 0, 
                                         nearestrd = 0, elk = 0, human = 0, bunny = 0, 
                                         cow = scaled_livestock[[1]][,2], moose = 0, 
                                         md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                         spp1 = "wolf", spp2 = "bobcat", cov = scaled_livestock[[1]][,1])
  #' 2021 significant covariate predictions
  wb_bunnies_21s <- spp_interactions_g(wb_21s_top, elev = 0, forest = 0, dist2sub = 0, 
                                       nearestrd = 0, elk = 0, human = 0, bunny = scaled_lagomorph[[2]][,2], 
                                       cow = 0, moose = 0, md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                       spp1 = "wolf", spp2 = "bobcat", cov = scaled_lagomorph[[2]][,1])
  
  #####  Wolf-Lion Predictions  #####
  #' 2020 significant covariate predictions
  wl_elev_20s <- spp_interactions_g(wl_20s_top, elev = scaled_elev[[1]][,2], forest = 0, dist2sub = 0, 
                                    nearestrd = 0, elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                    md = 0, wtd = 0, bigdeer = 0, smdeer = 0, spp1 = "wolf", 
                                    spp2 = "lion", cov = scaled_elev[[1]][,1])
  #' 2021 significant covariate predictions
  wl_nearestrd_21s <- spp_interactions_g(wl_21s_top, elev = 0, forest = 0, dist2sub = 0, 
                                       nearestrd = scaled_lognearestrd[[2]][,2], 
                                       elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                       md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                       spp1 = "wolf", spp2 = "lion", cov = scaled_lognearestrd[[2]][,1])
  wl_elev_21s <- spp_interactions_g(wl_21s_top, elev = scaled_elev[[2]][,2], forest = 0, dist2sub = 0, 
                                    nearestrd = 0, elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                    md = 0, wtd = 0, bigdeer = 0, smdeer = 0, spp1 = "wolf", 
                                    spp2 = "lion", cov = scaled_elev[[2]][,1])
  
  
  #####  Wolf-Coyote Predictions  #####
  #' 2020 significant covariate predictions
  wc_dist2burbs_20s <- spp_interactions_g(wc_20s_top, elev = 0, forest = 0,
                                          dist2sub = scaled_dist2burbs[[1]][,2], 
                                          nearestrd = 0, elk = 0, human = 0, bunny = 0, 
                                          cow = 0, moose = 0, md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                          spp1 = "wolf", spp2 = "coyote", cov = scaled_dist2burbs[[1]][,1])
  wc_elev_20s <- spp_interactions_g(wc_20s_top, elev = scaled_elev[[1]][,2], forest = 0, dist2sub = 0, 
                                    nearestrd = 0, elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                    md = 0, wtd = 0, bigdeer = 0, smdeer = 0, spp1 = "wolf", 
                                    spp2 = "coyote", cov = scaled_elev[[1]][,1])
  #' 2021 significant covariate predictions
  wc_elev_21s <- spp_interactions_g(wc_21s_top, elev = scaled_elev[[2]][,2], forest = 0, dist2sub = 0, 
                                    nearestrd = 0, elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                    md = 0, wtd = 0, bigdeer = 0, smdeer = 0, spp1 = "wolf", 
                                    spp2 = "coyote", cov = scaled_elev[[2]][,1])
  wc_bigdeer_21s <- spp_interactions_g(wc_21s_top, elev = 0, forest = 0, dist2sub = 0, 
                                       nearestrd = 0, elk = 0, human = 0, bunny = 0, cow = 0, 
                                       moose = 0, md = 0, wtd = 0, bigdeer = scaled_bigdeer[[2]][,2], smdeer = 0,
                                       spp1 = "wolf", spp2 = "coyote", cov = scaled_bigdeer[[2]][,1])
  wc_bunnies_21s <- spp_interactions_g(wc_21s_top, elev = 0, forest = 0, dist2sub = 0, 
                                       nearestrd = 0, elk = 0, human = 0, bunny = scaled_lagomorph[[2]][,2], 
                                       cow = 0, moose = 0, md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                       spp1 = "wolf", spp2 = "coyote", cov = scaled_lagomorph[[2]][,1])
  wc_nearestrd_21s <- spp_interactions_g(wc_21s_top, elev = 0, forest = 0, dist2sub = 0,
                                         nearestrd = scaled_lognearestrd[[2]][,2], 
                                         elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                         md = 0, wtd = 0, bigdeer = 0, smdeer = 0, spp1 = "wolf", 
                                         spp2 = "coyote", cov = scaled_lognearestrd[[2]][,1])
  
  #####  Lion-Coyote Predictions  #####
  #'  2020 significant covariate predictions
  lc_dist2burbs_20s <- spp_interactions_g(lc_20s_top, elev = 0, forest = 0,
                                          dist2sub = scaled_dist2burbs[[1]][,2], 
                                          nearestrd = 0, elk = 0, human = 0, bunny = 0, 
                                          cow = 0, moose = 0, md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                          spp1 = "lion", spp2 = "coyote", cov = scaled_dist2burbs[[1]][,1])
  lc_nearestrd_20s <- spp_interactions_g(lc_20s_top, elev = 0, forest = 0, dist2sub = 0,
                                         nearestrd = scaled_lognearestrd[[1]][,2], 
                                         elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                         md = 0, wtd = 0, bigdeer = 0, smdeer = 0, spp1 = "lion", 
                                         spp2 = "coyote", cov = scaled_lognearestrd[[1]][,1])
  lc_wtd_20s <- spp_interactions_g(lc_20s_top, elev = 0, forest = 0, dist2sub = 0, 
                                   nearestrd = 0, elk = 0, human = 0, bunny = 0, cow = 0, 
                                   moose = 0, md = 0, wtd = scaled_wtd[[1]][,2], bigdeer = 0, smdeer = 0,
                                   spp1 = "lion", spp2 = "coyote", cov = scaled_wtd[[1]][,1])
  lc_elev_20s <- spp_interactions_g(lc_20s_top, elev = scaled_elev[[1]][,2], forest = 0,
                                    dist2sub = 0, nearestrd = 0, elk = 0, human = 0, 
                                    bunny = 0, cow = 0, moose = 0, md = 0, wtd = 0, 
                                    bigdeer = 0, smdeer = 0, spp1 = "lion", spp2 = "coyote", 
                                    cov = scaled_elev[[1]][,1])
  
  #'  2021 significant covariate predictions
  lc_nearestrd_21s <- spp_interactions_g(lc_21s_top, elev = 0, forest = 0, dist2sub = 0,
                                         nearestrd = scaled_lognearestrd[[2]][,2], 
                                         elk = 0, human = 0, bunny = 0, cow = 0, moose = 0, 
                                         md = 0, wtd = 0, bigdeer = 0, smdeer = 0, spp1 = "lion", 
                                         spp2 = "coyote", cov = scaled_lognearestrd[[2]][,1])
  
  #####  Bobcat-Coyote Predictions  #####
  #'  2021 significant covariate predictions
  bc_human_21s <- spp_interactions_g(bc_21s_top, elev = 0, forest = 0, dist2sub = 0, 
                                   nearestrd = 0, elk = 0, human = scaled_human[[2]][,2], 
                                   bunny = 0, cow = 0, moose = 0, md = 0, wtd = 0, 
                                   bigdeer = 0, smdeer = 0, spp1 = "bobcat", spp2 = "coyote", 
                                   cov = scaled_human[[2]][,1])
  
  
  
  
  
  ####  Co-occurrence Figures  ####
  #####  Wolf-Bobcat Plots  #####
  #'  2020 co-occurrence 
  wb_dist2burbs_20s_predict <- wb_dist2burbs_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("bobcat absent", "bobcat present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "bobcat" = "Bobcat") 
  wb_dist2burbs_20s_facet <- ggplot(wb_dist2burbs_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Distance to urban/suburban area (m)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of proximity to human settelement on wolf - bobcat co-occurrence 2020")
  wb_dist2burbs_20s_facet
  
  wb_nearestrd_20s_predict <- wb_nearestrd_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("bobcat absent", "bobcat present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "bobcat" = "Bobcat") 
  wb_nearestrd_20s_facet <- ggplot(wb_nearestrd_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Log distance to nearest road") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of log distance to nearest road on wolf - bobcat co-occurrence 2020")
  wb_nearestrd_20s_facet
  
  wb_livestock_20s_predict <- wb_livestock_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("bobcat absent", "bobcat present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "bobcat" = "Bobcat") 
  wb_livestock_20s_facet <- ggplot(wb_livestock_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Livestock detections") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of livestock activity on wolf - bobcat co-occurrence 2020")
  wb_livestock_20s_facet
  
  #'  2021 co-occurrence
  wb_bunnies_21s_predict <- wb_bunnies_21s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("bobcat absent", "bobcat present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "bobcat" = "Bobcat") 
  wb_bunnies_21s_facet <- ggplot(wb_bunnies_21s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Lagomorph detections") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of lagomorph activity on wolf - bobcat co-occurrence 2021")
  wb_bunnies_21s_facet
  
  #####  Wolf-Lion Plots  #####
  #'  2020 co-occurrence
  wl_elev_20s_predict <- wl_elev_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("lion absent", "lion present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "lion" = "Lion") 
  wl_elev_20s_facet <- ggplot(wl_elev_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "lion absent" = "Lion absent", "lion present" = "Lion present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "lion absent" = "Lion absent", "lion present" = "Lion present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Elevation (m)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of elevation on wolf - lion co-occurrence 2020")
  wl_elev_20s_facet
  
  #'  2021 co-occurrence
  wl_nearestrd_21s_predict <- wl_nearestrd_21s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("lion absent", "lion present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "lion" = "Lion") 
  wl_nearestrd_21s_facet <- ggplot(wl_nearestrd_21s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "lion absent" = "Lion absent", "lion present" = "Lion present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "lion absent" = "Lion absent", "lion present" = "Lion present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Log distance to nearest road") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of log distance to nearest road on wolf - lion co-occurrence 2021")
  wl_nearestrd_21s_facet
  
  #####  Wolf-Coyote Plots  #####
  #'  2020 co-occurrence
  wc_dist2burbs_20s_predict <- wc_dist2burbs_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "coyote" = "Coyote") 
  wc_dist2burbs_20s_facet <- ggplot(wc_dist2burbs_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Distance to urban/suburban area (m)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of proximity to human settelement on wolf - coyote co-occurrence 2020")
  wc_dist2burbs_20s_facet
  
  wc_elev_20s_predict <- wc_elev_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "coyote" = "Coyote") 
  wc_elev_20s_facet <- ggplot(wc_elev_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Elevation (m)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of elevation on wolf - coyote co-occurrence 2020")
  wc_elev_20s_facet
  
  #'  2021 co-occurrence
  wc_bunnies_21s_predict <- wc_bunnies_21s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "coyote" = "Coyote") 
  wc_bunnies_21s_facet <- ggplot(wc_bunnies_21s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Lagomorph detections") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of lagomorph detections on wolf - coyote co-occurrence 2021")
  wc_bunnies_21s_facet
  
  wc_elev_21s_predict <- wc_elev_21s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "coyote" = "Coyote") 
  wc_elev_21s_facet <- ggplot(wc_elev_21s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Elevation (m)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of elevation on wolf - coyote co-occurrence 2021")
  wc_elev_21s_facet
  
  wc_bigdeer_21s_predict <- wc_bigdeer_21s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "coyote" = "Coyote") 
  wc_bigdeer_21s_facet <- ggplot(wc_bigdeer_21s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Elk and moose detections") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of elk and moose detections on wolf - coyote co-occurrence 2021")
  wc_bigdeer_21s_facet
  
  #####  Lion-Coyote Plots  #####
  #'  2020 co-occurrence
  lc_dist2burbs_20s_predict <- lc_dist2burbs_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "lion absent", "lion present")))
  newlabs <- c("lion" = "Lion", "coyote" = "Coyote") 
  lc_dist2burbs_20s_facet <- ggplot(lc_dist2burbs_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Distance to urban/suburban area (m)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of proximity to human settelement on wolf - coyote co-occurrence 2020")
  lc_dist2burbs_20s_facet
  
  lc_nearestrd_20s_predict <- lc_nearestrd_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "lion absent", "lion present")))
  newlabs <- c("lion" = "Lion", "coyote" = "Coyote") 
  lc_nearestrd_20s_facet <- ggplot(lc_nearestrd_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Log distance to nearest road") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of log distance to nearest road on lion - coyote co-occurrence 2020")
  lc_nearestrd_20s_facet
  
  lc_wtd_20s_predict <- lc_wtd_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "lion absent", "lion present")))
  newlabs <- c("lion" = "Lion", "coyote" = "Coyote") 
  lc_wtd_20s_facet <- ggplot(lc_wtd_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("White-tailed Deer Detections") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of white-tailed deer activity on wolf - coyote co-occurrence 2020")
  lc_wtd_20s_facet
  
  lc_elev_20s_predict <- lc_elev_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "lion absent", "lion present")))
  newlabs <- c("lion" = "Lion", "coyote" = "Coyote") 
  lc_elev_20s_facet <- ggplot(lc_elev_20s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Elevation (m)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of elevation on wolf - coyote co-occurrence 2020")
  lc_elev_20s_facet
  
  #'  2021 co-occurrence 
  lc_nearestrd_21s_predict <- lc_nearestrd_21s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "lion absent", "lion present")))
  newlabs <- c("lion" = "Lion", "coyote" = "Coyote") 
  lc_nearestrd_21s_facet <- ggplot(lc_nearestrd_21s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("lion absent" = "Lion absent", "lion present" = "Lion present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Log distance to nearest road") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of log distance to nearest road on lion - coyote co-occurrence 2021")
  lc_nearestrd_21s_facet
  
  #####  Bobcat-Coyote Plot  #####
  #'  2021 co-occurrence
  bc_human_21s_predict <- bc_human_21s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "bobcat absent", "bobcat present")))
  newlabs <- c("bobcat" = "Bobcat", "coyote" = "Coyote") 
  lc_human_20s_facet <- ggplot(bc_human_21s_predict, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present", 
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("bobcat absent" = "Bobcat absent", "bobcat present" = "Bobcat present", 
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Human Detections") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of human activity on bobcat - coyote co-occurrence 2020")
  lc_human_20s_facet
  