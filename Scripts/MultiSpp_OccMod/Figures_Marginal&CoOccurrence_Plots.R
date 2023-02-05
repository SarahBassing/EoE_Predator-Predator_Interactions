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
  #'  Predict cond. occupancy for each pairwise interaction over range of cattle activity
  #'  NOTE: setting all other continuous variables to 0 (their mean value) and 
  #'  setting categorical variables as either 0 or 1
  #'  For elk, study area = 0 (NE) because there were so few detections in OK
  #'  For md, wtd, and moose, study area = 1 (OK) b/c that is where the bulk of
  #'  cattle activity occurred and there were sufficient detections of these species
  #'  "pub" variable set to active grazing allotments (1)
  sppX_wolf_coy_20s <- spp_interactions_g(wc_20s_top, elev = 0, forest = 0, 
                                            dist2sub = 0, nearestrd = 0, elk = 0, 
                                            human = scaled_human[[1]][,2], bunny = 0, cow = 0, moose = 0, 
                                            md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                            spp1 = "wolf", spp2 = "coyote", cov = scaled_human[[1]][,1])
  
  sppX_wolf_coy_21s <- spp_interactions_g(wc_21s_top, elev = 0, forest = 0, 
                                          dist2sub = 0, nearestrd = 0, elk = 0, 
                                          human = scaled_human[[2]][,2], bunny = 0, cow = 0, moose = 0, 
                                          md = 0, wtd = 0, bigdeer = 0, smdeer = 0,
                                          spp1 = "wolf", spp2 = "coyote", cov = scaled_human[[2]][,1])
  
  
  #'  Plot co-occ
  tst <- sppX_wolf_coy_20s %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("coyote absent", "coyote present", "wolf absent", "wolf present")))
  newlabs <- c("wolf" = "Wolf", "coyote" = "Coyote") 
  tst_facet <- ggplot(tst, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
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
    xlab("Human activity (human detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of human activity on wolf - coyote co-occurrence")
  tst_facet
  