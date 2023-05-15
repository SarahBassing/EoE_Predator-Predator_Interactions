  #'  ---------------------------------
  #'  TBD figures and result tables
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  May 2023
  #'  ---------------------------------
  #'  Script to create figures and result tables based on results from analyses 
  #'  that estimated the effect of recent competitor presence and prey availability 
  #'  on the time between detections of predators.
  #'  ---------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(ggplot2)
  library(stringr)
  library(tidyverse)
  library(khroma)
  library(patchwork)
  library(grid)
  library(png)
  library(RCurl)
  library(rphylopic)
  
  #'  Load top models
  load("./Outputs/Time_btwn_Detections/tbd.comp.bear_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_X_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.lion_preyRAI.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.wolf_preyRAI.RData")
  
  #'  Load competitor models
  load("./Outputs/Time_btwn_Detections/tbd.comp.bear_competitor_detection.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.bob_competitor_detection.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.coy_competitor_detection.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.lion_competitor_detection.RData")
  load("./Outputs/Time_btwn_Detections/tbd.comp.wolf_competitor_detection.RData")
  
  #'  Load bundled data, including new covariate data
  load("./Data/Time_btwn_Detections/bear_bundled.RData")
  load("./Data/Time_btwn_Detections/bob_bundled.RData")
  load("./Data/Time_btwn_Detections/coy_bundled.RData")
  load("./Data/Time_btwn_Detections/lion_bundled.RData")
  load("./Data/Time_btwn_Detections/wolf_bundled.RData")
  
  #'  ------------------------
  ####  Species silhouettes  ####
  #'  ------------------------
  #'  Silhouettes for each species from PhyloPic in two different formats (PNG & rastergrob)
  #'  Cougar, mule deer, and white-tailed deer silhouettes created by the talented 
  #'  Gabriela Palomo-Munoz and uploaded to http://phylopic.org/
  cougurlGB <- "https://images.phylopic.org/images/cbe2a3c9-2c11-4f36-a51f-8a6c8de6a420/raster/1024x489.png?v=17cfd978c92.png"
  cougimgGB <- readPNG(getURLContent(cougurlGB), native = T)
  couggrid <- rasterGrob(cougimgGB, interpolate = TRUE)
  wolfurl <- "https://images.phylopic.org/images/8cad2b22-30d3-4cbd-86a3-a6d2d004b201/raster/1024x797.png?v=16ff245cc7c.png"
  wolfimg <- readPNG(getURLContent(wolfurl), native = T)
  wolfgrid <- rasterGrob(wolfimg, interpolate = TRUE)
  boburl <- "https://images.phylopic.org/images/ab6cfd4f-aef7-40fa-b5a5-1b79b7d112aa/raster/1024x740.png?v=17a30c18af1.png"
  bobimg <- readPNG(getURLContent(boburl), native = T)
  bobgrid <- rasterGrob(bobimg, interpolate = TRUE)
  coyurl <- "https://images.phylopic.org/images/5a0398e3-a455-4ca6-ba86-cf3f1b25977a/raster/1024x894.png?v=16fe8749858.png"
  coyimg <- readPNG(getURLContent(coyurl), native = T) 
  coygrid <- rasterGrob(coyimg, interpolate = TRUE)
  coyurlGB <- "https://images.phylopic.org/images/e6a2fa4b-85df-43b4-989c-34a65ba7eee3/raster/1024x911.png?v=17f2638df97.png"
  coyimgGB <- readPNG(getURLContent(coyurlGB), native = T)
  coygridGB <- rasterGrob(coyimg, interpolate = TRUE)
  bearurl <- "https://images.phylopic.org/images/369a7880-4798-41bf-851a-ec5da17fafa3/raster/1024x753.png?v=178afd80706.png"
  bearimg <- readPNG(getURLContent(bearurl), native = T)
  beargrid <- rasterGrob(bearimg, interpolate = TRUE)
  elkfurl <- "https://images.phylopic.org/images/97f83f5e-9afe-4ce8-812e-337f506ca841/raster/1024x1005.png?v=1402ea30c27.png"
  elkfimg <- readPNG(getURLContent(elkfurl), native = T)
  elkgrid <- rasterGrob(elkfimg, interpolate = TRUE)
  elkmurl <- "https://images.phylopic.org/images/72f2f997-e474-4caf-bbd5-72fc8dbcc40d/raster/866x1024.png?v=1356f9ea6de.png"
  elkmimg <- readPNG(getURLContent(elkmurl), native = T)
  elkmgrid <- rasterGrob(elkmimg, interpolate = TRUE)
  wtdurlGB1 <- "https://images.phylopic.org/images/8569838c-c725-4772-b0a3-b5eb04baaada/raster/1024x850.png?v=17cfdbaf920.png"
  wtdimgGB1 <- readPNG(getURLContent(wtdurlGB1), native = T)
  wtdgrid <- rasterGrob(wtdimgGB1, interpolate = TRUE)
  wtdurlGB2 <- "https://images.phylopic.org/images/6038e80c-398d-47b2-9a69-2b9edf436f64/raster/1023x1024.png?v=17cfdb9f8b6.png"
  wtdimgGB2 <- readPNG(getURLContent(wtdurlGB2), native = T)
  wtdgridGB2 <- rasterGrob(wtdimgGB2, interpolate = TRUE)
  bunnyurl <- "https://images.phylopic.org/images/f69eb95b-3d0d-491d-9a7f-acddd419afed/raster/925x1024.png?v=177f427b3d8.png"
  bunnyimg <- readPNG(getURLContent(bunnyurl), native = T)
  bunnygrid <- rasterGrob(bunnyimg, interpolate = TRUE)
  
  bear <- get_phylopic("369a7880-4798-41bf-851a-ec5da17fafa3")
  
  
  
  #'  ------------------
  ####  Format results  ####
  #'  ------------------
  #'  Snag and reformat coefficents and predictions from each top model
  coefs <- function(mod_out, spp, prey1, prey2, prey3, comp1, comp2, comp3, comp4) {
    Species <- spp
    Estimate <- round(unlist(mod_out$mean), 2)
    lci <- round(unlist(mod_out$q2.5), 2)
    uci <- round(unlist(mod_out$q97.5), 2)
    CI <- paste(" ", lci, "-", uci) # need that extra space in front b/c excel thinks this is an equation otherwise
    overlap0 <- unlist(mod_out$overlap0)
    out <- as.data.frame(cbind(Species, Estimate, CI, lci, uci, overlap0))
    out <- tibble::rownames_to_column(out, "row_names") %>%
      relocate(row_names, .after = Species)
    colnames(out) <- c("Species", "Parameter", "Estimate", "95% CI", "lci", "uci", "overlap0")
    renamed_out <- out %>%
      mutate(Parameter = ifelse(Parameter == "alpha0", "Intercept", Parameter),
             Parameter = ifelse(Parameter == "beta.prey1", paste("Prey RAI:", prey1), Parameter),
             Parameter = ifelse(Parameter == "beta.prey2", paste("Prey RAI:", prey2), Parameter),
             Parameter = ifelse(Parameter == "beta.prey3", paste("Prey RAI:", prey3), Parameter),
             Parameter = ifelse(Parameter == "beta.competitor1", paste("Competitor:", comp1), Parameter), 
             Parameter = ifelse(Parameter == "beta.competitor2", paste("Competitor:", comp2), Parameter), 
             Parameter = ifelse(Parameter == "beta.competitor3", paste("Competitor:", comp3), Parameter), 
             Parameter = ifelse(Parameter == "beta.competitor4", paste("Competitor:", comp4), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd1", paste("Mean TBD:", comp1), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd2", paste("Mean TBD:", comp2), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd3", paste("Mean TBD:", comp3), Parameter),
             Parameter = ifelse(Parameter == "spp.tbd4", paste("Mean TBD:", comp4), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd1", paste0("Competitor:Prey (", comp1, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd2", paste0("Competitor:Prey (", comp2, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd3", paste0("Competitor:Prey (", comp3, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.wtd4", paste0("Competitor:Prey (", comp4, " x ", prey1, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago1", paste0("Competitor:Prey (", comp1, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago2", paste0("Competitor:Prey (", comp2, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago3", paste0("Competitor:Prey (", comp3, " x ", prey2, ")"), Parameter),
             Parameter = ifelse(Parameter == "beta.interaction.lago4", paste0("Competitor:Prey (", comp4, " x ", prey2, ")"), Parameter),
             Parameter = gsub("spp.tbd.", "tbd ", Parameter),
             Parameter = ifelse(Parameter == "mu.tbd", "Mean TBD", Parameter)) %>%
      filter(lci != 0) %>%
      mutate(Estimate = as.numeric(Estimate),
             lci = as.numeric(lci),
             uci = as.numeric(uci))
    return(renamed_out)
  }
  tbd.bear.out <- coefs(tbd.bear.preyabund, spp = "Black bear", prey1 = "Elk", prey2 = "White-tailed deer")
  tbd.bob.out <- coefs(tbd.bob.compID.preyabund, spp = "Bobcat", prey1 = "White-tailed deer", prey2 = "Lagomorph", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.coy.out <- coefs(tbd.coy.compIDxpreyabund, spp = "Coyote", prey1 = "White-tailed deer", prey2 = "Lagomorph", comp1 = "Black bear", comp2 = "Bobcat", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.lion.out <- coefs(tbd.lion.preyabund, spp = "Mountain lion", prey1 = "Elk", prey2 = "White-tailed deer")
  tbd.wolf.out <- coefs(tbd.wolf.preyabund, spp = "Wolf", prey1 = "Elk", prey2 = "Moose", prey3 = "White-tailed deer")
  
  tbd.bear.comp <- coefs(tbd.bear.compID, spp = "Black bear", comp1 = "Coyote", comp2 = "Bobcat", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.bob.comp <- coefs(tbd.bob.compID, spp = "Bobcat", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.coy.comp <- coefs(tbd.coy.compID, spp = "Coyote", comp1 = "Black bear", comp2 = "Bobcat", comp3 = "Mountain lion", comp4 = "Wolf")
  tbd.lion.comp <- coefs(tbd.lion.compID, spp = "Mountain lion", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Bobcat", comp4 = "Wolf")
  tbd.wolf.comp <- coefs(tbd.wolf.compID, spp = "Wolf", comp1 = "Coyote", comp2 = "Black bear", comp3 = "Bobcat", comp4 = "Mountain lion")
  
  #'  Pull out just coefficient estimates
  bear.coefs <- tbd.bear.out[1:3,]
  bob.coefs <- tbd.bob.out[1:6,]
  coy.coefs <- tbd.coy.out[1:12,]
  lion.coefs <- tbd.lion.out[1:3,]
  wolf.coefs <- tbd.wolf.out[1:4,]
  
  #'  Pull out mean TBD estimates
  bear.mean.tbd <- tbd.bear.out[4,1:6]
  bob.mean.tbd <- tbd.bob.out[7:11,1:6]
  coy.mean.tbd <- tbd.coy.out[13:17,1:6]
  lion.mean.tbd <- tbd.lion.out[4,1:6]
  wolf.mean.tbd <- tbd.wolf.out[5,1:6]
  
  bear.comp.tbd <- tbd.bear.comp[5:9,1:6]
  bob.comp.tbd <- tbd.bob.out[7:11,1:6] #tbd.bob.comp[5:9,1:6]
  coy.comp.tbd <- tbd.coy.out[13:17,1:6] #tbd.coy.comp[5:9,1:6]
  lion.comp.tbd <- tbd.lion.comp[5:9,1:6]
  wolf.comp.tbd <- tbd.wolf.comp[5:9,1:6]
  
  #'  Pull out predicted TBD values
  bear.tbd.predictions <- tbd.bear.out[5:204,1:6]
  bob.tbd.predictions <- tbd.bob.out[12:811,1:6]
  coy.tbd.predictions <- tbd.coy.out[18:817,1:6]
  lion.tbd.predictions <- tbd.lion.out[5:204,1:6]
  wolf.tbd.predictions <- tbd.wolf.out[6:305,1:6]
  
  
  #'  ------------------
  ####  Result tables  ####
  #'  -----------------
  tbd.coefs <- rbind(bear.coefs, bob.coefs, coy.coefs, lion.coefs, wolf.coefs)
  mean.tbd <- rbind(bear.mean.tbd, bob.mean.tbd, coy.mean.tbd, lion.mean.tbd, wolf.mean.tbd)
  competitor.tbd <- rbind(bear.comp.tbd, bob.comp.tbd, coy.comp.tbd, lion.comp.tbd, wolf.comp.tbd)
  predicted.tbd <- rbind(bear.tbd.predictions, bob.tbd.predictions, coy.tbd.predictions, lion.tbd.predictions, wolf.tbd.predictions) %>%
    #'  Split out predictions based on categorical variable (mostly important for bobcat & coyote results)
    mutate(Prey_species = str_replace(Parameter, "tbd ", ""),
           Prey_species = str_extract(Prey_species, "[aA-zZ]+"),
           Obs_nmbr = as.numeric(str_extract(Parameter, "[0-9]+")),
           #'  If competitor ID had no effect, use reference category (coyote)
           Competitor_ID = ifelse(Species == "Black bear" | Species == "Mountain lion" | Species == "Wolf", "Coyote", NA),
           #'  If competitor ID had an effect, assign correct species to each data chunk
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr <101, "Coyote", Competitor_ID),
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr >100, "Black bear", Competitor_ID),
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr >200, "Mountain lion", Competitor_ID), 
           Competitor_ID = ifelse(Species == "Bobcat" & Obs_nmbr >300, "Wolf", Competitor_ID),
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr <101, "Black bear", Competitor_ID),
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr >100, "Bobcat", Competitor_ID),
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr >200, "Mountain lion", Competitor_ID), 
           Competitor_ID = ifelse(Species == "Coyote" & Obs_nmbr >300, "Wolf", Competitor_ID))
  
  #' #'  Save
  #' write.csv(tbd.coefs, "./Outputs/Tables/TBD_coefficient_estimates_allSpp.csv")
  #' write.csv(mean.tbd, "./Outputs/Tables/TBD_estimated_means_allSpp.csv")
  #' write.csv(competitor.tbd, "./Outputs/Tables/TBD_competitor_means_allSpp.csv")
  #' write.csv(predicted.tbd, "./Outputs/Tables/TBD_predicted_TBD_allSpp.csv")
  
  
  #'  ---------------------------------
  ####  Plot TBD by covariate effects  ####
  #'  ---------------------------------
  #'  Format new covariate data based on length of predictions for each predator species
  bear.covs <- c(bear_bundled$newcovs[,1], bear_bundled$newcovs[,3]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("elk", "wtd"), each = 100)) %>%
    rename("cov" = ".")
  bob.covs <- c(bob_bundled$newcovs[,3], bob_bundled$newcovs[,3], bob_bundled$newcovs[,3], bob_bundled$newcovs[,3], bob_bundled$newcovs[,4], bob_bundled$newcovs[,4], bob_bundled$newcovs[,4], bob_bundled$newcovs[,4]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("wtd", "lago"), each = 400)) %>%
    rename("cov" = ".")
  coy.covs <- c(coy_bundled$newcovs[,3], coy_bundled$newcovs[,3], coy_bundled$newcovs[,3], coy_bundled$newcovs[,3], coy_bundled$newcovs[,4], coy_bundled$newcovs[,4], coy_bundled$newcovs[,4], coy_bundled$newcovs[,4]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("wtd", "lago"), each = 400)) %>%
    rename("cov" = ".")
  lion.covs <- c(lion_bundled$newcovs[,1], lion_bundled$newcovs[,3]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("elk", "wtd"), each = 100)) %>%
    rename("cov" = ".")
  wolf.covs <- c(wolf_bundled$newcovs[,1], wolf_bundled$newcovs[,2], wolf_bundled$newcovs[,3]) %>%
    as.data.frame() %>%
    mutate(prey_spp = rep(c("elk", "moose", "wtd"), each = 100)) %>%
    rename("cov" = ".")
  
  #'  Review which coefficients did NOT overlap 0 - only retain predictions for
  #'  responses to prey RAI based on these relationships
  print(bear.coefs) # Prey RAI: Elk
  print(bob.coefs) # Prey RAI: White-tailed deer
  print(coy.coefs) # Competitor:Prey (Bobcat x White-tailed deer), (Mountain lion x White-tailed deer), (Wolf x White-tailed deer), (Wolf x Lagomorph)
  print(lion.coefs) # Prey RAI: Elk, Prey RAI: White-tailed deer
  print(wolf.coefs) # Prey RAI: White-tailed deer
  
  #'  Combine predictions and covariate data for plotting
  #'  Filter to only predictions derived from statistically meaningful relationships
  #'  Change Competitor_ID to "any" for species where previous competitor detection did not matter
  bear.predicted <- filter(predicted.tbd, Species == "Black bear") %>%
    bind_cols(bear.covs) %>%
    dplyr::select(-prey_spp) %>%
    filter(Prey_species == "elk") %>%
    mutate(Competitor_ID = "Any species")
  bob.predicted <- filter(predicted.tbd, Species == "Bobcat") %>%
    bind_cols(bob.covs) %>%
    dplyr::select(-prey_spp) %>%
    filter(Prey_species == "wtd")
  coy.predicted <- filter(predicted.tbd, Species == "Coyote") %>%
    bind_cols(coy.covs) %>%
    dplyr::select(-prey_spp) %>%
    filter(Prey_species != "lago" | Competitor_ID != "Bobcat") %>%
    filter(Prey_species != "lago" | Competitor_ID != "Mountain lion")
  lion.predicted <- filter(predicted.tbd, Species == "Mountain lion") %>%
    bind_cols(lion.covs) %>%
    dplyr::select(-prey_spp) %>%
    mutate(Competitor_ID = "Any species")
  wolf.predicted <- filter(predicted.tbd, Species == "Wolf") %>%
    bind_cols(wolf.covs) %>%
    dplyr::select(-prey_spp) %>%
    filter(Prey_species == "wtd") %>%
    mutate(Competitor_ID = "Any species")
  
  predicted.tbd.covs <- rbind(bear.predicted, bob.predicted, coy.predicted, lion.predicted, wolf.predicted)

  #'  Format results tables for plotting
  tbd_by_competitorID <- competitor.tbd %>%
    filter(Parameter != "Mean TBD") %>%
    arrange(Species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Previous_Species = str_replace(Parameter, "Mean TBD: ", ""),
           Previous_Species = factor(Previous_Species, levels = c("Coyote", "Bobcat", "Black bear", "Mountain lion", "Wolf")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) %>%
    relocate(Previous_Species, .after = Species)
  
  tbd_prey_prediction <- predicted.tbd.covs %>%
    # arrange(Species, Prey_species, Estimate) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Competitor_ID = factor(Competitor_ID, levels = c("Any species", "Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Prey_RAI = Prey_species,
           Prey_RAI = factor(Prey_RAI, levels = c("Elk" = "elk", "Moose" = "moose", "White-tailed deer" = "wtd", "Lagomorphs" = "lago")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) 
 
  #'  Choose colorblind-friendly scheme
  plot_scheme(colour("sunset")(11))
  colour("sunset")(11)
  # any.wolf.bear.lion.coy.bob_colors <- c("black", "#364B9A", "#98CAE1", "#DD3D2D", "#FDB366", "#A50026")
  any.bear.bob.coy.lion.wolf_colors <- c("black", "#98CAE1", "#A50026", "#DD3D2D", "#FDB366", "#364B9A")
  bear.wolf <- c("#98CAE1", "#364B9A")
  
  #####  Competitor ID effect  ####
  #'  Effect of previous competitor detection on latency of site use
  competitorID_plot <- ggplot(tbd_by_competitorID, aes(x = Previous_Species, y = (Estimate/60), group = Species)) +
    geom_errorbar(aes(ymin = (lci/60), ymax = (uci/60), color = Previous_Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Previous_Species), size = 2.5, position = position_dodge(width = 0.4)) +
    theme_bw() +
    guides(color = guide_legend(title = "Previously detected species")) +
    facet_wrap(~Species, scales = "free_y") + #, scales = "free_y"
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
    xlab("Previously detected species") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent competitor detection on latency of site use")
  competitorID_plot
  
  #' #'  Response of all species to the recent detection of a specific competitor
  #' competitorID_plot_flipped <- ggplot(tbd_by_competitorID, aes(x = Species, y = (Estimate/60), group = Previous_Species)) +
  #'   geom_errorbar(aes(ymin = (lci/60), ymax = (uci/60), color = Species), width = 0, position = position_dodge(width = 0.4)) +
  #'   geom_point(stat = "identity", aes(col = Species), size = 2.5, position = position_dodge(width = 0.4)) +
  #'   theme_bw() +
  #'   guides(color = guide_legend(title = "Following species detected")) +
  #'   facet_wrap(~Previous_Species) + #, scales = "free_y"
  #'   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  #'   theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
  #'   xlab("Species detected following detection of a competitor") +
  #'   ylab("Mean number of hours between detections") +
  #'   ggtitle("Competitor effect on latency of detection of following species")
  #' competitorID_plot_flipped
  
  
  #####  Prey relative abundance effects  ####
  tbd_wtdRAI_plot <- filter(tbd_prey_prediction, Prey_RAI == "wtd") %>%
    ggplot(aes(x = cov, y = (Estimate/60), group = Competitor_ID)) +
    geom_line(aes(color = Competitor_ID), lwd = 1.25) + 
    scale_color_manual(values = any.bear.bob.coy.lion.wolf_colors) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = (lci/60), ymax = (uci/60), fill = Competitor_ID), alpha = 0.3) +
    scale_fill_manual(values = any.bear.bob.coy.lion.wolf_colors) + 
    theme_bw() +
    xlab("Relative abundance of white-tailed deer (scaled)") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent competitor detection and white-tailed deer abundance on\nlatency of site use") +
    guides(color = guide_legend(title = "Previously detected\ncompetitor"),
           fill = guide_legend(title = "Previously detected\ncompetitor")) +
    facet_wrap(~Species, scale = "free_y")
  tbd_wtdRAI_plot
  
  tbd_lagoRAI_plot <- filter(tbd_prey_prediction, Prey_RAI == "lago") %>%
    ggplot(aes(x = cov, y = (Estimate/60), group = Competitor_ID)) +
    geom_line(aes(color = Competitor_ID), lwd = 1.25) + 
    scale_color_manual(values = bear.wolf) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = (lci/60), ymax = (uci/60), fill = Competitor_ID), alpha = 0.3) +
    scale_fill_manual(values = bear.wolf) + 
    theme_bw() +
    xlab("Relative abundance of lagomorphs (scaled)") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent competitor detection and\nlagomorph abundance on latency of site use") +
    guides(color = guide_legend(title = "Previously detected\ncompetitor"),
           fill = guide_legend(title = "Previously detected\ncompetitor")) +
    facet_wrap(~Species, scale = "free_y")
  tbd_lagoRAI_plot
  
  tbd_elkRAI_plot <- filter(tbd_prey_prediction, Prey_RAI == "elk") %>%
    ggplot(aes(x = cov, y = (Estimate/60), group = Competitor_ID)) +
    geom_line(aes(color = Competitor_ID), lwd = 1.25) + 
    scale_color_manual(values = any.bear.bob.coy.lion.wolf_colors) +
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = (lci/60), ymax = (uci/60), fill = Competitor_ID), alpha = 0.3) +
    scale_fill_manual(values = any.bear.bob.coy.lion.wolf_colors) + 
    theme_bw() +
    xlab("Relative abundance of elk (scaled)") +
    ylab("Mean number of hours between detections") +
    ggtitle("Effect of recent competitor detection and elk abundance on latency of site use") +
    # theme(legend.position = "none") +
    guides(color = guide_legend(title = "Previously detected\ncompetitor"),
           fill = guide_legend(title = "Previously detected\ncompetitor")) +
    facet_wrap(~Species, scale = "free_x")
  tbd_elkRAI_plot

  
  #'  Save
  ggsave("./Outputs/Time_btwn_Detections/Figures/Mean_TBD_by_competitorID.tiff", competitorID_plot, 
         units = "in", width = 7, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Time_btwn_Detections/Figures/Predicted_TBD_by_wtdRAI.tiff", tbd_wtdRAI_plot, 
         units = "in", width = 7, height = 7, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Time_btwn_Detections/Figures/Predicted_TBD_by_lagomorphRAI.tiff", tbd_lagoRAI_plot, 
         units = "in", width = 5, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Time_btwn_Detections/Figures/Predicted_TBD_by_elkRAI.tiff", tbd_elkRAI_plot, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  -----------------------------------
  #####  Plot coefficient effect sizes  ####
  #'  -----------------------------------
  #'  Format coefficients for plotting
  parameter_est <- tbd.coefs %>%
    arrange(Species) %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Bobcat", "Coyote", "Mountain lion", "Wolf")),
           Parameter = factor(Parameter, levels =  c("Competitor:Prey (Wolf x Lagomorph)", "Competitor:Prey (Wolf x White-tailed deer)", 
                                                     "Competitor:Prey (Mountain lion x Lagomorph)", "Competitor:Prey (Mountain lion x White-tailed deer)", 
                                                     "Competitor:Prey (Bobcat x Lagomorph)", "Competitor:Prey (Bobcat x White-tailed deer)", 
                                                     "Prey RAI: Lagomorph", "Prey RAI: White-tailed deer", "Prey RAI: Moose", "Prey RAI: Elk",
                                                     "Competitor: Wolf", "Competitor: Mountain lion", "Competitor: Coyote", "Competitor: Bobcat", 
                                                     "Competitor: Black bear", "Intercept")),
           Estimate = as.numeric(Estimate),
           lci = as.numeric(lci),
           uci = as.numeric(uci)) 
  
  #####  Coefficient effect sizes  ####
  meso_coef_plot <- filter(parameter_est, Species == "Bobcat" | Species == "Coyote") %>%
    mutate(Species = factor(Species, levels = c("Coyote", "Bobcat"))) %>%
    filter(Parameter != "Intercept") %>%  ############### NOTE: Excluding intercept so 95% CRIs of other coeffs are more visible!!!
    ggplot(aes(x = Parameter, y = (Estimate/60))) +
    geom_errorbar(aes(ymin = (lci/60), ymax = (uci/60), color = Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Species), size = 2.5, position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_bw() +
    coord_flip() +
    facet_wrap(~Species, scales = "free_x") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    xlab("Parameter") +
    ylab("Estimated coefficient (95% CRI)") +
    ggtitle("Coefficient estimates for bobcat and coyote analyses") #+
    #annotation_custom(beargrid, xmin=11.75, xmax=12.5) 
  meso_coef_plot
  
  coy_coef_plot <- filter(parameter_est, Species == "Coyote") %>%
    filter(Parameter != "Intercept") %>%  ############### NOTE: Excluding intercept so 95% CRIs of other coeffs are more visible!!!
    ggplot(aes(x = Parameter, y = (Estimate/60))) +
    geom_errorbar(aes(ymin = (lci/60), ymax = (uci/60), color = Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Species), size = 2.5, position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    xlab("Parameter") +
    ylab("Estimated coefficient (95% CRI)") +
    ggtitle("Coefficient estimates for bobcat and coyote analyses") +
    add_phylopic(wolfimg, x = 7.05, y = 0.3, ysize = 0.5, color = "black", alpha = 1)
  coy_coef_plot
    
    add_phylopic(wolfimg, x = 7.05, y = 0.8, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimgGB1, x = 6.1, y = 0.8, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(bunnyimg, x = 5.05, y = 0.8, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkfimg, x = 4.05, y = 0.8, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimgGB, x = 3.05, y = 0.8, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimgGB, x = 2, y = 0.8, ysize = 0.45, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.8, ysize = 0.4, color = "black", alpha = 1) +
    add_phylopic(bearimg, x = 1.05, y = 0.8, ysize = 0.4, color = "black", alpha = 1) +
    add_phylopic(wolfimg, x = 1.05, y = 0.8, ysize = 0.4, color = "black", alpha = 1)
  meso_coef_plot
  
  apex_coef_plot <- filter(parameter_est, Species == "Black bear" | Species == "Mountain lion" | Species == "Wolf") %>%
    mutate(Species = factor(Species, levels = c("Black bear", "Mountain lion", "Wolf"))) %>%
    filter(Parameter != "Intercept") %>%  ############### NOTE: Excluding intercept so 95% CRIs of other coeffs are more visible!!!
    ggplot(aes(x = Parameter, y = (Estimate/60))) +
    geom_errorbar(aes(ymin = (lci/60), ymax = (uci/60), color = Species), width = 0, position = position_dodge(width = 0.4)) +
    geom_point(stat = "identity", aes(col = Species), size = 2.5, position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_bw() +
    coord_flip() +
    facet_wrap(~Species, scales = "free_x") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none") +
    xlab("Parameter") +
    ylab("Estimated coefficient (95% CRI)") +
    ggtitle("Coefficient estimates for black bear, mountain lion, and wolf analyses")
  apex_coef_plot
  
  
  #'  Save
  ggsave("./Outputs/Time_btwn_Detections/Figures/Coefficient_size_mesopredators.tiff", meso_coef_plot, 
         units = "in", width = 7, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Time_btwn_Detections/Figures/Coefficient_size_largepredators.tiff", apex_coef_plot, 
         units = "in", width = 7, height = 7, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  
  