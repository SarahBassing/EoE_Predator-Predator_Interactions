  #'  ------------------------------
  #'  Plot co-detection results
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  -------------------------------
  #'  Script to create figures of predicted co-detection and covariate effects
  #'  from multispecies occupancy models run in JAGS.
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(ggplot2)
  library(khroma)
  library(stringr)
  library(patchwork)
  library(rphylopic)
  library(tidyverse)
  
  #'  Load covariate data
  # load("./Data/Covariates_extracted/Covariate_skinny_EoE20s21s.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Format_data_2spp_occmod_for_JAGS_img_updated_072924.RData")
  
  #' #'  Load covariate data
  #' load("./Data/Covariates_extracted/Covariate_skinny_EoE20s21s_updated_070824.RData")
  
  #'  Identify top models
  load("./Outputs/Tables/DIC_top_models.RData")
  print(topmodels)
  
  #'  -------------------------------------------------
  ####  Predict Pr(detection) across covariate values  ####
  #'  -------------------------------------------------
  #'  MASSIVE function the predict marginal and conditional detection for each species
  predict_detection <- function(mod, ncat, npoints, focal_cov, rho_cov, rho_inxs_cov, rho_cov_index, rho_inxs_cov_index) {
    #'  Rename model and review
    out1 <- mod
    #print(out1$summary[1:50,])
    #'  Number of states
    ncat <- ncat
    #'  Number of draws in the MCMC chain
    ndraws <- out1$mcmc.info$n.samples
    #'  Number of observations to predict
    npoints <- npoints
    
    #'  Create range of covariate values to predict across
    print(r <- range(focal_cov))
    print(mean.cov <- mean(focal_cov))
    print(sd.cov <- sd(focal_cov))
    cov.pred.orig <- seq(r[1], r[2], length.out = npoints)
    #'  Scale new data to be consistent with data used in model
    cov.pred <- (cov.pred.orig - mean.cov) / sd.cov
    
    #'  Create model matrices for prediction covariates on rho and rhox
    #'  Fill matrix with categorical or mean value of each covariate
    rho_cov <- rho_cov
    rho_covs <- matrix(rho_cov, nrow = npoints, ncol = length(rho_cov), byrow = TRUE)
    #'  Same thing but for covariates on interaction
    rho_inxs_cov <- rho_inxs_cov
    rho_inxs_covs <- matrix(rho_inxs_cov, nrow = npoints, ncol = length(rho_inxs_cov), byrow = TRUE)
    #'  Replace focal covariate column with new scaled data
    #'  New data must be in correct column based on order of covariates in original model
    rho_covs[,rho_cov_index] <- cov.pred
    print(head(rho_covs))
    rho_inxs_covs[,rho_inxs_cov_index] <- cov.pred
    print(head(rho_inxs_covs))
    
    #'  Assemble linear predictors across all iterations of the MCMC chains
    rhoSpp1 <- out1$sims.list$alphaSpp1 %*% t(rho_covs)
    rhoSpp2 <- out1$sims.list$alphaSpp2 %*% t(rho_covs)
    rhoSpp12 <- rhoSpp1 + rhoSpp2 + out1$sims.list$alphaSpp12 %*% t(rho_inxs_covs) 
    
    #'  Compute rho detection matrix (rdm) based on new data and means of other covariates
    rdm <- array(NA, dim = c(ndraws, npoints, ncat))
    rdm[,,1] <- 1
    rdm[,,2] <- exp(rhoSpp1)
    rdm[,,3] <- exp(rhoSpp2)
    rdm[,,4] <- exp(rhoSpp12)
    
    #'  Compute rho detection probability of each state from MCMC chains (rdp)
    #'  denominator = sum of rdm for each draw x point
    denom <- apply(rdm, 1:2, sum)
    rdpMC <- sweep(rdm, 1:2, denom, "/")
    dimnames(rdpMC)[[3]] <- c("U", "Spp1", "Spp2", "Both")
    str(rdpMC)
    
    #'  Fill in summary values for 'rdp'
    rdp <- array(NA, dim = c(npoints, ncat, 4))
    dimnames(rdp) <- list(NULL, c("U", "Spp1", "Spp2", "Both"),
                          c("mean", "SD", "lower", "upper"))
    rdp[,,1] <- apply(rdpMC, 2:3, mean)
    rdp[,,2] <- apply(rdpMC, 2:3, sd)
    rdp[,,3] <- apply(rdpMC, 2:3, quantile, probs = 0.025)
    rdp[,,4] <- apply(rdpMC, 2:3, quantile, probs = 0.975)
    
    #'  Compute and plot marginal probabilities
    #'  --------------------------------------- 
    #'  Occupancy probability for each species, regardless of presence or
    #'  absence of other species
    #'  Get MCMC chains
    marginalMC <- array(NA, dim = c(ndraws, npoints, 2))
    dimnames(marginalMC)[[3]] <- c("Spp1", "Spp2")
    marginalMC[,,"Spp1"] <- apply(rdpMC[,,c("Spp1", "Both")], 1:2, sum)
    marginalMC[,,"Spp2"] <- apply(rdpMC[,,c("Spp2", "Both")], 1:2, sum)
    
    #'  Fill in summary values for marginal probabilities
    marginal <- array(NA, dim = c(npoints, 2, 4))
    dimnames(marginal) <- list(NULL, c("Spp1", "Spp2"),
                               c("mean", "sd", "lower", "upper"))
    marginal[,,1] <- apply(marginalMC, 2:3, mean)
    marginal[,,2] <- apply(marginalMC, 2:3, sd)
    marginal[,,3] <- apply(marginalMC, 2:3, quantile, probs = 0.025)
    marginal[,,4] <- apply(marginalMC, 2:3, quantile, probs = 0.975)
    
    #'  Compute conditional detection for all species
    #'  --------------------------------------------- 
    #'  Detection probability for each species, conditional on the detection or
    #'  non-detection of the other species
    #'  Get MCMC chains
    conditionalMC <- array(NA, dim = c(ndraws, npoints, 4))
    dimnames(conditionalMC)[[3]] <- c("Spp1.alone", "Spp1.given.Spp2",
                                      "Spp2.alone", "Spp2.given.Spp1")
    
    conditionalMC[,,"Spp1.alone"] <- rdpMC[,,"Spp1"]/apply(rdpMC[,,c("U", "Spp1")], 1:2, sum)
    conditionalMC[,,"Spp1.given.Spp2"] <- rdpMC[,,"Both"]/marginalMC[,,"Spp2"]
    conditionalMC[,,"Spp2.alone"] <- rdpMC[,,"Spp2"]/apply(rdpMC[,,c("U", "Spp2")], 1:2, sum)
    conditionalMC[,,"Spp2.given.Spp1"] <- rdpMC[,,"Both"]/marginalMC[,,"Spp1"]
    
    #'  Fill in summary values
    conditional <- array(NA, dim = c(npoints, 4, 4))
    dimnames(conditional) <- list(NULL, c("Spp1.alone", "Spp1.given.Spp2",
                                          "Spp2.alone", "Spp2.given.Spp1"),
                                  c("mean", "sd", "lower", "upper"))
    conditional[,,1] <- apply(conditionalMC, 2:3, mean)
    conditional[,,2] <- apply(conditionalMC, 2:3, sd)
    conditional[,,3] <- apply(conditionalMC, 2:3, quantile, probs = 0.025)
    conditional[,,4] <- apply(conditionalMC, 2:3, quantile, probs = 0.975)
    
    #'  List predicted probabilities and new covariate data to plot
    predicted_probabilities <- list(marginal, conditional, cov.pred.orig, cov.pred)
    
    return(predicted_probabilities)
  }
  #####  Wolf-Black bear predictions  ####
  #'  Second most supported model (1.71 deltaDIC away from top)
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_px(.)_2024-08-29.RData")
  wolf.bear.eff.ung.yr1 <- predict_detection(mod = wolf.bear.hab.px, ncat = 4, npoints = 500,
                                             focal_cov = effort_eoe20s21s,
                                             rho_cov = c(1, 0, 0), rho_cov_index = 3,
                                             rho_inxs_cov = 0, rho_inxs_cov_index = 0)
  wolf.bear.eff.pred.yr1 <- predict_detection(mod = wolf.bear.hab.px, ncat = 4, npoints = 500,
                                              focal_cov = effort_eoe20s21s,
                                              rho_cov = c(1, 1, 0), rho_cov_index = 3,
                                              rho_inxs_cov = 0, rho_inxs_cov_index = 0)
  
  #####  Wolf-Lion predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2024-08-29.RData")  
  
  
  #####  Wolf-Coyote predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_px(.)_2024-08-29.RData") 
  wolf.coy.eff.ung.yr1 <- predict_detection(mod = wolf.coy.hab.px, ncat = 4, npoints = 500,
                                            focal_cov = effort_eoe20s21s,
                                            rho_cov = c(1, 0, 0), rho_cov_index = 3,
                                            rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  wolf.coy.eff.pred.yr1 <- predict_detection(mod = wolf.coy.hab.px, ncat = 4, npoints = 500,
                                             focal_cov = effort_eoe20s21s,
                                             rho_cov = c(1, 1, 0), rho_cov_index = 3,
                                             rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  wolf.coy.ung.yr1 <- predict_detection(mod = wolf.coy.hab.px, ncat = 4, npoints = 500,
                                        focal_cov = effort_eoe20s21s,
                                        rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                        rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  wolf.coy.pred.yr1 <- predict_detection(mod = wolf.coy.hab.px, ncat = 4, npoints = 500,
                                         focal_cov = effort_eoe20s21s,
                                         rho_cov = c(1, 1, 0), rho_cov_index = 0,
                                         rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  
  #####  Lion-Bear predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_px(.)_2024-08-30.RData")  
  
  
  #####  Lion-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_px(.)_2024-08-30.RData")  
  lion.bob.ung.yr1 <- predict_detection(mod = lion.bob.null.px, ncat = 4, npoints = 500,
                                        focal_cov = effort_eoe20s21s,
                                        rho_cov = c(1), rho_cov_index = 0,
                                        rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  lion.bob.pred.yr1 <- predict_detection(mod = lion.bob.null.px, ncat = 4, npoints = 500,
                                         focal_cov = effort_eoe20s21s,
                                         rho_cov = c(1), rho_cov_index = 0,
                                         rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  
  #####  Coyote-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_2024-08-30.RData")
  coy.bob.eff.ung.yr1 <- predict_detection(mod = coy.bob.habx.px, ncat = 4, npoints = 500,
                                           focal_cov = effort_eoe20s21s,
                                           rho_cov = c(1, 0, 0), rho_cov_index = 3,
                                           rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  coy.bob.eff.pred.yr1 <- predict_detection(mod = coy.bob.habx.px, ncat = 4, npoints = 500,
                                            focal_cov = effort_eoe20s21s,
                                            rho_cov = c(1, 1, 0), rho_cov_index = 3,
                                            rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  coy.bob.ung.yr1 <- predict_detection(mod = coy.bob.habx.px, ncat = 4, npoints = 500,
                                       focal_cov = effort_eoe20s21s,
                                       rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                       rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  coy.bob.pred.yr1 <- predict_detection(mod = coy.bob.habx.px, ncat = 4, npoints = 500,
                                        focal_cov = effort_eoe20s21s,
                                        rho_cov = c(1, 1, 0), rho_cov_index = 0,
                                        rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  
  #####  Black bear-Coyote predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_px(.)_2024-08-30.RData")
  bear.coy.eff.ung.yr1 <- predict_detection(mod = bear.coy.habx.px, ncat = 4, npoints = 500,
                                           focal_cov = effort_eoe20s21s,
                                           rho_cov = c(1, 0, 0), rho_cov_index = 3,
                                           rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  bear.coy.eff.pred.yr1 <- predict_detection(mod = bear.coy.habx.px, ncat = 4, npoints = 500,
                                            focal_cov = effort_eoe20s21s,
                                            rho_cov = c(1, 1, 0), rho_cov_index = 3,
                                            rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  bear.coy.ung.yr1 <- predict_detection(mod = bear.coy.habx.px, ncat = 4, npoints = 500,
                                       focal_cov = effort_eoe20s21s,
                                       rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                       rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  bear.coy.pred.yr1 <- predict_detection(mod = bear.coy.habx.px, ncat = 4, npoints = 500,
                                        focal_cov = effort_eoe20s21s,
                                        rho_cov = c(1, 1, 0), rho_cov_index = 0,
                                        rho_inxs_cov = c(1), rho_inxs_cov_index = 0)
  
  
  
  save.image(file = paste0("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/Predicted_rho-cov_relationships_", Sys.Date(), ".RData"))
  
  #'  --------------------------------
  ####  Plot marginal Pr(detection)  ####
  #'  --------------------------------
  load("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/Predicted_rho-cov_relationships_2024-08-30.RData")
  
  #'  Color-blind friendly color palette from Khroma
  plot_scheme(colour("sunset")(11))
  colour("sunset")(11)
  # two_colors <- c("#364B9A", "#A50026")
  two_colors <- c("#6EA6CD", "#DD3D2D")
  
  #'  Function to reformat data for easier use with ggplot2 and plot marginal detection
  plot_marginal_det <- function(predicted, spp1, spp2, covname, setup, spppair) {
    #'  Reformat data for ggplot
    #'  Snag marginal detection for each species, convert from wide to long format, 
    #'  and assign species name
    marg.det <- as.data.frame(predicted[[1]][,,"mean"]) %>%
      pivot_longer(cols = c(Spp1, Spp2), names_to = "Species") %>%
      arrange(Species) %>%
      mutate(Species = ifelse(Species == "Spp1", spp1, spp2),
             Species = factor(Species, levels = c(spp1, spp2)))
    #'  Snag lower 95% credible interval
    lower.marg <- as.data.frame(predicted[[1]][,,"lower"]) %>%
      pivot_longer(cols = c(Spp1, Spp2), names_to = "Species") %>%
      arrange(Species)
    #'  Snag upper 95% credible interval 
    upper.marg <- as.data.frame(predicted[[1]][,,"upper"]) %>%
      pivot_longer(cols = c(Spp1, Spp2), names_to = "Species") %>%
      arrange(Species)
    #'  Snag covariate values
    covs <- c(predicted[[3]], predicted[[3]]) 
    scaled.covs <- c(predicted[[4]], predicted[[4]])
    #'  Create single data frame
    predicted.marginal.det <- cbind(marg.det, lower.marg[,2], upper.marg[,2], covs, scaled.covs)
    names(predicted.marginal.det) <- c("Species", "marginal_det", "lowerCRI", "upperCRI", "covs", "scaled_covs")
    
    #'  Plot species-specific marginal detection probabilities & 95% CRI
    marg_det_plot <- ggplot(predicted.marginal.det, aes(x = covs, y = marginal_det, group = Species)) + 
      geom_line(aes(color = Species), lwd = 1.25) + 
      scale_color_manual(values = two_colors) + 
      #'  Add confidence intervals
      geom_ribbon(aes(ymin = lowerCRI, ymax = upperCRI, fill = Species), alpha = 0.3) +
      scale_fill_manual(values = two_colors) + 
      #'  Get rid of lines and gray background
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      #'  Force y-axis from 0 to 1
      ylim(0,1.0) +
      #'  Use list name as X-axis title
      xlab(covname) +
      ylab(paste("Marginal Detection Probability,", setup)) +
      labs(#title = paste(spppair, "Marginal Detection Probabilities"), 
        fill = "Species", color = "Species") +
      facet_wrap(~Species, scales = "free_y") +
      theme(legend.position="bottom")
    #'  Review figure
    plot(marg_det_plot)
    
    return(marg_det_plot)
  }
  #####  Wolf-Bear marginal detection  ####
  wolf.bear.marg.eff.ung <- plot_marginal_det(predicted = wolf.bear.eff.ung.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.eff.pred <- plot_marginal_det(predicted = wolf.bear.eff.pred.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.plots <- list(wolf.bear.marg.eff.ung, wolf.bear.marg.eff.pred)
  
  #####  Wolf-Lion marginal detection  ####
  #'  Nadda
  
  #####  Wolf-Coyote marginal detection  ####
  wolf.coy.marg.eff.ung <- plot_marginal_det(predicted = wolf.coy.eff.ung.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Sampling Effort", setup = "random sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.eff.pred <- plot_marginal_det(predicted = wolf.coy.eff.pred.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Sampling Effort", setup = "trail sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.plots <- list(wolf.coy.marg.eff.ung, wolf.coy.marg.eff.pred)
  
  #####  Lion-Bear marginal detection  ####
  #'  Nadda
  
  #####  Lion-Bobcat marginal detection  ####
  lion.bob.marg.ung <- plot_marginal_det(predicted = lion.bob.ung.yr1, spp1 = "Mountain Lion", spp2 = "Bobcat", covname = "Sampling Effort", setup = "random sites", spppair = "Mountain Lion - Bobcat")
  lion.bob.marg.pred <- plot_marginal_det(predicted = lion.bob.pred.yr1, spp1 = "Mountain Lion", spp2 = "Bobcat", covname = "Sampling Effort", setup = "trail sites", spppair = "Mountain Lion - Bobcat")
  lion.bob.marg.plots <- list(lion.bob.marg.ung, lion.bob.marg.pred)
  
  #####  Coyote-Bobcat marginal detection  ####
  coy.bob.marg.eff.ung <- plot_marginal_det(predicted = coy.bob.eff.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Sampling Effort", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.eff.pred <- plot_marginal_det(predicted = coy.bob.eff.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Sampling Effort", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.plots <- list(coy.bob.marg.eff.ung, coy.bob.marg.eff.pred)
  
  #####  Black bear-Coyote marginal detection  ####
  bear.coy.marg.eff.ung <- plot_marginal_det(predicted = bear.coy.eff.ung.yr1, spp1 = "Black bear", spp = "Coyote", covname = "Sampling Effort", setup = "random sites", spppair = "Black bear - Coyote")
  bear.coy.marg.eff.pred <- plot_marginal_det(predicted = bear.coy.eff.pred.yr1, spp1 = "Black bear", spp = "Coyote", covname = "Sampling Effort", setup = "trail sites", spppair = "Black bear - Coyote")
  bear.coy.marg.plots <- list(bear.coy.marg.eff.ung, bear.coy.marg.eff.pred)
  
  #'  -----------------------------------
  ####  Plot conditional Pr(detection)  ####
  #'  -----------------------------------
  #'  Color-blind friendly color palette from Khroma
  plot_scheme(colour("sunset")(11))
  colour("sunset")(11)
  # four_colors <- c("#364B9A", "#C2E4EF", "#FEDA8B", "#A50026")
  four_colors <- c("#364B9A", "#98CAE1", "#FDB366", "#A50026")
  
  #'  Function to reformat data for easier use with ggplot2 
  condish.det <- function(predicted, spp1, spp2) {
    condish.det <- as.data.frame(predicted[[2]][,,"mean"]) %>%
      #'  Snag conditional detection for each species, convert from wide to long format, 
      #'  and assign species name
      pivot_longer(cols = c(Spp1.alone, Spp1.given.Spp2, Spp2.alone, Spp2.given.Spp1), 
                   names_to = "Spp_Interaction") %>%
      arrange(Spp_Interaction) %>%
      #'  Make these labels meaningful with actual species names
      mutate(Spp_Interaction = factor(Spp_Interaction, levels = c("Spp1.alone", "Spp1.given.Spp2", "Spp2.alone", "Spp2.given.Spp1")),
             Spp_Interaction = gsub("Spp1", spp1, Spp_Interaction),
             Spp_Interaction = gsub("Spp2", spp2, Spp_Interaction),
             Spp_Interaction = gsub(".alone", " alone", Spp_Interaction),
             Spp_Interaction = gsub(".given.", " given ", Spp_Interaction), 
             Spp_Interaction = str_to_sentence(Spp_Interaction),
             Species = sub(" .*", "", Spp_Interaction),
             Species = ifelse(Species == "Mountain", "Mountain lion", Species)) %>%
      relocate(Species, .before = Spp_Interaction)
    #'  Snag lower 95% credible interval
    lower.condish <- as.data.frame(predicted[[2]][,,"lower"]) %>%
      pivot_longer(cols = c(Spp1.alone, Spp1.given.Spp2, Spp2.alone, Spp2.given.Spp1), names_to = "Spp_Interaction") %>%
      arrange(Spp_Interaction)
    #'  Snag upper 95% credible interval 
    upper.condish <- as.data.frame(predicted[[2]][,,"upper"]) %>%
      pivot_longer(cols = c(Spp1.alone, Spp1.given.Spp2, Spp2.alone, Spp2.given.Spp1), names_to = "Spp_Interaction") %>%
      arrange(Spp_Interaction)
    #'  Snag covariate values
    covs <- c(predicted[[3]], predicted[[3]]) 
    scaled.covs <- c(predicted[[4]], predicted[[4]])
    #'  Create single data frame
    predicted.conditional.det <- cbind(condish.det, lower.condish[,2], upper.condish[,2], covs, scaled.covs)
    names(predicted.conditional.det) <- c("Species", "Species_interaction", "conditional_det", "lowerCRI", "upperCRI", "covs", "scaled_covs")
    
    return(predicted.conditional.det)
  }
  
  #'  Function to plot conditional detection probability across range of values of a covariate
  plot_conditional_det_over_cov <- function(predicted, spp1, spp2, x, covname, setup, spppair) {
    #'  Call condish.det function to format data for ggplot
    predicted.conditional.det <- condish.det(predicted, spp1, spp2)
    
    #'  Plot species-specific conditional detection probabilities & 95% CRI
    condish_det_plot <- ggplot(predicted.conditional.det, aes(x = covs, y = conditional_det, group = Species_interaction)) + 
      geom_line(aes(color = Species_interaction), lwd = 1.25) + 
      scale_color_manual(values = four_colors, labels = c(paste(spp1, "not detected"), paste(spp1, "detected"), paste(spp2, "not detected"), paste(spp2, "detected"))) + 
      #'  Add confidence intervals
      geom_ribbon(aes(ymin = lowerCRI, ymax = upperCRI, fill = Species_interaction), alpha = 0.3) +
      scale_fill_manual(values = four_colors, labels = c(paste(spp1, "not detected"), paste(spp1, "detected"), paste(spp2, "not detected"), paste(spp2, "detected"))) + 
      #'  Get rid of lines and gray background
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      #'  Force y-axis from 0 to 1
      ylim(0,1.0) +
      #'  Use list name as X-axis title
      xlab(x) +
      ylab(paste("Conditional detection probability,", setup)) +
      labs(title = paste("Conditional detection in response to", covname),
           fill = "Species interaction", color = "Species interaction") +
      facet_wrap(~Species, scales = "free_y") +
      theme(legend.position="bottom")
    #'  Review figure
    plot(condish_det_plot)
    
    return(condish_det_plot)
  }
  
  #'  Function to plot mean conditional detection probability 
  plot_conditional_det_mean <- function(predicted, spp1, spp2, x, covname, setup, spppair, spp_order) {
    #'  Call condish.det function to format data for ggplot
    predicted.conditional.det <- condish.det(predicted, spp1, spp2)
    #'  Filter to a single observation per species and species interaction
    predicted.conditional.det <- predicted.conditional.det %>%
      group_by(Species, Species_interaction) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(Detection = ifelse(grepl("given", Species_interaction), "detected", "not detected"),
             Detection = factor(Detection, levels = c("detected", "not detected")),
             Species_detected = ifelse(Species == spp2, paste(spp1, Detection), paste(spp2, Detection)),
             Species_detected = str_to_sentence(Species_detected),
             Species_detected = factor(Species_detected, levels = spp_order))
    
    #'  Plot species-specific conditional detection probabilities & 95% CRI
    mean_condish_det_plot <- ggplot(predicted.conditional.det, aes(x = Detection, y = conditional_det, group = Species_detected)) + 
      geom_errorbar(aes(ymin = lowerCRI, ymax = upperCRI, color = Species_detected), width = 0, position = position_dodge(width = 0.4)) +
      scale_color_manual(values = four_colors) + 
      geom_point(stat = 'identity', aes(col = Species_detected), size = 2.5, position = position_dodge(width = 0.4)) +   
      #'  Get rid of lines and gray background
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      #'  Force y-axis from 0 to 1
      ylim(0,1.0) +
      #'  Use list name as X-axis title
      xlab(x) +
      ylab(paste("Conditional detection probability,", setup)) +
      labs(title = paste("Conditional detection in response to", covname),
           fill = "Species interaction", color = "Species interaction") +
      facet_wrap(~Species, scales = "free_y") +
      theme(legend.position="bottom")
    #'  Review figure
    plot(mean_condish_det_plot)
    
    return(mean_condish_det_plot)
  }
  #####  Wolf-Bear conditional detection  ####  
  #'  Nadda
  
  #####  Wolf-Lion conditional detection  ####
  #'  Nadda
  
  #####  Wolf-Coyote conditional detection  ####
  wolf.coy.condish.eff.ung <- plot_conditional_det_over_cov(predicted = wolf.coy.eff.ung.yr1, spp1 = "Wolf", spp2 = "Coyote", x = "Sampling Effort", covname = "sampling effort", setup = "random sites", spppair = "Wolf - Coyote")
  wolf.coy.condish.eff.pred <- plot_conditional_det_over_cov(predicted = wolf.coy.eff.pred.yr1, spp1 = "Wolf", spp2 = "Coyote", x = "Sampling Effort", covname = "sampling effort", setup = "trail sites", spppair = "Wolf - Coyote")
  wolf.coy.condish.ung <- plot_conditional_det_mean(predicted = wolf.coy.ung.yr1, spp1 = "Wolf", spp2 = "Coyote", x = "Competitor detection", covname = "competitor detection", setup = "random sites", spppair = "Wolf - Coyote", spp_order = c("Wolf not detected", "Wolf detected", "Coyote not detected", "Coyote detected"))
  wolf.coy.condish.pred <- plot_conditional_det_mean(predicted = wolf.coy.pred.yr1, spp1 = "Wolf", spp2 = "Coyote", x = "Competitor detection", covname = "competitor detection", setup = "trail sites", spppair = "Wolf - Coyote", spp_order = c("Wolf not detected", "Wolf detected", "Coyote not detected", "Coyote detected"))
  wolf.coy.condish.plots <- list(wolf.coy.condish.eff.ung, wolf.coy.condish.eff.pred, wolf.coy.condish.ung, wolf.coy.condish.pred)
  
  #####  Lion-Bear conditional detection  ####
  #'  Nadda
  
  #####  Lion-Bobcat conditional detection  ####
  lion.bob.condish.ung <- plot_conditional_det_mean(predicted = lion.bob.ung.yr1, spp1 = "Mountain lion", spp2 = "Bobcat", x = "Competitor detection", covname = "competitor detection", setup = "random sites", spppair = "Mountain lion - Bobcat", spp_order = c("Mountain lion not detected", "Mountain lion detected", "Bobcat not detected", "Bobcat detected"))
  lion.bob.condish.pred <- plot_conditional_det_mean(predicted = lion.bob.pred.yr1, spp1 = "Mountain lion", spp2 = "Bobcat", x = "Competitor detection", covname = "competitor detection", setup = "trail sites", spppair = "Mountain lion - Bobcat", spp_order = c("Mountain lion not detected", "Mountain lion detected", "Bobcat not detected", "Bobcat detected"))
  lion.bob.condish.plots <- list(lion.bob.condish.ung, lion.bob.condish.pred)
  
  #####  Coyote-Bobcat conditional detection  ####
  coy.bob.condish.eff.ung <- plot_conditional_det_over_cov(predicted = coy.bob.eff.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Sampling Effort", covname = "sampling effort", setup = "random sites")
  coy.bob.condish.eff.pred <- plot_conditional_det_over_cov(predicted = coy.bob.eff.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Sampling Effort", covname = "sampling effort", setup = "trail sites")
  coy.bob.condish.ung <- plot_conditional_det_mean(predicted = coy.bob.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Competitor detection", covname = "competitor detection", setup = "random sites", spppair = "Coyote - Bobcat", spp_order = c("Coyote not detected", "Coyote detected", "Bobcat not detected", "Bobcat detected"))
  coy.bob.condish.pred <- plot_conditional_det_mean(predicted = coy.bob.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Competitor detection", covname = "competitor detection", setup = "trail sites", spppair = "Coyote - Bobcat", spp_order = c("Coyote not detected", "Coyote detected", "Bobcat not detected", "Bobcat detected"))
  coy.bob.condish.plots <- list(coy.bob.condish.eff.ung, coy.bob.condish.eff.pred, coy.bob.condish.ung, coy.bob.condish.pred)
  
  #####  Black bear-Coyote conditional detection  ####
  bear.coy.condish.eff.ung <- plot_conditional_det_over_cov(predicted = bear.coy.eff.ung.yr1, spp1 = "Black bear", spp2 = "Coyote", x = "Sampling Effort", covname = "sampling effort", setup = "random sites")
  bear.coy.condish.eff.pred <- plot_conditional_det_over_cov(predicted = bear.coy.eff.pred.yr1, spp1 = "Black bear", spp2 = "Coyote", x = "Sampling Effort", covname = "sampling effort", setup = "trail sites")
  bear.coy.condish.ung <- plot_conditional_det_mean(predicted = bear.coy.ung.yr1, spp1 = "Black bear", spp2 = "Coyote", x = "Competitor detection", covname = "competitor detection", setup = "random sites", spppair = "Black bear - Coyote", spp_order = c("Black bear not detected", "Black bear detected", "Coyote not detected", "Coyote detected"))
  bear.coy.condish.pred <- plot_conditional_det_mean(predicted = bear.coy.pred.yr1, spp1 = "Black bear", spp2 = "Coyote", x = "Competitor detection", covname = "competitor detection", setup = "trail sites", spppair = "Black bear - Coyote", spp_order = c("Black bear not detected", "Black bear detected", "Coyote not detected", "Coyote detected"))
  bear.coy.condish.plots <- list(bear.coy.condish.eff.ung, bear.coy.condish.eff.pred, bear.coy.condish.ung, bear.coy.condish.pred)
  
  
  #'  Add species pair and camera placement to each data set
  add_info <- function(cond, spp_pair, cam_setup) {
    cond$Species_pair <- spp_pair
    cond$Setup <- cam_setup
    return(cond)
  }
  wolf.coy.condish.ung_new <- add_info(wolf.coy.condish.ung$data, spp_pair = "Wolf - Coyote", cam_setup = "Random sites")
  wolf.coy.condish.pred_new <- add_info(wolf.coy.condish.pred$data, spp_pair = "Wolf - Coyote", cam_setup = "Trail sites")
  lion.bob.condish.ung_new <- add_info(lion.bob.condish.ung$data, spp_pair = "Mountain lion - Bobcat", cam_setup = "Random sites")
  lion.bob.condish.pred_new <- add_info(lion.bob.condish.pred$data, spp_pair = "Mountain lion - Bobcat", cam_setup = "Trail sites")
  coy.bob.condish.ung_new <- add_info(coy.bob.condish.ung$data, spp_pair = "Coyote - Bobcat", cam_setup = "Random sites")
  coy.bob.condish.pred_new <- add_info(coy.bob.condish.pred$data, spp_pair = "Coyote - Bobcat", cam_setup = "Trail sites")
  
  conditional_det <- rbind(wolf.coy.condish.pred_new, lion.bob.condish.pred_new, coy.bob.condish.pred_new) 
  bob.coy.lion.wolf_colors <- c("#A50026", "#98CAE1", "#FDB366", "#364B9A")
  
  #'  Plot all conditional detection probabilities covariate together
  plot_all_condish_det <- function(predicted, x, ncolor) {
    predicted <- mutate(predicted, Species = factor(Species, levels = c("Bobcat", "Coyote", "Mountain lion", "Wolf")))
    cond_det_plot <- ggplot(predicted, aes(x = Detection, y = conditional_det, group = Species)) + 
      geom_errorbar(aes(ymin = lowerCRI, ymax = upperCRI, color = Species), width = 0, position = position_dodge(width = 0.4)) +
      scale_color_manual(values = ncolor) + 
      geom_point(stat = 'identity', aes(col = Species), size = 2.5, position = position_dodge(width = 0.4)) +   
      #'  Get rid of lines and gray background
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      #'  Force y-axis from 0 to 1
      ylim(0,1.0) +
      #'  Use list name as X-axis title
      xlab(x) +
      ylab("Conditional detection probability") +
      labs(title = "Co-detection probabilities for predator dyads", 
           fill = "Focal species", color = "Focal species") +
      facet_wrap(~Species_pair, scales = "free_y") +
      coord_cartesian(ylim = c(0, 0.70)) +
      theme(legend.position="bottom")
    #'  Review figure
    plot(cond_det_plot)
    
    return(cond_det_plot)
  }
  condish_plot <- plot_all_condish_det(predicted = conditional_det, x = "Whether competitor was detected at same site", 
                                       ncolor = bob.coy.lion.wolf_colors)
  
  #'  Plot each species separately while keeping pairings together
  cond_coybob_plot <- ggplot(conditional_det[conditional_det$Species_pair == "Coyote - Bobcat",], aes(x = Detection, y = conditional_det, group = Species)) + 
    geom_errorbar(aes(ymin = lowerCRI, ymax = upperCRI, color = Species), width = 0, position = position_dodge(width = 0.4)) +
    scale_color_manual(values = bob.coy.lion.wolf_colors) + 
    geom_point(stat = 'identity', aes(col = Species), size = 2.5, position = position_dodge(width = 0.4)) +   
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    #'  Force y-axis from 0 to 1
    ylim(0,1.0) +
    #'  Use list name as X-axis title
    xlab("") +
    ylab("Conditional detection probability (trail sites)") +
    facet_wrap(~Species, ncol = 1, scales = "free_y") +
    coord_cartesian(ylim = c(0, 0.70)) +
    #theme(strip.text = element_text(face = "bold", color = "black", hjust = 0, size = 12)) +
    theme(legend.position="none")
  cond_coybob_plot
  
  cond_apexmeso_plot <- ggplot(conditional_det[conditional_det$Species_pair != "Coyote - Bobcat",], aes(x = Detection, y = conditional_det, group = Species)) + 
    geom_errorbar(aes(ymin = lowerCRI, ymax = upperCRI, color = Species), width = 0, position = position_dodge(width = 0.4)) +
    scale_color_manual(values = bob.coy.lion.wolf_colors) + 
    geom_point(stat = 'identity', aes(col = Species), size = 2.5, position = position_dodge(width = 0.4)) +   
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    #'  Force y-axis from 0 to 1
    ylim(0,1.0) +
    #'  Use list name as X-axis title
    xlab("") +
    labs(fill = "Focal species in predator dyad", color = "Focal species in predator dyad") +
    facet_wrap(~Species, ncol = 2, scales = "free_y") +
    coord_cartesian(ylim = c(0, 0.70)) +
    #theme(strip.text = element_text(face = "bold", color = "black", hjust = 0, size = 12)) +
    theme(legend.position = "bottom") 
  cond_apexmeso_plot
  
  #'  remove legend from coy-bob plot
  coybob_guide <- cond_coybob_plot + guides(colour = "none")
  
  #'  Merge coybob and apexmeso plots into single figure
  condish_det_patchwork <- coybob_guide + cond_apexmeso_plot +
    plot_layout(widths = c(1,2)) + plot_annotation(title = 'Co-detection probabilities for predator dyads') +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  #'  Add a single unifying xaxis title to bottom of plot
  condish_det_patchwork <- wrap_elements(panel = condish_det_patchwork) +
    labs(tag = "Whether competitor was detected at same site") +
    theme(
      plot.tag = element_text(size = rel(1)),
      plot.tag.position = c(0.5, 0.10) #"bottom" # note: providing coordinates for tag position makes things a bit complicated when saving
    )
  condish_det_patchwork
  
  #'  -----------------------------
  ####  Save all the pretty plots  ####
  #'  -----------------------------
  #'  Marginal detection probabilities
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_eff_ung_marginal_det_plots.tiff", wolf.bear.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_eff_pred_marginal_det_plots.tiff", wolf.bear.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_eff_ung_marginal_det_plots.tiff", wolf.coy.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_eff_pred_marginal_det_plots.tiff", wolf.coy.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/lion-bob_eff_ung_marginal_det_plots.tiff", lion.bob.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/lion-bob_eff_pred_marginal_det_plots.tiff", lion.bob.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_eff_ung_marginal_det_plots.tiff", coy.bob.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_eff_pred_marginal_det_plots.tiff", coy.bob.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Conditional detection probabilities
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_eff_ung_conditional_det_plots.tiff", coy.bob.condish.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_eff_pred_conditional_det_plots.tiff", coy.bob.condish.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_mean_ung_conditional_det_plots.tiff", coy.bob.condish.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_mean_pred_conditional_det_plots.tiff", coy.bob.condish.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_eff_ung_conditional_det_plots.tiff", wolf.coy.condish.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_eff_pred_conditional_det_plots.tiff", wolf.coy.condish.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_mean_ung_conditional_det_plots.tiff", wolf.coy.condish.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_mean_pred_conditional_det_plots.tiff", wolf.coy.condish.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/lion-bob_mean_ung_conditional_det_plots.tiff", lion.bob.condish.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/lion-bob_mean_pred_conditional_det_plots.tiff", lion.bob.condish.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/all_signif_pairs_mean_pred_conditional_det_plots.tiff", condish_plot, 
         units = "in", width = 8, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/all_signif_pairs_mean_pred_conditional_det_plots_v2.tiff", condish_det_patchwork, 
         units = "in", width = 8, height = 7, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  
  #'  ----------------------------------------
  ####  Predict mean detection for Yr1 & Yr2  ####
  #'  ----------------------------------------
  #'  Load top models
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2024-07-17.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2024-08-29.RData")  
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2024-07-21.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_2024-08-29.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_2024-08-29.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData") 
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData") 
  
  #'  Add binary Year variable to covariate data frame
  stations_skinny_eoe20s21s <- mutate(stations_skinny_eoe20s21s, 
                                      Year = ifelse(Season == "Smr20", 0, 1))
  
  #'  Predict mean detection for Yr1 and Yr2
  wolf.bear.mean.yr1 <- predict_detection(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                          rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  wolf.bear.mean.yr2 <- predict_detection(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          rho_cov = c(1, 0, 1), rho_cov_index = 0,
                                          rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  wolf.coy.mean.yr1 <- predict_detection(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                         rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  wolf.coy.mean.yr2 <- predict_detection(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         rho_cov = c(1, 0, 1), rho_cov_index = 0,
                                         rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  wolf.lion.mean.yr1 <- predict_detection(mod = wolf.lion.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                          rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  wolf.lion.mean.yr2 <- predict_detection(mod = wolf.lion.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          rho_cov = c(1), rho_cov_index = 0,
                                          rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  lion.bear.mean.yr1 <- predict_detection(mod = lion.bear.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          rho_cov = c(1), rho_cov_index = 0,
                                          rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  lion.bear.mean.yr2 <- predict_detection(mod = lion.bear.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          rho_cov = c(1), rho_cov_index = 0,
                                          rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  lion.bob.mean.yr1 <- predict_detection(mod = lion.bob.null, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         rho_cov = c(1), rho_cov_index = 0,
                                         rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  lion.bob.mean.yr2 <- predict_detection(mod = lion.bob.null, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         rho_cov = c(1), rho_cov_index = 0,
                                         rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  coy.bob.mean.yr1 <- predict_detection(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                        focal_cov = stations_skinny_eoe20s21s$Year,
                                        rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                        rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  coy.bob.mean.yr2 <- predict_detection(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                        focal_cov = stations_skinny_eoe20s21s$Year,
                                        rho_cov = c(1, 0, 1), rho_cov_index = 0,
                                        rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  bear.coy.mean.yr1 <- predict_detection(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         rho_cov = c(1, 0, 0), rho_cov_index = 0,
                                         rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  bear.coy.mean.yr2 <- predict_detection(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         rho_cov = c(1, 0, 1), rho_cov_index = 0,
                                         rho_inxs_cov = c(0), rho_inxs_cov_index = 0)
  
  #'  Extract mean and 95% CRI for each observation (should be identical) and thin
  #'  to one value per species
  mean_det_prob <- function(predicted, spp1, spp2, yr) {
    #'  Snag mean
    marg.det <- as.data.frame(predicted[[1]][,,"mean"]) %>%
      pivot_longer(cols = c(Spp1, Spp2), names_to = "Species") %>%
      arrange(Species) %>%
      mutate(Species = ifelse(Species == "Spp1", spp1, spp2),
             Species = factor(Species, levels = c(spp1, spp2)))
    #'  Snag lower 95% credible interval
    lower.marg <- as.data.frame(predicted[[1]][,,"lower"]) %>%
      pivot_longer(cols = c(Spp1, Spp2), names_to = "Species") %>%
      arrange(Species)
    #'  Snag upper 95% credible interval 
    upper.marg <- as.data.frame(predicted[[1]][,,"upper"]) %>%
      pivot_longer(cols = c(Spp1, Spp2), names_to = "Species") %>%
      arrange(Species)
    #'  Create single data frame
    predicted.marginal.det <- cbind(marg.det, lower.marg[,2], upper.marg[,2])
    names(predicted.marginal.det) <- c("Species", "marginal_occ", "lowerCRI", "upperCRI")
    #'  Drop extra observations and round
    predicted.marginal.det <- as.data.frame(predicted.marginal.det) %>%
      group_by(Species) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(marginal_det = round(marginal_det, 2),
             lowerCRI = round(lowerCRI, 2),
             upperCRI = round(upperCRI, 2),
             Year = yr)
    return(predicted.marginal.det)
  }
  pmean_wolf.bear.yr1 <- mean_det_prob(wolf.bear.mean.yr1, spp1 = "Wolf", spp2 = "Black bear", yr = "2020")
  pmean_wolf.bear.yr2 <- mean_det_prob(wolf.bear.mean.yr2, spp1 = "Wolf", spp2 = "Black bear", yr = "2021")
  pmean_wolf.coy.yr1 <- mean_det_prob(wolf.coy.mean.yr1, spp1 = "Wolf", spp2 = "Coyote", yr = "2020")
  pmean_wolf.coy.yr2 <- mean_det_prob(wolf.coy.mean.yr2, spp1 = "Wolf", spp2 = "Coyote", yr = "2021")
  pmean_wolf.lion.yr1 <- mean_det_prob(wolf.lion.mean.yr1, spp1 = "Wolf", spp2 = "Mountain lion", yr = "2020")
  pmean_wolf.lion.yr2 <- mean_det_prob(wolf.lion.mean.yr2, spp1 = "Wolf", spp2 = "Mountain lion", yr = "2021")
  pmean_lion.bear.yr1 <- mean_det_prob(lion.bear.mean.yr1, spp1 = "Mountain lion", spp2 = "Black bear", yr = "2020")
  pmean_lion.bear.yr2 <- mean_det_prob(lion.bear.mean.yr2, spp1 = "Mountain lion", spp2 = "Black bear", yr = "2021")
  pmean_lion.bob.yr1 <- mean_det_prob(lion.bob.mean.yr1, spp1 = "Mountain lion", spp2 = "Bobcat", yr = "2020")
  pmean_lion.bob.yr2 <- mean_det_prob(lion.bob.mean.yr2, spp1 = "Mountain lion", spp2 = "Bobcat", yr = "2021")
  pmean_coy.bob.yr1 <- mean_det_prob(coy.bob.mean.yr1, spp1 = "Coyote", spp2 = "Bobcat", yr = "2020")
  pmean_coy.bob.yr2 <- mean_det_prob(coy.bob.mean.yr2, spp1 = "Coyote", spp2 = "Bobcat", yr = "2021")
  pmean_bear.coy.yr1 <- mean_det_prob(bear.coy.mean.yr1, spp1 = "Black bear", spp2 = "Coyote", yr = "2020")
  pmean_bear.coy.yr2 <- mean_det_prob(bear.coy.mean.yr2, spp1 = "Black bear", spp2 = "Coyote", yr = "2021")
  
  #'  Bind all species and years together
  mean_annual_p <- rbind(pmean_wolf.bear.yr1, pmean_wolf.bear.yr2, pmean_wolf.lion.yr1, 
                           pmean_wolf.lion.yr2, pmean_wolf.coy.yr1, pmean_wolf.coy.yr2, 
                           pmean_lion.bear.yr1, pmean_lion.bear.yr2, pmean_lion.bob.yr1,
                           pmean_lion.bob.yr2, pmean_coy.bob.yr1, pmean_coy.bob.yr2,
                           pmean_bear.coy.yr1, pmean_bear.coy.yr2) %>%
    #'  Reduce to one observation per species
    group_by(Species, Year) %>%
    slice(1L) %>%
    ungroup() %>%
    #'  Combine mean & CRI into single column
    mutate(CRI = paste0("(", lowerCRI, ", ", upperCRI, ")"),
           Mean = paste(marginal_det, CRI)) %>%
    dplyr::select(Species, Year, Mean) %>%
    spread(Year, Mean) %>%
    rename("Mean_CRI_2020" = "2020", "Mean_CRI_2021" = "2021")
  
  #'  Save!
  write.csv(mean_annual_p, "./Outputs/Tables/Summary_mean_p_yr1v2.csv")
  
  
