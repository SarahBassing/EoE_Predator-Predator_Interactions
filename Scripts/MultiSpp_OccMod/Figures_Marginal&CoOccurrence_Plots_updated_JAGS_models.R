  #'  ------------------------------
  #'  Plot co-occurrence results
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  -------------------------------
  #'  Script to create figures of predicted co-occurrence and covariate effects
  #'  from multispecies occupancy models run in JAGS. Also predict mean probabilties
  #'  of occupancy and detection for each species based on null model.
  #'  -------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(ggplot2)
  library(khroma)
  library(stringr)
  library(patchwork)
  # library(ggspatial)
  library(grid)
  library(png)
  library(RCurl)
  library(rphylopic)
  library(tidyverse)
  
  #'  Load covariate data
  load("./Data/Covariates_extracted/Covariate_skinny_EoE20s21s_updated_070824.RData")
  
  #'  Identify top models
  load("./Outputs/Tables/DIC_top_models.RData")
  print(topmodels)
  
  #'  -------------------------------------------------
  ####  Predict Pr(occupancy) across covariate values  ####
  #'  -------------------------------------------------
  #'  MASSIVE function the predict marginal and conditional occupancy for each species
  predict_occupancy <- function(mod, ncat, npoints, focal_cov, psi_cov, psi_inxs_cov, psi_cov_index, psi_inxs_cov_index) {
    #'  Rename model and review
    out1 <- mod
    # print(out1$summary[1:50,])
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
    
    #'  Create model matrices for prediction covariates on psi and psix
    #'  Fill matrix with categorical or mean value of each covariate
    psi_cov <- psi_cov
    psi_covs <- matrix(psi_cov, nrow = npoints, ncol = length(psi_cov), byrow = TRUE)
    #'  Same thing but for covariates on interaction
    psi_inxs_cov <- psi_inxs_cov
    psi_inxs_covs <- matrix(psi_inxs_cov, nrow = npoints, ncol = length(psi_inxs_cov), byrow = TRUE)
    #'  Replace focal covariate column with new scaled data
    #'  New data must be in correct column based on order of covariates in original model
    psi_covs[,psi_cov_index] <- cov.pred
    print(head(psi_covs))
    psi_inxs_covs[,psi_inxs_cov_index] <- cov.pred
    print(head(psi_inxs_covs))
    
    #'  Assemble linear predictors across all iterations of the MCMC chains
    psiSpp1 <- out1$sims.list$betaSpp1 %*% t(psi_covs)
    psiSpp2 <- out1$sims.list$betaSpp2 %*% t(psi_covs)
    psiSpp12 <- psiSpp1 + psiSpp2 + out1$sims.list$betaSpp12 %*% t(psi_inxs_covs) 
    
    #'  Compute latent state vector (lsv) based on new data and means of other covariates
    lsv <- array(NA, dim = c(ndraws, npoints, ncat))
    lsv[,,1] <- 1
    lsv[,,2] <- exp(psiSpp1)
    lsv[,,3] <- exp(psiSpp2)
    lsv[,,4] <- exp(psiSpp12)
    
    #'  Compute latent state probability of each state from MCMC chains (lsp)
    #'  denominator  = sum of lsv for each draw x point
    denom <- apply(lsv, 1:2, sum)
    lspMC <- sweep(lsv, 1:2, denom, "/")
    dimnames(lspMC)[[3]] <- c("U", "Spp1", "Spp2", "Both")
    str(lspMC)
    
    #'  Fill in summary values for 'lsp'
    lsp <- array(NA, dim = c(npoints, ncat, 4))
    dimnames(lsp) <- list(NULL, c("U", "Spp1", "Spp2", "Both"),
                          c("mean", "SD", "lower", "upper"))
    lsp[,,1] <- apply(lspMC, 2:3, mean)
    lsp[,,2] <- apply(lspMC, 2:3, sd)
    lsp[,,3] <- apply(lspMC, 2:3, quantile, probs = 0.025)
    lsp[,,4] <- apply(lspMC, 2:3, quantile, probs = 0.975)
    
    #'  Compute and plot marginal probabilities
    #'  --------------------------------------- 
    #'  Occupancy probability for each species, regardless of presence or
    #'  absence of other species
    #'  Get MCMC chains
    marginalMC <- array(NA, dim = c(ndraws, npoints, 2))
    dimnames(marginalMC)[[3]] <- c("Spp1", "Spp2")
    marginalMC[,,"Spp1"] <- apply(lspMC[,,c("Spp1", "Both")], 1:2, sum)
    marginalMC[,,"Spp2"] <- apply(lspMC[,,c("Spp2", "Both")], 1:2, sum)
    
    #'  Fill in summary values for marginal probabilities
    marginal <- array(NA, dim = c(npoints, 2, 4))
    dimnames(marginal) <- list(NULL, c("Spp1", "Spp2"),
                               c("mean", "sd", "lower", "upper"))
    marginal[,,1] <- apply(marginalMC, 2:3, mean)
    marginal[,,2] <- apply(marginalMC, 2:3, sd)
    marginal[,,3] <- apply(marginalMC, 2:3, quantile, probs = 0.025)
    marginal[,,4] <- apply(marginalMC, 2:3, quantile, probs = 0.975)
    
    #'  Compute conditional occupancy for all species
    #'  --------------------------------------------- 
    #'  Occupancy probability for each species, conditional on the presence or
    #'  absence of the other species
    #'  Get MCMC chains
    conditionalMC <- array(NA, dim = c(ndraws, npoints, 4))
    dimnames(conditionalMC)[[3]] <- c("Spp1.alone", "Spp1.given.Spp2",
                                      "Spp2.alone", "Spp2.given.Spp1")
    
    conditionalMC[,,"Spp1.alone"] <- lspMC[,,"Spp1"]/apply(lspMC[,,c("U", "Spp1")], 1:2, sum)
    conditionalMC[,,"Spp1.given.Spp2"] <- lspMC[,,"Both"]/marginalMC[,,"Spp2"]
    conditionalMC[,,"Spp2.alone"] <- lspMC[,,"Spp2"]/apply(lspMC[,,c("U", "Spp2")], 1:2, sum)
    conditionalMC[,,"Spp2.given.Spp1"] <- lspMC[,,"Both"]/marginalMC[,,"Spp1"]
    
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
  #'  Top model
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_2024-08-27.RData") 
  
  #'  Second most supported model (1.17 deltaDIC away from top)
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2024-07-17.RData")
  wolf.bear.for.ung.yr2 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$PercForest,
                                             psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 4,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.for.pred.yr2 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$PercForest,
                                              psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 4,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.elev.ung.yr2 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$Elev,
                                              psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 5,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.elev.pred.yr2 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                               focal_cov = stations_skinny_eoe20s21s$Elev,
                                               psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 5,
                                               psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.tri.ung.yr2 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$TRI,
                                             psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 6,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.tri.pred.yr2 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$TRI,
                                              psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 6,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  
  #####  Wolf-Lion predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2024-08-29.RData")  
  
  
  #####  Wolf-Coyote predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2024-07-21.RData") 
  
  wolf.coy.for.ung.yr2 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                            focal_cov = stations_skinny_eoe20s21s$PercForest,
                                            psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 4,
                                            psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.for.pred.yr2 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$PercForest,
                                             psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 4,
                                             psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.elev.ung.yr2 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 5,
                                             psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.elev.pred.yr2 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$Elev,
                                              psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 5,
                                              psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.tri.ung.yr2 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$TRI,
                                             psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 6,
                                             psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.tri.pred.yr2 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$TRI,
                                              psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 6,
                                              psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  
  #####  Lion-Bear predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2024-08-29.RData")  
  
  
  #####  Lion-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2024-08-29.RData")  
  
  
  #####  Coyote-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData")
  coy.bob.for.ung.yr2 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                           focal_cov = stations_skinny_eoe20s21s$PercForest,
                                           psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 4,
                                           psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  coy.bob.for.pred.yr2 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                            focal_cov = stations_skinny_eoe20s21s$PercForest,
                                            psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 4,
                                            psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  coy.bob.elev.ung.yr2 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500, 
                                            focal_cov = stations_skinny_eoe20s21s$Elev,
                                            psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 5,
                                            psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  coy.bob.elev.pred.yr2 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 5,
                                             psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  coy.bob.tri.ung.yr2 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500, 
                                           focal_cov = stations_skinny_eoe20s21s$TRI,
                                           psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 6,
                                           psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  coy.bob.tri.pred.yr2 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                            focal_cov = stations_skinny_eoe20s21s$TRI,
                                            psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 6,
                                            psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  
  
  #####  Bear-Coyote predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/bearcoy_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2024-07-24.RData")
  bear.coy.for.ung.yr2 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                           focal_cov = stations_skinny_eoe20s21s$PercForest,
                                           psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 4,
                                           psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  bear.coy.for.pred.yr2 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                            focal_cov = stations_skinny_eoe20s21s$PercForest,
                                            psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 4,
                                            psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  bear.coy.elev.ung.yr2 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500, 
                                            focal_cov = stations_skinny_eoe20s21s$Elev,
                                            psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 5,
                                            psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  bear.coy.elev.pred.yr2 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 5,
                                             psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  bear.coy.tri.ung.yr2 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500, 
                                           focal_cov = stations_skinny_eoe20s21s$TRI,
                                           psi_cov = c(1, 0, 1, 0, 0, 0), psi_cov_index = 6,
                                           psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  bear.coy.tri.pred.yr2 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                            focal_cov = stations_skinny_eoe20s21s$TRI,
                                            psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 6,
                                            psi_inxs_cov = c(1), psi_inxs_cov_index = 0)

  
  save.image(file = paste0("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/Predicted_psi-cov_relationships_", Sys.Date(), ".RData"))
  
  
  #'  --------------------------------
  ####  Plot marginal Pr(occupancy)  ####
  #'  --------------------------------
  load("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/Predicted_psi-cov_relationships_2024-08-29.RData")
  
  #'  Color-blind friendly color palette from Khroma
  plot_scheme(colour("sunset")(11))
  colour("sunset")(11)
  # two_colors <- c("#364B9A", "#A50026")
  two_colors <- c("#6EA6CD", "#DD3D2D")
  
  #'  Function to reformat data for easier use with ggplot2 and plot marginal occupancy
  plot_marginal_occ <- function(predicted, spp1, spp2, covname, setup, spppair) {
    #'  Reformat data for ggplot
    #'  Snag marginal occupancy for each species, convert from wide to long format, 
    #'  and assign species name
    marg.occ <- as.data.frame(predicted[[1]][,,"mean"]) %>%
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
    predicted.marginal.occ <- cbind(marg.occ, lower.marg[,2], upper.marg[,2], covs, scaled.covs)
    names(predicted.marginal.occ) <- c("Species", "marginal_occ", "lowerCRI", "upperCRI", "covs", "scaled_covs")
    
    #'  Plot species-specific marginal occupancy probabilities & 95% CRI
    marg_occ_plot <- ggplot(predicted.marginal.occ, aes(x = covs, y = marginal_occ, group = Species)) + 
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
      ylab(paste("Marginal Occupancy Probability,", setup)) +
      labs(#title = paste(spppair, "Marginal Occupancy Probabilities"), 
        fill = "Species", color = "Species") +
      facet_wrap(~Species, scales = "free_y") +
      theme(legend.position="bottom")
    #'  Review figure
    plot(marg_occ_plot)
    
    return(marg_occ_plot)
  }
  #####  Wolf-Bear marginal occupancy  ####
  #'  Second best model
  wolf.bear.marg.for.ung <- plot_marginal_occ(predicted = wolf.bear.for.ung.yr2, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.for.pred <- plot_marginal_occ(predicted = wolf.bear.for.pred.yr2, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "Trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.elev.ung <- plot_marginal_occ(predicted = wolf.bear.elev.ung.yr2, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.elev.pred <- plot_marginal_occ(predicted = wolf.bear.elev.pred.yr2, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.tri.ung <- plot_marginal_occ(predicted = wolf.bear.tri.ung.yr2, spp1 = "Wolf", spp2 = "Bear", covname = "Terrain ruggedness index", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.tri.pred <- plot_marginal_occ(predicted = wolf.bear.tri.pred.yr2, spp1 = "Wolf", spp2 = "Bear", covname = "Terrain ruggedness index", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.plots <- list(wolf.bear.marg.for.ung, wolf.bear.marg.for.pred, wolf.bear.marg.elev.ung, wolf.bear.marg.elev.pred, wolf.bear.marg.tri.ung, wolf.bear.marg.tri.pred)
  
  #####  Wolf-Coyote marginal occupancy  ####
  wolf.coy.marg.for.ung <- plot_marginal_occ(predicted = wolf.coy.for.ung.yr2, spp1 = "Wolf", spp2 = "Coyote", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Coyote") 
  wolf.coy.marg.for.pred <- plot_marginal_occ(predicted = wolf.coy.for.pred.yr2, spp1 = "Wolf", spp2 = "Coyote", covname = "Percent forest cover", setup = "trail sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.elev.ung <- plot_marginal_occ(predicted = wolf.coy.elev.ung.yr2, spp1 = "Wolf", spp2 = "Coyote", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.elev.pred <- plot_marginal_occ(predicted = wolf.coy.elev.pred.yr2, spp1 = "Wolf", spp2 = "Coyote", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.tri.ung <- plot_marginal_occ(predicted = wolf.coy.tri.ung.yr2, spp1 = "Wolf", spp2 = "Coyote", covname = "Terrain ruggedness index", setup = "random sites", spppair = "Wolf - Coyote") 
  wolf.coy.marg.tri.pred <- plot_marginal_occ(predicted = wolf.coy.tri.pred.yr2, spp1 = "Wolf", spp2 = "Coyote", covname = "Terrain ruggedness index", setup = "trail sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.plots <- list(wolf.coy.marg.for.ung, wolf.coy.marg.for.pred, wolf.coy.marg.elev.ung, wolf.coy.marg.elev.pred, wolf.coy.marg.tri.ung, wolf.coy.marg.tri.pred)
  
  #####  Wolf-Lion marginal occupancy  ####
  #'  Nadda
  
  #####  Lion-Bear marginal occupancy  ####
  #'  Nadda
  
  #####  Lion-Bobcat marginal occupancy  ####
  #'  Nadda
  
  #####  Coyote-Bobcat marginal occupancy  ####
  coy.bob.marg.for.ung <- plot_marginal_occ(predicted = coy.bob.for.ung.yr2, spp1 = "Coyote", spp2 = "Bobcat", covname = "Percent forest cover", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.for.pred <- plot_marginal_occ(predicted = coy.bob.for.pred.yr2, spp1 = "Coyote", spp2 = "Bobcat", covname = "Percent forest cover", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.elev.ung <- plot_marginal_occ(predicted = coy.bob.elev.ung.yr2, spp1 = "Coyote", spp2 = "Bobcat", covname = "Elevation (m)", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.elev.pred <- plot_marginal_occ(predicted = coy.bob.elev.pred.yr2, spp1 = "Coyote", spp2 = "Bobcat", covname = "Elevation (m)", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.tri.ung <- plot_marginal_occ(predicted = coy.bob.tri.ung.yr2, spp1 = "Coyote", spp2 = "Bobcat", covname = "Terrain ruggedness index", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.tri.pred <- plot_marginal_occ(predicted = coy.bob.tri.pred.yr2, spp1 = "Coyote", spp2 = "Bobcat", covname = "Terrain ruggedness index", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.plots <- list(coy.bob.marg.for.ung, coy.bob.marg.for.pred, coy.bob.marg.elev.ung, coy.bob.marg.elev.pred, coy.bob.marg.tri.ung, coy.bob.marg.tri.pred)
  
  #####  Bear-Coyote marginal occupancy  ####
  bear.coy.marg.for.ung <- plot_marginal_occ(predicted = bear.coy.for.ung.yr2, spp1 = "Bear", spp2 = "Coyote", covname = "Percent forest cover", setup = "random sites", spppair = "Black Bear - Coyote")
  bear.coy.marg.for.pred <- plot_marginal_occ(predicted = bear.coy.for.pred.yr2, spp1 = "Bear", spp2 = "Coyote", covname = "Percent forest cover", setup = "trail sites", spppair = "Black Bear - Coyote")
  bear.coy.marg.elev.ung <- plot_marginal_occ(predicted = bear.coy.elev.ung.yr2, spp1 = "Bear", spp2 = "Coyote", covname = "Elevation (m)", setup = "random sites", spppair = "Black Bear - Coyote")
  bear.coy.marg.elev.pred <- plot_marginal_occ(predicted = bear.coy.elev.pred.yr2, spp1 = "Bear", spp2 = "Coyote", covname = "Elevation (m)", setup = "trail sites", spppair = "Black Bear - Coyote")
  bear.coy.marg.tri.ung <- plot_marginal_occ(predicted = bear.coy.tri.ung.yr2, spp1 = "Bear", spp2 = "Coyote", covname = "Terrain ruggedness index", setup = "random sites", spppair = "Black Bear - Coyote")
  bear.coy.marg.tri.pred <- plot_marginal_occ(predicted = bear.coy.tri.pred.yr2, spp1 = "Bear", spp2 = "Coyote", covname = "Terrain ruggedness index", setup = "trail sites", spppair = "Black Bear - Coyote")
  bear.coy.marg.plots <- list(bear.coy.marg.for.ung, bear.coy.marg.for.pred, bear.coy.marg.elev.ung, bear.coy.marg.elev.pred, bear.coy.marg.tri.ung, bear.coy.marg.tri.pred) 
  
  #'  Add species pair and camera placement to each data set
  add_info <- function(marg, spp_pair, cam_setup) {
    marg$Species_pair <- spp_pair
    marg$Setup <- cam_setup
    return(marg)
  }
  wolf.bear.marg.for.ung_new <- add_info(wolf.bear.marg.for.ung$data, spp_pair = "wolf-bear", cam_setup = "Random sites")
  wolf.bear.marg.for.pred_new <- add_info(wolf.bear.marg.for.pred$data, spp_pair = "wolf-bear", cam_setup = "Trail sites")
  wolf.bear.marg.elev.ung_new <- add_info(wolf.bear.marg.elev.ung$data, spp_pair = "wolf-bear", cam_setup = "Random sites")
  wolf.bear.marg.elev.pred_new <- add_info(wolf.bear.marg.elev.pred$data, spp_pair = "wolf-bear", cam_setup = "Trail sites")
  wolf.bear.marg.tri.ung_new <- add_info(wolf.bear.marg.tri.ung$data, spp_pair = "wolf-bear", cam_setup = "Random sites")
  wolf.bear.marg.tri.pred_new <- add_info(wolf.bear.marg.tri.pred$data, spp_pair = "wolf-bear", cam_setup = "Trail sites")
  wolf.coy.marg.for.ung_new <- add_info(wolf.coy.marg.for.ung$data, spp_pair = "wolf-coy", cam_setup = "Random sites")
  wolf.coy.marg.for.pred_new <- add_info(wolf.coy.marg.for.pred$data, spp_pair = "wolf-coy", cam_setup = "Trail sites")
  wolf.coy.marg.elev.ung_new <- add_info(wolf.coy.marg.elev.ung$data, spp_pair = "wolf-coy", cam_setup = "Random sites")
  wolf.coy.marg.elev.pred_new <- add_info(wolf.coy.marg.elev.pred$data, spp_pair = "wolf-coy", cam_setup = "Trail sites")
  wolf.coy.marg.tri.ung_new <- add_info(wolf.coy.marg.tri.ung$data, spp_pair = "wolf-coy", cam_setup = "Random sites")
  wolf.coy.marg.tri.pred_new <- add_info(wolf.coy.marg.tri.pred$data, spp_pair = "wolf-coy", cam_setup = "Trail sites")
  coy.bob.marg.for.ung_new <- add_info(coy.bob.marg.for.ung$data, spp_pair = "coy-bob", cam_setup = "Random sites")
  coy.bob.marg.for.pred_new <- add_info(coy.bob.marg.for.pred$data, spp_pair = "coy-bob", cam_setup = "Trail sites")
  coy.bob.marg.elev.ung_new <- add_info(coy.bob.marg.elev.ung$data, spp_pair = "coy-bob", cam_setup = "Random sites")
  coy.bob.marg.elev.pred_new <- add_info(coy.bob.marg.elev.pred$data, spp_pair = "coy-bob", cam_setup = "Trail sites")
  coy.bob.marg.tri.ung_new <- add_info(coy.bob.marg.tri.ung$data, spp_pair = "coy-bob", cam_setup = "Random sites")
  coy.bob.marg.tri.pred_new <- add_info(coy.bob.marg.tri.pred$data, spp_pair = "coy-bob", cam_setup = "Trail sites")
  bear.coy.marg.for.ung_new <- add_info(bear.coy.marg.for.ung$data, spp_pair = "bear-coy", cam_setup = "Random sites")
  bear.coy.marg.for.pred_new <- add_info(bear.coy.marg.for.pred$data, spp_pair = "bear-coy", cam_setup = "Trail sites")
  bear.coy.marg.elev.ung_new <- add_info(bear.coy.marg.elev.ung$data, spp_pair = "bear-coy", cam_setup = "Random sites")
  bear.coy.marg.elev.pred_new <- add_info(bear.coy.marg.elev.pred$data, spp_pair = "bear-coy", cam_setup = "Trail sites")
  bear.coy.marg.tri.ung_new <- add_info(bear.coy.marg.tri.ung$data, spp_pair = "bear-coy", cam_setup = "Random sites")
  bear.coy.marg.tri.pred_new <- add_info(bear.coy.marg.tri.pred$data, spp_pair = "bear-coy", cam_setup = "Trail sites")
  
  
  #'  Combine all marginal probability estimates for each significant variable [and marginally significant]
  #'  wolf: elev, tri, [for]; bear: forest, tri, [elev]
  #'  wolf: elev, tri, [for]; coy: forest, tri
  #'  coy: forest, tri; bob: forest, [tri]
  #'  bear: forest, tri, [elev]; coy: forest, tri
  marginal_for <- rbind(wolf.bear.marg.for.ung_new, wolf.bear.marg.for.pred_new, 
                        wolf.coy.marg.for.ung_new, wolf.coy.marg.for.pred_new, 
                        coy.bob.marg.for.ung_new, coy.bob.marg.for.pred_new,
                        bear.coy.marg.for.ung_new, bear.coy.marg.for.pred_new)  %>%
    mutate(Species = factor(Species, levels = c("Bear", "Bobcat", "Coyote", "Wolf"))) %>%
    #'  Drop wolf-coy and bear-coy results b/c already have wolf & coy results included w/ other pairings
    filter(Species_pair != "wolf-coy") %>%
    filter(Species_pair != "bear-coy")
  marginal_elev <- rbind(wolf.bear.marg.elev.ung_new, wolf.bear.marg.elev.pred_new,
                         wolf.coy.marg.elev.ung_new, wolf.coy.marg.elev.pred_new,
                         bear.coy.marg.elev.ung_new, bear.coy.marg.elev.pred_new) %>%
    mutate(Species = factor(Species, levels = c("Bear", "Bobcat", "Coyote", "Wolf"))) %>%
    #'  Keep only wolf-bear results b/c elev not important for coy and wolf/bear responses already covered by wolf-bear
    filter(Species_pair == "wolf-bear")
  marginal_tri <- rbind(wolf.bear.marg.tri.ung_new, wolf.bear.marg.tri.pred_new, 
                        wolf.coy.marg.tri.ung_new, wolf.coy.marg.tri.pred_new, 
                        coy.bob.marg.tri.ung_new, coy.bob.marg.tri.pred_new,
                        bear.coy.marg.tri.ung_new, bear.coy.marg.tri.pred_new) %>%
    mutate(Species = factor(Species, levels = c("Bear", "Bobcat", "Coyote", "Wolf"))) %>%
    #'  Drop wolf-coy and bear-coy results b/c already have wolf & coy results included w/ other pairings
    filter(Species_pair != "wolf-coy") %>%
    filter(Species_pair != "bear-coy")
  
  #'  Set color combos for each species
  #'  wolf = "#364B9A", bear = "#98CAE1", lion = "#FDB366", coyote = "#DD3D2D", bobcat = "#A50026"
  bear_colors <- "#98CAE1"
  bear.bob.coy_colors <- c("#98CAE1", "#A50026", "#DD3D2D")
  bear.coy.wolf_colors <- c("#98CAE1", "#DD3D2D","#364B9A")
  bear.wolf_colors <- c("#98CAE1", "#364B9A")
  bear.bob.coy.wolf_colors <- c("#98CAE1", "#A50026", "#DD3D2D", "#364B9A")
  bear.bob.coy.lion.wolf_colors <- c("#98CAE1", "#A50026", "#DD3D2D", "#FDB366", "#364B9A")
  
  #'  Plot each species response to specific covariate together
  plot_margingal_occ_by_cov <- function(predicted, x, ylab, plottitle, ncolor) {
    marg_occ_plot <- ggplot(predicted, aes(x = covs, y = marginal_occ, group = Species)) + 
      geom_line(aes(color = Species), lwd = 1.25) + 
      scale_color_manual(values = ncolor) + 
      #'  Add confidence intervals
      geom_ribbon(aes(ymin = lowerCRI, ymax = upperCRI, fill = Species), alpha = 0.3) +
      scale_fill_manual(values = ncolor) + 
      #'  Get rid of lines and gray background
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      #'  Force y-axis from 0 to 1
      ylim(0,1.0) +
      #'  Use list name as X-axis title
      xlab(x) +
      ylab("Marginal occupancy probability") +
      labs(title = plottitle, 
           fill = "Species", color = "Species") +
      #facet_wrap(~Setup, scales = "free_y") +
      theme(legend.position="bottom")
    #'  Review figure
    plot(marg_occ_plot)
    
    return(marg_occ_plot)
  }
  marginal_for_pred_plot <- plot_margingal_occ_by_cov(marginal_for[marginal_for$Setup == "Trail sites",], x = "Percent forest cover", ylab = "", plottitle = "", ncolor = bear.bob.coy.wolf_colors)
  marginal_elev_pred_plot <- plot_margingal_occ_by_cov(marginal_elev[marginal_elev$Setup == "Trail sites",], x = "Elevation (m)", ylab = "Marginal occupancy probability", plottitle = "", ncolor = bear.wolf_colors) 
  marginal_tri_pred_plot <- plot_margingal_occ_by_cov(marginal_tri[marginal_tri$Setup == "Trail sites",], x = "Terrain ruggedness index", ylab = "", plottitle = "", ncolor = bear.bob.coy.wolf_colors) 
  
  marginal_for_ung_plot <- plot_margingal_occ_by_cov(marginal_for[marginal_for$Setup == "Random sites",], x = "Percent forest cover", ylab = "", plottitle = "", ncolor = bear.bob.coy.wolf_colors)
  marginal_elev_ung_plot <- plot_margingal_occ_by_cov(marginal_elev[marginal_elev$Setup == "Random sites",], x = "Elevation (m)", ylab = "Marginal occupancy probability", plottitle = "", ncolor = bear.wolf_colors) 
  marginal_tri_ung_plot <- plot_margingal_occ_by_cov(marginal_tri[marginal_tri$Setup == "Random sites",], x = "Terrain ruggedness index", ylab = "", plottitle = "", ncolor = bear.bob.coy.wolf_colors) 
  
  
  #'  Remove legend from elev and tri plots
  elev_pred_guide <- marginal_elev_pred_plot + theme(legend.title = element_blank()) + guides(colour = "none", fill = "none") 
  elev_ung_guide <- marginal_elev_ung_plot + theme(legend.title = element_blank()) + guides(colour = "none", fill = "none") 
  tri_pred_guide <- marginal_tri_pred_plot + theme(legend.title = element_blank()) + guides(colour = "none", fill = "none")
  tri_ung_guide <- marginal_tri_ung_plot + theme(legend.title = element_blank()) + guides(colour = "none", fill = "none")
  # for_pred_guide <- marginal_for_pred_plot + guides(title = "none")
  # for_ung_guide <- marginal_for_ung_plot + guides(title = "none")
  
  marginal_pred_patchwork <- marginal_for_pred_plot + elev_pred_guide + tri_pred_guide + 
    plot_annotation(title = "Species-specific marginal occupancy across percent forest cover, elevation, and terrain ruggedness gradients") +
    plot_annotation(tag_levels = 'a') + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  marginal_ung_patchwork <- marginal_for_ung_plot + elev_ung_guide + tri_ung_guide + 
    plot_annotation(title = "Species-specific marginal occupancy across percent forest cover, elevation, and terrain ruggedness gradients") +
    plot_annotation(tag_levels = 'a') + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_occ_plot_predcams.tiff", marginal_pred_patchwork, 
         units = "in", width = 11, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_occ_plot_ungcams.tiff", marginal_ung_patchwork, 
         units = "in", width = 11, height = 6, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  -----------------------------------
  ####  Plot conditional Pr(occupancy)  ####
  #'  -----------------------------------
  #'  Color-blind friendly color palette from Khroma
  plot_scheme(colour("sunset")(11))
  colour("sunset")(11)
  # four_colors <- c("#364B9A", "#C2E4EF", "#FEDA8B", "#A50026")
  four_colors <- c("#364B9A", "#98CAE1", "#FDB366", "#A50026")
  two_colots <- c("#FDB366", "#A50026")
  
  #'  Function to reformat data for easier use with ggplot2 and plot marginal occupancy
  plot_conditional_occ <- function(predicted, spp1, spp2, x, covname, setup, spppair) {
    #'  Reformat data for ggplot
    #'  Snag conditional occupancy for each species, convert from wide to long format, 
    #'  and assign species name
    condish.occ <- as.data.frame(predicted[[2]][,,"mean"]) %>%
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
             Species = sub(" .*", "", Spp_Interaction)) %>%
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
    predicted.conditional.occ <- cbind(condish.occ, lower.condish[,2], upper.condish[,2], covs, scaled.covs)
    names(predicted.conditional.occ) <- c("Species", "Species_interaction", "conditional_occ", "lowerCRI", "upperCRI", "covs", "scaled_covs")
    
    #'  Plot species-specific conditional occupancy probabilities & 95% CRI
    condish_occ_plot <- ggplot(predicted.conditional.occ, aes(x = covs, y = conditional_occ, group = Species_interaction)) + 
      geom_line(aes(color = Species_interaction), lwd = 1.25) + 
      scale_color_manual(values = four_colors, labels = c(paste(spp1, "absent"), paste(spp1, "present"), paste(spp2, "absent"), paste(spp2, "present"))) + 
      #'  Add confidence intervals
      geom_ribbon(aes(ymin = lowerCRI, ymax = upperCRI, fill = Species_interaction), alpha = 0.3) +
      scale_fill_manual(values = four_colors, labels = c(paste(spp1, "absent"), paste(spp1, "present"), paste(spp2, "absent"), paste(spp2, "present"))) + 
      #'  Get rid of lines and gray background
      theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      #'  Force y-axis from 0 to 1
      ylim(0,1.0) +
      #'  Use list name as X-axis title
      xlab(x) +
      ylab("Conditional probability of occupancy") +
      labs(#title = paste("Conditional occupancy in response to", covname),
           fill = "Species interaction", color = "Species interaction") +
      facet_wrap(~Species, scales = "free_y") +
      theme(legend.position="bottom")
    #'  Review figure
    plot(condish_occ_plot)
    
    return(condish_occ_plot)
  }
  #####  Wolf-Bear conditional occupancy  ####
  # wolf.bear.condish.elev.ung <- plot_conditional_occ(predicted = wolf.bear.elev.ung.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Black Bear")
  # wolf.bear.condish.elev.pred <- plot_conditional_occ(predicted = wolf.bear.elev.pred.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Black Bear")
  # wolf.bear.condish.for.ung <- plot_conditional_occ(predicted = wolf.bear.for.ung.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Black Bear")
  # wolf.bear.condish.for.pred <- plot_conditional_occ(predicted = wolf.bear.for.pred.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "trail sites", spppair = "Wolf - Black Bear")
  # wolf.bear.condish.plots <- list(wolf.bear.condish.elev.ung, wolf.bear.condish.elev.pred, wolf.bear.condish.for.ung, wolf.bear.condish.for.pred)
  
  #####  Wolf-Coyote conditional occupancy  ####
  # wolf.coy.condish.elev.ung <- plot_conditional_occ(predicted = wolf.coy.elev.ung.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Coyote")
  # wolf.coy.condish.elev.pred <- plot_conditional_occ(predicted = wolf.coy.elev.pred.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Coyote")
  # wolf.coy.condish.for.ung <- plot_conditional_occ(predicted = wolf.coy.for.ung.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Coyote")
  # wolf.coy.condish.for.pred <- plot_conditional_occ(predicted = wolf.coy.for.pred.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Percent forest cover", setup = "trail sites", spppair = "Wolf - Coyote")
  # wolf.coy.condish.plots <- list(wolf.coy.condish.elev.ung, wolf.coy.condish.elev.pred, wolf.coy.condish.for.ung, wolf.coy.condish.for.pred)
  
  #####  Wolf-Lion conditional occupancy  ####
  #'  Nadda
  
  #####  Lion-Bear conditional occupancy  ####
  #'  Nadda
  
  #####  Lion-Bobcat conditional occupancy  ####
  #'  Nadda
  
  #####  Coyote-Bobcat conditional occupancy  ####
  coy.bob.condish.for.ung <- plot_conditional_occ(predicted = coy.bob.for.ung.yr2, spp1 = "Coyote", spp2 = "Bobcat", x = "Percent forest cover", covname = "forest cover", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.for.pred <- plot_conditional_occ(predicted = coy.bob.for.pred.yr2, spp1 = "Coyote", spp2 = "Bobcat", x = "Percent forest cover", covname = "forest cover", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.elev.ung <- plot_conditional_occ(predicted = coy.bob.elev.ung.yr2, spp1 = "Coyote", spp2 = "Bobcat", x = "Elevation (m)", covname = "elevation", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.elev.pred <- plot_conditional_occ(predicted = coy.bob.elev.pred.yr2, spp1 = "Coyote", spp2 = "Bobcat", x = "Elevation (m)", covname = "elevation", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.tri.ung <- plot_conditional_occ(predicted = coy.bob.tri.ung.yr2, spp1 = "Coyote", spp2 = "Bobcat", x = "Terrain ruggedness (TRI)", covname = "terrain ruggedness", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.tri.pred <- plot_conditional_occ(predicted = coy.bob.tri.pred.yr2, spp1 = "Coyote", spp2 = "Bobcat", x = "Terrain ruggedness (TRI)", covname = "terrain ruggedness", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.plots <- list(coy.bob.condish.for.ung, coy.bob.condish.for.pred, coy.bob.condish.elev.ung, coy.bob.condish.elev.pred, coy.bob.condish.tri.ung, coy.bob.condish.tri.pred)
  # coy.bob.condish.wtd.ung <- plot_conditional_occ(predicted = coy.bob.wtd.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "White-tailed deer relative abundance (RAI)", covname = "white-tailed deer relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  # coy.bob.condish.wtd.pred <- plot_conditional_occ(predicted = coy.bob.wtd.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "White-tailed deer relative abundance (RAI)", covname = "white-tailed deer relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  # coy.bob.condish.lago.ung <- plot_conditional_occ(predicted = coy.bob.lago.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Lagomorph relative abundance (RAI)", covname = "lagomorph relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  # coy.bob.condish.lago.pred <- plot_conditional_occ(predicted = coy.bob.lago.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Lagomorph relative abundance (RAI) ", covname = "lagomorph relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  # coy.bob.condish.div.ung <- plot_conditional_occ(predicted = coy.bob.div.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Shannon's diversity index (H)", covname = "prey diversity",setup = "random sites", spppair = "Coyote - Bobcat")
  # coy.bob.condish.div.pred <- plot_conditional_occ(predicted = coy.bob.div.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", x = "Shannon's diversity index (H)", covname = "prey diversity",setup = "trail sites", spppair = "Coyote - Bobcat")
  # coy.bob.condish.plots <- list(coy.bob.condish.elev.ung, coy.bob.condish.elev.pred, coy.bob.condish.for.ung, coy.bob.condish.for.pred, coy.bob.condish.wtd.ung, coy.bob.condish.wtd.pred, 
  #                               coy.bob.condish.lago.ung, coy.bob.condish.lago.pred, coy.bob.condish.div.ung, coy.bob.condish.div.pred)
  
  #####  Bear-Coyote conditional occupancy  ####
  bear.coy.condish.for.ung <- plot_conditional_occ(predicted = bear.coy.for.ung.yr2, spp1 = "Black bear", spp2 = "Coyote", x = "Percent forest cover", covname = "forest cover", setup = "random sites", spppair = "Black bear - Coyote")
  bear.coy.condish.for.pred <- plot_conditional_occ(predicted = bear.coy.for.pred.yr2, spp1 = "Black bear", spp2 = "Coyote", x = "Percent forest cover", covname = "forest cover", setup = "trail sites", spppair = "Black bear - Coyote")
  bear.coy.condish.elev.ung <- plot_conditional_occ(predicted = bear.coy.elev.ung.yr2, spp1 = "Black bear", spp2 = "Coyote", x = "Elevation (m)", covname = "elevation", setup = "random sites", spppair = "Black bear - Coyote")
  bear.coy.condish.elev.pred <- plot_conditional_occ(predicted = bear.coy.elev.pred.yr2, spp1 = "Black bear", spp2 = "Coyote", x = "Elevation (m)", covname = "elevation", setup = "trail sites", spppair = "Black bear - Coyote")
  bear.coy.condish.tri.ung <- plot_conditional_occ(predicted = bear.coy.tri.ung.yr2, spp1 = "Black bear", spp2 = "Coyote", x = "Terrain ruggedness (TRI)", covname = "terrain ruggedness", setup = "random sites", spppair = "Black bear - Coyote")
  bear.coy.condish.tri.pred <- plot_conditional_occ(predicted = bear.coy.tri.pred.yr2, spp1 = "Black bear", spp2 = "Coyote", x = "Terrain ruggedness (TRI)", covname = "terrain ruggedness", setup = "trail sites", spppair = "Black bear - Coyote")
  bear.coy.condish.plots <- list(bear.coy.condish.for.ung, bear.coy.condish.for.pred, bear.coy.condish.elev.ung, bear.coy.condish.elev.pred, bear.coy.condish.tri.ung, bear.coy.condish.tri.pred)
  
  #'  Read in predator silhoettes
  boburl <- "https://images.phylopic.org/images/ab6cfd4f-aef7-40fa-b5a5-1b79b7d112aa/raster/1024x740.png?v=17a30c18af1.png"
  bobimg <- readPNG(getURLContent(boburl), native = T)
  bobgrid <- rasterGrob(bobimg, interpolate = TRUE)
  coyurl <- "https://images.phylopic.org/images/5a0398e3-a455-4ca6-ba86-cf3f1b25977a/raster/1024x894.png?v=16fe8749858.png"
  coyimg <- readPNG(getURLContent(coyurl), native = T) 
  coygrid <- rasterGrob(coyimg, interpolate = TRUE)
  coyurlGB <- "https://images.phylopic.org/images/e6a2fa4b-85df-43b4-989c-34a65ba7eee3/raster/1024x911.png?v=17f2638df97.png"
  coyimgGB <- readPNG(getURLContent(coyurlGB), native = T)
  coygridGB <- rasterGrob(coyimg, interpolate = TRUE)
  
  #'  Patchwork of conditional occupancy plots
  patchwork_conditional_coybob_pred <- coy.bob.condish.for.pred +
    inset_element(p = bobimg, left = 0.25, bottom = 0.91, right = 0.5, top = 1.17, ignore_tag = TRUE) +
    theme(rect = element_rect(fill = "transparent", linetype = "blank")) +
    inset_element(p = coyimgGB, left = 0.75, bottom = 0.9, right = 1.09, top = 1.18, ignore_tag = TRUE) +
    theme(rect = element_rect(fill = "transparent", linetype = "blank")) + #coy.bob.condish.elev.pred + 
    coy.bob.condish.tri.pred + 
    plot_layout(nrow = 2) + plot_annotation(tag_levels = 'a',
                                            title = "Habitat effects on co-occurrence probabilities") +
    plot_layout(guides = "collect") & theme(legend.position = "bottom") 
  patchwork_conditional_coybob_pred
  
  patchwork_conditional_coybob_ung <- coy.bob.condish.for.ung +
    inset_element(p = bobimg, left = 0.25, bottom = 0.91, right = 0.5, top = 1.17, ignore_tag = TRUE) +
    theme(rect = element_rect(fill = "transparent", linetype = "blank")) +
    inset_element(p = coyimgGB, left = 0.75, bottom = 0.9, right = 1.09, top = 1.18, ignore_tag = TRUE) +
    theme(rect = element_rect(fill = "transparent", linetype = "blank")) + #coy.bob.condish.elev.ung + 
    coy.bob.condish.tri.ung + 
    plot_layout(nrow = 2) + plot_annotation(tag_levels = 'a',
                                            title = "Habitat effects on co-occurrence probabilities") +
    plot_layout(guides = "collect") & theme(legend.position = "bottom") 
  patchwork_conditional_coybob_ung
  
  
  #'  Coyote-bobcat 
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_for_ungcam.tiff", coy.bob.condish.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_for_predcam.tiff", coy.bob.condish.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_elev_ungcam.tiff", coy.bob.condish.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_elev_predcam.tiff", coy.bob.condish.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_tri_ungcam.tiff", coy.bob.condish.plots[[5]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_tri_predcam.tiff", coy.bob.condish.plots[[6]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Co-occurrence plots
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_predcams.tiff", patchwork_conditional_coybob_pred, 
         units = "in", width = 6, height = 7, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/conditional_occ_plots_coybob_ungcams.tiff", patchwork_conditional_coybob_ung, 
         units = "in", width = 6, height = 7, dpi = 600, device = 'tiff', compression = 'lzw')
  
    #' #'  Effect of low and high WTD RAI values on coyote-bobcat conditional occupancy
  #' #'  with increasing lagomorph activity
  #' load("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/Predicted_psi-cov_coy-bob-lago_LowHiWTD_2023-07-13.RData")
  #' 
  #' #'  Function to reformat data for easier use with ggplot2 and plot conditional occupancy
  #' plot_conditional_occ_singlespp <- function(predicted, spp1, spp2, x, ylabtitle, spppair) {
  #'   #'  Reformat data for ggplot
  #'   #'  Snag conditional occupancy, convert from wide to long format, and assign species name
  #'   condish.occ <- as.data.frame(predicted[[2]][,,"mean"]) %>%
  #'     pivot_longer(cols = c(Spp1.alone, Spp1.given.Spp2, Spp2.alone, Spp2.given.Spp1), 
  #'                  names_to = "Spp_Interaction") %>%
  #'     arrange(Spp_Interaction) %>%
  #'     #'  Make these labels meaningful with actual species names
  #'     mutate(Spp_Interaction = factor(Spp_Interaction, levels = c("Spp1.alone", "Spp1.given.Spp2", "Spp2.alone", "Spp2.given.Spp1")),
  #'            Spp_Interaction = gsub("Spp1", spp1, Spp_Interaction),
  #'            Spp_Interaction = gsub("Spp2", spp2, Spp_Interaction),
  #'            Spp_Interaction = gsub(".alone", " alone", Spp_Interaction),
  #'            Spp_Interaction = gsub(".given.", " given ", Spp_Interaction), 
  #'            Spp_Interaction = str_to_sentence(Spp_Interaction),
  #'            Species = sub(" .*", "", Spp_Interaction)) %>%
  #'     relocate(Species, .before = Spp_Interaction)
  #'   #'  Snag lower 95% credible interval
  #'   lower.condish <- as.data.frame(predicted[[2]][,,"lower"]) %>%
  #'     pivot_longer(cols = c(Spp1.alone, Spp1.given.Spp2, Spp2.alone, Spp2.given.Spp1), names_to = "Spp_Interaction") %>%
  #'     arrange(Spp_Interaction)
  #'   #'  Snag upper 95% credible interval 
  #'   upper.condish <- as.data.frame(predicted[[2]][,,"upper"]) %>%
  #'     pivot_longer(cols = c(Spp1.alone, Spp1.given.Spp2, Spp2.alone, Spp2.given.Spp1), names_to = "Spp_Interaction") %>%
  #'     arrange(Spp_Interaction)
  #'   #'  Snag covariate values
  #'   covs <- c(predicted[[3]], predicted[[3]]) 
  #'   scaled.covs <- c(predicted[[4]], predicted[[4]])
  #'   #'  Create single data frame
  #'   predicted.conditional.occ <- cbind(condish.occ, lower.condish[,2], upper.condish[,2], covs, scaled.covs)
  #'   names(predicted.conditional.occ) <- c("Species", "Species_interaction", "conditional_occ", "lowerCRI", "upperCRI", "covs", "scaled_covs")
  #'   
  #'   #'  Filter to single species
  #'   predicted.conditional.occ <- filter(predicted.conditional.occ, Species == spp1)
  #'   
  #'   #'  Plot species-specific conditional occupancy probabilities & 95% CRI
  #'   condish_occ_plot <- ggplot(predicted.conditional.occ, aes(x = covs, y = conditional_occ, group = Species_interaction)) + 
  #'     geom_line(aes(color = Species_interaction), lwd = 1.25) + 
  #'     scale_color_manual(values = two_colots, labels = c(paste(spp2, "absent"), paste(spp2, "present"))) + 
  #'     #'  Add confidence intervals
  #'     geom_ribbon(aes(ymin = lowerCRI, ymax = upperCRI, fill = Species_interaction), alpha = 0.3) +
  #'     scale_fill_manual(values = two_colots, labels = c(paste(spp2, "absent"), paste(spp2, "present"))) + 
  #'     #'  Get rid of lines and gray background
  #'     theme_bw() +
  #'     theme(panel.border = element_blank()) +
  #'     theme(axis.line = element_line(color = 'black')) +
  #'     #'  Force y-axis from 0 to 1
  #'     ylim(0,1.0) +
  #'     #'  Use list name as X-axis title
  #'     xlab(x) +
  #'     ylab(ylabtitle) +
  #'     labs(#title = paste("Conditional occupancy in response to", covname),
  #'       fill = "Species interaction", color = "Species interaction") +
  #'     theme(legend.position="bottom")
  #'   #'  Review figure
  #'   plot(condish_occ_plot)
  #'   
  #'   return(condish_occ_plot)
  #' }
  #' coy.bob.condish.lago.pred.lowWTD <- plot_conditional_occ_singlespp(predicted = coy.bob.lago.pred.yr1_LowHiWTD[[1]], spp1 = "Coyote", spp2 = "Bobcat", x = "Lagomorph relative abundance (RAI), \nwhite-tailed deer RAI = 0", ylabtitle = "Conditional occupancy probability, trail sites", spppair = "Coyote - Bobcat")
  #' coy.bob.condish.lago.pred.medWTD <- plot_conditional_occ_singlespp(predicted = coy.bob.lago.pred.yr1_LowHiWTD[[2]], spp1 = "Coyote", spp2 = "Bobcat", x = "Lagomorph relative abundance (RAI), \nwhite-tailed deer RAI = 2", ylabtitle = "", spppair = "Coyote - Bobcat")
  #' coy.bob.condish.lago.pred.hiWTD <- plot_conditional_occ_singlespp(predicted = coy.bob.lago.pred.yr1_LowHiWTD[[3]], spp1 = "Coyote", spp2 = "Bobcat", x = "Lagomorph relative abundance (RAI), \nwhite-tailed deer RAI = 6", ylabtitle = "", spppair = "Coyote - Bobcat")
  #' 
  #' #'  Patchwork plots into single figure
  #' coy.bob.condish.lago.pred.wtd <- coy.bob.condish.lago.pred.lowWTD + coy.bob.condish.lago.pred.medWTD + coy.bob.condish.lago.pred.hiWTD +
  #'   plot_annotation(title = "Effect of white-tailed deer on coyote site use, conditional on bobcat site use and lagomorph abundance") +
  #'   plot_annotation(tag_levels = 'a') + 
  #'   plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  
  #' #'  -----------------------------
  #' ####  Save all the pretty plots  ####
  #' #'  -----------------------------
  #' #'  Marginal occupancy probabilities
  #' #'  Wolf-bear
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_ung_marginal_occ_plots.tiff", wolf.bear.marg.plots[[1]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_pred_marginal_occ_plots.tiff", wolf.bear.marg.plots[[2]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_ung_marginal_occ_plots.tiff", wolf.bear.marg.plots[[3]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_pred_marginal_occ_plots.tiff", wolf.bear.marg.plots[[4]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_div_ung_marginal_occ_plots.tiff", wolf.bear.marg.plots[[5]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_div_pred_marginal_occ_plots.tiff", wolf.bear.marg.plots[[6]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #'  Wolf-coyote
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_ung_marginal_occ_plots.tiff", wolf.coy.marg.plots[[1]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_pred_marginal_occ_plots.tiff", wolf.coy.marg.plots[[2]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_ung_marginal_occ_plots.tiff", wolf.coy.marg.plots[[3]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_pred_marginal_occ_plots.tiff", wolf.coy.marg.plots[[4]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #'  Coyote-bobcat
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_elev_ung_marginal_occ_plots.tiff", coy.bob.marg.plots[[1]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_elev_pred_marginal_occ_plots.tiff", coy.bob.marg.plots[[2]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_for_ung_marginal_occ_plots.tiff", coy.bob.marg.plots[[3]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_for_pred_marginal_occ_plots.tiff", coy.bob.marg.plots[[4]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #'  Combined by covariate
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_occ_elevation_plot.tiff", marginal_elev_plot, 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_occ_forestcov_plot.tiff", marginal_for_plot, 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_occ_preydiversity_plot.tiff", marginal_dif_plot, 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_occ_allcov_trailsites_plot.tiff", marginal_pred_patchwork, 
  #'        units = "in", width = 9, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_occ_allcov_randsites_plot.tiff", marginal_ung_patchwork, 
  #'        units = "in", width = 9, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' 
  #' 
  #' #' #'  Conditional occupancy probabilities
  #' #' #'  Wolf-bear  
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_ung_conditional_occ_plots.tiff", wolf.bear.condish.plots[[1]],
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_pred_conditional_occ_plots.tiff", wolf.bear.condish.plots[[2]],
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_ung_conditional_occ_plots.tiff", wolf.bear.condish.plots[[3]],
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_pred_conditional_occ_plots.tiff", wolf.bear.condish.plots[[4]],
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #' #'  Wolf-coyote 
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_ung_conditional_occ_plots.tiff", wolf.coy.condish.plots[[1]], 
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_pred_conditional_occ_plots.tiff", wolf.coy.condish.plots[[2]], 
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_ung_conditional_occ_plots.tiff", wolf.coy.condish.plots[[3]], 
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_pred_conditional_occ_plots.tiff", wolf.coy.condish.plots[[4]], 
  #' #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #'  Coyote-bobcat 
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_wtd_ung_conditional_occ_plots.tiff", coy.bob.condish.plots[[5]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_wtd_pred_conditional_occ_plots.tiff", coy.bob.condish.plots[[6]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_ung_conditional_occ_plots.tiff", coy.bob.condish.plots[[7]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_pred_conditional_occ_plots.tiff", coy.bob.condish.plots[[8]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_div_ung_conditional_occ_plots.tiff", coy.bob.condish.plots[[9]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_div_pred_conditional_occ_plots.tiff", coy.bob.condish.plots[[10]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' 
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_pred_conditional_occ_plots_wtd0.tiff", coy.bob.condish.lago.pred.lowWTD, 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_pred_conditional_occ_plots_wtd2.tiff", coy.bob.condish.lago.pred.medWTD, 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_pred_conditional_occ_plots_wtd6.tiff", coy.bob.condish.lago.pred.hiWTD, 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_pred_conditional_occ_plots_wtd0-6.tiff", coy.bob.condish.lago.pred.wtd, 
  #'        units = "in", width = 10, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  
  #'  ----------------------------------------
  ####  Predict mean occupancy for Yr1 & Yr2  ####
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
  
  #'  Predict mean occupancy for Yr1 and Yr2
  wolf.bear.mean.yr1 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.mean.yr2 <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.coy.mean.yr1 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 0,
                                         psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.coy.mean.yr2 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 0,
                                         psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.lion.mean.yr1 <- predict_occupancy(mod = wolf.lion.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 0), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.lion.mean.yr2 <- predict_occupancy(mod = wolf.lion.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 1), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  lion.bear.mean.yr1 <- predict_occupancy(mod = lion.bear.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 0), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  lion.bear.mean.yr2 <- predict_occupancy(mod = lion.bear.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 1), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  lion.bob.mean.yr1 <- predict_occupancy(mod = lion.bob.null, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         psi_cov = c(1, 0), psi_cov_index = 0,
                                         psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  lion.bob.mean.yr2 <- predict_occupancy(mod = lion.bob.null, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         psi_cov = c(1, 1), psi_cov_index = 0,
                                         psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  coy.bob.mean.yr1 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                        focal_cov = stations_skinny_eoe20s21s$Year,
                                        psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 0,
                                        psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  coy.bob.mean.yr2 <- predict_occupancy(mod = coy.bob.habx, ncat = 4, npoints = 500,
                                        focal_cov = stations_skinny_eoe20s21s$Year,
                                        psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 0,
                                        psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  bear.coy.mean.yr1 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                        focal_cov = stations_skinny_eoe20s21s$Year,
                                        psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 0,
                                        psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  bear.coy.mean.yr2 <- predict_occupancy(mod = bear.coy.habx, ncat = 4, npoints = 500,
                                        focal_cov = stations_skinny_eoe20s21s$Year,
                                        psi_cov = c(1, 1, 1, 0, 0, 0), psi_cov_index = 0,
                                        psi_inxs_cov = c(1), psi_inxs_cov_index = 0)
  
  #'  Extract mean and 95% CRI for each observation (should be identical) and thin
  #'  to one value per species
  mean_occ_prob <- function(predicted, spp1, spp2, yr) {
    #'  Snag mean
    marg.occ <- as.data.frame(predicted[[1]][,,"mean"]) %>%
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
    predicted.marginal.occ <- cbind(marg.occ, lower.marg[,2], upper.marg[,2])
    names(predicted.marginal.occ) <- c("Species", "marginal_occ", "lowerCRI", "upperCRI")
    #'  Drop extra observations and round
    predicted.marginal.occ <- as.data.frame(predicted.marginal.occ) %>%
      group_by(Species) %>%
      slice(1L) %>%
      ungroup() %>%
      mutate(marginal_occ = round(marginal_occ, 2),
             lowerCRI = round(lowerCRI, 2),
             upperCRI = round(upperCRI, 2),
             Year = yr)
    return(predicted.marginal.occ)
  }
  psimean_wolf.bear.yr1 <- mean_occ_prob(wolf.bear.mean.yr1, spp1 = "Wolf", spp2 = "Black bear", yr = "2020")
  psimean_wolf.bear.yr2 <- mean_occ_prob(wolf.bear.mean.yr2, spp1 = "Wolf", spp2 = "Black bear", yr = "2021")
  psimean_wolf.coy.yr1 <- mean_occ_prob(wolf.coy.mean.yr1, spp1 = "Wolf", spp2 = "Coyote", yr = "2020")
  psimean_wolf.coy.yr2 <- mean_occ_prob(wolf.coy.mean.yr2, spp1 = "Wolf", spp2 = "Coyote", yr = "2021")
  psimean_wolf.lion.yr1 <- mean_occ_prob(wolf.lion.mean.yr1, spp1 = "Wolf", spp2 = "Mountain lion", yr = "2020")
  psimean_wolf.lion.yr2 <- mean_occ_prob(wolf.lion.mean.yr2, spp1 = "Wolf", spp2 = "Mountain lion", yr = "2021")
  psimean_lion.bear.yr1 <- mean_occ_prob(lion.bear.mean.yr1, spp1 = "Mountain lion", spp2 = "Black bear", yr = "2020")
  psimean_lion.bear.yr2 <- mean_occ_prob(lion.bear.mean.yr2, spp1 = "Mountain lion", spp2 = "Black bear", yr = "2021")
  psimean_lion.bob.yr1 <- mean_occ_prob(lion.bob.mean.yr1, spp1 = "Mountain lion", spp2 = "Bobcat", yr = "2020")
  psimean_lion.bob.yr2 <- mean_occ_prob(lion.bob.mean.yr2, spp1 = "Mountain lion", spp2 = "Bobcat", yr = "2021")
  psimean_coy.bob.yr1 <- mean_occ_prob(coy.bob.mean.yr1, spp1 = "Coyote", spp2 = "Bobcat", yr = "2020")
  psimean_coy.bob.yr2 <- mean_occ_prob(coy.bob.mean.yr2, spp1 = "Coyote", spp2 = "Bobcat", yr = "2021")
  psimean_bear.coy.yr1 <- mean_occ_prob(bear.coy.mean.yr1, spp1 = "Black bear", spp2 = "Coyote", yr = "2020")
  psimean_bear.coy.yr2 <- mean_occ_prob(bear.coy.mean.yr2, spp1 = "Black bear", spp2 = "Coyote", yr = "2021")
  
  #'  Bind all species and years together
  mean_annual_psi <- rbind(psimean_wolf.bear.yr1, psimean_wolf.bear.yr2, psimean_wolf.lion.yr1, 
                           psimean_wolf.lion.yr2, psimean_wolf.coy.yr1, psimean_wolf.coy.yr2, 
                           psimean_lion.bear.yr1, psimean_lion.bear.yr2, psimean_lion.bob.yr1,
                           psimean_lion.bob.yr2, psimean_coy.bob.yr1, psimean_coy.bob.yr2,
                           psimean_bear.coy.yr1, psimean_bear.coy.yr2) %>%
    #'  Reduce to one observation per species
    group_by(Species, Year) %>%
    slice(1L) %>%
    ungroup() %>%
    #'  Combine mean & CRI into single column
    mutate(CRI = paste0("(", lowerCRI, ", ", upperCRI, ")"),
           Mean = paste(marginal_occ, CRI)) %>%
    dplyr::select(Species, Year, Mean) %>%
    spread(Year, Mean) %>%
    rename("Mean_CRI_2020" = "2020", "Mean_CRI_2021" = "2021")
  
  #'  Save!
  write.csv(mean_annual_psi, "./Outputs/Tables/Summary_mean_psi_yr1v2.csv")
  
  
  
