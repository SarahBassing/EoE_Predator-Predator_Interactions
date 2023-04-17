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
  # library(grid)
  # library(png)
  # library(RCurl)
  library(rphylopic)
  library(tidyverse)
  
  #'  Load covariate data
  load("./Data/Covariates_extracted/Covariate_skinny_EoE20s21s.RData")
  
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
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-11.RData")
  wolf.bear.elev.ung.yr1 <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                               focal_cov = stations_skinny_eoe20s21s$Elev,
                                               psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 4,
                                               psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.elev.pred.yr1 <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                                focal_cov = stations_skinny_eoe20s21s$Elev,
                                                psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 4,
                                                psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.for.ung.yr1 <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$PercForest,
                                              psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 5,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.for.pred.yr1 <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                               focal_cov = stations_skinny_eoe20s21s$PercForest,
                                               psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 5,
                                               psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.div.ung.yr1 <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$SppDiversity,
                                              psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 6,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.div.pred.yr1 <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                               focal_cov = stations_skinny_eoe20s21s$SppDiversity,
                                               psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 6,
                                               psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  
  #####  Wolf-Coyote predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2023-04-10.RData") 
  wolf.coy.elev.ung.yr1 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 0, 0, 0, 0), psi_cov_index = 4,
                                             psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.elev.pred.yr1 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 1, 0, 0, 0), psi_cov_index = 4,
                                             psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.for.ung.yr1 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                           focal_cov = stations_skinny_eoe20s21s$PercForest,
                                           psi_cov = c(1, 0, 0, 0, 0), psi_cov_index = 5,
                                           psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  wolf.coy.for.pred.yr1 <- predict_occupancy(mod = wolf.coy.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$PercForest,
                                             psi_cov = c(1, 1, 0, 0, 0), psi_cov_index = 5,
                                             psi_inxs_cov = 0, psi_inxs_cov_index = 0)
  
  #####  Wolf-Lion predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2023-04-06.RData")
  
  
  #####  Lion-Bear predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_2023-04-06.RData")
  
  
  #####  Lion-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_2023-04-06.RData")
  
  
  #####  Coyote-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_2023-04-11.RData")
  coy.bob.elev.ung.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500, 
                                            focal_cov = stations_skinny_eoe20s21s$Elev,
                                            psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), psi_cov_index = 4,
                                            psi_inxs_cov = c(1, 0, 0, 0, 0), psi_inxs_cov_index = 0)
  coy.bob.elev.pred.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), psi_cov_index = 4,
                                             psi_inxs_cov = c(1, 1, 0, 0, 0), psi_inxs_cov_index = 0)
  coy.bob.for.ung.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500,
                                           focal_cov = stations_skinny_eoe20s21s$PercForest,
                                           psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), psi_cov_index = 5,
                                           psi_inxs_cov = c(1, 0, 0, 0, 0), psi_inxs_cov_index = 0)
  coy.bob.for.pred.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500,
                                            focal_cov = stations_skinny_eoe20s21s$PercForest,
                                            psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), psi_cov_index = 5,
                                            psi_inxs_cov = c(1, 1, 0, 0, 0), psi_inxs_cov_index = 0)
  coy.bob.wtd.ung.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500, 
                                            focal_cov = stations_skinny_eoe20s21s$Nwtd,
                                            psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), psi_cov_index = 6,
                                            psi_inxs_cov = c(1, 0, 0, 0, 0), psi_inxs_cov_index = 3) # applying same cov to differ places in psi[1:8] & psix[1:5] covariate matrix
  coy.bob.wtd.pred.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Nwtd,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), psi_cov_index = 6,
                                             psi_inxs_cov = c(1, 1, 0, 0, 0), psi_inxs_cov_index = 3)
  coy.bob.lago.ung.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500, 
                                            focal_cov = stations_skinny_eoe20s21s$Nlagomorph,
                                            psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), psi_cov_index = 7,
                                            psi_inxs_cov = c(1, 0, 0, 0, 0), psi_inxs_cov_index = 4)
  coy.bob.lago.pred.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Nlagomorph,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), psi_cov_index = 7,
                                             psi_inxs_cov = c(1, 1, 0, 0, 0), psi_inxs_cov_index = 4)
  coy.bob.div.ung.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500, 
                                            focal_cov = stations_skinny_eoe20s21s$SppDiversity,
                                            psi_cov = c(1, 0, 0, 0, 0, 0, 0, 0), psi_cov_index = 8,
                                            psi_inxs_cov = c(1, 0, 0, 0, 0), psi_inxs_cov_index = 5)
  coy.bob.div.pred.yr1 <- predict_occupancy(mod = coy.bob.global, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$SppDiversity,
                                             psi_cov = c(1, 1, 0, 0, 0, 0, 0, 0), psi_cov_index = 8,
                                             psi_inxs_cov = c(1, 1, 0, 0, 0), psi_inxs_cov_index = 5)
  
  save.image(file = paste0("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/Predicted_psi-cov_relationships_", Sys.Date(), ".RData"))
  
  load("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/Predicted_psi-cov_relationships_2023-04-14.RData")
  
  #'  --------------------------------
  ####  Plot marginal Pr(occupancy)  ####
  #'  --------------------------------
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
  #'  Top model
  wolf.bear.marg.elev.ung <- plot_marginal_occ(predicted = wolf.bear.elev.ung.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.elev.pred <- plot_marginal_occ(predicted = wolf.bear.elev.pred.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.for.ung <- plot_marginal_occ(predicted = wolf.bear.for.ung.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.for.pred <- plot_marginal_occ(predicted = wolf.bear.for.pred.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "Trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.div.ung <- plot_marginal_occ(predicted = wolf.bear.div.ung.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Shannon's diversity index", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.div.pred <- plot_marginal_occ(predicted = wolf.bear.div.pred.yr1, spp1 = "Wolf", spp2 = "Bear", covname = "Shannon's diversity index", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.plots <- list(wolf.bear.marg.elev.ung, wolf.bear.marg.elev.pred, wolf.bear.marg.for.ung, wolf.bear.marg.for.pred, wolf.bear.marg.div.ung, wolf.bear.marg.div.pred)
  
  # wolf.bear.marg.elev.ung_a <- plot_marginal_occ(predicted = wolf.bear.elev.ung.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Black Bear")
  # wolf.bear.marg.elev.pred_a <- plot_marginal_occ(predicted = wolf.bear.elev.pred.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Black Bear")
  # wolf.bear.marg.for.ung_a <- plot_marginal_occ(predicted = wolf.bear.for.ung.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Black Bear")
  # wolf.bear.marg.for.pred_a <- plot_marginal_occ(predicted = wolf.bear.for.pred.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "Trail sites", spppair = "Wolf - Black Bear")
  # wolf.bear.marg.plots_altmod <- list(wolf.bear.marg.elev.ung_a, wolf.bear.marg.elev.pred_a, wolf.bear.marg.for.ung_a, wolf.bear.marg.for.pred_a)
  
  #####  Wolf-Coyote marginal occupancy  ####
  wolf.coy.marg.elev.ung <- plot_marginal_occ(predicted = wolf.coy.elev.ung.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.elev.pred <- plot_marginal_occ(predicted = wolf.coy.elev.pred.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.for.ung <- plot_marginal_occ(predicted = wolf.coy.for.ung.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Coyote") 
  wolf.coy.marg.for.pred <- plot_marginal_occ(predicted = wolf.coy.for.pred.yr1, spp1 = "Wolf", spp2 = "Coyote", covname = "Percent forest cover", setup = "trail sites", spppair = "Wolf - Coyote")
  wolf.coy.marg.plots <- list(wolf.coy.marg.elev.ung, wolf.coy.marg.elev.pred, wolf.coy.marg.for.ung, wolf.coy.marg.for.pred)
  
  #####  Wolf-Lion marginal occupancy  ####
  #'  Nadda
  
  #####  Lion-Bear marginal occupancy  ####
  #'  Nadda
  
  #####  Lion-Bobcat marginal occupancy  ####
  #'  Nadda
  
  #####  Coyote-Bobcat marginal occupancy  ####
  coy.bob.marg.elev.ung <- plot_marginal_occ(predicted = coy.bob.elev.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Elevation (m)", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.elev.pred <- plot_marginal_occ(predicted = coy.bob.elev.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Elevation (m)", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.for.ung <- plot_marginal_occ(predicted = coy.bob.for.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Percent forest cover", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.for.pred <- plot_marginal_occ(predicted = coy.bob.for.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Percent forest cover", setup = "trail sites", spppair = "Coyote - Bobcat")
  # coy.bob.marg.wtd.ung <- plot_marginal_occ(predicted = coy.bob.wtd.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "White-tailed deer relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  # coy.bob.marg.wtd.pred <- plot_marginal_occ(predicted = coy.bob.wtd.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "White-tailed deer relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  # coy.bob.marg.lago.ung <- plot_marginal_occ(predicted = coy.bob.lago.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Lagomorph relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  # coy.bob.marg.lago.pred <- plot_marginal_occ(predicted = coy.bob.lago.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Lagomorph relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  # coy.bob.marg.div.ung <- plot_marginal_occ(predicted = coy.bob.div.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Shannon's diversity index", setup = "random sites", spppair = "Coyote - Bobcat")
  # coy.bob.marg.div.pred <- plot_marginal_occ(predicted = coy.bob.div.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Shannon's diversity index", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.plots <- list(coy.bob.marg.elev.ung, coy.bob.marg.elev.pred, coy.bob.marg.for.ung, coy.bob.marg.for.pred)#, 
                             #coy.bob.marg.wtd.ung, coy.bob.marg.wtd.pred, coy.bob.marg.lago.ung, coy.bob.marg.lago.pred, coy.bob.marg.div.ung, coy.bob.marg.div.pred)
 
  
  #'  Add species pair and camera placement to each data set
  add_info <- function(marg, spp_pair, cam_setup) {
    marg$Species_pair <- spp_pair
    marg$Setup <- cam_setup
    return(marg)
  }
  wolf.bear.marg.elev.ung_new <- add_info(wolf.bear.marg.elev.ung$data, spp_pair = "wolf-bear", cam_setup = "Random sites")
  wolf.bear.marg.elev.pred_new <- add_info(wolf.bear.marg.elev.pred$data, spp_pair = "wolf-bear", cam_setup = "Trail sites")
  wolf.bear.marg.for.ung_new <- add_info(wolf.bear.marg.for.ung$data, spp_pair = "wolf-bear", cam_setup = "Random sites")
  wolf.bear.marg.for.pred_new <- add_info(wolf.bear.marg.for.pred$data, spp_pair = "wolf-bear", cam_setup = "Trail sites")
  wolf.bear.marg.div.ung_new <- add_info(wolf.bear.marg.div.ung$data, spp_pair = "wolf-bear", cam_setup = "Random sites")
  wolf.bear.marg.div.pred_new <- add_info(wolf.bear.marg.div.pred$data, spp_pair = "wolf-bear", cam_setup = "Trail sites")
  wolf.coy.marg.elev.ung_new <- add_info(wolf.coy.marg.elev.ung$data, spp_pair = "wolf-coy", cam_setup = "Random sites")
  wolf.coy.marg.elev.pred_new <- add_info(wolf.coy.marg.elev.pred$data, spp_pair = "wolf-coy", cam_setup = "Trail sites")
  wolf.coy.marg.for.ung_new <- add_info(wolf.coy.marg.for.ung$data, spp_pair = "wolf-coy", cam_setup = "Random sites")
  wolf.coy.marg.for.pred_new <- add_info(wolf.coy.marg.for.pred$data, spp_pair = "wolf-coy", cam_setup = "Trail sites")
  coy.bob.marg.elev.ung_new <- add_info(coy.bob.marg.elev.ung$data, spp_pair = "coy-bob", cam_setup = "Random sites")
  coy.bob.marg.elev.pred_new <- add_info(coy.bob.marg.elev.pred$data, spp_pair = "coy-bob", cam_setup = "Trail sites")
  coy.bob.marg.for.ung_new <- add_info(coy.bob.marg.for.ung$data, spp_pair = "coy-bob", cam_setup = "Random sites")
  coy.bob.marg.for.pred_new <- add_info(coy.bob.marg.for.pred$data, spp_pair = "coy-bob", cam_setup = "Trail sites")
  
  #'  Combine all marginal probability estimates for each significant variable
  marginal_elev <- rbind(wolf.bear.marg.elev.ung_new, wolf.bear.marg.elev.pred_new)
  marginal_for <- rbind(wolf.bear.marg.for.ung_new, wolf.bear.marg.for.pred_new, 
                        wolf.coy.marg.for.ung_new, wolf.coy.marg.for.pred_new, 
                        coy.bob.marg.for.ung_new, coy.bob.marg.for.pred_new) %>%
    #'  Drop wolf-coy results b/c already have wolf & coy results included w/ other pairings
    filter(Species_pair != "wolf-coy")
  marginal_div <- rbind(wolf.bear.marg.div.ung_new, wolf.bear.marg.div.pred_new) %>%
    #'  Drop wolf results since prey diversity not a meaningful relationship
    filter(Species != "Wolf")
  
  #'  Set color combos for each species
  #'  wolf = "#364B9A", bear = "#98CAE1", lion = "#FDB366", coyote = "#DD3D2D", bobcat = "#A50026"
  bear_colors <- "#98CAE1"
  wolf.bear_colors <- c("#364B9A", "#98CAE1")
  wolf.bear.coy_colors <- c("#364B9A", "#98CAE1", "#FDB366")
  wolf.bear.coy.bob_colors <- c("#364B9A", "#98CAE1", "#FDB366", "#A50026")
  wolf.bear.lion.coy.bob_colors <- c("#364B9A", "#98CAE1", "#DD3D2D", "#FDB366", "#A50026") 
  
  #'  Plot each species response to specific covariate together
  plot_margingal_occ_by_cov <- function(predicted, x, covname, ncolor) {
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
      labs(title = paste("Species-specific marginal occupancy", covname), 
           fill = "Species", color = "Species") +
      facet_wrap(~Setup, scales = "free_y") +
      theme(legend.position="bottom")
    #'  Review figure
    plot(marg_occ_plot)
    
    return(marg_occ_plot)
  }
  marginal_elev_plot <- plot_margingal_occ_by_cov(marginal_elev, x = "Elevation (m)", covname = "across elevation", ncolor = wolf.bear_colors)
  marginal_for_plot <- plot_margingal_occ_by_cov(marginal_for, x = "Percent forest cover", covname = "across percent forest cover", ncolor = wolf.bear.coy.bob_colors)
  marginal_dif_plot <- plot_margingal_occ_by_cov(marginal_div, x = "Shannon's diversity index (H)", covname = "and prey diversity", ncolor = bear_colors)
  
  
  #'  -----------------------------------
  ####  Plot conditional Pr(occupancy)  ####
  #'  -----------------------------------
  #'  Color-blind friendly color palette from Khroma
  plot_scheme(colour("sunset")(11))
  colour("sunset")(11)
  # four_colors <- c("#364B9A", "#C2E4EF", "#FEDA8B", "#A50026")
  four_colors <- c("#364B9A", "#98CAE1", "#FDB366", "#A50026")
  
  #'  Function to reformat data for easier use with ggplot2 and plot marginal occupancy
  plot_conditional_occ <- function(predicted, spp1, spp2, covname, setup, spppair) {
    #'  Reformat data for ggplot
    #'  Snag marginal occupancy for each species, convert from wide to long format, 
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
    
    #'  Plot species-specific marginal occupancy probabilities & 95% CRI
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
      xlab(covname) +
      ylab(paste("Conditional occupancy probability,", setup)) +
      labs(#title = paste(spppair, "Conditional Occupancy Probabilities"),
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
  coy.bob.condish.elev.ung <- plot_conditional_occ(predicted = coy.bob.elev.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Elevation (m)", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.elev.pred <- plot_conditional_occ(predicted = coy.bob.elev.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Elevation (m)", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.for.ung <- plot_conditional_occ(predicted = coy.bob.for.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Percent forest cover", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.for.pred <- plot_conditional_occ(predicted = coy.bob.for.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Percent forest cover", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.wtd.ung <- plot_conditional_occ(predicted = coy.bob.wtd.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "White-tailed deer relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.wtd.pred <- plot_conditional_occ(predicted = coy.bob.wtd.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "White-tailed deer relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.lago.ung <- plot_conditional_occ(predicted = coy.bob.lago.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Lagomorph relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.lago.pred <- plot_conditional_occ(predicted = coy.bob.lago.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Lagomorph relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.div.ung <- plot_conditional_occ(predicted = coy.bob.div.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Shannon's diversity index", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.div.pred <- plot_conditional_occ(predicted = coy.bob.div.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Shannon's diversity index", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.condish.plots <- list(coy.bob.condish.elev.ung, coy.bob.condish.elev.pred, coy.bob.condish.for.ung, coy.bob.condish.for.pred, coy.bob.condish.wtd.ung, coy.bob.condish.wtd.pred, 
                                coy.bob.condish.lago.ung, coy.bob.condish.lago.pred, coy.bob.condish.div.ung, coy.bob.condish.div.pred)
  
  
  #'  -----------------------------
  ####  Save all the pretty plots  ####
  #'  -----------------------------
  #'  Marginal occupancy probabilities
  #'  Wolf-bear
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_ung_marginal_plots.tiff", wolf.bear.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_pred_marginal_plots.tiff", wolf.bear.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_ung_marginal_plots.tiff", wolf.bear.marg.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_pred_marginal_plots.tiff", wolf.bear.marg.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_div_ung_marginal_plots.tiff", wolf.bear.marg.plots[[5]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_div_pred_marginal_plots.tiff", wolf.bear.marg.plots[[6]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Wolf-coyote
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_ung_marginal_plots.tiff", wolf.coy.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_pred_marginal_plots.tiff", wolf.coy.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_ung_marginal_plots.tiff", wolf.coy.marg.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_pred_marginal_plots.tiff", wolf.coy.marg.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Coyote-bobcat
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_elev_ung_marginal_plots.tiff", coy.bob.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_elev_pred_marginal_plots.tiff", coy.bob.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_for_ung_marginal_plots.tiff", coy.bob.marg.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_for_pred_marginal_plots.tiff", coy.bob.marg.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Combined by covariate
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_elevation_plot.tiff", marginal_elev_plot, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_forestcov_plot.tiff", marginal_for_plot, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/marginal_preydiversity_plot.tiff", marginal_dif_plot, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #' #'  Conditional occupancy probabilities
  #' #'  Wolf-bear  
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_ung_conditional_plots.tiff", wolf.bear.condish.plots[[1]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_pred_conditional_plots.tiff", wolf.bear.condish.plots[[2]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_ung_conditional_plots.tiff", wolf.bear.condish.plots[[3]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_pred_conditional_plots.tiff", wolf.bear.condish.plots[[4]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' #'  Wolf-coyote 
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_ung_conditional_plots.tiff", wolf.coy.condish.plots[[1]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_pred_conditional_plots.tiff", wolf.coy.condish.plots[[2]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_ung_conditional_plots.tiff", wolf.coy.condish.plots[[3]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_pred_conditional_plots.tiff", wolf.coy.condish.plots[[4]], 
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Coyote-bobcat 
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_wtd_ung_conditional_plots.tiff", coy.bob.condish.plots[[5]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_wtd_pred_conditional_plots.tiff", coy.bob.condish.plots[[6]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_ung_conditional_plots.tiff", coy.bob.condish.plots[[7]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_lago_pred_conditional_plots.tiff", coy.bob.condish.plots[[8]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_div_ung_conditional_plots.tiff", coy.bob.condish.plots[[9]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_div_pred_conditional_plots.tiff", coy.bob.condish.plots[[10]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  ----------------------------------------
  ####  Predict mean occupancy for Yr1 & Yr2  ####
  #'  ----------------------------------------
  #'  Load null models (includes year effect on marginal occupancy)
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(yr)_p(.)_2023-04-10.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(yr)_p(.)_2023-04-10.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(yr)_p(.)_2023-04-06.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(yr)_p(.)_2023-04-06.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(yr)_p(.)_2023-04-06.RData")
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(yr)_p(.)_2023-04-10.RData")
  
  #'  Add binary Year variable to covariate data fram
  stations_skinny_eoe20s21s <- mutate(stations_skinny_eoe20s21s, 
                                      Year = ifelse(Season == "Smr20", 0, 1))
  
  #'  Predict mean occupancy for Yr1 and Yr2
  wolf.bear.mean.yr1 <- predict_occupancy(mod = wolf.bear.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 0), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.mean.yr2 <- predict_occupancy(mod = wolf.bear.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 1), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.coy.mean.yr1 <- predict_occupancy(mod = wolf.coy.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 0), psi_cov_index = 0,
                                          psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.coy.mean.yr2 <- predict_occupancy(mod = wolf.coy.null, ncat = 4, npoints = 500,
                                          focal_cov = stations_skinny_eoe20s21s$Year,
                                          psi_cov = c(1, 1), psi_cov_index = 0,
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
  coy.bob.mean.yr1 <- predict_occupancy(mod = coy.bob.null, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         psi_cov = c(1, 0), psi_cov_index = 0,
                                         psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  coy.bob.mean.yr2 <- predict_occupancy(mod = coy.bob.null, ncat = 4, npoints = 500,
                                         focal_cov = stations_skinny_eoe20s21s$Year,
                                         psi_cov = c(1, 1), psi_cov_index = 0,
                                         psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  
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
  
  #'  Bind all species and years together
  mean_annual_psi <- rbind(psimean_wolf.bear.yr1, psimean_wolf.bear.yr2, psimean_wolf.coy.yr1, 
                           psimean_wolf.coy.yr2, psimean_wolf.lion.yr1, psimean_wolf.lion.yr2, 
                           psimean_lion.bear.yr1, psimean_lion.bear.yr2, psimean_lion.bob.yr1,
                           psimean_lion.bob.yr2, psimean_coy.bob.yr1, psimean_coy.bob.yr2) %>%
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
  
  
  
  