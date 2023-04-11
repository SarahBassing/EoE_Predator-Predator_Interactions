  #'  ------------------------------
  #'  Plot co-occurrence results
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  February 2023
  #'  -------------------------------
  #'  Script to create figures of predicted co-occurrence and covariate effects
  #'  from multispecies occupancy models run in JAGS.
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
  load("./Outputs/MultiSpp_OccMod_Outputs/DIC_top_models.RData")
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
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_habitat_yr)_p(setup_effort)_2023-04-08.RData")
  wolf.bear.elev.ung.yr1a <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 0, 0, 0, 0), psi_cov_index = 4,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.elev.pred.yr1a <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$Elev,
                                             psi_cov = c(1, 1, 0, 0, 0), psi_cov_index = 4,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.for.ung.yr1a <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$PercForest,
                                             psi_cov = c(1, 0, 0, 0, 0), psi_cov_index = 5,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.for.pred.yr1a <- predict_occupancy(mod = wolf.bear.hab, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$PercForest,
                                             psi_cov = c(1, 1, 0, 0, 0), psi_cov_index = 5,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  
  #'  Second top model (within 2 deltaDIC)
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfbear_psi(setup_preydiversity_yr)_p(setup_effort)_2023-04-09.RData")
  wolf.bear.elev.ung.yr1b <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$Elev,
                                              psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 4,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.elev.pred.yr1b <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                               focal_cov = stations_skinny_eoe20s21s$Elev,
                                               psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 4,
                                               psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.for.ung.yr1b <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$PercForest,
                                             psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 5,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.for.pred.yr1b <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$PercForest,
                                              psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 5,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.div.ung.yr1b <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                             focal_cov = stations_skinny_eoe20s21s$SppDiversity,
                                             psi_cov = c(1, 0, 0, 0, 0, 0), psi_cov_index = 6,
                                             psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  wolf.bear.div.pred.yr1b <- predict_occupancy(mod = wolf.bear.preydiv, ncat = 4, npoints = 500,
                                              focal_cov = stations_skinny_eoe20s21s$SppDiversity,
                                              psi_cov = c(1, 1, 0, 0, 0, 0), psi_cov_index = 6,
                                              psi_inxs_cov = c(0), psi_inxs_cov_index = 0)
  
  #####  Wolf-Coyote predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolfcoy_psi(setup_habitat_yr)_p(setup_effort)_2023-04-04.RData") 
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
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/wolflion_psi(.)_p(.)_2023-04-05.RData")
  
  
  #####  Lion-Bear predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbear_psi(.)_p(.)_2023-03-31.RData")
  
  
  #####  Lion-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/lionbob_psi(.)_p(.)_2023-04-04.RData")
  
  
  #####  Coyote-Bobcat predictions  ####
  load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(global)_psix(global)_p(setup_effort)_2023-04-07.RData")
  # load("./Outputs/MultiSpp_OccMod_Outputs/JAGS_output/coybob_psi(setup_habitat_yr)_psix(.)_p(setup_effort)_2023-04-04.RData")
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
      geom_ribbon(aes(ymin = lowerCRI, ymax = upperCRI, fill = Species), alpha = 0.2) +
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
  wolf.bear.marg.elev.ung_a <- plot_marginal_occ(predicted = wolf.bear.elev.ung.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.elev.pred_a <- plot_marginal_occ(predicted = wolf.bear.elev.pred.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.for.ung_a <- plot_marginal_occ(predicted = wolf.bear.for.ung.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.for.pred_a <- plot_marginal_occ(predicted = wolf.bear.for.pred.yr1a, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "Trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.plots_topmod <- list(wolf.bear.marg.elev.ung_a, wolf.bear.marg.elev.pred_a, wolf.bear.marg.for.ung_a, wolf.bear.marg.for.pred_a)
  
  #'  Second top model
  wolf.bear.marg.elev.ung_b <- plot_marginal_occ(predicted = wolf.bear.elev.ung.yr1b, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.elev.pred_b <- plot_marginal_occ(predicted = wolf.bear.elev.pred.yr1b, spp1 = "Wolf", spp2 = "Bear", covname = "Elevation (m)", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.for.ung_b <- plot_marginal_occ(predicted = wolf.bear.for.ung.yr1b, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.for.pred_b <- plot_marginal_occ(predicted = wolf.bear.for.pred.yr1b, spp1 = "Wolf", spp2 = "Bear", covname = "Percent forest cover", setup = "Trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.div.ung_b <- plot_marginal_occ(predicted = wolf.bear.div.ung.yr1b, spp1 = "Wolf", spp2 = "Bear", covname = "Shannon's diversity index", setup = "random sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.div.pred_b <- plot_marginal_occ(predicted = wolf.bear.div.pred.yr1b, spp1 = "Wolf", spp2 = "Bear", covname = "Shannon's diversity index", setup = "trail sites", spppair = "Wolf - Black Bear")
  wolf.bear.marg.plots_2ndmod <- list(wolf.bear.marg.elev.ung_b, wolf.bear.marg.elev.pred_b, wolf.bear.marg.for.ung_b, wolf.bear.marg.for.pred_b, wolf.bear.marg.div.ung_b, wolf.bear.marg.div.pred_b)
  
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
  coy.bob.marg.wtd.ung <- plot_marginal_occ(predicted = coy.bob.wtd.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "White-tailed deer relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.wtd.pred <- plot_marginal_occ(predicted = coy.bob.wtd.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "White-tailed deer relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.lago.ung <- plot_marginal_occ(predicted = coy.bob.lago.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Lagomorph relative abundance", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.lago.pred <- plot_marginal_occ(predicted = coy.bob.lago.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Lagomorph relative abundance", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.div.ung <- plot_marginal_occ(predicted = coy.bob.div.ung.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Shannon's diversity index", setup = "random sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.div.pred <- plot_marginal_occ(predicted = coy.bob.div.pred.yr1, spp1 = "Coyote", spp2 = "Bobcat", covname = "Shannon's diversity index", setup = "trail sites", spppair = "Coyote - Bobcat")
  coy.bob.marg.plots <- list(coy.bob.marg.elev.ung, coy.bob.marg.elev.pred, coy.bob.marg.for.ung, coy.bob.marg.for.pred, 
                             coy.bob.marg.wtd.ung, coy.bob.marg.wtd.pred, coy.bob.marg.lago.ung, coy.bob.marg.lago.pred, coy.bob.marg.div.ung, coy.bob.marg.div.pred)
 
  #' #'  Function to plot marginal probabilities
  #' marginal_prob_plots <- function(focal_cov, marg, covariate_name, spp1, spp2, plotit = T) {
  #'   #'  Plot marginal probability for each species
  #'   matplot(focal_cov, marg[,,"mean"],
  #'           type = "l", lty = 1:2, lwd = 3, col = c("black", "darkblue"), frame = FALSE, xlab = covariate_name, #col = 1:2
  #'           ylab = "Marginal occupancy probability", ylim = c(0, 1), #ylim = c(0, max(marg[,,"upper"]))
  #'           main = paste("Marginal Occupancy for", spp1, "and", spp2))
  #'   #'  Add CRIs
  #'   polygon(c(focal_cov, rev(focal_cov)), c(marg[,1,"lower"], rev(marg[,1,"upper"])),
  #'           border = NA, col = adjustcolor("black", alpha = 0.2))
  #'   polygon(c(focal_cov, rev(focal_cov)), c(marg[,2,"lower"], rev(marg[,2,"upper"])),
  #'           border = NA, col = adjustcolor("darkblue", alpha = 0.2))
  #'   # for(i in 1:2) {
  #'   #   polygon(c(focal_cov, rev(focal_cov)), c(marg[,i,"lower"], rev(marg[,i,"upper"])),
  #'   #           border = NA, col = adjustcolor(i, alpha = 0.2))
  #'   # }
  #'   legend("topright", c(spp1, spp2), lwd = 3, lty = 1:2, col = c("black", "darkblue"), horiz = TRUE, bty = "n") #col = 1:2
  #' }
  #' #####  Wolf-Bear marginal occupancy  ####
  #' (wolf.bear.marg.elev.ung <- marginal_prob_plots(focal_cov = wolf.bear.elev.ung.yr1b[[3]], marg = wolf.bear.elev.ung.yr1b[[1]],
  #'                                                 covariate_name = "Elevation (m)", spp1 = "Wolf", spp2 = "Bear"))
  #' (wolf.bear.marg.elev.pred <- marginal_prob_plots(focal_cov = wolf.bear.elev.pred.yr1b[[3]], marg = wolf.bear.elev.pred.yr1b[[1]],
  #'                                                  covariate_name = "Elevation (m)", spp1 = "Wolf", spp2 = "Bear"))
  #' (wolf.bear.marg.for.ung <- marginal_prob_plots(focal_cov = wolf.bear.for.ung.yr1b[[3]], marg = wolf.bear.for.ung.yr1b[[1]],
  #'                                                covariate_name = "Percent forest cover", spp1 = "Wolf", spp2 = "Bear"))
  #' (wolf.bear.marg.for.pred <- marginal_prob_plots(focal_cov = wolf.bear.for.pred.yr1b[[3]], marg = wolf.bear.for.pred.yr1b[[1]],
  #'                                                 covariate_name = "Percent forest cover", spp1 = "Wolf", spp2 = "Bear"))
  #' (wolf.bear.marg.div.ung <- marginal_prob_plots(focal_cov = wolf.bear.div.ung.yr1b[[3]], marg = wolf.bear.div.ung.yr1b[[1]],
  #'                                                covariate_name = "Shannon's diversity", spp1 = "Wolf", spp2 = "Bear"))
  #' (wolf.bear.marg.div.pred <- marginal_prob_plots(focal_cov = wolf.bear.div.pred.yr1b[[3]], marg = wolf.bear.div.pred.yr1b[[1]],
  #'                                                 covariate_name = "Shannon's diversity", spp1 = "Wolf", spp2 = "Bear"))
  #' 
  #' 
  #' #####  Wolf-Coyote marginal occupancy  ####
  #' (wolf.coy.marg.elev.ung <- marginal_prob_plots(focal_cov = wolf.coy.elev.ung.yr1[[3]], marg = wolf.coy.elev.ung.yr1[[1]], 
  #'                                                covariate_name = "Elevation (m)", spp1 = "Wolf", spp2 = "Coyote"))
  #' (wolf.coy.marg.elev.pred <- marginal_prob_plots(focal_cov = wolf.coy.elev.pred.yr1[[3]], marg = wolf.coy.elev.pred.yr1[[1]], 
  #'                                                 covariate_name = "Elevation (m)", spp1 = "Wolf", spp2 = "Coyote"))
  #' (wolf.coy.marg.for.ung <- marginal_prob_plots(focal_cov = wolf.coy.for.ung.yr1[[3]], marg = wolf.coy.for.ung.yr1[[1]], 
  #'                                               covariate_name = "Percent forest cover", spp1 = "Wolf", spp2 = "Coyote"))
  #' (wolf.coy.marg.for.pred <- marginal_prob_plots(focal_cov = wolf.coy.for.pred.yr1[[3]], marg = wolf.coy.for.pred.yr1[[1]], 
  #'                                                covariate_name = "Percent forest cover", spp1 = "Wolf", spp2 = "Coyote"))
  #' 
  #' #####  Wolf-Lion marginal occupancy  ####
  #' #'  Nadda
  #' 
  #' #####  Lion-Bear marginal occupancy  ####
  #' #'  Nadda
  #' 
  #' #####  Lion-Bobcat marginal occupancy  ####
  #' #'  Nadda
  #' 
  #' #####  Coyote-Bobcat marginal occupancy  ####
  #' (coy.bob.marg.elev.ung <- marginal_prob_plots(focal_cov = coy.bob.elev.ung.yr1[[3]], marg = coy.bob.elev.ung.yr1[[1]],
  #'                                               covariate_name = "Elevation (m)", spp1 = "Coyote", spp2 = "Bobcat"))
  #' (coy.bob.marg.elev.pred <- marginal_prob_plots(focal_cov = coy.bob.elev.pred.yr1[[3]], marg = coy.bob.elev.pred.yr1[[1]],
  #'                                                covariate_name = "Elevation (m)", spp1 = "Coyote", spp2 = "Bobcat"))
  #' (coy.bob.marg.for.ung <- marginal_prob_plots(focal_cov = coy.bob.for.ung.yr1[[3]], marg = coy.bob.for.ung.yr1[[1]],
  #'                                              covariate_name = "Percent forest cover", spp1 = "Coyote", spp2 = "Bobcat"))
  #' (coy.bob.marg.for.pred <- marginal_prob_plots(focal_cov = coy.bob.for.pred.yr1[[3]], marg = coy.bob.for.pred.yr1[[1]],
  #'                                               covariate_name = "Percent forest cover", spp1 = "Coyote", spp2 = "Bobcat"))
  #' (coy.bob.marg.div.ung <- marginal_prob_plots(focal_cov = coy.bob.div.ung.yr1[[3]], marg = coy.bob.div.ung.yr1[[1]],
  #'                                               covariate_name = "Shannon's diversity index", spp1 = "Coyote", spp2 = "Bobcat"))
  #' (coy.bob.marg.div.pred <- marginal_prob_plots(focal_cov = coy.bob.div.pred.yr1[[3]], marg = coy.bob.div.pred.yr1[[1]],
  #'                                                covariate_name = "Shannon's diversity index", spp1 = "Coyote", spp2 = "Bobcat"))
  
  
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
      geom_ribbon(aes(ymin = lowerCRI, ymax = upperCRI, fill = Species_interaction), alpha = 0.2) +
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
  
  #' #'  Function to plot conditional probabilities
  #' conditional_prob_plot <- function(focal_cov, cond, covariate_name, spp1, spp2, plotit = T) { 
  #'   #'  Plot conditional probabilities of occupancy
  #'   op <- par(mfrow = c(1,2))
  #'   matplot(focal_cov, cond[,1:2,"mean"],
  #'           type = "l", lty = 1:2, lwd = 3, col = c("black", "darkblue"), frame = FALSE, #col = 1:2
  #'           xlab = covariate_name, ylab = "Conditional occupancy probability",
  #'           ylim = c(0, 1), main = spp1)
  #'   polygon(c(focal_cov, rev(focal_cov)), c(cond[,1,"lower"], rev(cond[,1,"upper"])),
  #'           border = NA, col = adjustcolor("black", alpha = 0.2))
  #'   polygon(c(focal_cov, rev(focal_cov)), c(cond[,2,"lower"], rev(cond[,2,"upper"])),
  #'           border = NA, col = adjustcolor("darkblue", alpha = 0.2))
  #'   # for(i in 1:2){
  #'   #   polygon(c(focal_cov, rev(focal_cov)), c(cond[,i,"lower"], rev(cond[,i,"upper"])),
  #'   #           border = NA, col = adjustcolor(i, alpha = 0.2))
  #'   # }
  #'   legend("topright", c("alone", paste("given", spp2)), lwd = 3, lty = 1:3, col = c("black", "darkblue"), bty = 'n', cex = 1.5) #col = 1:2
  #'   
  #'   matplot(focal_cov, cond[,3:4,"mean"],
  #'           type = "l", lty = 1:2, lwd = 3, col = c("black", "darkblue"), frame = FALSE, #col = 1:2
  #'           xlab = covariate_name, ylab = "Conditional occupancy probability",
  #'           ylim = c(0, 1), main = spp2)
  #'   polygon(c(focal_cov, rev(focal_cov)), c(cond[,3,"lower"], rev(cond[,3,"upper"])),
  #'           border = NA, col = adjustcolor("black", alpha = 0.2))
  #'   polygon(c(focal_cov, rev(focal_cov)), c(cond[,4,"lower"], rev(cond[,4,"upper"])),
  #'           border = NA, col = adjustcolor("darkblue", alpha = 0.2))
  #'   # for(i in 1:2){
  #'   #   polygon(c(focal_cov, rev(focal_cov)), c(cond[,i+2,"lower"], rev(cond[,i+2,"upper"])),
  #'   #           border = NA, col = adjustcolor(i, alpha = 0.2))
  #'   # }
  #'   legend("topright", c("alone", paste("given", spp1)), lwd = 3, lty = 1:3, col = c("black", "darkblue"), bty = 'n', cex = 1.5) #col = 1:2
  #'   par(op)
  #' }
  #' #####  Wolf-Bear conditional occupancy  ####
  #' (wolf.bear.condish.elev.ung <- conditional_prob_plot(focal_cov = wolf.bear.elev.ung.yr1[[3]], cond = wolf.bear.elev.ung.yr1[[2]], 
  #'                                                      covariate_name = "Elevation (m)", spp1 = "wolf", spp2 = "bear"))
  #' (wolf.bear.condish.elev.pred <- conditional_prob_plot(focal_cov = wolf.bear.elev.pred.yr1[[3]], cond = wolf.bear.elev.pred.yr1[[2]],
  #'                                                       covariate_name = "Elevation (m)", spp1 = "wolf", spp2 = "bear"))
  #' (wolf.bear.condish.for.ung <- conditional_prob_plot(focal_cov = wolf.bear.for.ung.yr1[[3]], cond = wolf.bear.for.ung.yr1[[2]],
  #'                                                     covariate_name = "Percent forest cover", spp1 = "wolf", spp2 = "bear"))
  #' (wolf.bear.condish.for.pred <- conditional_prob_plot(focal_cov = wolf.bear.for.pred.yr1[[3]], cond = wolf.bear.for.pred.yr1[[2]],
  #'                                                      covariate_name = "Percent forest cover", spp1 = "wolf", spp2 = "bear"))
  #' 
  #' 
  #' #####  Wolf-Coyote conditional occupancy  ####
  #' (wolf.coy.condish.elev.ung <- conditional_prob_plot(focal_cov = wolf.coy.elev.ung.yr1[[3]], cond = wolf.coy.elev.ung.yr1[[2]],
  #'                                                     covariate_name = "Elevation (m)", spp1 = "wolf", spp2 = "coyote"))
  #' (wolf.coy.condish.elev.pred <- conditional_prob_plot(focal_cov = wolf.coy.elev.pred.yr1[[3]], cond = wolf.coy.elev.pred.yr1[[2]],
  #'                                                      covariate_name = "Elevation (m)", spp1 = "wolf", spp2 = "coyote"))
  #' (wolf.coy.condish.for.ung <- conditional_prob_plot(focal_cov = wolf.coy.for.ung.yr1[[3]], cond = wolf.coy.for.ung.yr1[[2]],
  #'                                                    covariate_name = "Percent forest cover", spp1 = "wolf", spp2 = "coyote"))
  #' (wolf.coy.condish.for.pred <- conditional_prob_plot(focal_cov = wolf.coy.for.pred.yr1[[3]], cond = wolf.coy.for.pred.yr1[[2]],
  #'                                                     covariate_name = "Percent forest cover", spp1 = "wolf", spp2 = "coyote"))
  #'                       
  #' #####  Wolf-Lion conditional occupancy  ####
  #' #'  Nadda
  #' 
  #' #####  Lion-Bear conditional occupancy  ####
  #' #'  Nadda
  #' 
  #' #####  Lion-Bobcat conditional occupancy  ####
  #' #'  Nadda
  #' 
  #' #####  Coyote-Bobcat conditional occupancy  ####
  #' (coy.bob.condish.elev.ung <- conditional_prob_plot(focal_cov = coy.bob.elev.ung.yr1[[3]], cond = coy.bob.elev.ung.yr1[[2]],
  #'                                                    covariate_name = "Elevation (m)", spp1 = "coyote", spp2 = "bobcat"))
  #' (coy.bob.condish.elev.pred <- conditional_prob_plot(focal_cov = coy.bob.elev.pred.yr1[[3]], cond = coy.bob.elev.pred.yr1[[2]],
  #'                                                     covariate_name = "Elevation (m)", spp1 = "coyote", spp2 = "bobcat"))
  #' (coy.bob.condish.for.ung <- conditional_prob_plot(focal_cov = coy.bob.for.ung.yr1[[3]], cond = coy.bob.for.ung.yr1[[2]],
  #'                                                   covariate_name = "Percent forest cover", spp1 = "coyote", spp2 = "bobcat"))
  #' (coy.bob.condish.for.pred <- conditional_prob_plot(focal_cov = coy.bob.for.pred.yr1[[3]], cond = coy.bob.for.pred.yr1[[2]],
  #'                                                    covariate_name = "Percent forest cover", spp1 = "coyote", spp2 = "bobcat"))
  
  
  
  #'  Save
  #'  Wolf-bear marginal
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_ung_marginal_plots.tiff", wolf.bear.marg.plots_2ndmod[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_pred_marginal_plots.tiff", wolf.bear.marg.plots_2ndmod[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_ung_marginal_plots.tiff", wolf.bear.marg.plots_2ndmod[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_pred_marginal_plots.tiff", wolf.bear.marg.plots_2ndmod[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_div_ung_marginal_plots.tiff", wolf.bear.marg.plots_2ndmod[[5]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_div_pred_marginal_plots.tiff", wolf.bear.marg.plots_2ndmod[[6]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Wolf-coyote marginal
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_ung_marginal_plots.tiff", wolf.coy.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_pred_marginal_plots.tiff", wolf.coy.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_ung_marginal_plots.tiff", wolf.coy.marg.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_pred_marginal_plots.tiff", wolf.coy.marg.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Coyote-bobcat marginal
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_elev_ung_marginal_plots.tiff", coy.bob.marg.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_elev_pred_marginal_plots.tiff", coy.bob.marg.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_for_ung_marginal_plots.tiff", coy.bob.marg.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/coy-bob_for_pred_marginal_plots.tiff", coy.bob.marg.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #' #'  Wolf-bear conditional 
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_ung_conditional_plots.tiff", wolf.bear.condish.plots[[1]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_elev_pred_conditional_plots.tiff", wolf.bear.condish.plots[[2]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_ung_conditional_plots.tiff", wolf.bear.condish.plots[[3]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #' ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-bear_for_pred_conditional_plots.tiff", wolf.bear.condish.plots[[4]],
  #'        units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Wolf-coyote conditional
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_ung_conditional_plots.tiff", wolf.coy.condish.plots[[1]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_elev_pred_conditional_plots.tiff", wolf.coy.condish.plots[[2]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_ung_conditional_plots.tiff", wolf.coy.condish.plots[[3]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/MultiSpp_OccMod_Outputs/Co-Occ_Plots/wolf-coy_for_pred_conditional_plots.tiff", wolf.coy.condish.plots[[4]], 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  #'  Coyote-bobcat conditional
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
  
  
  