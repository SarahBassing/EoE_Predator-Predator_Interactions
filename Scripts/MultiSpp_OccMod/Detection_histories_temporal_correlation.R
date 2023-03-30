  #'  ---------------------------------------
  #'  Test annual detection correlations
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  March 2023
  #'  ---------------------------------------
  #'  Assess temporal correlation between detections at same site between years.
  #'  Calculate correlation coefficient, run Chi-squared test, and run simple 
  #'  logistic regressions with year 2 detections as response variable and year 1 
  #'  detections as explanatory variable. Significant p-values for Chi-squared
  #'  test indicate we reject the null (that data are independent). Significant
  #'  effect of year 1 detections indicate year 1 data are predictive of year 2.
  #'  ---------------------------------------
  
  #'  Load libraries
  library(abind)
  library(tidyverse)
  
  #'  Load covariate and detection history data
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe20s_predators.RData")
  load("./Data/MultiSpp_OccMod_Outputs/Detection_Histories/DH_eoe21s_predators.RData")
  
  #'  Pair detection histories
  #'  Combine annual detection histories
  all_detections <- function(dh1, dh2) {
    dh1 <- dh1[[1]]
    dh2 <- dh2[[1]]
      
    #'  Rename columns in each detection history
    newcols <- c("occ1", "occ2", "occ3", "occ4", "occ5", "occ6", "occ7", "occ8", "occ9", "occ10", "occ11")
    colnames(dh1) <- newcols; colnames(dh2) <- newcols
      
    dh <- list(dh1, dh2)
      
    return(dh)
  }
  DH_bear <- all_detections(DH_eoe20s_predators[[1]], DH_eoe21s_predators[[1]])
  DH_bob <- all_detections(DH_eoe20s_predators[[2]], DH_eoe21s_predators[[2]])
  DH_coy <- all_detections(DH_eoe20s_predators[[3]], DH_eoe21s_predators[[3]])
  DH_lion <- all_detections(DH_eoe20s_predators[[4]], DH_eoe21s_predators[[4]])
  DH_wolf <- all_detections(DH_eoe20s_predators[[5]], DH_eoe21s_predators[[5]])
  
  #'  Combine species detection histories into an array (site x survey x species)
  detection_array <- function(spp1, spp2, name1, name2) {
    #'  List detection histories
    spp12_DH <- list(spp1, spp2)
    #'  Name lists based on specific species pairing
    names(spp12_DH) <- c(name1, name2)
    #'  Format list into a 3D array
    spp12_array <- abind(spp12_DH, along = 3)
    
    return(spp12_array)
  }
  wolf_bear_DH_20s <- detection_array(spp1 = DH_wolf[[1]], spp2 = DH_bear[[1]], name1 = "wolf", name2 = "bear")
  wolf_bear_DH_21s <- detection_array(spp1 = DH_wolf[[2]], spp2 = DH_bear[[2]], name1 = "wolf", name2 = "bear")
  wolf_coy_DH_20s <- detection_array(spp1 = DH_wolf[[1]], spp2 = DH_coy[[1]], name1 = "wolf", name2 = "coyote")
  wolf_coy_DH_21s <- detection_array(spp1 = DH_wolf[[2]], spp2 = DH_coy[[2]], name1 = "wolf", name2 = "coyote")
  wolf_lion_DH_20s <- detection_array(spp1 = DH_wolf[[1]], spp2 = DH_lion[[1]], name1 = "wolf", name2 = "lion")
  wolf_lion_DH_21s <- detection_array(spp1 = DH_wolf[[2]], spp2 = DH_lion[[2]], name1 = "wolf", name2 = "lion")
  lion_bear_DH_20s <- detection_array(spp1 = DH_lion[[1]], spp2 = DH_bear[[1]], name1 = "lion", name2 = "bear")
  lion_bear_DH_21s <- detection_array(spp1 = DH_lion[[2]], spp2 = DH_bear[[2]], name1 = "lion", name2 = "bear")
  lion_bob_DH_20s <- detection_array(spp1 = DH_lion[[1]], spp2 = DH_bob[[1]], name1 = "lion", name2 = "bobcat")
  lion_bob_DH_21s <- detection_array(spp1 = DH_lion[[2]], spp2 = DH_bob[[2]], name1 = "lion", name2 = "bobcat")
  coy_bob_DH_20s <- detection_array(spp1 = DH_coy[[1]], spp2 = DH_bob[[1]], name1 = "coyote", name2 = "bobcat")
  coy_bob_DH_21s <- detection_array(spp1 = DH_coy[[2]], spp2 = DH_bob[[2]], name1 = "coyote", name2 = "bobcat")
  
  #'  List 2-species detection arrays together for faster formatting below
  DH_array_list_20s <- list(wolf_bear_DH_20s, wolf_coy_DH_20s, wolf_lion_DH_20s, lion_bear_DH_20s, lion_bob_DH_20s, coy_bob_DH_20s)
  DH_array_list_21s <- list(wolf_bear_DH_21s, wolf_coy_DH_21s, wolf_lion_DH_21s, lion_bear_DH_21s, lion_bob_DH_21s, coy_bob_DH_21s)
  
  #'  Snag site names
  sites_20s <- rownames(wolf_bear_DH_20s)
  sites_21s <- rownames(wolf_bear_DH_21s)
  
  #'  Species-specific naive occupancy
  z_naive <- function(bundled_dh, loc) {
    #'  Naive occupancy for each species at each site (site x spp matrix)
    zinit <- apply(bundled_dh, c(1, 3), sum, na.rm = TRUE)
    zinit[zinit > 1] <- 1
    zinit <- as.data.frame(zinit)
    zinit <- rownames_to_column(zinit)
    return(zinit)
  }
  z_20s <- lapply(DH_array_list_20s, z_naive, loc = sites_20s)
  z_21s <- lapply(DH_array_list_21s, z_naive, loc = sites_21s)
  zcat_names <- c("wolf_bear_zcat", "wolf_coy_zcat", "wolf_lion_zcat",
                  "lion_bear_zcat", "lion_bob_zcat", "coy_bob_zcat")
  names(z_20s) <- zcat_names
  names(z_21s) <- zcat_names
  
  #'  Reduce detection data to just sites sampled >1 year
  resampled_sites <- function(naiveocc20, naiveocc21) {
    #'  Retain only smr21 sites that were surveyed the previous year
    resampled_sites_21s <- naiveocc21[(naiveocc21$rowname %in% naiveocc20$rowname),]
    #'  Retain only smr20 sites that were surveyed the next year
    resampled_sites_20s <- naiveocc20[(naiveocc20$rowname %in% resampled_sites_21s$rowname),]
    #'  Double check they are they same length
    print(nrow(resampled_sites_21s)); print(nrow(resampled_sites_20s))
    #'  Merge into single df
    resample_dets <- full_join(resampled_sites_20s, resampled_sites_21s, by = "rowname")
    names(resample_dets) <- c("loc", "zcat.x", "zcat.y")
    return(resample_dets)
  } 
  resamp.sites.bear <- resampled_sites(z_20s[[1]][,c(1,3)], z_21s[[1]][,c(1,3)])
  resamp.sites.bob <- resampled_sites(z_20s[[5]][,c(1,3)], z_21s[[5]][,c(1,3)])
  resamp.sites.coy <- resampled_sites(z_20s[[2]][,c(1,3)], z_21s[[2]][,c(1,3)])
  resamp.sites.lion <- resampled_sites(z_20s[[4]][,c(1:2)], z_21s[[4]][,c(1:2)])
  resamp.sites.wolf <- resampled_sites(z_20s[[1]][,c(1:2)], z_21s[[1]][,c(1:2)])
  resamp.sites.singlespp <- list(resamp.sites.bear, resamp.sites.bob, resamp.sites.coy, 
                       resamp.sites.lion, resamp.sites.wolf)

  #'  Test whether single-species detections are correlated across years
  #'  Using Spearman's rank-sum correlation test b/c categorical variables
  spp_cor <- function(dets) {
    det_cor <- cor(dets$zcat.x, dets$zcat.y) 
    det_cor <- round(det_cor, 3)
    print(det_cor)
    return(det_cor)
  }
  spp_det_cor <- lapply(resamp.sites.singlespp, spp_cor)
  spp_det_cor <- as.data.frame(unlist(spp_det_cor))
  spp <- c("Bear", "Bobcat", "Coyote", "Lion", "Wolf")
  spp_det_cor <- cbind(spp, spp_det_cor)
  names(spp_det_cor) <- c("Predators", "Pearson's correlation") 
  
  #'  Save
  write.csv(spp_det_cor, file = "./Outputs/Tables/Annual_detection_correlation.csv")
  
  #'  Chi-squared test of independence
  #'  Don't need Yates's correction for continuity b/c sufficient sample size
  chisq_test <- function(dets) {
    print(table(dets[,c(2:3)]))
    out <- chisq.test(dets$zcat.y, dets$zcat.z, correct = FALSE)
    print(out)
    return(out)
  }
  spp_chisq_tst <- lapply(resamp.sites.singlespp, chisq_test)
  
  #'  Logistic regression - can previous year predict current year
  logistic_reg <- function(dets) {
    dat <- dets %>%
      mutate(zcat.x = as.factor(zcat.x))
    mod <- glm(zcat.y ~ zcat.x, data = dat, family = "binomial")
    print(summary(mod))
    return(mod)
  }
  regressed <- lapply(resamp.sites.singlespp, logistic_reg)
  
  