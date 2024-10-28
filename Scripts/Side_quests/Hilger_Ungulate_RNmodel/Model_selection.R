  #'  -------------------------------------
  #'  Model selection
  #'  Northern Idaho Predator-Prey Project
  #'  Sarah B. Bassing
  #'  October 2024
  #'  -------------------------------------
  #'  Using model selection to identify best supported model for each species and month
  #'  -------------------------------------
  
  #'  Clear workspace
  rm(list = ls())

  #'  Load libraries
  library(AICcmodavg)
  library(tidyverse)
  
  #'  Load model outputs into a temporary working environment and list together
  elk_july_files <- list.files("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Elk_july_mods", pattern="*.RData", full.names=TRUE)
  elk_july_mods <- lapply(elk_july_files, load, temp_env <- new.env())  #.GlobalEnv
  elk_july_mods <- as.list(temp_env)
  
  elk_aug_files <- list.files("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/Elk_aug_mods", pattern="*.RData", full.names=TRUE)
  elk_aug_mods <- lapply(elk_aug_files, load, temp_env <- new.env())  
  elk_aug_mods <- as.list(temp_env)
  
  wtd_july_files <- list.files("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/WTD_july_mods", pattern="*.RData", full.names=TRUE)
  wtd_july_mods <- lapply(wtd_july_files, load, temp_env <- new.env())  
  wtd_july_mods <- as.list(temp_env)
  
  wtd_aug_files <- list.files("./Outputs/Hilger_RNmodel/JAGS_out/Fit_10.26.24/WTD_aug_mods", pattern="*.RData", full.names=TRUE)
  wtd_aug_mods <- lapply(wtd_aug_files, load, temp_env <- new.env())  
  wtd_aug_mods <- as.list(temp_env)
  
  #'  Function to calculate WAIC
  calc.waic <- function(x) {
    #'  log-likelihood for every iteration and camera
    like <- as.matrix(x$sims.list$log_N)
    #'  mean likelihood 
    fbar <- colMeans(exp(like)) 
    #'  mean variance in log-likelihood 
    Pw <- sum(apply(like,2,var)) 
    #'  WAIC
    WAIC<- -2*sum(log(fbar))+2*Pw
    
    return(WAIC)
  }
  elk_july_WAIC <- lapply(elk_july_mods, calc.waic) %>% do.call(rbind.data.frame, .) 
  elk_aug_WAIC <- lapply(elk_aug_mods, calc.waic) %>% do.call(rbind.data.frame, .) 
  wtd_july_WAIC <- lapply(wtd_july_mods, calc.waic) %>% do.call(rbind.data.frame, .) 
  wtd_aug_WAIC <- lapply(wtd_aug_mods, calc.waic) %>% do.call(rbind.data.frame, .) 
  
  #'  Function to calculate WAIC joint-likelihood approach (WAICj) based on Gaya & Ketz (2024)
  calc.waicj <- function(x) {
    #'  log-likelihood for every iteration and camera
    like_joint <- as.matrix(x$sims.list$loglike.new) 
    #'  mean likelihood 
    fbar_joint <- colMeans(exp(like_joint))
    #'  mean variance in log-likelihood 
    Pw_joint <- sum(apply(like_joint,2,var)) 
    #'  joint likelihood WAIC
    WAIC_joint<- -2*sum(log(fbar_joint))+2*Pw_joint
    
    return(WAIC_joint)
  }
  elk_july_WAICj <- lapply(elk_july_mods, calc.waicj) %>% do.call(rbind.data.frame, .) 
  elk_aug_WAICj <- lapply(elk_aug_mods, calc.waicj) %>% do.call(rbind.data.frame, .) 
  wtd_july_WAICj <- lapply(wtd_july_mods, calc.waicj) %>% do.call(rbind.data.frame, .) 
  wtd_aug_WAICj <- lapply(wtd_aug_mods, calc.waicj) %>% do.call(rbind.data.frame, .) 
  
  #'  Function to reformat and create model selection tables for WAIC and the joint-likelihood WAIC approach
  waic_modSelect <- function(x, mods, select_method) {
    names(x) <- "WAIC"
    modSel_tbl <- x %>%
      bind_cols(names(mods)) %>%
      mutate(deltaWAIC = min(WAIC) - WAIC,
             deltaWAIC = abs(deltaWAIC)) %>%
      relocate(WAIC, .before = deltaWAIC) %>%
      arrange(WAIC)
    names(modSel_tbl) <- c("Model Name", select_method, paste0("delta", select_method))
    return(modSel_tbl)
  }
  (modSel_elkJuly_WAIC <- waic_modSelect(elk_july_WAIC, elk_july_mods, "WAIC"))
  (modSel_elkJuly_WAICj <- waic_modSelect(elk_july_WAICj, elk_july_mods, "WAICj"))
  (modSel_elkAug_WAIC <- waic_modSelect(elk_aug_WAIC, elk_aug_mods, "WAIC"))
  (modSel_elkAug_WAICj <- waic_modSelect(elk_aug_WAICj, elk_aug_mods, "WAICj"))
  (modSel_wtdJuly_WAIC <- waic_modSelect(wtd_july_WAIC, wtd_july_mods, "WAIC"))
  (modSel_wtdJuly_WAICj <- waic_modSelect(wtd_july_WAICj, wtd_july_mods, "WAICj"))
  (modSel_wtdAug_WAIC <- waic_modSelect(wtd_aug_WAIC, wtd_aug_mods, "WAIC"))
  (modSel_wtdAug_WAICj <- waic_modSelect(wtd_aug_WAICj, wtd_aug_mods, "WAICj"))
  
  #'  Model selection using DIC
  (modSel_elkJuly_DIC <- dictab(cand.set = elk_july_mods, modnames = names(elk_july_mods), sort = TRUE))
  (modSel_elkAug_DIC <- dictab(cand.set = elk_aug_mods, modnames = names(elk_aug_mods), sort = TRUE))
  (modSel_wtdJuly_DIC <- dictab(cand.set = wtd_july_mods, modnames = names(wtd_july_mods), sort = TRUE))
  (modSel_wtdAug_DIC <- dictab(cand.set = wtd_aug_mods, modnames = names(wtd_aug_mods), sort = TRUE))
  
  
  