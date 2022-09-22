  #'  ---------------------------------
  #'  Cleaning raw detection data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ---------------------------------
  #'  Script to clean camea deteciton data, make sure annual datasets are consistent, 
  #'  filter into manageable data sets, and connect to camera locations.
  #'  ---------------------------------
  
  #'  Load libraries
  library(stringr)
  library(sf)
  library(raster)
  library(ggplot2)
  library(tidyverse)
  
  #'  Camera locations
  cams_eoe <- read.csv("./Data/IDFG camera data/cams_eoe_skinny.csv") %>%
    dplyr::select(-X)
  # cams_wolf <- read.csv("./Data/IDFG camera data/cams_wolf_skinny.csv") %>%
  #   dplyr::select(-X)
  
  #'  Detection data
  load("./Data/IDFG camera data/eoe2020s_detections.RData")
  # load("./Data/IDFG camera data/eoe2020w_detections.RData")
  # load("./Data/IDFG camera data/eoe2021s_detections.RData")
  
  # load("./Data/IDFG camera data/wolf2019s_detections.RData")
  # load("./Data/IDFG camera data/wolf2020s_detections.RData")
  # load("./Data/IDFG camera data/wolf2021s_detections.RData")
  
  #'  Filter HUGE data set into smaller chunks so data are more manageable
  filter_dets <- function(dets) {
    #'  Motion triggered detections
    Mt <- dets_s20[dets_s20$TriggerMode == "M",]
    wildM <- Mt[Mt$Wildlife != FALSE,]
    humanM <- Mt[Mt$Human != FALSE,]
    vehicleM <- Mt[Mt$Vehicle != FALSE,]
    livestockM <- Mt[Mt$Livestock != FALSE,]
    domesticM <- Mt[Mt$Pet_pack_horse != FALSE,]
    
    #'  Time triggered detections
    Tt <- dets_s20[dets_s20$TriggerMode == "T",]
    wildT <- Tt[Tt$Wildlife != FALSE,]
    humanT <- Tt[Tt$Human != FALSE,]
    vehicleT <- Tt[Tt$Vehicle != FALSE,]
    livestockT <- Tt[Tt$Livestock != FALSE,]
    domesticT <- Tt[Tt$Pet_pack_horse != FALSE,]
    
    #'  List smaller data sets together
    dets_list <- list(Mt, wildM, humanM, vehicleM, livestockM, domesticM, 
                      Tt, wildT, humanT, vehicleT, livestockT, domesticT)
    names(dets_list) <- c("Mt", "wildM", "humanM", "vehicleM", "livestockM", "domesticM", 
                          "Tt", "wildT", "humanT", "vehicleT", "livestockT", "domesticT")
    return(dets_list)
  }
  eoe20s_dets <- filter_dets(dets_s20)
  
  
  #'  Save smaller data sets as individual RData files
  save_Rdata <- function(dat, lname) {
    # assign(paste(lname), dat)
    save(dat, file = paste0("./Data/IDFG camera data/Split datasets/", lname, '.RData'))
  }
  dat_list <- list(eoe20s_dets[[11]], eoe20s_dets[[12]])
  names(dat_list) <- c("tst1", "tst2")
  lnames <- list("tst1", "tst2")
  tst <- mapply(save_Rdata, dat_list, lnames)
  
  hkguyfy <- dat_list[[1]]
  PLEASE <- dat_list
  saveit <- function(...) {
    x <- list(...)
    save(list=names(x), file = paste0("./Data/IDFG camera data/Split datasets/workdammit.RData"))
  }
  saveit(hkguyfy)
  saveit(dat_list)
  load("./Data/IDFG camera data/Split datasets/workdammit.RData")
  lapply(dat_list, saveit)
  foo <- 1
  saveit(bar=foo, file="hi.Rdata")
  
  
  
  library(parallel)
  #https://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html
  detectCores()
  start <- proc.time()
  cl <- makeCluster(3)
  clusterExport(cl, "eoe20s_dets")
  clusterEvalQ(cl, eoe20s_dets)
  save_dat <- parLapply(cl, eoe20s_dets, save_RData)
  stopCluster(cl)
  end <- proc.time()
  print(end - start)
  
  
  
  
  
  
  