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
    wildlifeT <- Tt[Tt$Wildlife != FALSE,]
    humanT <- Tt[Tt$Human != FALSE,]
    vehicleT <- Tt[Tt$Vehicle != FALSE,]
    livestockT <- Tt[Tt$Livestock != FALSE,]
    horseT <- Tt[Tt$Pet_pack_horse != FALSE,]
    
    #'  List smaller data sets together
    dets_list <- list(Mt, wildM, humanM, vehicleM, livestockM, domesticM, 
                      Tt, wildlifeT, humanT, vehicleT, livestockT, horseT)
    return(dets_list)
  }
  eoe20s_dets <- filter_dets(dets_s20)
  
  #' #'  Motion triggered detections
  #' eoe20s_M <- dets_s20[dets_s20$TriggerMode == "M",]
  #' eoe20s_wildM <- eoe20s_M[eoe20s_M$Wildlife != FALSE,]
  #' eoe20s_humanM <- eoe20s_M[eoe20s_M$Human != FALSE,]
  #' eoe20s_vehicleM <- eoe20s_M[eoe20s_M$Vehicle != FALSE,]
  #' eoe20s_livestockM <- eoe20s_M[eoe20s_M$Livestock != FALSE,]
  #' eoe20s_domesticM <- eoe20s_M[eoe20s_M$Pet_pack_horse != FALSE,]
  #' 
  #' #'  Time triggered detections
  #' eoe20s_T <- dets_s20[dets_s20$TriggerMode == "T",]
  #' eoe20s_wildlifeT <- eoe20s_T[eoe20s_T$Wildlife != FALSE,]
  #' eoe20s_humanT <- eoe20s_T[eoe20s_T$Human != FALSE,]
  #' eoe20s_vehicleT <- eoe20s_T[eoe20s_T$Vehicle != FALSE,]
  #' eoe20s_livestockT <- eoe20s_T[eoe20s_T$Livestock != FALSE,]
  #' eoe20s_horseT <- eoe20s_T[eoe20s_T$Pet_pack_horse != FALSE,]
  
  
  #'  Save smaller data sets as individual RData files
  #'  This takes a LONG time to do
  save_RData <- function(dat){
    save(dat, file = paste("./Data/IDFG camera data/Split datasets/", dat, ".RData"))
  }
  dat_list <- list(eoe20s_dets[[12]])
  lapply(dat_list, save_RData)
  # tst <- lapply(dat_list, save_RData)
  
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
  
  
  
  
  
  
  