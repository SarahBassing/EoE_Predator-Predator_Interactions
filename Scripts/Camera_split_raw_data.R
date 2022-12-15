  #'  ========================
  #'  Split raw data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ========================
  #'  Pull out, rename, and filter massive camera trap data sets into more 
  #'  manageable chunks and save for future analyses. Each data set is huge, 
  #'  especially the wolf ones. Likely need to load each data set and run script 
  #'  separately because data sets are so large.
  #'  ========================
  
  #'  Clean workspace
  rm(list = ls())  

  #'  Load libraries
  library(tidyverse)
  library(stringr)
  
  #'  =======================================================
  ####  Save deployment & detection data as separate files  ####
  #'  =======================================================
  #'  EoE 2020 summer
  load("./Data/IDFG camera data/eoe2020s_all_final.RData")
  cams_s20_eoe <- dep_dat_done; save(cams_s20_eoe, file = "./Data/IDFG camera data/eoe2020s_cameras.RData")
  dets_s20_eoe <- pic_dat_done; save(dets_s20_eoe, file = "./Data/IDFG camera data/eoe2020s_detections.RData")
  
  #'  EoE 2020 winter
  load("./Data/IDFG camera data/eoe2020w_pic_dep_loc.RData")
  cams_w20_eoe <- dep_dat; save(cams_w20_eoe, file = "./Data/IDFG camera data/eoe2020w_cameras.RData")
  dets_w20_eoe <- pic_dat; save(dets_w20_eoe, file = "./Data/IDFG camera data/eoe2020w_detections.RData")
  
  #'  EoE 2021 summer
  load("./Data/IDFG camera data/eoe2021s_dep_pic_loc.RData")
  cams_s21_eoe <- dep_dat; save(cams_s21_eoe, file = "./Data/IDFG camera data/eoe2021s_cameras.RData")
  dets_s21_eoe <- pic_dat; save(dets_s21_eoe, file = "./Data/IDFG camera data/eoe2021s_detections.RData")
  
  #'  Wolf 2019 summer/winter
  load("./Data/IDFG camera data/swwlf2019_dep_pic_loc.RData")
  cams_s19_wolf <- dep_dat; save(cams_s19_wolf, file = "./Data/IDFG camera data/wolf2019s_cameras.RData")
  dets_s19_wolf <- pic_dat; save(dets_s19_wolf, file = "./Data/IDFG camera data/wolf2019s_detections.RData")
  
  #'  Wolf 2020 summer/winter
  load("./Data/IDFG camera data/swwlf2020_dep_pic_loc.RData")
  cams_s20_wolf <- dep_dat; save(cams_s20_wolf, file = "./Data/IDFG camera data/wolf2020s_cameras.RData")
  dets_s20_wolf <- pic_dat; save(dets_s20_wolf, file = "./Data/IDFG camera data/wolf2020s_detections.RData")
  
  #'  Wolf 2021 summer/winter
  load("./Data/IDFG camera data/swwlf2021_dep_pic_loc.RData")
  cams_s21_wolf <- dep_dat; save(cams_s21_wolf, file = "./Data/IDFG camera data/wolf2021s_cameras.RData")
  dets_s21_wolf <- pic_dat; save(dets_s21_wolf, file = "./Data/IDFG camera data/wolf2021s_detections.RData")
  
  
  #'  ======================================================
  #####  Filter to more manageable data sets to work with  ####
  #'  ======================================================
  #'  Detection data
  load("./Data/IDFG camera data/eoe2020s_detections.RData")
  load("./Data/IDFG camera data/eoe2020w_detections.RData")
  load("./Data/IDFG camera data/eoe2021s_detections.RData")
  
  load("./Data/IDFG camera data/wolf2019s_detections.RData")
  load("./Data/IDFG camera data/wolf2020s_detections.RData")
  load("./Data/IDFG camera data/wolf2021s_detections.RData")
  
  #'  Filter HUGE data set into smaller chunks 
  filter_dets <- function(dets, prefixname) {
    #'  Make sure GMUs with "A"s are always uppercase
    dets <- mutate(dets, Gmu = toupper(Gmu))
    
    #'  Motion triggered detections
    allM <- dets[dets$TriggerMode == "M",]
    wildM <- allM[allM$Wildlife != FALSE,]
    humanM <- allM[allM$Human != FALSE,]
    vehicleM <- allM[allM$Vehicle != FALSE,]
    livestockM <- allM[allM$Livestock != FALSE,]
    domesticM <- allM[allM$Pet_pack_horse != FALSE,]
    
    #'  Time triggered detections
    allT <- dets[dets$TriggerMode == "T",]
    wildT <- allT[allT$Wildlife != FALSE,]
    humanT <- allT[allT$Human != FALSE,]
    vehicleT <- allT[allT$Vehicle != FALSE,]
    livestockT <- allT[allT$Livestock != FALSE,]
    domesticT <- allT[allT$Pet_pack_horse != FALSE,]
    
    #'  List smaller data sets together
    dets_list <- list(allM, wildM, humanM, vehicleM, livestockM, domesticM, 
                      allT, wildT, humanT, vehicleT, livestockT, domesticT)
    names(dets_list) <- c(paste0(prefixname,"allM"), paste0(prefixname,"wildM"), paste0(prefixname,"humanM"), paste0(prefixname,"vehicleM"), paste0(prefixname,"livestockM"), paste0(prefixname,"domesticM"), 
                          paste0(prefixname,"allT"), paste0(prefixname,"wildT"), paste0(prefixname,"humanT"), paste0(prefixname,"vehicleT"), paste0(prefixname,"livestockT"), paste0(prefixname,"domesticT"))
    return(dets_list)
  }
  eoe20s_dets <- filter_dets(dets_s20_eoe, prefixname = "eoe20s_")
  eoe20w_dets <- filter_dets(dets_w20_eoe, prefixname = "eoe20w_")
  eoe21s_dets <- filter_dets(dets_s21_eoe, prefixname = "eoe21s_")
  
  wolf19s_dets <- filter_dets(dets_s19_wolf, prefixname = "wolf19s_")
  wolf20s_dets <- filter_dets(dets_s20_wolf, prefixname = "wolf20s_")
  wolf21s_dets <- filter_dets(dets_s21_wolf, prefixname = "wolf21s_")
  
  
  #'  Save smaller data sets as individual RData files
  #'  https://stackoverflow.com/questions/30848091/saving-elements-of-a-list-as-data-frames-using-r
  save_Rdata <- function(dat) {
    invisible(lapply(names(dat), function(u) {
      assign(u, dat[[u]])
      save(list = u, file = paste0("./Data/IDFG camera data/Split datasets/", u, ".RData"))
    }))
  }
  save_Rdata(eoe20s_dets)
  save_Rdata(eoe20w_dets)
  save_Rdata(eoe21s_dets)
  
  save_Rdata(wolf19s_dets)
  save_Rdata(wolf20s_dets)
  save_Rdata(wolf21s_dets)
  
  
  #'  =====================================================================
  ####  Reduce data sets even further to just images flagged as "Keepers"  ####
  #'  =====================================================================
  #'  Keepers: flagged by S. Thompson as images with animals, some maintenance 
  #'  images, all images taken at noon, and a smattering of setup-retrieval day 
  #'  pics to help in case of ID errors.
  
  #'  Motion triggered images
  load("./Data/IDFG camera data/Split datasets/eoe20s_allM.RData")
  load("./Data/IDFG camera data/Split datasets/eoe20w_allM.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allM.RData")
  
  load("./Data/IDFG camera data/Split datasets/wolf19s_allM.RData")
  load("./Data/IDFG camera data/Split datasets/wolf20s_allM.RData")
  load("./Data/IDFG camera data/Split datasets/wolf21s_allM.RData")
  
  #'  Time triggered images
  load("./Data/IDFG camera data/Split datasets/eoe20s_allT.RData")
  load("./Data/IDFG camera data/Split datasets/eoe20w_allT.RData")
  load("./Data/IDFG camera data/Split datasets/eoe21s_allT.RData")
  
  load("./Data/IDFG camera data/Split datasets/wolf19s_allT.RData")
  load("./Data/IDFG camera data/Split datasets/wolf20s_allT.RData")
  load("./Data/IDFG camera data/Split datasets/wolf21s_allT.RData")
  
  
  #'  Filter down to "Keepers"
  keepers <- function(dets) {
    keep_dets <- filter(dets, SetLocation == "KeepSet")
    return(keep_dets)
  }
  eoe_motion <- list(eoe20s_allM, eoe20w_allM, eoe21s_allM)
  eoe_motion_skinny <- lapply(eoe_motion, keepers)
  save(eoe_motion_skinny, file = "./Data/IDFG camera data/Split datasets/eoe_motion_skinny.RData")
  
  eoe_time <- list(eoe20s_allT, eoe20w_allT, eoe21s_allT)
  eoe_time_skinny <- lapply(eoe_time, keepers)
  save(eoe_time_skinny, file = "./Data/IDFG camera data/Split datasets/eoe_time_skinny.RData")
  
  wolf_motion <- list(wolf19s_allM, wolf20s_allM, wolf21s_allM)
  wolf_motion_skinny <- lapply(wolf_motion, keepers)
  save(wolf_motion_skinny, file = "./Data/IDFG camera data/Split datasets/wolf_motion_skinny.RData")
  
  wolf_time <- list(wolf19s_allT, wolf20s_allT, wolf21s_allT)
  wolf_time_skinny <- lapply(wolf_time, keepers)
  save(wolf_time_skinny, file = "./Data/IDFG camera data/Split datasets/wolf_time_skinny.RData")
  
  