  #'  ========================
  #'  Split raw data
  #'  ICFWRU PredXPred Project
  #'  Sarah B. Bassing
  #'  September 2022
  #'  ========================
  #'  Pull out, rename, and filter massive camera trap data sets into more 
  #'  manageable chunks and save for future analyses. Each data set is huge, 
  #'  especially the wolf ones. Likely need to load and run script separately
  #'  for each original file.
  #'  ========================
  
  #'  Clean workspace
  rm(list = ls())  

  #'  Load libraries
  library(tidyverse)
  library(stringr)
  
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
  
  
  #'  Filter and save even further - need more manageable data sets to work with
  #'  Detection data
  load("./Data/IDFG camera data/eoe2020s_detections.RData")
  load("./Data/IDFG camera data/eoe2020w_detections.RData")
  load("./Data/IDFG camera data/eoe2021s_detections.RData")
  
  load("./Data/IDFG camera data/wolf2019s_detections.RData")
  load("./Data/IDFG camera data/wolf2020s_detections.RData")
  load("./Data/IDFG camera data/wolf2021s_detections.RData")
  
  #'  Filter HUGE data set into smaller chunks 
  filter_dets <- function(dets, prefixname) {
    #'  Motion triggered detections
    allM <- dets[dets$TriggerMode == "M",]
    wildM <- Mt[Mt$Wildlife != FALSE,]
    humanM <- Mt[Mt$Human != FALSE,]
    vehicleM <- Mt[Mt$Vehicle != FALSE,]
    livestockM <- Mt[Mt$Livestock != FALSE,]
    domesticM <- Mt[Mt$Pet_pack_horse != FALSE,]
    
    #'  Time triggered detections
    allT <- dets[dets$TriggerMode == "T",]
    wildT <- Tt[Tt$Wildlife != FALSE,]
    humanT <- Tt[Tt$Human != FALSE,]
    vehicleT <- Tt[Tt$Vehicle != FALSE,]
    livestockT <- Tt[Tt$Livestock != FALSE,]
    domesticT <- Tt[Tt$Pet_pack_horse != FALSE,]
    
    #'  List smaller data sets together
    dets_list <- list(Mt, wildM, humanM, vehicleM, livestockM, domesticM, 
                      Tt, wildT, humanT, vehicleT, livestockT, domesticT)
    names(dets_list) <- c(paste0(prefixname,"Mt"), paste0(prefixname,"wildM"), paste0(prefixname,"humanM"), paste0(prefixname,"vehicleM"), paste0(prefixname,"livestockM"), paste0(prefixname,"domesticM"), 
                          paste0(prefixname,"Tt"), paste0(prefixname,"wildT"), paste0(prefixname,"humanT"), paste0(prefixname,"vehicleT"), paste0(prefixname,"livestockT"), paste0(prefixname,"domesticT"))
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
  
  