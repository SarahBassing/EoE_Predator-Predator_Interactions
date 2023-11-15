  #'  ----------------------------------------
  #'  Spatial autocorrelations
  #'  ID CRU - Predator Interactions
  #'  Sarah Bassing
  #'  November 2023
  #'  ----------------------------------------
  #'  Attempting semi-variograms and Moran's I to assess spatial autocorrelation
  #'  in TIFC density index
  #'  https://gsp.humboldt.edu/olm/R/04_01_Variograms.html
  #'  ----------------------------------------
  
  #'  Load libraries
  library(usdm)
  library(gstat)
  library(sf)
  library(tidyverse)
  
  #'  Load spatial data
  load("./Shapefiles/IDFG spatial data/Camera_locations/spatial_tifc_list.R")
  wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
  
  #'  Spatial autocorrelation per species, year, & GMU
  spatial_autocor <- function(dat, sill, mod, rng, spp, yr) {
    #'  Filter and calculate spatial autocorrelation per GMU
    tifc_gmu1 <- dat %>% filter(Gmu == "1") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_gmu1) <- wgs84
    sv_gmu1 <- variogram(cpue_100km2 ~ 1, tifc_gmu1)
    sv_gmu1$gmu <- "GMU1"
    plot(sv_gmu1)
    
    tifc_gmu6 <- dat %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_gmu6) <- wgs84
    sv_gmu6 <- variogram(cpue_100km2 ~ 1, tifc_gmu6)
    sv_gmu6$gmu <- "GMU6"
    plot(sv_gmu6)
    
    tifc_gmu10a <- dat %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_gmu10a) <- wgs84
    sv_gmu10a <- variogram(cpue_100km2 ~ 1, tifc_gmu10a)
    sv_gmu10a$gmu <- "GMU10a"
    plot(sv_gmu10a)
    
    #'  Merge GMU-specific data together
    sv_all_gmu <- rbind(sv_gmu1, sv_gmu6, sv_gmu10a)
    
    #'  Fit model
    vg_mod <- vgm(psill = sill, model = mod, nugget = 1, range = rng)
    vg_modfitted <- fit.variogram(sv_all_gmu, model = vg_mod)
    print(plot(sv_all_gmu, model = vg_modfitted, main = paste(spp, yr, "- GMU1 (black), GMU6 (green), GMU10A (pink)"), pch = 16, col = factor(sv_all_gmu$gmu)))
    
    return(sv_all_gmu)
  }
  bear_vg_21s <- spatial_autocor(spatial_tifc_list[[1]][[2]], sill = 20000, mod = "Exp", rng  = 35, spp = "black bear", yr = "2021")
  bear_vg_22s <- spatial_autocor(spatial_tifc_list[[1]][[3]], sill = 5000, mod = "Exp", rng  = 35, spp = "black bear", yr = "2022")
  bob_vg_21s <- spatial_autocor(spatial_tifc_list[[2]][[2]], sill = 300, mod = "Gau", rng  = 35, spp = "bobcat", yr = "2021")
  bob_vg_22s <- spatial_autocor(spatial_tifc_list[[2]][[3]], sill = 150, mod = "Exp", rng  = 35, spp = "bobcat", yr = "2022")
  coy_vg_21s <- spatial_autocor(spatial_tifc_list[[3]][[2]], sill = 6000, mod = "Exp", rng  = 35, spp = "coyote", yr = "2021")
  coy_vg_22s <- spatial_autocor(spatial_tifc_list[[3]][[3]], sill = 1000, mod = "Exp", rng  = 35, spp = "coyote", yr = "2022")
  lion_vg_21s <- spatial_autocor(spatial_tifc_list[[4]][[2]], sill = 100, mod = "Exp", rng  = 35, spp = "mountain lion", yr = "2021")
  lion_vg_22s <- spatial_autocor(spatial_tifc_list[[4]][[3]], sill = 100, mod = "Exp", rng  = 35, spp = "mountain lion", yr = "2022")
  wolf_vg_21s <- spatial_autocor(spatial_tifc_list[[5]][[2]], sill = 100, mod = "Gau", rng  = 10, spp = "wolf", yr = "2021")
  wolf_vg_22s <- spatial_autocor(spatial_tifc_list[[5]][[3]], sill = 1000, mod = "Exp", rng  = 35, spp = "wolf", yr = "2022")
  
  
  #'  Spatial autocorrelation per species in GMU10A
  spatial_autocor_gmu10a <- function(dat, sill, mod, rng, spp) {
    #'  Filter and calculate spatial autocorrelation per year
    tifc_yr1 <- dat[[1]] %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr1) <- wgs84
    sv_yr1 <- variogram(cpue_100km2 ~ 1, tifc_yr1)
    sv_yr1$year <- "2020"
    plot(sv_yr1)
    
    tifc_yr2 <- dat[[2]] %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr2) <- wgs84
    sv_yr2 <- variogram(cpue_100km2 ~ 1, tifc_yr2)
    sv_yr2$year <- "2021"
    plot(sv_yr2)
    
    tifc_yr3 <- dat[[3]] %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr3) <- wgs84
    sv_yr3 <- variogram(cpue_100km2 ~ 1, tifc_yr3)
    sv_yr3$year <- "2022"
    plot(sv_yr3)
    
    #'  Merge GMU-specific data together
    sv_all_yr <- rbind(sv_yr1, sv_yr2, sv_yr3)
    plot(sv_all_yr, col = factor(sv_all_yr$year))
    
    #'  Fit model
    vg_mod <- vgm(psill = sill, model = mod, nugget = 1, range = rng)
    vg_modfitted <- fit.variogram(sv_all_yr, model = vg_mod)
    print(plot(sv_all_yr, model = vg_modfitted, main = paste(spp, "GMU10A, all years"), pch = 16, col = factor(sv_all_yr$year)))
 
    return(sv_all_yr)
  }
  bear_gmu10a_vg <- spatial_autocor_gmu10a(spatial_tifc_list[[1]], sill = 5000, mod = "Gau", rng = 40, spp = "black bear")
  bob_gmu10a_vg <- spatial_autocor_gmu10a(spatial_tifc_list[[2]], sill = 40, mod = "Sph", rng = 15, spp = "bobcat")
  coy_gmu10a_vg <- spatial_autocor_gmu10a(spatial_tifc_list[[3]], sill = 5000, mod = "Gau", rng = 10, spp = "coyote")
  lion_gmu10a_vg <- spatial_autocor_gmu10a(spatial_tifc_list[[4]], sill = 50, mod = "Exp", rng = 30, spp = "mountain lion")
  wolf_gmu10a_vg <- spatial_autocor_gmu10a(spatial_tifc_list[[5]], sill = 200, mod = "Gau", rng = 20, spp = "wolf")
  
  
  #'  Spatial autocorrelation per species in GMU10A
  spatial_autocor_gmu6 <- function(dat, sill, mod, rng, spp) {
    #'  Filter and calculate spatial autocorrelation per year
    tifc_yr1 <- dat[[1]] %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr1) <- wgs84
    sv_yr1 <- variogram(cpue_100km2 ~ 1, tifc_yr1)
    sv_yr1$year <- "2020"
    plot(sv_yr1)
    
    tifc_yr2 <- dat[[2]] %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr2) <- wgs84
    sv_yr2 <- variogram(cpue_100km2 ~ 1, tifc_yr2)
    sv_yr2$year <- "2021"
    plot(sv_yr2)
    
    tifc_yr3 <- dat[[3]] %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr3) <- wgs84
    sv_yr3 <- variogram(cpue_100km2 ~ 1, tifc_yr3)
    sv_yr3$year <- "2022"
    plot(sv_yr3)
    
    #'  Merge GMU-specific data together
    sv_all_yr <- rbind(sv_yr1, sv_yr2, sv_yr3)
    plot(sv_all_yr, col = factor(sv_all_yr$year))
    
    #'  Fit model
    vg_mod <- vgm(psill = sill, model = mod, nugget = 1, range = rng)
    vg_modfitted <- fit.variogram(sv_all_yr, model = vg_mod)
    print(plot(sv_all_yr, model = vg_modfitted, main = paste(spp, "GMU6, all years"), pch = 16, col = "black"))
    
    return(sv_all_yr)
  }
  bear_gmu6_vg <- spatial_autocor_gmu6(spatial_tifc_list[[1]], sill = 12000, mod = "Gau", rng = 40, spp = "black bear")
  bob_gmu6_vg <- spatial_autocor_gmu6(spatial_tifc_list[[2]], sill = 150, mod = "Gau", rng = 40, spp = "bobcat")
  coy_gmu6_vg <- spatial_autocor_gmu6(spatial_tifc_list[[3]], sill = 7000, mod = "Gau", rng = 5, spp = "coyote")
  lion_gmu6_vg <- spatial_autocor_gmu6(spatial_tifc_list[[4]], sill = 100, mod = "Exp", rng = 25, spp = "mountain lion")
  wolf_gmu6_vg <- spatial_autocor_gmu6(spatial_tifc_list[[5]], sill = 100, mod = "Exp", rng = 30, spp = "wolf")
  
  
  #'  Spatial autocorrelation per species in GMU10A
  spatial_autocor_gmu1 <- function(dat, sill, mod, rng, spp) {
    #'  Filter and calculate spatial autocorrelation per year
    tifc_yr2 <- dat[[2]] %>% filter(Gmu == "1") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr2) <- wgs84
    sv_yr2 <- variogram(cpue_100km2 ~ 1, tifc_yr2)
    sv_yr2$year <- "2021"
    plot(sv_yr2)
    
    tifc_yr3 <- dat[[3]] %>% filter(Gmu == "1") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr3) <- wgs84
    sv_yr3 <- variogram(cpue_100km2 ~ 1, tifc_yr3)
    sv_yr3$year <- "2022"
    plot(sv_yr3)
    
    #'  Merge GMU-specific data together
    sv_all_yr <- rbind(sv_yr1, sv_yr2, sv_yr3)
    plot(sv_all_yr, col = factor(sv_all_yr$year))
    
    #'  Fit model
    vg_mod <- vgm(psill = sill, model = mod, nugget = 1, range = rng)
    vg_modfitted <- fit.variogram(sv_all_yr, model = vg_mod)
    print(plot(sv_all_yr, model = vg_modfitted, main = paste(spp, "GMU1, all years"), pch = 16, col = "black"))
    
    return(sv_all_yr)
  }
  bear_gmu1_vg <- spatial_autocor_gmu1(spatial_tifc_list[[1]], sill = 12000, mod = "Gau", rng = 40, spp = "black bear")
  bob_gmu1_vg <- spatial_autocor_gmu1(spatial_tifc_list[[2]], sill = 150, mod = "Gau", rng = 40, spp = "bobcat")
  coy_gmu1_vg <- spatial_autocor_gmu1(spatial_tifc_list[[3]], sill = 7000, mod = "Gau", rng = 5, spp = "coyote")
  lion_gmu1_vg <- spatial_autocor_gmu1(spatial_tifc_list[[4]], sill = 100, mod = "Exp", rng = 25, spp = "mountain lion")
  wolf_gmu1_vg <- spatial_autocor_gmu1(spatial_tifc_list[[5]], sill = 100, mod = "Exp", rng = 30, spp = "wolf")
  
  
  
  #'  Create semi-variogram for each species for 2022
  annual_tifc_variogram <- function(dat, sill.1, sill.2, sill.3, mod.1, mod.2, mod.3, rng.1, rng.2, rng.3, spp) {
    #'  2022 semi-variograms
    tifc_yr3_gmu1 <- dat[[3]] %>% filter(Gmu == "1") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr3_gmu1) <- wgs84
    vg_yr3_gmu1 <- variogram(cpue_100km2 ~ 1, tifc_yr3_gmu1)
    vg_yr3_gmu1_mod <- vgm(psill = sill.1, model = mod.1, nugget = 1, range = rng.1)
    vg_yr3_gmu1_modfitted <- fit.variogram(vg_yr3_gmu1, model = vg_yr3_gmu1_mod)
    print(plot(vg_yr3_gmu1, model = vg_yr3_gmu1_modfitted, main = paste(spp, "GMU1, 2022"), pch = 16, col = "black"))
    
    tifc_yr3_gmu6 <- dat[[3]] %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr3_gmu6) <- wgs84
    vg_yr3_gmu6 <- variogram(cpue_100km2 ~ 1, tifc_yr3_gmu6)
    vg_yr3_gmu6_mod <- vgm(psill = sill.2, model = mod.2, nugget = 1, range = rng.2)
    vg_yr3_gmu6_modfitted <- fit.variogram(vg_yr3_gmu6, model = vg_yr3_gmu6_mod)
    print(plot(vg_yr3_gmu6, model = vg_yr3_gmu6_modfitted, main = paste(spp, "GMU6, 2022"), pch = 16, col = "black"))

    tifc_yr3_gmu10a <- dat[[3]] %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2)
    st_crs(tifc_yr3_gmu10a) <- wgs84
    vg_yr3_gmu10a <- variogram(cpue_100km2 ~ 1, tifc_yr3_gmu10a)
    vg_yr3_gmu10a_mod <- vgm(psill = sill.3, model = mod.3, nugget = 1, range = rng.3)
    vg_yr3_gmu10a_modfitted <- fit.variogram(vg_yr3_gmu10a, model = vg_yr3_gmu10a_mod)
    print(plot(vg_yr3_gmu10a, model = vg_yr3_gmu10a_modfitted, main = paste(spp, "GMU10A, 2022"), pch = 16, col = "black"))

    #'  List all
    yr3_vg <- list(vg_yr3_gmu1, vg_yr3_gmu6, vg_yr3_gmu10a)
    return(yr3_vg)
    
  }
  vg_bear <- annual_tifc_variogram(spatial_tifc_list[[1]], 
                                   sill.1 = 1200, sill.2 = 1500, sill.3 = 4500, 
                                   mod.1 = "Gau", mod.2 = "Exp", mod.3 = "Gau",
                                   rng.1 = 5, rng.2 = 7, rng.3 = 1, spp = "Black bear")
  vg_bob <- annual_tifc_variogram(spatial_tifc_list[[2]], 
                                  sill.1 = 80, sill.2 = 175, sill.3 = 40, 
                                  mod.1 = "Sph", mod.2 = "Sph", mod.3 = "Sph", 
                                  rng.1 = 20, rng.2 = 30, rng.3 = 30, spp = "Bobcat")
  vg_coy <- annual_tifc_variogram(spatial_tifc_list[[3]], 
                                  sill.1 = 600, sill.2 = 1000, sill.3 = 4500, 
                                  mod.1 = "Gau", mod.2 = "Nug", mod.3 = "Nug", 
                                  rng.1 = 10, rng.2 = 0, rng.3 = 0, spp = "Coyote") # gmu1 converges
  vg_lion <- annual_tifc_variogram(spatial_tifc_list[[4]], 
                                   sill.1 = 150, sill.2 = 50, sill.3 = 10, 
                                   mod.1 = "Sph", mod.2 = "Gau", mod.3 = "Sph", 
                                   rng.1 = 20, rng.2 = 15, rng.3 = 15, spp = "Mountain lion")
  vg_wolf <- annual_tifc_variogram(spatial_tifc_list[[5]], 
                                   sill.1 = 900, sill.2 = 70, sill.3 = 100, 
                                   mod.1 = "Gau", mod.2 = "Exp", mod.3 = "Exp", 
                                   rng.1 = 30, rng.2 = 15, rng.3 = 25, spp = "Wolf")
  
  
  
  
  
  
  
  
  
  
  
  #'  Create semi-variogram for each species and year
  annual_tifc_variogram <- function(dat, spp) {
    
    #'  Year 1 semi-variograms
    tifc_yr1_gmu6 <- dat[[1]] %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr1_gmu6) <- wgs84
    vg_yr1_gmu6 <- variogram(cpue_100km2 ~ 1, tifc_yr1_gmu6)
    # print(plot(vg_yr1_gmu6, main = paste(spp, "GMU6, 2020")))
    vg_yr1_gmu6_mod <- vgm(psill = max(vg_yr1_gmu6$gamma), model = "Gau", nugget = 1, range = max(vg_yr1_gmu6$gamma))
    # print(plot(vg_yr1_gmu6, model = vg_yr1_gmu6_mod))
    vg_yr1_gmu6_modfitted <- fit.variogram(vg_yr1_gmu6, model = vg_yr1_gmu6_mod)
    print(plot(vg_yr1_gmu6, model = vg_yr1_gmu6_modfitted, main = paste(spp, "GMU6, 2020")))
    
    tifc_yr1_gmu10a <- dat[[1]] %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr1_gmu10a) <- wgs84
    vg_yr1_gmu10a <- variogram(cpue_100km2 ~ 1, tifc_yr1_gmu10a)
    # print(plot(vg_yr1_gmu10a, main = paste(spp, "GMU10A, 2020")))
    vg_yr1_gmu10a_mod <- vgm(psill = max(vg_yr1_gmu10a$gamma), model = "Gau", nugget = 1, range = max(vg_yr1_gmu10a$gamma))
    # print(plot(vg_yr1_gmu10a, model = vg_yr1_gmu10a_mod))
    vg_yr1_gmu10a_modfitted <- fit.variogram(vg_yr1_gmu10a, model = vg_yr1_gmu10a_mod)
    print(plot(vg_yr1_gmu10a, model = vg_yr1_gmu10a_modfitted, main = paste(spp, "GMU10A, 2020")))
    
    #'  Year 2 semi-variograms
    tifc_yr2_gmu1 <- dat[[2]] %>% filter(Gmu == "1") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr2_gmu1) <- wgs84
    vg_yr2_gmu1 <- variogram(cpue_100km2 ~ 1, tifc_yr2_gmu1)
    # print(plot(vg_yr2_gmu1, main = paste(spp, "GMU1, 2021")))
    vg_yr2_gmu1_mod <- vgm(psill = max(vg_yr2_gmu1$gamma), model = "Gau", nugget = 1, range = max(vg_yr2_gmu1$gamma))
    # print(plot(vg_yr2_gmu1, model = vg_yr2_gmu1_mod))
    vg_yr2_gmu1_modfitted <- fit.variogram(vg_yr2_gmu1, model = vg_yr2_gmu1_mod)
    print(plot(vg_yr2_gmu1, model = vg_yr2_gmu1_modfitted, main = paste(spp, "GMU1, 2021")))
    
    tifc_yr2_gmu6 <- dat[[2]] %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr2_gmu6) <- wgs84
    vg_yr2_gmu6 <- variogram(cpue_100km2 ~ 1, tifc_yr2_gmu6)
    # print(plot(vg_yr2_gmu6, main = paste(spp, "GMU6, 2021")))
    vg_yr2_gmu6_mod <- vgm(psill = max(vg_yr2_gmu6$gamma), model = "Gau", nugget = 1, range = max(vg_yr2_gmu6$gamma))
    # print(plot(vg_yr2_gmu6, model = vg_yr2_gmu6_mod))
    vg_yr2_gmu6_modfitted <- fit.variogram(vg_yr2_gmu6, model = vg_yr2_gmu6_mod)
    print(plot(vg_yr2_gmu6, model = vg_yr2_gmu6_modfitted, main = paste(spp, "GMU6, 2021")))
    
    tifc_yr2_gmu10a <- dat[[2]] %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr2_gmu10a) <- wgs84
    vg_yr2_gmu10a <- variogram(cpue_100km2 ~ 1, tifc_yr2_gmu10a)
    # print(plot(vg_yr2_gmu10a, main = paste(spp, "GMU10A, 2021")))
    vg_yr2_gmu10a_mod <- vgm(psill = max(vg_yr2_gmu10a$gamma), model = "Gau", nugget = 1, range = max(vg_yr2_gmu10a$gamma))
    # print(plot(vg_yr2_gmu10a, model = vg_yr2_gmu10a_mod))
    vg_yr2_gmu10a_modfitted <- fit.variogram(vg_yr2_gmu10a, model = vg_yr2_gmu10a_mod)
    print(plot(vg_yr2_gmu10a, model = vg_yr2_gmu10a_modfitted, main = paste(spp, "GMU10A, 2021")))
    
    #'  Year 3 semi-variograms
    tifc_yr3_gmu1 <- dat[[3]] %>% filter(Gmu == "1") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr3_gmu1) <- wgs84
    vg_yr3_gmu1 <- variogram(cpue_100km2 ~ 1, tifc_yr3_gmu1)
    # print(plot(vg_yr3_gmu1, main = paste(spp, "GMU1, 2022")))
    vg_yr3_gmu1_mod <- vgm(psill = max(vg_yr3_gmu1$gamma), model = "Gau", nugget = 1, range = max(vg_yr3_gmu1$gamma))
    # print(plot(vg_yr3_gmu1, model = vg_yr3_gmu1_mod))
    vg_yr3_gmu1_modfitted <- fit.variogram(vg_yr3_gmu1, model = vg_yr3_gmu1_mod)
    print(plot(vg_yr3_gmu1, model = vg_yr3_gmu1_modfitted, main = paste(spp, "GMU1, 2022")))
    
    tifc_yr3_gmu6 <- dat[[3]] %>% filter(Gmu == "6") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr3_gmu6) <- wgs84
    vg_yr3_gmu6 <- variogram(cpue_100km2 ~ 1, tifc_yr3_gmu6)
    # print(plot(vg_yr3_gmu6, main = paste(spp, "GMU6, 2022")))
    vg_yr3_gmu6_mod <- vgm(psill = max(vg_yr3_gmu6$gamma), model = "Gau", nugget = 1, range = max(vg_yr3_gmu6$gamma))
    # print(plot(vg_yr3_gmu6, model = vg_yr3_gmu6_mod))
    vg_yr3_gmu6_modfitted <- fit.variogram(vg_yr3_gmu6, model = vg_yr3_gmu6_mod)
    print(plot(vg_yr3_gmu6, model = vg_yr3_gmu6_modfitted, main = paste(spp, "GMU6, 2022")))
    
    tifc_yr3_gmu10a <- dat[[3]] %>% filter(Gmu == "10A") %>% dplyr::select(cpue_100km2) 
    st_crs(tifc_yr3_gmu10a) <- wgs84
    vg_yr3_gmu10a <- variogram(cpue_100km2 ~ 1, tifc_yr3_gmu10a)
    # print(plot(vg_yr3_gmu10a, main = paste(spp, "GMU10A, 2022")))
    vg_yr3_gmu10a_mod <- vgm(psill = max(vg_yr3_gmu10a$gamma), model = "Gau", nugget = 1, range = max(vg_yr3_gmu10a$gamma))
    # print(plot(vg_yr3_gmu10a, model = vg_yr3_gmu10a_mod))
    vg_yr3_gmu10a_modfitted <- fit.variogram(vg_yr3_gmu10a, model = vg_yr3_gmu10a_mod)
    print(plot(vg_yr3_gmu10a, model = vg_yr3_gmu10a_modfitted, main = paste(spp, "GMU10A, 2022")))
    
    #'  List all
    all_vg <- list(vg_yr1_gmu6, vg_yr1_gmu10a, vg_yr2_gmu1, vg_yr2_gmu6, vg_yr2_gmu10a, 
                   vg_yr3_gmu1, vg_yr3_gmu6, vg_yr3_gmu10a)
    return(all_vg)
 
  }
  vg_bear <- annual_tifc_variogram(spatial_tifc_list[[1]], spp = "Black bear")
  vg_bob <- annual_tifc_variogram(spatial_tifc_list[[2]], spp = "Bobcat")
  vg_coy <- annual_tifc_variogram(spatial_tifc_list[[3]], spp = "Coyote")
  vg_lion <- annual_tifc_variogram(spatial_tifc_list[[4]], spp = "Mountain lion")
  vg_wolf <- annual_tifc_variogram(spatial_tifc_list[[5]], spp = "Wolf")
  
  
  
  
  