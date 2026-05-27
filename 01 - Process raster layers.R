# PROJECT: Habitat selection
# SCRIPT: 01 - Process raster layers
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 20 Apr 2026
# COMPLETED: 21 Apr 2026
# LAST MODIFIED: 27 May 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(terra)
library(landscapemetrics)

# ______________________________________________________________________________
# 2. Read rasters ----
# ______________________________________________________________________________

# raster directory
dir.rast <- "D:/hare_project/data_spatial/Rasters/"

# cover type
rast.cover.pre <- rast(paste0(dir.rast, "cover_type/cover_type_pre.tif"))
rast.cover.post <- rast(paste0(dir.rast, "cover_type/cover_type_post.tif"))

# distance 
rast.dEdge <- rast(paste0(dir.rast, "dist_rasters/dEdge.tif"))
rast.dOpen <- rast(paste0(dir.rast, "dist_rasters/dOpen.tif"))
rast.dDM <- rast(paste0(dir.rast, "dist_rasters/dDM.tif"))

# vegetation models
rast.stem.pre <- rast(paste0(dir.rast, "veg_pred/RF/pre_stem.tif"))
rast.stem.post <- rast(paste0(dir.rast, "veg_pred/RF/post_stem.tif"))
rast.vo.pre <- rast(paste0(dir.rast, "veg_pred/RF/pre_vo.tif"))
rast.vo.post <- rast(paste0(dir.rast, "veg_pred/RF/post_vo.tif"))
#rast.shrub <- rast(paste0(dir.rast, "veg_pred/RF/sr.tif"))

# topography
rast.twi <-  rast(paste0(dir.rast, "Topography/twi_10.tif"))
rast.vrm <-  rast(paste0(dir.rast, "Topography/vrmL_10.tif"))

# PCT
#rast.dPil <-  rast(paste0(dir.rast, "PCT/dPiles.tif"))
#rast.dRet <-  rast(paste0(dir.rast, "PCT/dRet.tif"))
#rast.dUnitInt <- rast(paste0(dir.rast, "PCT/dUnitInt.tif"))

# ______________________________________________________________________________
# 3. Mosaic canopy layers ----

# these will be the 10-m resolution ones

# ______________________________________________________________________________

# 2016 base layers
# merge as we read in
canopy.2016 <- merge(
  
  rast(paste0(dir.rast, "lidar_metrics/2016/1_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/2016/2_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/2016/3_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/2016/4_10.tif"))
  
)

# pre metrics
# merge as we read in
canopy.pre <- merge(
  
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_1A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_1B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_1C_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_2A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_2B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_2C_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_3A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_3B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_3C_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_4A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_4B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Pre/pre_4C_10.tif"))
  
)

# post metrics
# merge as we read in
canopy.post <- merge(
  
  rast(paste0(dir.rast, "lidar_metrics/Post/1A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/1B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/1C_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/2A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/2B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/2C_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/3A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/3B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/3C_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/4A_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/4B_10.tif")),
  rast(paste0(dir.rast, "lidar_metrics/Post/4C_10.tif"))
  
)

# merge and project
# pre
rast.ch.pre <- merge(canopy.pre$ch, canopy.2016$ch) |> project(rast.cover.pre)
rast.cc.pre <- merge(canopy.pre$cc, canopy.2016$cc) |> project(rast.cover.pre)

# post
rast.ch.post <- merge(canopy.post$ch, canopy.2016$ch) |> project(rast.cover.post)
rast.cc.post <- merge(canopy.post$cc, canopy.2016$cc) |> project(rast.cover.post)

# ______________________________________________________________________________
# 4. Landscape metrics on cover types ----
# ______________________________________________________________________________
# 4a. Convert to factors ---
# ______________________________________________________________________________

levels(rast.cover.pre) <- data.frame(
  
  id = 1:8,
  type = c("stand initiation",
           "doghair",
           "wetland",
           "open forest",
           "jackstraw",
           "spruce-fir",
           "other mature",
           "open")
  
)

levels(rast.cover.post) <- data.frame(
  
  id = 1:9,
  type = c("stand initiation",
           "doghair",
           "wetland",
           "open forest",
           "jackstraw",
           "spruce-fir",
           "other mature",
           "open",
           "thinned doghair")
  
)

# ______________________________________________________________________________
# 4b. Moving window calculations ---

# let's do patch diversity

# ______________________________________________________________________________

# define metrics to calculate
which.metrics <- c("lsm_l_shdi")

# define moving windows
# these must be matrices with odd rows/cols
# each pixel is 10 m
mw.50 <- matrix(1, nrow = 11, ncol = 11)
mw.100 <- matrix(1, nrow = 21, ncol = 21)

# calculate metric(s) within moving windows
#lsm.50.pre <- window_lsm(rast.cover.pre, window = mw.50, what = which.metrics)
#lsm.50.post <- window_lsm(rast.cover.post, window = mw.50, what = which.metrics)
#lsm.100.pre <- window_lsm(rast.cover.pre, window = mw.100, what = which.metrics)
#lsm.100.post <- window_lsm(rast.cover.post, window = mw.100, what = which.metrics)

# save to file
writeRaster(lsm.50.pre$layer_1$lsm_l_shdi, "data_raster/shdi_50_pre.tif", overwrite = T)
writeRaster(lsm.50.post$layer_1$lsm_l_shdi, "data_raster/shdi_50_post.tif", overwrite = T)
writeRaster(lsm.100.pre$layer_1$lsm_l_shdi, "data_raster/shdi_100_pre.tif", overwrite = T)
writeRaster(lsm.100.post$layer_1$lsm_l_shdi, "data_raster/shdi_100_post.tif", overwrite = T)

# read
rast.shdi.50.pre <- rast("data_raster/shdi_50_pre.tif")
rast.shdi.50.post <- rast("data_raster/shdi_50_post.tif")
rast.shdi.100.pre <- rast("data_raster/shdi_100_pre.tif")
rast.shdi.100.post <- rast("data_raster/shdi_100_post.tif")

# ______________________________________________________________________________
# 5. Resample as needed ----

# our target extent:
ext(rast.cover.pre)

# we'll also snap together

# ______________________________________________________________________________

rast.all <- c(
  
  # distance
  resample(rast.dEdge, rast.cover.pre),
  resample(rast.dOpen, rast.cover.pre),
  resample(rast.dDM, rast.cover.pre),
  
  # canopy
  resample(rast.ch.pre, rast.cover.pre),
  resample(rast.ch.post, rast.cover.pre),
  resample(rast.cc.pre, rast.cover.pre),
  resample(rast.cc.post, rast.cover.pre),
  
  # landscape
  resample(rast.shdi.50.pre, rast.cover.pre),
  resample(rast.shdi.50.post, rast.cover.pre),
  resample(rast.shdi.100.pre, rast.cover.pre),
  resample(rast.shdi.100.post, rast.cover.pre),
  
  # veg models
  resample(rast.stem.pre, rast.cover.pre),
  resample(rast.stem.post, rast.cover.pre),
  resample(rast.vo.pre, rast.cover.pre),
  resample(rast.vo.post, rast.cover.pre),
  
  # topography
  resample(rast.twi, rast.cover.pre),
  resample(rast.vrm, rast.cover.pre)
  
  # treatment
  #resample(rast.dPil, rast.cover.pre),
  #resample(rast.dRet, rast.cover.pre),
  #resample(rast.dUnitInt, rast.cover.pre)
  
)

# change names
names(rast.all) <- c("dEdge", "dOpen", "dDM",
                     "ch.pre", "ch.post", "cc.pre", "cc.post",
                     "shdi.50.pre", "shdi.50.post", "shdi.100.pre", "shdi.100.post",
                     "stem.pre", "stem.post", "vo.pre", "vo.post",
                     "twi", "vrm")

plot(rast.all)

# ______________________________________________________________________________
# 6. Write rasters ----
# ______________________________________________________________________________

writeRaster(rast.all, "data_raster/rast_all.tif", overwrite = T)
