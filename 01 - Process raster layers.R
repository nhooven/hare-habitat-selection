# PROJECT: Habitat selection
# SCRIPT: 01 - Process raster layers
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 20 Apr 2026
# COMPLETED: 21 Apr 2026
# LAST MODIFIED: 09 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(terra)
library(landscapemetrics)
library(sf)

# ______________________________________________________________________________
# 2. Read rasters ----
# ______________________________________________________________________________

# raster directory
dir.rast <- "D:/hare_project/data_spatial/Rasters/"

# cover type
rast.cover.pre <- rast(paste0(dir.rast, "cover_type/cover_type_pre.tif"))
rast.cover.post <- rast(paste0(dir.rast, "cover_type/cover_type_post.tif"))

# distance 
rast.dOM <- rast(paste0(dir.rast, "dist_rasters/dOM_new.tif"))
rast.dDM <- rast(paste0(dir.rast, "dist_rasters/dDM.tif"))

# vegetation models
rast.stem.pre <- rast(paste0(dir.rast, "veg_pred/RF/pre_stem.tif"))
rast.stem.post <- rast(paste0(dir.rast, "veg_pred/RF/post_stem.tif"))
rast.vo.pre <- rast(paste0(dir.rast, "veg_pred/RF/pre_vo.tif"))
rast.vo.post <- rast(paste0(dir.rast, "veg_pred/RF/post_vo.tif"))

# topography
rast.twi <-  rast(paste0(dir.rast, "Topography/twi_10.tif"))
rast.vrm <-  rast(paste0(dir.rast, "Topography/vrmL_10.tif"))

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

# fill in blanks (between lidar tiles)
rast.ch.pre <- focal(rast.ch.pre, w = 3, fun = mean, na.policy = "only", na.rm = T)
rast.cc.pre <- focal(rast.cc.pre, w = 3, fun = mean, na.policy = "only", na.rm = T)
rast.ch.post <- focal(rast.ch.post, w = 3, fun = mean, na.policy = "only", na.rm = T)
rast.cc.post <- focal(rast.cc.post, w = 3, fun = mean, na.policy = "only", na.rm = T)

# clamp rast.ch.post to highest real quantile
ch.post.quant <- quantile(values(rast.ch.post), na.rm = T, prob = 0.9999)

rast.ch.post <- clamp(rast.ch.post, upper = ch.post.quant)

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

# hard edge density
# mature forest vs all other classes (thus this won't change with treatment)
#rast.cover.edge <- rast.cover.pre

#mat.reclass.edge <- matrix(c(1, 1,
#                             2, 1,
#                             3, 1,
#                             4, 2,
#                             5, 2,
#                             6, 2,
#                             7, 2,
#                             8, 1),
#                           nrow = 8,
#                           byrow = T)

#rast.cover.edge <- classify(rast.cover.edge, rcl = mat.reclass.edge)

# ______________________________________________________________________________
# 4b. Moving window calculations ---
# ______________________________________________________________________________

# define moving windows
# these must be matrices with odd rows/cols
# each pixel is 10 m
#mw.50 <- matrix(1, nrow = 11, ncol = 11)
#mw.100 <- matrix(1, nrow = 21, ncol = 21)

# calculate metric(s) within moving windows
# patch diversity
#shdi.50.pre <- window_lsm(rast.cover.pre, window = mw.50, what = "lsm_l_shdi")
#shdi.50.post <- window_lsm(rast.cover.post, window = mw.50, what = "lsm_l_shdi")
#shdi.100.pre <- window_lsm(rast.cover.pre, window = mw.100, what = "lsm_l_shdi")
#shdi.100.post <- window_lsm(rast.cover.post, window = mw.100, what = "lsm_l_shdi")

# save to file
#writeRaster(lsm.50.pre$layer_1$lsm_l_shdi, "data_raster/shdi_50_pre.tif", overwrite = T)
#writeRaster(lsm.50.post$layer_1$lsm_l_shdi, "data_raster/shdi_50_post.tif", overwrite = T)
#writeRaster(lsm.100.pre$layer_1$lsm_l_shdi, "data_raster/shdi_100_pre.tif", overwrite = T)
#writeRaster(lsm.100.post$layer_1$lsm_l_shdi, "data_raster/shdi_100_post.tif", overwrite = T)

# read
#rast.shdi.50.pre <- rast("data_raster/shdi_50_pre.tif")
#rast.shdi.50.post <- rast("data_raster/shdi_50_post.tif")
rast.shdi.100.pre <- rast("data_raster/shdi_100_pre.tif")
rast.shdi.100.post <- rast("data_raster/shdi_100_post.tif")

# edge density (m/ha)
# write to raster to calculate in Fragstats
#writeRaster(rast.cover.edge, "data_raster/cover_edge.tif", overwrite = T)

# read (100 m)
rast.ed <- rast("data_raster/ed_100.tif")

# replace -999 with NA
rast.ed <- subst(rast.ed, -999, NA)

# ______________________________________________________________________________
# 5. Closest unit raster ----

# this is a sensible way to empirically determine which "unit" each animal
# is most closely aligned to

# first, we'll calculate distance rasters for every unit, then
# select the (integer) value for every LSF pixel for the smallest distance

# then later we can attribute a value to every available point and take the mode
# as the "corrected" unit

# ______________________________________________________________________________

# unit shapefile
units <- st_read("D:/hare_project/data_spatial/Units/units_fixed_utm/units_fixed_utm.shp") |>
  
  dplyr::select(name, geometry) |>
  
  arrange(name)

all.dists <- rast()

for (i in 1:12) {
  
  focal.unit <- units |> slice(i) |> vect()
  
  focal.dist <- distance(rast.cover.pre, focal.unit)
  
  all.dists <- c(all.dists, focal.dist)
  
}

# names
names(all.dists) <- units$name

# find which.min
# make matrix of values
val.mat <- matrix(NA, nrow = length(values(all.dists[[1]])), ncol = 12)

for (i in 1:12) {
  
  val.mat[ , i] <- values(all.dists[[i]])
  
}

min.dists <- apply(val.mat, 1, which.min)

rast.min.dists <- rast(all.dists[[1]], vals = min.dists)

plot(rast.min.dists)

# write
writeRaster(rast.min.dists, "data_raster/unit_prox.tif")

# ______________________________________________________________________________
# 5. Resample as needed ----

# our target extent:
ext(rast.cover.pre)

# we'll also snap together

# ______________________________________________________________________________

rast.all <- c(
  
  # landscape
  resample(rast.dOM, rast.cover.pre),
  resample(rast.dDM, rast.cover.pre),
  resample(rast.ed, rast.cover.pre),
  resample(rast.shdi.100.pre, rast.cover.pre),
  resample(rast.shdi.100.post, rast.cover.pre),
  
  # canopy
  resample(rast.ch.pre, rast.cover.pre),
  resample(rast.ch.post, rast.cover.pre),
  resample(rast.cc.pre, rast.cover.pre),
  resample(rast.cc.post, rast.cover.pre),
  
  # veg models
  resample(rast.stem.pre, rast.cover.pre),
  resample(rast.stem.post, rast.cover.pre),
  resample(rast.vo.pre, rast.cover.pre),
  resample(rast.vo.post, rast.cover.pre),
  
  # topography
  resample(rast.twi, rast.cover.pre),
  resample(rast.vrm, rast.cover.pre)
  
)

# change names
names(rast.all) <- c(
  
  "dOM", "dDM", "ed",
  "shdi.pre", "shdi.post",
  "ch.pre", "ch.post", "cc.pre", "cc.post",
  "stem.pre", "stem.post", "vo.pre", "vo.post",
  "twi", "vrm"
  
  )

plot(rast.all)

# ______________________________________________________________________________
# 6. Write rasters ----
# ______________________________________________________________________________

writeRaster(rast.all, "data_raster/rast_all.tif", overwrite = T)
