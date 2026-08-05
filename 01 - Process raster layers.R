# PROJECT: Habitat selection
# SCRIPT: 01 - Process raster layers
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 20 Apr 2026
# COMPLETED: 21 Apr 2026
# LAST MODIFIED: 05 Aug 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(terra)
library(landscapemetrics)
library(sf)

# ______________________________________________________________________________
# 2. Read spatial data ----
# ______________________________________________________________________________

# raster directory
dir.rast <- "D:/hare_project/data_spatial/Rasters/"

# cover type
rast.cover.pre <- rast(paste0(dir.rast, "cover_type/cover_type_pre.tif"))
rast.cover.post <- rast(paste0(dir.rast, "cover_type/cover_type_post.tif"))

# distance 
rast.dEdge <- rast(paste0(dir.rast, "dist_rasters/dEdge.tif"))

# vegetation models
rast.stem.pre <- rast(paste0(dir.rast, "veg_pred/RF/stem_pre_new.tif"))
rast.stem.post <- rast(paste0(dir.rast, "veg_pred/RF/stem_post_new.tif"))
rast.vo.pre <- rast(paste0(dir.rast, "veg_pred/RF/vo_pre_new.tif"))
rast.vo.post <- rast(paste0(dir.rast, "veg_pred/RF/vo_post_new.tif"))

# topography
rast.twi <-  rast(paste0(dir.rast, "Topography/twi_10.tif"))
rast.vrm <-  rast(paste0(dir.rast, "Topography/vrmL_10.tif"))

# 30 m canopy cover
rast.cc.pre <- rast(paste0(dir.rast, "TCC/tcc_2023_utm.tif"))
rast.cc.post <- rast(paste0(dir.rast, "TCC/tcc_2024_utm.tif"))

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

# units 
units <- st_read("D:/hare_project/data_spatial/Units/units_fixed_utm/units_fixed_utm.shp") |>
  
  dplyr::select(name, geometry) |>
  
  arrange(name) |>
  
  st_transform(crs(canopy.pre))

# ______________________________________________________________________________
# 3a. Canopy cover ----
# ______________________________________________________________________________

# change anything > 100 to NA
rast.cc.pre[rast.cc.pre > 100] <- NA
rast.cc.post[rast.cc.post > 100] <- NA

plot(rast.cc.pre)
plot(rast.cc.post)

# ______________________________________________________________________________
# 3b. Canopy height ----
# ______________________________________________________________________________

# mask to unit boundaries
canopy.pre.ch <- mask(crop(canopy.pre$ch, vect(units)), vect(units))
canopy.post.ch <- mask(crop(canopy.post$ch, vect(units)), vect(units))

# merge and project
rast.ch.pre <- merge(canopy.pre.ch, canopy.2016$ch) |> project(rast.cover.pre)
rast.ch.post <- merge(canopy.post.ch, canopy.2016$ch) |> project(rast.cover.post)

# fill in blanks (between lidar tiles)
rast.ch.pre <- focal(rast.ch.pre, w = 3, fun = mean, na.policy = "only", na.rm = T)
rast.ch.post <- focal(rast.ch.post, w = 3, fun = mean, na.policy = "only", na.rm = T)

# clamp rast.ch.post to highest real quantile
ch.post.quant <- quantile(values(rast.ch.post), na.rm = T, prob = 0.9999)

rast.ch.post <- clamp(rast.ch.post, upper = ch.post.quant)

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
# 6. Resample as needed ----

# our target extent:
ext(rast.cover.pre)

# we'll also snap together

# ______________________________________________________________________________

rast.all <- c(
  
  # conditions
  resample(rast.twi, rast.cover.pre),
  resample(rast.vrm, rast.cover.pre),
  
  # structure 
  resample(rast.stem.pre, rast.cover.pre),
  resample(rast.stem.post, rast.cover.pre),
  resample(rast.vo.pre, rast.cover.pre),
  resample(rast.vo.post, rast.cover.pre),
  resample(rast.ch.pre, rast.cover.pre),
  resample(rast.ch.post, rast.cover.pre),
  resample(rast.dEdge, rast.cover.pre),
  resample(rast.cc.pre, rast.cover.pre),
  resample(rast.cc.post, rast.cover.pre)
  
)

# change names
names(rast.all) <- c(
  
  "twi", "vrm",
  "stem.pre", "stem.post", 
  "vo.pre", "vo.post", 
  "ch.pre", "ch.post", 
  "dEdge", 
  "cc.pre", "cc.post"
  
  )

terra::plot(rast.all)

# ______________________________________________________________________________
# 6. Write rasters ----
# ______________________________________________________________________________

writeRaster(rast.all, "data_raster/rast_all.tif", overwrite = T)
