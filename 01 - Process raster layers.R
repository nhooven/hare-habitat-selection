# PROJECT: Habitat selection
# SCRIPT: 01 - Process raster layers
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 20 Apr 2026
# COMPLETED: 
# LAST MODIFIED: 
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(terra)

# ______________________________________________________________________________
# 2. Read rasters ----
# ______________________________________________________________________________

# raster directory
dir.rast <- "D:/hare_project/data_spatial/Rasters/"

# cover type
rast.cover.pre <- rast(paste0(dir.rast, "cover_type/cover_type_pre.tif"))
rast.cover.post <- rast(paste0(dir.rast, "cover_type/cover_type_post.tif"))

# distance to edge
rast.dEdge <- rast(paste0(dir.rast, "dist_rasters/dEdge.tif"))

# vegetation models
rast.stem.pre <- rast(paste0(dir.rast, "veg_pred/RF/pre_stem.tif"))
rast.stem.post <- rast(paste0(dir.rast, "veg_pred/RF/post_stem.tif"))
rast.vo.pre <- rast(paste0(dir.rast, "veg_pred/RF/pre_vo.tif"))
rast.vo.post <- rast(paste0(dir.rast, "veg_pred/RF/post_vo.tif"))

# topography
rast.twi <-  rast(paste0(dir.rast, "Topography/twi_10.tif"))
rast.slope <-  rast(paste0(dir.rast, "Topography/dtm_slope_deg_10.tif"))

# PCT
rast.dPil <-  rast(paste0(dir.rast, "PCT/dPiles.tif"))
rast.dRet <-  rast(paste0(dir.rast, "PCT/dRet.tif"))
rast.dUnitInt <- rast(paste0(dir.rast, "PCT/dUnitInt.tif"))

# ______________________________________________________________________________
# 3. Reproject as needed ----

# our target extent:
ext(rast.cover.pre)

# ______________________________________________________________________________