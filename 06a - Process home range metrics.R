# PROJECT: Habitat selection
# SCRIPT: 06a - Process home range metrics
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 11 Jun 2026
# COMPLETED: 11 Jun 2026
# LAST MODIFIED: 11 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(ctmm)
library(sf)
library(terra)
library(landscapemetrics)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# directory
gps.dir <- "D:/hare_project/data_analysis/General/hare-gps-processing-new/"

# AKDE models
all.akde <- readRDS(paste0(gps.dir, "data_cleaned/all_akde.rds"))

# model selection summaries
# we need these to attribute the correct AKDEs
model.select <- readRDS(paste0(gps.dir, "data_cleaned/all_model_select.rds"))

# raster directory
dir.rast <- "D:/hare_project/data_spatial/Rasters/"

# cover type
rast.cover.pre <- rast(paste0(dir.rast, "cover_type/cover_type_pre.tif"))
rast.cover.post <- rast(paste0(dir.rast, "cover_type/cover_type_post.tif"))

# ______________________________________________________________________________
# 3. Sample HR metrics ----

indivs.split <- model.select |> 
  
  dplyr::select(track_season_post, i) |>

  rename(TSPID = track_season_post) |>
  
  group_by(TSPID) |> 
  
  slice(1) |>
  
  split(~TSPID)

# ______________________________________________________________________________
# 3a. Function ----
# ______________________________________________________________________________

# this will take a list of TSPs
sample_hr_metrics <- function (x) {
  
  # extract correct HR
  focal.akde <- all.akde[[x$i]]
  
  # extract 99% contour
  focal.contour <- as.sf(focal.akde, level.UD = 0.99)[2, ] |>
    
    # convert to UTM
    st_transform(crs = "epsg:32611")
  
  # split string (easy way of choosing which raster to use)
  if("PRE" %in% str_split(x$TSPID, pattern = "_")[[1]]) { rast.focal <- rast.cover.pre }
  if("POST" %in% str_split(x$TSPID, pattern = "_")[[1]]) { rast.focal <- rast.cover.post }
  
  # crop and mask
  rast.mask <- mask(crop(rast.focal, vect(focal.contour)), vect(focal.contour))
  
  # calculate PLAND and SHDI
  pland <- lsm_c_pland(rast.mask)
  shdi <- lsm_l_shdi(rast.mask)
  
  # account for zeroes (9 doesn't matter here)
  p.OM4 <- ifelse(4 %in% pland$class, pland$value[pland$class == 4], 0)
  p.JsM <- ifelse(5 %in% pland$class, pland$value[pland$class == 5], 0)
  p.OM <- p.OM4 + p.JsM
  p.MsM <- ifelse(6 %in% pland$class, pland$value[pland$class == 6], 0)
  p.GM <- ifelse(7 %in% pland$class, pland$value[pland$class == 7], 0)
  p.DM <- p.MsM + p.GM
  
  # pack into df
  focal.metrics <- data.frame(TSPID = x$TSPID,
                              pOM = p.OM,
                              pDM = p.DM,
                              shdi = shdi$value)
  
  return(focal.metrics)
  
}

# ______________________________________________________________________________
# 3b. Use function ----
# ______________________________________________________________________________

all.hr.metrics <- do.call(rbind, lapply(indivs.split, sample_hr_metrics))

# ______________________________________________________________________________
# 4. Save to file ----
# ______________________________________________________________________________

saveRDS(all.hr.metrics, "data_cleaned/all_hr_metrics.rds")
