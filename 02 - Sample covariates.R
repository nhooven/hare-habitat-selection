# PROJECT: Habitat selection
# SCRIPT: 02 - Sample covariates
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 21 Apr 2026
# LAST MODIFIED: 29 May 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(terra)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# rasters
rast.all <- rast("data_raster/rast_all.tif")

# used/background points
pts <- readRDS("D:/hare_project/data_analysis/General/hare-gps-processing-new/data_cleaned/use_background.rds")

# ______________________________________________________________________________
# 3. Sample and split by season ----
# ______________________________________________________________________________

# promote to sf
pts.1 <- st_as_sf(pts, coords = c("x", "y"), crs = "epsg:32611") %>%
  
  # bind in extracted values
  bind_cols(
    
    .,
    
    # extract covariate values
    extract(x = rast.all, y = .) %>%
  
    # remove "ID"
    dplyr::select(-ID)
    
  ) |>
  
  as.data.frame() %>%
  
  # drop geometry
  dplyr::select(-geometry)
  
# split by season
pts.off <- pts.1 %>% filter(season == "off")
pts.on <- pts.1 %>% filter(season == "on")

# ______________________________________________________________________________
# 4. Write to file ----
# ______________________________________________________________________________

saveRDS(pts.off, "data_cleaned/data_off.rds")
saveRDS(pts.on, "data_cleaned/data_on.rds")
