# PROJECT: Habitat selection
# SCRIPT: 02 - Sample covariates
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 21 Apr 2026
# LAST MODIFIED: 09 Jun 2026
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
rast.prox <- rast("data_raster/unit_prox.tif")

names(rast.prox) <- "closest"

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
# 5. Unit proximity ----
# ______________________________________________________________________________

# function
add_corrected_trt <- function (.x) {
  
  .x.prox <- .x |> st_as_sf(coords = c("x", "y"), crs = "epsg:32611") %>%
    
    # bind in extracted values
    bind_cols(
      
      .,
      
      # extract covariate values
      extract(x = rast.prox, y = .) %>%
        
        # remove "ID"
        dplyr::select(-ID)
      
    )
  
  # group by TSP and pick the mode
  # Source - https://stackoverflow.com/a/25635740
  # Posted by jprockbelly, modified by community. See post 'Timeline' for change history
  # Retrieved 2026-06-02, License - CC BY-SA 3.0
  Mode <- function(x, na.rm = FALSE) {
    if(na.rm){
      x = x[!is.na(x)]
    }
    
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }

  .x.group <- .x.prox |> group_by(track_season_post) |>
    
    summarize(closest.mode = Mode(closest)) |>
    
    # corrected trt
    mutate(c.trt = case_when(
      
      closest.mode %in% c(1, 5, 8, 10) ~ "RET",  # 1A, 2B, 3B, 4A
      closest.mode %in% c(2, 4, 7, 11) ~ "PIL",  # 1B, 2A, 3A, 3B
      closest.mode %in% c(3, 6, 9, 12) ~ "CTRL"  # 1C, 2C, 3C, 4C
      
     )
    
    ) |>
    
    # drop columns
    dplyr::select(-c(geometry, closest.mode))
    
  # join in
  .x.1 <- .x |> left_join(.x.group) |> dplyr::select(-geometry)
  
  return(.x.1)
  
}

pts.corrected <- add_corrected_trt(pts)

# how many changed?
test <- pts.corrected |> 
  
  group_by(track_season_post) |>
  
  slice(1)

sum(test$trt != test$c.trt)  # only changed 3!

# add to sampled dfs
pts.off.1 <- cbind(pts.off, c.trt = pts.corrected$c.trt[1:nrow(pts.off)]) 
pts.on.1 <- cbind(pts.on, c.trt = pts.corrected$c.trt[(nrow(pts.off) + 1):(nrow(pts.off) + nrow(pts.on))]) 

# ______________________________________________________________________________
# 4. Write to file ----
# ______________________________________________________________________________

saveRDS(pts.off.1, "data_cleaned/data_off.rds")
saveRDS(pts.on.1, "data_cleaned/data_on.rds")
