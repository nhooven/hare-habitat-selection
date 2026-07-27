# PROJECT: Habitat selection
# SCRIPT: FIG - Example map
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 Jul 2026
# COMPLETED: 
# LAST MODIFIED: 27 Jul 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(cowplot)
library(spatialEco)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________
# 2a. Directories ----
# ______________________________________________________________________________

# spatial data
spat.dir <- "D:/hare_project/data_spatial/"

# ______________________________________________________________________________
# 2b. Vector data ----
# ______________________________________________________________________________

# PCT units
pct.units <- st_read(paste0(spat.dir, "Units/units_fixed_utm/units_fixed_utm.shp"))

# piles and retention
piles <- st_read(paste0(spat.dir, "Ground-truthing/all_piles.shp"))
ret <- st_read(paste0(spat.dir, "Ground-truthing/all_refugia.shp"))

# ______________________________________________________________________________
# 2c. Raster data ----
# ______________________________________________________________________________

cover <- rast(paste0(spat.dir, "Rasters/cover_type/cover_type_pre.tif"))
vo <- rast(paste0(spat.dir, "Rasters/veg_pred/RF/vo_post_new.tif"))

# ______________________________________________________________________________
# 3. Coordinate reference systems ----
# ______________________________________________________________________________

crs.utm <- "epsg:32611"

# ______________________________________________________________________________
# 4. Clean data ----
# ______________________________________________________________________________

# units
wr.units <- pct.units |> dplyr::select(name, geometry) |> 
  
  filter(name %in% c("1A", "1B", "1C"))

# piles
wr.piles <- piles |> select(name, geometry) |>
  
  mutate(unit = substr(name, 1, 2)) |>
  
  filter(unit == "1B")

# piles
wr.ret <- ret |> select(name, geometry) |>
  
  mutate(unit = substr(name, 1, 2)) |>
  
  filter(unit == "1A")

# ______________________________________________________________________________
# 5. Bounding box and crop rasters ----
# ______________________________________________________________________________

wr.bbox <- wr.units |> bbox_poly() |> st_buffer(dist = 100) |> bbox_poly()

cover.crop <- cover |> crop(wr.bbox)
vo.crop <- vo |> crop(wr.bbox)

# ______________________________________________________________________________
# 6. Plots ----
# ______________________________________________________________________________
# 6a. Basemap ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_spatraster(data = vo.crop) +
  
  # unit boundaries
  geom_sf(data = wr.units,
          fill = NA,
          color = "white") +
  
  # piles
  geom_sf(data = wr.piles,
          color = "white",
          shape = 3) +
  
  # patches
  geom_sf(data = wr.ret,
          fill = NA,
          color = "white") +
  
  # colors
  scale_fill_viridis_c() +
  
  # theme
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        
        # legend
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.2),
        legend.direction = "horizontal") -> map.base

# ______________________________________________________________________________
# 6b. Doghair (pre-thinned) ----

cover.dh <- cover.crop |> classify(rcl = matrix(data = c(1, NA,
                                                         2, 1,
                                                         3, NA,
                                                         4, NA,
                                                         5, NA,
                                                         6, NA,
                                                         7, NA,
                                                         8, NA),
                                                byrow = T,
                                                nrow = 8))

# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_spatraster(data = cover.dh) +
  
  # color
  scale_fill_continuous(palette = c("green3", "green3"),
                        na.value = "gray90") +
  
  # theme
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
