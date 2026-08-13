# PROJECT: Habitat selection
# SCRIPT: FIG - Covariate and relocation maps
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 13 Aug 2026
# COMPLETED: 13 Aug 2026
# LAST MODIFIED: 13 Aug 2026
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
library(ggspatial)       # scale bar

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
pct.units <- st_read(paste0(spat.dir, "Units/units_fixed_utm/units_fixed_utm.shp")) |>
  
  dplyr::select(name, geometry) |>
  
  arrange(name)

# piles and retention
piles <- st_read(paste0(spat.dir, "Ground-truthing/all_piles.shp")) |>
  
  select(name, geometry) |>
  
  mutate(unit = substr(name, 1, 2)) |>
  
  st_transform(crs = "epsg:32611")

ret <- st_read(paste0(spat.dir, "Ground-truthing/all_refugia.shp")) |> 
  
  select(name, geometry) |>
  
  mutate(unit = substr(name, 1, 2))

# relocations
pts <- readRDS("D:/hare_project/data_analysis/General/hare-gps-processing-new/data_cleaned/use_background.rds") |>
  
  filter(case == 1) |>
  
  st_as_sf(coords = c("x", "y"),
           crs = "epsg:32611") |>
  
  # treatment
  mutate(trt.1 = factor(case_when(year == "PRE" ~ "U",
                                  year %in% c("POST1", "POST2") & trt == "CTRL" ~ "U",
                                  year %in% c("POST1", "POST2") & trt == "RET" ~ "R",
                                  year %in% c("POST1", "POST2") & trt == "PIL" ~ "P"),
                        levels = c("U", "R", "P")))

# ______________________________________________________________________________
# 2c. Raster data ----
# ______________________________________________________________________________

rast.all <- rast(paste0(getwd(), "/data_raster/rast_all.tif"))

# ______________________________________________________________________________
# 3. Bounding boxes ----

# we'll use the same ones as the Ch 2 map

# ______________________________________________________________________________

# units
bbox.units <- rbind(
  
  bbox_poly(pct.units |> filter(name %in% c("1A", "1B", "1C"))),
  bbox_poly(pct.units |> filter(name %in% c("2A", "2B", "2C"))),
  bbox_poly(pct.units |> filter(name %in% c("3A", "3B", "3C"))),
  bbox_poly(pct.units |> filter(name %in% c("4A", "4B", "4C")))
  
)

# buffer them
bbox.units.buff <- rbind(
  
  bbox_poly(st_buffer(bbox.units[1, ], dist = 300)),
  bbox_poly(st_buffer(bbox.units[2, ], dist = 300)),
  bbox_poly(st_buffer(bbox.units[3, ], dist = 300)),
  bbox_poly(st_buffer(bbox.units[4, ], dist = 300))
  
)

# ______________________________________________________________________________
# 4. Function ----

# the left panel will be an enlarged map with unit boundaries,
# refugia, and all relocations (probably different colors)
# should there be a basemap? Maybe not...

# the right panel will be the same covariates as the WR example map

# ______________________________________________________________________________

plot_info_panel <- function (.cluster) {
  
  # correct subsets
  if (.cluster == 1) {
    
    clust.bbox <- bbox.units.buff[1, ]
    clust.units <- pct.units[1:3, ]
    clust.piles <- piles |> filter(unit == "1B")
    clust.ret <- ret |> filter(unit == "1A")
    clust.pts <- pts |> mutate(cluster = substr(site, 1, 1)) |> filter(cluster == 1)
    
  }
  
  if (.cluster == 2) {
    
    clust.bbox <- bbox.units.buff[2, ]
    clust.units <- pct.units[4:6, ]
    clust.piles <- piles |> filter(unit == "2A")
    clust.ret <- ret |> filter(unit == "2B")
    clust.pts <- pts |> mutate(cluster = substr(site, 1, 1)) |> filter(cluster == 2)
    
  }
  
  if (.cluster == 3) {
    
    clust.bbox <- bbox.units.buff[3, ]
    clust.units <- pct.units[7:9, ]
    clust.piles <- piles |> filter(unit == "3A")
    clust.ret <- ret |> filter(unit == "3B")
    clust.pts <- pts |> mutate(cluster = substr(site, 1, 1)) |> filter(cluster == 3)
    
  }
  
  if (.cluster == 4) {
    
    clust.bbox <- bbox.units.buff[4, ]
    clust.units <- pct.units[10:12, ]
    clust.piles <- piles |> filter(unit == "4B")
    clust.ret <- ret |> filter(unit == "4A")
    clust.pts <- pts |> mutate(cluster = substr(site, 1, 1)) |> filter(cluster == 4)
    
  }
  
  # ____________________________________________________________________________
  # RIGHT PANEL
  # ____________________________________________________________________________
  
  # plot
  ggplot() +
    
    theme_bw() +
    
    # bounding box
    geom_sf(data = clust.bbox,
            fill = NA,
            color = NA) +
    
    # unit boundaries
    geom_sf(data = clust.units,
            fill = NA,
            color = "black",
            linewidth = 0.3) +
    
    # piles
    geom_sf(data = clust.piles,
            color = "black",
            shape = 3,
            size = 0.6,
            stroke = 0.1) +
    
    # patches
    geom_sf(data = clust.ret,
            fill = NA,
            color = "gray15",
            linewidth = 0.3) +
    
    # relocations
    geom_sf(data = clust.pts,
            aes(color = season,
                shape = trt.1),
            alpha = 0.65,
            size = 0.4,
            stroke = 0.1) +
    
    coord_sf(expand = F) +
    
    # colors and shapes
    scale_color_manual(values = c("green4", "dodgerblue3")) +
    scale_shape_manual(values = c(1, 2, 5)) +
  
  # scale bar
  annotation_scale(aes(style = "ticks"),
                   location = "br",
                   line_width = 0.5,
                   height = unit(0.1, "cm"),
                   text_cex = 0.6,
                   pad_x = unit(0.7, "cm"),
                   pad_y = unit(0.5, "cm"),
                   line_col = "black",
                   text_col = "black") +
    
    # theme
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          
          # legend
          legend.position = "top",
          legend.direction = "horizontal",
          legend.background = element_rect(fill = NA),
          legend.frame = element_rect(color = "white"),
          legend.text.position = "top",
          legend.title = element_blank()) +
    
    # legend guides
    guides(color = guide_legend(override.aes = list(size = 2,
                                                    alpha = 1)),
           shape = guide_legend(override.aes = list(size = 2,
                                                    alpha = 1))) -> map.base
  
  # ____________________________________________________________________________
  # LEFT PANEL
  # ____________________________________________________________________________
  
  # crop rasters
  # crop to bbox
  rast.all.crop <- crop(rast.all, clust.bbox)
  
  # standardize and stack rasters 
  rast.focal.1 <- c(scale(rast.all.crop$vo.pre),
                    scale(rast.all.crop$vo.post),
                    scale(rast.all.crop$stem.pre),
                    scale(rast.all.crop$stem.post),
                    scale(rast.all.crop$ch.pre),
                    scale(rast.all.crop$ch.post),
                    scale(rast.all.crop$cc.pre),
                    scale(rast.all.crop$cc.post))
  
  rast.focal.2 <- c(scale(rast.all.crop$dEdge))
  
  rast.cond <- c(scale(rast.all.crop$twi),
                 scale(rast.all.crop$vrm))
  
  # names
  names(rast.focal.1) <- c("VO (pre)", "VO (post)", "STEM (pre)", "STEM (post)", 
                           "CH (pre)", "CH (post)", "CC (pre)", "CC (post)")
  names(rast.focal.2) <- c("dEdge")
  names(rast.cond) <- c("TWI", "VRM")
  
  # determine the appropriate range (quantiles)
  range.min <- min(quantile(values(rast.focal.1), 0.01, na.rm = T), 
                   quantile(values(rast.focal.2), 0.01, na.rm = T), 
                   quantile(values(rast.cond), 0.01, na.rm = T))
  
  range.max <- max(quantile(values(rast.focal.1), 0.99, na.rm = T), 
                   quantile(values(rast.focal.2), 0.99, na.rm = T), 
                   quantile(values(rast.cond), 0.99, na.rm = T)) 
  
  # clamp
  rast.focal.1 <- clamp(rast.focal.1, range.min, range.max)
  rast.focal.2 <- clamp(rast.focal.2, range.min, range.max)
  rast.cond <- clamp(rast.cond, range.min, range.max)
  
  # covariate plots
  ggplot() +
    
    theme_bw() +
    
    facet_wrap(~ lyr, 
               nrow = 4,
               strip.position = "right",
               labeller = labeller(lyr = c("VO (pre)" = "", 
                                           "VO (post)" = "VO", 
                                           "STEM (pre)" = "", 
                                           "STEM (post)" = "STEM", 
                                           "CH (pre)" = "", 
                                           "CH (post)" = "CH", 
                                           "CC (pre)" = "", 
                                           "CC (post)" = "CC"))) +
    
    # basemap
    geom_spatraster(data = rast.focal.1) +
    
    # unit boundaries
    geom_sf(data = clust.units,
            fill = NA,
            color = "white",
            linewidth = 0.2) +
    
    coord_sf(expand = F) +
    
    # theme
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          
          legend.position = "top",
          legend.key.width = unit(0.45, "cm"),
          legend.key.height = unit(0.2, "cm"),
          legend.title = element_text(size = 8),
          legend.title.position = "right",
          legend.text = element_text(size = 7,
                                     margin = margin(t = 3, "cm")),
          
          panel.spacing.x = unit(0.005, "cm"),
          strip.background = element_rect(color = NA,
                                          fill = NA),
          strip.text = element_text(size = 7,
                                    hjust = 0)) +
    
    # colors
    scale_fill_viridis_c(name = "SD") -> plot.focal.1
  
  # focal 2
  ggplot() +
    
    theme_bw() +
    
    facet_wrap(~ lyr, 
               nrow = 1,
               strip.position = "right") +
    
    # basemap
    geom_spatraster(data = rast.focal.2) +
    
    # unit boundaries
    geom_sf(data = clust.units,
            fill = NA,
            color = "white",
            linewidth = 0.2) +
    
    coord_sf(expand = F) +
    
    # theme
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          
          strip.background = element_rect(color = NA,
                                          fill = NA),
          strip.text = element_text(size = 7,
                                    hjust = 0)) +
    
    # colors
    scale_fill_viridis_c() -> plot.focal.2
  
  # cond
  ggplot() +
    
    theme_bw() +
    
    facet_wrap(~ lyr, 
               nrow = 1,
               strip.position = "right") +
    
    # basemap
    geom_spatraster(data = rast.cond) +
    
    # unit boundaries
    geom_sf(data = clust.units,
            fill = NA,
            color = "white",
            linewidth = 0.2) +
    
    coord_sf(expand = F) +
    
    # theme
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          
          strip.background = element_rect(color = NA,
                                          fill = NA),
          strip.text = element_text(size = 7,
                                    hjust = 0)) +
    
    # colors
    scale_fill_viridis_c() -> plot.cond
  
  # plot together
  maps.covar <- plot_grid(plot.focal.1,
                          plot.focal.2,
                          plot.cond,
                          nrow = 3,
                          rel_heights = c(3.05, 0.66, 0.675))
  
  # both big map and covars
  plot_grid(map.base, maps.covar, ncol = 2, rel_widths = c(1.5, 1))
  
}

# ______________________________________________________________________________
# 4b. Save plots ----
# ______________________________________________________________________________

plot_info_panel(1)

# save as high-res png
ggsave("fig_building/example_maps/info_panels/WR.png", 
       dpi = 600, 
       width = 5.98,
       height = 3.96,
       units = "in",
       bg = "white")

plot_info_panel(2)

ggsave("fig_building/example_maps/info_panels/CB.png", 
       dpi = 600, 
       width = 5.98,
       height = 3.96,
       units = "in",
       bg = "white")

plot_info_panel(3)

ggsave("fig_building/example_maps/info_panels/BB.png", 
       dpi = 600, 
       width = 5.98,
       height = 3.96,
       units = "in",
       bg = "white")

plot_info_panel(4)

ggsave("fig_building/example_maps/info_panels/CH.png", 
       dpi = 600, 
       width = 5.98,
       height = 3.96,
       units = "in",
       bg = "white")
