# PROJECT: Habitat selection
# SCRIPT: FIG - Example map
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 Jul 2026
# COMPLETED: 
# LAST MODIFIED: 07 Aug 2026
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
pct.units <- st_read(paste0(spat.dir, "Units/units_fixed_utm/units_fixed_utm.shp"))

# piles and retention
piles <- st_read(paste0(spat.dir, "Ground-truthing/all_piles.shp"))
ret <- st_read(paste0(spat.dir, "Ground-truthing/all_refugia.shp"))

# dissolved polygons
poly.dh <- st_read(paste0(spat.dir, "hand_digitization/dissolved_poly/2a_DhY_dis.shp"))
poly.om <- st_read(paste0(spat.dir, "hand_digitization/dissolved_poly/4_OM_dis.shp"))
poly.msm <- st_read(paste0(spat.dir, "hand_digitization/dissolved_poly/6_MsM_dis.shp"))
poly.jsm <- st_read(paste0(spat.dir, "hand_digitization/dissolved_poly/5_JsM_dis.shp"))

# relocations
pts <- readRDS("D:/hare_project/data_analysis/General/hare-gps-processing-new/data_cleaned/use_background.rds")

# ______________________________________________________________________________
# 2c. Raster data ----
# ______________________________________________________________________________

vo.pre <- rast(paste0(spat.dir, "Rasters/veg_pred/RF/vo_pre_new.tif"))
vo.post <- rast(paste0(spat.dir, "Rasters/veg_pred/RF/vo_post_new.tif"))

rast.all <- rast(paste0(getwd(), "/data_raster/rast_all.tif"))

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

# relocations
pts.1 <- pts |> 
  
  filter(case == 1 &
           site %in% c("1A", "1B", "1C") & 
           season == "off" &
           year %in% c("POST1", "POST2")) |>
  
  st_as_sf(coords = c("x", "y"),
           crs = crs.utm)

# ______________________________________________________________________________
# 5. Bounding box and crop rasters ----
# ______________________________________________________________________________

wr.bbox <- wr.units |> bbox_poly() |> st_buffer(dist = 100) |> bbox_poly()

vo.pre.crop <- vo.pre |> crop(wr.bbox)
vo.post.crop <- vo.post |> crop(wr.bbox)

# transform to UTM and return only those within the bbox
poly.dh <- poly.dh |> st_transform(crs = crs.utm) |> st_intersection(wr.bbox)
poly.om <- poly.om |> st_transform(crs = crs.utm) |> st_intersection(wr.bbox)
poly.msm <- poly.msm |> st_transform(crs = crs.utm) |> st_intersection(wr.bbox)
poly.jsm <- poly.jsm |> st_transform(crs = crs.utm) |> st_intersection(wr.bbox)

# intersect
pts.1 <- pts.1 |> 
  
  # add coordinates
  mutate(x = st_coordinates(pts.1)[, 1],
         y = st_coordinates(pts.1)[, 2]) |>
  
  st_intersection(wr.bbox)

# ______________________________________________________________________________
# 6. Plots ----
# ______________________________________________________________________________
# 6a. Basemap ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  geom_spatraster(data = vo.post.crop) +
  
  # unit boundaries
  geom_sf(data = wr.units,
          fill = NA,
          color = "white",
          linewidth = 0.5) +
  
  # piles
  geom_sf(data = wr.piles,
          color = "white",
          shape = 3,
          size = 0.9) +
  
  # patches
  geom_sf(data = wr.ret,
          fill = NA,
          color = "white",
          linewidth = 0.4) +
  
  # relocations
  #geom_path(data = pts.1,
   #         aes(x = x,
    #            y = y,
     #           group = track_season_post),
      #    color = "black",
       #   alpha = 0.5,
        # linewidth = 0.2) +
  
  # colors
  scale_fill_viridis_c() +
  
  # legend title
  labs(fill = "Visual obstruction") +
  
  # scale bar
  annotation_scale(aes(style = "ticks"),
                   location = "br",
                   line_width = 0.5,
                   height = unit(0.1, "cm"),
                   text_cex = 0.6,
                   pad_x = unit(0.7, "cm"),
                   pad_y = unit(0.5, "cm"),
                   line_col = "white",
                   text_col = "white") +
  
  # theme
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        
        legend.position = "none",
        # legend
        #legend.position = "inside",
        #legend.position.inside = c(0.83, 0.15),
        #legend.direction = "horizontal",
        #legend.background = element_rect(fill = NA),
        #legend.frame = element_rect(color = "white"),
        #legend.text = element_text(color = "white",
        #                           face = "bold"),
        #legend.text.position = "bottom",
        #legend.title = element_text(color = "white",
        #                            hjust = 0.5),
        #legend.title.position = "top",
        
        axis.title = element_blank()) -> map.base

map.base

# ______________________________________________________________________________
# 6b. Doghair (pre-thinned) ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # basemap
  geom_spatraster(data = vo.pre.crop) +
  
  # gray out the basemap
  geom_sf(data = st_difference(st_buffer(wr.bbox, 5), poly.dh),
          alpha = 0.65,
          color = NA) +
  
  geom_sf(data = poly.dh,
          color = "white",
          linewidth = 0.2,
          fill = NA) +
  
  # theme
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  
  # colors
  scale_fill_viridis_c() -> map.dh

# ______________________________________________________________________________
# 6c. Open mature ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # basemap
  geom_spatraster(data = vo.pre.crop) +
  
  # gray out the basemap
  geom_sf(data = st_difference(st_buffer(wr.bbox, 5), poly.om),
          alpha = 0.65,
          color = NA) +
  
  geom_sf(data = poly.om,
          color = "white",
          linewidth = 0.2,
          fill = NA) +
  
  # theme
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  
  # colors
  scale_fill_viridis_c() -> map.om

# ______________________________________________________________________________
# 6d. Dense mature ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # basemap
  geom_spatraster(data = vo.pre.crop) +
  
  # gray out the basemap
  geom_sf(data = st_difference(st_buffer(wr.bbox, 5), poly.msm),
          alpha = 0.65,
          color = NA) +
  
  geom_sf(data = poly.msm,
          color = "white",
          linewidth = 0.2,
          fill = NA) +
  
  # theme
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  
  # colors
  scale_fill_viridis_c() -> map.msm

# ______________________________________________________________________________
# 6e. Jackstraw ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  # basemap
  geom_spatraster(data = vo.pre.crop) +
  
  # gray out the basemap
  geom_sf(data = st_difference(st_buffer(wr.bbox, 5), poly.jsm),
          alpha = 0.65,
          color = NA) +
  
  geom_sf(data = poly.jsm,
          color = "white",
          linewidth = 0.2,
          fill = NA) +
  
  # theme
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  
  # colors
  scale_fill_viridis_c() -> map.jsm

# ______________________________________________________________________________
# 6f. Covariates ----

# I want to stack these so the stand variables are 2 x 2,
# dEdge is right below it and centered,
# and the conditions are below that, 1 x 2

# they should all be standardized and the legend should be shared
# and reflect the full range 

# we should use both the pre- and post- rasters

# ______________________________________________________________________________

# crop to bbox
rast.all.crop <- crop(rast.all, wr.bbox)

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
min(quantile(values(rast.focal.1), 0.01), 
    quantile(values(rast.focal.2), 0.01), 
    quantile(values(rast.cond), 0.01))   # -2.221115

max(quantile(values(rast.focal.1), 0.99), 
    quantile(values(rast.focal.2), 0.99), 
    quantile(values(rast.cond), 0.99))    # 3.806747

# clamp
rast.focal.1 <- clamp(rast.focal.1, -2.221115, 3.806747)
rast.focal.2 <- clamp(rast.focal.2, -2.221115, 3.806747)
rast.cond <- clamp(rast.cond, -2.221115, 3.806747)

# plots
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

# ______________________________________________________________________________
# 7. Plot together ----
# ______________________________________________________________________________

maps.covar <- plot_grid(plot.focal.1,
                        plot.focal.2,
                        plot.cond,
                        nrow = 3,
                        rel_heights = c(3.05, 0.66, 0.675))

# both big map and covars
plot_grid(map.base, maps.covar, ncol = 2, rel_widths = c(1.5, 1))

# save as high-res png
ggsave("fig_building/example_maps/example_fig/WR_example.png", 
       dpi = 600, 
       width = 5.98,
       height = 3.96,
       units = "in",
       bg = "white")



maps.cover <- plot_grid(map.dh, map.msm, map.om, map.jsm, 
                        ncol = 1,
                        align = "h",
                        axis = "l",
                        rel_heights = c(1.5, 1, 1, 1))

maps.cover

# NEXT: Standard plots for all clusters