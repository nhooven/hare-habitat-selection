# PROJECT: Habitat selection
# SCRIPT: FIG - Conceptual diagram
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 03 Aug 2026
# COMPLETED: 
# LAST MODIFIED: 03 Aug 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(cowplot)

# ______________________________________________________________________________
# 2. PANEL 1 - Alteration maps ----
# ______________________________________________________________________________

# read in shapefile
shp.dir <- "D:/hare_project/Maps/Ch 3 conceptual diagram/"

# read in
unit.sf <- st_read(paste0(shp.dir, "unit_polys.shp"))

# create buffered bounding box (square aspect ratio)
unit.bbox <- unit.sf |> 
  
  st_centroid() |>
  
  st_buffer(dist = 450) |> 
  
  spatialEco::bbox_poly()

plot(st_geometry(unit.bbox))
plot(st_geometry(unit.sf), add = T)

# rotate polygon for "altered"
unit.2.sf <- unit.sf |> spatialEco::rotate.polygon(angle = 191)

plot(st_geometry(unit.bbox))
plot(st_geometry(unit.2.sf), add = T)

# generate landscapes - "understory density"
library(NLMR)

ls.1 <- nlm_mosaicfield(ncol = 100,
                        nrow = 100,
                        resolution = 10,
                        mosaic_sd = 1.5)

ls.2 <- nlm_mosaicfield(ncol = 100,
                        nrow = 100,
                        resolution = 10,
                        mosaic_sd = 1.5)

# convert to rast, log-transform, and convert extent
ls.1.rast <- ls.1 |> rast()
ext(ls.1.rast) <- as.vector(st_bbox(unit.bbox))[c(1, 3, 2, 4)]

ls.2.rast <- ls.2 |> rast()
ext(ls.2.rast) <- as.vector(st_bbox(unit.bbox))[c(1, 3, 2, 4)]

# correct CRS
crs(ls.1.rast) <- crs(unit.sf)
crs(ls.2.rast) <- crs(unit.sf)
st_crs(unit.2.sf) <- crs(unit.sf)

# make a dense "unit"
unit.1 <- ls.1.rast |> crop(vect(unit.sf)) |> mask(vect(unit.sf))
unit.2 <- ls.2.rast |> crop(vect(unit.2.sf)) |> mask(vect(unit.2.sf))

# should have more understory than the surrounds
unit.1 <- unit.1 * 2
unit.2 <- unit.2 * 2

ls.1.rast.1 <- merge(unit.1, ls.1.rast)
ls.2.rast.1 <- merge(unit.2, ls.2.rast)

# alter it
unit.pct <- ls.2.rast.1 |> crop(vect(unit.2.sf)) |> mask(vect(unit.2.sf))
unit.pct.1 <- unit.pct * rast(unit.pct, vals = runif(n = length(values(unit.pct)),
                                                     0.2,
                                                     0.8))

# mosaic it
ls.3.rast.1 <- merge(unit.pct.1, ls.2.rast.1)

# standardize
#ls.1.rast.2 <- scale(ls.1.rast.1)
#ls.2.rast.2 <- scale(ls.2.rast.1)
#ls.3.rast.2 <- scale(ls.3.rast.1)

# plot
plot(ls.1.rast.1)
plot(st_geometry(unit.sf), add = T)

plot(ls.2.rast.1)
plot(st_geometry(unit.2.sf), add = T)

plot(ls.3.rast.1)
plot(st_geometry(unit.2.sf), add = T)

plot(c(ls.1.rast.1, ls.2.rast.1, ls.3.rast.1))

# save the rasters
writeRaster(ls.1.rast.1, paste0(shp.dir, "/ls_1.tif"), overwrite = T)
writeRaster(ls.2.rast.1, paste0(shp.dir, "/ls_2.tif"), overwrite = T)
writeRaster(ls.3.rast.1, paste0(shp.dir, "/ls_3.tif"), overwrite = T)

# determine ranges
min(c(min(values(ls.1.rast.1)), min(values(ls.2.rast.1)), min(values(ls.3.rast.1))))
max(c(max(values(ls.1.rast.1)), max(values(ls.2.rast.1)), max(values(ls.3.rast.1))))

# plots
plot_panel_1 <- function (.rast, .sf, .legend = F) {
  
  ggplot() + 
  
  theme_minimal() +
  
  geom_spatraster(data = .rast) +
  
  geom_sf(data = .sf,
          fill = NA,
          color = "white") + 
    
    coord_sf(expand = F) +
    
    scale_fill_viridis_c(name = "Understory",
                         limits = c(0, 1.89)) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          legend.position = ifelse(.legend == T,
                                   "right",
                                   "none"),
          legend.text = element_blank(),
          legend.title.position = "right",
          legend.title = element_text(hjust = 0.5,
                                      angle = 270))
  
}

# plot together
plot_grid(plot_panel_1(ls.1.rast.1, unit.sf, T),
          plot_grid(plot_panel_1(ls.2.rast.1, unit.2.sf, F),
                    plot_panel_1(ls.3.rast.1, unit.2.sf, F)),
          nrow = 2)

  
  
