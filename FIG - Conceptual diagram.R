# PROJECT: Habitat selection
# SCRIPT: FIG - Conceptual diagram
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 03 Aug 2026
# COMPLETED: 
# LAST MODIFIED: 11 Aug 2026
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

# read rasters
ls.1.rast.1 <- rast(paste0(shp.dir, "/ls_1.tif"))
ls.2.rast.1 <- rast(paste0(shp.dir, "/ls_2.tif"))
ls.3.rast.1 <- rast(paste0(shp.dir, "/ls_3.tif"))

# add CRS to unit 2
st_crs(unit.2.sf) <- st_crs(unit.sf)

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
                                      angle = 270),
          legend.key.width = unit(0.15, "cm"))
  
}

# plot together
plot_grid(plot_panel_1(ls.1.rast.1, unit.sf, T),
          plot_grid(plot_panel_1(ls.2.rast.1, unit.2.sf, F),
                    plot_panel_1(ls.3.rast.1, unit.2.sf, F)),
          nrow = 2)

# ______________________________________________________________________________
# 3. PANEL 2 - Telemetry tracks ----

library(ctmm)

n.indiv <- 5

# ______________________________________________________________________________

# define CTMM (we'll use a simple OUF process)
ctmm.mod <- ctmm(
  
  tau = c(3 %#% "hour", 0.1 %#% "hour"),
  isotropic = F,
  range = T,
  error = F,
  sigma = matrix(c(1600,
                   0,
                   0,
                   600),
                 nrow = 2)
  
)
  
# simulation parameters
sim.duration.sec <- 6 * 24 * 60 * 60

# sampling rate
sim.samp.rate <- 1800

# time steps
sim.timestep <- seq(1, sim.duration.sec, sim.samp.rate)

# define HRCs
# we want unaltered to maybe be more concentrated 
rast.ext <- as.numeric(st_bbox(ls.1.rast.1))

un.hrc <- matrix(data = c(runif(n.indiv, rast.ext[1] + 150, rast.ext[3] - 150),
                          runif(n.indiv, rast.ext[2] + 150, rast.ext[4] - 150)),
                 nrow = n.indiv)

alt.hrc <- matrix(data = c(runif(n.indiv, rast.ext[1] + 100, rast.ext[3] - 100),
                           runif(n.indiv, rast.ext[2] + 100, rast.ext[4] - 100)),
                  nrow = n.indiv)

# angles
un.angle <- matrix(data = c(runif(n.indiv, - pi / 2, pi / 2)),
                   nrow = n.indiv)

alt.angle <- matrix(data = c(runif(n.indiv, - pi / 2, pi / 2)),
                   nrow = n.indiv)

plot(un.hrc[ ,1], un.hrc[ ,2])
plot(alt.hrc[ ,1], alt.hrc[ ,2])

# simulate
un.tracks <- data.frame()
alt.tracks <- data.frame()

for (i in 1:n.indiv) {
  
  # UN
  # modify MM
  focal.mod <- ctmm.mod
  focal.mod$mu <- un.hrc[i,]
  focal.mod$sigma@par[3] <- un.angle[i]
  
  # simulate
  un.sim <- simulate(object = focal.mod,
                     t = sim.timestep,
                     complete = T)
    
  un.sim.df <- un.sim@.Data |> 
    
    do.call(cbind, args = _) |>
    
    as.data.frame() |>
    
    dplyr::select(V1, V2, V3) |>
    
    mutate(indiv = i)
  
  # bind in
  un.tracks <- rbind(un.tracks, un.sim.df)
  
  # ALT
  # modify MM
  focal.mod <- ctmm.mod
  focal.mod$mu <- alt.hrc[i, ]
  focal.mod$sigma@par[3] <- alt.angle[i]
  
  # simulate
  alt.sim <- simulate(object = focal.mod,
                     t = sim.timestep,
                     complete = T)
  
  alt.sim.df <- alt.sim@.Data |> 
    
    do.call(cbind, args = _) |>
    
    as.data.frame() |>
    
    dplyr::select(V1, V2, V3) |>
    
    mutate(indiv = i)
  
  # bind in
  alt.tracks <- rbind(alt.tracks, alt.sim.df)
  
}

# column names
colnames(un.tracks) <- c("t", "x", "y", "indiv")
colnames(alt.tracks) <- c("t", "x", "y", "indiv")

# plot function
plot_tracks <- function (.rast, .tracks) {
  
  ggplot() +
  
  theme_minimal() +
    
    geom_spatraster(data = .rast) +
    
    geom_path(data = .tracks,
               aes(x = x,
                   y = y,
                   group = indiv,
                   color = indiv)) + 
    
    coord_sf(expand = F) +
    
    scale_fill_viridis_c() +
    
    scale_color_continuous(palette = c("gray75", "white")) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          axis.title = element_blank())
  
}

# plot together
plot_grid(
  
  plot_tracks(ls.1.rast.1, un.tracks),
  plot_tracks(ls.3.rast.1, alt.tracks),
  nrow = 2
  
)

# HRs and sampled relocations
un.hr <- list()
alt.hr <- list()

un.avail <- data.frame()
alt.avail <- data.frame()

for (i in 1:n.indiv) {
  
  # UN
  focal.relocs <- un.tracks |> filter(indiv == i) |>
    
    st_as_sf(coords = c("x", "y"), 
             crs = "epsg:32611")
  
  # buffer and dissolve
  focal.hr <- st_buffer(focal.relocs, dist = 50) |> st_union() |> st_as_sf()
  
  # sample
  focal.avail <- st_sample(focal.hr, length(sim.timestep) * 10, type = "regular") |> 
    
    st_as_sf() |> mutate(indiv = i)
  
  # add in
  un.hr[[i]] <- focal.hr
  un.avail <- rbind(un.avail, focal.avail)
  
  # ALT
  focal.relocs <- alt.tracks |> filter(indiv == i) |>
    
    st_as_sf(coords = c("x", "y"), 
             crs = "epsg:32611")
  
  # buffer and dissolve
  focal.hr <- st_buffer(focal.relocs, dist = 50) |> st_union() |> st_as_sf()
  
  # sample
  focal.avail <- st_sample(focal.hr, length(sim.timestep) * 10, type = "regular") |> 
    
    st_as_sf() |> mutate(indiv = i)
  
  # add in
  alt.hr[[i]] <- focal.hr
  alt.avail <- rbind(alt.avail, focal.avail)
  
}

# plot
plot_points <- function (.used, .avail) {
  
  ggplot() +
    
    theme_minimal() +
    
    geom_sf(data = unit.bbox,
            fill = NA,
            color = NA) +
    
    # available
    geom_sf(data = .avail,
            aes(group = indiv),
            color = "black",
            size = 0.05,
            alpha = 0.1) + 
    
    # used
    geom_point(data = .used,
               aes(x = x,
                   y = y),
               size = 0.15) + 
    
    coord_sf(expand = F) +
    
    scale_color_continuous(palette = c("gray75", "white")) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          axis.title = element_blank())
  
}

# plot together
plot_grid(
  
  plot_points(un.tracks, un.avail),
  plot_points(alt.tracks, alt.avail),
  nrow = 2
  
)

# full panel
plot_grid(
  
  plot_grid(
    
    plot_tracks(ls.1.rast.1, un.tracks),
    plot_tracks(ls.3.rast.1, alt.tracks),
    nrow = 2
    
  ),
  
  plot_grid(
    
    plot_points(un.tracks, un.avail),
    plot_points(alt.tracks, alt.avail),
    nrow = 2
    
  ),
  
  ncol = 2
  
)
