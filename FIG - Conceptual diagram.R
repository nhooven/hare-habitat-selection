# PROJECT: Habitat selection
# SCRIPT: FIG - Conceptual diagram
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 03 Aug 2026
# COMPLETED: 20 Aug 2026
# LAST MODIFIED: 20 Aug 2026
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
    
    scale_fill_viridis_c(name = "X1",
                         limits = c(0, 1.89)) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          legend.position = ifelse(.legend == T,
                                   "right",
                                   "none"),
          legend.text = element_blank(),
          legend.title.position = "right",
          legend.title = element_text(hjust = 0.5,
                                      angle = 270,
                                      size = 8),
          legend.key.width = unit(0.40, "cm"),
          legend.key.height = unit(0.65, "cm"),
          legend.key.spacing = unit(0.1, "cm"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))
  
}

# plot together
plot_grid(plot_panel_1(ls.1.rast.1, unit.sf, T),
          plot_grid(plot_panel_1(ls.2.rast.1, unit.2.sf, F),
                    plot_panel_1(ls.3.rast.1, unit.2.sf, F)),
          nrow = 2)

# save as high-res png
ggsave("fig_building/conceptual/panel1.png", 
       dpi = 600, 
       width = 4,
       height = 4,
       units = "in",
       bg = "white")

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
# we'll do it manually
rast.ext <- as.numeric(st_bbox(ls.1.rast.1))

un.hrc <- matrix(data = c(rast.ext[1] + 450,
                          rast.ext[1] + 200,
                          rast.ext[1] + 350,
                          rast.ext[1] + 725,
                          rast.ext[1] + 625,
                          
                          rast.ext[2] + 625,
                          rast.ext[2] + 500,
                          rast.ext[2] + 350,
                          rast.ext[2] + 550,
                          rast.ext[2] + 275),
                 nrow = n.indiv)

plot(un.hrc[, 1], un.hrc[, 2])

alt.hrc <- matrix(data = c(rast.ext[1] + 750,
                           rast.ext[1] + 125,
                           rast.ext[1] + 350,
                           rast.ext[1] + 625,
                           rast.ext[1] + 275,
                           
                           rast.ext[2] + 625,
                           rast.ext[2] + 475,
                           rast.ext[2] + 250,
                           rast.ext[2] + 150,
                           rast.ext[2] + 675),
                  nrow = n.indiv)

plot(alt.hrc[, 1], alt.hrc[, 2])

# angles
un.angle <- matrix(data = c(runif(n.indiv, - pi / 2, pi / 2)),
                   nrow = n.indiv)

alt.angle <- matrix(data = c(runif(n.indiv, - pi / 2, pi / 2)),
                   nrow = n.indiv)

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
                   color = indiv),
              linewidth = 0.15) + 
    
    coord_sf(expand = F) +
    
    scale_fill_viridis_c() +
    
    scale_color_continuous(palette = c("gray75", "white")) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          axis.title = element_blank())
  
}

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
  focal.avail <- st_sample(focal.hr, length(sim.timestep), type = "regular") |> 
    
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
  focal.avail <- st_sample(focal.hr, length(sim.timestep), type = "regular") |> 
    
    st_as_sf() |> mutate(indiv = i)
  
  # add in
  alt.hr[[i]] <- focal.hr
  alt.avail <- rbind(alt.avail, focal.avail)
  
}

# bind HRs together
un.hr <- do.call(rbind, un.hr)
alt.hr <- do.call(rbind, alt.hr)

# plot
plot_points <- function (.used, .hr) {
  
  ggplot() +
    
    theme_minimal() +
    
    geom_sf(data = unit.bbox,
            fill = NA,
            color = NA) +
    
    geom_sf(data = .hr,
            fill = "gray75",
            alpha = 0.5,
            color = NA) +
    
    # used
    geom_point(data = .used,
               aes(x = x,
                   y = y),
               size = 0.005,
               alpha = 0.15,
               shape = 21) + 
    
    coord_sf(expand = F) +
    
    scale_color_continuous(palette = c("gray75", "white")) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          axis.title = element_blank())
  
}

# full panel
plot_grid(
  
  plot_grid(
    
  plot_tracks(ls.1.rast.1, un.tracks),
  plot_tracks(ls.3.rast.1, alt.tracks),
  nrow = 2
    
  ),
  
  plot_grid(
    
    plot_points(un.tracks, un.hr),
    plot_points(alt.tracks, alt.hr),
    nrow = 2
    
  ),
  
  nrow = 1
  
)

# save as high-res png
ggsave("fig_building/conceptual/panel2.png", 
       dpi = 600, 
       width = 4,
       height = 4,
       units = "in",
       bg = "white")

# ______________________________________________________________________________
# 4. PANEL 3 - Selection coefficients ----

# y-axis: parameter (1 focal, 2 others)
# x-axis: coefficient

# ______________________________________________________________________________

# create data
pop.data <- data.frame(var = c("X1", "X2", "X3"),
                       group = c("focal", "conditions", "conditions"),
                       mean = c(0.3, -0.35, -0.05),
                       lci = c(0.25, -0.40, -0.25),
                       uci = c(0.35, -0.30, 0.20))

set.seed(105)

ind.data <- data.frame(var = rep(c("X1", "X2", "X3"), each = 10),
                       group = rep(c("focal", "conditions", "conditions"), each = 10),
                       coef = c(rnorm(10, 0.3, 0.1),
                                rnorm(10, -0.35, 0.07),
                                rnorm(10, -0.05, 0.15))) |>
  
  mutate(group = factor(group, levels = c("focal", "conditions")))

# plot 
ggplot() +
  
  theme_classic() +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  geom_errorbar(data = pop.data,
                aes(x = mean,
                    y = var,
                    xmin = lci,
                    xmax = uci,
                    color = group),
                alpha = 0.25,
                height = 0,
                linewidth = 1.5) +
  
  geom_point(data = ind.data,
             aes(x = coef,
                 y = var,
                 color = group),
             alpha = 0.5) +
  
  geom_point(data = pop.data,
             aes(x = mean,
                 y = var,
                 color = group),
             shape = 21,
             fill = "white",
             size = 2.5,
             stroke = 0.8) +
  
  xlab(expression(beta)) +
  ylab("Habitat variable") +
  
  theme(legend.position = "none",
        strip.background = element_rect(color = NA),
        strip.text = element_text(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  
  scale_y_discrete(limits = rev) +
  
  scale_color_manual(values = c("gray", "darkblue"))

ggsave("fig_building/conceptual/panel3.png", 
       dpi = 600, 
       width = 4,
       height = 4,
       units = "in",
       bg = "white")

# ______________________________________________________________________________
# 5. PANEL 4 - Functional response predictions ----

# 2 x 2

# ______________________________________________________________________________
# 5a. H1 ----
# ______________________________________________________________________________

set.seed(427)

h1.data <- data.frame(coef = ind.data$coef[1:10],
                      avail = runif(10, 0.55, 0.85))

ggplot() +
  
  theme_classic() +
  
  geom_point(data = h1.data,
             aes(x = avail,
                 y = coef),
             shape = 21) +
  
  ylab(expression(beta)) +
  xlab("X1 availability") +
  
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  
  ggtitle("H1") -> h1.plot

# ______________________________________________________________________________
# 5b. H2 ----
# ______________________________________________________________________________

set.seed(334)

h2.data <- data.frame(coef = ind.data$coef[1:10],
                      avail = ((ind.data$coef[1:10] - 0.3) / 0.25) + rnorm(10, 0, 0.5))

ggplot() +
  
  theme_classic() +
  
  geom_point(data = h2.data,
             aes(x = avail,
                 y = coef),
             shape = 21) +
  
  geom_smooth(data = h2.data,
              aes(x = avail,
                  y = coef),
              method = "lm",
              color ="black") +
  
  ylab(expression(beta)) +
  xlab("X1 availability") +
  
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  
  ggtitle("H2") -> h2.plot

# ______________________________________________________________________________
# 5c. H3 ----
# ______________________________________________________________________________

h3.data <- data.frame(coef = ind.data$coef[1:10]) |>
  
  mutate(trt = c("UN", "UN", "ALT", "UN", "ALT", 
                 "ALT", "UN", "UN", "ALT", "ALT"))

h3.data.summ <- h3.data |>  

  group_by(trt) |>
  
  summarize(mean = mean(coef),
            se = sd(coef) / sqrt(n())) |>
  
  mutate(lci = mean - 1.645 * se,
         uci = mean + 1.645 * se)

ggplot() +
  
  theme_classic() +
  
  geom_point(data = h3.data,
             aes(x = trt,
                 y = coef,
                 color = trt,
                 shape = trt)) +
  
  geom_errorbar(data = h3.data.summ,
                aes(x = trt,
                    y = mean,
                    ymin = lci,
                    ymax = uci),
                width = 0,
                alpha = 0.25,
                linewidth = 1.1) +
  
  geom_point(data = h3.data.summ,
              aes(x = trt,
                  y = mean,
                  color = trt,
                  shape = trt),
             fill = "white",
             size = 2.2) +
  
  ylab(expression(beta)) +
  xlab("Treatment") +
  
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  
  scale_color_manual(values = c("darkorange3", "green4")) +
  scale_shape_manual(values = c(23, 21)) +
  
  scale_x_discrete(limits = rev) +
  
  ggtitle("H3") -> h3.plot

# ______________________________________________________________________________
# 5d. H4 ----
# ______________________________________________________________________________

h4.data <- data.frame(avail = runif(10, -2, 2),
                      trt = rep(c("UN", "ALT"), each = 5)) |>
  
  mutate(trt = factor(trt, levels = c("UN", "ALT")))

h4.data$coef <- c(rnorm(5, 0.3, 0.2),
                  0.2 + h4.data$avail[6:10] * 0.45 + rnorm(5, 0, 0.25))

ggplot() +
  
  theme_classic() +
  
  geom_point(data = h4.data,
             aes(x = avail,
                 y = coef,
                 color = trt,
                 shape = trt)) +
  
  geom_smooth(data = h4.data,
              aes(x = avail,
                  y = coef,
                  color = trt,
                  fill = trt,
                  linetype = trt),
              method = "lm",
              alpha = 0.15) +
  
  ylab(expression(beta)) +
  xlab("X1 availability") +
  
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +

  scale_color_manual(values = c("green4", "darkorange3")) +
  scale_fill_manual(values = c("green4", "darkorange3")) +
  scale_shape_manual(values = c(21, 23)) +
  
  ggtitle("H4") -> h4.plot

# ______________________________________________________________________________
# 5e. Plot together ----
# ______________________________________________________________________________

plot_grid(h1.plot, h2.plot, h3.plot, h4.plot, nrow = 2)

# save as high-res png
ggsave("fig_building/conceptual/panel4.png", 
       dpi = 600, 
       width = 4,
       height = 4,
       units = "in",
       bg = "white")
