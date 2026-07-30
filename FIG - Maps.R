# PROJECT: Habitat selection
# SCRIPT: FIG - Maps
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 30 Jul 2026
# COMPLETED: 
# LAST MODIFIED: 30 Jul 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(mgcv)
library(terra)
library(tidyterra)
library(patchwork)

# ______________________________________________________________________________
# 2. Read in models and data ----
# ______________________________________________________________________________

# HSF results
M.off <- readRDS("model_results/M_off.rds")
M.on <- readRDS("model_results/M_on.rds")

# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.cc <- readRDS("model_results/fr_models/off_cc.rds")
off.dEdge <- readRDS("model_results/fr_models/off_dEdge.rds")

# on
on.cc <- readRDS("model_results/fr_models/on_cc.rds")
on.dEdge <- readRDS("model_results/fr_models/on_dEdge.rds")

# means and SDs
mean.sd.off <- readRDS("data_for_model/mean_sd_off.rds")
mean.sd.on <- readRDS("data_for_model/mean_sd_on.rds")

# unit shapefile
units <- st_read("D:/hare_project/data_spatial/Units/units_fixed_utm/units_fixed_utm.shp") |>
  
  dplyr::select(name, geometry) |>
  
  arrange(name)

# ______________________________________________________________________________
# 4. Process raster stack ----
# ______________________________________________________________________________

# read
rast.all <- rast("data_raster/rast_all.tif")

# raster directory
dir.rast <- "D:/hare_project/data_spatial/Rasters/"

# functions
# standardize a raster layer function
stand_rast <- function (.layer, .var, .season = "off") {
  
  if (.season == "off") {
    
    layer.stand <- (.layer - mean.sd.off$mean[mean.sd.off$name == .var]) / mean.sd.off$sd[mean.sd.off$name == .var]
    
  } else {
    
    layer.stand <- (.layer - mean.sd.on$mean[mean.sd.on$name == .var]) / mean.sd.on$sd[mean.sd.on$name == .var]
    
  }
  
  return(layer.stand)
  
}

# prepare raster (availabilities and standardized)
prep_rast <- function (.site,
                       .radius = 300,          
                       .season = "off",
                       .year = "pre") {
  
  # which vars?
  const.vars <- c("dEdge", "twi", "vrm")
  
  if (.season == "off" & .year == "pre") { which.vars <- c("vo.pre", "ch.pre", "cc.pre", const.vars) }
  if (.season == "off" & .year == "post") { which.vars <- c("vo.post", "ch.post", "cc.post", const.vars) }
  if (.season == "on" & .year == "pre") { which.vars <- c("stem.pre", "ch.pre", "cc.pre", const.vars) }
  if (.season == "on" & .year == "post") { which.vars <- c("stem.post", "ch.post", "cc.post", const.vars) }
  
  # variable names
  if (.season == "off") { var.names <- c("vo", "ch", "cc", const.vars) }
  if (.season == "on") { var.names <- c("stem", "ch", "cc", const.vars) }
  
  # subset site
  focal.site <- units |> filter(name == .site)
  
  # extract centroid
  focal.site.cent <- st_centroid(focal.site)
  
  # buffered site bounding box
  # this differs from the script 08d function
  # we want all the site-specific maps to have the same dimensions
  # we'll need to tune this so it fits the largest on
  focal.bbox <- st_buffer(focal.site.cent, dist = 525) |>
    
    # create a bbox
    spatialEco::bbox_poly()
  
  # test plots
  #plot(st_geometry(focal.bbox))
  #plot(st_geometry(focal.site), add = T)
  
  # subset and crop rasts
  rast.crop <- rast.all |> subset(subset = which.vars) |> crop(vect(focal.bbox))
  
  # moving window analysis
  focal.mat <- focalMat(rast.crop, d = .radius, type = "circle") 
  focal.mat[focal.mat > 0] <- 1
  
  # apply over all layers
  rast.avail <- sapp(rast.crop, focal, w = focal.mat, "mean", na.rm = T)
  
  # names
  names(rast.avail) <- paste0("a.", var.names)
  names(rast.crop) <- var.names
  
  # standardize and bind together
  # AVAILABILITIES SHOULD NOT BE STANDARDIZED
  if (.season == "off") { 
    
    rast.stand.1 <- stand_rast(rast.crop$vo, "vo") 
    rast.stand.a1 <- rast.avail$a.vo 
    
  } else {
    
    rast.stand.1 <- stand_rast(rast.crop$stem, "stem")
    rast.stand.a1 <- rast.avail$a.stem
    
  }
  
  rast.stand <- c(
    
    # covariates
    rast.stand.1,
    stand_rast(rast.crop$ch, "ch", .season),
    stand_rast(rast.crop$cc, "cc", .season),
    stand_rast(rast.crop$dEdge, "dEdge", .season),
    stand_rast(rast.crop$twi, "twi", .season),
    stand_rast(rast.crop$twi^2, "twi2", .season),
    stand_rast(rast.crop$vrm, "vrm", .season),
    stand_rast(rast.crop$vrm^2, "vrm2", .season),
    
    # availabilities
    rast.stand.a1
    
  )
  
  # change names for squared terms
  names(rast.stand)[c(6, 8)] <- c("twi2", "vrm2")
  
  # return
  return(rast.stand)
  
}

# ______________________________________________________________________________
# 5. Calculate predictions by site ----
# ______________________________________________________________________________
# 5a. Mean only ----
# ______________________________________________________________________________

pred_hsf <- function (.site,
                      .season,
                      .year) {
  
  # ____________________________________________________________________________
  # ALWAYS calculate the "pre" raster to prevent artifacts
  
  # importantly, the TRT var should be the correct post one! otherwise the 
  # functional models assume unthinned (not what we want)
  
  # so, raster pre, model predictions post
  
  # ____________________________________________________________________________
  
  site.rast.pre <- prep_rast(.site = .site, .season = .season, .year = "pre")
  
  # add TRT for predictions
  TRT <- case_when(
    
    .year == "pre" ~ "UNTHIN",
    .year == "post" & .site %in% c("1C", "2C", "3C", "4C") ~ "UNTHIN",
    .year == "post" & .site %in% c("1A", "2B", "3B", "4A") ~ "RET",
    .year == "post" & .site %in% c("1B", "2A", "3A", "4B") ~ "PIL"
    
  )
  
  # add a TRT and clust raster
  site.rast.pre <- c(site.rast.pre,
                     rast(site.rast.pre, nlyrs = 1, names = "TRT", vals = TRT),
                     rast(site.rast.pre, nlyrs = 1, names = "cluster", vals = substr(.site, 1, 1)))
  
  if (.season == "off") {
    
    # base model
    hsf <- M.off[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.vo <- subset(site.rast.pre, c("a.vo", "TRT", "cluster"))
    
    # names
    names(rast.vo)[1] <- "avail"
    
    # FR predictions
    beta.vo <- predict(object = rast.vo, 
                       model = off.vo, 
                       fun = predict.gam, 
                       na.omit = T,
                       newdata.guaranteed = TRUE)
    
    beta.cc <- predict(object = rast.vo, 
                       model = off.cc, 
                       fun = predict.gam, 
                       na.omit = T,
                       newdata.guaranteed = TRUE)
    
    beta.dEdge <- predict(object = rast.vo, 
                          model = off.dEdge, 
                          fun = predict.gam, 
                          na.omit = T,
                          newdata.guaranteed = TRUE)
    
    # main coefs
    beta.ch <- hsf$mean[hsf$param == "ch"]
    beta.twi <- hsf$mean[hsf$param == "twi"]
    beta.twi2 <- hsf$mean[hsf$param == "twi2"]
    beta.vrm <- hsf$mean[hsf$param == "vrm"]
    beta.vrm2 <- hsf$mean[hsf$param == "vrm2"]
    
    # calculate log RSS prediction
    log.rss <- 
      
      # base
      beta.twi * site.rast.pre$twi +
      beta.twi2 * site.rast.pre$twi2 +
      beta.vrm * site.rast.pre$vrm +
      beta.vrm2 * site.rast.pre$vrm2 +
      beta.ch * site.rast.pre$ch +
      
      # functional responses
      beta.vo * site.rast.pre$vo +
      beta.cc * site.rast.pre$cc +
      beta.dEdge * site.rast.pre$dEdge
    
  } # season == "off"
  
  if (.season == "on") {
    
    # base model
    hsf <- M.on[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.stem <- subset(site.rast.pre, c("a.stem", "TRT", "cluster"))
    
    # names
    names(rast.stem)[1] <- "avail"
    
    # FR predictions
    beta.cc <- predict(object = rast.stem, 
                       model = on.cc, 
                       fun = predict.gam,
                       na.omit = T, 
                       newdata.guaranteed = TRUE)
    
    beta.dEdge <- predict(object = rast.stem, 
                          model = on.dEdge, 
                          fun = predict.gam,
                          na.omit = T,
                          newdata.guaranteed = TRUE)
    
    # main coefs
    beta.stem <- hsf$mean[hsf$param == "stem"]
    beta.ch <- hsf$mean[hsf$param == "ch"]
    beta.twi <- hsf$mean[hsf$param == "twi"]
    beta.twi2 <- hsf$mean[hsf$param == "twi2"]
    beta.vrm <- hsf$mean[hsf$param == "vrm"]
    beta.vrm2 <- hsf$mean[hsf$param == "vrm2"]
    
    # calculate log RSS prediction
    log.rss <- 
      
      # base
      beta.twi * site.rast.pre$twi +
      beta.twi2 * site.rast.pre$twi2 +
      beta.vrm * site.rast.pre$vrm +
      beta.vrm2 * site.rast.pre$vrm2 +
      beta.stem * site.rast.pre$stem +
      beta.ch * site.rast.pre$ch +
      
      # functional responses
      beta.cc * site.rast.pre$cc +
      beta.dEdge * site.rast.pre$dEdge
    
  } # season == "on"
  
  log.rss.pre <- log.rss
  
  # ____________________________________________________________________________
  
  # prepare rasters
  site.rast <- prep_rast(.site = .site, .season = .season, .year = .year)
  
  # add TRT for predictions
  TRT <- case_when(
    
    .year == "pre" ~ "UNTHIN",
    .year == "post" & .site %in% c("1C", "2C", "3C", "4C") ~ "UNTHIN",
    .year == "post" & .site %in% c("1A", "2B", "3B", "4A") ~ "RET",
    .year == "post" & .site %in% c("1B", "2A", "3A", "4B") ~ "PIL"
    
  )
  
  # add a TRT and clust raster
  site.rast <- c(site.rast,
                 rast(site.rast, nlyrs = 1, names = "TRT", vals = TRT),
                 rast(site.rast, nlyrs = 1, names = "cluster", vals = substr(.site, 1, 1)))
  
  if (.season == "off") {
    
    # base model
    hsf <- M.off[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.vo <- subset(site.rast, c("a.vo", "TRT", "cluster"))
    
    # names
    names(rast.vo)[1] <- "avail"
    
    # FR predictions
    beta.vo <- predict(object = rast.vo, 
                       model = off.vo, 
                       fun = predict.gam, 
                       na.omit = T,
                       newdata.guaranteed = TRUE)
    
    beta.cc <- predict(object = rast.vo, 
                       model = off.cc, 
                       fun = predict.gam, 
                       na.omit = T,
                       newdata.guaranteed = TRUE)
    
    beta.dEdge <- predict(object = rast.vo, 
                          model = off.dEdge, 
                          fun = predict.gam, 
                          na.omit = T,
                          newdata.guaranteed = TRUE)
    
    # main coefs
    beta.ch <- hsf$mean[hsf$param == "ch"]
    beta.twi <- hsf$mean[hsf$param == "twi"]
    beta.twi2 <- hsf$mean[hsf$param == "twi2"]
    beta.vrm <- hsf$mean[hsf$param == "vrm"]
    beta.vrm2 <- hsf$mean[hsf$param == "vrm2"]
    
    # calculate log RSS prediction
    log.rss <- 
      
      # base
      beta.twi * site.rast$twi +
      beta.twi2 * site.rast$twi2 +
      beta.vrm * site.rast$vrm +
      beta.vrm2 * site.rast$vrm2 +
      beta.ch * site.rast$ch +
      
      # functional responses
      beta.vo * site.rast$vo +
      beta.cc * site.rast$cc +
      beta.dEdge * site.rast$dEdge
    
  } # season == "off"
  
  if (.season == "on") {
    
    # base model
    hsf <- M.on[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.stem <- subset(site.rast, c("a.stem", "TRT", "cluster"))
    
    # names
    names(rast.stem)[1] <- "avail"
    
    # FR predictions
    beta.cc <- predict(object = rast.stem, 
                       model = on.cc, 
                       fun = predict.gam,
                       na.omit = T, 
                       newdata.guaranteed = TRUE)
    
    beta.dEdge <- predict(object = rast.stem, 
                          model = on.dEdge, 
                          fun = predict.gam,
                          na.omit = T,
                          newdata.guaranteed = TRUE)
    
    # main coefs
    beta.stem <- hsf$mean[hsf$param == "stem"]
    beta.ch <- hsf$mean[hsf$param == "ch"]
    beta.twi <- hsf$mean[hsf$param == "twi"]
    beta.twi2 <- hsf$mean[hsf$param == "twi2"]
    beta.vrm <- hsf$mean[hsf$param == "vrm"]
    beta.vrm2 <- hsf$mean[hsf$param == "vrm2"]
    
    # calculate log RSS prediction
    log.rss <- 
      
      # base
      beta.twi * site.rast$twi +
      beta.twi2 * site.rast$twi2 +
      beta.vrm * site.rast$vrm +
      beta.vrm2 * site.rast$vrm2 +
      beta.stem * site.rast$stem +
      beta.ch * site.rast$ch +
      
      # functional responses
      beta.cc * site.rast$cc +
      beta.dEdge * site.rast$dEdge
    
  } # season == "on"
  
  # mosaic in to prevent artifacts
  rss.pre <- exp(log.rss.pre)
  rss <- exp(log.rss)
  
  # slightly buffered unit 
  focal.unit <- units |> filter(name == .site) |> st_buffer(dist = 30)
  
  # mask
  rss.mask <- rss |> crop(focal.unit) |> mask(focal.unit)
  
  out.rast <- merge(rss.mask, rss.pre)
  
  return(out.rast)
  
} # f()

# ______________________________________________________________________________
# 7. Process rasters ----

# test each unit
plot(pred_hsf("2C", "off", "pre"))

# 525 m away from the unit centroid works

# ______________________________________________________________________________
# 7a. 1 - W Rabbit ----
# ______________________________________________________________________________

# 1A
off.pre.1A <- pred_hsf("1A", "off", "pre")
off.post.1A <- pred_hsf("1A", "off", "post")

on.pre.1A <- pred_hsf("1A", "on", "pre")
on.post.1A <- pred_hsf("1A", "on", "post")

# change
# this will be difference between the log(intensities)
off.1A.change <- log(off.post.1A) - log(off.pre.1A)
on.1A.change <- log(on.post.1A) - log(on.pre.1A)

# 1B
off.pre.1B <- pred_hsf("1B", "off", "pre")
off.post.1B <- pred_hsf("1B", "off", "post")

on.pre.1B <- pred_hsf("1B", "on", "pre")
on.post.1B <- pred_hsf("1B", "on", "post")

# change
off.1B.change <- log(off.post.1B) - log(off.pre.1B)
on.1B.change <- log(on.post.1B) - log(on.pre.1B)

# 1C (pre is fine)
off.1C <- pred_hsf("1C", "off", "pre")
on.1C <- pred_hsf("1C", "on", "pre")

# ______________________________________________________________________________
# 8. Build maps ----

# relocations 
pts <- readRDS("D:/hare_project/data_analysis/General/hare-gps-processing-new/data_cleaned/use_background.rds") |>
  
  filter(case == 1) |>
  
  st_as_sf(coords = c("x", "y"),
           crs = "epsg:32611")

# map theme
theme_hsf <- function () {
  
  theme_bw() +
    
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          
          legend.title = element_blank())
  
}

# functions
map_unthinned <- function (.unit, .season) {
  
  # choose correct raster
  if (.unit == "1C" & .season == "off") { .rast <- off.1C }
  if (.unit == "1C" & .season == "on") { .rast <- on.1C }
  if (.unit == "2C" & .season == "off") { .rast <- off.2C }
  if (.unit == "2C" & .season == "on") { .rast <- on.2C }
  if (.unit == "3C" & .season == "off") { .rast <- off.3C }
  if (.unit == "3C" & .season == "on") { .rast <- on.3C }
  if (.unit == "4C" & .season == "off") { .rast <- off.4C }
  if (.unit == "4C" & .season == "on") { .rast <- on.4C }
  
  # subset unit
  unit.poly <- units |> filter(name == .unit)
  
  # plot
  ggplot() +
    
    theme_hsf() +
    
    # intensity prediction
    geom_spatraster(data = .rast) +
    
    # unit boundary
    geom_sf(data = unit.poly,
            fill = NA,
            color = "white") +
    
    # correct viridis fill
    scale_fill_viridis_c(option = ifelse(.season == "off",
                                         "inferno",
                                         "mako"))
  
}

map_thinned <- function (.unit, .season) {
  
  # choose correct raster
  if (.unit == "1A" & .season == "off") { 
    
    .rast1 <- off.pre.1A
    .rast2 <- off.post.1A
    .rast3 <- off.1A.change 
    
    }
  
  if (.unit == "1A" & .season == "on") { 
    
    .rast1 <- on.pre.1A
    .rast2 <- on.post.1A 
    .rast3 <- on.1A.change 
    
  }
  
  if (.unit == "1B" & .season == "off") { 
    
    .rast1 <- off.pre.1B
    .rast2 <- off.post.1B
    .rast3 <- off.1B.change 
    
  }
  
  if (.unit == "1B" & .season == "on") { 
    
    .rast1 <- on.pre.1B
    .rast2 <- on.post.1B 
    .rast3 <- on.1B.change
    
  }
  
  # add other units as needed
  
  # stack intensity rasters
  rast.int <- c(.rast1, .rast2)
  
  # names
  names(rast.int) <- c("pre", "post")

  # subset unit
  unit.poly <- units |> filter(name == .unit)
  
  # pre-post plot
  ggplot() +
    
    theme_hsf() +
    
    facet_wrap(~ lyr) +
    
    # intensity prediction
    geom_spatraster(data = rast.int) +
    
    # unit boundary
    geom_sf(data = unit.poly,
            fill = NA,
            color = "white") +
    
    # correct viridis fill
    scale_fill_viridis_c(option = ifelse(.season == "off",
                                         "inferno",
                                         "mako")) +
    
    labs(fill = "intensity") +
    
    # legend position
    theme(legend.position = "bottom",
          legend.title.position = "top") +
    
    # strips
    theme(strip.background = element_blank(),
          strip.text = element_text(hjust = 0.05)) -> pp.plot
  
  # log-change plot
  ggplot() +
    
    theme_hsf() +
    
    # intensity prediction
    geom_spatraster(data = .rast3) +
    
    # unit boundary
    geom_sf(data = unit.poly,
            fill = NA,
            color = "black") +
    
    # diverging palette (continuous from colorbrewer)
    scale_fill_gradient2(low = "#8c510a",
                         mid = "white",
                         high = "#01665e") +
    
    # title
    theme(legend.position = "bottom") -> lc.plot
  
  # patchwork plots together
  pp.plot + lc.plot + plot_layout(widths = c(2, 1))
  
}

# ______________________________________________________________________________
# 8a. W Rabbit ----
# ______________________________________________________________________________

# 1C
map_unthinned("1C", "off")
map_unthinned("1C", "on")

# 1A
map_thinned("1A", "off")
map_thinned("1A", "on")

# 1B
map_thinned("1B", "off")
map_thinned("1B", "on")
