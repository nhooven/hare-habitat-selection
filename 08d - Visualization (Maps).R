# PROJECT: Habitat selection
# SCRIPT: 08d - Visualization (Maps)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 08 Jun 2026
# COMPLETED: 18 Jun 2026
# LAST MODIFIED: 18 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(sf)
library(mgcv)
library(terra)
library(tidyterra)

# ______________________________________________________________________________
# 2. Read in models and data ----
# ______________________________________________________________________________

# HSF results
M.off <- readRDS("model_results/M_off.rds")
M.on <- readRDS("model_results/M_on.rds")

# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.ch <- readRDS("model_results/fr_models/off_ch.rds")

# on
on.stem <- readRDS("model_results/fr_models/on_stem.rds")
on.ch <- readRDS("model_results/fr_models/on_ch.rds")
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
                       .radius = 150,          # ~ mean 99% HR radius
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
  
  # buffered site bounding box
  focal.bbox <- spatialEco::bbox_poly(focal.site) |>
    
    st_buffer(dist = 250) |>
    
    spatialEco::bbox_poly()
  
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
    stand_rast(rast.crop$cc^2, "cc2", .season),
    stand_rast(rast.crop$dEdge, "dEdge", .season),
    stand_rast(rast.crop$twi, "twi", .season),
    stand_rast(rast.crop$twi^2, "twi2", .season),
    stand_rast(rast.crop$vrm, "vrm", .season),
    stand_rast(rast.crop$vrm^2, "vrm2", .season),
    
    # availabilities
    rast.stand.a1
    
  )
  
  # change names for squared terms
  names(rast.stand)[c(4, 7, 9)] <- c("cc2", "twi2", "vrm2")
  
  # return
  return(rast.stand)
  
}

# use function
rast.test <- prep_rast("1A", .season = "off", .year = "pre")

# ______________________________________________________________________________
# 5. Calculate predictions by site ----
# ______________________________________________________________________________
# 5a. Mean only ----
# ______________________________________________________________________________

pred_hsf <- function (.site,
                      .season,
                      .year,
                      .log = F) {
  
  # prepare raster
  site.rast <- prep_rast(.site = .site, .season = .season, .year = .year)
  
  # add TRT for predictions
  TRT <- case_when(
    
    .year == "pre" ~ "UNTHIN",
    .year == "post" & .site %in% c("1C", "2C", "3C", "4C") ~ "UNTHIN",
    .year == "post" & .site %in% c("1A", "2B", "3B", "4A") ~ "RET",
    .year == "post" & .site %in% c("1B", "2A", "3A", "4B") ~ "PIL"
    
  )
  
  # add a TRT raster
  site.rast <- c(site.rast,
                 rast(site.rast, nlyrs = 1, names = "TRT", vals = TRT))
  
  if (.season == "off") {
    
    # base model
    hsf <- M.off[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.vo <- subset(site.rast, c("a.vo", "TRT"))
    
    # names
    names(rast.vo)[1] <- "avail"
    
    # FR predictions
    beta.vo <- predict(object = rast.vo, 
                       model = off.vo, 
                       fun = predict.gam, 
                       na.omit = T,
                       exclude = "s(cluster)", 
                       newdata.guaranteed = TRUE)
    
    beta.ch <- predict(object = rast.vo, 
                       model = off.ch, 
                       fun = predict.gam, 
                       na.omit = T,
                       exclude = "s(cluster)", 
                       newdata.guaranteed = TRUE)
    
    # main coefs
    beta.cc <- hsf$mean[hsf$param == "cc"]
    beta.cc2 <- hsf$mean[hsf$param == "cc2"]
    beta.twi <- hsf$mean[hsf$param == "twi"]
    beta.twi2 <- hsf$mean[hsf$param == "twi2"]
    beta.vrm <- hsf$mean[hsf$param == "vrm"]
    beta.vrm2 <- hsf$mean[hsf$param == "vrm2"]
    beta.dEdge <- hsf$mean[hsf$param == "dEdge"]
      
    # calculate log RSS prediction
    log.rss <- 
        
      # base
      beta.cc * site.rast$cc +
      beta.cc2 * site.rast$cc2 +
      beta.twi * site.rast$twi +
      beta.twi2 * site.rast$twi2 +
      beta.vrm * site.rast$vrm +
      beta.vrm2 * site.rast$vrm2 +
      beta.dEdge * site.rast$dEdge +
      
      # functional responses
      beta.vo * site.rast$vo +
      beta.ch * site.rast$ch
    
  } # season == "off"
  
  if (.season == "on") {
    
    # base model
    hsf <- M.on[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.stem <- subset(site.rast, c("a.stem", "TRT"))
    
    # names
    names(rast.stem)[1] <- "avail"
    
    # FR predictions
    beta.stem <- predict(object = rast.stem, 
                         model = on.stem, 
                         fun = predict.gam,
                         na.omit = T,
                         exclude = "s(cluster)", 
                         newdata.guaranteed = TRUE)
    
    beta.ch <- predict(object = rast.stem, 
                       model = on.ch, 
                       fun = predict.gam,
                       na.omit = T,
                       exclude = "s(cluster)", 
                       newdata.guaranteed = TRUE)
    
    beta.dEdge <- predict(object = rast.stem, 
                          model = on.dEdge, 
                          fun = predict.gam,
                          na.omit = T,
                          exclude = "s(cluster)", 
                          newdata.guaranteed = TRUE)
    
    # main coefs
    beta.cc <- hsf$mean[hsf$param == "cc"]
    beta.cc2 <- hsf$mean[hsf$param == "cc2"]
    beta.twi <- hsf$mean[hsf$param == "twi"]
    beta.twi2 <- hsf$mean[hsf$param == "twi2"]
    beta.vrm <- hsf$mean[hsf$param == "vrm"]
    beta.vrm2 <- hsf$mean[hsf$param == "vrm2"]
      
    # calculate log RSS prediction
    log.rss <- 
        
      # base
      beta.cc * site.rast$cc +
      beta.cc2 * site.rast$cc2 +
      beta.twi * site.rast$twi +
      beta.twi2 * site.rast$twi2 +
      beta.vrm * site.rast$vrm +
      beta.vrm2 * site.rast$vrm2 +
      
      # functional responses
      beta.stem * site.rast$stem +
      beta.ch * site.rast$ch +
      beta.dEdge * site.rast$dEdge
    
  } # season == "on"
  
  # log or exp?
  if (.log == TRUE) { out.rast <- log.rss } else { out.rast <- exp(log.rss) }
  
  return(out.rast)
  
}

# ______________________________________________________________________________
# 5b. Mean + SD ----

# returns a list of rasters:
  # [[1]]: mean predictions
  # [[2]]: SD predictions (from Monte Carlo resampling)

# ______________________________________________________________________________

pred_hsf_MC <- function (.site,
                         .season,
                         .year,
                         .log = F) {
  
  # prepare raster
  site.rast <- prep_rast(.site = .site, .season = .season, .year = .year)
  
  # add TRT for predictions
  TRT <- case_when(
    
    .year == "pre" ~ "UNTHIN",
    .year == "post" & .site %in% c("1C", "2C", "3C", "4C") ~ "UNTHIN",
    .year == "post" & .site %in% c("1A", "2B", "3B", "4A") ~ "RET",
    .year == "post" & .site %in% c("1B", "2A", "3A", "4B") ~ "PIL"
    
  )
  
  # add a TRT raster
  site.rast <- c(site.rast,
                 rast(site.rast, nlyrs = 1, names = "TRT", vals = TRT))
  
  # helper functions
  # sample from availability-weighted, spatially-varying betas
  sample.betas <- function (.rast) {
    
    rast.i <- rast(.rast$fit, vals = rnorm(length(values(.rast$fit)),
                                           values(.rast$fit),
                                           values(.rast$se.fit)))
    
    return(rast.i)
    
  }
  
  if (.season == "off") {
    
    # base model
    hsf <- M.off[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.vo <- subset(site.rast, c("a.vo", "TRT"))
    rast.dOM <- subset(site.rast, c("a.dOM", "TRT"))
    rast.dDM <- subset(site.rast, c("a.dDM", "TRT"))
    rast.shdi <- subset(site.rast, c("a.shdi", "TRT"))
    
    # names
    names(rast.vo)[1] <- "avail"
    names(rast.dOM)[1] <- "avail"
    names(rast.dDM)[1] <- "avail"
    names(rast.shdi)[1] <- "avail"
    
    # FR predictions
    beta.vo <- predict(object = rast.vo, model = off.vo, fun = predict.gam, se.fit = T, na.omit = T)
    beta.dOM <- predict(object = rast.dOM, model = off.dOM, fun = predict.gam, se.fit = T)
    beta.dDM <- predict(object = rast.dDM, model = off.dDM, fun = predict.gam, se.fit = T)
    beta.shdi <- predict(object = rast.shdi, model = off.shdi, fun = predict.gam, se.fit = T)
    
    # Monte Carlo predictions
    n.samp <- 100
    all.log.rss <- rast()
    
    for (i in 1:n.samp) {
      
      # sample from availability-weighted, spatially-varying betas
      beta.vo.i <- sample.betas(beta.vo)
      beta.dOM.i <- sample.betas(beta.dOM)
      beta.dDM.i <- sample.betas(beta.dDM)
      beta.shdi.i <- sample.betas(beta.shdi)
      
      # sample other coefs (non-spatially varying)
      beta.ch.i <- rnorm(1, hsf$mean[hsf$param == "ch"], hsf$sd[hsf$param == "ch"])
      beta.cc.i <- rnorm(1, hsf$mean[hsf$param == "cc"], hsf$sd[hsf$param == "cc"])
      beta.cc2.i <- rnorm(1, hsf$mean[hsf$param == "cc2"], hsf$sd[hsf$param == "cc2"])
      beta.twi.i <- rnorm(1, hsf$mean[hsf$param == "twi"], hsf$sd[hsf$param == "twi"])
      beta.twi2.i <- rnorm(1, hsf$mean[hsf$param == "twi2"], hsf$sd[hsf$param == "twi2"])
      beta.vrm.i <- rnorm(1, hsf$mean[hsf$param == "vrm"], hsf$sd[hsf$param == "vrm"])
      beta.vrm2.i <- rnorm(1, hsf$mean[hsf$param == "vrm2"], hsf$sd[hsf$param == "vrm2"])
      beta.ed.i <- rnorm(1, hsf$mean[hsf$param == "ed"], hsf$sd[hsf$param == "ed"])
      
      # calculate log RSS prediction
      log.rss <- 
        
        # base
        beta.ch.i * site.rast$ch +
        beta.cc.i * site.rast$cc +
        beta.cc2.i * site.rast$cc2 +
        beta.twi.i * site.rast$twi +
        beta.twi2.i * site.rast$twi2 +
        beta.vrm.i * site.rast$vrm +
        beta.vrm2.i * site.rast$vrm2 +
        beta.ed.i * site.rast$ed +
        
        # functional responses
        beta.vo.i * site.rast$vo +
        beta.dOM.i * site.rast$dOM +
        beta.dDM.i * site.rast$dDM +
        beta.shdi.i * site.rast$shdi
        
      # add
      add(all.log.rss) <- log.rss
      
    } # i
    
    # calculate mean and SD
    if (.log == TRUE) {
      
      out.rast <- c(app(all.log.rss, mean),
                    app(all.log.rss, sd))
      
    } else {
      
      out.rast <- c(app(exp(all.log.rss), mean),
                    app(exp(all.log.rss), sd))
      
    }
  
  } # season == "off"
  
  if (.season == "on") {
    
    # base model
    hsf <- M.on[[1]]
    
    # FR predictions + SEs
    # make rasters for prediction
    rast.stem <- subset(site.rast, c("a.stem", "TRT"))
    rast.ch <- subset(site.rast, c("a.ch", "TRT"))
    rast.cc2 <- subset(site.rast, c("a.cc", "TRT"))
    rast.dDM <- subset(site.rast, c("a.dDM", "TRT"))
    
    # names
    names(rast.stem)[1] <- "avail"
    names(rast.ch)[1] <- "avail"
    names(rast.cc2)[1] <- "avail"
    names(rast.dDM)[1] <- "avail"
    
    # FR predictions
    beta.stem <- predict(object = rast.stem, model = on.stem, fun = predict.gam, se.fit = T)
    beta.ch <- predict(object = rast.ch, model = on.ch, fun = predict.gam, se.fit = T)
    beta.cc2 <- predict(object = rast.cc2, model = on.cc2, fun = predict.gam, se.fit = T)
    beta.dDM <- predict(object = rast.dDM, model = on.dDM, fun = predict.gam, se.fit = T)
    
    # Monte Carlo predictions
    n.samp <- 100
    all.log.rss <- rast()
    
    for (i in 1:n.samp) {
      
      # sample from availability-weighted, spatially-varying betas
      beta.stem.i <- sample.betas(beta.stem)
      beta.ch.i <- sample.betas(beta.ch)
      beta.cc2.i <- sample.betas(beta.cc2)
      beta.dDM.i <- sample.betas(beta.dDM)
      
      # sample other coefs (non-spatially varying)
      beta.dOpen.i <- rnorm(1, hsf$mean[hsf$param == "dOpen"], hsf$sd[hsf$param == "dOpen"])
      beta.ed.i <- rnorm(1, hsf$mean[hsf$param == "ed"], hsf$sd[hsf$param == "ed"])
      beta.cc.i <- rnorm(1, hsf$mean[hsf$param == "cc"], hsf$sd[hsf$param == "cc"])
      beta.twi.i <- rnorm(1, hsf$mean[hsf$param == "twi"], hsf$sd[hsf$param == "twi"])
      beta.twi2.i <- rnorm(1, hsf$mean[hsf$param == "twi2"], hsf$sd[hsf$param == "twi2"])
      beta.vrm.i <- rnorm(1, hsf$mean[hsf$param == "vrm"], hsf$sd[hsf$param == "vrm"])
      beta.vrm2.i <- rnorm(1, hsf$mean[hsf$param == "vrm2"], hsf$sd[hsf$param == "vrm2"])
      
      # calculate log RSS prediction
      log.rss <- 
        
        # base
        beta.ed.i * site.rast$ed +
        beta.cc.i * site.rast$cc +
        beta.twi.i * site.rast$twi +
        beta.twi2.i * site.rast$twi2 +
        beta.vrm.i * site.rast$vrm +
        beta.vrm2.i * site.rast$vrm2 +
        
        # functional responses
        beta.stem.i * site.rast$stem +
        beta.ch.i * site.rast$ch +
        beta.cc2.i * site.rast$cc2 +
        beta.dOpen.i * site.rast$dOpen +
        beta.dDM.i * site.rast$dDM
      
      # add
      add(all.log.rss) <- log.rss
      
    } # i
    
    # calculate mean and SD
    if (.log == TRUE) {
      
      out.rast <- c(app(all.log.rss, mean),
                    app(all.log.rss, sd))
      
    } else {
      
      out.rast <- c(app(exp(all.log.rss), mean),
                    app(exp(all.log.rss), sd))
      
    }
    
  } # season == "on"
  
  return(out.rast)
  
}

# ______________________________________________________________________________
# 6. Map! ---- 
# ______________________________________________________________________________

# map just one
map_hsf <- function (.rast,
                     .site) {
  
  # subset site
  focal.site <- units |> filter(name == .site)
  
  # buffered site bounding box
  focal.bbox <- spatialEco::bbox_poly(focal.site) |>
    
    st_buffer(dist = 175) |>
    
    spatialEco::bbox_poly()
  
  # crop raster
  focal.rast <- crop(.rast, vect(focal.bbox)) 
    
  # clamp to reasonable value
  if (nlyr(focal.rast) == 2) {
    
    focal.rast$mean <- clamp(focal.rast$mean, upper = quantile(values(focal.rast$mean), prob = 0.99, na.rm = T))
    
  } else {
    
    focal.rast <- clamp(focal.rast, upper = quantile(values(focal.rast), prob = 0.99, na.rm = T))
    
  }
  
  names(focal.rast)[1] <- "mean"
    
  # plot
  ggplot() +
    
    theme_bw() +
    
    geom_raster(data = focal.rast,
                aes(fill = mean,
                    x = x,
                    y = y)) +
    
    geom_sf(data = focal.site,
            fill = NA,
            color = "white") +
    
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    
    scale_fill_viridis_c(name = "RSS")
  
}

# map two
map_hsf_2 <- function (.rast1,
                       .rast2,
                       .site) {
  
  # subset site
  focal.site <- units |> filter(name == .site)
  
  # buffered site bounding box
  focal.bbox <- spatialEco::bbox_poly(focal.site) |>
    
    st_buffer(dist = 175) |>
    
    spatialEco::bbox_poly()
  
  # crop raster
  focal.rast1 <- crop(.rast1, vect(focal.bbox))
  focal.rast2 <- crop(.rast2, vect(focal.bbox))
  
  # if multiple layers (i.e., output from pred_hsf_MC)
  
  if (nlyr(.rast1) == 2) {
    
   # clamp to reasonable value
   focal.rast1$mean <- clamp(focal.rast1$mean, upper = quantile(values(focal.rast1$mean), prob = 0.95, na.rm = T))
   focal.rast2$mean <- clamp(focal.rast2$mean, upper = quantile(values(focal.rast2$mean), prob = 0.95, na.rm = T))
  
   # add together
   focal.rasts <- c(focal.rast1$mean, focal.rast2$mean)
    
  } else {
    
    # clamp to reasonable value
    focal.rast1 <- clamp(focal.rast1, upper = quantile(values(focal.rast1), prob = 0.95, na.rm = T))
    focal.rast2 <- clamp(focal.rast2, upper = quantile(values(focal.rast2), prob = 0.95, na.rm = T))
    
    # add together
    focal.rasts <- c(focal.rast1, focal.rast2)
    
  }
  
  # change names
  names(focal.rasts) <- c("pre", "post")
  
  # plot
  ggplot() +
    
    theme_bw() +
    
    facet_wrap(~ lyr) +
    
    geom_spatraster(data = focal.rasts) +
    
    geom_sf(data = focal.site,
            fill = NA,
            color = "white") +
    
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(hjust = 0)) +
    
    scale_fill_viridis_c(name = "RSS")
  
}

# map change
map_hsf_change <- function (.rast1,
                            .rast2,
                            .site) {
  
  # subset site
  focal.site <- units |> filter(name == .site)
  
  # buffered site bounding box
  focal.bbox <- spatialEco::bbox_poly(focal.site) |>
    
    st_buffer(dist = 175) |>
    
    spatialEco::bbox_poly()
  
  # crop raster
  focal.rast1 <- crop(.rast1, vect(focal.bbox))
  focal.rast2 <- crop(.rast2, vect(focal.bbox))
  
  # difference - this will be the percent increase/decrease
  focal.rast.diff <- (1 - (focal.rast1 / focal.rast2)) * 10
  
  # which function output these?
  if (nlyr(focal.rast1) == 2) {
    
    # clamp to reasonable values
    focal.rast.diff$mean <- clamp(focal.rast.diff$mean, 
                                  lower = quantile(values(focal.rast.diff$mean), prob = 0.01, na.rm = T),
                                  upper = quantile(values(focal.rast.diff$mean), prob = 0.99, na.rm = T))
    
  } else {
    
    # clamp to reasonable values
    focal.rast.diff <- clamp(focal.rast.diff, 
                             lower = quantile(values(focal.rast.diff), prob = 0.01, na.rm = T),
                             upper = quantile(values(focal.rast.diff), prob = 0.99, na.rm = T))
    
  }
  
  names(focal.rast.diff)[1] <- "mean"
  
  # plot
  ggplot() +
    
    theme_bw() +
    
    geom_spatraster(data = focal.rast.diff$mean) +
    
    geom_sf(data = focal.site,
            fill = NA,
            color = "black") +
    
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "top") +
    
    scale_fill_gradient2(name = "% change", 
                         low = "#FF3300",
                         high = "purple",
                         midpoint = 0)
  
}

# ______________________________________________________________________________
# 7. Try these out ---- 
# ______________________________________________________________________________
# 7a. 1 - W Rabbit ----
# ______________________________________________________________________________

# 1A
off.pre.1A <- pred_hsf("1A", "off", "pre")
off.post.1A <- pred_hsf("1A", "off", "post")

on.pre.1A <- pred_hsf("1A", "on", "pre")
on.post.1A <- pred_hsf("1A", "on", "post")

map_hsf_2(off.pre.1A, off.post.1A, "1A")
map_hsf_2(on.pre.1A, on.post.1A, "1A")

map_hsf_change(off.pre.1A, off.post.1A, "1A")
map_hsf_change(on.pre.1A, on.post.1A, "1A")

# ______________________________________________________________________________

# 1B
off.pre.1B <- pred_hsf("1B", "off", "pre")
off.post.1B <- pred_hsf("1B", "off", "post")

on.pre.1B <- pred_hsf("1B", "on", "pre")
on.post.1B <- pred_hsf("1B", "on", "post")

# 1B
map_hsf_2(off.pre.1B, off.post.1B, "1B")
map_hsf_2(on.pre.1B, on.post.1B, "1B")

map_hsf_change(off.pre.1B, off.post.1B, "1B")
map_hsf_change(on.pre.1B, on.post.1B, "1B")

# ______________________________________________________________________________

# 1C
off.1C <- pred_hsf("1C", "off", "post")
on.1C <- pred_hsf("1C", "on", "post")

# 1C
map_hsf(off.1C, "1C")
map_hsf(on.1C, "1C")

# ______________________________________________________________________________
# 7b. 2 - Crazy Beetle Bug ----
# ______________________________________________________________________________

# 2A
off.pre.2A <- pred_hsf("2A", "off", "pre")
off.post.2A <- pred_hsf("2A", "off", "post")

on.pre.2A <- pred_hsf("2A", "on", "pre")
on.post.2A <- pred_hsf("2A", "on", "post")

map_hsf_2(off.pre.2A, off.post.2A, "2A")
map_hsf_2(on.pre.2A, on.post.2A, "2A")

map_hsf_change(off.pre.2A, off.post.2A, "2A")
map_hsf_change(on.pre.2A, on.post.2A, "2A")

# ______________________________________________________________________________

# 2B
off.pre.2B <- pred_hsf("2B", "off", "pre")
off.post.2B <- pred_hsf("2B", "off", "post")

on.pre.2B <- pred_hsf("2B", "on", "pre")
on.post.2B <- pred_hsf("2B", "on", "post")

map_hsf_2(off.pre.2B, off.post.2B, "2B")
map_hsf_2(on.pre.2B, on.post.2B, "2B")

map_hsf_change(off.pre.2B, off.post.2B, "2B")
map_hsf_change(on.pre.2B, on.post.2B, "2B")

# ______________________________________________________________________________

# 2C
off.2C <- pred_hsf("2C", "off", "post")
on.2C <- pred_hsf("2C", "on", "post")

map_hsf(off.2C, "2C")
map_hsf(on.2C, "2C")

# ______________________________________________________________________________
# 7c. 3 - Beetlejuice Bug ----
# ______________________________________________________________________________

# 3A
off.pre.3A <- pred_hsf("3A", "off", "pre")
off.post.3A <- pred_hsf("3A", "off", "post")

on.pre.3A <- pred_hsf("3A", "on", "pre")
on.post.3A <- pred_hsf("3A", "on", "post")

map_hsf_2(off.pre.3A, off.post.3A, "3A")
map_hsf_2(on.pre.3A, on.post.3A, "3A")

map_hsf_change(off.pre.3A, off.post.3A, "3A")
map_hsf_change(on.pre.3A, on.post.3A, "3A")

# ______________________________________________________________________________

# 3B
off.pre.3B <- pred_hsf("3B", "off", "pre")
off.post.3B <- pred_hsf("3B", "off", "post")

on.pre.3B <- pred_hsf("3B", "on", "pre")
on.post.3B <- pred_hsf("3B", "on", "post")

map_hsf_2(off.pre.3B, off.post.3B, "3B")
map_hsf_2(on.pre.3B, on.post.3B, "3B")

map_hsf_change(off.pre.3B, off.post.3B, "3B")
map_hsf_change(on.pre.3B, on.post.3B, "3B")

# ______________________________________________________________________________

# 3C
off.3C <- pred_hsf("3C", "off", "post")
on.3C <- pred_hsf("3C", "on", "post")

map_hsf(off.3C, "3C")
map_hsf(on.3C, "3C")

# ______________________________________________________________________________
# 7c. 4 - Chopaka ----
# ______________________________________________________________________________

# 4A
off.pre.4A <- pred_hsf("4A", "off", "pre")
off.post.4A <- pred_hsf("4A", "off", "post")

on.pre.4A <- pred_hsf("4A", "on", "pre")
on.post.4A <- pred_hsf("4A", "on", "post")

map_hsf_2(off.pre.4A, off.post.4A, "4A")
map_hsf_2(on.pre.4A, on.post.4A, "4A")

map_hsf_change(off.pre.4A, off.post.4A, "4A")
map_hsf_change(on.pre.4A, on.post.4A, "4A")

# ______________________________________________________________________________

# 4B
off.pre.4B <- pred_hsf("4B", "off", "pre")
off.post.4B <- pred_hsf("4B", "off", "post")

on.pre.4B <- pred_hsf("4B", "on", "pre")
on.post.4B <- pred_hsf("4B", "on", "post")

map_hsf_2(off.pre.4B, off.post.4B, "4B")
map_hsf_2(on.pre.4B, on.post.4B, "4B")

map_hsf_change(off.pre.4B, off.post.4B, "4B")
map_hsf_change(on.pre.4B, on.post.4B, "4B")

# ______________________________________________________________________________

# 4C
off.4C <- pred_hsf("4C", "off", "post")
on.4C <- pred_hsf("4C", "on", "post")

map_hsf(off.4C, "4C")
map_hsf(on.4C, "4C")

