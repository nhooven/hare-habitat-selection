# PROJECT: Habitat selection
# SCRIPT: 08c - Visualization (Mapping)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 01 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 03 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)

# ______________________________________________________________________________
# 2. Read in cleaned model results and shapefiles ----
# ______________________________________________________________________________

off.M6 <- readRDS("model_results/off_M6.rds")
off.M4 <- readRDS("model_results/off_M4.rds")

on.M6 <- readRDS("model_results/on_M6.rds")
on.M5 <- readRDS("model_results/on_M5.rds")

# unit shapefile
units <- st_read("D:/hare_project/data_spatial/Units/units_fixed_utm/units_fixed_utm.shp") |>
  
  dplyr::select(name, geometry) |>
  
  arrange(name)

# ______________________________________________________________________________
# 3. Covariate summaries ----
# ______________________________________________________________________________

# covariate means, SDs, and ranges
off.summ <- readRDS("data_for_model/mean_sd_off.rds")
on.summ <- readRDS("data_for_model/mean_sd_on.rds")

# ______________________________________________________________________________
# 4. Process raster stack ----
# ______________________________________________________________________________

# read
rast.all <- rast("data_raster/rast_all.tif")

# functions
# standardize a raster layer function
stand_rast <- function (.layer, .var, .season = "off") {
  
  if (.season == "off") {
    
    layer.stand <- (.layer - off.summ$mean[off.summ$name == .var]) / off.summ$sd[off.summ$name == .var]
    
  } else {
    
    layer.stand <- (.layer - on.summ$mean[on.summ$name == .var]) / on.summ$sd[on.summ$name == .var]
    
  }
  
  return(layer.stand)
  
}

# prepare raster (availabilities and standardized)
prep_rast <- function (.site,
                       .radius = 150,          # ~ mean 99% HR radius
                       .season = "off",
                       .year = "pre") {
  
  # which vars?
  const.vars <- c("dOpen", "dDM", "ed", "twi", "vrm")
  
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
    
    st_buffer(dist = 500) |>
    
    spatialEco::bbox_poly()
  
  # subset and crop rast
  rast.crop <- rast.all |>
    
    subset(subset = which.vars) |>
    
    crop(vect(focal.bbox))
  
  # moving window analysis
  focal.mat <- focalMat(rast.crop, d = .radius, type = "circle") 
  focal.mat[focal.mat > 0] <- 1
  
  # apply over all layers
  rast.avail <- sapp(rast.crop, focal, w = focal.mat, "mean", na.rm = T)
  
  # names
  names(rast.avail) <- paste0("a.", var.names)
  names(rast.crop) <- var.names
  
  # standardize and bind together
  if (.season == "off") { 
    
    rast.stand.1 <- stand_rast(rast.crop$vo, "vo") 
    rast.stand.a1 <- stand_rast(rast.avail$a.vo, "a.vo") 
    
  } else {
      
    rast.stand.1 <- stand_rast(rast.crop$stem, "stem")
    rast.stand.a1 <- stand_rast(rast.avail$a.stem, "stem")
    
  }
  
  rast.stand <- c(
    
    # covariates
    rast.stand.1,
    stand_rast(rast.crop$ch, "ch", .season),
    stand_rast(rast.crop$ch^2, "ch2", .season),
    stand_rast(rast.crop$cc, "cc", .season),
    stand_rast(rast.crop$cc^2, "cc2", .season),
    stand_rast(log(rast.crop$dOpen + 1), "dOpen", .season),
    stand_rast(log(rast.crop$dDM + 1), "dDM", .season),
    stand_rast(rast.crop$ed, "ed", .season),
    stand_rast(rast.crop$twi, "twi", .season),
    stand_rast(rast.crop$twi^2, "twi2", .season),
    stand_rast(rast.crop$vrm, "vrm", .season),
    stand_rast(rast.crop$vrm^2, "vrm2", .season),
    
    # availabilities
    rast.stand.a1,
    stand_rast(rast.avail$a.ch, "a.ch", .season),
    stand_rast(rast.avail$a.cc, "a.cc", .season),
    stand_rast(log(rast.avail$a.dOpen + 1), "a.dOpen", .season),
    stand_rast(log(rast.avail$a.dDM + 1), "a.dDM", .season),
    stand_rast(rast.avail$a.ed, "a.ed", .season)
    
  )
  
  # change names for squared terms
  names(rast.stand)[c(3, 5, 10, 12)] <- c("ch2", "cc2", "twi2", "vrm2")
  
  # return
  return(rast.stand)
  
}

# use function
rast.test <- prep_rast("1A", .season = "off", .year = "pre")

# ______________________________________________________________________________
# 6. Calculate predictions ----

# instead of calculating over the entire study area, this is made for examining
# a site at a time

# returns a raster

# ______________________________________________________________________________

pred_hsf <- function (.results,
                     .rast,
                     .season = "off",
                     .trt = "unthinned") {
  
  # extract means (no credible intervals for now)
  est.mean <- .results |> dplyr::select(param, mean) |> rename(est = mean)
  
  # treatment indicators
  RET <- ifelse(.trt == "retention", 1, 0)
  PIL <- ifelse(.trt == "piling", 1, 0)
  
  # add "zero" parameters if necessary
  # if not, add a "zero" column
  all.params.off <- off.M6[[1]]$param
  #all.params.off <- on.M6[[1]]$param
  
  # which one?
  all.params <- if (.season == "off") { all.params.off } else { all.params.on }
  
  # calculate RSS
  calc_rss <- function (.est) {
    
    # add zeroes if necessary
    # helper function
    add_zeroes <- function (.est, all.params) {
      
      # determine which covariates are not present
      which.not <- all.params[which(all.params %notin% .est$param)]
      
      # bind in zero columns
      .est.new <- rbind(.est,
                        data.frame(param = which.not,
                                   est = 0))
    
      return(.est.new)
      
    } # f()
    
    .est <- add_zeroes(.est, all.params)
    
    # correct season
    if (.season == "off") {
      
      log.rss <- 
        
        # main effects
        .est$est[.est$param == "vo"] * .rast$vo +
        .est$est[.est$param == "ch"] * .rast$ch +
        .est$est[.est$param == "ch2"] * .rast$ch2 +
        .est$est[.est$param == "cc"] * .rast$cc +
        .est$est[.est$param == "cc2"] * .rast$cc2 +
        .est$est[.est$param == "dOpen"] * .rast$dOpen +
        .est$est[.est$param == "dDM"] * .rast$dDM +
        .est$est[.est$param == "ed"] * .rast$ed +
        .est$est[.est$param == "twi"] * .rast$twi +
        .est$est[.est$param == "twi2"] * .rast$twi2 +
        .est$est[.est$param == "vrm"] * .rast$vrm +
        .est$est[.est$param == "vrm2"] * .rast$vrm2 +
        
        # stand-level functional responses (base)
        .est$est[.est$param == "vo:a.vo"] * .rast$vo * .rast$a.vo +
        .est$est[.est$param == "ch:a.ch"] * .rast$ch * .rast$a.ch +
        .est$est[.est$param == "ch2:a.ch"] * .rast$ch2 * .rast$a.ch +
        .est$est[.est$param == "cc:a.cc"] * .rast$cc * .rast$a.cc +
        .est$est[.est$param == "cc2:a.cc"] * .rast$cc2 * .rast$a.cc +
        
        # landscape-level functional responses (base)
        .est$est[.est$param == "dOpen:a.dOpen"] * .rast$dOpen * .rast$a.dOpen +
        .est$est[.est$param == "dDM:a.dDM"] * .rast$dDM * .rast$a.dDM +
        .est$est[.est$param == "ed:a.ed"] * .rast$ed * .rast$a.ed +
        
        # stand level TRT
        .est$est[.est$param == "vo:year.trt:ret"] * .rast$vo * RET +
        .est$est[.est$param == "vo:year.trt:pil"] * .rast$vo * PIL +
        .est$est[.est$param == "ch:year.trt:ret"] * .rast$ch * RET +
        .est$est[.est$param == "ch:year.trt:pil"] * .rast$ch * PIL +
        .est$est[.est$param == "ch2:year.trt:ret"] * .rast$ch2 * RET +
        .est$est[.est$param == "ch2:year.trt:pil"] * .rast$ch2 * PIL +
        .est$est[.est$param == "cc:year.trt:ret"] * .rast$cc * RET +
        .est$est[.est$param == "cc:year.trt:pil"] * .rast$cc * PIL +
        .est$est[.est$param == "cc2:year.trt:ret"] * .rast$cc2 * RET +
        .est$est[.est$param == "cc2:year.trt:pil"] * .rast$cc2 * PIL +
        
        # landscape-level TRT
        .est$est[.est$param == "dOpen:year.trt:ret"] * .rast$dOpen * RET +
        .est$est[.est$param == "dOpen:year.trt:pil"] * .rast$dOpen * PIL +
        .est$est[.est$param == "dDM:year.trt:ret"] * .rast$dDM * RET +
        .est$est[.est$param == "dDM:year.trt:pil"] * .rast$dDM * PIL +
        .est$est[.est$param == "ed:year.trt:ret"] * .rast$ed * RET +
        .est$est[.est$param == "ed:year.trt:pil"] * .rast$ed * PIL +
        
        # stand-level FR x TRT
        .est$est[.est$param == "vo:a.vo:year.trt:ret"] * .rast$vo * .rast$a.vo * RET +
        .est$est[.est$param == "vo:a.vo:year.trt:pil"] * .rast$vo * .rast$a.vo * PIL +
        .est$est[.est$param == "ch:a.ch:year.trt:ret"] * .rast$ch * .rast$a.ch * RET +
        .est$est[.est$param == "ch:a.ch:year.trt:pil"] * .rast$ch * .rast$a.ch * PIL +
        .est$est[.est$param == "ch2:a.ch:year.trt:ret"] * .rast$ch2 * .rast$a.ch * RET +
        .est$est[.est$param == "ch2:a.ch:year.trt:pil"] * .rast$ch2 * .rast$a.ch * PIL +
        .est$est[.est$param == "cc:a.cc:year.trt:ret"] * .rast$cc * .rast$a.cc * RET +
        .est$est[.est$param == "cc:a.cc:year.trt:pil"] * .rast$cc * .rast$a.cc * PIL +
        .est$est[.est$param == "cc2:a.cc:year.trt:ret"] * .rast$cc2 * .rast$a.cc * RET +
        .est$est[.est$param == "cc2:a.cc:year.trt:pil"] * .rast$cc2 * .rast$a.cc * PIL +
        
        # landscape-level FR x TRT
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:ret"] * .rast$dOpen * .rast$a.dOpen * RET +
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:pil"] * .rast$dOpen * .rast$a.dOpen * PIL +
        .est$est[.est$param == "dDM:a.dDM:year.trt:ret"] * .rast$dDM * .rast$a.dDM * RET +
        .est$est[.est$param == "dDM:a.dDM:year.trt:pil"] * .rast$dDM * .rast$a.dDM * PIL +
        .est$est[.est$param == "ed:a.ed:year.trt:ret"] * .rast$ed * .rast$a.ed * RET +
        .est$est[.est$param == "ed:a.ed:year.trt:pil"] * .rast$ed * .rast$a.ed * PIL
      
    } else {
      
      log.rss <-
      
        # main effects
        .est$est[.est$param == "stem"] * .rast$stem +
        .est$est[.est$param == "ch"] * .rast$ch +
        .est$est[.est$param == "ch2"] * .rast$ch2 +
        .est$est[.est$param == "cc"] * .rast$cc +
        .est$est[.est$param == "cc2"] * .rast$cc2 +
        .est$est[.est$param == "dOpen"] * .rast$dOpen +
        .est$est[.est$param == "dDM"] * .rast$dDM +
        .est$est[.est$param == "ed"] * .rast$ed +
        .est$est[.est$param == "twi"] * .rast$twi +
        .est$est[.est$param == "twi2"] * .rast$twi2 +
        .est$est[.est$param == "vrm"] * .rast$vrm +
        .est$est[.est$param == "vrm2"] * .rast$vrm2 +
        
        # stand-level functional responses (base)
        .est$est[.est$param == "stem:a.stem"] * .rast$stem * .rast$a.stem +
        .est$est[.est$param == "ch:a.ch"] * .rast$ch * .rast$a.ch +
        .est$est[.est$param == "ch2:a.ch"] * .rast$ch2 * .rast$a.ch +
        .est$est[.est$param == "cc:a.cc"] * .rast$cc * .rast$a.cc +
        .est$est[.est$param == "cc2:a.cc"] * .rast$cc2 * .rast$a.cc +
        
        # landscape-level functional responses (base)
        .est$est[.est$param == "dOpen:a.dOpen"] * .rast$dOpen * .rast$a.dOpen +
        .est$est[.est$param == "dDM:a.dDM"] * .rast$dDM * .rast$a.dDM +
        .est$est[.est$param == "ed:a.ed"] * .rast$ed * .rast$a.ed +
        
        # stand level TRT
        .est$est[.est$param == "stem:year.trt:ret"] * .rast$stem * RET +
        .est$est[.est$param == "stem:year.trt:pil"] * .rast$stem * PIL +
        .est$est[.est$param == "ch:year.trt:ret"] * .rast$ch * RET +
        .est$est[.est$param == "ch:year.trt:pil"] * .rast$ch * PIL +
        .est$est[.est$param == "ch2:year.trt:ret"] * .rast$ch2 * RET +
        .est$est[.est$param == "ch2:year.trt:pil"] * .rast$ch2 * PIL +
        .est$est[.est$param == "cc:year.trt:ret"] * .rast$cc * RET +
        .est$est[.est$param == "cc:year.trt:pil"] * .rast$cc * PIL +
        .est$est[.est$param == "cc2:year.trt:ret"] * .rast$cc2 * RET +
        .est$est[.est$param == "cc2:year.trt:pil"] * .rast$cc2 * PIL +
        
        # landscape-level TRT
        .est$est[.est$param == "dOpen:year.trt:ret"] * .rast$dOpen * RET +
        .est$est[.est$param == "dOpen:year.trt:pil"] * .rast$dOpen * PIL +
        .est$est[.est$param == "dDM:year.trt:ret"] * .rast$dDM * RET +
        .est$est[.est$param == "dDM:year.trt:pil"] * .rast$dDM * PIL +
        .est$est[.est$param == "ed:year.trt:ret"] * .rast$ed * RET +
        .est$est[.est$param == "ed:year.trt:pil"] * .rast$ed * PIL +
        
        # stand-level FR x TRT
        .est$est[.est$param == "stem:a.stem:year.trt:ret"] * .rast$stem * .rast$a.stem * RET +
        .est$est[.est$param == "stem:a.stem:year.trt:pil"] * .rast$stem * .rast$a.stem * PIL +
        .est$est[.est$param == "ch:a.ch:year.trt:ret"] * .rast$ch * .rast$a.ch * RET +
        .est$est[.est$param == "ch:a.ch:year.trt:pil"] * .rast$ch * .rast$a.ch * PIL +
        .est$est[.est$param == "ch2:a.ch:year.trt:ret"] * .rast$ch2 * .rast$a.ch * RET +
        .est$est[.est$param == "ch2:a.ch:year.trt:pil"] * .rast$ch2 * .rast$a.ch * PIL +
        .est$est[.est$param == "cc:a.cc:year.trt:ret"] * .rast$cc * .rast$a.cc * RET +
        .est$est[.est$param == "cc:a.cc:year.trt:pil"] * .rast$cc * .rast$a.cc * PIL +
        .est$est[.est$param == "cc2:a.cc:year.trt:ret"] * .rast$cc2 * .rast$a.cc * RET +
        .est$est[.est$param == "cc2:a.cc:year.trt:pil"] * .rast$cc2 * .rast$a.cc * PIL +
        
        # landscape-level FR x TRT
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:ret"] * .rast$dOpen * .rast$a.dOpen * RET +
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:pil"] * .rast$dOpen * .rast$a.dOpen * PIL +
        .est$est[.est$param == "dDM:a.dDM:year.trt:ret"] * .rast$dDM * .rast$a.dDM * RET +
        .est$est[.est$param == "dDM:a.dDM:year.trt:pil"] * .rast$dDM * .rast$a.dDM * PIL +
        .est$est[.est$param == "ed:a.ed:year.trt:ret"] * .rast$ed * .rast$a.ed * RET +
        .est$est[.est$param == "ed:a.ed:year.trt:pil"] * .rast$ed * .rast$a.ed * PIL
      
    }
    
    return(log.rss)
    
  } # f()
  
  # calculate RSS
  rast.rss <- exp(calc_rss(est.mean))
  
  # layer name
  names(rast.rss) <- "RSS"
  
  return(rast.rss)
  
} # f()

# ______________________________________________________________________________
# 7. Map predictions ----
# ______________________________________________________________________________

map_hsf <- function (.rast,
                     .site,
                     .log = FALSE) {
  
  # subset site
  focal.site <- units |> filter(name == .site)
  
  # buffered site bounding box
  focal.bbox <- spatialEco::bbox_poly(focal.site) |>
    
    st_buffer(dist = 200) |>
    
    spatialEco::bbox_poly()
  
  # crop raster
  focal.rast <- crop(.rast, vect(focal.bbox))
  
  # log transform
  if (.log == T) { focal.rast <- log(focal.rast) }
  
  # plot
  ggplot() +
    
    theme_bw() +
    
    geom_raster(data = focal.rast,
                aes(fill = RSS,
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
    
    scale_fill_viridis_c()
  
}

# ______________________________________________________________________________
# 8. Try functions out ----
# ______________________________________________________________________________

# before vs after
test.site <- "1B"

rast.pre <- prep_rast(.site = test.site, .season = "off", .year = "pre")
rast.post <- prep_rast(.site = test.site, .season = "off", .year = "post")

map_hsf(pred_hsf(.results = off.M4[[1]],
                 .rast = rast.pre,
                 .season = "off",
                 .trt = "unthinned"),
        .site = test.site,
        .log = F)

map_hsf(pred_hsf(.results = off.M4[[1]],
                 .rast = rast.post,
                 .season = "off",
                 .trt = "retention"),
        .site = test.site,
        .log = T)

# controls
test.site <- "4C"

rast.ctrl <- prep_rast(.site = test.site, .season = "off", .year = "pre")

map_hsf(pred_hsf(.results = off.M6[[1]],
                 .rast = rast.ctrl,
                 .season = "off",
                 .trt = "unthinned"),
        .site = test.site)
