# PROJECT: Habitat selection
# SCRIPT: 08e - Visualization (realized use)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 18 Jun 2026
# COMPLETED: 18 Jun 2026
# LAST MODIFIED: 23 Jul 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in models and data ----
# ______________________________________________________________________________

# results
M.off <- readRDS("model_results/M_off.rds")
M.on <- readRDS("model_results/M_on.rds")

# raw data
data.off <- readRDS("data_for_model/off_data.rds")
data.on <- readRDS("data_for_model/on_data.rds")

# means and SDs
mean.sd.off <- readRDS("data_for_model/mean_sd_off.rds")
mean.sd.on <- readRDS("data_for_model/mean_sd_on.rds")

mean.sd.off.trt <- readRDS("data_for_model/mean_sd_off_trt.rds")
mean.sd.on.trt <- readRDS("data_for_model/mean_sd_on_trt.rds")

# FR data
off.fr <- readRDS("data_for_model/off_fr.rds")
on.fr <- readRDS("data_for_model/on_fr.rds")

# ______________________________________________________________________________
# 3. Calculate realized intensity of use ----
# ______________________________________________________________________________
# 3a. Function ----
# ______________________________________________________________________________

calc_RIU <- function (.season) {
  
  if (.season == "off") {
    
    .model <- M.off
    .data <- data.off
    .mean.sd <- mean.sd.off
    
    # estimates
    pop.est <- .model[[1]]
    ind.est <- .model[[3]]
    
    # keep only available locations
    data.avail <- .data |> filter(case == 0)
    
    # split by TSPID
    data.split <- split(data.avail, ~track_season_post)
    
    # calculate population w(x)
    data.avail <- data.avail |>
      
      mutate(w.x.pop = exp(
        
          pop.est$mean[pop.est$param == "twi"] * twi +
          pop.est$mean[pop.est$param == "twi2"] * twi2 +
          pop.est$mean[pop.est$param == "vrm"] * vrm +
          pop.est$mean[pop.est$param == "vrm2"] * vrm2 +
          pop.est$mean[pop.est$param == "vo"] * vo +
          pop.est$mean[pop.est$param == "ch"] * ch +
          pop.est$mean[pop.est$param == "cc"] * cc +
          pop.est$mean[pop.est$param == "dEdge"] * dEdge
        
       )
      
      )
    
    # calculate individual-level w(x)
    # function
    calc_indiv_wx <- function (x) {
      
      ind.est.1 <- ind.est |> filter(TSPID == x$track_season_post[1])
      ind.avail <- data.avail |> filter(track_season_post == x$track_season_post[1])
      
      ind.avail <- ind.avail |> 
        
        mutate(w.x.ind = exp(
        
            ind.est.1$twi * twi +
            ind.est.1$twi2 * twi2 +
            ind.est.1$vrm * vrm +
            ind.est.1$vrm2 * vrm2 +
            ind.est.1$vo * vo +
            ind.est.1$ch * ch +
            ind.est.1$cc * cc +
            ind.est.1$dEdge * dEdge
        
         )
        
        ) |>
        
        dplyr::select(w.x.ind)
      
      return(ind.avail)
      
    }
    
    # apply function
    w.x.ind <- do.call(rbind, lapply(data.split, calc_indiv_wx)) 
    
    # final cleaning
    data.avail.1 <- data.avail |>
      
      # bind in
      bind_cols(w.x.ind) |>
      
      # keep only necessary columns
      dplyr::select(track_season_post,
                    cc:w.x.ind)
    
    return(data.avail.1)
    
    # .season == "off
  } else {
    
    .model <- M.on
    .data <- data.on
    .mean.sd <- mean.sd.on
    
    # estimates
    pop.est <- .model[[1]]
    ind.est <- .model[[3]]
    
    # keep only available locations
    data.avail <- .data |> filter(case == 0)
    
    # split by TSPID
    data.split <- split(data.avail, ~track_season_post)
    
    # calculate population w(x)
    data.avail <- data.avail |>
      
      mutate(w.x.pop = exp(
        
          pop.est$mean[pop.est$param == "twi"] * twi +
          pop.est$mean[pop.est$param == "twi2"] * twi2 +
          pop.est$mean[pop.est$param == "vrm"] * vrm +
          pop.est$mean[pop.est$param == "vrm2"] * vrm2 +
          pop.est$mean[pop.est$param == "stem"] * stem +
          pop.est$mean[pop.est$param == "ch"] * ch +
          pop.est$mean[pop.est$param == "cc"] * cc +
          pop.est$mean[pop.est$param == "dEdge"] * dEdge
        
      )
      
      )
    
    # calculate individual-level w(x)
    # function
    calc_indiv_wx <- function (x) {
      
      ind.est.1 <- ind.est |> filter(TSPID == x$track_season_post[1])
      ind.avail <- data.avail |> filter(track_season_post == x$track_season_post[1])
      
      ind.avail <- ind.avail |> 
        
        mutate(w.x.ind = exp(
          
            ind.est.1$twi * twi +
            ind.est.1$twi2 * twi2 +
            ind.est.1$vrm * vrm +
            ind.est.1$vrm2 * vrm2 +
            ind.est.1$stem * stem +
            ind.est.1$ch * ch +
            ind.est.1$cc * cc +
            ind.est.1$dEdge * dEdge
          
        )
        
        ) |>
        
        dplyr::select(w.x.ind)
      
      return(ind.avail)
      
    }
    
    # apply function
    w.x.ind <- do.call(rbind, lapply(data.split, calc_indiv_wx)) 
    
    # final cleaning
    data.avail.1 <- data.avail |>
      
      # bind in
      bind_cols(w.x.ind) |>
      
      # keep only necessary columns
      dplyr::select(track_season_post,
                    cc:w.x.ind)
    
    return(data.avail.1)
    
  } # .season == on
  
} # f()

# ______________________________________________________________________________
# 3b. Use function ----
# ______________________________________________________________________________

off.RIU <- calc_RIU("off")
on.RIU <- calc_RIU("on")

# ______________________________________________________________________________
# 4. Process - unstandardize and add useful identifiers ----
# ______________________________________________________________________________
# 4a. Function ----
# ______________________________________________________________________________

process_RIU <- function (.RIU, .season) {
  
  if (.season == "off") {
    
    .mean.sd <- mean.sd.off
    .fr <- off.fr
    
    # .season == "off
  } else {
    
    .mean.sd <- mean.sd.on
    .fr <- on.fr
    
  } # .season == "on
  
  # standardize
  .RIU.1 <- .RIU |>
    
    mutate(
      
      twi = (twi * .mean.sd$sd[.mean.sd$name == "twi"]) + .mean.sd$mean[.mean.sd$name == "twi"],
      vrm = (vrm * .mean.sd$sd[.mean.sd$name == "vrm"]) + .mean.sd$mean[.mean.sd$name == "vrm"],
      vo = (vo * .mean.sd$sd[.mean.sd$name == "vo"]) + .mean.sd$mean[.mean.sd$name == "vo"],
      stem = (stem * .mean.sd$sd[.mean.sd$name == "stem"]) + .mean.sd$mean[.mean.sd$name == "stem"],
      ch = (ch * .mean.sd$sd[.mean.sd$name == "ch"]) + .mean.sd$mean[.mean.sd$name == "ch"],
      cc = (cc * .mean.sd$sd[.mean.sd$name == "cc"]) + .mean.sd$mean[.mean.sd$name == "cc"],
      dEdge = (dEdge * .mean.sd$sd[.mean.sd$name == "dEdge"]) + .mean.sd$mean[.mean.sd$name == "dEdge"]
      
    ) |>
    
    # change canopy height to m
    mutate(ch = ch / 3.2) |>
    
    # change TSPID
    rename(TSPID = track_season_post) |>
    
    # keep columns we need
    dplyr::select(TSPID, twi, vrm, stem, vo, ch, cc, dEdge, w.x.pop, w.x.ind)
  
  # add in identifiers
  .fr.1 <- .fr |>
    
    group_by(TSPID) |>
    
    slice(1) |>
    
    # keep relevant columns
    dplyr::select(TSPID, sex, TRT, MRID, cluster)
  
  # join in
  .RIU.2 <- .RIU.1 |> left_join(.fr.1)
  
  return(.RIU.2)
  
} # f()

# ______________________________________________________________________________
# 4b. Use function ----
# ______________________________________________________________________________

off.RIU.1 <- process_RIU(off.RIU, "off")
on.RIU.1 <- process_RIU(on.RIU, "on")

# ______________________________________________________________________________
# 5. Plots ----

# bind together
all.RIU <- rbind(off.RIU.1 |> mutate(season = "snow-off"),
                 on.RIU.1 |> mutate(season = "snow-on")) |>
  
  # factor levels
  mutate(TRT = factor(TRT,
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("unthinned", "retention", "piling")))

# means
all.means.trt <- bind_rows(
  
  mean.sd.off.trt |> 
    
    dplyr::select(name, mean, TRT) |> 
    filter(name %in% c("ch", "cc", "dEdge", "stem", "twi", "vo", "vrm")) |>
    pivot_wider(names_from = name, values_from = mean) |>
    mutate(season = "snow-off",
           TRT = factor(TRT,
                        levels = c("UNTHIN", "RET", "PIL"),
                        labels = c("unthinned", "retention", "piling"))),
  
  mean.sd.on.trt |> 
    
    dplyr::select(name, mean, TRT) |> 
    filter(name %in% c("ch", "cc", "dEdge", "stem", "twi", "vo", "vrm")) |>
    pivot_wider(names_from = name, values_from = mean) |>
    mutate(season = "snow-on",
           TRT = factor(TRT,
                        levels = c("UNTHIN", "RET", "PIL"),
                        labels = c("unthinned", "retention", "piling")))
  
)

# ______________________________________________________________________________
# 5a. Wetness ----
# ______________________________________________________________________________

ggplot(data = all.RIU) +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  # mean lines
  geom_vline(data = all.means.trt,
             aes(xintercept = twi),
             linetype = "dashed") +
  
  # individual level effects
  geom_smooth(aes(x = twi,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = twi,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Topographic wetness (index)") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# ______________________________________________________________________________
# 5b. Ruggedness ----
# ______________________________________________________________________________

ggplot(data = all.RIU) +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  # mean lines
  geom_vline(data = all.means.trt,
             aes(xintercept = vrm),
             linetype = "dashed") +
  
  # individual level effects
  geom_smooth(aes(x = vrm,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = vrm,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Ruggedness (index)") +
  ylab("Realized intensity of use") +
  
  # axis scales
  scale_x_continuous(breaks = c(0, 0.03, 0.06, 0.09)) +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# ______________________________________________________________________________
# 5c. Canopy height ----
# ______________________________________________________________________________

ggplot(data = all.RIU) +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  # mean lines
  geom_vline(data = all.means.trt,
             aes(xintercept = ch / 3.2),
             linetype = "dashed") +
  
  # individual level effects
  geom_smooth(aes(x = ch,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = ch,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Canopy height (m)") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  # axis scales
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# ______________________________________________________________________________
# 5d. Canopy cover ----
# ______________________________________________________________________________

ggplot(data = all.RIU) +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  # mean lines
  geom_vline(data = all.means.trt,
             aes(xintercept = cc),
             linetype = "dashed") +
  
  # individual level effects
  geom_smooth(aes(x = cc,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = cc,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Canopy cover (%)") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  # axis scales
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# ______________________________________________________________________________
# 5d. Distance to edge ----
# ______________________________________________________________________________

ggplot(data = all.RIU) +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  # mean lines
  geom_vline(data = all.means.trt,
             aes(xintercept = dEdge),
             linetype = "dashed") +
  
  # individual level effects
  geom_smooth(aes(x = dEdge,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = dEdge,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Distance to edge (m)") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  # axis scales
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# ______________________________________________________________________________
# 5e. Visual obstruction ----
# ______________________________________________________________________________

ggplot(data = all.RIU |> filter(season == "snow-off")) +
  
  theme_bw() +
  
  facet_wrap(~ TRT) +
  
  # mean lines
  geom_vline(data = all.means.trt |> filter(season == "snow-off"),
             aes(xintercept = vo),
             linetype = "dashed") +
  
  # individual level effects
  geom_smooth(aes(x = vo,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = vo,
                  y = w.x.pop),
              color = "green4",
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Visual obstruction (%)") +
  ylab("Realized intensity of use") +
  
  # axis scales
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8),
                     labels = c(40, 60, 80)) +
  
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# 619 x 258

# ______________________________________________________________________________
# 5g. Stem density ----
# ______________________________________________________________________________

ggplot(data = all.RIU |> filter(season == "snow-on")) +
  
  theme_bw() +
  
  facet_wrap(~ TRT) +
  
  # mean lines
  geom_vline(data = all.means.trt |> filter(season == "snow-on"),
             aes(xintercept = stem),
             linetype = "dashed") +
  
  # individual level effects
  geom_smooth(aes(x = stem,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = stem,
                  y = w.x.pop),
              color = "dodgerblue3",
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Conifer stem density (stems/ha x 0.01)") +
  ylab("Realized intensity of use") +
  
  # axis scales
  scale_x_continuous(breaks = c(5000, 10000, 15000, 20000),
                     labels = c(50, 100, 150, 200)) +
  
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# 619 x 258

# ______________________________________________________________________________
# 6. Realized use of refugia, unit, and cover types ----
# ______________________________________________________________________________
# 6a. Rasters and all points ----

library(terra)
library(sf)

# ______________________________________________________________________________

# raster directory
dir.rast <- "D:/hare_project/data_spatial/Rasters/"

rast.cover.pre <- rast(paste0(dir.rast, "cover_type/cover_type_pre.tif"))
rast.cover.post <- rast(paste0(dir.rast, "cover_type/cover_type_post.tif"))
rast.dRet <- rast(paste0(dir.rast, "PCT/dRet.tif"))
rast.dPil <- rast(paste0(dir.rast, "PCT/dPiles.tif"))
rast.dUnitInt <- rast(paste0(dir.rast, "PCT/dUnitInt.tif"))
rast.dUnitExt <- rast(paste0(dir.rast, "PCT/dUnitExt.tif"))

pts <- readRDS("D:/hare_project/data_analysis/General/hare-gps-processing-new/data_cleaned/use_background.rds")

# since we dropped some NAs, we need to know which points were which
rast.all <- rast("data_raster/rast_all.tif")

# ______________________________________________________________________________
# 6b. Clean points ----
# ______________________________________________________________________________

pts.1 <- st_as_sf(pts, coords = c("x", "y"), crs = "epsg:32611") %>%
  
  # bind in extracted values
  bind_cols(
    
    .,
    
    # extract covariate values
    extract(x = rast.all, y = .) %>%
      
      # remove "ID"
      dplyr::select(-ID)
    
  ) |>
  
  filter(case == 0)

# correct values
pts.2 <- pts.1 |> 
    
    # attribute correct values
    mutate(
      
      # CONDITIONS
      cc = case_when(year == "PRE" ~ cc.pre,
                     year %in% c("POST1", "POST2") ~ cc.post),
      wsr = case_when(year == "PRE" ~ wsr.pre,
                      year %in% c("POST1", "POST2") ~ wsr.post),
      
      # LOCAL
      stem = case_when(year == "PRE" ~ stem.pre,
                       year %in% c("POST1", "POST2") ~ stem.post),
      vo = case_when(year == "PRE" ~ vo.pre,
                     year %in% c("POST1", "POST2") ~ vo.post),
      ch = case_when(year == "PRE" ~ ch.pre,
                     year %in% c("POST1", "POST2") ~ ch.post)
      
    ) |>
    
    # drop variables
    dplyr::select(
      
      -c(trt,
         season,
         cc.pre,
         cc.post,
         stem.pre,
         stem.post,
         vo.pre,
         vo.post,
         ch.pre,
         ch.post)
      
    ) |>
    
    # drop NAs
    drop_na(c(cc, twi, vrm, stem, vo, ch, dEdge)) |>
    
    # case weights
    mutate(w = ifelse(case == 0, 5000, 1)) |>
    
    # keep correct variables
    dplyr::select(geometry)

# sample the focal rasters and add to the RIU df
# combine rasters
rast.focal <- c(rast.cover.pre, rast.cover.post, rast.dRet, rast.dPil, rast.dUnitInt, rast.dUnitExt)

# extract
rast.focal.extract <- extract(rast.focal, pts.2)

colnames(rast.focal.extract) <- c("ID", "cover.pre", "cover.post", "dRet", "dPil", "dUnitInt", "dUnitExt")

# add in
all.RIU.1 <- all.RIU |>
  
  bind_cols(
    
    rast.focal.extract |>
      
      dplyr::select(-ID) |>
      
      mutate(in.out = ifelse(dUnitExt == 0, "in", "out"),
             dUnitEdge = ifelse(dUnitExt == 0,
                                -dUnitInt,
                                dUnitExt))
      
    
  ) |>
  
  # correct cover type
  mutate(cover = ifelse(TRT == "unthinned", cover.pre, cover.post))

# ______________________________________________________________________________
# 6c. Plots ----
# ______________________________________________________________________________

# in/out
ggplot(data = all.RIU.1) +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  # individual
  stat_summary(aes(x = in.out,
                   y = w.x.ind,
                   group = TSPID,
                   color = season),
               fun = mean,
               geom = "point",
               color = "gray",
               alpha = 0.15) +
  
  # population
  stat_summary(aes(x = in.out,
                   y = w.x.pop,
                   color = season),
               fun = mean,
               geom = "point") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Inside/outside unit") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  # axis scales
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# unit edge
ggplot(data = all.RIU.1) +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  geom_vline(xintercept = 0) +
  
  # individual level effects
  geom_smooth(aes(x = dUnitEdge,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = dUnitEdge,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Distance to unit edge (m)") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# cover
all.RIU.1 |>
  
  mutate(cover = factor(cover,
                        levels = 1:9,
                        labels = c("other",
                                   "unthinned doghair",
                                   "other",
                                   "open mature",
                                   "jackstraw",
                                   "dense mature",
                                   "dense mature",
                                   "other",
                                   "thinned doghair"))) |>
  
  drop_na(cover) |>

ggplot() +
  
  theme_bw() +
  
  facet_grid(season ~ TRT) +
  
  # individual
  stat_summary(aes(x = cover,
                   y = w.x.ind,
                   group = TSPID,
                   color = season),
               fun = mean,
               geom = "point",
               color = "gray",
               alpha = 0.15) +
  
  # population
  stat_summary(aes(x = cover,
                   y = w.x.pop,
                   color = season),
               fun = mean,
               geom = "point") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0),
        axis.title.x = element_blank()) +
  
  # axis titles
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  # axis scales
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

# distance to refugia
all.RIU.ret <- all.RIU.1 |> filter(TRT == "retention")
all.RIU.pil <- all.RIU.1 |> filter(TRT == "piling")

plot.ret <- ggplot(data = all.RIU.ret) +
  
  theme_bw() +
  
  facet_grid(~ season) +
  
  # individual level effects
  geom_smooth(aes(x = dRet,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = dRet,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Distance to retention patches (m)") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

plot.pil <- ggplot(data = all.RIU.pil) +
  
  theme_bw() +
  
  facet_grid(~ season) +
  
  # individual level effects
  geom_smooth(aes(x = dPil,
                  y = w.x.ind,
                  group = TSPID),
              color = "gray",
              se = F,
              method = "gam",
              linewidth = 0.1,
              alpha = 0.05) +
  
  # population level effects
  geom_smooth(aes(x = dPil,
                  y = w.x.pop,
                  color = season,
                  fill = season),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_rect(color = NA),
        strip.text = element_text(hjust = 0)) +
  
  # axis titles
  xlab("Distance to piles (m)") +
  ylab("Realized intensity of use") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue3")) +
  scale_fill_manual(values = c("green4", "dodgerblue3")) +
  
  scale_y_continuous(limits = c(0, 15)) +   # ensure that no values are dropped
  
  coord_cartesian(ylim = c(0.22, 4.5))

cowplot::plot_grid(plot.ret, plot.pil, nrow = 2)
