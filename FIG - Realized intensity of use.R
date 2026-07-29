# PROJECT: Habitat selection
# SCRIPT: FIG - Realized intensity of use
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 28 Jul 2026
# COMPLETED: 28 Jul 2026
# LAST MODIFIED: 29 Jul 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(cowplot)

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

# FR models
# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.cc <- readRDS("model_results/fr_models/off_cc.rds")
off.dEdge <- readRDS("model_results/fr_models/off_dEdge.rds")

# on
on.cc <- readRDS("model_results/fr_models/on_cc.rds")
on.dEdge <- readRDS("model_results/fr_models/on_dEdge.rds")

# data
off.fr.data <- readRDS("data_for_model/off_fr.rds")
on.fr.data <- readRDS("data_for_model/on_fr.rds")

# HS data (bind in site column)
off.hs.data <- readRDS("data_for_model/off_data.rds") 
on.hs.data <- readRDS("data_for_model/on_data.rds")

# bind in site and sex columns
off.hs.data <- off.hs.data |>
  
  left_join(
    
    readRDS("data_cleaned/data_off.rds") |> 
      
      dplyr::select(site, sex, track_season_post) |>
      
      group_by(track_season_post) |>
      
      slice(1)
    
  )

on.hs.data <- on.hs.data |>
  
  left_join(
    
    readRDS("data_cleaned/data_on.rds") |> 
      
      dplyr::select(site, sex, track_season_post) |>
      
      group_by(track_season_post) |>
      
      slice(1)
    
  )

# ______________________________________________________________________________
# 3. Calculate realized intensity of use ----
# ______________________________________________________________________________
# 3a. Function ----
# ______________________________________________________________________________

calc_RIU <- function (.season) {
  
  if (.season == "off") {
    
    # use/availability data
    hs.data <- off.hs.data |>
      
      # add TRT and cluster
      mutate(TRT = c.trt,
             cluster = substr(site, 1, 1)) |>
      
      dplyr::select(track_season_post, case, TRT, cluster, 
                    vo, ch, cc, twi, twi2, vrm, vrm2, dEdge) |>
      
      rename(TSPID = track_season_post) |>
      
      # keep only available locations
      filter(case == 0)
    
    # base coefficients
    coef.base <- M.off[[1]] |> dplyr::select(param, mean)
    
    # FR coefficients
    coef.fr <- off.fr.data |> dplyr::select(TSPID, TRT, cluster, a.vo) |>
      
      group_by(TSPID) |> slice(1) |> ungroup() |>
      
      mutate(
        
        beta.vo = predict(off.vo, data.frame(TRT = TRT, cluster = cluster, avail = a.vo)),
        beta.cc = predict(off.cc, data.frame(TRT = TRT, cluster = cluster, avail = a.vo)),
        beta.dEdge = predict(off.dEdge, data.frame(TRT = TRT, cluster = cluster, avail = a.vo))
        
        
      ) |>
      
      dplyr::select(TSPID, beta.vo, beta.cc, beta.dEdge)
    
    # add implied coefs to all use/avail
    hs.data.fr <- hs.data |>
      
      dplyr::select(TSPID) |>
      
      left_join(coef.fr)
    
    # calculate w(x)
    hs.data$pred.base <- exp(
      
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$vo * coef.base$mean[coef.base$param == "vo"] +
        hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
        hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
        hs.data$dEdge * coef.base$mean[coef.base$param == "dEdge"]
      
    )
    
    hs.data$pred.fr <- exp(
      
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$vo * hs.data.fr$beta.vo +
        hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
        hs.data$cc * hs.data.fr$beta.cc +
        hs.data$dEdge * hs.data.fr$beta.dEdge
      
    )
    
    # .season == "off
  } else {
    
    # use/availability data
    hs.data <- on.hs.data |>
      
      # add TRT and cluster
      mutate(TRT = c.trt,
             cluster = substr(site, 1, 1)) |>
      
      dplyr::select(track_season_post, case, TRT, cluster, 
                    stem, ch, cc, twi, twi2, vrm, vrm2, dEdge) |>
      
      rename(TSPID = track_season_post) |>
      
      # keep only available locations
      filter(case == 0)
    
    # base coefficients
    coef.base <- M.on[[1]] |> dplyr::select(param, mean)
    
    # FR coefficients
    coef.fr <- on.fr.data |> dplyr::select(TSPID, TRT, cluster, a.stem) |>
      
      group_by(TSPID) |> slice(1) |> ungroup() |>
      
      mutate(
        
        beta.cc = predict(on.cc, data.frame(TRT = TRT, cluster = cluster, avail = a.stem)),
        beta.dEdge = predict(on.dEdge, data.frame(TRT = TRT, cluster = cluster, avail = a.stem))
        
        
      ) |>
      
      dplyr::select(TSPID, beta.cc, beta.dEdge)
    
    # add implied coefs to all use/avail
    hs.data.fr <- hs.data |>
      
      dplyr::select(TSPID) |>
      
      left_join(coef.fr)
    
    # calculate w(x)
    hs.data$pred.base <- exp(
      
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$stem * coef.base$mean[coef.base$param == "stem"] +
        hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
        hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
        hs.data$dEdge * coef.base$mean[coef.base$param == "dEdge"]
      
    )
    
    hs.data$pred.fr <- exp(
      
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$stem * coef.base$mean[coef.base$param == "stem"] +
        hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
        hs.data$cc * hs.data.fr$beta.cc +
        hs.data$dEdge * hs.data.fr$beta.dEdge
      
    )
    
  } # .season == on
  
  return(hs.data)
  
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
    .fr <- off.fr.data
    
    # standardize
    .RIU.1 <- .RIU |>
      
      mutate(
        
        twi = (twi * .mean.sd$sd[.mean.sd$name == "twi"]) + .mean.sd$mean[.mean.sd$name == "twi"],
        vrm = (vrm * .mean.sd$sd[.mean.sd$name == "vrm"]) + .mean.sd$mean[.mean.sd$name == "vrm"],
        vo = (vo * .mean.sd$sd[.mean.sd$name == "vo"]) + .mean.sd$mean[.mean.sd$name == "vo"],
        #stem = (stem * .mean.sd$sd[.mean.sd$name == "stem"]) + .mean.sd$mean[.mean.sd$name == "stem"],
        ch = (ch * .mean.sd$sd[.mean.sd$name == "ch"]) + .mean.sd$mean[.mean.sd$name == "ch"],
        cc = (cc * .mean.sd$sd[.mean.sd$name == "cc"]) + .mean.sd$mean[.mean.sd$name == "cc"],
        dEdge = (dEdge * .mean.sd$sd[.mean.sd$name == "dEdge"]) + .mean.sd$mean[.mean.sd$name == "dEdge"]
        
      ) |>
      
      # change canopy height to m
      mutate(ch = ch / 3.2) |>
      
      # keep columns we need
      dplyr::select(TSPID, twi, vrm, vo, ch, cc, dEdge, pred.base, pred.fr)
    
    # .season == "off
  } else {
    
    .mean.sd <- mean.sd.on
    .fr <- on.fr.data
    
    # standardize
    .RIU.1 <- .RIU |>
      
      mutate(
        
        twi = (twi * .mean.sd$sd[.mean.sd$name == "twi"]) + .mean.sd$mean[.mean.sd$name == "twi"],
        vrm = (vrm * .mean.sd$sd[.mean.sd$name == "vrm"]) + .mean.sd$mean[.mean.sd$name == "vrm"],
        #vo = (vo * .mean.sd$sd[.mean.sd$name == "vo"]) + .mean.sd$mean[.mean.sd$name == "vo"],
        stem = (stem * .mean.sd$sd[.mean.sd$name == "stem"]) + .mean.sd$mean[.mean.sd$name == "stem"],
        ch = (ch * .mean.sd$sd[.mean.sd$name == "ch"]) + .mean.sd$mean[.mean.sd$name == "ch"],
        cc = (cc * .mean.sd$sd[.mean.sd$name == "cc"]) + .mean.sd$mean[.mean.sd$name == "cc"],
        dEdge = (dEdge * .mean.sd$sd[.mean.sd$name == "dEdge"]) + .mean.sd$mean[.mean.sd$name == "dEdge"]
        
      ) |>
      
      # change canopy height to m
      mutate(ch = ch / 3.2) |>
      
      # keep columns we need
      dplyr::select(TSPID, twi, vrm, stem, ch, cc, dEdge, pred.base, pred.fr)
    
  } # .season == "on
  
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
# 5. Process for plotting ----
# ______________________________________________________________________________

# bind together
off.RIU.2 <- off.RIU.1 |> mutate(season = "snow-off") |>
  
  # factor levels
  mutate(TRT = factor(TRT,
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P")))

on.RIU.2 <- on.RIU.1 |> mutate(season = "snow-on") |>
  
  # factor levels
  mutate(TRT = factor(TRT,
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P")))

# ______________________________________________________________________________
# 6. Plots ----
# ______________________________________________________________________________
# 6a. Understory ----
# ______________________________________________________________________________

# off - vo
ggplot(data = off.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.off |> filter(name == "vo"),
             aes(xintercept = mean),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = vo,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.title.x = element_text(size = 9)) +
  
  # axis titles
  xlab("Visual obstruction (%)") +
  ylab("Realized intensity of use") +
  
  # axis scales
  scale_x_continuous(breaks = c(0.4, 0.6, 0.8),
                     labels = c(40, 60, 80)) +
  
  coord_cartesian(ylim = c(0.22, 4)) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> off.vo.plot

# on - stem
ggplot(data = on.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.on |> filter(name == "stem"),
             aes(xintercept = mean),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = stem,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.title.x = element_text(size = 9)) +
  
  # axis titles
  xlab("Stem density (stems/ha x 0.01)") +
  ylab("Realized intensity of use") +
  
  # axis scales
  scale_x_continuous(breaks = c(2500, 7500, 12500, 17500),
                     labels = c(25, 75, 125, 175)) +
  
  coord_cartesian(ylim = c(0.22, 4)) +
  
  scale_color_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_fill_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> on.stem.plot

# ______________________________________________________________________________
# 6b. Canopy height ----
# ______________________________________________________________________________

# off
ggplot(data = off.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.off |> filter(name == "ch"),
             aes(xintercept = mean / 3.2),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = ch,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.title = element_blank()) +
  
  # axis titles
  xlab("Canopy height (m)") +
  
  # axis limits
  coord_cartesian(ylim = c(0.22, 4),
                  xlim = c(1, 33)) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> off.ch.plot

# on - ch
ggplot(data = on.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.on |> filter(name == "ch"),
             aes(xintercept = mean / 3.2),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = ch,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9)) +
  
  # axis titles
  xlab("Canopy height (m)") +
  
  # axis limits
  coord_cartesian(ylim = c(0.22, 4),
                  xlim = c(1, 33)) +
  
  scale_color_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_fill_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> on.ch.plot

# ______________________________________________________________________________
# 6c. Canopy cover ----
# ______________________________________________________________________________

# off
ggplot(data = off.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.off |> filter(name == "cc"),
             aes(xintercept = mean),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = cc,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.title = element_blank()) +
  
  # axis titles
  xlab("Canopy cover (%)") +
  
  scale_x_continuous(breaks = seq(15, 75, 15)) +
  
  coord_cartesian(ylim = c(0.22, 4)) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> off.cc.plot

# on
ggplot(data = on.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.on |> filter(name == "cc"),
             aes(xintercept = mean),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = cc,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9)) +
  
  # axis titles
  xlab("Canopy cover (%)") +
  
  scale_x_continuous(breaks = seq(15, 75, 15)) +
  
  coord_cartesian(ylim = c(0.22, 4)) +
  
  scale_color_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_fill_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> on.cc.plot

# ______________________________________________________________________________
# 6d. Distance to edge ----
# ______________________________________________________________________________

# off
ggplot(data = off.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.off |> filter(name == "dEdge"),
             aes(xintercept = mean),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = dEdge,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = c(0.80, 0.75),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_blank()) +
  
  coord_cartesian(ylim = c(0.22, 4),
                  xlim = c(0, 275)) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> off.dEdge.plot

# on
ggplot(data = on.RIU.2) +
  
  theme_classic() +
  
  # mean lines
  geom_vline(data = mean.sd.on |> filter(name == "dEdge"),
             aes(xintercept = mean),
             linetype = "dashed") +
  
  # population level effects
  geom_smooth(aes(x = dEdge,
                  y = pred.fr,
                  color = TRT,
                  fill = TRT,
                  linetype = TRT),
              linewidth = 0.75,
              alpha = 0.25,
              method = "gam") +
  
  # theme arguments
  theme(legend.position = c(0.80, 0.75),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9)) +
  
  # axis titles
  xlab("Distance from mature edge (m)") +

  coord_cartesian(ylim = c(0.22, 4),
                  xlim = c(0, 275)) +
  
  scale_color_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_fill_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) -> on.dEdge.plot

# ______________________________________________________________________________
# 7. Plot together ----
# ______________________________________________________________________________

off.plots <- plot_grid(off.vo.plot, off.ch.plot, off.cc.plot, off.dEdge.plot,
                       nrow = 1)

on.plots <- plot_grid(on.stem.plot, on.ch.plot, on.cc.plot, on.dEdge.plot,
                       nrow = 1)

plot_grid(off.plots, on.plots, nrow = 2)

# 837 x 443

# NEXT: TWI and VRM plots?