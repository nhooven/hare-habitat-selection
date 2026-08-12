# PROJECT: Habitat selection
# SCRIPT: FIG - Base RSS and RIU
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 12 Aug 2026
# COMPLETED: 12 Aug 2026
# LAST MODIFIED: 12 Aug 2026
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

# ______________________________________________________________________________
# 3. RSS function ----
# ______________________________________________________________________________

plot_rss <- function (.var, .season) {
  
  # use correct model and means
  if (.season == "off") { .mod <- M.off } else { .mod <- M.on }
  if (.season == "off") { .mean.sd <- mean.sd.off } else { .mean.sd <- mean.sd.on }

  # linear coefficients
  if (.var %in% c("vo", "stem", "ch", "cc", "dEdge")) {
    
    # mean.sd (range will be in SDs)
    mean.sd.pop <- .mean.sd |> filter(name == .var) |>
      
      mutate(min.s = (min - mean) / sd,
             max.s = (max - mean) / sd)
    
    # availability sequence
    seq.avail <- seq(mean.sd.pop$min.s, mean.sd.pop$max.s, length.out = 100) 
  
    # pop beta
    beta.pop <- .mod[[1]] |> filter(param == .var)
    
    # indiv betas
    beta.ind <- .mod[[3]] |> dplyr::select(.var)
    
    # newdata
    # pop
    new.data.pop <- data.frame(x.s = seq.avail,
                               x = (seq.avail * mean.sd.pop$sd) + mean.sd.pop$mean,
                               y = seq.avail * beta.pop$mean)
    
    # ind
    new.data.ind <- data.frame()
    
    for (i in 1:nrow(beta.ind)) {
      
      new.data.focal <- data.frame(x.s = seq.avail,
                                   x = (seq.avail * mean.sd.pop$sd) + mean.sd.pop$mean,
                                   y = seq.avail * beta.ind[i, ],
                                   ind = i)
      
      new.data.ind <- rbind(new.data.ind, new.data.focal)
      
    } # i
      
      # otherwise, quadratic effects
    } else {
      
      # mean.sd (range will be in SDs)
      mean.sd.pop <- .mean.sd |> filter(name %in% c(.var, paste0(.var, 2))) |>
        
        mutate(min.s = (min - mean) / sd,
               max.s = (max - mean) / sd)
      
      # availability sequences
      seq.avail1 <- seq(mean.sd.pop$min.s[1], mean.sd.pop$max.s[1], length.out = 100) 
      seq.avail2 <- ((((seq.avail1 * mean.sd.pop$sd[1]) + mean.sd.pop$mean[1])^2) - 
                    mean.sd.pop$mean[2]) / mean.sd.pop$sd[2]
      
      # pop beta
      beta.pop <- .mod[[1]] |> filter(param %in% c(.var, paste0(.var, 2)))
      
      # indiv betas
      beta.ind <- .mod[[3]] |> dplyr::select(.var, paste0(.var, 2))
      
      # newdata
      # pop
      new.data.pop <- data.frame(x.s = seq.avail1,
                                 x = (seq.avail1 * mean.sd.pop$sd[1]) + mean.sd.pop$mean[1],
                                 y = seq.avail1 * beta.pop$mean[1] + seq.avail2 * beta.pop$mean[2])
      
      # ind
      new.data.ind <- data.frame()
      
      for (i in 1:nrow(beta.ind)) {
        
        new.data.focal <- data.frame(x.s = seq.avail1,
                                     x = (seq.avail1 * mean.sd.pop$sd[1]) + mean.sd.pop$mean[1],
                                     y = seq.avail1 * beta.ind[i, 1] + seq.avail2 * beta.ind[i, 2],
                                     ind = i)
        
        new.data.ind <- rbind(new.data.ind, new.data.focal)
      
    } # i
  
  }
  
  # prepare title for plot
  x.title <- case_when(
    
    .var == "vo" ~ "Visual obstruction (%)",
    .var == "stem" ~ "Stem density (stems/ha)",
    .var == "ch" ~ "Canopy height (ft)",
    .var == "cc" ~ "Canopy cover (%)",
    .var == "dEdge" ~ "Distance to mature edge (m)",
    .var == "twi" ~ "Wetness",
    .var == "vrm" ~ "Ruggedness"
    
  )
  
  # plot
  ggplot() +
    
    theme_classic() +
    
    # indifference line
    geom_hline(yintercept = 0) +
    
    # mean
    geom_vline(xintercept = mean.sd.pop$mean[1],
               linetype = "dashed") +
    
    # individual responses
    geom_line(data = new.data.ind,
              aes(x = x,
                  y = y,
                  group = ind),
              color = "gray",
              alpha = 0.35) +
    
    # population response
    geom_line(data = new.data.pop,
              aes(x = x,
                  y = y),
              color = ifelse(.season == "off",
                             "green4",
                             "dodgerblue3"),
              linewidth = 0.9) +
    
    # theme
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10)) +
    
    # limits
    coord_cartesian(ylim = c(-3, 3)) +
    
    # titles
    xlab(x.title)
  
}

# test
plot_rss("vrm", "on")

# ______________________________________________________________________________
# 3. RIU ----
# ______________________________________________________________________________
# 3a. Functions ----
# ______________________________________________________________________________

# population
calc_RIU_pop <- function (.season) {
  
  if (.season == "off") {
    
    # use/availability data
    hs.data <- data.off |>
      
      # keep only available locations
      filter(case == 0) |>
    
      dplyr::select(vo, ch, cc, twi, twi2, vrm, vrm2, dEdge)
    
    # base coefficients
    coef.base <- M.off[[1]] |> dplyr::select(param, mean)
    
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
    
    # .season == "off
  } else {
    
    # use/availability data
    hs.data <- data.on |>
      
      filter(case == 0) |>
      
      dplyr::select(stem, ch, cc, twi, twi2, vrm, vrm2, dEdge)
    
    # base coefficients
    coef.base <- M.on[[1]] |> dplyr::select(param, mean)
    
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
    
  } # .season == on
  
  return(hs.data)
  
} # f()

# use
off.RIU.pop <- calc_RIU_pop("off")
on.RIU.pop <- calc_RIU_pop("on")

# individual
calc_RIU_ind <- function (.season) {
  
  if (.season == "off") {
    
    # use/availability data
    hs.data <- data.off |>
      
      # keep only available locations
      filter(case == 0) |>
      
      dplyr::select(track_season_post, vo, ch, cc, twi, twi2, vrm, vrm2, dEdge)
    
    # calculate by individual
    data.ind.all <- data.frame()
    
    for (i in 1:length(unique(hs.data$track_season_post))) {
      
      # data
      data.ind <- hs.data |> filter(track_season_post == unique(hs.data$track_season_post)[i])
      
      # coefs
      coef.ind <- M.off[[3]][i, ]
      
      data.ind$pred.base <- exp(
        
        data.ind$twi * coef.ind$twi +
          data.ind$twi2 * coef.ind$twi2 +
          data.ind$vrm * coef.ind$vrm +
          data.ind$vrm2 * coef.ind$vrm2 +
          data.ind$vo * coef.ind$vo +
          data.ind$ch * coef.ind$ch +
          data.ind$cc * coef.ind$cc +
          data.ind$dEdge * coef.ind$dEdge
        
      )
      
      data.ind.all <- rbind(data.ind.all, data.ind)
      
    } # i
    
    # .season == "off
  } else {
    
    # use/availability data
    hs.data <- data.on |>
      
      # keep only available locations
      filter(case == 0) |>
      
      dplyr::select(track_season_post, stem, ch, cc, twi, twi2, vrm, vrm2, dEdge)
    
    # calculate by individual
    data.ind.all <- data.frame()
    
    for (i in 1:length(unique(hs.data$track_season_post))) {
      
      # data
      data.ind <- hs.data |> filter(track_season_post == unique(hs.data$track_season_post)[i])
      
      # coefs
      coef.ind <- M.on[[3]][i, ]
      
      data.ind$pred.base <- exp(
        
        data.ind$twi * coef.ind$twi +
          data.ind$twi2 * coef.ind$twi2 +
          data.ind$vrm * coef.ind$vrm +
          data.ind$vrm2 * coef.ind$vrm2 +
          data.ind$stem * coef.ind$stem +
          data.ind$ch * coef.ind$ch +
          data.ind$cc * coef.ind$cc +
          data.ind$dEdge * coef.ind$dEdge
        
      )
      
      data.ind.all <- rbind(data.ind.all, data.ind)
      
    } # i
    
  } # .season == on
  
  return(data.ind.all)
  
} # f()

# use
off.RIU.ind <- calc_RIU_ind("off")
on.RIU.ind <- calc_RIU_ind("on")

# ______________________________________________________________________________
# 3b. Process - unstandardize and add useful identifiers ----
# ______________________________________________________________________________

# function
process_RIU <- function (.RIU, .season) {
  
  if (.season == "off") {
    
    .mean.sd <- mean.sd.off
    
    # standardize
    .RIU.1 <- .RIU |>
      
      mutate(
        
        twi = (twi * .mean.sd$sd[.mean.sd$name == "twi"]) + .mean.sd$mean[.mean.sd$name == "twi"],
        vrm = (vrm * .mean.sd$sd[.mean.sd$name == "vrm"]) + .mean.sd$mean[.mean.sd$name == "vrm"],
        vo = (vo * .mean.sd$sd[.mean.sd$name == "vo"]) + .mean.sd$mean[.mean.sd$name == "vo"],
        ch = (ch * .mean.sd$sd[.mean.sd$name == "ch"]) + .mean.sd$mean[.mean.sd$name == "ch"],
        cc = (cc * .mean.sd$sd[.mean.sd$name == "cc"]) + .mean.sd$mean[.mean.sd$name == "cc"],
        dEdge = (dEdge * .mean.sd$sd[.mean.sd$name == "dEdge"]) + .mean.sd$mean[.mean.sd$name == "dEdge"]
        
      )
    
    # .season == "off
  } else {
    
    .mean.sd <- mean.sd.on
    
    # standardize
    .RIU.1 <- .RIU |>
      
      mutate(
        
        twi = (twi * .mean.sd$sd[.mean.sd$name == "twi"]) + .mean.sd$mean[.mean.sd$name == "twi"],
        vrm = (vrm * .mean.sd$sd[.mean.sd$name == "vrm"]) + .mean.sd$mean[.mean.sd$name == "vrm"],
        stem = (stem * .mean.sd$sd[.mean.sd$name == "stem"]) + .mean.sd$mean[.mean.sd$name == "stem"],
        ch = (ch * .mean.sd$sd[.mean.sd$name == "ch"]) + .mean.sd$mean[.mean.sd$name == "ch"],
        cc = (cc * .mean.sd$sd[.mean.sd$name == "cc"]) + .mean.sd$mean[.mean.sd$name == "cc"],
        dEdge = (dEdge * .mean.sd$sd[.mean.sd$name == "dEdge"]) + .mean.sd$mean[.mean.sd$name == "dEdge"]
        
      )
    
  } # .season == "on
  
  
  return(.RIU.1)
  
} # f()

# use
off.RIU.pop.1 <- process_RIU(off.RIU.pop, "off")
on.RIU.pop.1 <- process_RIU(on.RIU.pop, "on")

off.RIU.ind.1 <- process_RIU(off.RIU.ind, "off")
on.RIU.ind.1 <- process_RIU(on.RIU.ind, "on")

# ______________________________________________________________________________
# 3c. Plot ----
# ______________________________________________________________________________

plot_riu <- function (.var, .season) {
  
  if (.season == "off") { .mean.sd <- mean.sd.off } else { .mean.sd <- mean.sd.on }
  if (.season == "off") { .RIU.pop <- off.RIU.pop.1 } else { .RIU.pop <- on.RIU.pop.1 }
  if (.season == "off") { .RIU.ind <- off.RIU.ind.1 } else { .RIU.ind <- on.RIU.ind.1 }
  
  # correct variable
  .mean.sd <- .mean.sd |> filter(name == .var)
  .RIU.pop <- .RIU.pop |> dplyr::select(c(pred.base, .var))
  .RIU.ind <- .RIU.ind |> dplyr::select(c(pred.base, .var, track_season_post))
  
  # names
  colnames(.RIU.pop) <- c("pred.base", "var")
  colnames(.RIU.ind) <- c("pred.base", "var", "ind")
  
  # prepare title for plot
  x.title <- case_when(
    
    .var == "vo" ~ "Visual obstruction (%)",
    .var == "stem" ~ "Stem density (stems/ha)",
    .var == "ch" ~ "Canopy height (ft)",
    .var == "cc" ~ "Canopy cover (%)",
    .var == "dEdge" ~ "Distance to mature edge (m)",
    .var == "twi" ~ "Wetness",
    .var == "vrm" ~ "Ruggedness"
    
  )
  
  # plot
  ggplot() +
    
    theme_classic() +
    
    # mean
    geom_vline(xintercept = .mean.sd$mean,
               linetype = "dashed") +
    
    # individual responses
    geom_smooth(data = .RIU.ind,
              aes(x = var,
                  y = pred.base,
                  group = factor(ind)),
              color = "gray",
              alpha = 0.15,
              method = "gam",
              se = F,
              linewidth = 0.4) +
    
    # population response
    geom_smooth(data = .RIU.pop,
              aes(x = var,
                  y = pred.base),
              color = ifelse(.season == "off",
                             "green4",
                             "dodgerblue3"),
              method = "gam",
              se = F,
              linewidth = 0.9) +
    
    # theme
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10)) +
    
    # limits
    coord_cartesian(ylim = c(0, 7)) +
    
    # titles
    xlab(x.title)
    
} # f()

# ______________________________________________________________________________
# 4. Plot together ----
# ______________________________________________________________________________

plot_grid(
  
  plot_rss("vo", "off"), plot_rss("stem", "on"), plot_riu("vo", "off"), plot_riu("stem", "on"),
  plot_rss("ch", "off"), plot_rss("ch", "on"), plot_riu("ch", "off"), plot_riu("ch", "on"),
  plot_rss("cc", "off"), plot_rss("cc", "on"), plot_riu("cc", "off"), plot_riu("cc", "on"),
  plot_rss("dEdge", "off"), plot_rss("dEdge", "on"), plot_riu("dEdge", "off"), plot_riu("dEdge", "on"),
  
  plot_rss("twi", "off"), plot_rss("twi", "on"), plot_riu("twi", "off"), plot_riu("twi", "on"),
  plot_rss("vrm", "off"), plot_rss("vrm", "on"), plot_riu("vrm", "off"), plot_riu("vrm", "on"),
  
  nrow = 6
  
)

# 800 x 1350