# PROJECT: Habitat selection
# SCRIPT: 04a - Data exploration (stand-level)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 
# LAST MODIFIED: 29 May 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(hoove)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

data.off <- readRDS("data_for_model/off_data.rds")
data.on <- readRDS("data_for_model/on_data.rds")

# means/SDs
mean.sd.off <- readRDS("data_for_model/mean_sd_off.rds")
mean.sd.on <- readRDS("data_for_model/mean_sd_on.rds")

# ______________________________________________________________________________
# 3. Use/availability distributions ----

# generic function
plot_use_avail <- function (.covar) {
  
  # subset data
  focal.data.off <- data.off |> dplyr::select(c(case, .covar)) 
  focal.data.on <- data.on |> dplyr::select(c(case, .covar))
  
  # means/SDs
  focal.mean.sd.off <- mean.sd.off |> filter(name == .covar)
  focal.mean.sd.on <- mean.sd.on |> filter(name == .covar)
  
  # rename
  colnames(focal.data.off)[2] <- "var"
  colnames(focal.data.on)[2] <- "var"
  
  # back-transform
  focal.data.off$var <- (focal.data.off$var * focal.mean.sd.off$sd) + focal.mean.sd.off$mean
  focal.data.on$var <- (focal.data.on$var * focal.mean.sd.on$sd) + focal.mean.sd.on$mean
  
  # plots
  off.plot <- ggplot() +
    
    theme_hoove() +
    
    # mean
    geom_vline(xintercept = focal.mean.sd.off$mean,
               linetype = "dashed") +
    
    geom_density(data = focal.data.off,
                 aes(x = var,
                     fill = as.factor(case)),
                 alpha = 0.5) +
    
    scale_fill_manual(values = c("gray", "#FF3300")) +
    
    theme(legend.position = "none",
          axis.title = element_blank()) +
    
    ggtitle(paste0(.covar, " (snow-off)"))
  
  on.plot <- ggplot() +
    
    theme_hoove() +
    
    # mean
    geom_vline(xintercept = focal.mean.sd.on$mean,
               linetype = "dashed") +
    
    geom_density(data = focal.data.on,
                 aes(x = var,
                     fill = as.factor(case)),
                 alpha = 0.5) +
    
    scale_fill_manual(values = c("gray", "lightblue")) +
    
    theme(legend.position = "none",
          axis.title = element_blank()) +
    
    ggtitle(paste0(.covar, " (snow-on)"))
  
  return(cowplot::plot_grid(off.plot, on.plot))
  
}

# ______________________________________________________________________________
# 3a. Plots ----
# ______________________________________________________________________________

plot_use_avail("stem")   # positive
plot_use_avail("vo")     # positive
plot_use_avail("ch")     # intermediate peak
plot_use_avail("cc")     # intermediate peak (below mean off, above mean on)

# let's keep these parameterizations the same
