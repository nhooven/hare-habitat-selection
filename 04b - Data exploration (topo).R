# PROJECT: Habitat selection
# SCRIPT: 04b - Data exploration (topo)
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
library(INLA)

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
plot_use_avail <- function (.covar,
                            .xlim = NULL) {
  
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
    
    ggtitle(paste0(.covar, " (snow-off)")) +
    
    coord_cartesian(xlim = .xlim)
  
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
    
    ggtitle(paste0(.covar, " (snow-on)")) +
    
    coord_cartesian(xlim = .xlim)
  
  return(cowplot::plot_grid(off.plot, on.plot))
  
}

# ______________________________________________________________________________
# 3a. Plots ----
# ______________________________________________________________________________

plot_use_avail("twi")    # intermediate peak, maybe selection during winter
plot_use_avail("vrm", .xlim = c(0, 0.015))   # slight positive
plot_use_avail("north")  # positive
plot_use_avail("east")   # positive

# ______________________________________________________________________________
# 4. Random slope parameterization ----

# while individuals could vary, I have good reason to think that site cluster
# is more important. We could drop the individual-level slopes if the models
# do better

# ______________________________________________________________________________
# 4a. Prepare data ----

hyper.list <- list(theta = list(initial = log(1),
                                fixed = F,
                                prior = "pc.prec",
                                param = c(3, 0.05)))

# ______________________________________________________________________________

# residuals for log(AKDE)
data.off$g.s <- residuals(lm(log(akde) ~ twi + twi2 + east + east2 + north + north2,
                             data = data.off))

data.on$g.s <- residuals(lm(log(akde) ~ twi + twi2 + east + east2 + north + north2,
                             data = data.on))

# correlation
covs.lin <- c("twi", "east", "north", "g.s")

cor(data.off |> dplyr::select(covs.lin), method = "spearman")
cor(data.on |> dplyr::select(covs.lin), method = "spearman")

# copy MRID and cluster
data.off <- data.off |>
  
  mutate(
    
    MRID1 = MRID,
    MRID2 = MRID,
    MRID3 = MRID,
    MRID4 = MRID,
    cluster = cluster,
    cluster1 = cluster,
    cluster2 = cluster,
    cluster3 = cluster
    
  )

data.on <- data.on |>
  
  mutate(
    
    MRID1 = MRID,
    MRID2 = MRID,
    MRID3 = MRID,
    MRID4 = MRID,
    cluster = cluster,
    cluster1 = cluster,
    cluster2 = cluster,
    cluster3 = cluster
    
  )

# ______________________________________________________________________________
# 4b. Model formulae ----
# ______________________________________________________________________________

M.indiv.form <- case ~ 
  
  # population-level effects
  g.s +
  twi + twi2 +
  north + 
  east + 
  
  # random intercepts
  f(MRID, model = "iid", hyper = list(theta = list(initial = log(1/1e6), fixed = T))) +
  
  # random slopes
  f(MRID1, twi, model = "iid", hyper = hyper.list) +
  f(MRID2, twi2, model = "iid", hyper = hyper.list) +
  f(MRID3, north, model = "iid", hyper = hyper.list) + 
  f(MRID4, east, model = "iid", hyper = hyper.list)

M.clust.form <- case ~ 
  
  # population-level effects
  g.s +
  twi + twi2 +
  north + 
  east + 
  
  # random intercepts
  f(MRID, model = "iid", hyper = list(theta = list(initial = log(1/1e6), fixed = T))) +
  
  # random slopes
  f(cluster, twi, model = "iid", hyper = hyper.list) +
  f(cluster1, twi2, model = "iid", hyper = hyper.list) +
  f(cluster2, north, model = "iid", hyper = hyper.list) + 
  f(cluster3, east, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 4c. Fit models ----

compute.list <- list(cpo = T,      # CPO/PIT (similar to PPPvals)
                     dic = T)

# ______________________________________________________________________________

# off
M.indiv.off.fit <- inla(M.indiv.form,
                        weights = data.off$w,
                        family = "binomial",
                        data = data.off,
                        control.compute = compute.list)  

M.clust.off.fit <- inla(M.clust.form,
                        weights = data.off$w,
                        family = "binomial",
                        data = data.off,
                        control.compute = compute.list) 

summary(M.indiv.off.fit)  # 148, 194919 ***
summary(M.clust.off.fit)  # 82,  1948987

# on
M.indiv.on.fit <- inla(M.indiv.form,
                        weights = data.on$w,
                        family = "binomial",
                        data = data.on,
                        control.compute = compute.list)  

M.clust.on.fit <- inla(M.clust.form,
                        weights = data.on$w,
                        family = "binomial",
                        data = data.on,
                        control.compute = compute.list) 

summary(M.indiv.on.fit)  
summary(M.clust.on.fit)  
