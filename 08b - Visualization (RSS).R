# PROJECT: Habitat selection
# SCRIPT: 08b - Visualization (RSS)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 08 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 08 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(mgcv)

# ______________________________________________________________________________
# 2. Read in models and data ----
# ______________________________________________________________________________

# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.dOpen <- readRDS("model_results/fr_models/off_dOpen.rds")
off.dDM <- readRDS("model_results/fr_models/off_dDM.rds")
off.ed <- readRDS("model_results/fr_models/off_ed.rds")

# on
on.stem <- readRDS("model_results/fr_models/on_stem.rds")
on.ch <- readRDS("model_results/fr_models/on_ch.rds")
on.cc2 <- readRDS("model_results/fr_models/on_cc2.rds")
on.dOpen <- readRDS("model_results/fr_models/on_dOpen.rds")
on.dDM <- readRDS("model_results/fr_models/on_dDM.rds")

# data
off.data <- readRDS("data_for_model/off_fr.rds")
on.data <- readRDS("data_for_model/on_fr.rds")

# means and SDs
mean.sd.off <- readRDS("data_for_model/mean_sd_off.rds")
mean.sd.on <- readRDS("data_for_model/mean_sd_on.rds")

# ______________________________________________________________________________
# 3. Predict and plot function ----
# ______________________________________________________________________________

rss_pred_plot <- function (
    
  .model,
  .which.model,
  .param,
  .avail,
  .season
  
) {
  
  # choose correct data, treatment color values, and means/SDs
  if (.season == "off") { data.1 <- off.data
                          col.values <- c("gray", "green4", "green4")
                          mean.sd <- mean.sd.off }
  
  if (.season == "on") { data.1 <- on.data 
                         col.values <- c("gray", "dodgerblue3", "dodgerblue3")
                         mean.sd <- mean.sd.on }
  
  # choose correct parameter and availability
  data.1 <- data.1 |> filter(param == .param) |> rename("avail" = .avail)
  mean.sd.1 <- mean.sd |> filter(name == .param) 
  
  # x-axis title
  x.title <- case_when(
    
    .param == "vo" ~ "Visual obstruction (%)",
    .param == "stem" ~ "Conifer stem density (stems/ha)",
    .param == "ch" ~ "Canopy height (ft)",
    .param == "dOpen" ~ "Distance from open (m)",
    .param == "dDM" ~ "Distance from dense mature (m)",
    .param == "ed" ~ "Edge density (m/ha)"
    
  )
  
  # variable range
  var.range <- seq(mean.sd.1$min, mean.sd.1$max, length.out = 100)
  var.range.s <- (var.range - mean.sd.1$mean) / mean.sd.1$sd
  
  # levels of availability (3 for now)
  levels.avail <- c(quantile(data.1$avail, prob = 0.05, na.rm = T),
                    quantile(data.1$avail, prob = 0.50, na.rm = T),
                    quantile(data.1$avail, prob = 0.95, na.rm = T))
  
  # M2 - FR 
  if (.which.model == 2) {
    
    # predict on functional response model
    new.data.fr <- data.frame(avail = levels.avail)
    
    # add predictions
    new.data.fr <- new.data.fr |>
      
      mutate(fit = predict(.model, new.data.fr, se.fit = T)$fit,
             se.fit = predict(.model, new.data.fr, se.fit = T)$se.fit) |>
      
      # CIs
      mutate(low = fit - 1.645 * se.fit,
             upp = fit + 1.645 * se.fit) |>
      
      # convert avail to factor for joining
      mutate(avail = round(avail, digits = 2) |> factor())
    
    # RSS prediction
    new.data.rss <- data.frame(x = rep(var.range.s, times = 3),
                               avail = rep(levels.avail, times = length(var.range.s)) |>
                                       round(digits = 2) |>
                                       factor()) |>
      
      left_join(new.data.fr |> dplyr::select(avail, fit, low, upp)) |>
      
      # calculate log RSS
      mutate(rss.est = x * fit,
             rss.low = x * low,
             rss.upp = x * upp) |>
      
      # back-transform x
      mutate(x.1 = (x * mean.sd.1$sd) + mean.sd.1$mean) |>
      
      # keep only relevant columns for plotting
      dplyr::select(x.1, avail, rss.est, rss.low, rss.upp) |>
      
      # label factors
      mutate(avail = factor(avail, labels = paste0("availability = ", levels(avail))))
    
    # plot
    out.plot <- ggplot(data = new.data.rss) +
      
      theme_bw() +
      
      facet_wrap(~ avail) +
      
      geom_hline(yintercept = 0,
                 linetype = "dashed") +
      
      geom_ribbon(aes(x = x.1,
                      y = rss.est,
                      ymin = rss.low,
                      ymax = rss.upp),
                  fill = col.values[2],
                  alpha = 0.25) +
      
      geom_line(aes(x = x.1,
                    y = rss.est),
                color = col.values[2],
                linewidth = 0.9) +
      
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            legend.title = element_blank(),
            legend.position = "none",
            strip.background = element_rect(color = NA),
            strip.text = element_text(hjust = 0)) +
      
      xlab(x.title) +
      ylab("log(Relative selection strength)")
    
  }
  
  # M3 - TRT
  if (.which.model == 3) {
    
    # predict on functional response model
    new.data.fr <- data.frame(TRT = c("UNTHIN", "RET", "PIL") |>
                              factor(levels = c("UNTHIN", "RET", "PIL")))
    
    # add predictions
    new.data.fr <- new.data.fr |>
      
      mutate(fit = predict(.model, new.data.fr, se.fit = T)$fit,
             se.fit = predict(.model, new.data.fr, se.fit = T)$se.fit) |>
      
      # CIs
      mutate(low = fit - 1.645 * se.fit,
             upp = fit + 1.645 * se.fit)
    
    # RSS prediction
    new.data.rss <- data.frame(x = rep(var.range.s, times = 3),
                               TRT  = c("UNTHIN", "RET", "PIL") |>
                                 rep(each = length(var.range.s)) |>
                                 factor(levels = c("UNTHIN", "RET", "PIL"))) |>
      
      left_join(new.data.fr |> dplyr::select(TRT, fit, low, upp)) |>
      
      # calculate log RSS
      mutate(rss.est = x * fit,
             rss.low = x * low,
             rss.upp = x * upp) |>
      
      # back-transform x
      mutate(x.1 = (x * mean.sd.1$sd) + mean.sd.1$mean) |>
      
      # keep only relevant columns for plotting
      dplyr::select(x.1, TRT, rss.est, rss.low, rss.upp) |>
      
      # label factors
      mutate(TRT = factor(TRT, labels = c("unthinned", "retention", "piling")))
    
    # plot
    out.plot <- ggplot(data = new.data.rss) +
      
      theme_bw() +
      
      facet_grid(~ TRT) +
      
      geom_hline(yintercept = 0,
                 linetype = "dashed") +
      
      geom_ribbon(aes(x = x.1,
                      y = rss.est,
                      ymin = rss.low,
                      ymax = rss.upp,
                      fill = TRT),
                  alpha = 0.25) +
      
      geom_line(aes(x = x.1,
                    y = rss.est,
                    color = TRT,
                    linetype = TRT),
                linewidth = 0.9) +
      
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            legend.title = element_blank(),
            legend.position = "none",
            strip.background = element_rect(color = NA),
            strip.text = element_text(hjust = 0)) +
      
      xlab(x.title) +
      ylab("log(Relative selection strength)") +
      
      scale_linetype_manual(values = c("solid", "solid", "dashed")) +
      scale_color_manual(values = col.values) +
      scale_fill_manual(values = col.values)
    
  }
  
  # M4 - FR x TRT
  if (.which.model == 4) {
    
    # predict on functional response model
    new.data.fr <- data.frame(avail = rep(levels.avail, times = 3),
                              TRT = c("UNTHIN", "RET", "PIL") |>
                                rep(each = 3) |>
                                factor(levels = c("UNTHIN", "RET", "PIL")))
    
    # add predictions
    new.data.fr <- new.data.fr |>
      
      mutate(fit = predict(.model, new.data.fr, se.fit = T)$fit,
             se.fit = predict(.model, new.data.fr, se.fit = T)$se.fit) |>
      
      # CIs
      mutate(low = fit - 1.645 * se.fit,
             upp = fit + 1.645 * se.fit) |>
      
      # convert avail to factor for joining
      mutate(avail = round(avail, digits = 2) |> factor())
    
    # RSS prediction
    new.data.rss <- data.frame(x = rep(var.range.s, times = 3 * 3),
                               avail = rep(levels.avail, times = length(var.range.s) * 3 * 3) |>
                                       round(digits = 2) |>
                                       factor(),
                               TRT  = c("UNTHIN", "RET", "PIL") |>
                                      rep(each = length(var.range.s) * 3 * 3) |>
                                      factor(levels = c("UNTHIN", "RET", "PIL"))) |>
      
      left_join(new.data.fr |> dplyr::select(avail, TRT, fit, low, upp)) |>
      
      # calculate log RSS
      mutate(rss.est = x * fit,
             rss.low = x * low,
             rss.upp = x * upp) |>
      
      # back-transform x
      mutate(x.1 = (x * mean.sd.1$sd) + mean.sd.1$mean) |>
      
      # keep only relevant columns for plotting
      dplyr::select(x.1, avail, TRT, rss.est, rss.low, rss.upp) |>
      
      # label factors
      mutate(TRT = factor(TRT, labels = c("unthinned", "retention", "piling")),
             avail = factor(avail, labels = paste0("availability = ", levels(avail))))
    
    # plot
    out.plot <- ggplot(data = new.data.rss) +
      
      theme_bw() +
      
      facet_grid(avail ~ TRT) +
      
      geom_hline(yintercept = 0,
                 linetype = "dashed") +
      
      geom_ribbon(aes(x = x.1,
                      y = rss.est,
                      ymin = rss.low,
                      ymax = rss.upp,
                      fill = TRT),
                  alpha = 0.25) +
      
      geom_line(aes(x = x.1,
                    y = rss.est,
                    color = TRT,
                    linetype = TRT),
                linewidth = 0.9) +
      
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            legend.title = element_blank(),
            legend.position = "none",
            strip.background = element_rect(color = NA),
            strip.text = element_text(hjust = 0)) +
      
      xlab(x.title) +
      ylab("log(Relative selection strength)") +
      
      scale_linetype_manual(values = c("solid", "solid", "dashed")) +
      scale_color_manual(values = col.values) +
      scale_fill_manual(values = col.values)
    
  }
  
  return(out.plot)
  
  }

# ______________________________________________________________________________
# 4. Use function ----
# ______________________________________________________________________________

rss_pred_plot(off.vo, 4, "vo", "a.vo", "off")
rss_pred_plot(off.dOpen, 2, "dOpen", "a.dOpen", "off")
rss_pred_plot(off.dDM, 3, "dDM", "a.dDM", "off")
rss_pred_plot(off.ed, 4, "ed", "a.ed", "off")

rss_pred_plot(on.stem, 4, "stem", "a.stem", "on")
rss_pred_plot(on.ch, 4, "ch", "a.ch", "on")
rss_pred_plot(on.dOpen, 4, "dOpen", "a.dOpen", "on")
rss_pred_plot(on.dDM, 4, "dDM", "a.dDM", "on")

# ______________________________________________________________________________
# 5. Canopy cover ----
# ______________________________________________________________________________

# add on model results
M.on <- readRDS("model_results/M_on.rds")

fr_pred_plot_cc <- function (.model = on.cc2) {
  
  # choose correct data, treatment color values, and means/SDs
  data.1 <- on.data 
  col.values <- c("gray", "dodgerblue3", "dodgerblue3")
  mean.sd <- mean.sd.on 
  
  # choose correct parameter and availability
  # cc
  data.1.cc <- data.1 |> filter(param == "cc") |> rename("avail" = "a.cc")
  mean.sd.1.cc <- mean.sd |> filter(name == "cc") 
  
  data.1.cc2 <- data.1 |> filter(param == "cc2") |> rename("avail" = "a.cc")
  mean.sd.1.cc2 <- mean.sd |> filter(name == "cc2") 
  
  # x-axis title
  x.title <- "Canopy cover (%)"
  
  # variable range
  var.range <- seq(mean.sd.1.cc$min, mean.sd.1.cc$max, length.out = 100)
  var.range2 <- var.range^2
  
  # standardize
  var.range.s <- (var.range - mean.sd.1.cc$mean) / mean.sd.1.cc$sd
  var.range2.s <- (var.range2 - mean.sd.1.cc2$mean) / mean.sd.1.cc2$sd
  
  # levels of availability (3 for now)
  levels.avail <- c(quantile(data.1.cc$avail, prob = 0.05, na.rm = T),
                    quantile(data.1.cc$avail, prob = 0.50, na.rm = T),
                    quantile(data.1.cc$avail, prob = 0.95, na.rm = T))
  
  # M2 - FR 
  # predict on functional response model
  # we'll construct CIs using Monte Carlo resampling
  new.data.fr <- data.frame(avail = levels.avail)
  
  # add predictions
  new.data.fr <- new.data.fr |>
    
    mutate(fit = predict(.model, new.data.fr, se.fit = T)$fit,
           se.fit = predict(.model, new.data.fr, se.fit = T)$se.fit) |>
    
    # convert avail to factor for joining
    mutate(avail = round(avail, digits = 2) |> factor())
  
  # RSS prediction
  new.data.rss <- data.frame(x = rep(var.range.s, times = 3),
                             x2 = rep(var.range2.s, times = 3),
                             avail = rep(levels.avail, times = length(var.range.s)) |>
                               round(digits = 2) |>
                               factor()) |>
    
    left_join(new.data.fr |> dplyr::select(avail, fit, se.fit)) 
  
  # resample for CIs
  # assume Gaussian***
  pred.matrix <- matrix(data = NA, nrow = nrow(new.data.rss), ncol = 500)
  
  for (i in 1:nrow(new.data.rss)) {
    
    for (j in 1:500 ) {
      
      cc.draw <- rnorm(1, M.on[[1]]$mean[M.on[[1]]$param == "cc"], M.on[[1]]$sd[M.on[[1]]$param == "cc"])
      cc2.draw <- rnorm(1, new.data.rss$fit[i], new.data.rss$se.fit[i])
      
      pred.matrix[i, j] <- new.data.rss$x[i] * cc.draw + new.data.rss$x2[i] * cc2.draw
      
    } # j
    
  } # i
  
  # add summaries to df
  new.data.rss <- new.data.rss |>
    
    mutate(rss.est = apply(pred.matrix, 1, mean),
           rss.low = apply(pred.matrix, 1, quantile, prob = 0.05),
           rss.upp = apply(pred.matrix, 1, quantile, prob = 0.95)) |>
    
    # back-transform x
    mutate(x.1 = (x * mean.sd.1.cc$sd) + mean.sd.1.cc$mean) |>
    
    # keep only relevant columns for plotting
    dplyr::select(x.1, avail, rss.est, rss.low, rss.upp) |>
    
    # label factors
    mutate(avail = factor(avail, labels = paste0("availability = ", levels(avail))))
  
  # plot
  out.plot <- ggplot(data = new.data.rss) +
    
    theme_bw() +
    
    facet_wrap(~ avail) +
    
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    
    geom_ribbon(aes(x = x.1,
                    y = rss.est,
                    ymin = rss.low,
                    ymax = rss.upp),
                fill = col.values[2],
                alpha = 0.25) +
    
    geom_line(aes(x = x.1,
                  y = rss.est),
              color = col.values[2],
              linewidth = 0.9) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.title = element_blank(),
          legend.position = "none",
          strip.background = element_rect(color = NA),
          strip.text = element_text(hjust = 0)) +
    
    xlab(x.title) +
    ylab("log(Relative selection strength)")
  
  return(out.plot)
  
}

fr_pred_plot_cc()
