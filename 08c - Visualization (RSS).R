# PROJECT: Habitat selection
# SCRIPT: 08c - Visualization (RSS)
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
library(mgcv)

# ______________________________________________________________________________
# 2. Read in models and data ----
# ______________________________________________________________________________

# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.ch <- readRDS("model_results/fr_models/off_ch.rds")

# on
on.stem <- readRDS("model_results/fr_models/on_stem.rds")
on.ch <- readRDS("model_results/fr_models/on_ch.rds")
on.dEdge <- readRDS("model_results/fr_models/on_dEdge.rds")

# data
off.data <- readRDS("data_for_model/off_fr.rds")
on.data <- readRDS("data_for_model/on_fr.rds")

# means and SDs
mean.sd.off <- readRDS("data_for_model/mean_sd_off.rds")
mean.sd.on <- readRDS("data_for_model/mean_sd_on.rds")

mean.sd.off.trt <- readRDS("data_for_model/mean_sd_off_trt.rds")
mean.sd.on.trt <- readRDS("data_for_model/mean_sd_on_trt.rds")

# ______________________________________________________________________________
# 3. Predict 

# using treatment-specific means and ranges

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
                          mean.sd <- mean.sd.off
                          mean.sd.trt <- mean.sd.off.trt }
  
  if (.season == "on") { data.1 <- on.data 
                         col.values <- c("gray", "dodgerblue3", "dodgerblue3")
                         mean.sd <- mean.sd.on
                         mean.sd.trt <- mean.sd.on.trt }
  
  # choose correct parameter and availability
  data.1 <- data.1 |> filter(param == .param) |> rename("avail" = .avail)
  mean.sd.1 <- mean.sd |> filter(name == .param) 
  mean.sd.trt.1 <- mean.sd.trt |> filter(name == .param) 
  
  # x-axis title
  x.title <- case_when(
    
    .param == "vo" ~ "Visual obstruction (%)",
    .param == "stem" ~ "Conifer stem density (stems/ha)",
    .param == "ch" ~ "Canopy height (ft)",
    .param == "dEdge" ~ "Distance from edge (m)"
    
  )
  
  # variable ranges (UNTHIN, RET, PIL)
  var.range <- seq(mean.sd.1$min, 
                   mean.sd.1$max, 
                   length.out = 100)
  
  var.range.trt <- list(seq(mean.sd.trt.1$min[mean.sd.trt.1$TRT == "UNTHIN"], 
                            mean.sd.trt.1$max[mean.sd.trt.1$TRT == "UNTHIN"], 
                            length.out = 100),
                        seq(mean.sd.trt.1$min[mean.sd.trt.1$TRT == "RET"], 
                            mean.sd.trt.1$max[mean.sd.trt.1$TRT == "RET"], 
                            length.out = 100),
                        seq(mean.sd.trt.1$min[mean.sd.trt.1$TRT == "PIL"], 
                            mean.sd.trt.1$max[mean.sd.trt.1$TRT == "PIL"], 
                            length.out = 100))
  
  var.range.s <- (var.range - mean.sd.1$mean) / mean.sd.1$sd
    
  var.range.trt.s <- list((var.range.trt[[1]] - mean.sd.1$mean) / mean.sd.1$sd,
                          (var.range.trt[[2]] - mean.sd.1$mean) / mean.sd.1$sd,
                          (var.range.trt[[3]] - mean.sd.1$mean) / mean.sd.1$sd)
  
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
      
      mutate(fit = predict(.model, 
                           new.data.fr, 
                           se.fit = T, 
                           exclude = "s(cluster)",
                           newdata.guaranteed = TRUE)$fit,
             se.fit = predict(.model, 
                              new.data.fr, 
                              se.fit = T, 
                              exclude = "s(cluster)",
                              newdata.guaranteed = TRUE)$se.fit) |>
      
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
      
      mutate(fit = predict(.model, 
                           new.data.fr, 
                           se.fit = T, 
                           exclude = "s(cluster)",
                           newdata.guaranteed = TRUE)$fit,
             se.fit = predict(.model, 
                              new.data.fr, 
                              se.fit = T, 
                              exclude = "s(cluster)",
                              newdata.guaranteed = TRUE)$se.fit) |>
      
      # CIs
      mutate(low = fit - 1.645 * se.fit,
             upp = fit + 1.645 * se.fit)
    
    # RSS prediction
    new.data.rss <- bind_rows(
      
      data.frame(x = var.range.trt.s[[1]],
                 TRT = "UNTHIN"),
      data.frame(x = var.range.trt.s[[2]],
                 TRT = "RET"),
      data.frame(x = var.range.trt.s[[3]],
                 TRT = "PIL"),
      
    ) |>
      
      mutate(TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL"))) |>
      
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
      
      mutate(fit = predict(.model, 
                           new.data.fr, 
                           se.fit = T, 
                           exclude = "s(cluster)",
                           newdata.guaranteed = TRUE)$fit,
             se.fit = predict(.model, 
                              new.data.fr, 
                              se.fit = T, 
                              exclude = "s(cluster)",
                              newdata.guaranteed = TRUE)$se.fit) |>
      
      # CIs
      mutate(low = fit - 1.645 * se.fit,
             upp = fit + 1.645 * se.fit) |>
      
      # convert avail to factor for joining
      mutate(avail = round(avail, digits = 2) |> factor())
    
    # RSS prediction
    new.data.rss <- bind_rows(
      
      # UNTHIN
      data.frame(x = rep(var.range.trt.s[[1]], times = 3),
                 avail = rep(levels.avail, times = length(var.range.trt.s[[1]]) * 3) |>
                   round(digits = 2) |>
                   factor(),
                 TRT = "UNTHIN"),
      
      # RET
      data.frame(x = rep(var.range.trt.s[[2]], times = 3),
                 avail = rep(levels.avail, times = length(var.range.trt.s[[2]]) * 3) |>
                   round(digits = 2) |>
                   factor(),
                 TRT = "RET"),
      
      # PIL
      data.frame(x = rep(var.range.trt.s[[3]], times = 3),
                 avail = rep(levels.avail, times = length(var.range.trt.s[[3]]) * 3) |>
                   round(digits = 2) |>
                   factor(),
                 TRT = "PIL")
      
    ) |>
      
      mutate(TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL"))) |>
      
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
    
    if (.param == "vo") {
      
      out.plot <- out.plot + scale_x_continuous(breaks = c(0.4, 0.6, 0.8),
                                                labels = c(40, 60, 80))
      
    }
    
  }
  
  return(out.plot)
  
}

# use function
rss_pred_plot(off.vo, 4, "vo", "a.vo", "off") # ***
rss_pred_plot(off.ch, 4, "ch", "a.vo", "off")

rss_pred_plot(on.stem, 2, "stem", "a.stem", "on")
rss_pred_plot(on.ch, 3, "ch", "a.stem", "on")
rss_pred_plot(on.dEdge, 4, "dEdge", "a.stem", "on")
