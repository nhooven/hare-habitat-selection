# PROJECT: Habitat selection
# SCRIPT: 08a - Visualization (functional response predictions)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 05 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 05 Jun 2026
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

# ______________________________________________________________________________
# 3. Predict and plot function ----
# ______________________________________________________________________________

fr_pred_plot <- function (
    
  .model,
  .which.model,
  .param,
  .avail,
  .season
  
) {
  
  # choose correct data
  if (.season == "off") { data.1 <- off.data }
  if (.season == "on") { data.1 <- on.data }
  
  # treatment color values
  if (.season == "off") { col.values <- c("gray", "green4", "green4") }
  if (.season == "on") { col.values <- c("gray", "dodgerblue3", "dodgerblue3") }
  
  # choose correct parameter and availability
  data.1 <- data.1 |> filter(param == .param) |> rename("avail" = .avail)
  
  # x-axis title
  x.title <- case_when(
    
    .avail == "a.vo" ~ "Available visual obstruction (%)",
    .avail == "a.stem" ~ "Available conifer stem density (stems/ha)",
    .avail == "a.ch" ~ "Available canopy height (ft)",
    .avail == "a.dOpen" ~ "Available distance from open (m)",
    .avail == "a.dDM" ~ "Available distance from dense mature (m)",
    .avail == "a.ed" ~ "Available edge density (m/ha)"
    
  )
  
  # sequence of availability
  seq.avail <- seq(min(data.1$avail, na.rm = T), 
                   max(data.1$avail, na.rm = T), 
                   length.out = 100)
  
  # M2 - FR 
  if (.which.model == 2) {
    
    # newdata
    new.data <- data.frame(avail = seq.avail)
    
    # predict
    pred <- predict(.model, new.data, se.fit = T)
    
    # predict.df
    pred.df <- cbind(new.data,
                     data.frame(beta = pred$fit,
                                low = pred$fit - 1.645 * pred$se.fit,
                                upp = pred$fit + 1.645 * pred$se.fit))
    
    # plot
    out.plot <- ggplot(data = pred.df) +
      
      theme_classic() +
      
      geom_hline(yintercept = 0,
                 linetype = "dashed") +
      
      geom_ribbon(aes(x = avail,
                      y = beta,
                      ymin = low,
                      ymax = upp),
                  fill = col.values[2],
                  alpha = 0.25) +
      
      geom_line(aes(x = avail,
                    y = beta),
                color = col.values[2],
                linewidth = 0.9) +
      
      theme(legend.title = element_blank(),
            legend.position = "top") +
      
      xlab(x.title) +
      ylab("Selection coefficient")
    
  }
  
  # M3 - TRT
  if (.which.model == 3) {
    
    # newdata
    new.data <- data.frame(TRT = c("UNTHIN", "RET", "PIL") |>
                             factor(levels = c("UNTHIN", "RET", "PIL")))
    
    # predict
    pred <- predict(.model, new.data, se.fit = T)
    
    # predict.df
    pred.df <- cbind(new.data,
                     data.frame(beta = pred$fit,
                                low = pred$fit - 1.645 * pred$se.fit,
                                upp = pred$fit + 1.645 * pred$se.fit))
    
    # plot
    out.plot <- ggplot(data = pred.df) +
      
      theme_classic() +
      
      geom_hline(yintercept = 0,
                 linetype = "dashed") +
      
      geom_errorbar(aes(x = TRT,
                        y = beta,
                        ymin = low,
                        ymax = upp,
                        color = TRT),
                    width = 0,
                    linewidth = 1.2) +
      
      geom_point(aes(x = TRT,
                     y = beta,
                     fill = TRT),
                 shape = 21,
                 size = 2.5) +
      
      theme(legend.position = "none") +
      
      xlab("Treatment") +
      ylab("Selection coefficient") +
      
      scale_color_manual(values = col.values) +
      scale_fill_manual(values = col.values)
    
  }
  
  # M4 - FR x TRT
  if (.which.model == 4) {
    
    # newdata
    new.data <- data.frame(avail = rep(seq.avail, 3),
                           TRT = c("UNTHIN", "RET", "PIL") |>
                             rep(each = length(seq.avail)) |>
                             factor(levels = c("UNTHIN", "RET", "PIL")))
    
    # predict
    pred <- predict(.model, new.data, se.fit = T)
    
    # predict.df
    pred.df <- cbind(new.data,
                     data.frame(beta = pred$fit,
                                low = pred$fit - 1.645 * pred$se.fit,
                                upp = pred$fit + 1.645 * pred$se.fit))
    
    # plot
    out.plot <- ggplot(data = pred.df) +
      
      theme_classic() +
      
      geom_hline(yintercept = 0,
                 linetype = "dashed") +
      
      geom_ribbon(aes(x = avail,
                      y = beta,
                      ymin = low,
                      ymax = upp,
                      fill = TRT),
                  alpha = 0.25) +
      
      geom_line(aes(x = avail,
                    y = beta,
                    color = TRT,
                    linetype = TRT),
                linewidth = 0.9) +
      
      theme(legend.title = element_blank(),
            legend.position = "top") +
      
      xlab(x.title) +
      ylab("Selection coefficient") +
      
      scale_linetype_manual(values = c("solid", "solid", "dashed")) +
      scale_color_manual(values = col.values) +
      scale_fill_manual(values = col.values)
    
  }
  
  return(out.plot)
  
}

# ______________________________________________________________________________
# 4. Plots ----
# ______________________________________________________________________________

fr_pred_plot(off.vo, 4, "vo", "a.vo", "off")
fr_pred_plot(off.dOpen, 2, "dOpen", "a.dOpen", "off")
fr_pred_plot(off.dDM, 3, "dDM", "a.dDM", "off")
fr_pred_plot(off.ed, 4, "ed", "a.ed", "off")

fr_pred_plot(on.stem, 4, "stem", "a.stem", "on")
fr_pred_plot(on.ch, 4, "ch", "a.ch", "on")
fr_pred_plot(on.dOpen, 4, "dOpen", "a.dOpen", "on")
fr_pred_plot(on.dDM, 4, "dDM", "a.dDM", "on")

cowplot::plot_grid(fr_pred_plot(off.dOpen, 2, "dOpen", "a.dOpen", "off"),
                   fr_pred_plot(on.dOpen, 4, "dOpen", "a.dOpen", "on"))

cowplot::plot_grid(fr_pred_plot(off.dDM, 3, "dDM", "a.dDM", "off"),
                   fr_pred_plot(on.dDM, 4, "dDM", "a.dDM", "on"))
