# PROJECT: Habitat selection
# SCRIPT: FIG - Functional response predictions
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 28 Jul 2026
# COMPLETED: 29 Jul 2026
# LAST MODIFIED: 29 Jul 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(cowplot)
library(patchwork)

# ______________________________________________________________________________
# 2. Read in models and data ----
# ______________________________________________________________________________

# model results
M.off <- readRDS("model_results/M_off.rds")[[1]]
M.on <- readRDS("model_results/M_on.rds")[[1]]

# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.cc <- readRDS("model_results/fr_models/off_cc.rds")
off.dEdge <- readRDS("model_results/fr_models/off_dEdge.rds")

# on
on.cc <- readRDS("model_results/fr_models/on_cc.rds")
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
# Theme ----
# ______________________________________________________________________________

theme_fr_rss <- function () {
  
  theme_classic() +
    
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          strip.text = element_text(size = 7))
  
}

# ______________________________________________________________________________
# 3. OFF - Visual obstruction ----

# M3

# ______________________________________________________________________________
# 3a. FR predictions ----
# ______________________________________________________________________________

# filter data
off.vo.data <- off.data |> filter(param == "vo")

# predictions
# newdata
off.vo.fr.newdata <- data.frame(TRT = c("UNTHIN", "RET", "PIL"))

# predict
off.vo.fr.pred <- predict(off.vo, 
                          off.vo.fr.newdata, 
                          exclude = "s(cluster)",
                          se.fit = T,
                          newdata.guaranteed = T)

# dfs for plotting
off.vo.fr.forPlot <- off.vo.fr.newdata |>
  
  mutate(pred = off.vo.fr.pred$fit,
         l90 = off.vo.fr.pred$fit - 1.645 * off.vo.fr.pred$se.fit,
         u90 = off.vo.fr.pred$fit + 1.645 * off.vo.fr.pred$se.fit,
         
         TRT = factor(TRT, 
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P")))

# plot
ggplot() +
  
  theme_fr_rss() +
  
  # population-level effect
  geom_rect(data = M.off |> filter(param == "vo"),
            aes(xmin = 0.4,
                xmax = 3.6,
                ymax = upp,
                ymin = low),
            fill = "gray97") +
  
  geom_hline(data = M.off |> filter(param == "vo"),
             aes(yintercept = mean),
             color = "darkgray") +
  
  # treatment-specific predictions
  geom_errorbar(data = off.vo.fr.forPlot,
                aes(x = TRT,
                    y = pred,
                    ymin = l90,
                    ymax = u90,
                    color = TRT),
                width = 0,
                linewidth = 1.2) +
  
  geom_point(data = off.vo.fr.forPlot,
             aes(x = TRT,
                 y = pred,
                 fill = TRT),
             shape = 21,
             color = "black",
             size = 2) +
  
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  
  xlab("Treatment") +
  ylab(expression(beta)) -> off.vo.fr.plot

# ______________________________________________________________________________
# 3b. RSS predictions ----
# ______________________________________________________________________________

# filter data
off.vo.mean.sd <- mean.sd.off |> filter(name == "vo")
off.vo.mean.sd.trt <- mean.sd.off.trt |> filter(name == "vo")

# availability ranges
off.vo.range.trt <- list(seq(off.vo.mean.sd.trt$min[off.vo.mean.sd.trt$TRT == "UNTHIN"], 
                             off.vo.mean.sd.trt$max[off.vo.mean.sd.trt$TRT == "UNTHIN"], 
                             length.out = 100),
                         seq(off.vo.mean.sd.trt$min[off.vo.mean.sd.trt$TRT == "RET"], 
                             off.vo.mean.sd.trt$max[off.vo.mean.sd.trt$TRT == "RET"], 
                             length.out = 100),
                         seq(off.vo.mean.sd.trt$min[off.vo.mean.sd.trt$TRT == "PIL"], 
                             off.vo.mean.sd.trt$max[off.vo.mean.sd.trt$TRT == "PIL"], 
                             length.out = 100))

off.vo.range.trt.s <- list((off.vo.range.trt[[1]] - off.vo.mean.sd$mean) / off.vo.mean.sd$sd,
                           (off.vo.range.trt[[2]] - off.vo.mean.sd$mean) / off.vo.mean.sd$sd,
                           (off.vo.range.trt[[3]] - off.vo.mean.sd$mean) / off.vo.mean.sd$sd)

# RSS prediction
off.vo.rss.newdata <- bind_rows(
  
  data.frame(x = off.vo.range.trt.s[[1]],
             TRT = "UNTHIN"),
  data.frame(x = off.vo.range.trt.s[[2]],
             TRT = "RET"),
  data.frame(x = off.vo.range.trt.s[[3]],
             TRT = "PIL"),
  
) |>
  
  mutate(TRT = factor(TRT, 
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P"))) |>
  
  left_join(off.vo.fr.forPlot |> dplyr::select(TRT, pred, l90, u90)) |>
  
  # calculate log RSS
  mutate(rss.est = x * pred,
         rss.low = x * l90,
         rss.upp = x * u90) |>
  
  # back-transform x
  mutate(x.1 = (x * off.vo.mean.sd$sd) + off.vo.mean.sd$mean) |>
  
  # keep only relevant columns for plotting
  dplyr::select(x.1, TRT, rss.est, rss.low, rss.upp)

# plot
ggplot(data = off.vo.rss.newdata) +
  
  theme_fr_rss() +
  
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
  
  theme(axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        
        legend.position = c(0.8, 0.29),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 5),
        legend.background = element_rect(fill = NA)) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  
  xlab("Visual obstruction (%)") +
  ylab("log(RSS)") -> off.vo.rss.plot

# ______________________________________________________________________________
# 4. OFF - Canopy cover ----

# M3

# ______________________________________________________________________________
# 4a. FR predictions ----
# ______________________________________________________________________________

# filter data
off.cc.data <- off.data |> filter(param == "cc")

# predictions
# newdata
off.cc.fr.newdata <- data.frame(TRT = c("UNTHIN", "RET", "PIL"))

# predict
off.cc.fr.pred <- predict(off.cc, 
                          off.cc.fr.newdata, 
                          exclude = "s(cluster)",
                          se.fit = T,
                          newdata.guaranteed = T)

# dfs for plotting
off.cc.fr.forPlot <- off.cc.fr.newdata |>
  
  mutate(pred = off.cc.fr.pred$fit,
         l90 = off.cc.fr.pred$fit - 1.645 * off.cc.fr.pred$se.fit,
         u90 = off.cc.fr.pred$fit + 1.645 * off.cc.fr.pred$se.fit,
         
         TRT = factor(TRT, 
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P")))

# plot
ggplot() +
  
  theme_fr_rss() +
  
  # population-level effect
  geom_rect(data = M.off |> filter(param == "cc"),
            aes(xmin = 0.4,
                xmax = 3.6,
                ymax = upp,
                ymin = low),
            fill = "gray97") +
  
  geom_hline(data = M.off |> filter(param == "cc"),
             aes(yintercept = mean),
             color = "darkgray") +
  
  # treatment-specific predictions
  geom_errorbar(data = off.cc.fr.forPlot,
                aes(x = TRT,
                    y = pred,
                    ymin = l90,
                    ymax = u90,
                    color = TRT),
                width = 0,
                linewidth = 1.2) +
  
  geom_point(data = off.cc.fr.forPlot,
             aes(x = TRT,
                 y = pred,
                 fill = TRT),
             shape = 21,
             color = "black",
             size = 2) +
  
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  
  xlab("Treatment") +
  ylab(expression(beta)) -> off.cc.fr.plot

# ______________________________________________________________________________
# 4b. RSS predictions ----
# ______________________________________________________________________________

# filter data
off.cc.mean.sd <- mean.sd.off |> filter(name == "cc")
off.cc.mean.sd.trt <- mean.sd.off.trt |> filter(name == "cc")

# availability ranges
off.cc.range.trt <- list(seq(off.cc.mean.sd.trt$min[off.cc.mean.sd.trt$TRT == "UNTHIN"], 
                             off.cc.mean.sd.trt$max[off.cc.mean.sd.trt$TRT == "UNTHIN"], 
                             length.out = 100),
                         seq(off.cc.mean.sd.trt$min[off.cc.mean.sd.trt$TRT == "RET"], 
                             off.cc.mean.sd.trt$max[off.cc.mean.sd.trt$TRT == "RET"], 
                             length.out = 100),
                         seq(off.cc.mean.sd.trt$min[off.cc.mean.sd.trt$TRT == "PIL"], 
                             off.cc.mean.sd.trt$max[off.cc.mean.sd.trt$TRT == "PIL"], 
                             length.out = 100))

off.cc.range.trt.s <- list((off.cc.range.trt[[1]] - off.cc.mean.sd$mean) / off.cc.mean.sd$sd,
                           (off.cc.range.trt[[2]] - off.cc.mean.sd$mean) / off.cc.mean.sd$sd,
                           (off.cc.range.trt[[3]] - off.cc.mean.sd$mean) / off.cc.mean.sd$sd)

# RSS prediction
off.cc.rss.newdata <- bind_rows(
  
  data.frame(x = off.cc.range.trt.s[[1]],
             TRT = "UNTHIN"),
  data.frame(x = off.cc.range.trt.s[[2]],
             TRT = "RET"),
  data.frame(x = off.cc.range.trt.s[[3]],
             TRT = "PIL"),
  
) |>
  
  mutate(TRT = factor(TRT, 
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P"))) |>
  
  left_join(off.cc.fr.forPlot |> dplyr::select(TRT, pred, l90, u90)) |>
  
  # calculate log RSS
  mutate(rss.est = x * pred,
         rss.low = x * l90,
         rss.upp = x * u90) |>
  
  # back-transform x
  mutate(x.1 = (x * off.cc.mean.sd$sd) + off.cc.mean.sd$mean) |>
  
  # keep only relevant columns for plotting
  dplyr::select(x.1, TRT, rss.est, rss.low, rss.upp)

# plot
ggplot(data = off.cc.rss.newdata) +
  
  theme_fr_rss() +
  
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
        legend.position = c(0.8, 0.3),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 5),
        legend.background = element_rect(fill = NA)) +
  
  scale_color_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_fill_manual(values = c("gray65", "green4", "darkorange3")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  
  scale_x_continuous(breaks = seq(15, 75, 15)) +
  
  xlab("Canopy cover (%)") +
  ylab("log(RSS)") -> off.cc.rss.plot

# ______________________________________________________________________________
# 5. OFF - dEdge ----

# M2

# ______________________________________________________________________________
# 5a. FR predictions ----
# ______________________________________________________________________________

# filter data
off.dEdge.data <- off.data |> filter(param == "dEdge")

# predictions
# availability sequence
off.dEdge.avail <- seq(quantile(off.dEdge.data$a.vo, prob = 0.05, na.rm = T), 
                       quantile(off.dEdge.data$a.vo, prob = 0.95, na.rm = T), 
                       length.out = 100)

# newdata
off.dEdge.fr.newdata <- data.frame(avail = off.dEdge.avail)

# predict
off.dEdge.fr.pred <- predict(off.dEdge, 
                             off.dEdge.fr.newdata, 
                             exclude = "s(cluster)",
                             se.fit = T,
                             newdata.guaranteed = T)

# dfs for plotting
off.dEdge.fr.forPlot <- off.dEdge.fr.newdata |>
  
  mutate(pred = off.dEdge.fr.pred$fit,
         l90 = off.dEdge.fr.pred$fit - 1.645 * off.dEdge.fr.pred$se.fit,
         u90 = off.dEdge.fr.pred$fit + 1.645 * off.dEdge.fr.pred$se.fit)

# plot
ggplot() +
  
  theme_fr_rss() +
  
  # population-level effect
  geom_rect(data = M.off |> filter(param == "dEdge"),
            aes(xmin = 0.51,
                xmax = 0.75,
                ymax = upp,
                ymin = low),
            fill = "gray97") +
  
  geom_hline(data = M.off |> filter(param == "dEdge"),
             aes(yintercept = mean),
             color = "darkgray") +
  
  # treatment-specific predictions
  geom_ribbon(data = off.dEdge.fr.forPlot,
                aes(x = avail,
                    y = pred,
                    ymin = l90,
                    ymax = u90),
              fill = "green4",
              alpha = 0.25) +
  
  geom_line(data = off.dEdge.fr.forPlot,
            aes(x = avail,
                y = pred),
             color = "green4",
             linewidth = 1.2) +
  
  theme(axis.title.y = element_blank()) +
  
  coord_cartesian(xlim = c(0.52, 0.72)) +
  
  xlab("Mean % visual obstruction") +
  ylab(expression(beta)) -> off.dEdge.fr.plot

# ______________________________________________________________________________
# 5b. RSS predictions ----
# ______________________________________________________________________________

# filter data
off.dEdge.mean.sd <- mean.sd.off |> filter(name == "dEdge")

# understory availability levels (just pull from fr df to join)
off.dEdge.avail.levels <- c(off.dEdge.fr.forPlot$avail[1],
                            off.dEdge.fr.forPlot$avail[50],
                            off.dEdge.fr.forPlot$avail[100])

# availability ranges
off.dEdge.range <- seq(off.dEdge.mean.sd$min, 
                       off.dEdge.mean.sd$max, 
                       length.out = 100)

off.dEdge.range.s <- (off.dEdge.range - off.dEdge.mean.sd$mean) / off.dEdge.mean.sd$sd

# RSS prediction
off.dEdge.rss.newdata <- data.frame(x = rep(off.dEdge.range.s, 3)) |>

  # add availability levels
  mutate(avail = c(rep(off.dEdge.avail.levels[1], 100),
                   rep(off.dEdge.avail.levels[2], 100),
                   rep(off.dEdge.avail.levels[3], 100))) |>
         
  left_join(off.dEdge.fr.forPlot |> dplyr::select(avail, pred, l90, u90)) |>
  
  # calculate log RSS
  mutate(rss.est = x * pred,
         rss.low = x * l90,
         rss.upp = x * u90) |>
  
  # back-transform x
  mutate(x.1 = (x * off.dEdge.mean.sd$sd) + off.dEdge.mean.sd$mean) |>
  
  # availability to sensible factor
  mutate(avail = factor(avail,
                        levels = unique(factor(avail)),
                        labels = c("Mean VO = 0.51",
                                   "Mean VO = 0.63",
                                   "Mean VO = 0.75"))) |>
  
  # keep only relevant columns for plotting
  dplyr::select(x.1, avail, rss.est, rss.low, rss.upp)

# plot
ggplot(data = off.dEdge.rss.newdata) +
  
  theme_fr_rss() +
  
  facet_wrap(~ avail) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_ribbon(aes(x = x.1,
                  y = rss.est,
                  ymin = rss.low,
                  ymax = rss.upp),
              fill = "green4",
              alpha = 0.25) +
  
  geom_line(aes(x = x.1,
                y = rss.est),
            color = "green4",
            linewidth = 0.9) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  
  geom_text(aes(x = 2,
                y = -0.82,
                label = avail),
            size = 2,
            hjust = 0,
            check_overlap = T) +
  
  xlab("Distance to mature edge (m)") +
  ylab("log(RSS)") -> off.dEdge.rss.plot

# ______________________________________________________________________________
# 6. ON - Canopy cover ----

# M2

# ______________________________________________________________________________
# 6a. FR predictions ----
# ______________________________________________________________________________

# filter data
on.cc.data <- on.data |> filter(param == "cc")

# predictions
# availability sequence
on.cc.avail <- seq(quantile(on.cc.data$a.stem, prob = 0.05, na.rm = T), 
                   quantile(on.cc.data$a.stem, prob = 0.95, na.rm = T), 
                   length.out = 100)

# newdata
on.cc.fr.newdata <- data.frame(avail = on.cc.avail)

# predict
on.cc.fr.pred <- predict(on.cc, 
                         on.cc.fr.newdata, 
                         exclude = "s(cluster)",
                         se.fit = T,
                         newdata.guaranteed = T)

# dfs for plotting
on.cc.fr.forPlot <- on.cc.fr.newdata |>
  
  mutate(pred = on.cc.fr.pred$fit,
         l90 = on.cc.fr.pred$fit - 1.645 * on.cc.fr.pred$se.fit,
         u90 = on.cc.fr.pred$fit + 1.645 * on.cc.fr.pred$se.fit)

# plot
ggplot() +
  
  theme_fr_rss() +
  
  # population-level effect
  geom_rect(data = M.on |> filter(param == "cc"),
            aes(xmin = 3500,
                xmax = 10500,
                ymax = upp,
                ymin = low),
            fill = "gray97") +
  
  geom_hline(data = M.on |> filter(param == "cc"),
             aes(yintercept = mean),
             color = "darkgray") +
  
  # treatment-specific predictions
  geom_ribbon(data = on.cc.fr.forPlot,
              aes(x = avail,
                  y = pred,
                  ymin = l90,
                  ymax = u90),
              fill = "dodgerblue3",
              alpha = 0.25) +
  
  geom_line(data = on.cc.fr.forPlot,
            aes(x = avail,
                y = pred),
            color = "dodgerblue3",
            linewidth = 1.2) +
  
  coord_cartesian(xlim = c(3900, 10000)) +
  
  theme(axis.title.y = element_blank()) +
  
  xlab("Mean stems/ha") +
  ylab(expression(beta)) -> on.cc.fr.plot

# ______________________________________________________________________________
# 6b. RSS predictions ----
# ______________________________________________________________________________

# filter data
on.cc.mean.sd <- mean.sd.on |> filter(name == "cc")

# understory availability levels (just pull from fr df to join)
on.cc.avail.levels <- c(on.cc.fr.forPlot$avail[1],
                            on.cc.fr.forPlot$avail[50],
                            on.cc.fr.forPlot$avail[100])

# availability ranges
on.cc.range <- seq(on.cc.mean.sd$min, 
                       on.cc.mean.sd$max, 
                       length.out = 100)

on.cc.range.s <- (on.cc.range - on.cc.mean.sd$mean) / on.cc.mean.sd$sd

# RSS prediction
on.cc.rss.newdata <- data.frame(x = rep(on.cc.range.s, 3)) |>
  
  # add availability levels
  mutate(avail = c(rep(on.cc.avail.levels[1], 100),
                   rep(on.cc.avail.levels[2], 100),
                   rep(on.cc.avail.levels[3], 100))) |>
  
  left_join(on.cc.fr.forPlot |> dplyr::select(avail, pred, l90, u90)) |>
  
  # calculate log RSS
  mutate(rss.est = x * pred,
         rss.low = x * l90,
         rss.upp = x * u90) |>
  
  # back-transform x
  mutate(x.1 = (x * on.cc.mean.sd$sd) + on.cc.mean.sd$mean) |>
  
  # availability to sensible factor
  mutate(avail = factor(avail,
                        levels = unique(factor(avail)),
                        labels = c("Mean stems = 3580",
                                   "Mean stems = 6982",
                                   "Mean stems = 10454"))) |>
  
  # keep only relevant columns for plotting
  dplyr::select(x.1, avail, rss.est, rss.low, rss.upp)

# plot
ggplot(data = on.cc.rss.newdata) +
  
  theme_fr_rss() +
  
  facet_wrap(~ avail) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  
  geom_ribbon(aes(x = x.1,
                  y = rss.est,
                  ymin = rss.low,
                  ymax = rss.upp),
              fill = "dodgerblue3",
              alpha = 0.25) +
  
  geom_line(aes(x = x.1,
                y = rss.est),
            color = "dodgerblue3",
            linewidth = 0.9) +
  
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  
  geom_text(aes(x = 1,
                y = -1.27,
                label = avail),
            size = 2,
            hjust = -0.3,
            check_overlap = T) +
  
  scale_x_continuous(breaks = seq(15, 75, 15)) +
  
  xlab("Canopy cover (%)") +
  ylab("log(RSS)") -> on.cc.rss.plot

# ______________________________________________________________________________
# 7. ON - dEdge ----

# M4

# ______________________________________________________________________________
# 7a. FR predictions ----
# ______________________________________________________________________________

# filter data
on.dEdge.data <- on.data |> filter(param == "dEdge")

# predictions
# availability sequence
on.dEdge.avail <- seq(quantile(on.dEdge.data$a.stem, prob = 0.05, na.rm = T), 
                      quantile(on.dEdge.data$a.stem, prob = 0.95, na.rm = T), 
                      length.out = 100)

# newdata
on.dEdge.fr.newdata <- data.frame(avail = rep(on.dEdge.avail, 3),
                                  TRT = rep(c("UNTHIN", "RET", "PIL"),
                                            each = length(on.dEdge.avail)))

# predict
on.dEdge.fr.pred <- predict(on.dEdge, 
                            on.dEdge.fr.newdata, 
                            exclude = "s(cluster)",
                            se.fit = T,
                            newdata.guaranteed = T)

# dfs for plotting
on.dEdge.fr.forPlot <- on.dEdge.fr.newdata |>
  
  mutate(pred = on.dEdge.fr.pred$fit,
         l90 = on.dEdge.fr.pred$fit - 1.645 * on.dEdge.fr.pred$se.fit,
         u90 = on.dEdge.fr.pred$fit + 1.645 * on.dEdge.fr.pred$se.fit,
         
         TRT = factor(TRT, 
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P")))
# plot
ggplot() +
  
  theme_fr_rss() +
  
  # population-level effect
  geom_rect(data = M.on |> filter(param == "dEdge"),
            aes(xmin = 3500,
                xmax = 10500,
                ymax = upp,
                ymin = low),
            fill = "gray97") +
  
  geom_hline(data = M.on |> filter(param == "dEdge"),
             aes(yintercept = mean),
             color = "darkgray") +
  
  # treatment-specific predictions
  geom_ribbon(data = on.dEdge.fr.forPlot,
              aes(x = avail,
                  y = pred,
                  fill = TRT,
                  ymin = l90,
                  ymax = u90),
              alpha = 0.25) +
  
  geom_line(data = on.dEdge.fr.forPlot,
            aes(x = avail,
                y = pred,
                color = TRT,
                linetype = TRT),
            linewidth = 1.2) +
  
  coord_cartesian(xlim = c(3900, 10000)) +
  
  theme(legend.position = "none",
        #legend.title = element_blank(),
        #legend.key.size = unit(0.4, "cm"),
        #legend.text = element_text(size = 6),
        #legend.background = element_rect(fill = NA),
        axis.title.y = element_blank()) +
  
  scale_color_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_fill_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  
  xlab("Mean stems/ha") +
  ylab(expression(beta)) -> on.dEdge.fr.plot

# ______________________________________________________________________________
# 7b. RSS predictions ----
# ______________________________________________________________________________

# filter data
on.dEdge.mean.sd <- mean.sd.on |> filter(name == "dEdge")
on.dEdge.mean.sd.trt <- mean.sd.on.trt |> filter(name == "dEdge")

# understory availability levels (just pull from fr df to join)
on.dEdge.avail.levels <- c(on.dEdge.fr.forPlot$avail[1],
                            on.dEdge.fr.forPlot$avail[50],
                            on.dEdge.fr.forPlot$avail[100])

# availability ranges
on.dEdge.range.trt <- list(seq(on.dEdge.mean.sd.trt$min[on.dEdge.mean.sd.trt$TRT == "UNTHIN"], 
                               on.dEdge.mean.sd.trt$max[on.dEdge.mean.sd.trt$TRT == "UNTHIN"], 
                               length.out = 100),
                           seq(on.dEdge.mean.sd.trt$min[on.dEdge.mean.sd.trt$TRT == "RET"], 
                               on.dEdge.mean.sd.trt$max[on.dEdge.mean.sd.trt$TRT == "RET"], 
                               length.out = 100),
                           seq(on.dEdge.mean.sd.trt$min[on.dEdge.mean.sd.trt$TRT == "PIL"], 
                               on.dEdge.mean.sd.trt$max[on.dEdge.mean.sd.trt$TRT == "PIL"], 
                               length.out = 100))

on.dEdge.range.trt.s <- list((on.dEdge.range.trt[[1]] - on.dEdge.mean.sd$mean) / on.dEdge.mean.sd$sd,
                             (on.dEdge.range.trt[[2]] - on.dEdge.mean.sd$mean) / on.dEdge.mean.sd$sd,
                             (on.dEdge.range.trt[[3]] - on.dEdge.mean.sd$mean) / on.dEdge.mean.sd$sd)

# RSS prediction
on.dEdge.rss.newdata <- bind_rows(
  
  data.frame(x = on.dEdge.range.trt.s[[1]],
             TRT = "UNTHIN"),
  data.frame(x = on.dEdge.range.trt.s[[2]],
             TRT = "RET"),
  data.frame(x = on.dEdge.range.trt.s[[3]],
             TRT = "PIL")
  
)

# repeat and add availability levels
on.dEdge.rss.newdata.1 <- bind_rows(on.dEdge.rss.newdata,
                                    on.dEdge.rss.newdata,
                                    on.dEdge.rss.newdata) |>
  
  # add availability levels
  mutate(avail = c(rep(on.dEdge.avail.levels[1], 300),
                   rep(on.dEdge.avail.levels[2], 300),
                   rep(on.dEdge.avail.levels[3], 300))) |>
  
  # treatment factor
  mutate(TRT = factor(TRT, 
                      levels = c("UNTHIN", "RET", "PIL"),
                      labels = c("U", "R", "P"))) |>
  
  left_join(on.dEdge.fr.forPlot |> dplyr::select(avail, TRT, pred, l90, u90)) |>
  
  # calculate log RSS
  mutate(rss.est = x * pred,
         rss.low = x * l90,
         rss.upp = x * u90) |>
  
  # back-transform x
  mutate(x.1 = (x * on.dEdge.mean.sd$sd) + on.dEdge.mean.sd$mean) |>
  
  # availability to sensible factor
  mutate(avail = factor(avail,
                        levels = unique(factor(avail)),
                        labels = c("Mean stems = 3580",
                                   "Mean stems = 6982",
                                   "Mean stems = 10454"))) |>
  
  # keep only relevant columns for plotting
  dplyr::select(x.1, TRT, avail, rss.est, rss.low, rss.upp)

# plot
ggplot(data = on.dEdge.rss.newdata.1) +
  
  theme_fr_rss() +
  
  facet_wrap(~ avail) +
  
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
        strip.background = element_blank(),
        strip.text = element_blank(),
        
        legend.title = element_blank(),
        legend.position = c(0.08, 0.75),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 6)) +
  
  geom_text(aes(x = 2,
                y = -0.25,
                label = avail),
            size = 2,
            hjust = 0,
            check_overlap = T) + 
  
  scale_color_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_fill_manual(values = c("gray65", "dodgerblue3", "dodgerblue4")) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  
  xlab("Distance to mature edge (m)") +
  ylab("log(RSS)") -> on.dEdge.rss.plot

# ______________________________________________________________________________
# 8. Plot together ----
# ______________________________________________________________________________

# define design
# T L B R
plot.design <- c(
  
  area(1, 1, 1, 1), area(1, 2, 1, 2.5),
  area(2, 1, 2, 1), area(2, 2, 2, 2.5),
  area(3, 1, 3, 1), area(3, 2, 3, 4.8),
  area(4, 1, 4, 1), area(4, 2, 4, 4.8),
  area(5, 1, 5, 1), area(5, 2, 5, 4.8)
  
)

off.vo.fr.plot + off.vo.rss.plot + 
off.cc.fr.plot + off.cc.rss.plot + 
off.dEdge.fr.plot + off.dEdge.rss.plot + 
on.cc.fr.plot + on.cc.rss.plot +
on.dEdge.fr.plot + on.dEdge.rss.plot +
plot_layout(design = plot.design)

# 520 x 788