# PROJECT: Habitat selection
# SCRIPT: FIG - Functional response predictions
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 28 Jul 2026
# COMPLETED: 
# LAST MODIFIED: 28 Jul 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

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
  
  theme_classic() +
  
  # population-level effect
  geom_rect(data = M.off |> filter(param == "vo"),
            aes(xmin = 0.4,
                xmax = 3.6,
                ymax = upp,
                ymin = low),
            fill = "gray95") +
  
  geom_hline(data = M.off |> filter(param == "vo"),
             aes(yintercept = mean)) +
  
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
             size = 2.5) +
  
  theme(legend.position = "none") +
  
  scale_color_manual(values = c("gray65", "green4", "green4")) +
  scale_fill_manual(values = c("gray65", "green4", "green4")) +
  
  xlab("Treatment") +
  ylab(expression(beta)) -> off.vo.fr.plot

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
  
  theme_classic() +
  
  # population-level effect
  geom_rect(data = M.off |> filter(param == "cc"),
            aes(xmin = 0.4,
                xmax = 3.6,
                ymax = upp,
                ymin = low),
            fill = "gray95") +
  
  geom_hline(data = M.off |> filter(param == "cc"),
             aes(yintercept = mean)) +
  
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
             size = 2.5) +
  
  theme(legend.position = "none") +
  
  scale_color_manual(values = c("gray65", "green4", "green4")) +
  scale_fill_manual(values = c("gray65", "green4", "green4")) +
  
  xlab("Treatment") +
  ylab(expression(beta)) -> off.cc.fr.plot

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
  
  theme_classic() +
  
  # population-level effect
  geom_rect(data = M.off |> filter(param == "dEdge"),
            aes(xmin = 0.51,
                xmax = 0.75,
                ymax = upp,
                ymin = low),
            fill = "gray95") +
  
  geom_hline(data = M.off |> filter(param == "dEdge"),
             aes(yintercept = mean)) +
  
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
  
  coord_cartesian(xlim = c(0.52, 0.72)) +
  
  xlab("Mean visual obstruction in UD") +
  ylab(expression(beta)) -> off.dEdge.fr.plot

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
  
  theme_classic() +
  
  # population-level effect
  geom_rect(data = M.on |> filter(param == "cc"),
            aes(xmin = 3500,
                xmax = 10500,
                ymax = upp,
                ymin = low),
            fill = "gray95") +
  
  geom_hline(data = M.on |> filter(param == "cc"),
             aes(yintercept = mean)) +
  
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
  
  xlab("Mean stem density in UD") +
  ylab(expression(beta)) -> on.cc.fr.plot

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
  
  theme_classic() +
  
  # population-level effect
  geom_rect(data = M.on |> filter(param == "dEdge"),
            aes(xmin = 3500,
                xmax = 10500,
                ymax = upp,
                ymin = low),
            fill = "gray95") +
  
  geom_hline(data = M.on |> filter(param == "dEdge"),
             aes(yintercept = mean)) +
  
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
  
  theme(legend.position = c(0.77, 0.87),
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  
  scale_color_manual(values = c("gray65", "dodgerblue3", "dodgerblue3")) +
  scale_fill_manual(values = c("gray65", "dodgerblue3", "dodgerblue3")) +
  scale_linetype_manual(values = c("dotted", "solid", "solid")) +
  
  xlab("Mean stem density in UD") +
  ylab(expression(beta)) -> on.dEdge.fr.plot
