# PROJECT: Habitat selection
# SCRIPT: 08a - Visualization (parameter estimates)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 
# LAST MODIFIED: 27 May 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(hoove)
library(cowplot)

# ______________________________________________________________________________
# 2. Read in model results ----
# ______________________________________________________________________________

M1.off <- readRDS("model_results/M1_off.rds")

# ______________________________________________________________________________
# 3. Process ----
# ______________________________________________________________________________

# reorder/rename factors
param.levels <- M1.off[[1]]$param
param.labels <- c("visual obstruction", "canopy height", "canopy height2",
                  "canopy cover", "canopy cover2",
                  "wetness", "wetness2", "ruggedness", "ruggedness2",
                  "log(distance open)", "log(distance dense mature)")

M1.off[[1]] <- M1.off[[1]] |>
  
  mutate(param = factor(param,
                        levels = param.levels,
                        labels = param.labels))

M1.off[[2]] <- M1.off[[2]] |>
  
  mutate(param = factor(param,
                        levels = param.levels,
                        labels = param.labels))

# add variable class

# ______________________________________________________________________________
# 4. Dual parameter estimate and RS variance plot ----
# ______________________________________________________________________________

plot.param <- ggplot(data = M1.off[[1]]) +
  
  theme_hoove() +
  
  # vline
  geom_vline(xintercept = 0,
             linetype = "dashed") + 
  
  # CIs
  geom_errorbar(aes(x = mean,
                    y = param,
                    xmin = low,
                    xmax = upp),
                width = 0,
                linewidth = 2,
                color = "gray50") +
  
  # means
  geom_point(aes(x = mean,
                 y = param),
             shape = 21,
             color = "black",
             fill = "white",
             size = 2,
             stroke = 1.5) +
  
  theme(axis.title.y = element_blank()) +
  
  xlab("Standardized coefficient") +
  
  scale_y_discrete(limits = rev)

plot.param

plot.var <- ggplot(M1.off[[2]]) +
  
  theme_hoove() +
  
  geom_col(aes(x = variance,
               y = param),
           color = "black",
           fill = "#FF3300",
           alpha = 0.5) +
  
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  
  xlab("Random slope variance") +
  
  scale_y_discrete(limits = rev)

plot.var

# plot together
plot_grid(plot.param, plot.var,
          rel_widths = c(1, 0.6))
