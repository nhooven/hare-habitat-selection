# PROJECT: Habitat selection
# SCRIPT: 08a - Visualization (parameter estimates)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 11 Jun 2026
# LAST MODIFIED: 18 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(hoove)
library(cowplot)

# ______________________________________________________________________________
# 2. Read in cleaned model results ----
# ______________________________________________________________________________

M.off <- readRDS("model_results/M_off.rds")
M.on <- readRDS("model_results/M_on.rds")

# ______________________________________________________________________________
# 3. Re-order and rename factor levels ----
# ______________________________________________________________________________
# 3a. Define order and labels ----
# ______________________________________________________________________________

# levels 
levels.off <- M.off[[1]]$param
levels.on <- M.on[[1]]$param

# labels
labels.off <- c("log(G[s])", 
                "canopy cover", "canopy cover2",
                "wetness", "wetness2", "ruggedness", "ruggedness2",
                "visual obstruction", "canopy height", "distance to edge")

labels.on <- c("log(G[s])", 
               "canopy cover", "canopy cover2",
               "wetness", "wetness2", "ruggedness", "ruggedness2",
               "stem density", "canopy height", "distance to edge")

# ______________________________________________________________________________
# 3b. Function ----
# ______________________________________________________________________________

apply_cov_labels <- function (.results,
                              .levels,
                              .labels) {
  
  .results.1 <- .results |>
    
    mutate(param = factor(param,
                          levels = .levels,
                          labels = .labels))
  
  return(.results.1)
  
}

# use on population-level parameter estimates
off.pop <- apply_cov_labels(M.off[[1]], levels.off, labels.off)
on.pop <- apply_cov_labels(M.on[[1]], levels.on, labels.on)

# use on random slope variances
off.rs.var <- apply_cov_labels(M.off[[2]], levels.off, labels.off)
on.rs.var <- apply_cov_labels(M.on[[2]], levels.on, labels.on)

# ______________________________________________________________________________
# 4. Population-level effect plot ----
# ______________________________________________________________________________

# function
plot_pop <- function (.pop) {
  
  ggplot(data = .pop) +
    
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
  
}

# use function
plot_pop(off.pop)

# ______________________________________________________________________________
# 5. Random slope variance plot ----
# ______________________________________________________________________________

plot_var <- function (.rs.var) {
  
  ggplot(data = .rs.var) +
    
    theme_hoove() +
    
    geom_col(aes(x = variance,
                 y = param),
             color = "black",
             fill = "#FF3300",
             alpha = 0.5) +
    
    theme(axis.title.y = element_blank()) +
    
    xlab("Random slope variance") +
    
    scale_y_discrete(limits = rev)
  
}

# use
plot_var(off.rs.var)

# ______________________________________________________________________________
# 6. Both seasons together ----
# ______________________________________________________________________________
# 6a. Prep dfs ----
# ______________________________________________________________________________

# population level effects
both.pop <- bind_rows(
  
  off.pop |> mutate(season = "snow-off"),
  on.pop |> mutate(season = "snow-on")
  
) |>
  
  # remove log(g(s)) for now
  filter(param != "log(G[s])") |>
  
  # groupings
  mutate(
    
    group = case_when(
      
      param %in% c("visual obstruction",
                   "stem density",
                   "canopy height",
                   "distance to edge") ~ "focal",
      param %in% c("canopy cover",
                   "canopy cover2",
                   "wetness",
                   "wetness2",
                   "ruggedness",
                   "ruggedness2") ~ "conditions"
      
    )
    
  ) |>
  
  # factor order
  mutate(
    
    param = factor(param, levels = c("canopy cover",
                                     "canopy cover2",
                                     "wetness",
                                     "wetness2",
                                     "ruggedness",
                                     "ruggedness2",
                                     "visual obstruction",
                                     "stem density",
                                     "canopy height",
                                     "distance to edge"),
                   labels = c("canopy cover",
                              "canopy cover (sq)",
                              "wetness",
                              "wetness (sq)",
                              "ruggedness",
                              "ruggedness (sq)",
                              "visual obstruction",
                              "stem density",
                              "canopy height",
                              "distance to edge")),
    
    group = factor(group, levels = c("focal", "conditions"))
    
  )

# random slopes
both.rs <- bind_rows(
  
  M.off[[3]] |> pivot_longer(cols = cc:dEdge) |>
    
    rename(param = name, mean = value) |>
    
    mutate(season = "snow-off") |>
    
    apply_cov_labels(levels.off, labels.off),
  
  M.on[[3]] |> pivot_longer(cols = cc:dEdge) |>
    
    rename(param = name, mean = value) |>
    
    mutate(season = "snow-on") |>
    
    apply_cov_labels(levels.on, labels.on)
  
) |>
  
  # drop log(g(s))
  filter(param != "log(G[s])") |>
  
  # groupings
  mutate(
    
    group = case_when(
      
      param %in% c("canopy cover",
                   "canopy cover2",
                   "wetness",
                   "wetness2",
                   "ruggedness",
                   "ruggedness2") ~ "conditions",
      param %in% c("visual obstruction",
                   "stem density",
                   "canopy height",
                   "distance to edge") ~ "focal"
      
    )
    
  ) |>
  
  # factor order
  mutate(
    
    param = factor(param, levels = c("canopy cover",
                                     "canopy cover2",
                                     "wetness",
                                     "wetness2",
                                     "ruggedness",
                                     "ruggedness2",
                                     "visual obstruction",
                                     "stem density",
                                     "canopy height",
                                     "distance to edge"),
                   labels = c("canopy cover",
                              "canopy cover (sq)",
                              "wetness",
                              "wetness (sq)",
                              "ruggedness",
                              "ruggedness (sq)",
                              "visual obstruction",
                              "stem density",
                              "canopy height",
                              "distance to edge")),
    
    group = factor(group, levels = c("focal", "conditions"))
    
  )

# ______________________________________________________________________________
# 6b. Plot ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  facet_wrap(~ group,
             scales = "free_y",
             nrow = 2,
             strip.position = "right",
             space = "free_y") +      # ggplot says I can't use this, but it works anyway
  
  # vertical line
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  # random slopes
  geom_point(data = both.rs,
             aes(x = mean,
                 y = param,
                 color = season,
                 group = season),
             size = 2,  
             position = position_dodge(width = 0.75),
             alpha = 0.1) +
  
  # population-level effects
  geom_errorbar(data = both.pop,
                aes(x = mean,
                    xmin = low,
                    xmax = upp,
                    y = param,
                    group = season),
                position = position_dodge(width = 0.75),
                width = 0,
                linewidth = 1) +
  
  geom_point(data = both.pop,
             aes(x = mean,
                 y = param,
                 group = season,
                 fill = season),
             size = 1.4,   # slightly smaller
             shape = 23,
             color = "white",
             position = position_dodge(width = 0.75)) +
  
  # theme arguments
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        
        legend.position = c(0.82, 0.15),
        legend.background = element_rect(fill = "white",
                                         color = "gray"),
        legend.title = element_blank(),
        legend.margin = margin(1, 10, 1, 0.5),  # TRBL
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black",
                                 size = 8),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA)) +
  
  # legend guide
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  
  # axis range
  coord_cartesian(xlim = c(-2.2, 2.2)) +
  
  # labels
  xlab("Selection coefficient") +
  
  # colors
  scale_color_manual(values = c("green4", "dodgerblue2")) +
  scale_fill_manual(values = c("green4", "dodgerblue2")) +
  
  scale_y_discrete(limits = rev)

# 465 x 481
