# PROJECT: Habitat selection
# SCRIPT: 11 - Multivariate habitat selection
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 23 Jul 2026
# COMPLETED: 23 Jul 2026
# LAST MODIFIED: 23 Jul 2026
# R VERSION: 4.5.2

# instead of testing relationships for single coefficients, it might be useful
# to show multivariate patterns (i.e., an ordination of all selection coefficients)
# by season, and just see if we can pick out patterns between the treatments

# PREDICTION: you won't be able to tell in two dimensions. Differences are 
# more or less random when you consider all responses together

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(vegan)

# ______________________________________________________________________________
# 2. Read individual slopes ----
# ______________________________________________________________________________

fr.data.off <- readRDS("data_for_model/off_fr.rds")
fr.data.on <- readRDS("data_for_model/on_fr.rds")

# ______________________________________________________________________________
# 3. Clean ----

# we need a n x beta "community data matrix" for vegan

# ______________________________________________________________________________

comm.off <- fr.data.off |>
  
  dplyr::select(TSPID, param, beta) |>
  
  pivot_wider(names_from = param, values_from = beta) |>
  
  # drop TSPID and G.s
  dplyr::select(-c(TSPID, g.s))

comm.on <- fr.data.on |>
  
  dplyr::select(TSPID, param, beta) |>
  
  pivot_wider(names_from = param, values_from = beta) |>
  
  # drop TSPID and G.s
  dplyr::select(-c(TSPID, g.s))

# ______________________________________________________________________________
# 4. Fit NMDS ----
# ______________________________________________________________________________
# 4a. Off ----
# ______________________________________________________________________________

nmds.off <- metaMDS(
  
  comm.off,
  distance = "euclidean",  # Bray doesn't work with negative values
  k = 2,
  autotransform = F,
  try = 40
  
)

plot(nmds.off)

# ______________________________________________________________________________
# 4b. On ----
# ______________________________________________________________________________

nmds.on <- metaMDS(
  
  comm.on,
  distance = "euclidean",  # Bray doesn't work with negative values
  k = 2,
  autotransform = F,
  try = 40
  
)

plot(nmds.on)

# ______________________________________________________________________________
# 5. Evaluate NMDS ----
# ______________________________________________________________________________

# stress
nmds.off$stress
nmds.on$stress

# goodness of fit
hist(goodness(nmds.off))
hist(goodness(nmds.on))

# distance vs dissimilarity
plot(nmds.off$diss, nmds.off$dist)
plot(nmds.on$diss, nmds.on$dist)

# stressplot
stressplot(nmds.off)
stressplot(nmds.on)

# ______________________________________________________________________________
# 6. Extract scores for plotting ----
# ______________________________________________________________________________

scores.off <- fr.data.off |>
  
  group_by(TSPID) |>
  
  slice(1) |>
  
  dplyr::select(TSPID, sex, TRT, cluster) |>
  
  ungroup() |>
  
  mutate(NMDS1 = scores(nmds.off)[ , 1],
         NMDS2 = scores(nmds.off)[ , 2],
         TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")))
  
scores.on <- fr.data.on |>
  
  group_by(TSPID) |>
  
  slice(1) |>
  
  dplyr::select(TSPID, sex, TRT, cluster) |>
  
  ungroup() |>
  
  mutate(NMDS1 = scores(nmds.on)[ , 1],
         NMDS2 = scores(nmds.on)[ , 2],
         TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")))

# ______________________________________________________________________________
# 7. Plot function ----
# ______________________________________________________________________________

plot_nmds <- function (.scores) {
  
  if (nrow(.scores) == 97) { 
    
    .col <- c("gray", "green3", "green4") 
    
    } else { 
    
    .col <- c("gray", "dodgerblue3", "dodgerblue4")
    
    }
  
  ggplot(.scores) +
    
    theme_classic() +
    
    stat_ellipse(aes(x = NMDS1,
                     y = NMDS2,
                     linetype = TRT,
                     color = TRT),
                 linewidth = 1.1) +
    
    geom_point(aes(x = NMDS1,
                   y = NMDS2,
                   shape = TRT,
                   fill = TRT),
               color = "black") +
    
    scale_fill_manual(values = .col) +
    scale_color_manual(values = .col) +
    scale_shape_manual(values = c(21, 22, 23)) +
    
    theme(legend.position = "top",
          legend.title = element_blank())
  
}

# plot together
cowplot::plot_grid(plot_nmds(scores.off), plot_nmds(scores.on))
