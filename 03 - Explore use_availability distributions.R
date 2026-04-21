# PROJECT: Habitat selection
# SCRIPT: 03 - Explore use/availability distributions
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 
# LAST MODIFIED: 21 Apr 2026
# R VERSION: 4.4.3

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----

# RDS for now

# ______________________________________________________________________________

data.off <- readRDS("data_cleaned/data_off.rds")
data.on <- readRDS("data_cleaned/data_on.rds")

# ______________________________________________________________________________
# 3. Distribution plot function ----
# ______________________________________________________________________________

dist_plot <- function (df, var) {
  
  # subset
  df.1 <- df %>%
    
    dplyr::select(case,
                  all_of(var)) %>%
    
    # rename for generality
    rename(covariate = var)
  
  # plot
  out.plot <- ggplot() +
    
    theme_classic() +
    
    geom_density(data = df.1,
                 aes(x = covariate,
                     color = as.factor(case),
                     fill = as.factor(case),
                     linetype = as.factor(case)),
                 linewidth = 0.9,
                 alpha = 0.25) + 
    
    scale_color_manual(values = c("lightgray", "darkgreen")) +
    scale_fill_manual(values = c("lightgray", "darkgreen")) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    
    labs(x = var)
  
  return(out.plot)
  
}

# ______________________________________________________________________________
# 4. Snow off ----
# ______________________________________________________________________________

# distance
dist_plot(data.off, "dEdge")
dist_plot(data.off, "dDh")
dist_plot(data.off, "dMsM")
dist_plot(data.off, "dJsM")
dist_plot(data.off, "dOM")

# canopy
dist_plot(data.off, "ch.pre")
dist_plot(data.off, "ch.post")
dist_plot(data.off, "cc.pre")
dist_plot(data.off, "cc.post")

# veg models 
dist_plot(data.off, "stem.pre")
dist_plot(data.off, "stem.post")
dist_plot(data.off, "vo.pre")
dist_plot(data.off, "vo.post")
dist_plot(data.off, "shrub")

# topography
dist_plot(data.off, "twi")
dist_plot(data.off, "slope")
dist_plot(data.off, "tpi")
dist_plot(data.off, "vrm")

# ______________________________________________________________________________
# 5. Snow on ----
# ______________________________________________________________________________

# distance
dist_plot(data.on, "dEdge")
dist_plot(data.on, "dDh")
dist_plot(data.on, "dMsM")
dist_plot(data.on, "dJsM")
dist_plot(data.on, "dOM")

# canopy
dist_plot(data.on, "ch.pre")
dist_plot(data.on, "ch.post")
dist_plot(data.on, "cc.pre")
dist_plot(data.on, "cc.post")

# veg models 
dist_plot(data.on, "stem.pre")
dist_plot(data.on, "stem.post")
dist_plot(data.on, "vo.pre")
dist_plot(data.on, "vo.post")
dist_plot(data.on, "shrub")

# topography
dist_plot(data.on, "twi")
dist_plot(data.on, "slope")
dist_plot(data.on, "tpi")
dist_plot(data.on, "vrm")

# ______________________________________________________________________________
# 6. Treatment variables ----
# ______________________________________________________________________________
# 6a. Snow-off ----
# ______________________________________________________________________________

# dPil
ggplot(
  
  data = data.off %>% filter(trt == "PIL" & year %in% c("POST1", "POST2")) 
  
  ) +
  
  facet_wrap(~ year) +
  
  theme_classic() +
  
  geom_density(aes(x = dPil,
                   color = as.factor(case),
                   fill = as.factor(case),
                   linetype = as.factor(case)),
               linewidth = 0.6,
               alpha = 0.25) + 
  
  scale_color_manual(values = c("lightgray", "darkgreen")) +
  scale_fill_manual(values = c("lightgray", "darkgreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  
  coord_cartesian(xlim = c(0, 400))

# dRet
ggplot(
  
  data = data.off %>% filter(trt == "RET" & year %in% c("POST1", "POST2")) 
  
) +
  
  facet_wrap(~ year) +
  
  theme_classic() +
  
  geom_density(aes(x = dRet,
                   color = as.factor(case),
                   fill = as.factor(case),
                   linetype = as.factor(case)),
               linewidth = 0.6,
               alpha = 0.25) + 
  
  scale_color_manual(values = c("lightgray", "darkgreen")) +
  scale_fill_manual(values = c("lightgray", "darkgreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  
  coord_cartesian(xlim = c(0, 450))

# ______________________________________________________________________________
# 6b. Snow-on ----
# ______________________________________________________________________________

# dPil
ggplot(
  
  data = data.on %>% filter(trt == "PIL" & year %in% c("POST1", "POST2")) 
  
) +
  
  facet_wrap(~ year) +
  
  theme_classic() +
  
  geom_density(aes(x = dPil,
                   color = as.factor(case),
                   fill = as.factor(case),
                   linetype = as.factor(case)),
               linewidth = 0.6,
               alpha = 0.25) + 
  
  scale_color_manual(values = c("lightgray", "darkgreen")) +
  scale_fill_manual(values = c("lightgray", "darkgreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  
  coord_cartesian(xlim = c(0, 400))

# dRet
ggplot(
  
  data = data.on %>% filter(trt == "RET" & year %in% c("POST1", "POST2")) 
  
) +
  
  facet_wrap(~ year) +
  
  theme_classic() +
  
  geom_density(aes(x = dRet,
                   color = as.factor(case),
                   fill = as.factor(case),
                   linetype = as.factor(case)),
               linewidth = 0.6,
               alpha = 0.25) + 
  
  scale_color_manual(values = c("lightgray", "darkgreen")) +
  scale_fill_manual(values = c("lightgray", "darkgreen")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  
  coord_cartesian(xlim = c(0, 450))
