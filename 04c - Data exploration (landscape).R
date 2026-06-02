# PROJECT: Habitat selection
# SCRIPT: 04c - Data exploration (landscape)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 
# LAST MODIFIED: 02 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(hoove)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

data.off <- readRDS("data_for_model/off_data.rds")
data.on <- readRDS("data_for_model/on_data.rds")

# ______________________________________________________________________________
# 3. "Additive habitat use" plots ----

# Holbrook et al 2019
# doi.org/10.1002/eap.1852

# function
# generic function
plot_hab_use <- function (.covar) {
  
  # subset data and calculate individual-specific values
  focal.data.off <- data.off |> dplyr::select(c(track_season_post, case, .covar, year.trt, ret, pil))
  focal.data.on <- data.on |> dplyr::select(c(track_season_post, case, .covar, year.trt, ret, pil))
  
  colnames(focal.data.off)[3] <- "var"
  colnames(focal.data.on)[3] <- "var"
  
  # calculate individual-specific values
  focal.data <- focal.data.off |>
    
    # add season
    mutate(season = "off") |>
    
    # add treatment
    mutate(trt = case_when(
      
      year.trt == 0 ~ "unthinned",
      ret == 0 & pil == 0 & year.trt == 1 ~ "unthinned",
      ret == 1 & year.trt == 1 ~ "ret",
      pil == 1 & year.trt == 1 ~ "pil"
      
    )
    
    ) |>
    
    # bind
    bind_rows(
      
      focal.data.on |>
        
        # add season
        mutate(season = "on") |>
        
        # add treatment
        mutate(trt = case_when(
          
          year.trt == 0 ~ "unthinned",
          ret == 0 & pil == 0 & year.trt == 1 ~ "unthinned",
          ret == 1 & year.trt == 1 ~ "ret",
          pil == 1 & year.trt == 1 ~ "pil"
          
        )
        
        )
      
    ) |>
    
    group_by(track_season_post, trt, season, case) |>
    
    summarize(mean = mean(var),
              sd = sd(var)) |>
    
    pivot_wider(names_from = case,
                values_from = c(mean, sd)) |>
    
    rename("available" = mean_0,
           "used" = mean_1,
           "available_sd" = sd_0,
           "used_sd" = sd_1)

  # plot
  out.plot <- ggplot(focal.data) +
    
    theme_bw() +
    
    facet_grid(trt ~ season) +
    
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed") +
    
    geom_point(aes(x = available,
                   y = used,
                   color = season),
               size = 0.5) +
    
    geom_smooth(aes(x = available,
                    y = used,
                    color = season,
                    fill = season),
                method = "gam") +
    
    theme(panel.grid = element_blank(),
          strip.background = element_rect(color = NA),
          strip.text = element_text(hjust = 0),
          legend.position = "none") +
    
    scale_color_manual(values = c("#FF3300", "lightblue4")) +
    
    ggtitle(.covar)
  
  out.plot
  
}

# ______________________________________________________________________________

plot_hab_use("dOpen")      # relaxed selection (i.e., relax avoidance of cover)
plot_hab_use("dDM")        # relaxed selection (i.e., relax avoidance of cover)
plot_hab_use("ed")         # trade-off

# looks like we could include FRs for all of these