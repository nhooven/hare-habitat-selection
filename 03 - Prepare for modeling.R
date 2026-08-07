# PROJECT: Habitat selection
# SCRIPT: 03 - Prepare for modeling
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 01 Jun 2026
# LAST MODIFIED: 06 Aug 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(mefa4)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

data.off <- readRDS("data_cleaned/data_off.rds")
data.on <- readRDS("data_cleaned/data_on.rds")

# ______________________________________________________________________________
# 3. Prepare for modeling ----
# ______________________________________________________________________________
# 3a. Correct values for each treatment, drop NAs ----
# ______________________________________________________________________________

prep_1 <- function (x) {
  
  x.1 <- x |> 
    
    # attribute correct values
    mutate(
      
      # CONDITIONS
      cc = case_when(year == "PRE" ~ cc.pre,
                     year %in% c("POST1", "POST2") ~ cc.post),
      
      # LOCAL
      stem = case_when(year == "PRE" ~ stem.pre,
                       year %in% c("POST1", "POST2") ~ stem.post),
      vo = case_when(year == "PRE" ~ vo.pre,
                     year %in% c("POST1", "POST2") ~ vo.post),
      ch = case_when(year == "PRE" ~ ch.pre,
                     year %in% c("POST1", "POST2") ~ ch.post)
      
    ) |>
    
    # drop variables
    dplyr::select(
      
      -c(trt,
         season,
         cc.pre,
         cc.post,
         stem.pre,
         stem.post,
         vo.pre,
         vo.post,
         ch.pre,
         ch.post)
      
    ) |>
    
    # drop NAs
    drop_na(c(cc, twi, vrm,
              stem, vo, ch, dEdge)) |>
    
    # case weights
    mutate(w = ifelse(case == 0, 5000, 1)) |>
    
    # keep correct variables
    dplyr::select(
      
      track_season_post,   # for calculating availability
      year,
      c.trt,
      case,
      w,
      
      # predictors
      akde,
      
      # CONDITIONS
      cc, twi, vrm,
      
      # FOCAL
      stem, vo, ch, dEdge
      
    )
  
  return(x.1)
  
}

# use
data.off.1 <- prep_1(data.off)
data.on.1 <- prep_1(data.on)

# ______________________________________________________________________________
# 4. Remove any tracks with < 10 relocations ----
# ______________________________________________________________________________

n.used.off <- data.off.1 |> group_by(track_season_post) |>
  
  filter(case == 1) |>
  
  summarize(n.used = n())

n.used.on <- data.on.1 |> group_by(track_season_post) |>
  
  filter(case == 1) |>
  
  summarize(n.used = n())

# remove anything with < 10 relocations
off.less10 <- n.used.off$track_season_post[n.used.off$n.used < 10]  # n = 8
on.less10 <- n.used.on$track_season_post[n.used.on$n.used < 10] # n = 2

data.off.2 <- data.off.1 |> filter(track_season_post %notin% off.less10)
data.on.2 <- data.on.1 |> filter(track_season_post %notin% on.less10)

# ______________________________________________________________________________
# 5. Extract attributes to bind in later ----
# ______________________________________________________________________________

data.off.attr <- data.off.2 |> dplyr::select(track_season_post:akde)
data.on.attr <- data.on.2 |> dplyr::select(track_season_post:akde)

data.off.cov <- data.off.2 |> dplyr::select(cc:dEdge)
data.on.cov <- data.on.2 |> dplyr::select(cc:dEdge)

# ______________________________________________________________________________
# 6. Transformations ----

# function
transform_covs <- function (x) {
  
  x.1 <- x |>
    
    # covariate transformations
    # squared
    mutate(twi2 = twi^2,
           vrm2 = vrm^2)
  
  return(x.1)
  
}

# ______________________________________________________________________________

data.off.3 <- transform_covs(data.off.cov)
data.on.3 <- transform_covs(data.on.cov)

# ______________________________________________________________________________
# 7. Save means, SDs, and ranges ----

# function
mean_sds <- function (x) {
  
  x.1 <- x |>
    
    pivot_longer(cols = everything()) |>
    
    group_by(name) |>
    
    summarize(
      
      mean = mean(value),
      sd = sd(value),
      min = min(value),
      max = max(value)
      
    )
  
  return(x.1)
  
}

mean_sds_trt <- function (x, y) {
  
  x.1 <- x |>
    
    bind_cols(y |> dplyr::select(year, c.trt)) |>
    
    mutate(TRT = case_when(
      
      year == "PRE" ~ "UNTHIN",
      year %in% c("POST1", "POST2") & c.trt == "CTRL" ~ "UNTHIN",
      year %in% c("POST1", "POST2") & c.trt == "RET" ~ "RET",
      year %in% c("POST1", "POST2") & c.trt == "PIL" ~ "PIL"
      
    )
    
    ) |>
    
    dplyr::select(-c(year, c.trt)) |>
    
    pivot_longer(cols = cc:vrm2) |>
    
    group_by(name, TRT) |>
    
    summarize(
      
      mean = mean(value),
      sd = sd(value),
      min = min(value),
      max = max(value)
      
    ) |>
    
    ungroup()
  
  return(x.1)
  
}

# ______________________________________________________________________________

# use
mean.sd.off <- mean_sds(data.off.3)
mean.sd.on <- mean_sds(data.on.3)

mean.sd.off.trt <- mean_sds_trt(data.off.3, data.off.attr)
mean.sd.on.trt <- mean_sds_trt(data.on.3, data.on.attr)

# ______________________________________________________________________________
# 8. Standardize ----

# function
standardize_across <- function (x) {
  
  # standardize function
  standardize <- function (x) {
    
    x.1 <- (x - mean(x)) / sd(x)
    
    return(x.1)
    
  }
  
  x.1 <- x |>
    
    mutate(
      
      across(everything(),
             standardize)
      
    )
  
  return(x.1)
  
}

# ______________________________________________________________________________

data.off.4 <- standardize_across(data.off.3)
data.on.4 <- standardize_across(data.on.3)

# ______________________________________________________________________________
# 9. Bind back in ----
# ______________________________________________________________________________

data.off.5 <- cbind(data.off.attr, data.off.4)
data.on.5 <- cbind(data.on.attr, data.on.4)

# ______________________________________________________________________________
# 10. Save to files ----
# ______________________________________________________________________________

saveRDS(data.off.5, "data_for_model/off_data.rds")
saveRDS(data.on.5, "data_for_model/on_data.rds")

saveRDS(mean.sd.off, "data_for_model/mean_sd_off.rds")
saveRDS(mean.sd.on, "data_for_model/mean_sd_on.rds")

saveRDS(mean.sd.off.trt, "data_for_model/mean_sd_off_trt.rds")
saveRDS(mean.sd.on.trt, "data_for_model/mean_sd_on_trt.rds")

# ______________________________________________________________________________
# 11. Summaries for results ----
# ______________________________________________________________________________

# n tracks
length(unique(data.off.5$track_season_post))
length(unique(data.on.5$track_season_post))

# total relocations
data.off.relocs <- data.off.5 |> filter(case == 1) |>
  
  group_by(track_season_post) |>
  
  summarize(relocs = n())

data.on.relocs <- data.on.5 |> filter(case == 1) |>
  
  group_by(track_season_post) |>
  
  summarize(relocs = n())

all.relocs <- rbind(data.off.relocs, data.on.relocs)

# summaries
min(all.relocs$relocs)
max(all.relocs$relocs)
median(all.relocs$relocs)
sd(all.relocs$relocs)

# total relocs and background
nrow(data.off.5)
nrow(data.off.5 |> filter(case == 1))

nrow(data.on.5)
nrow(data.on.5 |> filter(case == 1))
