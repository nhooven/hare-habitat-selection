# PROJECT: Habitat selection
# SCRIPT: 03 - Prepare for modeling
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 01 Jun 2026
# LAST MODIFIED: 09 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

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
      
      # CH
      ch = case_when(
        
        year == "PRE" ~ ch.pre,
        year %in% c("POST1", "POST2") ~ ch.post
        
      ),
      
      # CC
      cc = case_when(
        
        year == "PRE" ~ cc.pre,
        year %in% c("POST1", "POST2") ~ cc.post
        
      ),
      
      # VO
      vo = case_when(
        
        year == "PRE" ~ vo.pre,
        year %in% c("POST1", "POST2") ~ vo.post
        
      ),
      
      # stem
      stem = case_when(
        
        year == "PRE" ~ stem.pre,
        year %in% c("POST1", "POST2") ~ stem.post
        
      )
      
    ) |>
    
    # drop variables
    dplyr::select(
      
      -c(trt,
         season,
         ch.pre,
         ch.post,
         cc.pre,
         cc.post,
         vo.pre,
         vo.post)
      
    ) |>
    
    # drop NAs
    drop_na(
      
      c(
        
        dOpen,
        dDM,
        ch,
        cc,
        ed,
        vo,
        stem,
        twi,
        vrm
        
      )
      
    ) |>
    
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
      
      # landscape covariates
      vo,
      stem,
      ch,
      cc,
      twi,
      vrm,
      dOpen,
      dDM,
      ed
      
    )
  
  return(x.1)
  
}

# use
data.off.1 <- prep_1(data.off)
data.on.1 <- prep_1(data.on)

# ______________________________________________________________________________
# 5. Extract attributes to bind in later ----
# ______________________________________________________________________________

data.off.attr <- data.off.1 |> dplyr::select(track_season_post:akde)
data.on.attr <- data.on.1 |> dplyr::select(track_season_post:akde)

data.off.cov <- data.off.1 |> dplyr::select(vo:ed)
data.on.cov <- data.on.1 |> dplyr::select(vo:ed)

# ______________________________________________________________________________
# 6. Transformations ----

# function
transform_covs <- function (x) {
  
  x.1 <- x |>
    
    # covariate transformations
    # squared
    mutate(
      
      cc2 = cc^2,
      twi2 = twi^2,
      vrm2 = vrm^2
      
    ) 
  
  return(x.1)
  
}

# ______________________________________________________________________________

data.off.2 <- transform_covs(data.off.cov)
data.on.2 <- transform_covs(data.on.cov)

# ______________________________________________________________________________
# 6. Save means, SDs, and ranges ----

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
    
    pivot_longer(cols = vo:vrm2) |>
    
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
mean.sd.off <- mean_sds(data.off.2)
mean.sd.on <- mean_sds(data.on.2)

mean.sd.off.trt <- mean_sds_trt(data.off.2, data.off.attr)
mean.sd.on.trt <- mean_sds_trt(data.on.2, data.on.attr)

# ______________________________________________________________________________
# 7. Standardize ----

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

data.off.3 <- standardize_across(data.off.2)
data.on.3 <- standardize_across(data.on.2)

# ______________________________________________________________________________
# 8. Bind back in ----
# ______________________________________________________________________________

data.off.4 <- cbind(data.off.attr, data.off.3)
data.on.4 <- cbind(data.on.attr, data.on.3)

# ______________________________________________________________________________
# 9. Examine n used locations per TSP ----
# ______________________________________________________________________________

n.used.off <- data.off.4 |> group_by(track_season_post) |>
  
  filter(case == 1) |>
  
  summarize(n.used = n())

n.used.on <- data.on.4 |> group_by(track_season_post) |>
  
  filter(case == 1) |>
  
  summarize(n.used = n())

# some of these would never work with individual-level models

# ______________________________________________________________________________
# 10. Save to files ----
# ______________________________________________________________________________

saveRDS(data.off.4, "data_for_model/off_data.rds")
saveRDS(data.on.4, "data_for_model/on_data.rds")

saveRDS(mean.sd.off, "data_for_model/mean_sd_off.rds")
saveRDS(mean.sd.on, "data_for_model/mean_sd_on.rds")

saveRDS(mean.sd.off.trt, "data_for_model/mean_sd_off_trt.rds")
saveRDS(mean.sd.on.trt, "data_for_model/mean_sd_on_trt.rds")
