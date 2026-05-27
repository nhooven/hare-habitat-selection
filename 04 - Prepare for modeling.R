# PROJECT: Habitat selection
# SCRIPT: 04 - Prepare for modeling
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 
# LAST MODIFIED: 27 May 2026
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
# 4. Subsample available points ----

# right now we're using ~100 available per used
# let's go down to 20

# ______________________________________________________________________________

# split by track_season_post
data.off.split <- split(data.off, ~ track_season_post)
data.on.split <- split(data.on, ~ track_season_post)

# function to subsample
subsample_avail <- function (x) {
  
  # split 
  x.1 <- x %>% filter(case == 1)
  
  # how many avail?
  total.avail <- nrow(x) - nrow(x.1)
  
  x.0 <- x %>% filter(case == 0) %>%
    
    # subsample
    slice(sample(1:total.avail, size = round(total.avail / 5, digits = 0)))
  
  # bind together
  x.all <- rbind(x.1, x.0)
  
  return(x.all)
  
}

# apply function
data.off <- do.call(rbind, lapply(data.off.split, subsample_avail))
data.on <- do.call(rbind, lapply(data.on.split, subsample_avail))

# remove split lists
rm(data.off.split)
rm(data.on.split)

# ______________________________________________________________________________
# 5. Prepare for modeling ----
# ______________________________________________________________________________
# 5a. Correct values for each treatment, drop NAs ----
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
      
      # SHDI
      shdi = case_when(
        
        year == "PRE" ~ shdi.100.pre,
        year %in% c("POST1", "POST2") ~ shdi.100.post
        
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
    
    # indicators for year (just pre-post treatment for now)
    # and each treatment
    mutate(
      
      year.trt = case_when(
        
        year == "PRE" ~ 0,
        year %in% c("POST1", "POST2") ~ 1
        
      ),
      
      ret = ifelse(trt == "RET", 1, 0),
      pil = ifelse(trt == "PIL", 1, 0)
      
    ) |>
    
    # drop variables
    dplyr::select(
      
      -c(year,
         trt,
         season,
         ch.pre,
         ch.post,
         cc.pre,
         cc.post,
         shdi.100.pre,
         shdi.100.post,
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
        shdi,
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
      sex,
      MRID,
      year.trt,
      ret,
      pil,
      case,
      w,
      
      # predictors
      akde,
      
      # landscape covariates
      dOpen,
      dDM,
      ch,
      cc,
      shdi,
      vo,
      stem,
      twi,
      vrm
      
    )
  
  return(x.1)
  
}

# use
data.off.1 <- prep_1(data.off)
data.on.1 <- prep_1(data.on)

# ______________________________________________________________________________
# 5b. Mean availability ----
# ______________________________________________________________________________

calc_avail <- function (x) {
  
  x.1 <- x |>
    
    # keep identifier and each landscape covariate
    dplyr::select(track_season_post,
                  dOpen:vrm) |>
    
    # pivot longer
    pivot_longer(dOpen:vrm) |>
    
    # group and summarize
    group_by(track_season_post, name) |>
    summarize(mean = mean(value)) |>
    
    # remove groups
    ungroup() |>
    
    # pivot wider 
    pivot_wider(names_from = name,
                names_prefix = "a.",
                values_from = mean)
  
  return(x.1)
  
}

# use
data.off.2 <- data.off.1 |> left_join(calc_avail(data.off.1))
data.on.2 <- data.on.1 |> left_join(calc_avail(data.on.1))

# ______________________________________________________________________________
# 5c. Transformations ----

# function
transform_covs <- function (x) {
  
  x.1 <- x |>
    
    # covariate transformations
    # log
    mutate(
      
      dOpen = log(dOpen + 1),
      dDM = log(dDM + 1)
      
      ) |>
    
    # squared
    mutate(
      
      across(
        
        ch:vrm,
        function (x) x * x,
        .names = "{.col}2"
        
      )
      
    ) |>
    
    # availability transformations
    # these will all be log for now
    mutate(
      
      across(
        
        a.cc:a.vrm,
        log
        
      )
      
    )
  
  return(x.1)
  
}

# ______________________________________________________________________________

data.off.3 <- transform_covs(data.off.2)
data.on.3 <- transform_covs(data.on.2)

# ______________________________________________________________________________
# 5c. Save means/SDs ----

# function
mean_sds <- function (x) {
  
  x.1 <- x |>
    
    pivot_longer(cols = c(dOpen:vrm2)) |>
    
    group_by(name) |>
    
    summarize(
      
      mean = mean(value),
      sd = sd(value)
      
    )
  
  return(x.1)
  
}

# ______________________________________________________________________________

# use
mean.sd.off <- mean_sds(data.off.3)
mean.sd.on <- mean_sds(data.on.3)

# ______________________________________________________________________________
# 5d. Standardize ----

# function
standardize_across <- function (x) {
  
  # standardize function
  standardize <- function (x) {
    
    x.1 <- (x - mean(x)) / sd(x)
    
    return(x.1)
    
  }
  
  x.1 <- x |>
    
    mutate(
      
      across(dOpen:vrm2,
             standardize)
      
    )
  
  return(x.1)
  
}

# ______________________________________________________________________________

data.off.4 <- standardize_across(data.off.3)
data.on.4 <- standardize_across(data.on.3)

# ______________________________________________________________________________
# 6. Correlation ----
# ______________________________________________________________________________

cor(data.off.4 %>% dplyr::select(c(dOM:vrm2)), method = "spearman")
cor(data.on.4 %>% dplyr::select(c(dOM:vrm2)), method = "spearman")

# ______________________________________________________________________________
# 7. Save to files ----
# ______________________________________________________________________________

saveRDS(data.off.4, "data_for_model/off_data.rds")
saveRDS(data.on.4, "data_for_model/on_data.rds")

saveRDS(mean.sd.off, "data_for_model/mean_sd_off.rds")
saveRDS(mean.sd.on, "data_for_model/mean_sd_on.rds")
