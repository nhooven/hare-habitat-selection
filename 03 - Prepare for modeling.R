# PROJECT: Habitat selection
# SCRIPT: 03 - Prepare for modeling
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 21 Apr 2026
# COMPLETED: 01 Jun 2026
# LAST MODIFIED: 01 Jun 2026
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
      
      # SHDI
      #shdi = case_when(
        
      #  year == "PRE" ~ shdi.pre,
      #  year %in% c("POST1", "POST2") ~ shdi.post
        
      #),
      
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
        vrm,
        north,
        east
        
      )
      
    ) |>
    
    # case weights
    mutate(w = ifelse(case == 0, 5000, 1)) |>
    
    # cluster
    mutate(cluster = substr(site, 1, 1)) |>
    
    # keep correct variables
    dplyr::select(
      
      track_season_post,   # for calculating availability
      cluster,
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
      vo,
      stem,
      ch,
      cc,
      twi,
      vrm,
      north,
      east,
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
# 4. Mean availability ----
# ______________________________________________________________________________

calc_avail <- function (x) {
  
  x.1 <- x |>
    
    # keep identifier and each landscape covariate
    dplyr::select(track_season_post,
                  vo:ed) |>
    
    # pivot longer
    pivot_longer(vo:ed) |>
    
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
# 5. Extract attributes to bind in later ----
# ______________________________________________________________________________

data.off.attr <- data.off.2 |> dplyr::select(track_season_post:akde)
data.on.attr <- data.on.2 |> dplyr::select(track_season_post:akde)

data.off.cov <- data.off.2 |> dplyr::select(vo:a.vrm)
data.on.cov <- data.on.2 |> dplyr::select(vo:a.vrm)

# ______________________________________________________________________________
# 6. Transformations ----

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
      
      ch2 = ch^2,
      cc2 = cc^2,
      twi2 = twi^2,
      vrm2 = vrm^2,
      north2 = north^2,
      east2 = east^2
        
      ) 
  
  #|>
    
    # availability transformations
    # these will all be log for now
    #mutate(
      
      #across(
        
        #a.cc:a.vrm,
        #log,
        #.names = "l.{.col}"
        
      #)
      
    #)
  
  return(x.1)
  
}

# ______________________________________________________________________________

data.off.3 <- transform_covs(data.off.cov)
data.on.3 <- transform_covs(data.on.cov)

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

# ______________________________________________________________________________

# use
mean.sd.off <- mean_sds(data.off.3)
mean.sd.on <- mean_sds(data.on.3)

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

data.off.4 <- standardize_across(data.off.3)
data.on.4 <- standardize_across(data.on.3)

# ______________________________________________________________________________
# 8. Bind back in ----
# ______________________________________________________________________________

data.off.5 <- cbind(data.off.attr, data.off.4)
data.on.5 <- cbind(data.on.attr, data.on.4)

# ______________________________________________________________________________
# 9. Save to files ----
# ______________________________________________________________________________

saveRDS(data.off.5, "data_for_model/off_data.rds")
saveRDS(data.on.5, "data_for_model/on_data.rds")

saveRDS(mean.sd.off, "data_for_model/mean_sd_off.rds")
saveRDS(mean.sd.on, "data_for_model/mean_sd_on.rds")
