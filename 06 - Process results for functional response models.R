# PROJECT: Habitat selection
# SCRIPT: 06 - Process results for functional response models
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 04 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 05 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# results
M.off <- readRDS("model_results/M_off.rds")
M.on <- readRDS("model_results/M_on.rds")

# raw data
data.off <- readRDS("data_cleaned/data_off.rds")
data.on <- readRDS("data_cleaned/data_on.rds")

# ______________________________________________________________________________
# 3. Join slopes and (deviation) SEs together ----
# ______________________________________________________________________________

# function
join_slopes <- function (.slopes, .sds, .season) {
  
  if (.season == "off") {
    
    # pivot
  slopes.joined <- .slopes |> 
    
    # pivot slopes
    pivot_longer(cols = vo:ed) |>
    
    rename(param = name,
           beta = value) |>
    
    # join
    left_join(
      
      .sds |> 
        
        pivot_longer(cols = vo:ed) |>
        
        rename(param = name,
               sd = value)
      
    )
    
  } else {
    
    # pivot
    slopes.joined <- .slopes |> 
      
      # pivot slopes
      pivot_longer(cols = stem:ed) |>
      
      rename(param = name,
             beta = value) |>
      
      # join
      left_join(
        
        .sds |> 
          
          pivot_longer(cols = stem:ed) |>
          
          rename(param = name,
                 sd = value)
        
      )
    
  }
  
  return(slopes.joined)
  
}

# use
off.slopes <- join_slopes(M.off[[3]], M.off[[4]], "off")
on.slopes <- join_slopes(M.on[[3]], M.on[[4]], "on")

# ______________________________________________________________________________
# 4. FR predictors ----
# ______________________________________________________________________________

# function
calc_fr <- function (.data) {
  
  # calculate availability
  data.avail <- .data |>
    
    filter(case == 0) |>
    
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
    
    group_by(track_season_post) |>
    
    summarize(a.vo = mean(vo),
              a.stem = mean(stem),
              a.ch = mean(ch),
              a.cc = mean(cc),
              a.dOpen = mean(dOpen),
              a.dDM = mean(dDM),
              a.ed = mean(ed)) |>
    
    rename(TSPID = track_season_post)
  
  # factors for modeling
  data.fact <- .data |>
    
    group_by(track_season_post) |>
    
    slice(1) |>
    
    ungroup() |>
    
    # three-level factor for treatment
    mutate(TRT = case_when(
      
      year == "PRE" ~ "UNTHIN",
      year %in% c("POST1", "POST2") & c.trt == "CTRL" ~ "UNTHIN",
      year %in% c("POST1", "POST2") & c.trt == "RET" ~ "RET",
      year %in% c("POST1", "POST2") & c.trt == "PIL" ~ "PIL"
      
    )
    
    ) |>
    
    rename(TSPID = track_season_post) |>
    
    dplyr::select(TSPID, sex, TRT) |>
    
    left_join(data.avail)
    
    return(data.fact)
  
}

# use
off.fr <- left_join(off.slopes, calc_fr(data.off))
on.fr <- left_join(on.slopes, calc_fr(data.on))

# ______________________________________________________________________________
# 5. Save to file ----
# ______________________________________________________________________________

saveRDS(off.fr, "data_for_model/off_fr.rds")
saveRDS(on.fr, "data_for_model/on_fr.rds")

