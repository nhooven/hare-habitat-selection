# PROJECT: Habitat selection
# SCRIPT: 08b - Visualization (RSS)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 01 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 03 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(hoove)
library(cowplot)
library(mefa4)

# ______________________________________________________________________________
# 2. Read in cleaned model results ----
# ______________________________________________________________________________

off.M6 <- readRDS("model_results/off_M6.rds")
off.M4 <- readRDS("model_results/off_M4.rds")

on.M6 <- readRDS("model_results/on_M6.rds")
on.M5 <- readRDS("model_results/on_M5.rds")

# ______________________________________________________________________________
# 3. Covariate summaries and functions to deal with them ----
# ______________________________________________________________________________

# covariate means, SDs, and ranges
off.summ <- readRDS("data_for_model/mean_sd_off.rds")
on.summ <- readRDS("data_for_model/mean_sd_on.rds")

# helper functions
# calculate the correct squared covariate values
correct_sq <- function (.range, .summ, .var) {
  
  # unstandardize
  range.un <- (.range * .summ$sd[.summ$name == .var]) + .summ$mean[.summ$name == .var]
    
  # square
  range.un2 <- range.un^2
  
  # standardize
  var2 <- paste0(.var, 2)
  
  range.s2 <- (range.un2 - .summ$mean[.summ$name == var2]) / .summ$sd[.summ$name == var2]
  
  return(range.s2)
  
}

# return the standardized range of a covariate
range_s <- function (.summ, .var, .length = 100) {
  
  cov.min <- (.summ$min[.summ$name == .var] - .summ$mean[.summ$name == .var]) / .summ$sd[.summ$name == .var]
  cov.max <- (.summ$max[.summ$name == .var] - .summ$mean[.summ$name == .var]) / .summ$sd[.summ$name == .var]
  
  cov.seq <- seq(cov.min, cov.max, length.out = .length)
  
  return(cov.seq)
  
}

# unstandardize values of a covariate
unstand <- function (.df, .summ, .var) {
  
  .vals <- .df[.var]
  
  .df[.var] <- (.vals * .summ$sd[.summ$name == .var]) + .summ$mean[.summ$name == .var]
  
  return(.df)
  
}

# calculate the un-log transformed values
# notably, these are also unstandardized
unlog <- function (.df, .summ, .var) {
 
  .vals <- .df[.var]
  
  .vals.unstand <- (.vals * .summ$sd[.summ$name == .var]) + .summ$mean[.summ$name == .var]
  .df[.var] <- exp(.vals.unstand)
  
  return(.df)
   
}

# ______________________________________________________________________________
# 4. RSS prediction ----

# this function will calculate RSS predictions generally for a given model and a 
# data.frame of values

# importantly, the .pred.df must include values for ALL covariates
# to begin with, we will assume x2 is always an average location (all vals = 0)

# NOTE: this cannot currently handle TRT-only responses, just the FR x TRT

# ______________________________________________________________________________

predict_rss <- function (.results,
                         .pred.df,
                         .season = "off",
                         .log = TRUE) {
  
  # extract means and credible limits
  est.mean <- .results |> dplyr::select(param, mean) |> rename(est = mean)
  est.low <- .results |> dplyr::select(param, low) |> rename(est = low)
  est.upp <- .results |> dplyr::select(param, upp) |> rename(est = upp)
  
  # for flexibility, check if a variable is in .pred.df
  # if not, add a "zero" column
  all.covars.off <- c("vo", "ch", "ch2", "cc", "cc2", "dOpen", "dDM", "ed", "twi", "twi2", "vrm", "vrm2",
                      "a.vo", "a.ch", "a.cc", "a.dOpen", "a.dDM", "a.ed",
                      "year.trt", "ret", "pil")
  
  all.covars.on <- c("stem", "ch", "ch2", "cc", "cc2", "dOpen", "dDM", "ed", "twi", "twi2", "vrm", "vrm2",
                     "a.stem", "a.ch", "a.cc", "a.dOpen", "a.dDM", "a.ed",
                     "year.trt", "ret", "pil")
  
  # which one?
  all.covars <- if (.season == "off") { all.covars.off } else { all.covars.on }
  
  # add columns if necessary
  # helper function
  add_columns <- function (.pred.df, all.covars) {
    
    # determine which covariates are not present
    which.not <- all.covars[which(all.covars %notin% colnames(.pred.df))]
    
    # bind in zero columns
    .pred.df.new <- .pred.df
    
    for (i in 1:length(which.not)) { .pred.df.new <- cbind(.pred.df.new, data.frame(0)) }
    
    # and add names
    colnames(.pred.df.new) <- c(colnames(.pred.df), which.not)
    
    return(.pred.df.new)
    
  } # f()
  
  pred.df.new <- add_columns(.pred.df, all.covars)
  
  # calculate RSS
  # for convenience, we're going to calculate the whole model structure here
  # as above, this will NOT accommodate the TRT-only interaction
  calc_rss <- function (.est) {
    
    # correct season
    if (.season == "off") {
      
      log.rss <- 
        
        # main effects
        .est$est[.est$param == "vo"] * pred.df.new$vo +
        .est$est[.est$param == "ch"] * pred.df.new$ch +
        .est$est[.est$param == "ch2"] * pred.df.new$ch2 +
        .est$est[.est$param == "cc"] * pred.df.new$cc +
        .est$est[.est$param == "cc2"] * pred.df.new$cc2 +
        .est$est[.est$param == "dOpen"] * pred.df.new$dOpen +
        .est$est[.est$param == "dDM"] * pred.df.new$dDM +
        .est$est[.est$param == "ed"] * pred.df.new$ed +
        .est$est[.est$param == "twi"] * pred.df.new$twi +
        .est$est[.est$param == "twi2"] * pred.df.new$twi2 +
        .est$est[.est$param == "vrm"] * pred.df.new$vrm +
        .est$est[.est$param == "vrm2"] * pred.df.new$vrm2 +
        
        # stand-level functional responses (base)
        .est$est[.est$param == "vo:a.vo"] * pred.df.new$vo * pred.df.new$a.vo +
        .est$est[.est$param == "ch:a.ch"] * pred.df.new$ch * pred.df.new$a.ch +
        .est$est[.est$param == "ch2:a.ch"] * pred.df.new$ch2 * pred.df.new$a.ch +
        .est$est[.est$param == "cc:a.cc"] * pred.df.new$cc * pred.df.new$a.cc +
        .est$est[.est$param == "cc2:a.cc"] * pred.df.new$cc2 * pred.df.new$a.cc +
        
        # landscape-level functional responses (base)
        .est$est[.est$param == "dOpen:a.dOpen"] * pred.df.new$dOpen * pred.df.new$a.dOpen +
        .est$est[.est$param == "dDM:a.dDM"] * pred.df.new$dDM * pred.df.new$a.dDM +
        .est$est[.est$param == "ed:a.ed"] * pred.df.new$ed * pred.df.new$a.ed +
        
        # stand level TRT
        .est$est[.est$param == "vo:year.trt:ret"] * pred.df.new$vo * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "vo:year.trt:pil"] * pred.df.new$vo * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch:year.trt:ret"] * pred.df.new$ch * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch:year.trt:pil"] * pred.df.new$ch * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch2:year.trt:ret"] * pred.df.new$ch2 * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch2:year.trt:pil"] * pred.df.new$ch2 * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc:year.trt:ret"] * pred.df.new$cc * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc:year.trt:pil"] * pred.df.new$cc * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc2:year.trt:ret"] * pred.df.new$cc2 * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc2:year.trt:pil"] * pred.df.new$cc2 * pred.df.new$year.trt * pred.df.new$pil +
        
        # landscape-level TRT
        .est$est[.est$param == "dOpen:year.trt:ret"] * pred.df.new$dOpen * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dOpen:year.trt:pil"] * pred.df.new$dOpen * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "dDM:year.trt:ret"] * pred.df.new$dDM * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dDM:year.trt:pil"] * pred.df.new$dDM * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ed:year.trt:ret"] * pred.df.new$ed * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ed:year.trt:pil"] * pred.df.new$ed * pred.df.new$year.trt * pred.df.new$pil +
        
        # stand-level FR x TRT
        .est$est[.est$param == "vo:a.vo:year.trt:ret"] * pred.df.new$vo * pred.df.new$a.vo * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "vo:a.vo:year.trt:pil"] * pred.df.new$vo * pred.df.new$a.vo * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch:a.ch:year.trt:ret"] * pred.df.new$ch * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch:a.ch:year.trt:pil"] * pred.df.new$ch * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch2:a.ch:year.trt:ret"] * pred.df.new$ch2 * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch2:a.ch:year.trt:pil"] * pred.df.new$ch2 * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc:a.cc:year.trt:ret"] * pred.df.new$cc * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc:a.cc:year.trt:pil"] * pred.df.new$cc * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc2:a.cc:year.trt:ret"] * pred.df.new$cc2 * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc2:a.cc:year.trt:pil"] * pred.df.new$cc2 * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$pil +
        
        # landscape-level FR x TRT
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:ret"] * pred.df.new$dOpen * pred.df.new$a.dOpen * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:pil"] * pred.df.new$dOpen * pred.df.new$a.dOpen * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "dDM:a.dDM:year.trt:ret"] * pred.df.new$dDM * pred.df.new$a.dDM * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dDM:a.dDM:year.trt:pil"] * pred.df.new$dDM * pred.df.new$a.dDM * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ed:a.ed:year.trt:ret"] * pred.df.new$ed * pred.df.new$a.ed * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ed:a.ed:year.trt:pil"] * pred.df.new$ed * pred.df.new$a.ed * pred.df.new$year.trt * pred.df.new$pil
        
    } else {
      
      # main effects
      .est$est[.est$param == "stem"] * pred.df.new$stem +
        .est$est[.est$param == "ch"] * pred.df.new$ch +
        .est$est[.est$param == "ch2"] * pred.df.new$ch2 +
        .est$est[.est$param == "cc"] * pred.df.new$cc +
        .est$est[.est$param == "cc2"] * pred.df.new$cc2 +
        .est$est[.est$param == "dOpen"] * pred.df.new$dOpen +
        .est$est[.est$param == "dDM"] * pred.df.new$dDM +
        .est$est[.est$param == "ed"] * pred.df.new$ed +
        .est$est[.est$param == "twi"] * pred.df.new$twi +
        .est$est[.est$param == "twi2"] * pred.df.new$twi2 +
        .est$est[.est$param == "vrm"] * pred.df.new$vrm +
        .est$est[.est$param == "vrm2"] * pred.df.new$vrm2 +
        
        # stand-level functional responses (base)
        .est$est[.est$param == "stem:a.stem"] * pred.df.new$stem * pred.df.new$a.stem +
        .est$est[.est$param == "ch:a.ch"] * pred.df.new$ch * pred.df.new$a.ch +
        .est$est[.est$param == "ch2:a.ch"] * pred.df.new$ch2 * pred.df.new$a.ch +
        .est$est[.est$param == "cc:a.cc"] * pred.df.new$cc * pred.df.new$a.cc +
        .est$est[.est$param == "cc2:a.cc"] * pred.df.new$cc2 * pred.df.new$a.cc +
        
        # landscape-level functional responses (base)
        .est$est[.est$param == "dOpen:a.dOpen"] * pred.df.new$dOpen * pred.df.new$a.dOpen +
        .est$est[.est$param == "dDM:a.dDM"] * pred.df.new$dDM * pred.df.new$a.dDM +
        .est$est[.est$param == "ed:a.ed"] * pred.df.new$ed * pred.df.new$a.ed +
        
        # stand level TRT
        .est$est[.est$param == "stem:year.trt:ret"] * pred.df.new$stem * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "stem:year.trt:pil"] * pred.df.new$stem * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch:year.trt:ret"] * pred.df.new$ch * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch:year.trt:pil"] * pred.df.new$ch * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch2:year.trt:ret"] * pred.df.new$ch2 * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch2:year.trt:pil"] * pred.df.new$ch2 * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc:year.trt:ret"] * pred.df.new$cc * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc:year.trt:pil"] * pred.df.new$cc * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc2:year.trt:ret"] * pred.df.new$cc2 * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc2:year.trt:pil"] * pred.df.new$cc2 * pred.df.new$year.trt * pred.df.new$pil +
        
        # landscape-level TRT
        .est$est[.est$param == "dOpen:year.trt:ret"] * pred.df.new$dOpen * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dOpen:year.trt:pil"] * pred.df.new$dOpen * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "dDM:year.trt:ret"] * pred.df.new$dDM * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dDM:year.trt:pil"] * pred.df.new$dDM * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ed:year.trt:ret"] * pred.df.new$ed * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ed:year.trt:pil"] * pred.df.new$ed * pred.df.new$year.trt * pred.df.new$pil +
        
        # stand-level FR x TRT
        .est$est[.est$param == "stem:a.stem:year.trt:ret"] * pred.df.new$stem * pred.df.new$a.stem * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "stem:a.stem:year.trt:pil"] * pred.df.new$stem * pred.df.new$a.stem * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch:a.ch:year.trt:ret"] * pred.df.new$ch * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch:a.ch:year.trt:pil"] * pred.df.new$ch * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ch2:a.ch:year.trt:ret"] * pred.df.new$ch2 * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ch2:a.ch:year.trt:pil"] * pred.df.new$ch2 * pred.df.new$a.ch * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc:a.cc:year.trt:ret"] * pred.df.new$cc * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc:a.cc:year.trt:pil"] * pred.df.new$cc * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "cc2:a.cc:year.trt:ret"] * pred.df.new$cc2 * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "cc2:a.cc:year.trt:pil"] * pred.df.new$cc2 * pred.df.new$a.cc * pred.df.new$year.trt * pred.df.new$pil +
        
        # landscape-level FR x TRT
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:ret"] * pred.df.new$dOpen * pred.df.new$a.dOpen * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dOpen:a.dOpen:year.trt:pil"] * pred.df.new$dOpen * pred.df.new$a.dOpen * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "dDM:a.dDM:year.trt:ret"] * pred.df.new$dDM * pred.df.new$a.dDM * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "dDM:a.dDM:year.trt:pil"] * pred.df.new$dDM * pred.df.new$a.dDM * pred.df.new$year.trt * pred.df.new$pil +
        .est$est[.est$param == "ed:a.ed:year.trt:ret"] * pred.df.new$ed * pred.df.new$a.ed * pred.df.new$year.trt * pred.df.new$ret +
        .est$est[.est$param == "ed:a.ed:year.trt:pil"] * pred.df.new$ed * pred.df.new$a.ed * pred.df.new$year.trt * pred.df.new$pil
      
    }
    
    # exponentiate if necessary
    if (.log == TRUE) { rss <- log.rss } else { rss <- exp(log.rss) } 
    
    return(rss)
    
  } # f()
  
  # add RSS columns
  pred.df.new$rss.est <- calc_rss(est.mean)
  pred.df.new$rss.low <- calc_rss(est.low)
  pred.df.new$rss.upp <- calc_rss(est.upp)
  
  return(pred.df.new)
  
} # f()

# ______________________________________________________________________________
# 5. RSS plotting ----

# this function will return a clean plot given an x, y, and (possibly) a group

# ______________________________________________________________________________

plot_rss <- function (.pred,
                      .x,
                      .group = NULL,
                      .log = TRUE) {
  
  # change x name
  colnames(.pred)[colnames(.pred) == .x] <- "x"
  
  # change group name if necessary, if not, include a placeholder
  if (is.null(.group) == FALSE) { 
    
    colnames(.pred)[colnames(.pred) == .group] <- "group" 
    
    .pred$group <- factor(.pred$group)
    
  } else {
      
    .pred$group <- "1"
    
  }
  
  # informative x-axis title
  x.titles <- data.frame(var = c("vo", "stem", "ch", "cc", "dOpen", "dDM", "ed", "twi", "vrm"),
                         title = c("Visual obstruction (%)",
                                   "Conifer stem density (stems/ha)",
                                   "Canopy height (m)",
                                   "Canopy cover (%)",
                                   "Distance to open (m)",
                                   "Distance to dense mature (m)",
                                   "Edge density (m/ha)",
                                   "Wetness",
                                   "Ruggedness"))
  
  # plot
  ggplot(data = .pred) +
    
    theme_hoove() +
    
    # indifference line
    geom_hline(yintercept = ifelse(.log == T, 0, 1)) +
    
    # credible interval
    geom_ribbon(aes(x = x,
                    ymin = rss.low,
                    ymax = rss.upp,
                    group = group,
                    fill = group),
                alpha = 0.25) +
    
    # mean
    geom_line(aes(x = x,
                  y = rss.est,
                  color = group),
              linewidth = 1.2) +
    
    # axis titles
    ylab(ifelse(.log == T, 
                "log(Relative selection strength)", 
                "Relative selection strength")) +
    
    xlab(x.titles$title[x.titles$var == .x])
  
}

# ______________________________________________________________________________
# 6. Test plots ----
# ______________________________________________________________________________

# linear
off.M6[[1]] |> 
  predict_rss(data.frame(vo = range_s(off.summ, "vo"))) |> 
  unstand(off.summ, "vo") |> 
  plot_rss("vo")

# squared
off.M6[[1]] |> 
  predict_rss(data.frame(twi = range_s(off.summ, "twi"),
                         twi2 = correct_sq(range_s(off.summ, "twi"), off.summ, "twi"))) |> 
  unstand(off.summ, "twi") |> 
  plot_rss("twi")

# distance
off.M6[[1]] |> 
  predict_rss(data.frame(dDM = range_s(off.summ, "dDM"))) |> 
  unstand(off.summ, "dDM") |> 
  plot_rss("dDM")

# functional responses
fr.pred.df <- data.frame(cc = rep(range_s(off.summ, "cc"), 3),
                         cc2 = rep(correct_sq(range_s(off.summ, "cc"), off.summ, "cc"), 3),
                         a.cc = c(rep(-3, 100), rep(0, 100), rep(3, 100),
                                  rep(-3, 100), rep(0, 100), rep(3, 100),
                                  rep(-3, 100), rep(0, 100), rep(3, 100)),
                         year.trt = 1,
                         ret = c(rep(0, 300), rep(1, 300), rep(0, 300)),
                         pil = c(rep(0, 300), rep(0, 300), rep(1, 300)),
                         z = c(rep("low.unt", 100), rep("mid.unt", 100), rep("hi.unt", 100),
                               rep("low.ret", 100), rep("mid.ret", 100), rep("hi.ret", 100),
                               rep("low.pil", 100), rep("mid.pil", 100), rep("hi.pil", 100)))

off.M6[[1]] |> 
  predict_rss(fr.pred.df) |> 
  unstand(off.summ, "cc") |> 
  plot_rss("cc", .group = "z")

# ______________________________________________________________________________
# 7. Interactions plots ----
# ______________________________________________________________________________

plot_rss_int <- function (.pred,
                          .x,
                          .group,
                          .log = TRUE) {
  
  # change x name
  colnames(.pred)[colnames(.pred) == .x] <- "x"
  
  # add treatment variable
  .pred$trt <- case_when(
    
    .pred$year.trt == 0 | (.pred$ret == 0 & .pred$pil == 0) ~ "unthinned",
    .pred$year.trt == 1 & .pred$ret == 1 & .pred$pil == 0 ~ "retention",
    .pred$year.trt == 1 & .pred$ret == 0 & .pred$pil == 1 ~ "piling"
    
  )
  
  # group name
  colnames(.pred)[colnames(.pred) == .group] <- "group" 
  
  .pred$group <- factor(.pred$group)
  
  # informative x-axis title
  x.titles <- data.frame(var = c("vo", "stem", "ch", "cc", "dOpen", "dDM", "ed", "twi", "vrm"),
                         title = c("Visual obstruction (%)",
                                   "Conifer stem density (stems/ha)",
                                   "Canopy height (m)",
                                   "Canopy cover (%)",
                                   "Distance to open (m)",
                                   "Distance to dense mature (m)",
                                   "Edge density (m/ha)",
                                   "Wetness",
                                   "Ruggedness"))
  
  # plot
  ggplot(data = .pred) +
    
    theme_bw() +
    
    facet_wrap(~ trt) +
    
    # indifference line
    geom_hline(yintercept = ifelse(.log == T, 0, 1)) +
    
    # credible interval
    geom_ribbon(aes(x = x,
                    ymin = rss.low,
                    ymax = rss.upp,
                    group = group,
                    fill = group),
                alpha = 0.25) +
    
    # mean
    geom_line(aes(x = x,
                  y = rss.est,
                  color = group),
              linewidth = 1.2) +
    
    # axis titles
    ylab(ifelse(.log == T, 
                "log(Relative selection strength)", 
                "Relative selection strength")) +
    
    xlab(x.titles$title[x.titles$var == .x]) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          strip.background = element_rect(color = NA)) +
    
    scale_color_viridis_d() +
    scale_fill_viridis_d()
  
} # f()

# landscape
fr.pred.df <- data.frame(ch = rep(range_s(off.summ, "ch"), 3),
                         ch2 = correct_sq(rep(range_s(off.summ, "ch"), 3), off.summ, "ch"),
                         a.ch = c(rep(-3, 100), rep(0, 100), rep(3, 100),
                                     rep(-3, 100), rep(0, 100), rep(3, 100),
                                     rep(-3, 100), rep(0, 100), rep(3, 100)),
                         year.trt = 1,
                         ret = c(rep(0, 300), rep(1, 300), rep(0, 300)),
                         pil = c(rep(0, 300), rep(0, 300), rep(1, 300)),
                         z = c(rep("low.unt", 100), rep("mid.unt", 100), rep("hi.unt", 100),
                               rep("low.ret", 100), rep("mid.ret", 100), rep("hi.ret", 100),
                               rep("low.pil", 100), rep("mid.pil", 100), rep("hi.pil", 100)))

off.M6[[1]] |> 
  predict_rss(fr.pred.df) |> 
  unstand(off.summ, "ch") |> 
  plot_rss_int("ch", .group = "a.ch")

off.M6[[1]] |> 
  predict_rss(fr.pred.df) |> 
  unstand(off.summ, "vo") |> 
  plot_rss_int("vo", "a.vo")

fr.pred.df.1 <- data.frame(ed = rep(range_s(off.summ, "ed"), 3),
                         a.ed = c(rep(-3, 100), rep(0, 100), rep(3, 100),
                                  rep(-3, 100), rep(0, 100), rep(3, 100),
                                  rep(-3, 100), rep(0, 100), rep(3, 100)),
                         year.trt = 1,
                         ret = c(rep(0, 300), rep(1, 300), rep(0, 300)),
                         pil = c(rep(0, 300), rep(0, 300), rep(1, 300)),
                         z = c(rep("low.unt", 100), rep("mid.unt", 100), rep("hi.unt", 100),
                               rep("low.ret", 100), rep("mid.ret", 100), rep("hi.ret", 100),
                               rep("low.pil", 100), rep("mid.pil", 100), rep("hi.pil", 100)))

off.M6[[1]] |> 
  predict_rss(fr.pred.df.1) |> 
  unstand(off.summ, "ed") |> 
  plot_rss_int("ed", "a.ed")
