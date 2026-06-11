# PROJECT: Habitat selection
# SCRIPT: 09 - Assess predictive ability
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 09 Jun 2026
# R VERSION: 4.5.2

# here we want to look at overall predictive ability
# as well as by-individual

# we can compare the base HSF and the implied FR HSF for this

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(mgcv)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

# HSF results
M.off <- readRDS("model_results/M_off.rds")
M.on <- readRDS("model_results/M_on.rds")

# fr models
# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.dOpen <- readRDS("model_results/fr_models/off_dOpen.rds")
off.dDM <- readRDS("model_results/fr_models/off_dDM.rds")
off.ed <- readRDS("model_results/fr_models/off_ed.rds")

# on
on.stem <- readRDS("model_results/fr_models/on_stem.rds")
on.ch <- readRDS("model_results/fr_models/on_ch.rds")
on.cc2 <- readRDS("model_results/fr_models/on_cc2.rds")
on.dOpen <- readRDS("model_results/fr_models/on_dOpen.rds")
on.dDM <- readRDS("model_results/fr_models/on_dDM.rds")

# data
off.fr.data <- readRDS("data_for_model/off_fr.rds")
on.fr.data <- readRDS("data_for_model/on_fr.rds")

# HS data
off.hs.data <- readRDS("data_for_model/off_data.rds")
on.hs.data <- readRDS("data_for_model/on_data.rds")

# ______________________________________________________________________________
# 3. Boyce index ----
# ______________________________________________________________________________
# 3a. Generic function ---
# ______________________________________________________________________________

boyce <- function (.hsf,
                   .which = "pop",      # population vs individual level
                   .pred = "base") {
  
  if (.hsf == "M.off") {
    
    # use/availability data
    hs.data <- off.hs.data |>
      
      # add TRT
      mutate(TRT = case_when(
        
        year == "PRE" ~ "UNTHIN",
        year %in% c("POST1", "POST2") & c.trt == "CTRL" ~ "UNTHIN",
        year %in% c("POST1", "POST2") & c.trt == "RET" ~ "RET",
        year %in% c("POST1", "POST2") & c.trt == "PIL" ~ "PIL"
        
      )
      
    ) |>
      
      dplyr::select(track_season_post, case, TRT, 
                    vo, ch, cc, cc2, twi, twi2, vrm, vrm2, dOpen, dDM, ed) |>
      
      rename(TSPID = track_season_post)
    
    # base coefficients
    coef.base <- M.off[[1]] |> dplyr::select(param, mean)
    
    # FR coefficients
    coef.fr <- off.fr.data |> dplyr::select(TSPID, TRT, a.vo:a.ed) |>
      
      group_by(TSPID) |> slice(1) |> ungroup() |>
      
      mutate(
        
        beta.vo = predict(off.vo, data.frame(TRT = TRT, avail = a.vo)),
        beta.dOpen = predict(off.dOpen, data.frame(TRT = TRT, avail = a.dOpen)),
        beta.dDM = predict(off.dDM, data.frame(TRT = TRT, avail = a.dDM)),
        beta.ed = predict(off.ed, data.frame(TRT = TRT, avail = a.ed))
        
      ) |>
      
      dplyr::select(TSPID, beta.vo:beta.ed)
    
    # add implied coefs to all use/avail
    hs.data.fr <- hs.data |>
      
      dplyr::select(TSPID) |>
      
      left_join(coef.fr)
    
    # calculate w(x)
    hs.data$pred.base <- exp(
      
      hs.data$vo * coef.base$mean[coef.base$param == "vo"] +
      hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
      hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
      hs.data$cc2 * coef.base$mean[coef.base$param == "cc2"] +
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
      hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
      hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
      hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
      hs.data$dOpen * coef.base$mean[coef.base$param == "dOpen"] +
      hs.data$dDM * coef.base$mean[coef.base$param == "dDM"] +
      hs.data$ed * coef.base$mean[coef.base$param == "ed"]
      
    )
    
    hs.data$pred.fr <- exp(
      
      hs.data$vo * hs.data.fr$beta.vo +
      hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
      hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
      hs.data$cc2 * coef.base$mean[coef.base$param == "cc2"] +
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
      hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
      hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
      hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
      hs.data$dOpen * hs.data.fr$beta.dOpen +
      hs.data$dDM * hs.data.fr$beta.dDM +
      hs.data$ed * hs.data.fr$beta.ed
      
    )
    
  } # off
  
  if (.hsf == "M.on") {
    
    # use/availability data
    hs.data <- on.hs.data |>
      
      # add TRT
      mutate(TRT = case_when(
        
        year == "PRE" ~ "UNTHIN",
        year %in% c("POST1", "POST2") & c.trt == "CTRL" ~ "UNTHIN",
        year %in% c("POST1", "POST2") & c.trt == "RET" ~ "RET",
        year %in% c("POST1", "POST2") & c.trt == "PIL" ~ "PIL"
        
      )
      
      ) |>
      
      dplyr::select(track_season_post, case, TRT, 
                    stem, ch, cc, cc2, twi, twi2, vrm, vrm2, dOpen, dDM, ed) |>
      
      rename(TSPID = track_season_post)
    
    # base coefficients
    coef.base <- M.on[[1]] |> dplyr::select(param, mean)
    
    # FR coefficients
    coef.fr <- on.fr.data |> dplyr::select(TSPID, TRT, a.stem:a.ed) |>
      
      group_by(TSPID) |> slice(1) |> ungroup() |>
      
      mutate(
        
        beta.stem = predict(on.stem, data.frame(TRT = TRT, avail = a.stem)),
        beta.ch = predict(on.ch, data.frame(TRT = TRT, avail = a.ch)),
        beta.cc2 = predict(on.cc2, data.frame(TRT = TRT, avail = a.cc)),
        beta.dOpen = predict(on.dOpen, data.frame(TRT = TRT, avail = a.dOpen)),
        beta.dDM = predict(on.dDM, data.frame(TRT = TRT, avail = a.dDM))
        
      ) |>
      
      dplyr::select(TSPID, beta.stem:beta.dDM)
    
    # add implied coefs to all use/avail
    hs.data.fr <- hs.data |>
      
      dplyr::select(TSPID) |>
      
      left_join(coef.fr)
    
    # calculate w(x)
    hs.data$pred.base <- exp(
      
      hs.data$stem * coef.base$mean[coef.base$param == "stem"] +
        hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
        hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
        hs.data$cc2 * coef.base$mean[coef.base$param == "cc2"] +
        hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$dOpen * coef.base$mean[coef.base$param == "dOpen"] +
        hs.data$dDM * coef.base$mean[coef.base$param == "dDM"] +
        hs.data$ed * coef.base$mean[coef.base$param == "ed"]
      
    )
    
    hs.data$pred.fr <- exp(
      
      hs.data$stem * hs.data.fr$beta.stem +
        hs.data$ch * hs.data.fr$beta.ch +
        hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
        hs.data$cc2 * hs.data.fr$beta.cc2 +
        hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$dOpen * hs.data.fr$beta.dOpen +
        hs.data$dDM * hs.data.fr$beta.dDM +
        hs.data$ed * coef.base$mean[coef.base$param == "ed"]
      
    )
    
  } # on
      
    # data.frames for Boyce index
    # helper function
    bin_hsf <- function (x) {
      
      # create quantiles
      quant <- quantile(x, seq(0, 1, 0.1), na.rm = T)
      
      # fill bin variable
      bin <- rep(NA, length(x))
      
      for (j in 1:10) {
        
        bin[x >= quant[j] & x < quant[j + 1]] <- j
        
      }
      
      # switch NAs to 10
      bin[is.na(bin)] <- 10
      
      return(bin)
      
    }
    
    # population-level
    if (.which == "pop") {
      
      # base
      if (.pred == "base") {
        
        boyce.data.base <- hs.data |>
          
          dplyr::select(case, pred.base) |>
          
          mutate(bins = bin_hsf(pred.base)) |>
          
          # group by bin and compute U/A ratio
          group_by(case, bins) |>
          
          summarize(total = n()) |>
          
          pivot_wider(names_from = case, values_from = total) |>
          
          ungroup() |>
          
          mutate(u.a.ratio = `1` / `0`) |>
          
          dplyr::select(bins, u.a.ratio)
        
        return(boyce.data.base)
        
      }
      
      if (.pred == "fr") {
        
        boyce.data.fr <- hs.data |>
          
          dplyr::select(case, pred.fr) |>
          
          mutate(bins = bin_hsf(pred.fr)) |>
          
          # group by bin and compute U/A ratio
          group_by(case, bins) |>
          
          summarize(total = n()) |>
          
          pivot_wider(names_from = case, values_from = total) |>
          
          ungroup() |>
          
          mutate(u.a.ratio = `1` / `0`) |>
          
          dplyr::select(bins, u.a.ratio)
        
        return(boyce.data.fr)
        
      }
      
    } # pop
    
    if (.which == "indiv") {
      
      # split hs.data
      hs.data.split <- split(hs.data, ~TSPID)
      
      # base
      if (.pred == "base") {
        
        # function to apply
        indiv_boyce_data_base <- function (x) {
          
          suppressMessages(
          
          x.1 <- x |>
            
            dplyr::select(case, pred.base) |>
            
            mutate(bins = bin_hsf(pred.base)) |>
            
            # group by bin and compute U/A ratio
            group_by(case, bins) |>
            
            summarize(total = n()) |>
            
            pivot_wider(names_from = case, values_from = total) |>
            
            ungroup() |>
            
            mutate(u.a.ratio = `1` / `0`) |>
            
            dplyr::select(bins, u.a.ratio) |>
            
            drop_na(u.a.ratio) |>
            
            mutate(TSPID = x$TSPID[1])
          
          )
          
        }
        
        indiv.boyce.data.base.list <- lapply(hs.data.split, indiv_boyce_data_base)
        
        # bind together for plotting
        indiv.boyce.data.base <- do.call(rbind, indiv.boyce.data.base.list)
        
        return(indiv.boyce.data.base)
      
    } # base
      
    if (.pred == "fr") {
      
      # function to apply
      indiv_boyce_data_fr <- function (x) {
        
        suppressMessages(
          
          x.1 <- x |>
            
            dplyr::select(case, pred.fr) |>
            
            mutate(bins = bin_hsf(pred.fr)) |>
            
            # group by bin and compute U/A ratio
            group_by(case, bins) |>
            
            summarize(total = n()) |>
            
            pivot_wider(names_from = case, values_from = total) |>
            
            ungroup() |>
            
            mutate(u.a.ratio = `1` / `0`) |>
            
            dplyr::select(bins, u.a.ratio) |>
            
            drop_na(u.a.ratio) |>
            
            mutate(TSPID = x$TSPID[1])
          
        )
        
      }
      
      indiv.boyce.data.fr.list <- lapply(hs.data.split, indiv_boyce_data_fr)
      
      # bind together for plotting
      indiv.boyce.data.fr <- do.call(rbind, indiv.boyce.data.fr.list)
      
      return(indiv.boyce.data.fr)
      
    } # fr
    
  } # indiv
  
} # f()

# ______________________________________________________________________________
# 3b. Use ---
# ______________________________________________________________________________

off.pop <- rbind(
  
  cbind(boyce("M.off", "pop", "base"),
        predictions = "base"),
  
  cbind(boyce("M.off", "pop", "fr"),
        predictions = "fr")
  
)

on.pop <- rbind(
  
  cbind(boyce("M.on", "pop", "base"),
        predictions = "base"),
  
  cbind(boyce("M.on", "pop", "fr"),
        predictions = "fr")
  
)

# individual-level
off.indiv.base <- boyce("M.off", "indiv", "base")
off.indiv.fr <- boyce("M.off", "indiv", "fr")

on.indiv.base <- boyce("M.on", "indiv", "base")
on.indiv.fr <- boyce("M.on", "indiv", "fr")

# ______________________________________________________________________________
# 4. Correlations ---
# ______________________________________________________________________________
# 4a. Population-level ----
# ______________________________________________________________________________

# off
off.pop.base <- off.pop |> filter(predictions == "base")
off.pop.fr <- off.pop |> filter(predictions == "fr")

cor.test(off.pop.base$bins, off.pop.base$u.a.ratio)$estimate
cor.test(off.pop.fr$bins, off.pop.fr$u.a.ratio)$estimate

# on
on.pop.base <- on.pop |> filter(predictions == "base")
on.pop.fr <- on.pop |> filter(predictions == "fr")

cor.test(on.pop.base$bins, on.pop.base$u.a.ratio)$estimate
cor.test(on.pop.fr$bins, on.pop.fr$u.a.ratio)$estimate

# ______________________________________________________________________________
# 4b. Individual-level ----
# ______________________________________________________________________________

indiv_corr <- function (x) {
  
  if (nrow(x) > 1) {
  
  corr <- cor.test(x$bins, x$u.a.ratio)$estimate
  
  return(corr)
  
  }
  
}

# off
off.indiv.base.corr <- lapply(split(off.indiv.base, ~TSPID), indiv_corr) |>
  
  do.call(rbind, args = _) |>
  
  as.data.frame()

off.indiv.fr.corr <- lapply(split(off.indiv.fr, ~TSPID), indiv_corr)  |>
  
  do.call(rbind, args = _) |>
  
  as.data.frame()

# on
on.indiv.base.corr <- lapply(split(on.indiv.base, ~TSPID), indiv_corr) |>
  
  do.call(rbind, args = _) |>
  
  as.data.frame()

on.indiv.fr.corr <- lapply(split(on.indiv.fr, ~TSPID), indiv_corr)  |>
  
  do.call(rbind, args = _) |>
  
  as.data.frame()

# ______________________________________________________________________________
# 5. Plots ---
# ______________________________________________________________________________
# 5a. U/A vs bin ----
# ______________________________________________________________________________

bin_plot <- function (.pred,
                      .which = "pop") {
  
  if (.which == "pop") {
  
  # labels
  .pred$predictions <- factor(.pred$predictions,
                              labels = c("base", "functional response"))
  
  ggplot(data = .pred) +
    
    theme_classic() +
    
    geom_line(aes(x = bins,
                  y = u.a.ratio,
                  color = predictions),
              linewidth = 1.1) +
    
    geom_point(aes(x = bins,
                   y = u.a.ratio,
                   color = predictions),
               size = 2.25,
               shape = 21,
               fill = "white") +
    
    scale_x_continuous(breaks = 1:10) +
    
    theme(legend.title = element_blank(),
          legend.position = c(0.7, 0.2)) +
    
    xlab("HSF bin") +
    ylab("Use/availability ratio")
  
  # individual plots
  } else {
    
    
    
    
  }
  
}

# ______________________________________________________________________________
# 5b. Distribution of correlations ----
# ______________________________________________________________________________

corr_dist_plot <- function (x) {
  
  ggplot(data = x) +
    
    theme_classic() +
    
    geom_histogram(aes(x = cor),
                   fill = "aquamarine3",
                   color = "white",
                   bins = 25) +
    
    # median and mean
    geom_vline(xintercept = median(x$cor),
               linewidth = 0.9) +
    geom_vline(xintercept = mean(x$cor), 
               linetype = "dashed",
               linewidth = 0.9) +
    
    coord_cartesian(xlim = c(-1.0, 1.0),
                    ylim = c(0, 10)) +
    
    scale_y_continuous(breaks = seq(0, 10, 2)) +
    
    xlab("Spearman's correlation") +
    ylab("Count of individuals")
  
}

# use
corr_dist_plot(off.indiv.base.corr)
corr_dist_plot(off.indiv.fr.corr)

corr_dist_plot(on.indiv.base.corr)
corr_dist_plot(on.indiv.fr.corr)

# ______________________________________________________________________________
# 6. Correlation distributions by TRT and sex ----
# ______________________________________________________________________________
# 6a. Bind in indices ----
# ______________________________________________________________________________

# off
off.indiv.df <- bind_rows(
  
  data.frame(TSPID = rownames(off.indiv.base.corr),
             cor = off.indiv.base.corr$cor,
             season = "off",
             model = "base"),
  
  data.frame(TSPID = rownames(off.indiv.fr.corr),
             cor = off.indiv.fr.corr$cor,
             season = "off",
             model = "fr")
) |>
  
  # join in
  left_join(
    
    off.fr.data |> group_by(TSPID) |> slice(1) |> dplyr::select(TSPID, sex, TRT)
    
  )

# on
on.indiv.df <- bind_rows(
  
  data.frame(TSPID = rownames(on.indiv.base.corr),
             cor = on.indiv.base.corr$cor,
             season = "on",
             model = "base"),
  
  data.frame(TSPID = rownames(on.indiv.fr.corr),
             cor = on.indiv.fr.corr$cor,
             season = "on",
             model = "fr")
) |>
  
  # join in
  left_join(
    
    on.fr.data |> group_by(TSPID) |> slice(1) |> dplyr::select(TSPID, sex, TRT)
    
  )

# ______________________________________________________________________________
# 6b. Plot ----
# ______________________________________________________________________________

plot_corr_by_fact <- function (.df,
                               .season) {
  
  col.value <- ifelse(.season == "off", "green4", "dodgerblue3")
  
  # factor levels
  .df <- .df |>
    
    mutate(TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")),
           model = factor(model, labels = c("base", "functional response")))
  
  ggplot(.df) +
    
    theme_bw() +
    
    facet_grid(model ~ TRT) +
    
    geom_density(aes(x = cor,
                     linetype = sex),
                 fill = col.value,
                 color = col.value,
                 linewidth = 0.7,
                 alpha = 0.15) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          strip.background = element_rect(color = NA),
          strip.text = element_text(hjust = 0),
          axis.title.y = element_blank(),
          legend.position = c(0.15, 0.85)) +
    
    coord_cartesian(xlim = c(-0.8, 1),
                    ylim = c(0, 3)) +
    
    xlab("Spearman's correlation")
  
}

plot_corr_by_fact(off.indiv.df, "off")
plot_corr_by_fact(on.indiv.df, "on")

# sex only
plot_corr_by_sex <- function (.df,
                              .season) {
  
  col.value <- ifelse(.season == "off", "green4", "dodgerblue3")
  
  # factor levels
  .df <- .df |>
    
    mutate(TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")),
           model = factor(model, labels = c("base", "functional response")))
  
  ggplot(.df) +
    
    theme_bw() +
    
    facet_grid(~ model) +
    
    geom_density(aes(x = cor,
                     linetype = sex),
                 fill = col.value,
                 color = col.value,
                 linewidth = 0.7,
                 alpha = 0.15) +
    
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          strip.background = element_rect(color = NA),
          strip.text = element_text(hjust = 0),
          axis.title.y = element_blank(),
          legend.position = c(0.15, 0.85)) +
    
    coord_cartesian(xlim = c(-0.8, 1),
                    ylim = c(0, 3)) +
    
    xlab("Spearman's correlation")
  
}

plot_corr_by_sex(off.indiv.df, "off")
plot_corr_by_sex(on.indiv.df, "on")

# median correlation
off.indiv.df |> group_by(sex, model) |>  summarize(median.cor = median(cor))
on.indiv.df |> group_by(sex, model) |>  summarize(median.cor = median(cor))
