# PROJECT: Habitat selection
# SCRIPT: 09 - Assess predictive ability
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 09 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 20 Jul 2026
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

# FR models
# off 
off.vo <- readRDS("model_results/fr_models/off_vo.rds")
off.cc <- readRDS("model_results/fr_models/off_cc.rds")
off.dEdge <- readRDS("model_results/fr_models/off_dEdge.rds")

# on
on.cc <- readRDS("model_results/fr_models/on_cc.rds")
on.dEdge <- readRDS("model_results/fr_models/on_dEdge.rds")

# data
off.fr.data <- readRDS("data_for_model/off_fr.rds")
on.fr.data <- readRDS("data_for_model/on_fr.rds")

# HS data (bind in site column)
off.hs.data <- readRDS("data_for_model/off_data.rds") 
on.hs.data <- readRDS("data_for_model/on_data.rds")

# bind in site and sex columns
off.hs.data <- off.hs.data |>
  
  left_join(
  
    readRDS("data_cleaned/data_off.rds") |> 
    
    dplyr::select(site, sex, track_season_post) |>
    
    group_by(track_season_post) |>
    
    slice(1)
  
  )

on.hs.data <- on.hs.data |>
  
  left_join(
    
    readRDS("data_cleaned/data_on.rds") |> 
      
      dplyr::select(site, sex, track_season_post) |>
      
      group_by(track_season_post) |>
      
      slice(1)
    
  )

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
      
      # add TRT and cluster
      mutate(TRT = c.trt,
             cluster = substr(site, 1, 1)) |>
      
      dplyr::select(track_season_post, case, TRT, cluster, 
                    vo, ch, cc, twi, twi2, vrm, vrm2, dEdge) |>
      
      rename(TSPID = track_season_post)
    
    # add residuals
    hs.data$g.s <- residuals(lm(log(akde) ~ 
                                  twi + twi2 + vrm + vrm2 + 
                                  vo + ch + cc + dEdge,
                                data = off.hs.data))
    
    # base coefficients
    coef.base <- M.off[[1]] |> dplyr::select(param, mean)
    
    # FR coefficients
    coef.fr <- off.fr.data |> dplyr::select(TSPID, TRT, cluster, a.vo) |>
      
      group_by(TSPID) |> slice(1) |> ungroup() |>
      
      mutate(
        
        beta.vo = predict(off.vo, data.frame(TRT = TRT, cluster = cluster, avail = a.vo)),
        beta.cc = predict(off.cc, data.frame(TRT = TRT, cluster = cluster, avail = a.vo)),
        beta.dEdge = predict(off.dEdge, data.frame(TRT = TRT, cluster = cluster, avail = a.vo))

        
      ) |>
      
      dplyr::select(TSPID, beta.vo, beta.cc, beta.dEdge)
    
    # add implied coefs to all use/avail
    hs.data.fr <- hs.data |>
      
      dplyr::select(TSPID) |>
      
      left_join(coef.fr)
    
    # calculate w(x)
    hs.data$pred.base <- exp(
      
      hs.data$g.s * coef.base$mean[coef.base$param == "g.s"] +
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
      hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
      hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
      hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
      hs.data$vo * coef.base$mean[coef.base$param == "vo"] +
      hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
      hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
      hs.data$dEdge * coef.base$mean[coef.base$param == "dEdge"]
      
    )
    
    hs.data$pred.fr <- exp(
      
      hs.data$g.s * coef.base$mean[coef.base$param == "g.s"] +
      hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
      hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
      hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
      hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
      hs.data$vo * hs.data.fr$beta.vo +
      hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
      hs.data$cc * hs.data.fr$beta.cc +
      hs.data$dEdge * hs.data.fr$beta.dEdge
      
    )
    
  } # off
  
  if (.hsf == "M.on") {
    
    # use/availability data
    hs.data <- on.hs.data |>
      
      # add TRT and cluster
      mutate(TRT = c.trt,
             cluster = substr(site, 1, 1)) |>
      
      dplyr::select(track_season_post, case, TRT, cluster, 
                    stem, ch, cc, twi, twi2, vrm, vrm2, dEdge) |>
      
      rename(TSPID = track_season_post)
    
    # add residuals
    hs.data$g.s <- residuals(lm(log(akde) ~ 
                                  twi + twi2 + vrm + vrm2 + 
                                  vo + ch + cc + dEdge,
                                data = on.hs.data))
    
    # base coefficients
    coef.base <- M.on[[1]] |> dplyr::select(param, mean)
    
    # FR coefficients
    coef.fr <- on.fr.data |> dplyr::select(TSPID, TRT, cluster, a.stem) |>
      
      group_by(TSPID) |> slice(1) |> ungroup() |>
      
      mutate(
        
        beta.cc = predict(on.cc, data.frame(TRT = TRT, cluster = cluster, avail = a.stem)),
        beta.dEdge = predict(on.dEdge, data.frame(TRT = TRT, cluster = cluster, avail = a.stem))
        
        
      ) |>
      
      dplyr::select(TSPID, beta.cc, beta.dEdge)
    
    # add implied coefs to all use/avail
    hs.data.fr <- hs.data |>
      
      dplyr::select(TSPID) |>
      
      left_join(coef.fr)
    
    # calculate w(x)
    hs.data$pred.base <- exp(
      
      hs.data$g.s * coef.base$mean[coef.base$param == "g.s"] +
        hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$stem * coef.base$mean[coef.base$param == "stem"] +
        hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
        hs.data$cc * coef.base$mean[coef.base$param == "cc"] +
        hs.data$dEdge * coef.base$mean[coef.base$param == "dEdge"]
      
    )
    
    hs.data$pred.fr <- exp(
      
      hs.data$g.s * coef.base$mean[coef.base$param == "g.s"] +
        hs.data$twi * coef.base$mean[coef.base$param == "twi"] +
        hs.data$twi2 * coef.base$mean[coef.base$param == "twi2"] +
        hs.data$vrm * coef.base$mean[coef.base$param == "vrm"] +
        hs.data$vrm2 * coef.base$mean[coef.base$param == "vrm2"] +
        hs.data$stem * coef.base$mean[coef.base$param == "stem"] +
        hs.data$ch * coef.base$mean[coef.base$param == "ch"] +
        hs.data$cc * hs.data.fr$beta.cc +
        hs.data$dEdge * hs.data.fr$beta.dEdge
      
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
        model = "base"),
  
  cbind(boyce("M.off", "pop", "fr"),
        model = "fr")
  
)

on.pop <- rbind(
  
  cbind(boyce("M.on", "pop", "base"),
        model = "base"),
  
  cbind(boyce("M.on", "pop", "fr"),
        model = "fr")
  
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
off.pop.base <- off.pop |> filter(model == "base")
off.pop.fr <- off.pop |> filter(model == "fr")

cor.test(off.pop.base$bins, off.pop.base$u.a.ratio)$estimate
cor.test(off.pop.fr$bins, off.pop.fr$u.a.ratio)$estimate

# on
on.pop.base <- on.pop |> filter(model == "base")
on.pop.fr <- on.pop |> filter(model == "fr")

cor.test(on.pop.base$bins, on.pop.base$u.a.ratio)$estimate
cor.test(on.pop.fr$bins, on.pop.fr$u.a.ratio)$estimate

# ______________________________________________________________________________
# 4b. Individual-level ----
# ______________________________________________________________________________

# helper function - calculate correlation
indiv_corr <- function (x) {
  
  if (nrow(x) > 1) {
  
  corr <- cor(x$bins, x$u.a.ratio)
  
  return(corr)
  
  }
  
}

# function to do it all self-contained
indiv_corr_df <- function (.boyce,
                           .season,
                           .model) {
  
  .boyce.df <- lapply(split(.boyce, ~TSPID), indiv_corr) |>
    
    do.call(rbind, args = _) |>
    
    as.data.frame() |>
    
    rename(cor = V1) |>
    
    mutate(season = .season,
           model = .model)
  
  # add rownames as TSP
  .boyce.df$TSPID <- rownames(.boyce.df)
  
  return(.boyce.df)
  
}

# off
off.indiv.base.corr <- indiv_corr_df(off.indiv.base, "off", "base")
off.indiv.fr.corr <- indiv_corr_df(off.indiv.fr, "off", "fr")
on.indiv.base.corr <- indiv_corr_df(on.indiv.base, "on", "base")
on.indiv.fr.corr <- indiv_corr_df(on.indiv.fr, "on", "fr")

# ______________________________________________________________________________
# 5. U/A vs bin plots ---
# ______________________________________________________________________________
# 5a. Prepare data ----
# ______________________________________________________________________________

all.pop <- bind_rows(
  
  off.pop |> mutate(season = "off"),
  on.pop |> mutate(season = "on")
  
) |>
  
  # factors
  mutate(
    
    model = factor(model, labels = c("base", "functional")),
    season = factor(season, labels = c("snow-off", "snow-on"))
    
  )

all.indiv <- bind_rows(
  
  off.indiv.base |> left_join(off.indiv.base.corr),
  off.indiv.fr |> left_join(off.indiv.fr.corr),
  on.indiv.base |> left_join(on.indiv.base.corr),
  on.indiv.fr |> left_join(on.indiv.fr.corr)
  
) |>
  
  # factors
  mutate(
    
    model = factor(model, labels = c("base", "functional")),
    season = factor(season, labels = c("snow-off", "snow-on"))
    
  )

# ______________________________________________________________________________
# 5b. Plot ----
# ______________________________________________________________________________

ggplot() +
  
  theme_bw() +
  
  facet_grid(season ~ model) +
  
  # individual
  geom_line(data = all.indiv,
            aes(x = bins,
                y = u.a.ratio,
                group = TSPID,
                color = cor),
            linewidth = 0.2,
            alpha = 0.5) +
  
  # population
  geom_line(data = all.pop,
            aes(x = bins,
                y = u.a.ratio),
            color = "black",
            linewidth = 0.9) +
  
  theme(panel.grid = element_blank(),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA),
        axis.text = element_text(color = "black")) +
  
  scale_color_viridis_c("Correlation") +
  
  scale_x_continuous(breaks = c(1:10)) +
  
  xlab("HSF score bin") +
  ylab("Used/available ratio")

# ______________________________________________________________________________
# 6. Individual comparison plots ----
# ______________________________________________________________________________
# 6a. Prep data ----
# ______________________________________________________________________________

# add identifiers with all indivs
indiv.corr.id <- all.indiv |>
  
  group_by(TSPID, model) |>
  
  slice(1) |>
  
  dplyr::select(TSPID, cor, season, model) |>
  
  left_join(
    
    bind_rows(
      
      readRDS("data_cleaned/data_off.rds") |> 
        
        dplyr::select(site, sex, track_season_post, c.trt) |>
        
        group_by(track_season_post) |>
        
        slice(1) |>
        
        rename(TSPID = track_season_post) |>
        
        mutate(season = "snow-off"),
      
      readRDS("data_cleaned/data_on.rds") |> 
        
        dplyr::select(site, sex, track_season_post, c.trt) |>
        
        group_by(track_season_post) |>
        
        slice(1) |>
        
        rename(TSPID = track_season_post) |>
        
        mutate(season = "snow-on")
      
    )
    
  ) |>
  
  # add cluster
  mutate(cluster = substr(site, 1, 1))

# ______________________________________________________________________________
# 6b. Histograms - season x model ----

# calculate mean correlations
indiv.corr.means.season.model <- indiv.corr.id |>
  
  group_by(season, model) |>
  
  summarize(mean.cor = mean(cor))

# ______________________________________________________________________________

ggplot(data = indiv.corr.id) +
  
  theme_bw() +
  
  facet_grid(season ~ model) +
  
  geom_histogram(aes(x = cor,
                     fill = season),
                 color = "white") +
  
  # means
  geom_vline(data = indiv.corr.means.season.model,
             aes(xintercept = mean.cor),
             linetype = "dashed",
             linewidth = 0.8) +
  
  theme(panel.grid = element_blank(),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  
  scale_fill_manual(values = c("green4", "dodgerblue2")) +
  
  xlab("Spearman's correlation") +
  ylab("Individual tracks")

# ______________________________________________________________________________
# 6c. Performance comparison - season x model ----
# ______________________________________________________________________________

ggplot(data = indiv.corr.id |> pivot_wider(names_from = model,
                                           values_from = cor)) +
  
  theme_bw() +
  
  facet_grid(~ season) +
  
  # 1:1 line
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +

  geom_point(aes(x = base,
                 y = functional,
                 color = season),
             alpha = 0.5) +
  
  theme(panel.grid = element_blank(),
        strip.text = element_text(hjust = 0),
        strip.background = element_rect(color = NA),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  
  scale_color_manual(values = c("green4", "dodgerblue2")) +
  
  coord_cartesian(xlim = c(-0.1, 1),
                  ylim = c(-0.1, 1)) +
  
  xlab("Base correlation") +
  ylab("Functional correlation")
