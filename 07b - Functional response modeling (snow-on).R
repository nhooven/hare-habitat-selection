# PROJECT: Habitat selection
# SCRIPT: 07b - Functional response modeling (snow-on)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 05 Jun 2026
# COMPLETED: 18 Jun 2026
# LAST MODIFIED: 18 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(mgcv)
library(gratia)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

fr.data <- readRDS("data_for_model/on_fr.rds")

# ______________________________________________________________________________
# 3. Function - Fit models ----
# ______________________________________________________________________________

fr_model <- function (.param, .avail, .smooth = "cs", .k = 10) {
  
  # clean data
  focal.data <- fr.data |> 
    
    filter(param == .param) |>
    
    rename("avail" = .avail) |>
    
    # factors and weights
    mutate(cluster = factor(cluster),
           TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")),
           w = 1 / sd^2)
  
  # initialize model list
  model.list <- list()
  
  # models
  # M1 - NULL
  model.list[[1]] <- gam(
    
    beta ~ s(cluster, bs = "re"),
    data = focal.data,
    family = "gaussian",
    method = "REML",
    weights = w
    
  )
  
  # M2 - FR
  model.list[[2]] <- gam(
    
    beta ~ s(cluster, bs = "re") +
      s(avail, k = .k, bs = .smooth),
    data = focal.data,
    family = "gaussian",
    method = "REML",
    weights = w
    
  )
  
  # M3 - TRT
  model.list[[3]] <- gam(
    
    beta ~ TRT + s(cluster, bs = "re"),
    data = focal.data,
    family = "gaussian",
    method = "REML",
    weights = w
    
  )
  
  # M4 - FR x TRT
  model.list[[4]] <- gam(
    
    beta ~ s(cluster, bs = "re") +
      s(avail, k = .k, by = TRT, bs = .smooth),
    data = focal.data,
    family = "gaussian",
    method = "REML",
    weights = w
    
  )
  
  return(model.list)
  
}

# ______________________________________________________________________________
# 4. Function - AIC table ----
# ______________________________________________________________________________

aic_tab <- function(.models) {
  
  # model names
  model.names <- c("M1 - NULL",
                   "M2 - FR",
                   "M3 - TRT",
                   "M4 - FR x TRT")
  
  # loop through each model
  aic.tab <- data.frame()
  
  for (i in 1:length(.models)) {
    
    focal.model <- .models[[i]]
    
    # extract AIC, estimated degrees of freedom, and number of observations
    aic <- round(focal.model$aic, digits = 1)
    edf <- round(sum(focal.model$edf), digits = 1)
    
    # create df to hold everything
    focal.row <- data.frame("model" = model.names[i],
                            "edf" = edf,
                            "AIC" = aic,
                            "LL" = as.numeric(logLik.gam(focal.model)))
    
    # bind into table
    aic.tab <- rbind(aic.tab, focal.row)
    
  }
  
  # sort table and calculate delta AIC
  aic.tab.1 <- aic.tab %>%
    
    arrange(AIC) %>%
    
    mutate(dAIC = AIC - min(AIC)) %>%
    
    dplyr::select(model, edf, AIC, dAIC, LL)
  
  # return
  return(aic.tab.1)
  
}

# ______________________________________________________________________________
# 5. Fit models ----
# ______________________________________________________________________________
# 5a. Stem ----
# ______________________________________________________________________________

models.stem <- fr_model("stem", "a.stem")

aic_tab(models.stem)  # M2

summary(models.stem[[2]])
plot(models.stem[[2]])
appraise(models.stem[[2]])

# ______________________________________________________________________________
# 5b. CH ----
# ______________________________________________________________________________

models.ch <- fr_model("ch", "a.stem")

aic_tab(models.ch)  # M3

summary(models.ch[[3]])
appraise(models.ch[[3]])

# ______________________________________________________________________________
# 5c. dEdge ----
# ______________________________________________________________________________

models.dEdge <- fr_model("dEdge", "a.stem")

aic_tab(models.dEdge)  # M4

summary(models.dEdge[[4]])
plot(models.dEdge[[4]])
appraise(models.dEdge[[4]])

# ______________________________________________________________________________
# 6. Save top models (if not null) ----
# ______________________________________________________________________________

saveRDS(models.stem[[2]], "model_results/fr_models/on_stem.rds")
saveRDS(models.ch[[3]], "model_results/fr_models/on_ch.rds")
saveRDS(models.dEdge[[4]], "model_results/fr_models/on_dEdge.rds")
