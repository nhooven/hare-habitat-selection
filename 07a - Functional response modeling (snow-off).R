# PROJECT: Habitat selection
# SCRIPT: 07a - Functional response modeling (snow-off)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 05 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 09 Jun 2026
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

fr.data <- readRDS("data_for_model/off_fr.rds")

# ______________________________________________________________________________
# 3. Function - Fit models ----
# ______________________________________________________________________________

fr_model <- function (.param, .avail, .smooth = "cs") {
  
  # clean data
  focal.data <- fr.data |> 
    
    filter(param == .param) |>
    
    rename("avail" = .avail) |>
    
    mutate(TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")),
           w = 1 / sd^2)
  
  # initialize model list
  model.list <- list()
  
  # models
  # M1 - NULL
  model.list[[1]] <- gam(
    
    beta ~ 1,
    data = focal.data,
    family = "gaussian",
    method = "REML",
    weights = w
    
  )
  
  # M2 - FR
  model.list[[2]] <- gam(
    
    beta ~ s(avail, bs = .smooth),
    data = focal.data,
    family = "gaussian",
    method = "REML",
    weights = w
    
  )
  
  # M3 - TRT
  model.list[[3]] <- gam(
    
    beta ~ TRT,
    data = focal.data,
    family = "gaussian",
    method = "REML",
    weights = w
    
  )
  
  # M4 - FR x TRT
  model.list[[4]] <- gam(
    
    beta ~ s(avail, by = TRT, bs = .smooth),
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
# 5a. VO ----
# ______________________________________________________________________________

models.vo <- fr_model("vo", "a.vo")

aic_tab(models.vo)  # M4

summary(models.vo[[4]])
plot(models.vo[[4]])
appraise(models.vo[[4]])

# ______________________________________________________________________________
# 5b. CH ----
# ______________________________________________________________________________

models.ch <- fr_model("ch", "a.ch")

aic_tab(models.ch)  # M1

summary(models.ch[[1]])
#plot(models.ch[[1]])
appraise(models.ch[[1]])

# ______________________________________________________________________________
# 5c. CC and CC2 ----
# ______________________________________________________________________________

models.cc <- fr_model("cc", "a.cc")

aic_tab(models.cc)  # M1

models.cc2 <- fr_model("cc2", "a.cc")

aic_tab(models.cc2)  # M1

# ______________________________________________________________________________
# 5d. dOM ----
# ______________________________________________________________________________

models.dOM <- fr_model("dOM", "a.dOM")

aic_tab(models.dOM)  # M4

summary(models.dOM[[4]])
plot(models.dOM[[4]])
appraise(models.dOM[[4]])

# ______________________________________________________________________________
# 5e. dDM ----
# ______________________________________________________________________________

models.dDM <- fr_model("dDM", "a.dDM")

aic_tab(models.dDM)  # M3

summary(models.dDM[[3]])
#plot(models.dDM[[3]])
appraise(models.dDM[[3]])

# ______________________________________________________________________________
# 5f. ed ----
# ______________________________________________________________________________

models.ed <- fr_model("ed", "a.ed")

aic_tab(models.ed)  # M1

#summary(models.ed[[1]])
#plot(models.ed[[1]])
#appraise(models.ed[[1]])

# ______________________________________________________________________________
# 5g. shdi ----
# ______________________________________________________________________________

models.shdi <- fr_model("shdi", "a.shdi")

aic_tab(models.shdi)  # M4

summary(models.shdi[[4]])
plot(models.shdi[[4]])
appraise(models.shdi[[4]])

# ______________________________________________________________________________
# 6. Save top models (if not null) ----
# ______________________________________________________________________________

saveRDS(models.vo[[4]], "model_results/fr_models/off_vo.rds")
saveRDS(models.dOM[[4]], "model_results/fr_models/off_dOM.rds")
saveRDS(models.dDM[[3]], "model_results/fr_models/off_dDM.rds")
saveRDS(models.shdi[[4]], "model_results/fr_models/off_shdi.rds")
