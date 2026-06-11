# PROJECT: Habitat selection
# SCRIPT: 07b - Functional response modeling (snow-on)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 05 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 11 Jun 2026
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
# 5a. Stem ----
# ______________________________________________________________________________

models.stem <- fr_model("stem", "a.stem")

aic_tab(models.stem)  # M4

summary(models.stem[[4]])
plot(models.stem[[4]])
appraise(models.stem[[4]])

# ______________________________________________________________________________
# 5b. CH ----
# ______________________________________________________________________________

models.ch <- fr_model("ch", "a.ch")

aic_tab(models.ch)  # M4

summary(models.ch[[4]])
plot(models.ch[[4]])
appraise(models.ch[[4]])

# ______________________________________________________________________________
# 5c. CC and CC2 ----
# ______________________________________________________________________________

models.cc <- fr_model("cc", "a.cc")

aic_tab(models.cc)  # M1

models.cc2 <- fr_model("cc2", "a.cc")

aic_tab(models.cc2)  # M4

summary(models.cc2[[4]])
plot(models.cc2[[4]])
appraise(models.cc2[[4]])

# ______________________________________________________________________________
# 5d. dOM ----

# weird outlier leading to a ridiculous relationship
test <- fr.data |> filter(param == "dOM")

plot(test$a.dOM, test$beta)

# let's remove it for the model
fr.data <- fr.data |> filter(beta < 2.5)

# ______________________________________________________________________________

models.dOM <- fr_model("dOM", "a.dOM")

aic_tab(models.dOM)  # M2

summary(models.dOM[[2]])
plot(models.dOM[[2]])
appraise(models.dOM[[2]])

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

summary(models.ed[[1]])
#plot(models.ed[[1]])
appraise(models.ed[[1]])

# ______________________________________________________________________________
# 5g. shdi ----
# ______________________________________________________________________________

models.shdi <- fr_model("shdi", "a.shdi")

aic_tab(models.shdi)  # M2

summary(models.shdi[[2]])
plot(models.shdi[[2]])
appraise(models.shdi[[2]])

# ______________________________________________________________________________
# 6. Save top models (if not null) ----
# ______________________________________________________________________________

saveRDS(models.stem[[4]], "model_results/fr_models/on_stem.rds")
saveRDS(models.ch[[4]], "model_results/fr_models/on_ch.rds")
saveRDS(models.cc2[[4]], "model_results/fr_models/on_cc2.rds")
saveRDS(models.dOM[[2]], "model_results/fr_models/on_dOM.rds")
saveRDS(models.dDM[[3]], "model_results/fr_models/on_dDM.rds")
saveRDS(models.shdi[[2]], "model_results/fr_models/on_shdi.rds")
