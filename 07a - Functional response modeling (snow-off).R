# PROJECT: Habitat selection
# SCRIPT: 07a - Functional response modeling (snow-off)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 05 Jun 2026
# COMPLETED: 
# LAST MODIFIED: 16 Jun 2026
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
# 3. Examine availability distributions ----
# ______________________________________________________________________________

tsp.data <- fr.data |> group_by(TSPID) |> slice(1) 

hist(tsp.data$a.vo)
hist(tsp.data$a.stem)
hist(tsp.data$a.ch)
hist(tsp.data$a.cc)
hist(tsp.data$pOM)   
hist(tsp.data$pDM)   

# ______________________________________________________________________________
# 3. Function - Fit models ----
# ______________________________________________________________________________

fr_model <- function (.param, .avail, .smooth = "cr") {
  
  # clean data
  focal.data <- fr.data |> 
    
    filter(param == .param) |>
    
    rename("avail" = .avail) |>
    
    mutate(TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")),
           w = 1 / sd^2) |>
    
    # cluster factor
    mutate(cluster = factor(cluster))
  
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
      s(avail, k = 5, m = 1, bs = .smooth),
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
      s(avail, k = 5, m = 1, by = TRT, bs = .smooth),
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

# we need to think of these terms as part of the same parabola
test <- fr.data |> filter(param %in% c("cc", "cc2")) |> 
  
  dplyr::select(TSPID, param, beta) |>
  
  pivot_wider(names_from = param,
              values_from = beta)

plot(test$cc, test$cc2)

# STRONG correlation here
# let's try a tensor product spline

# visualize the correlation
cc.spline <- gam(
  
  cc2 ~ s(cc, bs = "cr"),
  data = test,
  family = "gaussian",
  method = "REML"
  
) |>
  
  plot()

# ______________________________________________________________________________

# clean data
cc.data <- fr.data |> 
  
  filter(param == "cc") |>  
  
  dplyr::select(TSPID, param, beta, sd, a.cc, TRT, cluster) |>
  
  mutate(TRT = factor(TRT, levels = c("UNTHIN", "RET", "PIL")),
         w = 1 / sd^2) |>
  
  # cluster factor
  mutate(cluster = factor(cluster)) |>
  
  bind_cols(
    
    fr.data |> filter(param == "cc2") |> dplyr::select(beta) |> rename(cc2 = beta)
    
  ) |>
  
  # remove outlier?
  filter(beta > 1.4)

# initialize model list
cc.model.list <- list()

# models
# M1 - NULL
cc.model.list[[1]] <- gam(
  
  beta ~ s(cluster, bs = "re") + 
    s(cc2, bs = "cr"),
  data = cc.data,
  family = "gaussian",
  method = "REML",
  weights = w
  
)

# M2 - FR
cc.model.list[[2]] <- gam(
  
  beta ~ s(cluster, bs = "re") +
    te(cc2, a.cc, bs = "cr"),
  data = cc.data,
  family = "gaussian",
  method = "REML",
  weights = w
  
)

# M3 - TRT
cc.model.list[[3]] <- gam(
  
  beta ~ s(cluster, bs = "re") +
    s(cc2, by = TRT, bs = "cr"),
  data = cc.data,
  family = "gaussian",
  method = "REML",
  weights = w
  
)

# M4 - FR x TRT
cc.model.list[[4]] <- gam(
  
  beta ~ s(cluster, bs = "re") +
    te(cc2, a.cc, by = TRT, bs = "cr"),
  data = cc.data,
  family = "gaussian",
  method = "REML",
  weights = w
  
)

# AIC table
aic_tab(cc.model.list) # M3

summary(cc.model.list[[3]])
plot(cc.model.list[[3]])
appraise(cc.model.list[[3]])

# ______________________________________________________________________________
# 5d. dOM ----
# ______________________________________________________________________________

models.dOM <- fr_model("dOM", "pOM")

aic_tab(models.dOM)  # M3

summary(models.dOM[[3]])
#plot(models.dOM[[3]])
appraise(models.dOM[[3]])

# ______________________________________________________________________________
# 5e. dDM ----
# ______________________________________________________________________________

models.dDM <- fr_model("dDM", "pDM")

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
#appraise(models.ed[[1]])

# ______________________________________________________________________________
# 6. Save top models (if not null) ----
# ______________________________________________________________________________

saveRDS(models.vo[[4]], "model_results/fr_models/off_vo.rds")
saveRDS(models.dOM[[3]], "model_results/fr_models/off_dOM.rds")
saveRDS(models.dDM[[3]], "model_results/fr_models/off_dDM.rds")

saveRDS(cc.model.list[[3]], "model_results/fr_models/off_cc.rds")
