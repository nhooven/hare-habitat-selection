# PROJECT: Habitat selection
# SCRIPT: 06a - Modeling (snow-off)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 
# LAST MODIFIED: 27 May 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(INLA)
library(tictoc)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

data.off <- readRDS("data_for_model/off_data.rds")

# ______________________________________________________________________________
# 3. Correlation ----
# ______________________________________________________________________________

covs.lin <- c("dOpen", "dDM", "ch", "cc", "vo", "twi", "vrm")

cor(data.off |> dplyr::select(covs.lin), method = "spearman")

# ______________________________________________________________________________
# 3. Setup ----

# we'll fit 11 random slopes

# ______________________________________________________________________________

# copy MRID
data.off <- data.off |>
  
  mutate(
    
    MRID1 = MRID,
    MRID2 = MRID,
    MRID3 = MRID,
    MRID4 = MRID,
    MRID5 = MRID,
    MRID6 = MRID,
    MRID7 = MRID,
    MRID8 = MRID,
    MRID9 = MRID,
    MRID10 = MRID,
    MRID11 = MRID
    
  )

# ______________________________________________________________________________
# 3. Model formulae ----

hyper.list <- list(theta = list(initial = log(1),
                                fixed = F,
                                prior = "pc.prec",
                                param = c(1, 0.05)))

# ______________________________________________________________________________
# 3a. M1 - Base model ----
# ______________________________________________________________________________

M1.form <- case ~ 
  
  # population-level effects
  # structural (n = 5)
  vo +
  ch + ch2 +
  cc + cc2 +
  
  # conditions (n = 4)
  twi + twi2 +
  vrm + vrm2 +
  
  # landscape (n = 2)
  dOpen +
  dDM +
  
  # random intercepts
  f(MRID, model = "iid", hyper = list(theta = list(initial = log(1/1e6), fixed = T))) +
  
  # random slopes
  f(MRID1, vo, values = unique(MRID), model = "iid", hyper = hyper.list) +
  f(MRID2, ch, model = "iid", hyper = hyper.list) +
  f(MRID3, ch2, model = "iid", hyper = hyper.list) + 
  f(MRID4, cc, model = "iid", hyper = hyper.list) + 
  f(MRID5, cc2, model = "iid", hyper = hyper.list) +
  f(MRID6, twi, model = "iid", hyper = hyper.list) +
  f(MRID7, twi2, model = "iid", hyper = hyper.list) +
  f(MRID8, vrm, model = "iid", hyper = hyper.list) +
  f(MRID9, vrm2, model = "iid", hyper = hyper.list) +
  f(MRID10, dOpen, model = "iid", hyper = hyper.list) +
  f(MRID11, dDM, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3b. M2 - Functional response model ----
# ______________________________________________________________________________

M2.form <- case ~ 
  
  # population-level effects
  # structural (n = 5)
  vo +
  ch + ch2 +
  cc + cc2 +
  
  # conditions (n = 4)
  twi + twi2 +
  vrm + vrm2 +
  
  # landscape (n = 2)
  dOpen +
  dDM +
  
  # functional responses (n = 5)
  vo : a.vo +
  ch : a.ch + ch2 : a.ch +
  cc : a.cc + cc2 : a.cc +
  
  # random intercepts
  f(MRID, model = "iid", hyper = list(theta = list(initial = log(1/1e6), fixed = T))) +
  
  # random slopes
  f(MRID1, vo, model = "iid", hyper = hyper.list) +
  f(MRID2, ch, model = "iid", hyper = hyper.list) +
  f(MRID3, ch2, model = "iid", hyper = hyper.list) + 
  f(MRID4, cc, model = "iid", hyper = hyper.list) + 
  f(MRID5, cc2, model = "iid", hyper = hyper.list) +
  f(MRID6, twi, model = "iid", hyper = hyper.list) +
  f(MRID7, twi2, model = "iid", hyper = hyper.list) +
  f(MRID8, vrm, model = "iid", hyper = hyper.list) +
  f(MRID9, vrm2, model = "iid", hyper = hyper.list) +
  f(MRID10, dOpen, model = "iid", hyper = hyper.list) +
  f(MRID11, dDM, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3c. M3 - Treatment model ----
# ______________________________________________________________________________

M3.form <- case ~ 
  
  # population-level effects
  # structural (n = 5)
  vo +
  ch + ch2 +
  cc + cc2 +
  
  # conditions (n = 4)
  twi + twi2 +
  vrm + vrm2 +
  
  # landscape (n = 2)
  dOpen +
  dDM +
  
  # treatment (n = 6)
  vo : year.trt : ret +
  vo : year.trt : pil +
  ch : year.trt : ret +
  ch : year.trt : pil +
  ch2 : year.trt : ret +
  ch2 : year.trt : pil +
  cc : year.trt : ret +
  cc : year.trt : pil +
  cc2 : year.trt : ret +
  cc2 : year.trt : pil +
  
  # random intercepts
  f(MRID, model = "iid", hyper = list(theta = list(initial = log(1/1e6), fixed = T))) +
  
  # random slopes
  f(MRID1, vo, model = "iid", hyper = hyper.list) +
  f(MRID2, ch, model = "iid", hyper = hyper.list) +
  f(MRID3, ch2, model = "iid", hyper = hyper.list) + 
  f(MRID4, cc, model = "iid", hyper = hyper.list) + 
  f(MRID5, cc2, model = "iid", hyper = hyper.list) +
  f(MRID6, twi, model = "iid", hyper = hyper.list) +
  f(MRID7, twi2, model = "iid", hyper = hyper.list) +
  f(MRID8, vrm, model = "iid", hyper = hyper.list) +
  f(MRID9, vrm2, model = "iid", hyper = hyper.list) +
  f(MRID10, dOpen, model = "iid", hyper = hyper.list) +
  f(MRID11, dDM, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3d. M4 - Functional x treatment responses model ----
# ______________________________________________________________________________

M4.form <- case ~ 
  
  # population-level effects
  # structural (n = 5)
  vo +
  ch + ch2 +
  cc + cc2 +
  
  # conditions (n = 4)
  twi + twi2 +
  vrm + vrm2 +
  
  # landscape (n = 2)
  dOpen +
  dDM +
  
  # treatment (n = 6)
  vo : year.trt : ret : a.vo +
  vo : year.trt : pil : a.vo +
  ch : year.trt : ret : a.ch +
  ch : year.trt : pil : a.ch +
  ch2 : year.trt : ret : a.ch +
  ch2 : year.trt : pil : a.ch +
  cc : year.trt : ret : a.cc +
  cc : year.trt : pil : a.cc +
  cc2 : year.trt : ret : a.cc +
  cc2 : year.trt : pil : a.cc +
  
  # random intercepts
  f(MRID, model = "iid", hyper = list(theta = list(initial = log(1/1e6), fixed = T))) +
  
  # random slopes
  f(MRID1, vo, model = "iid", hyper = hyper.list) +
  f(MRID2, ch, model = "iid", hyper = hyper.list) +
  f(MRID3, ch2, model = "iid", hyper = hyper.list) + 
  f(MRID4, cc, model = "iid", hyper = hyper.list) + 
  f(MRID5, cc2, model = "iid", hyper = hyper.list) +
  f(MRID6, twi, model = "iid", hyper = hyper.list) +
  f(MRID7, twi2, model = "iid", hyper = hyper.list) +
  f(MRID8, vrm, model = "iid", hyper = hyper.list) +
  f(MRID9, vrm2, model = "iid", hyper = hyper.list) +
  f(MRID10, dOpen, model = "iid", hyper = hyper.list) +
  f(MRID11, dDM, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 4. Fit models ----
# ______________________________________________________________________________
# 4a. M1 - Base model ----

# 151 s

# ______________________________________________________________________________

tic()
M1.fit <- inla(M1.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = list(cpo = T,
                                      dic = T))
toc()
  
  
summary(M1.fit)

# 27 May 2026
# Fit all models with CPO/PIT (adds a little bit of time)

# save base model for visualization
saveRDS(M1.fit, "model_tests/M1_fit.rds")

# ______________________________________________________________________________
# 4b. M2 - Functional responses ----

# 193 s

# ______________________________________________________________________________

tic()
M2.fit <- inla(M2.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = list(cpo = T,
                                      dic = T))
toc()


summary(M2.fit)

# ______________________________________________________________________________
# 4c. M3 - Treatment responses ----

# 180 s

# ______________________________________________________________________________

tic()
M3.fit <- inla(M3.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = list(cpo = T,
                                      dic = T))
toc()


summary(M3.fit)

# ______________________________________________________________________________
# 4d. M4 - Functional x treatment responses ----

# 209 s

# ______________________________________________________________________________

tic()
M4.fit <- inla(M4.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = list(cpo = T,
                                      dic = T))
toc()


summary(M4.fit)

# ______________________________________________________________________________
# 05. Model selection with DIC ----

# function
inla_dic_table <- function (x) {
  
  # x is a list of models
  dic.table <- data.frame()
  
  for (i in 1:length(x)) {
    
    model.focal <- x[[i]]
    
    dic.table.focal <- data.frame(
      
      Model = case_when(i == 1 ~ "base",
                        i == 2 ~ "fr",
                        i == 3 ~ "trt",
                        i == 4 ~ "fr x trt"),
      k.eff = model.focal$dic$p.eff,
      DIC = model.focal$dic$dic,
      mLL = model.focal$mlik[2]
      
    )
    
    # bind in 
    dic.table <- rbind(dic.table, dic.table.focal)
    
  }
  
  # compute dWAIC and arrange
  dic.table <- dic.table |>
    
    mutate(dDIC = DIC - min(DIC)) |>
    
    # arrange
    arrange(dDIC) |>
    
    dplyr::select(Model, k.eff, DIC, dDIC, mLL)
  
  return(dic.table)
  
}

# ______________________________________________________________________________

inla_dic_table(list(M1.fit,
                    M2.fit,
                    M3.fit,
                    M4.fit))
