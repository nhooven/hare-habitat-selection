# PROJECT: Habitat selection
# SCRIPT: 06a - Modeling (snow-off)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 
# LAST MODIFIED: 02 Jun 2026
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

# residuals for log(AKDE)
data.off$g.s <- residuals(lm(log(akde) ~ vo + ch + ch2 + cc + cc2 + 
                               twi + twi2 + vrm + vrm2 + 
                               dOpen + dDM + ed,
                             data = data.off))

# ______________________________________________________________________________
# 3. Correlation ----
# ______________________________________________________________________________

covs.lin <- c("vo", "ch", "cc", "twi", "vrm", "dOpen", "dDM", "ed", "g.s")

cor(data.off |> dplyr::select(covs.lin), method = "spearman")

# ______________________________________________________________________________
# 3. Setup ----

# we'll fit 12 random slopes

# ______________________________________________________________________________

# copy MRID
data.off <- data.off |>
  
  mutate(
    
    MRID = MRID,
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
    MRID11 = MRID,
    MRID12 = MRID
    
  )

# ______________________________________________________________________________
# 3. Model formulae ----

hyper.list <- list(theta = list(initial = log(1),
                                fixed = F,
                                prior = "pc.prec",
                                param = c(3, 0.05)))

# RS structure
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
  f(MRID11, dDM, model = "iid", hyper = hyper.list) +
  f(MRID12, ed, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3a. M1 - STAND + (LAND x AVAIL) + TOPO ----
# ______________________________________________________________________________

M1.form <- case ~ 
  
  # population-level effects
  g.s +
  
  # STAND
  vo +
  ch + ch2 +
  cc + cc2 +
  
  # LAND
  dOpen + dOpen : a.dOpen +
  dDM + dDM : a.dDM +
  ed + ed : a.ed +
  
  # TOPO
  twi + twi2 + vrm + vrm2 +
  
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
  f(MRID11, dDM, model = "iid", hyper = hyper.list) +
  f(MRID12, ed, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3b. M2 - (STAND x AVAIL) + (LAND x AVAIL) + TOPO ----
# ______________________________________________________________________________

M2.form <- case ~ 
  
  # population-level effects
  g.s +
  
  # STAND
  vo + vo : a.vo +
  ch + ch : a.ch + ch2 + ch2 : a.ch +
  ch + cc : a.cc + cc2 + cc2 : a.cc +
  
  # LAND
  dOpen + dOpen : a.dOpen +
  dDM + dDM : a.dDM +
  ed + ed : a.ed +
  
  # TOPO
  twi + twi2 + vrm + vrm2 +
  
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
  f(MRID11, dDM, model = "iid", hyper = hyper.list) +
  f(MRID12, ed, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3c. M3 - (STAND X TRT) + (LAND x AVAIL) + TOPO ----
# ______________________________________________________________________________

M3.form <- case ~ 
  
  # population-level effects
  g.s +
  
  # STAND
  vo + 
    vo : year.trt : ret + 
    vo : year.trt : pil +
  ch + 
    ch : year.trt : ret + 
    ch : year.trt : pil + 
  ch2 + 
    ch2 : year.trt : ret + 
    ch2 : year.trt : pil + 
  cc + 
    cc : year.trt : ret + 
    cc : year.trt : pil + 
  cc2 + 
    cc2 : year.trt : ret + 
    cc2 : year.trt : pil + 
  
  # LAND
  dOpen + dOpen : a.dOpen +
  dDM + dDM : a.dDM +
  ed + ed : a.ed +
  
  # TOPO
  twi + twi2 + vrm + vrm2 +
  
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
  f(MRID11, dDM, model = "iid", hyper = hyper.list) +
  f(MRID12, ed, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3d. M4 - STAND + (LAND x AVAIL x TRT) + TOPO ----
# ______________________________________________________________________________

M4.form <- case ~ 
  
  # population-level effects
  g.s +
  
  # STAND
  vo +
  ch + ch2 +
  cc + cc2 +
  
  # LAND
  dOpen + dOpen : a.dOpen + 
    dOpen : year.trt : ret +
    dOpen : year.trt : pil +
    dOpen : a.dOpen : year.trt : ret +
    dOpen : a.dOpen : year.trt : pil +
  dDM + dDM : a.dDM +
    dDM : year.trt : ret +
    dDM : year.trt : pil +
    dDM : a.dDM : year.trt : ret +
    dDM : a.dDM : year.trt : pil +
  ed + ed : a.ed +
    ed : year.trt : ret +
    ed : year.trt : pil +
    ed : a.ed : year.trt : ret +
    ed : a.ed : year.trt : pil +
  
  # TOPO
  twi + twi2 + vrm + vrm2 +
  
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
  f(MRID11, dDM, model = "iid", hyper = hyper.list) +
  f(MRID12, ed, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3e. M5 - (STAND x AVAIL x TRT) + (LAND x AVAIL) + TOPO ----
# ______________________________________________________________________________

M5.form <- case ~ 
  
  # population-level effects
  g.s +
  
  # STAND
  vo + vo : a.vo +
    vo : year.trt : ret +
    vo : year.trt : pil +
    vo : a.vo : year.trt : ret + 
    vo : a.vo : year.trt : pil +
  ch + ch : a.ch +
    ch : year.trt : ret +
    ch : year.trt : pil +
    ch : a.ch : year.trt : ret + 
    ch : a.ch : year.trt : pil + 
  ch2 + ch2 : a.ch +
    ch2 : year.trt : ret +
    ch2 : year.trt : pil +
    ch2 : a.ch : year.trt : ret + 
    ch2 : a.ch : year.trt : pil + 
  cc + cc : a.cc +
    cc : year.trt : ret +
    cc : year.trt : pil +
    cc : a.cc : year.trt : ret + 
    cc : a.cc : year.trt : pil + 
  cc2 + cc2 : a.cc +
    cc2 : year.trt : ret +
    cc2 : year.trt : pil +
    cc2 : a.cc : year.trt : ret + 
    cc2 : a.cc : year.trt : pil + 
  
  # LAND
  dOpen + dOpen : a.dOpen + 
  dDM + dDM : a.dDM +
  ed + ed : a.ed +
  
  # TOPO
  twi + twi2 + vrm + vrm2 +
  
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
  f(MRID11, dDM, model = "iid", hyper = hyper.list) +
  f(MRID12, ed, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 3f. M6 - (STAND x AVAIL x TRT) + (LAND x AVAIL x TRT) + TOPO ----
# ______________________________________________________________________________

M6.form <- case ~ 
  
  # population-level effects
  g.s +
  
  # STAND
  vo + vo : a.vo +
    vo : year.trt : ret +
    vo : year.trt : pil +
    vo : a.vo : year.trt : ret + 
    vo : a.vo : year.trt : pil +
  ch + ch : a.ch +
    ch : year.trt : ret +
    ch : year.trt : pil +
    ch : a.ch : year.trt : ret + 
    ch : a.ch : year.trt : pil + 
  ch2 + ch2 : a.ch +
    ch2 : year.trt : ret +
    ch2 : year.trt : pil +
    ch2 : a.ch : year.trt : ret + 
    ch2 : a.ch : year.trt : pil + 
  cc + cc : a.cc +
    cc : year.trt : ret +
    cc : year.trt : pil +
    cc : a.cc : year.trt : ret + 
    cc : a.cc : year.trt : pil + 
  cc2 + cc2 : a.cc +
    cc2 : year.trt : ret +
    cc2 : year.trt : pil +
    cc2 : a.cc : year.trt : ret + 
    cc2 : a.cc : year.trt : pil + 
  
  # LAND
  dOpen + dOpen : a.dOpen + 
    dOpen : year.trt : ret +
    dOpen : year.trt : pil +
    dOpen : a.dOpen : year.trt : ret +
    dOpen : a.dOpen : year.trt : pil +
  dDM + dDM : a.dDM +
    dDM : year.trt : ret +
    dDM : year.trt : pil +
    dDM : a.dDM : year.trt : ret +
    dDM : a.dDM : year.trt : pil +
  ed + ed : a.ed +
    ed : year.trt : ret +
    ed : year.trt : pil +
    ed : a.ed : year.trt : ret +
    ed : a.ed : year.trt : pil +
  
  # TOPO
  twi + twi2 + vrm + vrm2 +
  
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
  f(MRID11, dDM, model = "iid", hyper = hyper.list) +
  f(MRID12, ed, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 4. Fit models ----

compute.list <- list(cpo = T,      # CPO/PIT (similar to PPPvals)
                     dic = T)      # more appropriate than WAIC

# ______________________________________________________________________________

M1.fit <- inla(M1.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = compute.list)  

M2.fit <- inla(M2.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = compute.list)  

M3.fit <- inla(M3.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = compute.list)  

M4.fit <- inla(M4.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = compute.list)  

M5.fit <- inla(M5.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = compute.list)  

M6.fit <- inla(M6.form,
               weights = data.off$w,
               family = "binomial",
               data = data.off,
               control.compute = compute.list)  

# ______________________________________________________________________________
# 5. Model comparison ----

# function
inla_cpo_table <- function (x) {
  
  # x is a list of models
  comp.table <- data.frame()
  
  for (i in 1:length(x)) {
    
    model.focal <- x[[i]]
    
    comp.table.focal <- data.frame(
      
      Model = case_when(i == 1 ~ "Land FR",
                        i == 2 ~ "Stand FR",
                        i == 3 ~ "Stand TRT",
                        i == 4 ~ "Land FR x TRT",
                        i == 5 ~ "Stand FR x TRT",
                        i == 6 ~ "Overall FR x TRT"),
      k.fixed = nrow(model.focal$summary.fixed),
      CPO = -sum(log(model.focal$cpo$cpo)),
      mLL = model.focal$mlik[2]
      
    )
    
    # bind in 
    comp.table <- rbind(comp.table, comp.table.focal)
    
  }
  
  # compute dDIC and arrange
  comp.table <- comp.table |>
    
    mutate(dCPO = max(CPO) - CPO) |>
    
    # arrange
    arrange(dCPO) |>
    
    dplyr::select(Model, k.fixed, CPO, dCPO, mLL)
  
  return(comp.table)
  
}

# ______________________________________________________________________________

inla_cpo_table(list(M1.fit,
                     M2.fit,
                     M3.fit,
                     M4.fit,
                     M5.fit,
                     M6.fit))

# weird issues with effective parameters for M4 (DIC problem for IPP models?)
# M4 performs the best here

summary(M4.fit)

# ______________________________________________________________________________
# 6. Save to file ----
# ______________________________________________________________________________

saveRDS(M4.fit, "model_tests/off_M4.rds")
