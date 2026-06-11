# PROJECT: Habitat selection
# SCRIPT: 06a - HSF modeling (snow-off)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 
# LAST MODIFIED: 11 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(INLA)

# ______________________________________________________________________________
# 2. Read in data ----
# ______________________________________________________________________________

hs.data <- readRDS("data_for_model/off_data.rds")

# residuals for log(AKDE)
hs.data$g.s <- residuals(lm(log(akde) ~ vo + ch + cc + cc2 + 
                               twi + twi2 + vrm + vrm2 + 
                               dOM + dDM + ed + shdi,
                             data = hs.data))

# ______________________________________________________________________________
# 3. Correlation and VIF ----

# VIF function
# input is a df of covariates
calc_vif <- function (x) {
  
  all.vif <- data.frame()
  
  for (j in 1:ncol(x)) {
    
    x.1 <- cbind(x[ ,j], x[ ,-j])
    
    # rename response
    names(x.1)[1] <- "response"
    
    x.lm <- lm(response ~ ., data = x.1)
    
    # extract r2
    r2 <- summary(x.lm)$r.squared
    
    # calculate VIF
    vif.1 <- 1 / (1 - r2)
    
    # bind in
    all.vif <- rbind(all.vif,
                     data.frame(cov = names(x)[j],
                                vif = vif.1))
    
  }
  
  return(all.vif)
  
}

# ______________________________________________________________________________

# subset just the linear coefficients
# assume squared terms will be highly correlated to their linear terms
covs.lin <- hs.data |> dplyr::select(vo, ch, cc, twi, vrm, dOM, dDM, ed, shdi, g.s)

# correlation
cor(covs.lin, method = "spearman") |> round(2)   # nothing over 0.65

# VIF
calc_vif(covs.lin)  # all < 2.5

# ______________________________________________________________________________
# 3. Setup ----
# ______________________________________________________________________________
# 3a. Copy TSP for random slopes ----
# ______________________________________________________________________________

hs.data <- hs.data |>
  
  rename(TSPID = track_season_post) |>
  
  mutate(
    
    TSPID = TSPID,
    TSPID1 = TSPID,
    TSPID2 = TSPID,
    TSPID3 = TSPID,
    TSPID4 = TSPID,
    TSPID5 = TSPID,
    TSPID6 = TSPID,
    TSPID7 = TSPID,
    TSPID8 = TSPID,
    TSPID9 = TSPID,
    TSPID10 = TSPID,
    TSPID11 = TSPID,
    TSPID12 = TSPID,
    TSPID13 = TSPID
    
  )

# ______________________________________________________________________________
# 3b. Hyperparameter prior list ----

hyper.list <- list(theta = list(initial = log(1),
                                fixed = F,
                                prior = "pc.prec",
                                param = c(3, 0.05)))

# ______________________________________________________________________________
# 3c. Fixed effect prior list ----

# fixed effect priors (Eisaguirre et al. 2025)
fixed.list <- list(mean = 0,
                   prec = 1)

# ______________________________________________________________________________
# 3d. Compute control parameters ----

compute.list <- list(config = T)

# ______________________________________________________________________________
# 4. Model formula ----
# ______________________________________________________________________________

M.form <- case ~ 
  
  # population-level effects
  g.s +
  
  # STAND
  vo +
  ch + 
  cc + cc2 +
  
  # TOPO
  twi + twi2 + vrm + vrm2 +
  
  # LAND
  dOM +
  dDM + 
  ed +
  shdi +
  
  # random intercepts
  f(TSPID, model = "iid", hyper = list(theta = list(initial = log(1/1e6), fixed = T))) +
  
  # random slopes
  f(TSPID1, g.s, model = "iid", hyper = hyper.list) +
  f(TSPID2, vo, model = "iid", hyper = hyper.list) +
  f(TSPID3, ch, model = "iid", hyper = hyper.list) +
  f(TSPID4, cc, model = "iid", hyper = hyper.list) + 
  f(TSPID5, cc2, model = "iid", hyper = hyper.list) +
  f(TSPID6, twi, model = "iid", hyper = hyper.list) +
  f(TSPID7, twi2, model = "iid", hyper = hyper.list) +
  f(TSPID8, vrm, model = "iid", hyper = hyper.list) +
  f(TSPID9, vrm2, model = "iid", hyper = hyper.list) +
  f(TSPID10, dOM, model = "iid", hyper = hyper.list) +
  f(TSPID11, dDM, model = "iid", hyper = hyper.list) +
  f(TSPID12, ed, model = "iid", hyper = hyper.list) +
  f(TSPID13, shdi, model = "iid", hyper = hyper.list)

# ______________________________________________________________________________
# 5. Fit model ----
# ______________________________________________________________________________

M.fit <- inla(M.form,
              weights = hs.data$w,
              family = "binomial",
              data = hs.data,
              control.fixed = fixed.list,
              control.compute = compute.list) 

# ______________________________________________________________________________
# 6. Summary ----
# ______________________________________________________________________________

summary(M.fit)

# ______________________________________________________________________________
# 7. Save model ----
# ______________________________________________________________________________

saveRDS(M.fit, "model_tests/off.rds")
