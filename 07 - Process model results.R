# PROJECT: Habitat selection
# SCRIPT: 07 - Process model results
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 01 Jun 2026
# LAST MODIFIED: 01 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(INLA)

# ______________________________________________________________________________
# 2. Load models ----
# ______________________________________________________________________________

# off
off.M6 <- readRDS("model_tests/off_M6.rds")
off.M4 <- readRDS("model_tests/off_M4.rds")

on.M6 <- readRDS("model_tests/on_M6.rds")
on.M5 <- readRDS("model_tests/on_M5.rds")

# ______________________________________________________________________________
# 3. Population-level effects ----

# function
extract_pop_effects <- function (.model,
                                 .ci = 0.9) {
  
  # population-level marginals
  marg.pop <- .model$marginals.fixed
  
  # effects (drop intercepts)
  pop.effects <- cbind(
      
      # posterior means 
      do.call(rbind, lapply(marg.pop, FUN = inla.emarginal, fun = mean))[-1, ],
      
      # posterior sd
      .model$summary.fixed[-1 , 2],
  
      # posterior HPD intervals
      do.call(rbind, lapply(marg.pop, FUN = inla.hpdmarginal, p = .ci))[-1, ]
     
    ) |>
    
    as.data.frame()
  
  # variable names
  pop.effects <- cbind(names(.model$marginals.fixed)[-1], pop.effects)
  
  # names
  colnames(pop.effects) <- c("param", "mean", "sd", "low", "upp")
  rownames(pop.effects) <- NULL
  
  return(pop.effects)
  
}

# ______________________________________________________________________________

off.M6.pop.eff <- extract_pop_effects(off.M6, 0.9)
off.M4.pop.eff <- extract_pop_effects(off.M4, 0.9)

on.M6.pop.eff <- extract_pop_effects(on.M6, 0.9)
on.M5.pop.eff <- extract_pop_effects(on.M5, 0.9)

# ______________________________________________________________________________
# 4. Random effect variances ----

# function
extract_RS_variance <- function (.model) {
  
  # hyperparameter marginals
  marg.hyp <- .model$marginals.hyperpar
  
  # random slope variances (from Muff et al)
  hyp.var <- sapply(marg.hyp,
                     
                     function(y)
                       
                       # must transform from precision
                       inla.emarginal(function (y) y,
                                      inla.tmarginal(function (y) 1 / y, y))
                     
  ) |>
    
    unname()
  
  # add names
  # ignore intercept and g.s
  hyp.var.df <- data.frame(param = names(.model$marginals.fixed)[3:(length(hyp.var) + 2)],
                           variance = hyp.var,
                           sd = sqrt(hyp.var))
  
  return(hyp.var.df)
  
}

# ______________________________________________________________________________

off.M6.rs.var <- extract_RS_variance(off.M6)
off.M4.rs.var <- extract_RS_variance(off.M4)

on.M6.rs.var <- extract_RS_variance(on.M6)
on.M5.rs.var <- extract_RS_variance(on.M5)

# ______________________________________________________________________________
# 5. Individual-level slopes ----

# function
extract_RS <- function (.model) {
  
  # random slope adjustment marginals (remove intercept)
  marg.rs <- .model$marginals.random[-1]
  
  # param names (remove intercept and g.s)
  param.names <- rownames(.model$summary.fixed)[3:(length(marg.rs) + 2)]
  
  # internal function - to apply across list elements (i.e., parameters)
  extract_RS_int <- function (.param) {
    
    focal.param.deviation <- do.call(rbind,
                                     lapply(.param,
                                            inla.emarginal,
                                            fun = mean))
    
    return(focal.param.deviation)
    
  }
  
  all.deviations <- do.call(cbind, lapply(marg.rs, extract_RS_int))
  
  # add column-wise
  # population-level effects
  pop.means <- do.call(cbind, lapply(.model$marginals.fixed, 
                                     inla.emarginal, 
                                     fun = mean))[, -c(1, 2)]
  
  cond.slopes <- pop.means + all.deviations |>
    
    as.data.frame()
  
  # parameter names
  colnames(cond.slopes) <- param.names
  
  # add MRIDs (double check that these are in the right order!!)
  cond.slopes <- cbind(MRID = .model$summary.random$MRID$ID,
                       cond.slopes)
  
  # remove rownmaes for cleanness
  rownames(cond.slopes) <- NULL
  
  # return
  return(cond.slopes)
  
}

# ______________________________________________________________________________

off.M6.rs <- extract_RS(off.M6)
off.M4.rs <- extract_RS(off.M4)

on.M6.rs <- extract_RS(on.M6)
on.M5.rs <- extract_RS(on.M5)

# ______________________________________________________________________________
# 6. Write to file ----
# ______________________________________________________________________________
# 6a. Lists ----
# ______________________________________________________________________________

off.M6.list <- list(off.M6.pop.eff,
                    off.M6.rs.var,
                    off.M6.rs)

off.M4.list <- list(off.M4.pop.eff,
                    off.M4.rs.var,
                    off.M4.rs)

on.M6.list <- list(on.M6.pop.eff,
                   on.M6.rs.var,
                   on.M6.rs)

on.M5.list <- list(on.M5.pop.eff,
                   on.M5.rs.var,
                   on.M5.rs)

# ______________________________________________________________________________
# 6b. Write ----
# ______________________________________________________________________________

saveRDS(off.M6.list, "model_results/off_M6.rds")
saveRDS(off.M4.list, "model_results/off_M4.rds")

saveRDS(on.M6.list, "model_results/on_M6.rds")
saveRDS(on.M5.list, "model_results/on_M5.rds")
