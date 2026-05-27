# PROJECT: Habitat selection
# SCRIPT: 07 - Process model results
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

# ______________________________________________________________________________
# 2. Load models ----
# ______________________________________________________________________________

M1.off <- readRDS("model_tests/M1_fit.rds")

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

M1.off.pop.effects <- extract_pop_effects(M1.off, 0.9)

# ______________________________________________________________________________
# 4. Random effect variances ----

# function
extract_RS_variance <- function (.model) {
  
  # hyperparameter marginals
  marg.hyp <- .model$marginals.hyperpar
  
  # random slope variances (from Muff et al)
  hyp.var <- sapply(marg.hyp,
                     
                     function(y)
                       
                       inla.emarginal(function (y) y,
                                      inla.tmarginal(function (y) 1 / y, y))
                     
  ) |>
    
    unname()
  
  # add names
  hyp.var.df <- data.frame(param = names(.model$marginals.fixed)[-1],
                           variance = hyp.var,
                           sd = sqrt(hyp.var))
  
  return(hyp.var.df)
  
}

# ______________________________________________________________________________

M1.off.rs.var <- extract_RS_variance(M1.off)

# ______________________________________________________________________________
# 5. Individual-level slopes ----
# ______________________________________________________________________________

# TBD

# ______________________________________________________________________________
# 6. Write to file ----
# ______________________________________________________________________________
# 6a. Lists ----
# ______________________________________________________________________________

M1.off.list <- list(M1.off.pop.effects,
                    M1.off.rs.var)

# ______________________________________________________________________________
# 6b. Write ----
# ______________________________________________________________________________

saveRDS(M1.off.list, "model_results/M1_off.rds")
