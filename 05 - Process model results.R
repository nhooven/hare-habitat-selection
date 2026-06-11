# PROJECT: Habitat selection
# SCRIPT: 05 - Process model results
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 01 Jun 2026
# LAST MODIFIED: 11 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(INLA)

# ______________________________________________________________________________
# 2. Load models ----
# ______________________________________________________________________________

M.off <- readRDS("model_tests/off.rds")
M.on <- readRDS("model_tests/on.rds")

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

M.off.pop.eff <- extract_pop_effects(M.off, 0.9)
M.on.pop.eff <- extract_pop_effects(M.on, 0.9)

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
  # ignore intercept
  hyp.var.df <- data.frame(param = names(.model$marginals.fixed)[2:(length(hyp.var) + 1)],
                           variance = hyp.var,
                           sd = sqrt(hyp.var))
  
  return(hyp.var.df)
  
}

# ______________________________________________________________________________

M.off.rs.var <- extract_RS_variance(M.off)
M.on.rs.var <- extract_RS_variance(M.on)

# ______________________________________________________________________________
# 5. Individual-level slopes ----

# function
extract_RS <- function (.model) {
  
  # random slope adjustment marginals (remove intercept and g.s)
  marg.rs <- .model$marginals.random[-c(1, 2)]
  
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
  
  cond.slopes <- sweep(all.deviations, 2, pop.means, FUN = "+") |> as.data.frame()
  
  # extract SDs
  # importantly, these are the SDs for the DEVIATIONS, not the conditional estimates
  # however we can still use them for IVW since they're all on the same scale
  cond.slopes.sd <- cbind(.model$summary.random$TSPID2$sd,
                          .model$summary.random$TSPID3$sd,
                          .model$summary.random$TSPID4$sd,
                          .model$summary.random$TSPID5$sd,
                          .model$summary.random$TSPID6$sd,
                          .model$summary.random$TSPID7$sd,
                          .model$summary.random$TSPID8$sd,
                          .model$summary.random$TSPID9$sd,
                          .model$summary.random$TSPID10$sd,
                          .model$summary.random$TSPID11$sd,
                          .model$summary.random$TSPID12$sd,
                          .model$summary.random$TSPID13$sd) |> as.data.frame()
  
  # parameter names
  colnames(cond.slopes) <- param.names
  colnames(cond.slopes.sd) <- param.names
  
  # add MRIDs (double check that these are in the right order!!)
  cond.slopes <- cbind(TSPID = .model$summary.random$TSPID$ID,
                       cond.slopes)
  
  cond.slopes.sd <- cbind(TSPID = .model$summary.random$TSPID$ID,
                          cond.slopes.sd)
  
  # remove rownames for cleanness
  rownames(cond.slopes) <- NULL
  rownames(cond.slopes.sd) <- NULL
  
  # return
  return(list(cond.slopes, cond.slopes.sd))
  
}

# ______________________________________________________________________________

M.off.rs <- extract_RS(M.off)
M.on.rs <- extract_RS(M.on)

# ______________________________________________________________________________
# 6. Write to file ----
# ______________________________________________________________________________
# 6a. Lists ----
# ______________________________________________________________________________

M.off.list <- list(M.off.pop.eff,
                   M.off.rs.var,
                   M.off.rs[[1]],
                   M.off.rs[[2]])

M.on.list <- list(M.on.pop.eff,
                   M.on.rs.var,
                   M.on.rs[[1]],
                   M.on.rs[[2]])

# ______________________________________________________________________________
# 6b. Write ----
# ______________________________________________________________________________

saveRDS(M.off.list, "model_results/M_off.rds")
saveRDS(M.on.list, "model_results/M_on.rds")
