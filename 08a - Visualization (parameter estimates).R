# PROJECT: Habitat selection
# SCRIPT: 08a - Visualization (parameter estimates)
# AUTHOR: Nate Hooven
# EMAIL: nathan.d.hooven@gmail.com
# BEGAN: 27 May 2026
# COMPLETED: 
# LAST MODIFIED: 03 Jun 2026
# R VERSION: 4.5.2

# ______________________________________________________________________________
# 1. Load packages ----
# ______________________________________________________________________________

library(tidyverse)
library(hoove)
library(cowplot)

# ______________________________________________________________________________
# 2. Read in cleaned model results ----
# ______________________________________________________________________________

off.M6 <- readRDS("model_results/off_M6.rds")
off.M4 <- readRDS("model_results/off_M4.rds")

on.M6 <- readRDS("model_results/on_M6.rds")
on.M5 <- readRDS("model_results/on_M5.rds")

# ______________________________________________________________________________
# 3. Re-order and rename factor levels ----
# ______________________________________________________________________________
# 3a. Define order and labels ----
# ______________________________________________________________________________

# levels 
levels.off.M6 <- off.M6[[1]]$param
levels.off.M4 <- off.M4[[1]]$param
levels.on.M6 <- on.M6[[1]]$param
levels.on.M5 <- on.M5[[1]]$param

# base levels
levels.base.off <- levels.off.M6[1:13]
levels.base.on <- levels.on.M6[1:13]

# labels
labels.base.off <- c("log(G[s])", 
                     "visual obstruction", "canopy height", "canopy height2",
                     "canopy cover", "canopy cover2",
                     "dOpen", "dDense Mature", "edge density",
                     "wetness", "wetness2", "ruggedness", "ruggedness2")

labels.base.on <- c("log(G[s])", 
                    "stem density", "canopy height", "canopy height2",
                    "canopy cover", "canopy cover2",
                    "dOpen", "dDense Mature", "edge density",
                    "wetness", "wetness2", "ruggedness", "ruggedness2")

labels.off.M6 <- c(labels.base.off,
                   "visual obstruction x A",
                   "canopy height x A", "canopy height2 x A",
                   "canopy cover x A", "canopy cover2 x A",
                   "dOpen x A", "dDense Mature x A", "edge density x A",
                   "visual obstruction x RET", "visual obstruction x PIL",
                   "canopy height x RET", "canopy height x PIL", 
                   "canopy height2 x RET", "canopy height2 x PIL",
                   "canopy cover x RET", "canopy cover x PIL", 
                   "canopy cover2 x RET", "canopy cover2 x PIL",
                   "dOpen x RET", "dOpen x PIL",
                   "dDense Mature x RET", "dDense Mature x PIL",
                   "edge density x RET", "edge density x PIL",
                   "visual obstruction x A x RET", "visual obstruction x A x PIL",
                   "canopy height x A x RET", "canopy height x A x PIL", 
                   "canopy height2 x A x RET", "canopy height2 x A x PIL",
                   "canopy cover x A x RET", "canopy cover x A x PIL", 
                   "canopy cover2 x A x RET", "canopy cover2 x A x PIL",
                   "dOpen x A x RET", "dOpen x A x PIL", 
                   "dDense Mature x A x RET", "dDense Mature x A x PIL", 
                   "edge density x A x RET", "edge density x A x PIL")

labels.off.M4 <- c(labels.base.off,
                   "dOpen x A", "dDense Mature x A", "edge density x A",
                   "dOpen x RET", "dOpen x PIL",
                   "dDense Mature x RET", "dDense Mature x PIL",
                   "edge density x RET", "edge density x PIL",
                   "dOpen x A x RET", "dOpen x A x PIL", 
                   "dDense Mature x A x RET", "dDense Mature x A x PIL", 
                   "edge density x A x RET", "edge density x A x PIL")

labels.on.M6 <- c(labels.base.on,
                  "stem density x A",
                  "canopy height x A", "canopy height2 x A",
                  "canopy cover x A", "canopy cover2 x A",
                  "log(dOpen) x A", "log(dDense Mature) x A", "edge density x A",
                  "stem density x A x RET", "stem density x A x PIL",
                  "canopy height x A x RET", "canopy height x A x PIL", 
                  "canopy height2 x A x RET", "canopy height2 x A x PIL",
                  "canopy cover x A x RET", "canopy cover x A x PIL", 
                  "canopy cover2 x A x RET", "canopy cover2 x A x PIL",
                  "log(dOpen) x A x RET", "log(dOpen) x A x PIL", 
                  "log(dDense Mature) x A x RET", "log(dDense Mature) x A x PIL", 
                  "edge density x A x RET", "edge density x PIL")

labels.on.M5 <- c(labels.base.on,
                  "stem density x A",
                  "canopy height x A", "canopy height2 x A",
                  "canopy cover x A", "canopy cover2 x A",
                  "log(dOpen) x A", "log(dDense Mature) x A", "edge density x A",
                  "stem density x A x RET", "stem density x A x PIL",
                  "canopy height x A x RET", "canopy height x A x PIL", 
                  "canopy height2 x A x RET", "canopy height2 x A x PIL",
                  "canopy cover x A x RET", "canopy cover x A x PIL", 
                  "canopy cover2 x A x RET", "canopy cover2 x A x PIL")

# double check
cbind(levels.off.M6, labels.off.M6)
cbind(levels.off.M4, labels.off.M4)
cbind(levels.on.M6, labels.on.M6)
cbind(levels.on.M5, labels.on.M5)

# ______________________________________________________________________________
# 3b. Function ----
# ______________________________________________________________________________

apply_cov_labels <- function (.results,
                              .levels,
                              .labels) {
  
  .results.1 <- .results |>
    
    mutate(param = factor(param,
                          levels = .levels,
                          labels = .labels))
  
  return(.results.1)
  
}

# use on population-level parameter estimates
off.M6.pop <- apply_cov_labels(off.M6[[1]], levels.off.M6, labels.off.M6)
off.M4.pop <- apply_cov_labels(off.M4[[1]], levels.off.M4, labels.off.M4)
on.M6.pop <- apply_cov_labels(on.M6[[1]], levels.on.M6, labels.on.M6)
on.M5.pop <- apply_cov_labels(on.M5[[1]], levels.on.M5, labels.on.M5)

# use on random slope variances
off.M6.rs.var <- apply_cov_labels(off.M6[[2]], levels.base.off, labels.base.off)
off.M4.rs.var <- apply_cov_labels(off.M4[[2]], levels.base.off, labels.base.off)
on.M6.rs.var <- apply_cov_labels(on.M6[[2]], levels.base.on, labels.base.on)
on.M5.rs.var <- apply_cov_labels(on.M5[[2]], levels.base.on, labels.base.on)

# add variable class?

# ______________________________________________________________________________
# 4. Population-level effect plot ----
# ______________________________________________________________________________

# function
plot_pop <- function (.pop) {
  
  ggplot(data = .pop) +
    
    theme_hoove() +
    
    # vline
    geom_vline(xintercept = 0,
               linetype = "dashed") + 
    
    # CIs
    geom_errorbar(aes(x = mean,
                      y = param,
                      xmin = low,
                      xmax = upp),
                  width = 0,
                  linewidth = 2,
                  color = "gray50") +
    
    # means
    geom_point(aes(x = mean,
                   y = param),
               shape = 21,
               color = "black",
               fill = "white",
               size = 2,
               stroke = 1.5) +
    
    theme(axis.title.y = element_blank()) +
    
    xlab("Standardized coefficient") +
    
    scale_y_discrete(limits = rev)
  
}

# use function
plot_pop(off.M6.pop)
plot_pop(off.M4.pop)
plot_pop(on.M6.pop)
plot_pop(on.M5.pop)

# ______________________________________________________________________________
# 5. Random slope variance plot ----
# ______________________________________________________________________________

plot_var <- function (.rs.var) {
  
  ggplot(data = .rs.var) +
    
    theme_hoove() +
    
    geom_col(aes(x = variance,
                 y = param),
             color = "black",
             fill = "#FF3300",
             alpha = 0.5) +
    
    theme(axis.title.y = element_blank()) +
    
    xlab("Random slope variance") +
    
    scale_y_discrete(limits = rev)
  
}

# use
plot_var(off.M6.rs.var)
plot_var(off.M4.rs.var)
plot_var(on.M6.rs.var)
plot_var(on.M5.rs.var)
