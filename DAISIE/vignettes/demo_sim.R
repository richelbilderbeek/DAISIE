## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
 collapse = TRUE,
 comment = "#>"
)

## ------------------------------------------------------------------------
library(DAISIE)

## ------------------------------------------------------------------------
n_replicates <- 2

## ----fig.width=7, fig.height=7-------------------------------------------
pars <- c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)

island_replicates <- DAISIE_sim( 
  time = 4, 
  M = 1000, 
  pars = pars, 
  replicates = n_replicates,
  plot_sims = FALSE,
  verbose = FALSE
)
DAISIE_plot_sims(island_replicates, use_dev_new = FALSE)

## ------------------------------------------------------------------------
island_replicates[[1]]

## ----fig.width=7, fig.height=7-------------------------------------------
pars <- c(2.550687345, 2.683454548, 10,0.00933207, 1.010073119) 

island_replicates_K <- DAISIE_sim( 
  time = 4, 
  M = 1000, 
  pars = pars, 
  replicates = n_replicates,
  plot_sims = FALSE,
  verbose = FALSE
) 
DAISIE_plot_sims(island_replicates_K, use_dev_new = FALSE)

## ----fig.width=7, fig.height=7-------------------------------------------
pars_type1 <- c(0.38, 0.55, Inf, 0.004, 1.10) 
pars_type2 <- c(2.28, 0.55, Inf, 0.004, 1.10) 

island_replicates_2types <- DAISIE_sim( 
  time = 4, 
  M = 1000, 
  pars = c(pars_type1,pars_type2), 
  replicates = n_replicates, 
  prop_type2_pool = 0.163,
  plot_sims = FALSE,
  verbose = FALSE
)
DAISIE_plot_sims(island_replicates_2types, use_dev_new = FALSE)

