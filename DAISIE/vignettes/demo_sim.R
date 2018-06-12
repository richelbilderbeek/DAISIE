## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
 collapse = TRUE,
 comment = "#>"
)

## ------------------------------------------------------------------------
library(DAISIE)

## ------------------------------------------------------------------------
pars <- c(2.550687345,2.683454548,Inf,0.00933207,1.010073119)

island_replicates <- DAISIE_sim( 
  time = 4, 
  M = 1000, 
  pars = pars, 
  replicates = 100,
  verbose = FALSE
)

