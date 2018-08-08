context("DAISIE_sim")

test_that("A clean run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  expect_silent(
    DAISIE_sim( 
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1,
      plot_sims = FALSE,
      verbose = FALSE
    )
  )
})
