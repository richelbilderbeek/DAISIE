context("DAISIE_sim")

test_that("A clean classic run should produce no output", {
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

test_that("A clean ontogeny run should produce no output", {
  n_mainland_species <- 1000
  island_age <- 4
  clado_rate <- 2.550687345 # cladogenesis rate
  ext_rate <- 2.683454548 # extinction rate
  clade_carr_cap <- 10.0  # clade-level carrying capacity
  imm_rate <- 0.00933207 # immigration rate
  ana_rate <- 1.010073119 # anagenesis rate
  max_area <- 1000
  peak_time <- 0.5
  sharpness <- 1
  total_island_age <- 5
  mu_min <- 0.5
  mu_max <- 100
  island_ontogeny <- "quadratic"
  
  expect_silent(
    DAISIE_sim(
      time = island_age, 
      M = n_mainland_species, 
      pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
      replicates = 1, 
      island_ontogeny = island_ontogeny,
      Apars = c(max_area, peak_time, sharpness, total_island_age),
      Epars = c(mu_min, mu_max),
      plot_sims = FALSE,
      verbose = FALSE
    )
  )