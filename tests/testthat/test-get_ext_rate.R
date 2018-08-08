context("get_ext_rate")

test_that("classic behaviour", {
  
  
  carr_cap <- 10
  ps_ext_rate <- 2
  n_species <- 5
  n_mainland_species <- 2
  expected <- DAISIE_calc_clade_ext_rate(
    ps_ext_rate = ps_ext_rate,
    n_species = n_species
  )
  created <- get_ext_rate(
    timeval = 1.0,
    totaltime = 10.0,
    mu = ps_ext_rate,
    Apars =  NULL,
    Epars = NULL,
    island_ontogeny = NULL,
    island_spec = matrix(data = NA, nrow = n_species, ncol = 1),
    K = carr_cap,
    extcutoff = 1000
  )
  expect_equal(expected, created)
})
