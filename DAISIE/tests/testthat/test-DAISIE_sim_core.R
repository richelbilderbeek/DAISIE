context("DAISIE_sim_core")

test_that("Clean run should be silent", {

  set.seed(42)
  n_mainland_species <- 1
  sim_time <- 10
  clado_rate <- 1.0 
  ext_rate <- 0.1
  carr_cap <- 4
  imm_rate <- 1.0
  ana_rate <- 1.0
  
  expect_silent(
    DAISIE::DAISIE_sim_core(
      time = sim_time,
      mainland_n = n_mainland_species,
      pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
    )
  )

})
