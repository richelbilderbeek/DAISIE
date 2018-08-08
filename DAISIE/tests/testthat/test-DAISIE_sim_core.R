context("DAISIE_sim_core")

test_that("new and v1.4 should give same results", {
  
  sim_time <- 10
  n_mainland_species <- 1 
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10 
  imm_rate <- 1.0
  ana_rate <- 1.0
  pars <- c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
  rng_seed <- 42
  set.seed(rng_seed)
  new <- DAISIE_sim_core(
    time = sim_time, 
    mainland_n = n_mainland_species, 
    pars = pars
  )
  set.seed(rng_seed)
  old <- DAISIE:::DAISIE_sim_core_1_4(
    time = sim_time, 
    mainland_n = n_mainland_species, 
    pars = pars
  )
  
  expect_true(all(names(new) == names(old)))
  # stt_table has different content
  expect_true(nrow(new$stt_table) == nrow(old$stt_table))
  # different branching times
  expect_equal(length(new$branching_times), length(old$branching_times))
  expect_true(new$stac == old$stac)
  expect_true(new$missing_species == old$missing_species)
  expect_true(length(new$other_clades_same_ancestor) == length(old$other_clades_same_ancestor))
  expect_true(new$other_clades_same_ancestor[[1]]$species_type == old$other_clades_same_ancestor[[1]]$species_type)

  # stt_table, branching_times, missing species are different
  expect_equal(new$stt_table, old$stt_table)
  expect_equal(new$branching_times, old$branching_times)
  expect_true(new$other_clades_same_ancestor[[1]]$brts_miss == old$other_clades_same_ancestor[[1]]$brts_miss)
})

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

test_that("Pedro's should run silent", {
  
  set.seed(234567890)
  DAISIE_sim_core(
    time = 10, 
    mainland_n = 1, 
    pars = c(2.5, 2.2, 10, 0.009, 1.01), 
    Apars = create_area_params(
                               max_area = 10,
                               proportional_peak_t = 0.5,
                               peak_sharpness = 1,
                               total_island_age = 5), 
    Epars = c(1.7, 100), 
    island_ontogeny = "quadratic"
  )
  
  expect_silent(
    DAISIE_sim_core(
      time = 10, 
      mainland_n = 1, 
      pars = c(2.5, 2.2, 10, 0.009, 1.01), 
      Apars = c(5000, 0.2, 1, 15), 
      Epars = c(1.7, 100), 
      island_ontogeny = "quadratic"
    )
  )
  
})

