context("DAISIE_sim_core_1_4")

test_that("should give same results", {
  
  sim_time <- 10
  n_mainland_species <- 1 
  clado_rate <- 1.0
  ext_rate <- 0.5
  carr_cap <- 10 
  imm_rate <- 1.0
  ana_rate <- 1.0
  rng_seed <- 42
  set.seed(rng_seed)
  new <- DAISIE:::DAISIE_sim_core_checked(
    sim_time = sim_time, 
    n_mainland_species = n_mainland_species, 
    clado_rate = clado_rate, 
    ext_rate = ext_rate,
    carr_cap = carr_cap,
    imm_rate = imm_rate,
    ana_rate = ana_rate,
    island_ontogeny = NULL
  )
  set.seed(rng_seed)
  old <- DAISIE:::DAISIE_sim_core_checked_1_4(
    sim_time = sim_time, 
    n_mainland_species = n_mainland_species, 
    clado_rate = clado_rate, 
    ext_rate = ext_rate,
    carr_cap = carr_cap,
    imm_rate = imm_rate,
    ana_rate = ana_rate
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

  # skip("For Pedro")
  expect_true(new$stt_table == old$stt_table)
  expect_true(new$branching_times == old$branching_times)
  expect_true(new$other_clades_same_ancestor[[1]]$brts_miss == old$other_clades_same_ancestor[[1]]$brts_miss)
})