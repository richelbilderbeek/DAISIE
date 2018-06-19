context("DAISIE_sim_core_checked")

test_that("use", {

  set.seed(42)
  result <- NULL
  expect_silent(
    result <- DAISIE::DAISIE_sim_core_checked(
      sim_time = 4, 
      n_mainland_species = 1, 
      clado_rate = 1.0, 
      ext_rate = 0.1,
      carr_cap = 10,
      imm_rate = 1.0,
      ana_rate = 1.0
    )
  )
  result
  expect_true("stt_table" %in% names(result))
  expect_true("branching_times" %in% names(result))
  expect_true("stac" %in% names(result))
  expect_true("missing_species" %in% names(result))
  # May also return 'other_clades_same_ancestor'
})
