context("DAISIE_sim_core_checked")

test_that("use", {

  result <- NULL
  expect_silent(
    result <- DAISIE::DAISIE_sim_core_checked(
      sim_time = 1, 
      n_mainland_species = 1, 
      clado_rate = 0.1, 
      ext_rate = 0.1,
      carr_cap = 10,
      imm_rate = 0.1,
      ana_rate = 1.0
    )
  )
  expect_true("stt_table" %in% names(result))
  expect_true("branching_times" %in% names(result))
  expect_true("stac" %in% names(result))
  expect_true("missing_species" %in% names(result))
})
