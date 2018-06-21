context("DAISIE_ONEcolonist")

test_that("use", {

  skip("WIP richelbilderbeek")
  sim_time <- 10

  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  expect_equal(2 * 2, 4)
})
