context("get_thor")

test_that("minimal use", {
  
  expect_silent(
    get_thor(
      timeval = 1,
      totaltime = 5,
      Apars = create_area_params(max_area = 10,
                                 proportional_peak_t = 0.5,
                                 peak_sharpness = 1,
                                 total_island_age = 5),
      ext_multiplier = 0.5,
      island_ontogeny = "quadratic",
      thor = NULL
    )
  )
})

test_that("classic behavior", {
  
  # That is, a simulation without island ontogeny
  total_time <- 12.34
  expected <- total_time
  created <- get_thor(
    timeval = 1,
    totaltime = total_time,
    Apars = NULL,
    ext_multiplier = 0.5,
    island_ontogeny = NULL,
    thor = NULL
  )
  expect_equal(created, expected)
  
})
