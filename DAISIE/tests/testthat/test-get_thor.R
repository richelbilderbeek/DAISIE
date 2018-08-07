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
