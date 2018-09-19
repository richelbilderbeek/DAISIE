context("test-update_thor_timeval")

test_that("use", {
  thor_timeval <- update_thor_timeval(timeval = 7,
                      totaltime = 10, 
                      Apars = create_area_params(
                        max_area = 10,
                        proportional_peak_t = 0.5,
                        peak_sharpness = 1,
                        total_island_age = 5),
                      ext_multiplier = 0.5,
                      island_ontogeny = "quadratic",
                      thor = 6)
  
  expect_true(if ("thor" %in% names(thor_timeval)) TRUE)
  expect_true(if ("timeval" %in% names(thor_timeval)) TRUE)
  expect_equal(thor_timeval$thor, 8.5)
  expect_equal(thor_timeval$timeval, 6)
  
})
