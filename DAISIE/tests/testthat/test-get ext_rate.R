context("get ext_rate")

test_that("ext rate is a number", {

  
  get_ext_rate(
    timeval = 0,
    totaltime = 1,
    mu = 2,
    Apars = c(10, 0.5, 1, 15),
    Epars = c(1, 10), 
    island_ontogeny = "quadratic", 
    extcutoff = 1000,
    island_spec = c(),
    K = 10
  )
  
  expect_silent(
    is.numeric(
      get_ext_rate(
        timeval = 0,
        totaltime = 1,
        mu = 2,
        Apars = c(10, 0.5, 1, 15),
        Epars = c(1, 10), 
        island_ontogeny = "quadratic", 
        extcutoff = 1000,
        island_spec = c(),
        K = 10
      )
    )
  )
})
