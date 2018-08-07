context("update_rates")

test_that("update rates use", {

  #testit::assert(is.matrix(c()))
  # Does not give errors. One day, it can be checked to be silent  
  set.seed(42)
  update_rates(
    timeval = 0, 
    totaltime = 1, 
    gam = 0.009, 
    mu = 2.0, 
    laa = 1.0, 
    lac = 2.5, 
    Apars = c(1.0, 0.5, 1.0, 1.0), 
    Epars = c(0.5, 10.0),
    island_ontogeny = "quadratic", 
    extcutoff = 1000.0, 
    K = 3, 
    island_spec = c(), 
    mainland_n = 1, 
    thor_ext = 0.5, 
    thor_c_i = 0.25
  )

})
