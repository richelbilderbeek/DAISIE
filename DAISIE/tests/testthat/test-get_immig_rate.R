context("get_immig_rate")

test_that("multiplication works", {
  immig <- c()
  timepoints <- seq(0, 10, by = 0.01)
  
  for (i in 1:1000) {
    immig[i] <- get_immig_rate(timepoints[i], totaltime = 10, gam = 0.001,
                               Apars = c(5000, 0.2, 1, 15), 
                               island_spec = matrix(ncol = 1,), 
                               island_ontogeny = "quadratic", 
                               mainland_n = 1000)                        
  }
  return()
  plot(immig)
  
})
