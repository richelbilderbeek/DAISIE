context("ddd_sample2")

test_that("sample 1:4 using prob 0 in 4 is not the same as sampling 1:3", {
  
  set.seed(42)
  event_1 <- DDD::sample2(1:4, 1e4, replace = T, prob = c(1, 1, 1, 0))
  set.seed(42)
  event_2 <- DDD::sample2(1:3, 1e4, replace = T, prob = c(1, 1, 1))
  skip("For Richel")
  expect_equal(event_1, event_2)
  
})
