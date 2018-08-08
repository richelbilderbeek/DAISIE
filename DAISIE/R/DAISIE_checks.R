# Integrate test by Giovanni Laudanno
island_area_for_test <- function(timeval, totaltime, Apars, island_function_shape){
  testit::assert(are_area_params(Apars))
  Tmax <- Apars$total_island_age # total time A PARS 1
  Amax <- Apars$max_area # maximum area
  Topt <- Apars$proportional_peak_t # peak position in %
  peak <- Apars$peak_sharpness # peakiness - we specify a value of 1 but this is flexible.
  proptime <- timeval/Tmax
  # Constant
  if (is.null(island_function_shape)){
    return(Apars$max_area)
  }
  # Beta function
  if(island_function_shape == "quadratic") {
    
    f <- Topt / (1 - Topt)
    a <- f * peak/ ( 1 + f)
    b <- peak / (1 + f)
    At <- Amax * proptime^a * (1 - proptime)^ b/ ((a / (a + b))^a * (b
                                                                     / (a + b))^b)
    return(At)}
  
  #Linear decline
  if(island_function_shape == "linear") {
    b <- Amax # intercept (peak area)
    m <- -(b / Topt) # slope
    At <- m * timeval + b
    return(At)
  }
}


# Function to describe changes in extinction rate through time. From
# Valente et al 2014 ProcB
get_ext_rate_for_test <- function(timeval, totaltime, mu,
                         Apars, Epars,
                         island_function_shape,
                         extcutoff, N,
                         K){
  # Epars[1] and Epars[2] (mu_min, mu_p) must be user specified
  testit::assert(are_area_params(Apars))
  if (is.null(island_function_shape)){
    extrate <- mu * N
    
  } else {
    
    X <- log(Epars[1] / Epars[2]) / log(0.1)
    extrate <- Epars[1]/((island_area_for_test(timeval, totaltime, Apars,
                                      island_function_shape) / Apars$max_area)^X)
    extrate[which(extrate > extcutoff)] <- extcutoff
    extrate <- extrate * N
    return(extrate)
  }
}

MU   <- function(s, N = 1000){
  totaltime <- 10
  Apars <- create_area_params(1000, 0.2, 1, totaltime * 1.5)
  Epars <- c(1.7, 20)
  island_function_shape <- 'quadratic'
  extcutoff <- 1000
  # N <- 1000
  K <- Inf
  
  MU.out <- get_ext_rate_for_test(timeval = s,
                         totaltime = totaltime,
                         mu = 0,
                         Apars = Apars,
                         Epars = Epars,
                         island_function_shape = island_function_shape,
                         extcutoff = extcutoff,
                         N = N,
                         K = K)
  return(MU.out)
}
RHO0 <- function(t, t0 = 0, P0 = 1000){
  RHO.out <- integrate(f = MU, lower = t0, upper = t, N = P0)$value
  return(RHO.out)
}
RHO  <- function(t, t0 = 0, P0 = 1000){
  return(Vectorize(RHO0(t = t, t0 = t0, P0 = P0)))
}
Pt   <- function(t, t0 = 0, P0 = 1000){
  out <- P0 * exp(-RHO(t = t, P0 = P0, t0 = t0))
  return(out)
}
Pt(0, 10, 1000)


PPt <- tt <- seq(0.04,10,0.02);
for (i in 1:length(tt))
{
  PPt[i] <- Pt(t = tt[i])
}

