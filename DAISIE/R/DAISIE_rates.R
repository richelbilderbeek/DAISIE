# Function to describe changes in area through time. Adapted from
# Valente et al 2014 ProcB

island_area <- function(timeval, totaltime, Apars, island_function_shape){

  Tmax <- Apars[4] # total time A PARS 1
  Amax <- Apars[1] # maximum area
  Topt <- Apars[2] # peak position in %
  peak <- Apars[3] # peakiness - we specify a value of 1 but this is flexible.
  proptime<- timeval/Tmax	
  # Constant
  if (is.null(island_function_shape)){
    return(Apars[1])
  }	
  if(island_function_shape == "quadratic"){

    f <- Topt / (1 - Topt)
    a <- f * peak/ ( 1 + f)
    b <- peak / (1 + f) 
    At <- Amax * proptime^a * (1 - proptime)^ b/ ((a / (a + b))^a * (b / (a + b))^b)
    return(At)}
  
  #Linear decline
  if(island_function_shape == "linear"){
    b <- Amax # intercept (peak area)
    m <- -(b / Topt) # slope
    At <- m * timeval + b
    return(At)
  }
}


# Function to describe changes in extinction rate through time. From
# Valente et al 2014 ProcB
get_ext_rate <- function(timeval, totaltime, mu,
                         Apars, Epars, 
                         island_function_shape, 
                         extcutoff, island_spec,
                         K){
  # Epars[1] and Epars[2] (mu_min, mu_p) must be user specified
  if (is.null(island_function_shape)){
    extrate <- mu * length(island_spec[,1])
  
    } else {
      
      
    X <- log(Epars[1] / Epars[2]) / log(0.1)
    extrate <- Epars[1]/((island_area(timeval, totaltime, Apars, island_function_shape) / Apars[1])^X)
    extrate[which(extrate > extcutoff)] <- extcutoff
    extrate[which(extrate > extcutoff)] <- extcutoff
    extrate <- extrate * length(island_spec[,1])
    return(extrate)
  }
}

# Function to calculate anagenesis rate given number of immigrant species
get_ana_rate <- function(timeval, totaltime, laa,
                         Apars, Epars, island_function_shape,
                         extcutoff, island_spec, K) {
  ana_rate = laa * length(which(island_spec[,4] == "I"))
  return(ana_rate)
} 

# Function to calculate cladogenesis rate given number of island species
get_clado_rate <- function(timeval, totaltime,
                           lac, Apars, Epars,
                           island_function_shape, 
                           extcutoff, island_spec,
                           K) {
  clado_rate = max(c(length(island_spec[,1]) * (lac * (1 - length(island_spec[, 1]) / K)), 0), na.rm = T)
  return(clado_rate)
}

get_immig_rate <- function(timeval, totaltime,
                           gam, Apars, Epars,
                           island_function_shape, 
                           extcutoff, island_spec,
                           K, mainland_n) {
  immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K), 0), na.rm = T)
  return(immig_rate)
}

get_thor <- function(timeval, totaltime, Apars, ext_multiplier, island_ontogeny, thor) {

  if (is.null(island_ontogeny)) {
    thor <- totaltime
    
  } else {
    
    if (is.null(thor)){
      thor <- Apars[2] * Apars[4]
      
    } else if (timeval >= thor) {
      
      thor <- timeval + ext_multiplier * (totaltime - timeval)
      thor <- min(totaltime, thor)
    } 
  }

  return(thor)
}

