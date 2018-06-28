#' Function to describe changes in area through time. Adapted from
#' Valente et al 2014 ProcB
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param Apars a numeric vector:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @family rates calculation
island_area <- function(timeval, totaltime, Apars, island_ontogeny){

  Tmax <- Apars[4] # total time used to be Apars[4] in ProcB paper
  Amax <- Apars[1] # maximum area
  Topt <- Apars[2] # peak position in %
  peak <- Apars[3] # peakiness - we specify a value of 1 but this is flexible.
  proptime<- timeval/Tmax	
  # Constant
  if (is.null(island_ontogeny)){
    return(Apars[1])
  }	
  if(island_ontogeny == "quadratic"){

    f <- Topt / (1 - Topt)
    a <- f * peak/ ( 1 + f)
    b <- peak / (1 + f) 
    At <- Amax * proptime^a * (1 - proptime)^ b/ ((a / (a + b))^a * (b / (a + b))^b)
    return(At)}
  
  #Linear decline
  if(island_ontogeny == "linear"){
    b <- Amax # intercept (peak area)
    m <- -(b / Topt) # slope
    At <- m * timeval + b
    return(At)
  }
}


#' Function to describe changes in extinction rate through time. From
#' Valente et al 2014 ProcB
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param mu per capita extinction rate in no ontogeny model
#' @param Apars a numeric vector:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars a numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. 
#' Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param extcutoff cutoff for extinction rate preventing it from being too 
#' large and slowing down simulation
#' @param island_spec matrix containing state of system
#' @param K carrying capacity
#' @seealso Does the same as \link{DAISIE_calc_clade_ext_rate}
#' @family rates calculation
#' @author Pedro Neves
get_ext_rate <- function(timeval, totaltime, mu,
                         Apars, Epars, 
                         island_ontogeny, 
                         extcutoff, island_spec,
                         K){
  # Epars[1] and Epars[2] (mu_min, mu_p) must be user specified
  if (is.null(island_ontogeny)){
    extrate <- mu * length(island_spec[,1])
    return(extrate)
    
    } else {
      
      
    X <- log(Epars[1] / Epars[2]) / log(0.1)
    extrate <- Epars[1]/((island_area(timeval, totaltime, Apars, island_ontogeny) / Apars[1])^X)
    extrate[which(extrate > extcutoff)] <- extcutoff
    extrate[which(extrate > extcutoff)] <- extcutoff
    extrate <- extrate * length(island_spec[,1])
    extrate
  }
}

#' Calculate anagenesis rate
#' @description Internal function. 
#' Calculates the anagenesis rate given the current number of
#' immigrant species and the per capita rate.
#' @param laa per capita anagenesis rate
#' @param island_spec matrix with current state of system
#' @seealso Does the same as \link{DAISIE_calc_clade_ana_rate}
#' @family rates calculation
#' @author Pedro Neves
get_ana_rate <- function(laa, island_spec) {
  ana_rate = laa * length(which(island_spec[,4] == "I"))
  ana_rate
} 

#' Calculate cladogenesis rate
#' @description Internal function. 
#' Calculates the cladogenesis rate given the current number of
#' species in the system, the carrying capacity and the per capita cladogenesis
#' rate
#' @param lac per capita cladogenesis rate
#' @param island_ontogeny a string describing the type of island ontogeny. 
#' Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param island_spec matrix with current state of system
#' @param K carrying capacity
#' @seealso Does the same as \link{DAISIE_calc_clade_clado_rate}
#' @author Pedro Neves
get_clado_rate <- function(lac, island_ontogeny, island_spec, K) {
  clado_rate = max(c(length(island_spec[,1])
                     * (lac * (1 - length(island_spec[, 1]) / K)),
                     0),
                   na.rm = T)
  clado_rate
}

#' Calculate immigration rate
#' @description Internal function. 
#' Calculates the immigration rate given the current number of
#' species in the system, the carrying capacity
#' @param gam per capita immigration rate
#' @param island_spec matrix with current state of system
#' @param K carrying capacity
#' @param mainland_n total number of species present in the mainland
#' @seealso Does the same as \link{DAISIE_calc_clade_imm_rate}
#' @family rates calculation
#' @author Pedro Neves
get_immig_rate <- function(gam,
                           island_spec,
                           K, mainland_n) {
  immig_rate = max(c(mainland_n 
                     * gam * (1 - length(island_spec[,1]) / K), 0), na.rm = T)
  immig_rate
}

#' Function to calculate and update horizon for maximum extinction rate
#' @description Internal function. 
#' Calculates when the next horizon for maximum extinction will be in the 
#' simulation
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param Apars a numeric vector:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param ext_multiplier reduces or increases distance of horizon to current
#' simulation time
#' @param island_ontogeny a string describing the type of island ontogeny.
#'  Can be \code{NULL}, \code{"quadratic"} for a beta function
#'   describing area through time, or \code{"linear"} for a linear function
#' @param thor time of horizon for max extinction
#' @family rates calculation
#' @author Pedro Neves
get_thor <- function(timeval,
                     totaltime,
                     Apars,
                     ext_multiplier,
                     island_ontogeny,
                     thor) {
  # Function calculates where the horizon for max(ext_rate) is.
  if (is.null(island_ontogeny)) {
    thor <- totaltime
    return(thor)
  } else {
    
    if (is.null(thor)){
      thor <- Apars[2] * Apars[4]
      return(thor)
      
    } else if (timeval >= thor) {
      
      thor <- timeval + ext_multiplier * (totaltime - timeval)
      thor <- min(totaltime, thor)
      thor
    } 
  }
}

