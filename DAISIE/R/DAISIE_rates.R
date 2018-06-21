#' Function to describe changes in area through time. Adapted from
#' Valente et al 2014 ProcB
#' @param timeval something
#' @param totaltime something
#' @param Apars something
#' @param island_function_shape something
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


#' Function to describe changes in extinction rate through time. From
#' Valente et al 2014 ProcB
#' @param timeval something
#' @param totaltime something
#' @param mu something
#' @param Apars something
#' @param Epars something
#' @param island_function_shape something
#' @param extcutoff something
#' @param island_spec something
#' @param K something
#' @seealso Does the same as \link{DAISIE_calc_clade_ext_rate}
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

#' Function to calculate anagenesis rate given number of immigrant species
#' @param timeval something
#' @param totaltime something
#' @param laa something
#' @param Apars something
#' @param Epars something
#' @param island_function_shape something
#' @param extcutoff something
#' @param island_spec something
#' @param K something
#' @seealso Does the same as \link{DAISIE_calc_clade_ana_rate}
get_ana_rate <- function(timeval, totaltime, laa,
                         Apars, Epars, island_function_shape,
                         extcutoff, island_spec, K) {
  ana_rate = laa * length(which(island_spec[,4] == "I"))
  return(ana_rate)
} 

#' Function to calculate cladogenesis rate given number of island species
#' @param timeval something
#' @param totaltime something
#' @param lac something
#' @param Apars something
#' @param Epars something
#' @param island_function_shape something
#' @param extcutoff something
#' @param island_spec something
#' @param K something
#' @seealso Does the same as \link{DAISIE_calc_clade_clado_rate}
get_clado_rate <- function(timeval, totaltime,
                           lac, Apars, Epars,
                           island_function_shape, 
                           extcutoff, island_spec,
                           K) {
  clado_rate = max(c(length(island_spec[,1]) * (lac * (1 - length(island_spec[, 1]) / K)), 0), na.rm = T)
  return(clado_rate)
}

#' Something
#' @param timeval something
#' @param totaltime something
#' @param gam something
#' @param Apars something
#' @param Epars something
#' @param island_function_shape something
#' @param extcutoff something
#' @param island_spec something
#' @param K something
#' @param mainland_n something
#' @seealso Does the same as \link{DAISIE_calc_clade_imm_rate}
get_immig_rate <- function(timeval, totaltime,
                           gam, Apars, Epars,
                           island_function_shape, 
                           extcutoff, island_spec,
                           K, mainland_n) {
  immig_rate = max(c(mainland_n * gam * (1 - length(island_spec[,1])/K), 0), na.rm = T)
  return(immig_rate)
}

#' Something
#' @param timeval something
#' @param totaltime something
#' @param Apars something
#' @param ext_multiplier something
#' @param island_ontogeny something
#' @param thor something
get_thor <- function(timeval, totaltime, Apars, ext_multiplier, island_ontogeny, thor) {
  # Function calculates where the horizon for max(ext_rate) is.
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

#' Calculate the clade-wide extinction rate
#' @param ps_ext_rate per species extinction rate
#' @param n_species number of species in that clade
#' @return the clade's extinction rate
#' @author Richel J.C. Bilderbeek
#' @examples 
#'   testit::assert(
#'     DAISIE_calc_clade_ext_rate(
#'       ps_ext_rate = 0.2, 
#'       n_species = 4
#'     ) == 0.8
#'   )
#' @export
DAISIE_calc_clade_ext_rate <- function(ps_ext_rate, n_species) {
  testit::assert(ps_ext_rate >= 0.0)
  testit::assert(n_species >= 0)
  ps_ext_rate * n_species
}

#' Calculate the clade-wide effective anagenesis rate.
#' With 'effective', this means that if an immigrant
#' undergoes anagenesis, it will become a new species.
#' Would such a species undergo anagenesis again, no net new
#' species is created; the species only gets renamed
#' @param ps_ana_rate per species anagensis rate
#' @param n_immigrants number of immigrants in that clade
#' @return the clade's effective anagenesis rate
#' @author Richel J.C. Bilderbeek
#' @examples 
#'   testit::assert(
#'     DAISIE_calc_clade_ana_rate(
#'       ps_ana_rate = 0.3,  
#'       n_immigrants = 5
#'     ) == 1.5
#'   )
#' @export
DAISIE_calc_clade_ana_rate <- function(ps_ana_rate, n_immigrants) {
  testit::assert(ps_ana_rate >= 0.0)
  testit::assert(n_immigrants >= 0)
  ps_ana_rate * n_immigrants
}

#' Calculate the clade-wide cladogenesis rate.
#' @param ps_clado_rate per species cladogenesis rate
#' @param n_species number of species in that clade
#' @param carr_cap carrying capacity, number of species this clade will
#'   grow to
#' @return the clade's cladogenesis rate, which is at least zero. This
#'   rate will be zero if there are more species than the carrying capacity
#'   allows for
#' @note For clade-specific carrying capacity, 
#'   each clade is simulated seperately in \code{\link{DAISIE_sim}}
#' @author Richel J.C. Bilderbeek
#' @examples 
#'   testit::assert(
#'     DAISIE_calc_clade_clado_rate(
#'       ps_clado_rate = 0.2,  
#'       n_species = 5,
#'       carr_cap = 10
#'     ) == 0.5
#'   )
#'   testit::assert(
#'     DAISIE_calc_clade_clado_rate(
#'       ps_clado_rate = 0.2,  
#'       n_species = 2,
#'       carr_cap = 1
#'     ) == 0.0
#'   )
#' @export
DAISIE_calc_clade_clado_rate <- function(ps_clado_rate, n_species, carr_cap) {
  testit::assert(ps_clado_rate >= 0.0)
  testit::assert(n_species >= 0)
  testit::assert(carr_cap >= 0)
  max(
    0.0,
    n_species * ps_clado_rate * (1.0 - (n_species / carr_cap))
  )
}

#' Calculate the clade-wide immigration rate.
#' @param ps_imm_rate per species immigration rate
#' @param n_island_species number of species in that clade on the island
#' @param n_mainland_species number of species in that clade on the mainland
#' @param carr_cap carrying capacity, number of species this clade will
#'   grow to
#' @return the clade's immigration rate, which is at least zero. This
#'   rate will be zero if there are more species than the carrying capacity
#'   allows for
#' @author Richel J.C. Bilderbeek
#' @examples 
#'   testit::assert(
#'     DAISIE_calc_clade_imm_rate(
#'       ps_imm_rate = 0.1, 
#'       n_island_species = 5, 
#'       n_mainland_species = 2, 
#'       carr_cap = 10
#'     ) == 0.1
#'   )
#'   testit::assert(
#'     DAISIE_calc_clade_imm_rate(
#'       ps_imm_rate = 0.1, 
#'       n_island_species = 5, 
#'       n_mainland_species = 2, 
#'       carr_cap = 1
#'     ) == 0.0
#'   )
#' @export
DAISIE_calc_clade_imm_rate <- function(
  ps_imm_rate, 
  n_island_species, 
  n_mainland_species, 
  carr_cap
) {
  testit::assert(ps_imm_rate >= 0.0)
  testit::assert(n_island_species >= 0)
  testit::assert(n_mainland_species >= 0)
  testit::assert(carr_cap >= 0)
  max(
    0.0,
     n_mainland_species * ps_imm_rate * (1.0 - (n_island_species / carr_cap))
  )
}
