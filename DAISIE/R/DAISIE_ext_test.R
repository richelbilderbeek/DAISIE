#' Internal function of the DAISIE simulation
#' @param time simulated amount of time
#' @param mainland_n number of mainland species, that
#'   is, the number of species that can potentially colonize the island.
#'   If \code{\link{DAISIE_sim}} uses a clade-specific diversity dependence,
#'   this value is set to 1. 
#'   If \code{\link{DAISIE_sim}} uses an island-specific diversity dependence,
#'   this value is set to the number of mainland species.
#' @param pars a numeric vector:
#' \itemize{
#'   \item{[1]: cladogenesis rate}
#'   \item{[2]: extinction rate}
#'   \item{[3]: carrying capacity}
#'   \item{[4]: immigration rate}
#'   \item{[5]: anagenesis rate}
#' }
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
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
DAISIE_exinction_test <- function(
  time,
  mainland_n,
  pars,
  Apars = NULL,
  Epars = NULL,
  island_ontogeny = NULL
) {
  timeval <- 0
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  
  extcutoff <- max(1000, 1000 * (laa + lac + gam))
  ext_multiplier <- 0.5
  stt <- matrix(ncol = 2)
  # if(pars[4] == 0) 
  # {
  #   stop('Rate of colonisation is zero. Island cannot be colonised.')
  # }  
  
  if (are_area_params(Apars) && is.null(island_ontogeny)){
    stop("Apars specified for contant island_ontogeny. Set Apars to NULL")
  }
  
  if (!is.null(island_ontogeny) && island_ontogeny != "linear" && island_ontogeny != "quadratic") {
    stop("Please select valid island ontogeny model. Options are no ontogeny: NULL, 'linear' or 'quadratic'.")
  }
  
  mainland_spec <- seq(1, mainland_n, 1)
  maxspecID <- mainland_n
  
  island_spec = matrix(ncol = 7, nrow = 1000)
  island_spec[,4] = "I"
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time","nI","nA","nC")
  stt_table[1,] <- c(totaltime,0,0,0)
  
  # Pick thor (before timeval, to set Amax thor)
  thor <- get_thor(0, totaltime, Apars, ext_multiplier, island_ontogeny, thor = NULL)
  
  #### Start Gillespie ####
  while(timeval < totaltime) {
    if (timeval < thor) {
      rates <- update_rates(timeval = timeval, totaltime = totaltime, gam = gam,
                            mu = mu, laa = laa, lac = lac, Apars = Apars,
                            Epars = Epars, island_ontogeny = island_ontogeny,
                            extcutoff = extcutoff, K = K,
                            island_spec = island_spec, mainland_n, thor)
      if (is.na(timeval) == T) {
        timeval <- totaltime
      } else {
      timeval <- pick_timeval(rates, timeval)
      }
      # Determine event
      # If statement prevents odd behaviour of sample when rates are 0
      if (is.null(island_ontogeny)) {
        possible_event <- sample(
          1:4, 1, 
          prob = c(rates[[1]], rates[[2]], 
                   rates[[3]], rates[[4]]), 
                                 replace = FALSE)
      } else if (sum(rates[[1]], rates[[2]], 
                     rates[[3]], rates[[4]], 
                     rates[[5]]) > 0) {
        possible_event <- sample(1:5, 1, prob = c(rates[[1]], rates[[2]], 
                                                  rates[[3]], rates[[4]], 
                                                  (rates[[5]] - rates[[2]])),
                                 replace = FALSE)
      }
      if (is.nan(timeval) == T) {
        timeval <- totaltime
      }

      if (timeval < totaltime) {
        # Run event
        
        
        new_state <- DAISIE_sim_update_state(timeval = timeval,
                                             possible_event = possible_event,
                                             maxspecID = maxspecID,
                                             mainland_spec = mainland_spec,
                                             island_spec = island_spec)
        
        island_spec <- new_state$island_spec
        maxspecID <- new_state$maxspecID
        nspec <- nrow(island_spec)
        stt <- rbind(stt, c(nspec, timeval))
      }
      stt_table <- rbind(stt_table,
                         c(totaltime - timeval,
                           length(which(island_spec[,4] == "I")),
                           length(which(island_spec[,4] == "A")),
                           length(which(island_spec[,4] == "C"))))

    } else {
      ##### After thor is reached ####
      # Recalculate thor
      thor <- get_thor(timeval = timeval, totaltime = totaltime, Apars = Apars,
                       ext_multiplier = ext_multiplier,
                       island_ontogeny = island_ontogeny, thor = thor)
    }
  }
  
  return(stt)
}

