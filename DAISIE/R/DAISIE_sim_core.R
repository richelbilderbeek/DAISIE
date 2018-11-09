#' Internal function of the DAISIE simulation
#' @param time Simulated amount of time
#' @param mainland_n A numeric stating the number of mainland species, that
#'   is, the number of species that can potentially colonize the island.
#'   If \code{\link{DAISIE_sim}} uses a clade-specific diversity dependence,
#'   this value is set to 1. 
#'   If \code{\link{DAISIE_sim}} uses an island-specific diversity dependence,
#'   this value is set to the number of mainland species.
#' @param pars A numeric vector:
#' \itemize{
#'   \item{[1]: cladogenesis rate}
#'   \item{[2]: extinction rate}
#'   \item{[3]: carrying capacity}
#'   \item{[4]: immigration rate}
#'   \item{[5]: anagenesis rate}
#' }
DAISIE_sim_core <- function(
  time,
  mainland_n,
  pars) {
  testit::assert(length(pars) == 5)
  
  if (pars[4] == 0) {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  timeval <- 0
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  
  mainland_spec <- seq(1, mainland_n, 1)
  maxspecID <- mainland_n
  
  island_spec = c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time","nI","nA","nC")
  stt_table[1,] <- c(totaltime,0,0,0)
  
  #### Start Gillespie ####
  while (timeval < totaltime) {
    # Calculate rates
    rates <- update_rates(tiImeval = timeval,
                          totaltime = totaltime,
                          gam = gam,
                          mu = mu,
                          laa = laa,
                          lac = lac,
                          K = K, 
                          extcutoff = NULL,
                          island_spec = island_spec,
                          mainland_n = mainland_n)
    
    old_timeval <- timeval
    timeval_and_dt <- calc_next_timeval(rates, timeval)
    timeval <- timeval_and_dt$timeval
    dt <- timeval_and_dt$dt
    
    testit::assert(are_rates(rates))
    
    # Determine event
    possible_event <- DAISIE_sample_event(
      rates = rates)
    
    updated_state <- DAISIE_sim_update_state(
      timeval = timeval, 
      totaltime = totaltime,
      possible_event = possible_event,
      maxspecID = maxspecID,
      mainland_spec = mainland_spec,
      island_spec = island_spec,
      stt_table = stt_table
    )
    
    island_spec <- updated_state$island_spec
    maxspecID <- updated_state$maxspecID
    stt_table <- updated_state$stt_table
    
  } 
  
  #### Finalize stt_table ####
  stt_table <- rbind(stt_table, 
                     c(0, 
                       stt_table[nrow(stt_table), 2],
                       stt_table[nrow(stt_table), 3],
                       stt_table[nrow(stt_table), 4]))
  
  island <- DAISIE_create_island(stt_table = stt_table,
                                 totaltime = totaltime,
                                 island_spec = island_spec,
                                 mainland_n = mainland_n)
  return(island)
}

