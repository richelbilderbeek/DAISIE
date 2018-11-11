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
    rates <- update_rates(timeval = timeval,
                          totaltime = totaltime,
                          gam = gam,
                          mu = mu,
                          laa = laa,
                          lac = lac,
                          K = K, 
                          extcutoff = NULL,
                          island_spec = island_spec,
                          mainland_n = mainland_n)
    
    timeval_and_dt <- calc_next_timeval(rates, timeval)
    timeval <- timeval_and_dt$timeval
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

DAISIE_ONEcolonist <- function(time,island_spec,stt_table)
{
  totaltime <- time
  ### number of independent colonisations
  uniquecolonisation <- as.numeric(unique(island_spec[,"Colonisation time (BP)"]))
  number_colonisations <- length(uniquecolonisation)

  ### if there is only one independent colonisation - anagenetic and cladogenetic
  #species are classed as stac=2; immigrant classed as stac=4:
  if(number_colonisations == 1)
  {
    if(island_spec[1,"Species type"] == "I")
    {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(totaltime,as.numeric(island_spec[1,"Colonisation time (BP)"])),
                          stac = 4,
                          missing_species = 0)
    }
    if(island_spec[1,"Species type"] == "A")
    {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(totaltime,as.numeric(island_spec[1,"Colonisation time (BP)"])),
                          stac = 2,
                          missing_species = 0)
    }
    if(island_spec[1,"Species type"] == "C")
    {
      descendants <- list(stt_table = stt_table,
                          branching_times = c(totaltime,rev(sort(as.numeric(island_spec[,"branching time (BP)"])))),
                          stac = 2,
                          missing_species = 0)
    }
  }

  ### if there are two or more independent colonisations, all species are classed as stac=3 and put within same list item:
  else if(number_colonisations > 1)
  {
    print("trigger1")
    descendants <- list(stt_table = stt_table,
                        branching_times = NA,stac = 2,missing_species = 0,
                        other_clades_same_ancestor = list())
    ### create table with information on other clades with same ancestor, but this information is not used in DAISIE_ML
    oldest <- which(as.numeric(island_spec[,"Colonisation time (BP)"]) == max(as.numeric(island_spec[,"Colonisation time (BP)"])))

    oldest_table <- island_spec[oldest,]
    if(class(oldest_table) == 'character')
    {
      oldest_table <- t(as.matrix(oldest_table))
    }
    if(oldest_table[1,'Species type'] == 'A')
    {
      descendants$branching_times <- c(totaltime, as.numeric(oldest_table[1,"Colonisation time (BP)"]))
    } else if(oldest_table[1,'Species type'] == 'C')
    {
      descendants$branching_times <- c(totaltime, rev(sort(as.numeric(oldest_table[,'branching time (BP)']))))
    }

    youngest_table = island_spec[-oldest,]
    if(class(youngest_table) == 'character')
    {
      youngest_table <- t(as.matrix(youngest_table))
    }

    uniquecol <- as.numeric(unique(youngest_table[,"Colonisation time (BP)"]))

    descendants$missing_species <- length(which(youngest_table[,"Species type"]!='I'))
    for(colonisation in 1:length(uniquecol))
    {
      descendants$other_clades_same_ancestor[[colonisation]] <- list(brts_miss = NA,species_type = NA)

      samecolonisation <- which(as.numeric(youngest_table[,"Colonisation time (BP)"]) == uniquecol[colonisation])

      if(youngest_table[samecolonisation[1],"Species type"] == "I")
      {
        descendants$stac <- 3
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "I"
      } else if(youngest_table[samecolonisation[1],"Species type"] == "A")
      {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- as.numeric(youngest_table[samecolonisation,"Colonisation time (BP)"])
        descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "A"
      } else if (youngest_table[samecolonisation[1],"Species type"] == "C")
      {
        descendants$other_clades_same_ancestor[[colonisation]]$brts_miss <- rev(sort(as.numeric(youngest_table[samecolonisation,"branching time (BP)"])))
        descendants$other_clades_same_ancestor[[colonisation]]$species_type <- "C"
      }
    }
  }
  descendants
}
