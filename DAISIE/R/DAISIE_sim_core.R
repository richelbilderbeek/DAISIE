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
#' @param Apars A named list containing area parameters as created by create_area_params:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param Epars A numeric vector:
#' \itemize{
#'   \item{[1]: minimum extinction when area is at peak}
#'   \item{[2]: extinction rate when current area is 0.10 of maximum area}
#' }
#' @param island_ontogeny A string describing the type of island ontogeny. Can be \code{NULL},
#' \code{quadratic} for a beta function describing area through time,
#'  or \code{linear} for a linear function
DAISIE_sim_core <- function(
  time,
  mainland_n,
  pars,
  Apars = NULL,
  Epars = NULL,
  island_ontogeny = NULL,
  keep_final_state = FALSE,
  island_spec = NULL,
  stt_table = NULL
) {
  testit::assert(is.logical(keep_final_state))
  testit::assert(length(pars) == 5)
  testit::assert(is.null(Apars) || are_area_params(Apars))
  testit::assert(is.null(island_spec) || is.matrix(island_spec))
  
  if (pars[4] == 0) {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  if (!is.null(Apars) && is.null(island_ontogeny)) {
    stop("Apars specified for constant island_ontogeny. Set Apars to NULL.")
  }
  
  if ((is.null(Epars) || is.null(Apars)) && !is.null(island_ontogeny)) {
    stop("Island ontogeny specified but Area parameters and/or extinction 
         parameters not available. Please either set island_ontogeny to NULL, or 
         specify Apars and Epars.")
  }
  
  if (!is.null(island_ontogeny) && island_ontogeny != "linear" && island_ontogeny != "quadratic") {
    stop("Please select valid island ontogeny model. Options are NULL, 'linear' or 'quadratic'.")
  }
  
  timeval <- 0
  totaltime <- time
  lac <- pars[1]
  mu <- pars[2]
  K <- pars[3]
  gam <- pars[4]
  laa <- pars[5]
  extcutoff <- max(1000, 1000 * (laa + lac + gam))
  testit::assert(is.numeric(extcutoff))
  ext_multiplier <- 0.5
  testit::assert((totaltime <= Apars$total_island_age) || is.null(Apars))
  
  #### Start Gillespie ####
  if (is.null(stt_table)) {
    island_spec = c()
    stt_table <- matrix(ncol = 4)
    colnames(stt_table) <- c("Time","nI","nA","nC")
    stt_table[1,] <- c(totaltime,0,0,0)
  } else {
    stt_table <- stt_table # Remove this later
    colnames(stt_table) <- c("Time","nI","nA","nC")
    stt_table[1,1] <- c(totaltime)
    stt_table <- stt_table[nrow(stt_table), ]
  }
  
  mainland_spec <- seq(1, mainland_n, 1)
  maxspecID <- mainland_n
  
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # Pick t_hor (before timeval, to set Amax t_hor)
  t_hor <- get_t_hor(
    timeval = 0,
    totaltime = totaltime,
    Apars = Apars,
    ext = 0,
    ext_multiplier = ext_multiplier,
    island_ontogeny = island_ontogeny, 
    t_hor = NULL,
    dt = 0,
    old_timeval = 0
  )
  
  while (timeval < totaltime) {
    # Calculate rates
    rates <- update_rates(timeval = timeval,
                          totaltime = totaltime,
                          gam = gam,
                          mu = mu,
                          laa = laa,
                          lac = lac,
                          Apars = Apars,
                          Epars = Epars,
                          island_ontogeny = island_ontogeny,
                          extcutoff = extcutoff,
                          K = K,
                          island_spec = island_spec,
                          mainland_n = mainland_n,
                          t_hor = t_hor)
    
    old_timeval <- timeval
    timeval_and_dt <- calc_next_timeval(rates, timeval)
    timeval <- timeval_and_dt$timeval
    dt <- timeval_and_dt$dt
    
    if (timeval <= t_hor) {
      testit::assert(are_rates(rates))
      
      # Determine event
      possible_event <- DAISIE_sample_event(
        rates = rates,
        island_ontogeny = island_ontogeny
      )
      
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
    } else {
      #### After t_hor is reached ####
      
      timeval <- t_hor
      t_hor <- get_t_hor(
        timeval = timeval,
        totaltime = totaltime,
        Apars = Apars,
        ext = rates$ext_rate,
        ext_multiplier = ext_multiplier,
        island_ontogeny = island_ontogeny, 
        t_hor = t_hor,
        dt = dt,
        old_timeval = old_timeval
      )
    }
    # TODO Check if this is redundant, or a good idea
    if (rates$ext_rate_max >= extcutoff && length(island_spec[,1]) == 0) {
      timeval <- totaltime
    }
  }
  
  # Finalize stt_table 
  stt_table <- rbind(stt_table, 
                     c(0, 
                       stt_table[nrow(stt_table), 2],
                       stt_table[nrow(stt_table), 3],
                       stt_table[nrow(stt_table), 4]))
  
  island <- DAISIE_create_island(stt_table = stt_table,
                                 totaltime = totaltime,
                                 island_spec = island_spec,
                                 mainland_n = mainland_n,
                                 keep_final_state = keep_final_state)
  return(island)
}





#' Calculates when the next timestep will be.
#'
#' @param rates list of numeric with probabilities of each event
#' @param timeval current time of simulation
calc_next_timeval <- function(rates, timeval) {
  # Calculates when next event will happen
  testit::assert(are_rates(rates))
  testit::assert(timeval >= 0)
  totalrate <- rates$immig_rate_max + rates$ana_rate + rates$clado_rate_max + rates$ext_rate_max
  dt <- rexp(1, totalrate)
  timeval <- timeval + dt
  return(list(timeval = timeval, dt = dt))
}

#' Updates state of island given sampled event
#' 
#' Makes the event happen by updating island species matrix and species IDs.
#' What event happens is determined by the sampling in the algorithm.
#' !!!!!!!!!!!!THIS DOCUMENTATION MUST BE CONFIRMED!!!!!!!!!!!!
#' @param timeval current time of simulation
#' @param totaltime simulated amount of time
#' @param possible_event numeric indicating what event will happen.
#' @param maxspecID current species IDs
#' @param mainland_spec number of mainland species
#' @param island_spec A Matrix with species on island (state of system at each time point)
#' @param stt_table A species-through-time table
DAISIE_sim_update_state <- function(timeval,
                                    totaltime,
                                    possible_event,
                                    maxspecID,
                                    mainland_spec,
                                    island_spec,
                                    stt_table)
{ 
  if (possible_event > 4) {
    # Nothing happens
  }
  
  ##########################################
  #IMMIGRATION
  if (possible_event == 1)
  {
    colonist = DDD::sample2(mainland_spec,1)
    
    if (length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }
    
    if (length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA))
    }
    
    if (length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
    }
  }
  
  ##########################################
  #EXTINCTION
  if (possible_event == 2)
  { 	
    extinct = DDD::sample2(1:length(island_spec[,1]),1)
    #this chooses the row of species data to remove
    
    typeofspecies = island_spec[extinct,4]
    
    if(typeofspecies == "I")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove immigrant
    
    if(typeofspecies == "A")
    {
      island_spec = island_spec[-extinct,]
    }
    #remove anagenetic
    
    if(typeofspecies == "C")
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival totaltime
      sisters = intersect(which(island_spec[,2] == island_spec[extinct,2]),which(island_spec[,3] == island_spec[extinct,3]))
      survivors = sisters[which(sisters != extinct)]
      
      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic	
        island_spec[survivors,4] = "A"
        island_spec[survivors,c(5,6)] = c(NA,NA)
        island_spec[survivors,7] = "Clado_extinct"
        island_spec = island_spec[-extinct,]
      }
      
      if(length(sisters) >= 3)
      {		
        numberofsplits = nchar(island_spec[extinct,5])
        
        mostrecentspl = substring(island_spec[extinct,5],numberofsplits)
        
        if(mostrecentspl=="B")
        { 
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }
        
        motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")
        
        possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]
        
        #different rules depending on whether a B or A is removed. B going extinct is simpler because it only 
        #carries a record of the most recent speciation			
        if(mostrecentspl == "A")
        {								
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec[possiblesister,6] == max(as.numeric(island_spec[possiblesister,6])))]
          island_spec[tochange,6] = island_spec[extinct,6]	
        }
        
        #remove the offending A/B from these species
        island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),
                                              substring(island_spec[possiblesister,5],numberofsplits + 1,
                                                        nchar(island_spec[possiblesister,5])),sep = "")	
        island_spec = island_spec[-extinct,]
      }
    }
    island_spec = rbind(island_spec)	
  }
  
  ##########################################
  #ANAGENESIS
  if(possible_event == 3)
  {    
    immi_specs = which(island_spec[,4] == "I")
    
    #we only allow immigrants to undergo anagenesis
    if(length(immi_specs) == 1)
    {
      anagenesis = immi_specs
    }
    if(length(immi_specs) > 1)
    {
      anagenesis = DDD::sample2(immi_specs,1)
    }
    
    maxspecID = maxspecID + 1
    island_spec[anagenesis,4] = "A"
    island_spec[anagenesis,1] = maxspecID
    island_spec[anagenesis,7] = "Immig_parent"
  }
  
  ##########################################
  #CLADOGENESIS - this splits species into two new species - both of which receive 
  if(possible_event == 4)
  { 		
    tosplit = DDD::sample2(1:length(island_spec[,1]),1)
    
    #if the species that speciates is cladogenetic
    if(island_spec[tosplit,4] == "C")
    {
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      oldstatus = island_spec[tosplit,5]
      island_spec[tosplit,5] = paste(oldstatus,"A",sep = "")
      #island_spec[tosplit,6] = timeval
      island_spec[tosplit,7] = NA
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],
                                        "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      
      maxspecID = maxspecID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      
      island_spec[tosplit,4] = "C"
      island_spec[tosplit,1] = maxspecID + 1
      island_spec[tosplit,5] = "A"
      island_spec[tosplit,6] = island_spec[tosplit,3]
      island_spec[tosplit,7] = NA
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplit,2],island_spec[tosplit,3],"C","B",timeval,NA))
      
      maxspecID = maxspecID + 2
    } 
  }
  
  
  if (possible_event <= 4) {
    stt_table <- rbind(stt_table,
                       c(totaltime - timeval,
                         length(which(island_spec[,4] == "I")),
                         length(which(island_spec[,4] == "A")),
                         length(which(island_spec[,4] == "C"))))
  }
  
  updated_state <- list(island_spec = island_spec, 
                        maxspecID = maxspecID, 
                        stt_table = stt_table)
  updated_state
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