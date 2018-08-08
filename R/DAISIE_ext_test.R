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
  stt[1,] <- c(1000, 0)
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
        possible_event <- sample(1:4, 1, prob = c(rates[[1]], rates[[2]], 
                                                  rates[[3]], rates[[4]]), 
                                 replace = FALSE)
      } else if (sum(rates[[1]], rates[[2]], 
                     rates[[3]], rates[[4]], 
                     rates[[5]]) > 0){
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

#' Calculates algorithm rates
#' @description Internal function that updates the all the rates and 
#' max extinction horizon at time t.
#' @family rates calculation
#' @param timeval current time of simulation
#' @param totaltime total time of simulation
#' @param gam per capita immigration rate
#' @param mu per capita extinction rate in no ontogeny model
#' @param laa per capita anagenesis rate
#' @param lac per capita cladogenesis rate
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
#' @param K carrying capacity
#' @param island_spec matrix containing state of system
#' @param mainland_n total number of species present in the mainland
#' @param thor time of horizon for max extinction
update_rates <- function(timeval, totaltime,
                         gam, mu, laa, lac, Apars, Epars,
                         island_ontogeny, 
                         extcutoff,
                         K, 
                         island_spec, mainland_n, thor) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  
  
  immig_rate <- get_immig_rate(gam = gam,
                               island_spec = island_spec,
                               K = K,
                               mainland_n = mainland_n)
  
  ext_rate <- get_ext_rate(timeval = timeval,
                           totaltime = totaltime,
                           mu = mu,
                           Apars = Apars,
                           Epars = Epars, 
                           island_ontogeny = island_ontogeny, 
                           extcutoff = extcutoff,
                           island_spec = island_spec,
                           K = K)
  
  ana_rate <- get_ana_rate(laa = laa, island_spec = island_spec)
  
  clado_rate <- get_clado_rate(lac = lac,
                               island_ontogeny = island_ontogeny, 
                               island_spec = island_spec,
                               K = K)
  
  if (is.null(island_ontogeny)) {
    
    ext_rate_max <- ext_rate
    
  } else if ((Apars$proportional_peak_t * Apars$total_island_age) > timeval) {
    
    ext_rate_max <- ext_rate
    
  } else {
    
    ext_rate_max <- get_ext_rate(timeval = thor, totaltime = totaltime, mu = mu,
                                 Apars = Apars, Epars = Epars,
                                 island_ontogeny = island_ontogeny, 
                                 extcutoff = extcutoff, island_spec = island_spec,
                                 K = K)
  }
  
  rates <- list(immig_rate, ext_rate, ana_rate, clado_rate, ext_rate_max)
  return(rates)
}

#' Something
#' @param rates something
#' @param timeval something
pick_timeval <- function(rates, timeval) {
  # Calculates when next event will happen
  totalrate <- rates[[1]] + rates[[3]] + rates[[4]] + rates[[5]]
  dt <- rexp(1, totalrate)
  timeval <- timeval + dt
  return(timeval)
}

#' Does something
#' @param timeval something
#' @param possible_event something
#' @param maxspecID something
#' @param mainland_spec something
#' @param island_spec something
DAISIE_sim_update_state <- function(timeval, possible_event,maxspecID,mainland_spec,island_spec)
{  
  if (possible_event > 4) {
    # Nothing happens
  }
  ##########################################
  #IMMIGRATION
  if(possible_event == 1)
  {  	
    colonist = DDD::sample2(mainland_spec,1)
    
    if(length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
    } else
    {
      isitthere = c()
    }
    
    if(length(isitthere) == 0)
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA))
    }
    
    if(length(isitthere) != 0)
    {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
    }
  }
  
  ##########################################
  #EXTINCTION
  if(possible_event == 2)
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
  return(list(island_spec = island_spec, maxspecID = maxspecID))
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
