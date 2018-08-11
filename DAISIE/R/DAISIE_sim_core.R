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
DAISIE_sim_core <- function(
  time,
  mainland_n,
  pars,
  Apars = NULL,
  Epars = NULL,
  island_ontogeny = NULL
) {
  testit::assert(length(pars) == 5)
  testit::assert(is.null(Apars) || are_area_params(Apars))
  
  if(pars[4] == 0) 
  {
    stop('Rate of colonisation is zero. Island cannot be colonised.')
  }  
  
  if (are_area_params(Apars) && is.null(island_ontogeny)){
    stop("Apars specified for contant island_ontogeny. Set Apars to NULL")
  }
  
  if (!is.null(island_ontogeny) && island_ontogeny != "linear" && island_ontogeny != "quadratic") {
    stop("Please select valid island ontogeny model. Options are no ontogeny: NULL, 'linear' or 'quadratic'.")
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
  
  
  mainland_spec <- seq(1, mainland_n, 1)
  maxspecID <- mainland_n
  
  island_spec = c()
  stt_table <- matrix(ncol = 4)
  colnames(stt_table) <- c("Time","nI","nA","nC")
  stt_table[1,] <- c(totaltime,0,0,0)
  testit::assert(is.null(Apars) || are_area_params(Apars))
  # Pick thor (before timeval, to set Amax thor)
  thor_ext <- get_thor(
    0,
    totaltime,
    Apars,
    ext_multiplier,
    island_ontogeny, 
    thor = NULL
  )
  thor_c_i <- get_thor_half(0,
                          totaltime,
                          Apars,
                          ext_multiplier,
                          island_ontogeny, 
                          thor_c_i = NULL
                          )
  

  #### Start Gillespie ####
  while (timeval <= totaltime) {
    if (timeval < thor_ext) {
      if (timeval < thor_c_i) {
        
        

      rates <- update_rates(timeval = timeval, totaltime = totaltime, gam = gam,
                            mu = mu, laa = laa, lac = lac, Apars = Apars,
                            Epars = Epars, island_ontogeny = island_ontogeny,
                            extcutoff = extcutoff, K = K,
                            island_spec = island_spec, mainland_n, thor_ext, thor_c_i)
      testit::assert(are_rates(rates))

      timeval <- calc_next_timeval(rates, timeval)

      # Determine event
      # If statement prevents odd behaviour of sample when rates are 0
      if (is.null(island_ontogeny)) {
        possible_event <- sample(1:4, 1, prob = c(rates$immig_rate,
                                                  rates$ext_rate,
                                                  rates$ana_rate,
                                                  rates$clado_rate), 
                                 replace = FALSE)
      } else {

      testit::assert(are_rates(rates))  
        possible_event <- sample(1:7, 1, prob = c(
          rates$immig_rate,
          rates$ext_rate,
          rates$ana_rate,
          rates$clado_rate,
          (rates$ext_rate_max - rates$ext_rate),
          (rates$immig_rate_max - rates$immig_rate),
          (rates$clado_rate_max - rates$clado_rate)),
          replace = FALSE)
        
      }

      if (timeval <= totaltime) {
        # Run event

        new_state <- DAISIE_sim_update_state(timeval = timeval,
                                             possible_event = possible_event,
                                             maxspecID = maxspecID,
                                             mainland_spec = mainland_spec,
                                             island_spec = island_spec)
        
        island_spec <- new_state$island_spec
        maxspecID <- new_state$maxspecID

      }
      stt_table <- rbind(stt_table,
                         c(totaltime - timeval,
                           length(which(island_spec[,4] == "I")),
                           length(which(island_spec[,4] == "A")),
                           length(which(island_spec[,4] == "C"))))
      } else {
        thor_c_i <- get_thor_half(
          timeval = timeval,
          totaltime = totaltime,
          Apars = Apars,
          ext_multiplier = ext_multiplier,
          island_ontogeny = island_ontogeny, 
          thor_c_i = thor_c_i
          )
      }
    } else {
      #### After thor is reached ####
      # Recalculate thor
      testit::assert(are_area_params(Apars))
      thor_ext <- get_thor(timeval = timeval,
                       totaltime = totaltime,
                       Apars = Apars,
                       ext_multiplier = ext_multiplier,
                       island_ontogeny = island_ontogeny, 
                       thor = thor_ext)
    }
  }
  stt_table[nrow(stt_table),1] <- 0
  
  ############# 
  ### if there are no species on the island branching_times = island_age, stac = 0, missing_species = 0 
  if(length(island_spec[,1]) == 0)
  {
    island <- list(stt_table = stt_table, branching_times = totaltime, stac = 0, missing_species = 0)
  } else
  {
    cnames <- c("Species","Mainland Ancestor","Colonisation time (BP)",
                "Species type","branch_code","branching time (BP)","Anagenetic_origin")
    colnames(island_spec) <- cnames
    
    ### set ages as counting backwards from present
    island_spec[,"branching time (BP)"] <- totaltime - as.numeric(island_spec[,"branching time (BP)"])
    island_spec[,"Colonisation time (BP)"] <- totaltime - as.numeric(island_spec[,"Colonisation time (BP)"])
    if(mainland_n == 1)
    {
      island <- DAISIE_ONEcolonist(totaltime,island_spec,stt_table)
    } else if (mainland_n > 1) 
    {  
      ### number of colonists present
      colonists_present <- sort(as.numeric(unique(island_spec[,'Mainland Ancestor'])))
      number_colonists_present <- length(colonists_present) 
      
      island_clades_info <- list()  
      for(i in 1:number_colonists_present)
      {
        subset_island <- island_spec[which(island_spec[,'Mainland Ancestor']==colonists_present[i]),] 
        if(class(subset_island) != 'matrix')
        {
          subset_island <- rbind(subset_island[1:7])
          colnames(subset_island) <- cnames
        }
        island_clades_info[[i]] <- DAISIE_ONEcolonist(totaltime,island_spec=subset_island,stt_table=NULL)
        island_clades_info[[i]]$stt_table <- NULL
      }
      island <- list(stt_table = stt_table, taxon_list = island_clades_info)
    }
  }
  return(island)
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
#' large and slowing down simulation. Should be big
#' @param K carrying capacity
#' @param island_spec matrix containing state of system
#' @param mainland_n total number of species present in the mainland
#' @param thor_ext time of horizon for max extinction
#' @param thor_c_i time of horizon for max cladogenesis and immigration
update_rates <- function(timeval, totaltime,
                         gam, mu, laa, lac, Apars, Epars,
                         island_ontogeny, 
                         extcutoff,
                         K, 
                         island_spec, mainland_n, thor_ext, thor_c_i) {
  # Function to calculate rates at time = timeval. Returns list with each rate.
  testit::assert(is.numeric(timeval))
  testit::assert(is.numeric(totaltime))
  testit::assert(is.numeric(gam))
  testit::assert(is.numeric(mu))
  testit::assert(is.numeric(laa))
  testit::assert(is.numeric(lac))
  testit::assert(is.null(Apars) || are_area_params(Apars))
  testit::assert(is.null(Epars) || is.numeric(Epars))
  testit::assert(is.character(island_ontogeny) || is.null(island_ontogeny))
  testit::assert(is.numeric(extcutoff))
  testit::assert(is.numeric(K))
  testit::assert(is.matrix(island_spec) || is.null(island_spec))
  testit::assert(is.numeric(mainland_n))
  testit::assert(is.numeric(thor_ext))
  testit::assert(is.numeric(thor_c_i))
  
  immig_rate <- get_immig_rate(timeval = timeval,
                               totaltime = totaltime,
                               gam = gam,
                               Apars = Apars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K,
                               mainland_n = mainland_n)
  testit::assert(is.numeric(immig_rate))
  ext_rate <- get_ext_rate(timeval = timeval,
                           totaltime = totaltime,
                           mu = mu,
                           Apars = Apars,
                           Epars = Epars, 
                           island_ontogeny = island_ontogeny, 
                           extcutoff = extcutoff,
                           island_spec = island_spec,
                           K = K)
  testit::assert(is.numeric(ext_rate))
  
  ana_rate <- get_ana_rate(laa = laa,
                           island_spec = island_spec)
  testit::assert(is.numeric(ana_rate))
  clado_rate <- get_clado_rate(timeval = timeval, 
                               totaltime = totaltime,
                               lac = lac,
                               Apars = Apars,
                               island_ontogeny = island_ontogeny,
                               island_spec = island_spec,
                               K = K)
  testit::assert(is.numeric(clado_rate))
  if (is.null(island_ontogeny)) {
    
    immig_rate_max <- immig_rate
    testit::assert(is.numeric(immig_rate_max))
    ext_rate_max <- ext_rate
    testit::assert(is.numeric(ext_rate_max))
    clado_rate_max <- clado_rate
    testit::assert(is.numeric(clado_rate_max))
  } else if ((Apars$proportional_peak_t * Apars$total_island_age) > timeval) {

    ext_rate_max <- ext_rate
    testit::assert(is.numeric(ext_rate_max))
  } else {
    # No ontogeny, max rate is thor, which in this case is totaltime (from get_thor)
    immig_rate_max <- get_immig_rate(timeval = thor_c_i,
                                     totaltime = totaltime,
                                     gam = gam,
                                     Apars = Apars,
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = 0.05, 
                                     mainland_n = mainland_n)
    testit::assert(is.numeric(immig_rate_max))

    ext_rate_max <- get_ext_rate(timeval = thor_ext,
                                 totaltime = totaltime,
                                 mu = mu,
                                 Apars = Apars, 
                                 Epars = Epars,
                                 island_ontogeny = island_ontogeny, 
                                 extcutoff = extcutoff, 
                                 island_spec = island_spec,
                                 K = K)
    testit::assert(is.numeric(ext_rate_max) && ext_rate_max >= 0.0)
    clado_rate_max <- get_clado_rate(timeval = thor_c_i, 
                                     totaltime = totaltime,
                                     lac = lac,
                                     Apars = Apars,
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = 0.05)
    testit::assert(is.numeric(clado_rate_max))
  }
  
  if ((((Apars$proportional_peak_t * Apars$total_island_age) / 2) > timeval) && !is.null(island_ontogeny)) {
    clado_rate_max <- get_clado_rate(((Apars$proportional_peak_t * Apars$total_island_age) / 2),
                                     totaltime = totaltime, 
                                     lac = lac,
                                     Apars = Apars, 
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = 0.05)
    testit::assert(is.numeric(clado_rate_max))
    
    immig_rate_max <- get_immig_rate(((Apars$proportional_peak_t * Apars$total_island_age) / 2),
                                     totaltime = totaltime, 
                                     gam = gam,
                                     Apars = Apars, 
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = 0.05,
                                     mainland_n = mainland_n)
    testit::assert(is.numeric(immig_rate_max))
  } else {

    clado_rate_max <- get_clado_rate(thor_c_i,
                                     totaltime = totaltime, 
                                     lac = lac,
                                     Apars = Apars, 
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = K)
    testit::assert(is.numeric(clado_rate_max))
    immig_rate_max <- get_immig_rate(thor_c_i,
                                     totaltime = totaltime, 
                                     gam = gam,
                                     Apars = Apars, 
                                     island_ontogeny = island_ontogeny,
                                     island_spec = island_spec,
                                     K = K,
                                     mainland_n = mainland_n)
    testit::assert(is.numeric(immig_rate_max))
    }


  rates <- create_rates(
    immig_rate = immig_rate,
    ext_rate = ext_rate,
    ana_rate = ana_rate,
    clado_rate = clado_rate,
    ext_rate_max = ext_rate_max,
    immig_rate_max = immig_rate_max,
    clado_rate_max = clado_rate_max
  )

  return(rates)
}

#' Calculates when the next timestep will be.
#' @param rates list of numeric with probabilities of each event
#' @param timeval current time of simulation
calc_next_timeval <- function(rates, timeval) {
  # Calculates when next event will happen
  testit::assert(are_rates(rates))
  testit::assert(timeval >= 0)
  totalrate <- rates$immig_rate_max + rates$ana_rate + rates$clado_rate_max + rates$ext_rate_max
  dt <- rexp(1, totalrate)
  timeval <- timeval + dt
  return(timeval)
}

#' Updates state of island given sampled event
#' 
#' Makes the event happen by updating island species matrix and species IDs.
#' What event happens is determined by the sampling in the algorithm.
#' !!!!!!!!!!!!THIS DOCUMENTATION MUST BE CONFIRMED!!!!!!!!!!!!
#' 
#' 
#' @param timeval current time of simulation
#' @param possible_event numeric indicating what event will happen.
#' @param maxspecID current species IDs
#' @param mainland_spec number of mainland species
#' @param island_spec matrix with species on island (state of system at each time point)
DAISIE_sim_update_state <- function(timeval, possible_event,maxspecID,mainland_spec,island_spec)
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
