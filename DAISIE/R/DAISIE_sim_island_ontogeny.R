#' Simulate islands with given parameters.
#' @description This function simulates islands with given cladogenesis,
#'  extinction, Kprime, immigration and anagenesis parameters.
#'   If a single parameter set is provided (5 parameters) it simulates islands
#'    where all species have the same macro-evolutionary process.
#'     If two paramater sets (10 parameters) are provided, it simulates islands
#'      where two different macro-evolutionary processes operate, one applying
#'       to type 1 species and other to type 2 species.
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
DAISIE_sim_island_ontogeny <- function(
  time,
  M,
  pars,
  replicates,
  divdepmodel = 'CS',
  prop_type2_pool = NA,
  replicates_apply_type2 = TRUE,
  sample_freq = 25,
  plot_sims = TRUE,
  island_ontogeny = NULL, # NULL = no effect; "quadratic" = quadratic function; "linear" = linear function
  Apars = NULL,
  Epars = NULL,
  verbose = TRUE,
  ...
) {
  totaltime <- time
  island_replicates  = list()
  
  if(divdepmodel =='IW')
  {
    if(length(pars) > 5)
    {
       stop('Island-wide carrying capacity model not yet implemented for two types of mainland species')
    }
    
    for(rep in 1:replicates)
    {
      island_replicates[[rep]] <- DAISIE_sim_core_island_ontogeny(
        time = totaltime,
        mainland_n = M,
        pars = pars,
        island_ontogeny = island_ontogeny,
        Apars = Apars,
        Epars = Epars
      )
      if (verbose == TRUE) {
        print(paste("Island replicate ",rep,sep = ""))
      }
    } 
    island_replicates = DAISIE_format_IW(island_replicates = island_replicates,
                                         time = totaltime,M = M,sample_freq = sample_freq)
  }
      
  if(divdepmodel == 'CS')
  {
    if(length(pars) == 5)
    { 
      for(rep in 1:replicates)
      {
        island_replicates[[rep]] = list() 
        # Run each clade seperately
        full_list = list()
        for(m_spec in 1:M) 
        { 	
          full_list[[m_spec]] <- DAISIE_sim_core_island_ontogeny(
            time = totaltime,
            mainland_n = 1,
            pars = pars,
            island_ontogeny = island_ontogeny,
            Apars = Apars,
            Epars = Epars
          )
          # print(full_list)
        }
        
        island_replicates[[rep]] = full_list
        if (verbose == TRUE) {
          print(paste("Island replicate ",rep,sep = ""))
        }
      } 
    }
    
    if(length(pars) == 10)
    {
      if(is.na(prop_type2_pool))
      {
        stop('prop_type2_pool (fraction of mainland species that belongs to the second subset of species) must be specified when running model with two species types')
      }
      
      if(replicates_apply_type2 == TRUE)
      {
        island_replicates = DAISIE_sim_min_type2(time = totaltime,
                                                 M = M,
                                                 pars = pars,
                                                 replicates = replicates, 
                                                 prop_type2_pool = prop_type2_pool)
      } else
      {
        for(rep in 1:replicates)
        {
          pool2 = DDD::roundn(M * prop_type2_pool)
          pool1 = M - pool2
          
          lac_1 = pars[1]
          mu_1 = pars[2]
          K_1 = pars[3]
          gam_1 = pars[4]
          laa_1 = pars[5]
          
          lac_2 = pars[6]
          mu_2 = pars[7]
          K_2 = pars[8]
          gam_2 = pars[9]
          laa_2 = pars[10]
          
          full_list = list()
          
          #### species of pool1
          for(m_spec in 1:pool1) 
          { 	
            full_list[[m_spec]] = DAISIE_sim_core_island_ontogeny(time = totaltime,mainland_n = 1,pars = c(lac_1,mu_1,K_1,gam_1,laa_1))
            full_list[[m_spec]]$type1or2  = 1
          }
          
          #### species of pool2
          for(m_spec in (pool1 + 1):(pool1 + pool2)) 
          { 	
            full_list[[m_spec]] = DAISIE_sim_core_island_ontogeny(time = totaltime,mainland_n = 1,pars = c(lac_2,mu_2,K_2,gam_2,laa_2))
            full_list[[m_spec]]$type1or2 = 2
          }
          island_replicates[[rep]] = full_list
          if (verbose == TRUE) {
            print(paste("Island replicate ",rep,sep = ""))
          }
        }
      }
    }

    island_replicates = DAISIE_format_CS(island_replicates = island_replicates,
                                         time = totaltime,
                                         M = M,
                                         sample_freq = sample_freq)
  }
  if(plot_sims == TRUE)
  { 
    DAISIE_plot_sims(island_replicates)
  }
  return(island_replicates)
}
