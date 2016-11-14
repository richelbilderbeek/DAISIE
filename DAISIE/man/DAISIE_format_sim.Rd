\name{DAISIE_format_sim}
\alias{DAISIE_format_sim}
\title{Converting simulation output from RAW to DAISIE format}
\description{Formats islands in RAW format (simulated in DAISIE_sim with format=FALSE) to a DAISIE format}
\usage{
DAISIE_format_sim(
   island_replicates,
   time,
   M,
   sample_freq
   )
}

\arguments{
\item{island_replicates}{Island replicates in RAW format produced in DAISIE_sim with format=FALSE option}
\item{time}{ Length of time the islands were simulated with}
\item{M}{The size of the mainland pool, i.e the number of species that can
 potentially colonize the island }
\item{sample_freq}{Specifies the number of units time should be divided by for 
plotting purposes. Larger values will lead to 
plots with higher definition, but will also run slower.}
}
\details{Returns R list object that contains the simulated islands in DAISIE format that can be run on other
 DAISIE functions. }
\value{

Each simulated dataset is an element of the list, which can be called using [[x]]. For example if the object is called island_replicates, the 1st replicate can be called using island_replicates[[1]]
Each of the island replicates is a list in itself. The first (e.g. island_replicates[[x]][[1]]) element of that list has the following components: \cr
  \code{$island_age} - the island age \cr
  Then, depending on whether a distinction between types is made, we have:\cr
  \code{$not_present} - the number of mainland lineages that are not present on the island \cr  
  or:\cr
  \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr  
  \code{$not_present_type2} - the number of mainland lineages of type 2 that are not present on the island \cr 
  \code{$stt_all} - STT table for all species on the island (nI - number of non-endemic species; nA - number of anagenetic species, nC - number of cladogenetic species, present - number of independent colonisations present )\cr
  \code{$stt_stt_type1} - STT table for type 1 species on the island - only if 2 types of species were simulated (nI - number of non-endemic species; nA - number of anagenetic species, nC - number of cladogenetic species, present - number of independent colonisations present )\cr   
  \code{$stt_stt_type2} - STT table for type 2 species on the island - only if 2 types of species were simulated (nI - number of non-endemic species; nA - number of anagenetic species, nC - number of cladogenetic species, present - number of independent colonisations present )\cr   

  The subsequent elements of the list each contain information on a single colonist lineage on the island and has 4 components:\cr
  \code{$branching_times} - island age and stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species these should be island age and branching times of the radiation including the stem age of the radiation.\cr                          
  \code{$stac} - the status of the colonist \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr
  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$type_1or2} - whether the colonist belongs to type 1 or type 2 \cr
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852.}
\author{
Luis Valente
}
\seealso{
\code{\link{DAISIE_plot_sims}}
\code{\link{DAISIE_sim}}
}
\examples{ 

data(islands_10reps_RAW)
island_replicates = DAISIE_format_sim(
    island_replicates = islands_10reps_RAW,
    time = 4,
    M = 1000,
    sample_freq = 25
    )

}
\keyword{ models }