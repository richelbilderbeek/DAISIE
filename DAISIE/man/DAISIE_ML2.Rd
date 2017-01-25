\name{DAISIE_ML2}
\alias{DAISIE_ML2}
\title{Maximization of the loglikelihood under the DAISIE model, when there are multiple island systems}
\description{
This function computes the maximum likelihood estimates of the parameters of the
DAISIE model for data from lineages colonizing several island systems. It also outputs the
corresponding loglikelihood that can be used in model comparisons.
}
\usage{
DAISIE_ML2(
   datalist,
   idparsmat,
   initparsopt,
   idparsopt,
   parsfix,
   idparsfix,
   res = 100,
   ddmodel = 0,
   cond = 0,
   tol = c(1e-04, 1e-05, 1e-07),
   maxiter = 1000 * round((1.25)^length(idparsopt)),
   methode = "lsodes",
   optimmethod = 'subplex'
   )
}
\arguments{
  \item{datalist}{A list of data objects containing information on colonisation and branching times for each island system. This object can be generated using the DAISIE_dataprep function, which converts a user-specified data table into a data object, but the object can of course also be entered directly. Each data object is an R list object with the following elements.\cr
  The first element of the list has two three components: \cr \cr
  \code{$island_age} - the island age \cr
  Then, depending on whether a distinction between types is made, we have:\cr
  \code{$not_present} - the number of mainland lineages that are not present on the island \cr  
  or:\cr
  \code{$not_present_type1} - the number of mainland lineages of type 1 that are not present on the island \cr  
  \code{$not_present_type2} - the number of mainland lineages of type 2 that are not present on the island \cr \cr  
  The remaining elements of the list each contains information on a single colonist lineage on the island and has 5 components:\cr \cr
  \code{$colonist_name} - the name of the species or clade that colonized the island \cr  
  \code{$branching_times} - island age and stem age of the population/species in the case of Non-endemic, Non-endemic_MaxAge and Endemic anagenetic species. For cladogenetic species these should be island age and branching times of the radiation including the stem age of the radiation.\cr                          
  \code{$stac} - the status of the colonist \cr \cr
  * Non_endemic_MaxAge: 1 \cr
  * Endemic: 2 \cr
  * Endemic&Non_Endemic: 3 \cr
  * Non_endemic: 4 \cr \cr
  \code{$missing_species} - number of island species that were not sampled for particular clade (only applicable for endemic clades) \cr
  \code{$type1or2} - whether the colonist belongs to type 1 or type 2. Currently only one type is implemented. \cr
  }
  \item{idparsmat}{ Matrix containing the ids of the parameters, linking them to initparsopt and parsfix. Per island system we use the following order: \cr \cr
  * lac = (initial) cladogenesis rate \cr
  * mu = extinction rate \cr
  * K = maximum number of species possible in the clade \cr
  * gam = (initial) immigration rate \cr
  * laa = (initial) anagenesis rate \cr
  Example: idparsmat = rbind(c(1,2,3,4,5),c(1,2,3,6,7)) has different rates of immigration and anagenesis for the two islands.} 
  \item{initparsopt}{The initial values of the parameters that must be optimized}
  \item{idparsopt}{
   The ids of the parameters that must be optimized.}
  \item{idparsfix}{The ids of the parameters that should not be optimized.}
  \item{parsfix}{The values of the parameters that should not be optimized }
  \item{res}{Sets the maximum number of species for which a probability must be computed, must be larger than the size of the largest clade }
  \item{ddmodel}{Sets the model of diversity-dependence: \cr \cr
    ddmodel = 0 : no diversity dependence \cr
    ddmodel = 1 : linear dependence in speciation rate \cr
    ddmodel = 11: linear dependence in speciation rate and in immigration rate \cr
    ddmodel = 2 : exponential dependence in speciation rate\cr
    ddmodel = 21: exponential dependence in speciation rate and in immigration rate\cr
  }
  \item{cond}{
    cond = 0 : conditioning on island age \cr
    cond = 1 : conditioning on island age and non-extinction of the island biota \cr
  }
  \item{tol}{Sets the tolerances in the optimization. Consists of: \cr
   reltolx = relative tolerance of parameter values in optimization \cr
   reltolf = relative tolerance of function value in optimization \cr
   abstolx = absolute tolerance of parameter values in optimization
  }
  \item{maxiter}{Sets the maximum number of iterations in the optimization}
  \item{methode}{Method of the ODE-solver. See package deSolve for details. Default is "lsodes"}
  \item{optimmethod}{Method used in likelihood optimization. Default is "subplex" (see subplex package). Alternative is 'simplex' which was the method in previous versions.}  
}
\details{The result of sort(c(idparsopt, idparsfix, idparsnoshift)) should be identical to c(1:10). If not, an error is reported that the input is incoherent. The same happens when the length of initparsopt is different from the length of idparsopt, and the length of parsfix is different from the length of idparsfix.\cr
Including the 11th parameter (p_f) in either idparsopt or idparsfix (and therefore initparsopt or parsfix) is optional. If this parameter is not specified, then the information in the data is used, otherwise the information in the data is overruled.
}
\value{The output is a dataframe containing estimated parameters and maximum loglikelihood.
  \item{lambda_c}{ gives the maximum likelihood estimate of lambda^c, the rate of cladogenesis for each island system}
  \item{mu}{ gives the maximum likelihood estimate of mu, the extinction rate for each island system}
  \item{K}{ gives the maximum likelihood estimate of K, the carrying-capacity for each island system}
  \item{gamma}{ gives the maximum likelihood estimate of gamma, the immigration rate for each island system}
  \item{lambda_a}{ gives the maximum likelihood estimate of lambda^a, the rate of anagenesis for each island system}
  \item{loglik}{ gives the maximum loglikelihood}
  \item{df}{ gives the number of estimated parameters, i.e. degrees of feedom}
  \item{conv}{ gives a message on convergence of optimization; conv = 0 means convergence}
}
\references{
Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852. <DOI:10.1111/ele.12461>.}
\author{
Rampal S. Etienne
}
\seealso{
 \code{\link{DAISIE_ML}}
 \code{\link{DAISIE_loglik_all}},
 \code{\link{DAISIE_sim}}
}
\examples{
cat("
# When both islands have the same rates except for immigration and anagensis rate,
# and we want to optimize all 5 parameters, we use:

data(Galapagos_datalist)
DAISIE_ML2(
   datalist = list(Galapagos_datalist,Galapagos_datalist)
   initparsopt = c(2.5,2.7,20,0.009,1.01),
   idparsmat = rbind(1:5,c(1:3,6,7))
   idparsopt = 1:5,
   parsfix = NULL,
   idparsfix = NULL
)
")
}
\keyword{models}