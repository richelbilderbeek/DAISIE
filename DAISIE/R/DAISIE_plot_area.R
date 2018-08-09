#### THESE FUNCTIONS NEED MORE TESTING


#' Plots island area function through time
#'
#' @param totaltime total time of simulation
#' @param Apars a numeric vector:
#' \itemize{
#'   \item{[1]: maximum area}
#'   \item{[2]: value from 0 to 1 indicating where in the island's history the 
#'   peak area is achieved}
#'   \item{[3]: sharpness of peak}
#'   \item{[4]: total island age}
#' }
#' @param island_ontogeny a string describing the type of island ontogeny. Can be \code{NULL},
#' \code{"quadratic"} for a beta function describing area through time,
#'  or \code{"linear"} for a linear function
#' @param resolution numeric indicating resolution of plot. Should be < 0.
#' @family rates calculation
#'
#' @return a plot with the area size through time
#' @export
DAISIE_plot_area <- function(totaltime,
                             Apars,
                             island_ontogeny = "quadratic",
                             resolution
) {
  testit::assert(are_area_params(Apars))
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  
  axis <- seq(0, totaltime, by = resolution)
  area <- c()
  for (i in seq_along(axis)) {
    testit::assert(are_area_params(Apars))
    area[i] <- DAISIE::island_area(timeval = axis[i],
                                totaltime = totaltime,
                                Apars = Apars,
                                island_ontogeny = island_ontogeny
                                )
    
  }
  island_area_time <- data.frame(Area = area, Time = axis, Totaltime = totaltime)
  
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Area <- NULL; rm(Area) # nolint, fixes warning: no visible binding for global variable
  ggplot2::ggplot(data = island_area_time, ggplot2::aes(x = Time, y = Area)) +
    ggplot2::geom_line(size = 1.5)
  invisible(island_area_time)
}

#' Plots extinction rate function through time
#'
#' @param island_area_time something
#' @param totaltime something
#' @param K something
#' @param Apars something
#' @param Epars something
#' @param island_ontogeny something
#' @param removed_timepoints something
#'
#' @return extinction rate through time plot
#' @export
DAISIE_plot_extinction <- function(island_area_time,
                                   totaltime,
                                   K, 
                                   Apars, 
                                   Epars, 
                                   island_ontogeny = "quadratic", 
                                   removed_timepoints) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  ext_rate <- c()
  for (i in seq_along(island_area_time$Time)) {
    ext_rate[i] <- DAISIE::get_ext_rate(timeval = island_area_time$Time[i],
                                   totaltime = totaltime,
                                   Apars = Apars,
                                   Epars = Epars,
                                   mu = NA, 
                                   K = K, 
                                   extcutoff = 1100, 
                                   island_spec = matrix(ncol = 1),
                                   island_ontogeny = island_ontogeny
                                   )
  }
  
  ext_rate_time <- data.frame(Extinction = ext_rate[removed_timepoints:length(ext_rate)], Time = island_area_time$Time[removed_timepoints:length(island_area_time$Time)])
  
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Extinction <- NULL; rm(Extinction) # nolint, fixes warning: no visible binding for global variable
  ggplot2::ggplot(data = ext_rate_time, ggplot2::aes(x = Time, y = Extinction)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::ylim(0, 10)
  invisible(ext_rate_time)
}

#' Plot immigration rate through time
#'
#' @param island_area_time Something
#' @param totaltime Something
#' @param K Something
#' @param Apars Something
#' @param gam Something
#' @param mainland_n Something
#' @param island_ontogeny Something
#' @param removed_timepoints Something
#' @param immig_rate Something
#'
#' @return a plot with immigration rate through time
#' @export
#'
DAISIE_plot_immigration <- function(island_area_time,
                                   totaltime,
                                   K, 
                                   Apars, 
                                   gam,
                                   mainland_n,
                                   island_ontogeny = "quadratic", 
                                   removed_timepoints,
                                   immig_rate) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  ext_rate <- c()
  for (i in seq_along(island_area_time$Time)) {
    ext_rate[i] <- get_immig_rate(timeval = island_area_time$Time[i],
                                  totaltime = totaltime,
                                  Apars = Apars,
                                  gam = NA, 
                                  K = K, 
                                  mainland_n = 1000, 
                                  island_spec = matrix(ncol = 1),
                                  island_ontogeny = island_ontogeny
    )
  }
  
  immig_rate_time <- data.frame(Immigration = immig_rate[removed_timepoints:length(immig_rate)], Time = island_area_time$Time[removed_timepoints:length(island_area_time$Time)])
  
  Time <- NULL; rm(Time) # nolint, fixes warning: no visible binding for global variable
  Immigration <- NULL; rm(Immigration) # nolint, fixes warning: no visible binding for global variable
  ggplot2::ggplot(data = immig_rate_time, ggplot2::aes(x = Time, y = Immigration)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::ylim(0, 10)
  invisible(immig_rate_time)
}

