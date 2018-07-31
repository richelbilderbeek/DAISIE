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
                             resolution) {
  axis <- seq(0, totaltime, by = resolution)
  area <- c()
  for (i in seq_along(axis)) {
    area[i] <- DAISIE::island_area(timeval = axis[i],
                                totaltime = totaltime,
                                Apars = Apars,
                                island_ontogeny = island_ontogeny
                                )
    
  }
  island_area_time <- data.frame(Area = area, Time = axis, Totaltime = totaltime)
  ggplot2::ggplot(data = island_area_time, ggplot2::aes(x = Time, y = Area)) +
    ggplot2::geom_line(size = 1.5)
  invisible(island_area_time)
}

#' Plots extinction rate function through time
#'
#' @param island_area_time something
#' @param K something
#' @param Apars something
#' @param Epars something
#' @param island_ontogeny something
#' @param removed_timepoints something
#'
#' @return extinction rate through time plot
#' @export
DAISIE_plot_extinction <- function(island_area_time, 
                                   K, 
                                   Apars, 
                                   Epars, 
                                   island_ontogeny = "quadratic", 
                                   removed_timepoints) {
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
  # str(ext_rate_time)
  ggplot2::ggplot(data = ext_rate_time, ggplot2::aes(x = Time, y = Extinction)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::ylim(0, 10)
  invisible(ext_rate_time)
}
