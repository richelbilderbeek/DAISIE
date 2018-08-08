are_area_params <- function(area_params) {
  if (!"max_area" %in% names(area_params)) return(FALSE)
  if (!"proportional_peak_t" %in% names(area_params)) return(FALSE)
  if (!"peak_sharpness" %in% names(area_params)) return(FALSE)
  if (!"total_island_age" %in% names(area_params)) return(FALSE)
  if (area_params$max_area < 0.0) return(FALSE)
  if (area_params$proportional_peak_t <= 0.0) return(FALSE)
  if (area_params$proportional_peak_t >= 1.0) return(FALSE)
  if (area_params$peak_sharpness <= 0) return(FALSE)
  if (area_params$total_island_age <= 0.0) return(FALSE)
  TRUE
}
create_area_params <- function(max_area,
                               proportional_peak_t,
                               peak_sharpness,
                               total_island_age) {
  testit::assert(max_area > 0.0)
  testit::assert(proportional_peak_t >= 0.0)
  testit::assert(proportional_peak_t <= 1.0)
  testit::assert(peak_sharpness >= 0)
  testit::assert(total_island_age >= 0.0)
  list(max_area = max_area, 
       proportional_peak_t = proportional_peak_t,
       peak_sharpness = peak_sharpness,
       total_island_age = total_island_age)
}