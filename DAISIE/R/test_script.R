n_mainland_species <- 1000
island_age <- 10
clado_rate <- 0.0001 # cladogenesis rate
ext_rate <- 0.5 # extinction rate (not used)
clade_carr_cap <- 0.5  # clade-level carrying capacity
imm_rate <- 0.001 # immigration rate
ana_rate <- 0.2 # anagenesis rate
max_area <- 1000
peak_time <- 0.1
sharpness <- 1
total_island_age <- 12
mu_min <- 0.05
mu_max <- 100
island_ontogeny <- "quadratic"
extcutoff <- 1000
keep_final_state = TRUE


out <- DAISIE_sim(
  time = island_age, 
  M = n_mainland_species, 
  pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
  replicates = 1, 
  island_ontogeny = island_ontogeny,
  Apars = create_area_params(max_area, peak_time, sharpness, total_island_age),
  Epars = c(mu_min, mu_max),
  extcutoff = extcutoff,
  plot_sims = TRUE,
  verbose = FALSE,
  keep_final_state = keep_final_state
)



n_mainland_species <- 1000
island_age <- 4
clado_rate <- 0.7 # cladogenesis rate
ext_rate <- 0.5 # extinction rate (not used)
clade_carr_cap <- 20  # clade-level carrying capacity
imm_rate <- 0.001 # immigration rate
ana_rate <- 0.4 # anagenesis rate
max_area <- 4000
peak_time <- 0.1
sharpness <- 1
total_island_age <- 10
mu_min <- 0.05
mu_max <- 100
island_ontogeny <- "quadratic"
extcutoff <- 1000
keep_final_state = TRUE




midway_out <- DAISIE_sim(
  time = 11, 
  M = n_mainland_species, 
  pars = c(clado_rate, ext_rate, clade_carr_cap, imm_rate, ana_rate),
  replicates = 1, 
  island_ontogeny = NULL,
  extcutoff = extcutoff,
  plot_sims = FALSE,
  verbose = FALSE,
  stored_data = out
)
