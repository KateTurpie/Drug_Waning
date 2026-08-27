# ---------------------------------------------------------------------------- #

#### ---- Calibrate Model ---- ####
#
# Purpose : Calibrate alpha, offset and immunity
#
# By : Kate Turpie
# Updated : 2026-05-16

# ---------------------------------------------------------------------------- #

burn_n <- 200 #idk why i chose this, just spiritually decided lol

calibration_grid <- expand.grid(
  alpha = c(0.20, 0.35, 0.50),
  offset = c(30, 60, 90, 120, 150),
  immunity = c(0.25, 0.50, 0.75)
)

# Total annual RSV infections among infants aged under 1
infant_case_target <- 469871

# Day 0 is 1 September
peak_target <- 105   # approximately mid-December
width_target <- 135  # days above 10% of peak infant incidence

base <- input_all()
summaries <- vector("list", nrow(calibration_grid))
curves <- vector("list", nrow(calibration_grid))

for (i in seq_len(nrow(calibration_grid))) {
  
  inputs <- base$inputs
  inputs$alpha <- calibration_grid$alpha[i]
  inputs$offset <- calibration_grid$offset[i]
  inputs$immunity <- calibration_grid$immunity[i]
  
  # Switch prophylaxis off while generating the background epidemic
  inputs$rho[] <- 0
  inputs$start_S <- inputs$start_S + inputs$start_P
  inputs$start_P[] <- 0
  inputs$tau_start <- burn_n + 366
  inputs$tau_end <- burn_n + 367
  
  system <- dust2::dust_system_create(
    model_sirs,
    inputs,
    n_particles = 1,
    deterministic = TRUE
  )
  
  dust2::dust_system_set_state_initial(system)
  
  raw <- dust2::dust_system_simulate(
    system,
    0:(burn_n + 365)
  )
  
  state_names <- names(unlist(dust2::dust_unpack_index(system)))
  
  if (nrow(raw) == length(state_names)) {
    rownames(raw) <- state_names
  }
  
  results <- output_tidyup(raw)
  
  # Find and correctly order cumulative total-infection columns
  infection_cols <- grep(
    "^TotalCases_pop[0-9]+$",
    names(results),
    value = TRUE
  )
  
  infection_cols <- infection_cols[
    order(as.integer(sub("TotalCases_pop", "", infection_cols)))
  ]
  
  # Final model year: t = 365 to t = 730
  year <- results[
    results$t >= burn_n & results$t <= burn_n + 365,
  ]
  
  if (nrow(year) != 366) {
    stop("Expected 366 daily rows in the final model year.")
  }
  
  infant_cumulative <- year[[infection_cols[1]]]
  daily_infant_infections <- diff(infant_cumulative)
  
  infant_cases <-
    infant_cumulative[nrow(year)] -
    infant_cumulative[1]
  
  peak_day <- which.max(daily_infant_infections)
  
  season_width <- sum(
    daily_infant_infections >
      0.10 * max(daily_infant_infections)
  )
  
  # Lower is better: infant burden, peak timing and epidemic duration
  score <-
    log(infant_cases / infant_case_target)^2 +
    ((peak_day - peak_target) / 30)^2 + # measured error in months
    ((season_width - width_target) / 45)^2 # measured error by 6 weeks 
                                          #  (1/3 of ideal season length)
  
  summaries[[i]] <- data.frame(
    calibration_grid[i, ],
    score = score,
    peak_day = peak_day,
    season_width = season_width,
    cases_0_1 = infant_cases
  )
  
  curves[[i]] <- data.frame(
    day = 1:365,
    daily_infections = daily_infant_infections
  )
}

# Rank candidates and retain the best infant epidemic
calibration_results <- do.call(rbind, summaries)

best_row <- which.min(calibration_results$score)
best_params <- calibration_results[best_row, ]
best_curve <- curves[[best_row]]

calibration_results <- calibration_results[
  order(calibration_results$score),
]

best_params
head(calibration_results, 10)

plot(
  best_curve$day,
  best_curve$daily_infections,
  type = "l",
  lwd = 2,
  xlab = "Day since 1 September",
  ylab = "Daily infant infections",
  main = "Best UK-like infant RSV epidemic"
)
