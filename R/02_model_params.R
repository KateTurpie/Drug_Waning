# ---------------------------------------------------------------------------- #

#### ---- Parameter Functions ---- ####
#
# Purpose : 
#
# By : Kate Turpie
# Updated : 2026-05-16

# ---------------------------------------------------------------------------- #

# ---- Set Model Parameter ----

input_contacts <- function(country = "United Kingdom",
                           age_limits = c(0, 1, 2, 3, 4, 5, 18, 26, 40, 60, 75),
                           age_groups = c("0-1", "1", "2", "3", "4", "5-17", "18-25", "26-39", "40-59", "60-74", "75+")) {
  # Check
  if (length(age_limits) != (length(age_groups))) {
    stop("age_limits and age_groups must the be the same length.")
  }
  
  # Load contact matrix
  data(polymod)
  contacts <- contact_matrix(
    polymod,
    countries = country,
    age.limits = age_limits,
    symmetric = TRUE,
    return.demography = TRUE
  )
  
  age_number = length(age_groups)
  
  # Output list with contact matrix, demography stats, number of age groups in analysis, and age group labels
  out <- list(
    contacts = contacts$matrix,
    demog = contacts$demography$proportion,
    n_age_groups = age_number,
    age_groups = age_groups
  )
  
  return(out)
}

input_infection <- function(R0 = 3,
                            gamma = 7,
                            delta = 275,
                            immunity = 0.5,
                            contacts = NULL) {
  
  # Checks 
  if (is.null(contacts) == TRUE) {
    stop("Please input contact matrix from function : input_contacts()")
  }
  if (gamma <= 0 || R0 <= 0) {
    stop("gamma and R) must be greater than 0.")
  }
  
  # Days to rates
  if (delta != 0) {delta = 1 / delta}
  if (gamma != 0) {gamma = 1 / gamma}
  
  # Build Next Generation Matrix and get Dominant Eigenvalue
  NGM <- (t(contacts)) * (1 / gamma)
  max_eig <- max(Re(eigen(NGM)$values))
  
  # Calculate Beta
  beta = R0 / max_eig
  
  out <- list(
    beta = beta,
    gamma = gamma,
    delta = delta,
    immunity = immunity
  )
  
  return(out)
  
}

input_prophylactic <- function(efficacy = 0.9, 
                               duration = 90, 
                               protection = c(0,0,0,0,0,0,0,0,0,0,0)) {
  
  # Days to rates
  if (duration != 0) {duration_rate = 1 / duration} else {duration_rate = 0} 
  
  out <- list(chi = efficacy, omega = duration_rate, rho = protection)
  
  return(out)
}

input_epi <- function(p_ma = c(0.2,0.2,0.2,0.2,0.2,0.2,1/6,1/6,1/6,1/6,1/6),
                      alpha = 0, 
                      offset = 0) {
  out <- list(p_ma = p_ma,
              alpha = alpha,
              offset = offset)
  
  return(out)
}

input_trial <- function(trial_duration = 274, #9 months
                        trial_start_time = 0 ) { #default to beginning of run
  
  tau_end = trial_duration + trial_start_time
  
  if(is.null(trial_start_time) == TRUE) {
    tau_start <- 1
    message(paste0("Trial starting from Day 1 of model run & running until Day ", tau_end))
  }
  else {
    tau_start <- trial_start_time
    message(paste0("Trial starting from Day ", trial_start_time, " of model run & running until Day ", tau_end))
  }
  
  out <- list(tau_end = tau_end, tau_start = tau_start)
  
  return(out)
}

# ---- Set Initial Conditions ----

input_initial_conditions <- function(population = 68350000,
                                     demog,
                                     set.seed_true = TRUE,
                                     seed = 99,
                                     start_P_seed = FALSE,
                                     n_trial_groups = 1,
                                     n_trial_pop = 1000) {
  if (is.null(demog) == TRUE) {
    stop("Please input demog from function : input_contacts()")
  }
  
  #divide up susceptibles by age group
  start_S = demog * population
  
  #remove older adults
  start_R <- rep(0, length(start_S))
  start_R = round(start_S * c(0.0, 0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 0.7, 0.7, 0.6, 0.6))
  start_S = start_S - start_R
  
  # seed infection
  if (set.seed_true == TRUE) {
    set.seed(seed)
  }
  
  start_Iu <- rep(1, length(start_S))
  start_S <- start_S - start_Iu
  
  # seed prophylaxis (if required)
  if (start_P_seed == FALSE) {
    start_P = rep(0, length(start_S))
  }
  else {
    start_P <- round(start_S * c(0.9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    start_S <- start_S - start_P
  }
  
  # initialize other compartments to 0 across all age groups
  
  start_Ip <- rep(0, length(start_S))
  
  return(
    list(
      start_S = start_S,
      start_Iu = start_Iu,
      start_Ip = start_Ip,
      start_P = start_P,
      start_R = start_R,
      start_Sc = n_trial_pop / 2,
      start_Stx = n_trial_pop / 2,
      n_age_groups_trial = n_trial_groups
    )
  )
}

# ---- Set Model Meta Parameters ----

input_model_meta <- function(burn_in = TRUE,
                             burn_n = 200,
                             days_n = 365,
                             stochastic = FALSE,
                             model_n_run = NULL) {
  
  if (burn_in == TRUE) {end_time = burn_n + days_n} else {end_time = days_n}
  
  if(stochastic == FALSE) {n_runs = 1} 
  else if (is.null(model_n_run) == FALSE) {n_runs = model_n_run} 
  else {stop("Enter n_run")}
  
  out <- list(end_time = end_time, 
              n_runs = n_runs,
              burn = burn_n)
  
  return(out)
}

# ---- Set defaults for all params ----

input_all <- function(defaults = TRUE) {
  
  if (defaults == FALSE) {
    stop("Please set parameters manually using the following functions : \n
          input_contacts() : retrieves polymod contact matrix\n
          input_infection() : sets beta, gamma, delta, & omega\n
          input_prophylactic() : sets chi, omega & rho\n
          input_epi() : sets p_ma, alpha, offset\n
          input_parameters() : sets trial duration & start time\n
          input_initial_conditions() : sets starting population in each compartment\n
         input_model_meta() : sets burn in, duration, & number of runs")
  }
  
  c_params <- input_contacts()
  i_params <- input_infection(contacts = c_params$contacts)
  p_params <- input_prophylactic()
  e_params <- input_epi()
  initial_conditions <- input_initial_conditions(demog = c_params$demog)
  meta_params <- input_model_meta()
  t_params <- input_trial(trial_start_time = meta_params$burn)
  
  inputs = c(c_params, i_params, p_params, e_params, t_params, initial_conditions)
  
  out <- list(inputs = inputs, meta_params = meta_params)
  
  return(out)
}

# ---- Model Calibration ----
