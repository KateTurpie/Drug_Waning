# ---------------------------------------------------------------------------- #

#### ---- Parameter Functions ---- ####
#
# Purpose : 
#
# By : Kate Turpie
# Updated : 2026-05-16

# ---------------------------------------------------------------------------- #

# ---- Run Model Functions ----

run_model <- function(model,
                      params,
                      dar,
                      deterministic = TRUE,
                      param_one = NULL,
                      one_values = NULL,
                      param_two = NULL,
                      two_values = NULL,
                      all_combos = NULL,
                      quiet = TRUE) {
  
  # ---- Running Functions ----
  
  ## ---- Run Once ----
  run_once <- function(model, params, deterministic) {
    simulator <- dust2::dust_system_create(model,
                                           params,
                                           n_particles=dar$n_runs,
                                           deterministic = deterministic)
    dust2::dust_system_set_state_initial(simulator)
    run_output <- dust2::dust_system_simulate(simulator, 1:dar$end_time)
    rownames(run_output) <- names(unlist(dust2::dust_unpack_index(simulator)))
    run_output <- output_tidyup(run_output)
    return(run_output)
  }
  ## ---- Vary One Parameter ----
  run_one_param <- function(param_name, values) {
    grid <- data.frame(value = values)
    results <- lapply(values, function(v) {
      p <- params
      p[param_name] <- v
      run_once(model, p, deterministic)
    })
    grid$result <- results
    grid %>% unnest(result)
  }
  ## ---- Vary Two Parameters ----
  run_two_params <- function(p1, v1, p2, v2, all_combos) {
    
    if (p2 == "omega") {
      v2 <- 1/v2 
    }
    
    grid <- if (all_combos) expand_grid(one = v1, two = v2)
    else data.frame(one = v1, two = v2)
    results <- pmap(grid, function(one, two){
      p <- params
      p[p1] <- one
      p[p2] <- two
      message("Running the model")
      message(paste0("Param One: ", one, " && Param Two: ", two))
      run_once(model, p, deterministic)
    })
    grid$results <- results
    grid %>% unnest(results)
  }
  
  # ---- Run Model ----
  
  # Run the model either:
  # (1) with a single set of parameters,
  # (2) varying a single parameter,
  # (3) varying two parameters
  
  # -- 1 --
  if(is.null(param_one) && is.null(param_two)) {
    if (!quiet){
      mode <- ifelse(deterministic, "deterministic", "stochastic")
      message(paste("Running", mode, "simulation..."))}
    result <- run_once(model, params, deterministic)
    if (!quiet){
      message("Done 🎉")}
    return(result)
  }
  # -- 2 --
  if(!is.null(param_one) && is.null(param_two)) {
    result <- run_one_param(param_one,one_values)
    return(result)
  }
  if(is.null(param_one) && !is.null(param_two)) {
    message(paste("🪄Please use param_one to run single-variations✨"))
    message(paste("It's wierd you used param_two :D"))
  }
  # -- 3 --
  if(!is.null(param_one) && !is.null(param_two)){
    result <- run_two_params(param_one, one_values, param_two, two_values, all_combos)
    return(result)
  }
}

run_model_from_contours <- function(values, save = FALSE) {
  
  update_chi_omega <- output_contours(values, save)
  model_population_analysis <- list()
  
  for (i in 1:length(update_chi_omega)) {
    chi_values <- update_chi_omega[[i]]$chi
    omega_values <- update_chi_omega[[i]]$omega
    
    message(paste0("Running model for ", names(update_chi_omega)[i], " population 🪼"))
    model_population_analysis[[names(update_chi_omega)[i]]] <- suppressMessages(
      run_model(
        model = model_sirs,
        params = params$inputs,
        dar = params$meta_params,
        param_one = "chi", 
        one_values = chi_values,
        param_two = "omega", 
        two_values = omega_values, 
        all_combos = FALSE
      ))
  }
  
  message("ALL DONE WOOHOO 🎉")
  return(model_population_analysis)
  
}


# ---- Main Figures & Analysis ----

# Figure One : Trial Efficacy Heat Map
#' @input model_run : multiple model trial runs
#' @input trial_len : length of trial (days)

figure_1 <- function(model_run, trial_len, trial_start = 0) {
  
  # Trial length label and absolute model endpoint
  len_label <- paste0(
    round(trial_len / 30, 1),
    " months"
  )
  
  trial_endpoint <- trial_start + trial_len
  
  # Prepare results
  runs <- model_run %>%
    mutate(
      # 'one' contains initial efficacy
      initial_efficacy = one,
      
      # 'two' contains the waning rate omega
      duration_protection = 1/two
      ) %>%
    filter(
      dplyr::near(t, trial_endpoint),
      is.finite(duration_protection)
    ) %>%
    select(
      t,
      initial_efficacy,
      duration_protection,
      efficacy_measured
    )
  
  # Give an informative error if the endpoint was not simulated
  if (nrow(runs) == 0) {
    stop(
      "No model results found at t = ",
      trial_endpoint,
      ". Model output ranges from t = ",
      min(model_run$t),
      " to t = ",
      max(model_run$t),
      "."
    )
  }
  
  # Create graph
  graph <- ggplot(
    runs,
    aes(
      x = duration_protection,
      y = initial_efficacy,
      fill = efficacy_measured
    )
  ) +
    geom_tile() +
    
    metR::geom_contour2(
      aes(
        z = efficacy_measured,
        label = after_stat(
          scales::percent(level, accuracy = 1)
        )
      ),
      breaks = c(0.25, 0.50, 0.75, 0.90), #(uncomment to get points for fig2+)
      label.placer =
        metR::label_placer_fraction(0.35),
      colour = "white",
      skip = 0
    ) +
    
    labs(
      title = paste(
        "RCT Duration:",
        len_label
      ),
      x = "Average Duration of Protection (Days)",
      y = "Initial Efficacy",
      fill = "RCT-Measured Efficacy"
    ) +
    
    scale_y_continuous(
      labels = scales::label_percent(),
      limits = c(0.15, 1.05),
      expand = c(0, 0)
    ) +
    
    scale_x_continuous(
      limits = c(30, 1081),
      breaks = seq(
        30,
        1080,
        by = 150
      ),
      expand = expansion(
        mult = c(0, 0.019)
      )
    ) +
    
    scale_fill_gradient2(
      low = "mistyrose2",
      mid = "#EA4E82",
      high = "#8B1134",
      midpoint = 0.50,
      limits = c(0, 1),
      labels = scales::label_percent(accuracy = 1),
      name = "RCT-Measured Efficacy",
      oob = scales::squish, 
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = unit(9, "cm"),
        barheight = unit(0.45, "cm")
      )
    ) +
    
    theme_clean() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom",
      panel.spacing = unit(
        1.25,
        "lines"
      ), 
      panel.border = element_blank(),
      plot.background = element_rect(colour = NA, fill = NA),
      panel.background = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      legend.box.background = element_blank()
    )
  
  print(graph)
  
  return(graph)
}

# Figure Two : Efficacy Decay Curves
#' @input model_run : multiple model trial runs
#' @input heatmap : result from figure_1()
#' @input multi : boolean, TRUE for supplemental figure with multiple trial lengths

figure_2 <- function(model_run, heatmap, multi = FALSE){
  
  # -- Extract Heat Map Lines --
  plot_pawprint <- ggplot_build(heatmap)
  contour_data <- plot_pawprint$data[[2]]
  
  # -- Calculate Efficacy Decay Over Time --
  t <- seq(0, 365, by = 1)  # adjust time range as needed
  
  efficacy_decay <- contour_data %>%
    select(label, order, x, y)%>%
    mutate(x = 1/x) %>%
    mutate(id = row_number()) %>%        
    tidyr::crossing(time = t) %>%
    mutate(true_efficacy = y * exp(-x * time)) %>%  
    mutate(VE = as.numeric(order))
  
  # -- Build Line Graphs -- 
  graph <- ggplot(efficacy_decay, aes(
    x = time,
    y = true_efficacy,
    group = id,
    color = 1/x
  )) +
    geom_line(alpha = 0.65) +
    scale_color_gradient2(
      midpoint = 500,
      high = "#8B1134",
      mid = "#EA4E82",
      low = "thistle1", 
      limits = c(0, 1080),
      breaks = c(250, 500, 750, 1000),
      labels = c("250", "500", "750", "1000"),
      name = "Duration of Protection",
      oob = scales::squish,
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(6, "cm"),
        barheight = unit(0.45, "cm"))
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    geom_vline(xintercept = 365, linetype = "dashed") +
    annotate("text",
             x = 365, y = 0.95,
             label = "365 days",
             angle = 90,
             vjust = -0.5,
             hjust = 0.75,
             size = 3
    ) +
    labs(
      x = "Days Since Administration",
      y = "Efficacy",
      color = "Duration of Protection"
    ) +
    theme_clean() +
    theme(
      legend.position = "bottom",
      legend.box = element_blank(),
      legend.box.margin = margin(5, 5, 5, 5),
      legend.margin = margin(10, 15, 10, 10),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"), 
      panel.spacing = unit(
        1.25,
        "lines"
      ), 
      panel.border = element_blank(),
      plot.background = element_rect(colour = NA, fill = NA),
      panel.background = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      legend.box.background = element_blank()
    )
  
  # Facet Wrap: Multi = TRUE for supplemental figure to include 
  if (multi) {
    graph <- graph + facet_grid(trial_len ~ VE, labeller = labeller(
      VE = function(x) paste0("RCT VE: ", as.numeric(x) * 100, "%"),
      trial_len = function(x) paste0(as.numeric(x) / 30, " Months")))
  } else {
    graph <- graph + facet_wrap(~ VE, nrow = 2, labeller = labeller(
      VE = function(x) paste0("RCT-measured VE: ", 
                              as.numeric(x) * 100, "%"))) +
      geom_vline(xintercept = 150, linetype = "dashed") +
      annotate("text",
               x = 150, y = 0.95,
               label = "150 days",
               angle = 90,
               vjust = -0.5,
               hjust = 0.75,
               size = 3)
  }
  print(graph)
  return(list(graph = graph, values = contour_data))
}

# Figure Three : Population-level (Infant) Impact of RCT-measured VE
#' @input model_runs : multiple model population runs
#' @input baseline : baseline model population run
#' @input metric : metric to use for calculating population-level impact
#' @input x_breaks : breaks for x-axis
#' @input y_label : label for y-axis
#' @input y_max : maximum value for y-axis

figure_3 <- function(
    model_run,
    baseline = NULL,
    metric = "I_ma1",
    x_breaks = c(90, 180, 270, 360),
    y_label = "Number of Medically-Attended RSV infections in infants under one", 
    y_max = 50000
) {
  
  plot_data <- model_run %>%
    select(t, VE, infant_infections = all_of(metric)) %>%
    group_by(t, VE) %>%
    summarise(
      median = median(infant_infections, na.rm = TRUE),
      lower = min(infant_infections, na.rm = TRUE),
      upper = max(infant_infections, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Flexible colours for however many VEs exist
  ve_colours <- grDevices::colorRampPalette(
    c("#F6A6B2", "#9B8FE7", "#7ED8C6", "#FFBE85", "#D01856")
  )(length(unique(plot_data$VE)))
  
  names(ve_colours) <- (unique(plot_data$VE))
  
  # Base plot
  graph <- ggplot(
    plot_data,
    aes(
      x = t,
      y = median,
      colour = VE,
      fill = VE,
      group = VE
    )
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.2,
      colour = NA
    ) +
    geom_line(linewidth = 0.5) +
    scale_colour_manual(
      values = ve_colours,
      name = "RCT-measured VE"
    ) +
    scale_fill_manual(
      values = ve_colours,
      name = "RCT-measured VE"
    ) +
    guides(fill = "none") +
    scale_x_continuous(
      breaks = x_breaks,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      labels = scales::label_comma(),
      expand = c(0, 0), 
      limits = c(0, y_max)
    ) +
    labs(
      x = "Time (days)",
      y = y_label
    ) +
    theme_clean() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10)
    )
  
  # Add baseline if supplied
  if (!is.null(baseline)) {
    
    baseline_data <- baseline %>%
      select(t, infant_infections = all_of(metric)) %>%
      group_by(t) %>%
      summarise(
        infant_infections = median(infant_infections, na.rm = TRUE),
        .groups = "drop"
      )
    
    graph <- graph +
      geom_line(
        data = baseline_data,
        aes(x = t, y = infant_infections),
        inherit.aes = FALSE,
        colour = "grey30",
        linewidth = 0.5
      )
  }
  
  print(graph)
  return(graph)
  
}

# Figure Four : Range of Impact by RCT-measured VE
#' @input model_population_all : list of model population runs for each trial length
#' @input baseline : baseline model population run
#' @input metric : metric to use for calculating population-level impact
#' @input y_label : label for y-axis

figure_4 <- function(
    model_population_all,
    baseline,
    metric = "I_ma1",
    y_label = "Proportion of medically-attended RSV cases averted"
) {
  
  # ---- Clean baseline metric ----
  if (!metric %in% names(baseline)) {
    stop("Metric column '", metric, "' not found in baseline.")
  }
  
  baseline_total <- baseline %>%
    summarise(
      baseline_cases = sum(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pull(baseline_cases)
  
  # ---- Bind all trial length datasets together ----
  fig4_data <- purrr::imap_dfr(
    model_population_all,
    function(df, trial_length_name) {
      
      if (!metric %in% names(df)) {
        stop("Metric column '", metric, "' not found in model output for trial length ", trial_length_name)
      }
      
      # Add run_id if somehow missing
      if (!"run_id" %in% names(df)) {
        df <- df %>%
          group_by(VE, t) %>%
          mutate(run_id = row_number()) %>%
          ungroup()
      }
      
      df %>%
        select(
          t,
          VE,
          run_id,
          metric_value = all_of(metric)
        ) %>%
        group_by(VE, run_id) %>%
        summarise(
          intervention_cases = sum(metric_value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          baseline_cases = baseline_total,
          cases_averted = baseline_cases - intervention_cases,
          proportion_cases_averted = cases_averted / baseline_cases,
          trial_length = as.numeric(trial_length_name)
        )
    }
  )
  
  # ---- Preserve sensible ordering ----
  fig4_data <- fig4_data %>%
    mutate(
      trial_length = factor(
        trial_length,
        levels = sort(unique(as.numeric(trial_length)))
      )
    )
  
  if (is.factor(fig4_data$VE)) {
    ve_labels <- levels(fig4_data$VE)
  } else {
    ve_labels <- unique(fig4_data$VE)
  }
  
  fig4_data <- fig4_data %>%
    mutate(
      VE = factor(VE, levels = ve_labels)
    )
  
  # ---- Same colour scheme as figure 3 ----
  ve_colours <- grDevices::colorRampPalette(
    c("#F6A6B2", "#9B8FE7", "#7ED8C6", "#FFBE85", "#D01856")
  )(length(ve_labels))
  
  names(ve_colours) <- ve_labels
  
  # ---- Plot ----
  graph <- ggplot(
    fig4_data,
    aes(
      x = trial_length,
      y = proportion_cases_averted,
      fill = VE
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      alpha = 0.75,
      outlier.alpha = 0.35,
      outlier.size = 0.7
    ) +
    scale_fill_manual(
      values = ve_colours,
      name = "RCT-measured VE"
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = c(0, 0)
    ) +
    coord_cartesian(
      ylim = c(0, NA)
    ) +
    labs(
      x = "Trial length (days)",
      y = y_label
    ) +
    theme_clean() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold")
    )
    
  print(graph)
  
  return(
    list(
      graph = graph,
      data = fig4_data
    )
  )
}

# Figure Five : Impact of Mis-parameterisation on Population-level Impact
#' @input model_run : multiple model trial runs
#' @input baseline : baseline model population run
#' @input params : list of model parameters
#' @input model : model function to use for running mis-parameterised scenarios
#' @input metric : metric to use for calculating population-level impact
#' @input waning_durations : vector of waning durations to test
#' @input y_label : label for y-axis
#' @input x_label : label for x-axis

figure_5 <- function(
    model_run,
    baseline,
    params,
    model = model_sirs,
    metric = "I_ma1",
    waning_durations = c(150, 500, 1000),
    y_label = "RCT-measured efficacy",
    x_label = "Error in predicted\npopulation-level impact"
) {
  
  # ---- Baseline medically-attended burden ----
  baseline_total <- baseline %>%
    summarise(
      baseline_cases = sum(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pull(baseline_cases)
  
  # ---- Existing RCT-compatible contour runs ----
  compatible_data <- model_run %>%
    select(t, VE, run_id, metric_value = all_of(metric)) %>%
    group_by(VE, run_id) %>%
    summarise(
      intervention_cases = sum(metric_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      baseline_cases = baseline_total,
      proportion_cases_averted = (baseline_cases - intervention_cases) / baseline_cases
    )
  
  compatible_summary <- compatible_data %>%
    group_by(VE) %>%
    summarise(
      compatible_median = median(proportion_cases_averted, na.rm = TRUE),
      compatible_min = min(proportion_cases_averted, na.rm = TRUE),
      compatible_max = max(proportion_cases_averted, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Extract VE values as numbers ----
  ve_values <- compatible_summary %>%
    mutate(
      VE_numeric = as.numeric(gsub("%", "", as.character(VE))) / 100
    )
  
  # ---- Additional mis-parameterised runs ----
  misparam_outputs <- list()
  
  for (duration in waning_durations) {
    
    message("Running mis-parameterised VE scenarios with waning duration ", duration, " days...")
    
    scenario_outputs <- list()
    
    for (i in seq_len(nrow(ve_values))) {
      
      ve_label <- as.character(ve_values$VE[i])
      ve_numeric <- ve_values$VE_numeric[i]
      
      scenario_outputs[[ve_label]] <- suppressMessages(run_model(
        model = model,
        params = params$inputs,
        dar = params$meta_params,
        param_one = "chi",
        one_values = ve_numeric,
        param_two = "omega",
        two_values = duration,
        all_combos = FALSE
      ) %>%
        output_removeburn(burn_in = params$meta_params$burn) %>%
        output_selectinfantsonly() %>%
        mutate(
          VE = ve_label,
          VE = factor(VE, levels = levels(factor(model_run$VE))),
          waning_duration = duration
        ))
    }
    
    misparam_outputs[[as.character(duration)]] <- bind_rows(scenario_outputs)
  }
  
  misparam_data <- bind_rows(misparam_outputs)
  
  # ---- Population impact from mis-parameterised runs ----
  misparam_impact <- misparam_data %>%
    group_by(VE, waning_duration) %>%
    summarise(
      intervention_cases = sum(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      baseline_cases = baseline_total,
      proportion_cases_averted = (baseline_cases - intervention_cases) / baseline_cases
    )
  
  # ---- Error relative to compatible median ----
  misparam_error <- misparam_impact %>%
    left_join(compatible_summary, by = "VE") %>%
    mutate(
      error = proportion_cases_averted - compatible_median,
      parameterisation = paste0("Waning efficacy (", waning_duration, " days)")
    )
  
  # ---- Compatible uncertainty range around median ----
  compatible_error <- compatible_summary %>%
    mutate(
      xmin = compatible_min - compatible_median,
      xmax = compatible_max - compatible_median,
      parameterisation = "RCT-compatible efficacy & duration"
    )
  
  # ---- Plot data for the fixed-duration runs ----
  duration_error <- misparam_error %>%
    mutate(
      xmin = pmin(0, error),
      xmax = pmax(0, error)
    )
  
  # ---- Combine for output ----
  output_data <- bind_rows(
    compatible_error %>%
      select(VE, parameterisation, xmin, xmax),
    duration_error %>%
      select(VE, parameterisation, xmin, xmax)
  )
  
  parameterisation_levels <- c(
    "RCT-compatible efficacy & duration",
    paste0("Waning efficacy (", waning_durations, " days)")
  )
  
  output_data <- output_data %>%
    mutate(
      VE = factor(VE, levels = rev(levels(factor(model_run$VE)))),
      parameterisation = factor(parameterisation, levels = parameterisation_levels)
    )
  
  # ---- Dynamic colours so names actually match parameterisation levels ----
  duration_colours <- c("#F6A6B2", "#E25178", "#D01856")[seq_along(waning_durations)]
  
  parameterisation_colours <- c(
    "RCT-compatible efficacy & duration" = "#FFA600",
    setNames(
      duration_colours,
      paste0("Waning efficacy (", waning_durations, " days)")
    )
  )
  
  # ---- Make VE numeric positions so we can manually dodge horizontal bars ----
  output_data <- output_data %>%
    mutate(
      VE_num = as.numeric(VE),
      parameterisation_num = as.numeric(parameterisation),
      dodge_offset = (parameterisation_num - mean(seq_along(levels(parameterisation)))) * 0.16,
      y_mid = VE_num + dodge_offset,
      y_min = y_mid - 0.07,
      y_max = y_mid + 0.07
    )
  
  # ---- Background stripes like old version ----
  stripes <- tibble::tibble(
    VE_num = seq_along(levels(output_data$VE)),
    ystart = VE_num - 0.5,
    yend = VE_num + 0.5
  ) %>%
    filter(VE_num %% 2 == 1)
  
  # ---- Plot ----
  graph <- ggplot() +
    
    geom_rect(
      data = stripes,
      aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = ystart,
        ymax = yend
      ),
      fill = "lavenderblush3",
      alpha = 0.20
    ) +
    
    geom_vline(
      xintercept = 0,
      colour = "black",
      linewidth = 0.5
    ) +
    
    geom_rect(
      data = output_data,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = y_min,
        ymax = y_max,
        fill = parameterisation
      ),
      colour = NA
    ) +
    
    scale_x_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    
    scale_y_continuous(
      breaks = seq_along(levels(output_data$VE)),
      labels = levels(output_data$VE),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    
    scale_fill_manual(
      values = parameterisation_colours,
      name = "Model parameterisation"
    ) +
    
    labs(
      x = x_label,
      y = y_label
    ) +
    
    theme_clean() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold")
    )
  
  print(graph)
  
  return(
    list(
      graph = graph,
      data = output_data,
      compatible_summary = compatible_summary,
      misparam_impact = misparam_impact
    )
  )
}

# ---- Supplemental Figures & Analysis ----

# Contact Matrix Visualisation
figure_s1 <- function(cm) {
  
  cm <- as.data.frame(
    as.table(cm),
    responseName = "contacts"
  )
  
  age_labels <- c(
    "[0,1)"   = "<1",
    "[1,2)"   = "1",
    "[2,3)"   = "2",
    "[3,4)"   = "3",
    "[4,5)"   = "4",
    "[5,18)"  = "5–17",
    "[18,26)" = "18–25",
    "[26,40)" = "26–39",
    "[40,60)" = "40–59",
    "[60,75)" = "60–74",
    "75+"     = "75+"
  )
  
  graph <- ggplot(cm,
    aes(x = age.group, y = contact.age.group, fill = contacts)) +
    geom_tile(
      colour = "white",
      linewidth = 0.5
    ) +
    scale_x_discrete(labels = age_labels) +
    scale_y_discrete(labels = age_labels) +
    scale_fill_gradientn(
      colours =  c(
        "#FCE6EE",  # pale blush
        "#F4A1C1",  # soft pink
        "#D85AAF",  # raspberry
        "#8B3FA5",  # plum
        "#29306F"   # muted indigo
      )) +
    coord_equal() +
    guides(
      fill = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(7, "cm"),
        barheight = unit(0.4, "cm")
      )
    ) +
    labs(
      fill = "Mean daily contacts"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(),
    )
  
  print(graph)
  return(graph)
}

# Population Dynamics Over Time
figure_s2 <- function(model_run){
  
  sirs <- model_run %>%
    group_by(t) %>%
    summarise(
      S_Total  = sum(across(starts_with("S") & !matches("Sc|Stx"))),
      I_Total  = sum(across(starts_with("I") & !matches("Iu|Ip|Ic|Itx"))), 
      R_Total  = sum(across(starts_with("R"))),
      P_Total  = sum(across(starts_with("P"))),
      .groups = 'drop'
    ) %>% 
    select(t, S_Total, I_Total, R_Total, P_Total) %>%
    pivot_longer(
      cols = -c(t), 
      names_to = "Compartment",
      values_to = "Population"
    ) %>%
    mutate(
      Compartment = case_when(
        Compartment == "S_Total" ~ "Susceptible",
        Compartment == "R_Total" ~ "Recovered",
        Compartment == "P_Total" ~ "Prophylaxis",
        Compartment == "I_Total" ~ "Infected",
        TRUE ~ Compartment
      )
    )
  
  compartment_colors <- c(
    "Susceptible" = "#F1A0BA",  # soft pink
    "Recovered"   = "#CC4778",  # raspberry
    "Infected"    = "#65AFAD",  # muted teal
    "Prophylaxis" = "#713B91"   # purple
  )
  
  graph <- ggplot(sirs, aes(x = t, y = Population, color = Compartment, group = Compartment)) +
    geom_line(linewidth = 1) +
    labs(
      title = "Population-Level Model Dynamics Over Time",
      x = "Time (Days)",
      y = "Number of People",
      color = "Compartment"
    ) +
    scale_color_manual(values = compartment_colors) +
    scale_y_continuous(
      labels = label_number(scale = 1e-6, suffix = "M"))+
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  print(graph)
  return(graph)
}

# Trial Dynamics Over Time
figure_s3 <- function(model_run){
  
  run <- model_run %>% select(t, Sc, Stx, Ic, Itx) %>%
    pivot_longer(
      cols = -c(t), 
      names_to = "Compartment",
      values_to = "Population"
    ) %>%
    mutate(
      Compartment = case_when(
        Compartment == "Sc" ~ "S (Control)",
        Compartment == "Stx" ~ "S (Treatment)",
        Compartment == "Ic" ~ "I (Control)",
        Compartment == "Itx" ~ "I (Treatment)",
        TRUE ~ Compartment
      )
    )
  
  compartment_colors <- c(
    "S (Control)" = "#F1A0BA",  # soft pink
    "S (Treatment)" = "#9299C1", #pale indigo
    "I (Control)" = "#CC4778",  # raspberry
    "I (Treatment)" = "#29306F"  # muted indigo
  )
  
  graph <- ggplot(run, aes(x = t, y = Population, color = Compartment, group = Compartment)) +
    geom_line(linewidth = 1) +
    labs(
      title = "Model Dynamics Over Time",
      x = "Time (Days)",
      y = "Number of People",
      color = "Compartment"
    ) +
    scale_color_manual(values = compartment_colors) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  print(graph)
  return(graph)
  
}

# ---- Odin Output Tidying ----

output_tidyup <- function(results) {
  
  #when running the model one time only 
  if (length(dim(results)) == 2){
    
    #convert to data frame
    df <- expand.grid(state = rownames(results), t = seq_len(ncol(results)))
    df$value <- mapply(function(state, t) results[state, t], df$state, df$t)
    
    #switch to wide
    df <- df %>%
      pivot_wider(
        names_from = state, 
        values_from = value
      )
    
    return(df)
  }
  
  #when running the model more than once
  else if (length(dim(results)) == 3){ 
    
    #convert to data frame
    df <- expand.grid(
      state = rownames(results), 
      particle = seq_len(ncol(results)), 
      t = seq_len(dim(results)[3])
    )
    df$value <- mapply(function(state,particle, t) results[state, particle, t], 
                       df$state, df$particle, df$t)
    
    #switch to wide
    df <- df %>%
      pivot_wider(
        names_from = c(state),
        values_from = value
      )
    
    return(df)
  }
  
  #if something random happens i guess
  else {
    stop("unsupported number of dimensions for results")
  }
  
}

output_removeburn <- function(model_run, burn_in = 200, follow_up = 365) { 
  output <- model_run %>%
    filter(t > burn_in) %>%
    mutate(t = t - burn_in)
  return(output)
}

output_contours <- function(values, save = FALSE) {
  
  # Cut off unneccessary gunk
  dataset <- values %>%
    mutate(omega = x) %>%
    mutate(chi = y) %>%
    select(order, chi, omega)
  
  if (save) {
    write_csv(dataset, "contour_data.csv")
    message("Contour data saved to contour_data.csv in working directory, ur welc 🫦")
  }
  
  # Output list
  output <- list()
  
  for (i in unique(dataset$order)) {
    output[[paste0("VE_", i)]] <- dataset %>%
      filter(order == i) %>%
      select(-order)
  }
  
  return(output)
}

output_selectinfantsonly <- function(model_run) {
  output <- model_run %>%
    select(t, S1, P1, Iu1, Ip1, I1, Iu_ma1, Ip_ma1, I_ma1, R1, S21) %>%
    mutate(infant_popcheck = S1+P1+I1 + R1 + S21)
  return(output)
}

output_cleanandmerge <- function(input) {
  
  labels <- names(input)
  
  lab_clean <- gsub("VE_", "", labels)
  lab_clean <- gsub("%", "", lab_clean)
  
  lab_num <- suppressWarnings(as.numeric(lab_clean))
  
  clean_labels <- ifelse(
    is.na(lab_num),
    labels,
    ifelse(
      lab_num <= 1,
      paste0(lab_num * 100, "%"),
      paste0(lab_num, "%")
    )
  )
  
  names(input) <- clean_labels
  
  output <- purrr::imap_dfr(
    input,
    function(df, ve_name) {
      df %>%
        group_by(t) %>%
        mutate(run_id = row_number()) %>%
        ungroup() %>%
        mutate(
          VE = ve_name,
          VE = factor(VE, levels = clean_labels)
        )
    }
  )
  
  return(output)
}

output_dynamicfiguring <- function(months = c(90, 180, 270, 360)) {
  
  # --- Get all the contours (aka figure 2 / supplemental figure 4) --- 
  
  # Efficacy Decay Curves for all Trial Lengths
  fig_s4a <- figure_2(model_trial_analysis, fig_1a)
  fig_s4b <- figure_2(model_trial_analysis, fig_1b)
  fig_s4c <- figure_2(model_trial_analysis, fig_1c)
  fig_s4d <- figure_2(model_trial_analysis, fig_1d)
  fig_s4 <- patchwork::wrap_plots(fig_s4a$graph, fig_s4b$graph, fig_s4c$graph, fig_s4d$graph, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  fig_s4
  
  # Save this as a nice graph
  ggsave(filename = "figure_s4.png", plot = fig_s4, width = 9, height = 7, dpi = 300)
  
  # --- U1 Infections for all Decay Curves ---
  fig_s4_allcontours <- list(s4a = fig_s4a$values, s4b = fig_s4b$values, s4c = fig_s4c$values, s4d = fig_s4d$values)
  
  output <- list() # datasets used for supplemental figure 5 + other main figures
  
  for (fig_name in names(fig_s4_allcontours)) {
    
    message("Running population analysis for ", fig_name, "...")
    
    # Run population model from this contour dataset
    model_population_analysis <- run_model_from_contours(
      fig_s4_allcontours[[fig_name]],
      FALSE
    )
    
    # Remove burn-in + select infant-only outputs
    for (i in seq_along(model_population_analysis)) {
      model_population_analysis[[i]] <- output_removeburn(
        model_population_analysis[[i]],
        burn_in = params$meta_params$burn
      )
      
      model_population_analysis[[i]] <- output_selectinfantsonly(
        model_population_analysis[[i]]
      )
    }
    
    # Clean labels + merge list into one dataframe
    model_population_analysis_merged <- output_cleanandmerge(
      model_population_analysis
    )
    
    # Save output
    output[[fig_name]] <- model_population_analysis_merged
    
    # Make graph
    fig_s5 <- figure_3(
      model_population_analysis_merged,
      baseline,
      y_max = 20000
    )
    
    # Save graph
    ggsave(
      filename = paste0("figure_s5_", fig_name, ".png"),
      plot = fig_s5,
      width = 8,
      height = 5,
      dpi = 300
    )
    
    message("Graph saved as figure_s5_", fig_name, ".png 🫡")
  }
  
  names(output) <- months
  
  return(output)
}



