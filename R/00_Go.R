# ---- Click Thru ----
source("R/00_model_set_up.R")

# ---- Baseline Model Run ----
params <- input_all()

baseline <- run_model(
  model = model_sirs,
  params = params$input,
  dar = params$meta_params,
  deterministic = TRUE
)

# ---- Model Runs - Trial ---- 

# Set Parameters
params$meta_params$trial_duration <- 90 
chi_values <- seq(0.2, 1, by = .05)
omega_values <- seq(90, 1080, by = 15)

# Run Model
model_trial_analysis <- run_model(
  model = model_sirs,
  params = params$input,
  dar = params$meta_params,
  deterministic = TRUE, 
  param_one = "chi", 
  one_values = chi_values,
  param_two = "omega", 
  two_values = omega_values, 
  all_combos = TRUE
)

# Tidy Output
model_trial_analysis <- output_removeburn(model_trial_analysis, burn_in = params$meta_params$burn)
model_trial_analysis <- model_trial_analysis %>% mutate(efficacy_measured = ((Ic/500 - (Itx/500))/(Ic/500)))

# Figures
fig_1a <- figure_1(model_trial_analysis, trial_len = 90)
fig_1b <- figure_1(model_trial_analysis, trial_len = 180)
fig_1c <- figure_1(model_trial_analysis, trial_len = 270)
fig_1d <- figure_1(model_trial_analysis, trial_len = 365)

fig_1 <- patchwork::wrap_plots(fig_1a, fig_1b, fig_1c, fig_1d, ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")
fig_1

fig2 <- figure_2(model_trial_analysis, fig_1b, multi = FALSE)

# Save EVERYTHING (except the subfigures, I don't want those 💅)
# ((But you can save them if you want, just un-comment them 💖))
write_csv(model_trial_analysis, "model_trial_analysis.csv")
ggsave(filename = "figure_1.png", plot = fig_1, width = 9, height = 7, dpi = 300)
ggsave(filename = "figure_2.png", plot = fig_2, width = 9, height = 7, dpi = 300)

#ggsave(filename = "figure_1a.png", plot = fig_1a, width = 9, height = 7, dpi = 300)
#ggsave(filename = "figure_1b.png", plot = fig_1b, width = 9, height = 7, dpi = 300)
#ggsave(filename = "figure_1c.png", plot = fig_1c, width = 9, height = 7, dpi = 300)
#ggsave(filename = "figure_1d.png", plot = fig_1d, width = 9, height = 7, dpi = 300)

# ---- Model Runs - Population ----

# Set Parameters
update_initial_conditions <- input_initial_conditions(demog = params$inputs$demog, start_P_seed = TRUE)
params$inputs$start_S <- update_initial_conditions$start_S
params$inputs$start_P <- update_initial_conditions$start_P

# Run the Model for Each Efficacy
model_population_analysis <- run_model_from_contours(fig2$values, TRUE)

# Tidy Output
for (i in 1:length(model_population_analysis)) {
  model_population_analysis[[i]] <- output_removeburn(model_population_analysis[[i]], burn_in = params$meta_params$burn)
  model_population_analysis[[i]] <- output_selectinfantsonly(model_population_analysis[[i]])
}
model_population_analysis_merged <- output_cleanandmerge(model_population_analysis)
baseline <- output_removeburn(baseline, burn_in = params$meta_params$burn)

# Additional Model Running (for Figures 4 - 5)
model_population_all <- output_dynamicfiguring()

# Figures 
fig_3 <- figure_3(model_population_analysis_merged, baseline, y_max = 20000)
fig_4 <- figure_4(
  model_population_all = model_population_all,
  baseline = baseline,
  metric = "I_ma1"
)

fig_4$graph
fig_5 <- figure_5(
  model_run = model_population_all[["180"]],
  baseline = baseline,
  params = params
)

fig_5$graph

ggsave(filename = "figure_3.png", plot = fig_3, width = 9, height = 7, dpi = 300)
ggsave(filename = "figure_4.png", plot = fig_4, width = 9, height = 7, dpi = 300)
ggsave(filename = "figure_5.png", plot = fig_5, width = 9, height = 7, dpi = 300)

# ---- Supplemental Figures that are Just for Vibes ----

# Contact Matrix
fig_s1 <- figure_s1(cm = params$inputs$contacts) 

# Model Baseline SIR Graph for Dynamic Transmission Model & Trial Model 
fig_s2 <- figure_s2(baseline)
fig_s3 <- figure_s3(baseline)

# Save Supplements
ggsave(filename = "figure_s1.png", plot = fig_s1, width = 9, height = 7, dpi = 300)
ggsave(filename = "figure_s2.png", plot = fig_s2, width = 9, height = 7, dpi = 300)
ggsave(filename = "figure_s3.png", plot = fig_s3, width = 9, height = 7, dpi = 300)

