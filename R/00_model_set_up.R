# ---------------------------------------------------------------------------- #

# ---- RSV Model Set Up File ---- #
# By : Kate Turpie
# Updated : 07-Jul-2026

# ---------------------------------------------------------------------------- #

# ---- Load All Necessary Packages ----

library(odin2)
library(dust2)
library(monty)
library(tidyverse)
library(socialmixr)
library(ggthemes)
library(scales)
library(viridis)
library(purrr)
library(metR)

# ---- Source All Relevant Files & Functions ----

source("./R/01_model_sirs_odin.R")
source("./R/02_model_params.R")
source("./R/03_model_calibrate.R")
source("./R/04_model_outputs.R")