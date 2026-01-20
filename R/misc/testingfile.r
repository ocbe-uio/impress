


  library(tibble)
  library(labelled)
  library(dplyr)
  library(purrr)
  library(readxl)
  library(glue)
  library(haven)
  library(readr)
  library(stringr)
  library(tidyr)
  library(yaml)
  library(fs)
  library(here)
  library(rlang)
  
library(tidyverse)
source("R/external/functions.R")
source("R/make_ad/helpers_ad.R")
source("R/make_ad/make_adsl.R")
source("R/make_raw/make_shamrand.R")
source("R/make_raw/make_raw.R")
source("R/make_rd/helpers_rd.R")
cfg <- yaml::read_yaml("config/cfg.yml")
cfg_adam <- yaml::read_yaml("config/adam.yml")
#raw <- make_raw(cfg)
rawsham <- make_shamrand( raw, cfg)
adsl <- make_adsl(rawsham, cfg_adam, cfg, identifiers = "subjectid")

paste0(rantrt, if_else(starts_with(rantrt, 'A'), paste0('_', DOSE02P), ''))

# R console

# Restart your R session.
rstudioapi::restartSession()

library(targets)
library(tarchetypes)
tar_make()

# Loads globals like tar_option_set() packages, simulate_data(), and analyze_data():
tar_load_globals()

# Load the data that the target depends on.
stems <- targets::tar_meta()[['name']][targets::tar_meta()[['type']] == "stem"]
targets::tar_load(all_of(stems))

# Run the command of the errored target.
make_adsl(rawsham, cfg_adam, cfg, identifiers = "subjectid")


pick(raw, "event_dates") %>%
  select(eventname) %>%
  distinct() %>%
  writexl::write_xlsx("metadata/visits2.xlsx")
