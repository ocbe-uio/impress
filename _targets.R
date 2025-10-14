# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) 
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "tibble", "labelled", "dplyr", "purrr", "readxl", "glue",
    "haven", "readr", "stringr", "tidyr", "yaml", "fs", "here",
    "rlang"
  ) # Packages that your targets need for their tasks.
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/external/functions.R")
tar_source("R/make_ad/make_adsl.R")
tar_source("R/make_raw/make_shamrand.R")
tar_source("R/make_raw/make_raw.R")
tar_source("R/make_ad/make_adeff.R")


# Replace the target list below with your own:
list(
  tar_target(
    cfg_adam_file,
    "config/adam.yml",
    format = "file"
  ),
  tar_target(
    cfg_adam,
    yaml::read_yaml(cfg_adam_file)
  ),
  tar_target(
    cfg_file,
    "config/cfg.yml",
    format = "file"
  ),
  tar_target(
    cfg,
    yaml::read_yaml(cfg_file)
  ),
  tar_target(
    raw,
    make_raw(cfg)
  ),
  tar_target(
    shamraw,
    make_shamrand(raw, cfg)
  ),
  tar_target(
    adsl,
    make_adsl(shamraw, cfg_adam, cfg, identifiers = "subjectid")
  ), 
  tar_target(
    raw_domains,
    as_raw_domains(shamraw)
  ),
  tar_target(
    adeff,
    make_adeff(
      lineage_yml  = cfg_adam_file,
      raw          = raw_domains,
      adsl         = adsl,
      cfg          = cfg  
    )
  )
)
