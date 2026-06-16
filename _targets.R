# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(quarto)
# library(tarchetypes) # Load other packages as needed.


# Set target options:
tar_option_set(
  packages = c(
    "tibble", "labelled", "dplyr", "purrr", "readxl", "glue",
    "haven", "readr", "stringr", "tidyr", "yaml", "fs", "here",
    "rlang", "admiral", "lubridate", "quarto", "tern", "rtables.officer",
    "ggplot2", "lme4", "emmeans", "broom", "DoseFinding", "survival", "survminer",
    "patchwork"
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
tar_source("R/make_ad/helpers_ad.R")
tar_source("R/make_ad/make_adsl.R")
tar_source("R/make_raw/make_shamrand.R")
tar_source("R/make_raw/make_raw.R")
tar_source("R/make_ad/make_admri.R")
tar_source("R/make_ad/make_adeff.R")
tar_source("R/make_ad/make_adbl.R")
tar_source("R/make_ad/make_adtte.R")
tar_source("R/make_rd/make_rdbl.R")
tar_source("R/make_rd/helpers_rd.R")
tar_source("R/make_rd/make_rdmri.R")
tar_source("R/make_rd/make_rdeff.R")
tar_source("R/make_rd/make_rdtte.R")
tar_source("R/make_ad/make_adlb.R")
tar_source("R/make_rd/make_rdlb.R")
# Safety (SAP §7)
tar_source("R/make_ad/make_adae.R")
tar_source("R/make_ad/make_advs.R")
tar_source("R/make_ad/make_adex.R")
tar_source("R/make_ad/make_adlbsaf.R")
tar_source("R/make_rd/make_rdae.R")
tar_source("R/make_rd/make_rdvs.R")
tar_source("R/make_rd/make_rdlbsaf.R")
tar_source("R/make_rd/make_rdex.R")

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
    dat,
    effective_raw(raw, shamraw, cfg)
  ),
  tar_target(
    adsl,
    make_adsl(dat, cfg)
  ),
  tar_target(
    admri,
    make_admri(dat, adsl, cfg)
  ),
  tar_target(
    adeff,
    make_adeff(dat, adsl, cfg)
  ),
  tar_target(
    adbl,
    make_adbl(dat, cfg, adsl, admri)
  ),
  tar_target(
    adtte,
    make_adtte(dat, adsl, cfg)
  ),
  tar_target(
    adlb,
    make_adlb(dat, adsl, cfg)
  ),
  tar_target(
    lb_vars,
    rlang::set_names(unique(as.character(adlb$paramcd)))
  ),
  tar_target(
    lb_summaries,
    make_rdlb_batch(adlb, lb_vars, cfg)
  ),
  tar_target(
    lb_late_summaries,
    make_rdlb_late_batch(adlb, lb_vars, cfg, lb_summaries)
  ),
  tarchetypes::tar_render(
    lb_modelfit_report,
    path       = "reports/model fit/impress_lb_model_fit.Rmd",
    output_format = "pdf_document",
    output_file   = "impress_lb_model_fit.pdf",
    output_dir    = "reports/model fit"
  ),
  tarchetypes::tar_render(
    tte_ph_report,
    path       = "reports/model fit/impress_tte_ph_check.Rmd",
    output_format = "pdf_document",
    output_file   = "impress_tte_ph_check.pdf",
    output_dir    = "reports/model fit"
  ),
  tar_target(
    tte_vars,
    rlang::set_names(unique(as.character(adtte$paramcd)))
  ),
  tar_target(
    tte_sections,
    purrr::set_names(tte_vars) |>
      purrr::map(~ make_tte_section(adtte, .x, cfg))
  ),
  tar_target(
    primary_mri_vars,
    rlang::set_names(c("mrt1cf", "mrt2cf", "mrt1ss", "mrt2ss"))
  ),
  tar_target(
    primary_mri_labels,
    c(
      mrt1cf = "Mean relative cerebral blood flow (PTA-ROI)",
      mrt2cf = "Mean relative cerebral blood flow (ROI2)",
      mrt1ss = "Mean relative solid stress (PTA-ROI)",
      mrt2ss = "Mean relative solid stress (ROI2)"
    )
  ),
  tar_target(
    primary_mri_sections,
    make_rdmri_batch(admri, primary_mri_vars, cfg)
  ),
  tar_target(
    other_mri_vars,
    setdiff(unique(as.character(admri$paramcd)), primary_mri_vars)
  ),
  tar_target(
    other_mri_summaries,
    purrr::set_names(other_mri_vars) |>
      purrr::map(~ summarize_mri_endpoint(admri, .x, cfg))
  ),
  tar_target(
    neuro_vars,
    rlang::set_names(c("kps", "ecog", "nano_tot"))
  ),
  tar_target(
    neuro_summaries,
    purrr::set_names(neuro_vars) |>
      purrr::map(~ make_cont_section(adeff, .x, cfg))
  ),
  tar_target(
    qol_vars,
    rlang::set_names(c("qlq_c30", "qlq_bn20"))
  ),
  tar_target(
    qol_summaries,
    purrr::set_names(qol_vars) |>
      purrr::map(~ make_cont_section(adeff, .x, cfg))
  ),
  tar_target(
    steroid_summaries,
    make_steroid_section(adeff, cfg)
  ),
  tar_target(
    rano_vars,
    rlang::set_names(c("trsdisea", "trorresp"))
  ),
  tar_target(
    rano_summaries,
    purrr::set_names(rano_vars) |>
      purrr::map(~ summarise_rano_polr(adeff, .x, cfg))
  ),
  tar_target(
    cycle11_mri_vars,
    rlang::set_names(unique(as.character(admri$paramcd)))
  ),
  tar_target(
    cycle11_mri_summaries,
    purrr::set_names(cycle11_mri_vars) |>
      purrr::map(~ summarize_mri_cycle11(admri, .x, cfg))
  ),
  # ---- Safety (SAP §7) ----------------------------------------------------
  tar_target(
    adae,
    make_adae(dat, adsl, cfg)
  ),
  tar_target(
    advs,
    make_advs(dat, adsl, cfg)
  ),
  tar_target(
    adex,
    make_adex(dat, adsl, cfg)
  ),
  tar_target(
    adlbsaf,
    make_adlbsaf(dat, adsl, cfg)
  ),
  tar_target(
    ae_summaries,
    make_rdae(adae, adsl, cfg)
  ),
  tar_target(
    vs_vars,
    rlang::set_names(c("SBP", "DBP", "PULSE"))
  ),
  tar_target(
    vs_summaries,
    make_rdvs_batch(advs, vs_vars, cfg)
  ),
  tar_target(
    lbsaf_summaries,
    make_rdlbsaf(adlbsaf, cfg)
  ),
  tar_target(
    ex_summary,
    make_rdex(adex, cfg)
  ),
  tar_target(
    deaths_summary,
    make_deaths_section(adtte, cfg)
  ),
  tar_target(
    tbl_baseline,
    make_rdbl(adbl, cfg),
    format = "rds"
  ),
  tarchetypes::tar_render(
    report_docx,
    path = "reports/impress_statistical_analysis.Rmd",
    output_format = "word_document",
    output_file = "impress_statistical_analysis.docx",
    output_dir = "reports"
  ),
  tarchetypes::tar_render(
    report_pdf,
    path = "reports/impress_statistical_analysis.Rmd",
    output_format = "pdf_document",
    output_file = "impress_statistical_analysis.pdf",
    output_dir = "reports"
  )

)
