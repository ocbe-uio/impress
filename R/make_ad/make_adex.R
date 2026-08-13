# ADEX dataset creation — extent of exposure (SAP §7.7)
#


make_adex <- function(raw, adsl, cfg) {

  visits <- visits_from_yaml(cfg)
  ex <- raw |> get_raw("ex")

adex_base <-
    adsl %>%
    left_join(
      ex %>%
        transmute(
          subjid    = subjectid,
          eventname = eventname,
          adt       = as_date(eventdate),
          avisit    = lookup_avisit(eventname, visits$map, fallback_title = TRUE),
          avisitn   = lookup_avisitn(eventname, visits$map, unmapped = visits$defaults$avisitn_unmapped),
          across(all_of(starts_with("ex")))
        ),
      by = "subjid"
    ) |>
    mutate(
      randdt = as_date(randdt),
      ablfl = if_else(avisit == "Randomization And Baseline MRI", "Y", "N", missing = "N")
    ) |> 
    mutate(exstdat = as_date(exstdat),
            exendat = as_date(exendat))

 
  
  # Summarise the duration, dose and daily dose per participant
  exsummary <- adex_base |> 
    mutate(exstdat = as_date(exstdat),
            exendat = as_date(exendat)) |> 
    mutate(trtdur = exendat - exstdat,
         trtdur = as.numeric(trtdur),
         exdosstrcd = str_replace(exdosstrcd, ",", "."),
         exdosstrcd = as.numeric(exdosstrcd)) |> 
    mutate(trtdose1 = excons*exdosstrcd,
        trtdose2 = ex2cons * 12.5,
        trtdose = rowSums(cbind(trtdose1, trtdose2), na.rm = TRUE)) |> 
    group_by(subjid) |> 
    summarise(totdur = max(trtdur, na.rm = TRUE), 
            totdose = sum(trtdose, na.rm = TRUE), 
            moddose = any(exmodyncd == 1, na.rm = TRUE),
            .groups = "drop_last") |> 
    mutate(dailydose = totdose/totdur)
  
  adex <- adsl |> 
    left_join(exsummary, by = "subjid")
  
  
  


}

# Standalone execution -------------------------------------------------------
if (sys.nframe() == 0L) {
  suppressMessages({
    library(dplyr); library(tidyr); library(stringr); library(purrr)
    library(tibble); library(lubridate); library(yaml); library(targets)
  })
  source("R/external/functions.R")
  source("R/make_ad/helpers_ad.R")

  #cfg <- yaml::read_yaml("config/cfg.yml")
  tar_load(raw); tar_load(adsl); tar_load(cfg)
  adex <- make_adex(raw, adsl, cfg)
  message("adex: ", nrow(adex), " subjects")
  print(adex %>% dplyr::count(trt02p))
}

  
  