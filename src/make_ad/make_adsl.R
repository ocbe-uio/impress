# This code makes the adsl.rds file for the Impress trial

library(tidyverse)
source("src/external/functions.R")

library(purrr)
library(dplyr)

raw <- read_rds("data/raw/raw.rds")

dm <- raw %>%
    pick("dm")

ran <- raw %>% 
    pick("ran")

if (nrow(ran) == 0) {
  source("src/make_raw/make_pseudorand.R")
  ran <- ran0
}



adsl0 <- dm |>
    mutate(
        studyid = "ImPRESS",
        subjid = subjectid,
        siteid = "OUS",
        age = dmage,
        ageu = "years"
    ) |>
    select(studyid, subjid, siteid, age, ageu, sex)

