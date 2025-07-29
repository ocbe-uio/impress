# Code to make the ran dataset using a mock-up randomisation.
# Used until the real randomisation is available
#

library(tidyverse)
library(labelled)
source("src/external/functions.R")

# Pick the dm table from the raw data

dm <- raw |> pick("dm")

# Pick the rand table from the raw data

ran <- raw |> pick("ran")
elig <- raw |> pick("elig")



# Create a new dataset with mock-up randomisation data

seed <- 1234

set.seed(seed) 


alloc <- tibble(diagno = 1:3,
                Diagnosis = c("Recurrent glioblastoma",
                              "Newly diagnosed glioblastoma",
                              "Brain metastases")
                ) %>% 
  group_by(diagno) %>%
  crossing(blockno = 1:20) %>%
  group_by(diagno, blockno) %>%
  mutate(blocksize = if_else(diagno == 3, sample(c(3,6),1), 18)) %>%
  mutate(allocation = if_else(diagno == 3, 
                              map(blocksize, ~rep(1:3,.x/3)), 
                              list(1:18))) %>%
  unnest(cols = allocation) %>%
  mutate(rand = runif(n())) %>%
  arrange(diagno, blockno,rand) %>%
  group_by(diagno) %>%
  mutate(randno = diagno * 1000 + 1:n())


  randgrp <- tibble(diagno = 1:3) %>%
    crossing(allocation = 1:18) %>%
    filter(diagno != 3 | (diagno == 3 & allocation < 4)) %>%
    mutate(dose = if_else(diagno == 3, 50,
                            case_when(
                              between(allocation, 1, 6) ~ 25,
                              between(allocation, 7, 12) ~ 50,
                              between(allocation, 13, 18) ~ 100
                            ))) %>%
    mutate(startup = if_else(diagno == 3, 
                             case_when(
                               allocation == 1 ~ 1,
                               allocation == 2 ~ 91, 
                               allocation == 3 ~ 181
                             ), 
                             case_when(
                               allocation %in% c(1, 4, 7, 10, 13, 16) ~ 1,
                               allocation %in% c(2, 5, 8, 11, 14, 17) ~ 15,
                               allocation %in% c(3, 6, 9, 12, 15, 18) ~ 29
                             ))) %>%
    mutate(fu_dose = if_else(diagno == 3, NA_real_,
                               case_when(
                                 between(allocation, 1, 3) ~ 25,
                                 between(allocation, 7, 9) ~ 50,
                                 between(allocation, 13, 18) ~ 100, 
                                 TRUE ~ 0
                               ))) %>%
    mutate(Group =if_else(diagno == 3, paste0(allocation), 
                          case_when(
                            allocation %in% c(1,4) ~ "1_1",
                            allocation %in% c(2,5) ~ "1_2",
                            allocation %in% c(3,6) ~ "1_3",
                            allocation %in% c(7,10) ~ "2_4",
                            allocation %in% c(8,11) ~ "2_5",
                            allocation %in% c(9,12) ~ "2_6",
                            allocation %in% c(13,16) ~ "3_1",
                            allocation %in% c(14,17) ~ "3_2",
                            allocation %in% c(15,18) ~ "3_3"
                          ))) %>%
    mutate(Group = case_when(
      diagno == 1 ~ paste0("AR", Group),
      diagno == 2 ~ paste0("AN", Group),
      diagno == 3 ~ paste0("BM", Group),
    ))

    
  alloc <- alloc %>%
    left_join(randgrp, by = c("diagno","allocation")) %>%
    ungroup %>% 
    rename(randiacd = diagno, 
           randia = Diagnosis, 
           ranno = randno, 
           rantrt = Group, 
           randos = dose, 
           ranst = startup,
           ranfudos = fu_dose)


ran0 <- dm %>% 
  select(siteseq, sitename, subjectseq, subjectid, cohort) %>% 
  left_join(elig %>% select(subjectid, eventdate), by = "subjectid") %>%
  rename(randat = eventdate, randia = cohort) %>% 
  mutate(diagno = case_when(
    randia == "Recurrent glioblastoma" ~ 1,
    randia == "Newly diagnosed glioblastoma" ~ 2,
    randia == "Brain metastases" ~ 3
  )) %>%
  group_by(randia) %>% 
  mutate(ranno = diagno * 1000 + 1:n()) %>% 
  left_join(alloc, by = c("randia", "ranno")) %>% 
  select(-blockno, -blocksize, -allocation, -rand, -diagno) %>% 
  mutate(randos = paste0(randos, "mg"),
          ranst = case_when(
            ranst == 1 ~ "Cycle 1 - day 1",
            ranst == 15 ~ "Cycle 2 - day 15",
            ranst == 29 ~ "Cycle 3 - day 29",
            ranst == 91 ~ "Cycle 2 - day 91",
            ranst == 181 ~ "Cycle 3 - day 181",
            TRUE ~ "NA"
          ),
         ranfudos = paste0(ranfudos, "mg")) %>%
  mutate(randia = factor(randia, levels = levels(ran$randia)),
         rantrt = factor(rantrt, levels = levels(ran$rantrt)),
         randos = factor(randos, levels = levels(ran$randos)),
         ranst = factor(ranst, levels = levels(ran$ranst)),
         ranfudos = factor(ranfudos, levels = levels(ran$ranfudos))) %>%
  copy_labels_from(ran)
  
  
  
