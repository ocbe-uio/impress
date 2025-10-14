# Code to import data from Viedoc into labelled datasets

# The input is a folder with the exported data from Viedoc. The folder should contain
# exported datasets including the MedDRA coding and the CodeLists and Items files.
# The format should be long (one row per activity)

# The output is a data list (raw) where the datasets are stored in the column "data"
# To retreive a data set use the following code 
# 
# dm <- pick(raw, "dm")
# 
# This retrieves the dm dataset
# 
# Note that the numerical category-variables (ending with cd) are kept because they 
# may contain information (i.e. that the numbers have an interpretation). In R all 
# categorical variables (factors) are coded 1, 2, 3 etc regardless of the value in 
#the "cd" variable. If not needed they can be removed by adding the line
# %>% mutate(data = map(data,remove_cd))
# or do it by each dataset like
# dm <- pick(raw, "dm") %>% 
#   remove_cd()

# Both the "pick" and "remove_cd" are functions defined in "functions.R" in the
# external folder

# The "data_lbl" datasets are ready for export to SAS, Stata or SPSS. Here the 
# categorical variables have been reduced to one variable with value labels 
# according to the -labelled- package. Just use e.g. the write_dta-function. 

# Version:
# 2024-08-19: Updated the code according to updates from the MPX-RESPONSE repository. 

library(tidyverse)
library(glue)
library(haven)
library(labelled)
library(readxl)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  export_name <- "_20250729_133232" #use the last exported data
} else if (length(args) != 0) {
  export_name <- args[1]
}


export_folder <- glue("data/raw/{export_name}")
prefix <- export_name # needs to be adjusted according to the specifications of the output file name
   
source(file.path("R/external/functions.R"))


# Compile the codelists from the CodeList file and from the CodeList sheet in the MRI excel files

cl <- read_csv(file.path(export_folder,paste0(export_name,"_CodeLists.csv")), skip = 1, show_col_types = FALSE) %>%
  rename_all(tolower) %>%
  #mutate(datatype = if_else(formatname == "ICMFFMT", "integer", datatype)) %>% 
  group_by(formatname) %>%
  nest(value_labels = c(datatype, codevalue, codetext))

# Add the codelists from the MRI excel files 

my_read_cl <- function(file){
  data <- read_excel(file, skip = 1, sheet = "CodeLists") %>%
    rename_all(tolower) 
  return(data)
}

files <- list.files("data/raw/_20250729_133232/FileData/", recursive = TRUE)

mricl <- tibble(file = paste0("data/raw/_20250729_133232/FileData/", files)) %>% 
  mutate(raw = map(file, my_read_cl)) %>%
  unnest(raw) %>% 
  select(-file) %>% 
  distinct(formatname, datatype, codevalue, codetext) %>% 
  group_by(formatname) %>% 
  nest(value_labels = c(datatype, codevalue, codetext)) 

cl <- bind_rows(cl, mricl) %>% 
  distinct(formatname, .keep_all = TRUE) %>% 
  arrange(formatname) %>% 
  ungroup()

#############################
# Compile the items from the Items file and from the Items sheet in the MRI excel files
#############################

items <- read_csv(file.path(export_folder, paste0(export_name, "_Items.csv")), skip = 1, show_col_types = FALSE) %>%
  rename_all(tolower) %>%
  mutate(id_ = id,
         id = tolower(id)) %>%
  mutate(categorical = if_else(!is.na(formatname), 2,
                               if_else(paste0(id, "cd") == lead(id), 1, 0)
  )) %>%
  mutate(formatname = if_else(categorical == 1, lead(formatname), formatname))

my_read_excel <- function(file){
  my_cols <- read_excel(file, skip = 1, n_max = 1, sheet = "MRI") %>% 
    colnames() %>% 
    as_tibble_col(column_name = "id_")
  
  my_items <- read_excel(file, skip = 1,  sheet = "Items") %>%
    rename_all(tolower) %>%
    mutate(id_ = id,
           id = tolower(id)) %>%
    mutate(categorical = if_else(!is.na(formatname), 2,
                                 if_else(paste0(id, "cd") == lead(id), 1, 0)
    )) %>%
    mutate(formatname = if_else(categorical == 1, lead(formatname), formatname)) %>%
    rename_all(tolower) %>% 
    mutate(cols_abb = case_when(
      datatype == "date" ~ "date",
      datatype == "datetime" ~ "text",
      datatype == "double" ~ "numeric",
      datatype == "integer" ~ "numeric",
      datatype == "string" ~ "text",
      datatype == "text" ~ "text"
    ))
  
  my_col_types <- my_cols %>% 
    left_join(my_items %>% select(id_, cols_abb), by = "id_") %>% 
    mutate(cols_abb = if_else(is.na(cols_abb), "guess", cols_abb)) %>% 
    select(cols_abb) %>% 
    as_vector()
  
  data <- read_excel(file, col_types = my_col_types, skip = 1, sheet = "MRI") %>% 
    rename_all(tolower)
  return(data)
}


my_read_items <- function(file){
  data <- read_excel(file, skip = 1, sheet = "Items") %>%
    rename_all(tolower) %>%
    mutate(id_ = id,
           id = tolower(id)) %>%
    mutate(categorical = if_else(!is.na(formatname), 2,
                                 if_else(paste0(id, "cd") == lead(id), 1, 0)
    )) %>%
    mutate(formatname = if_else(categorical == 1, lead(formatname), formatname)) %>%
    rename_all(tolower)
  return(data)
}

mriitems <- tibble(file = paste0("data/raw/_20250729_133232/FileData/", files)) %>% 
  mutate(raw = map(file, my_read_items)) %>%
  unnest(raw) %>%
  select(-file) %>% 
  distinct(id, label, datatype, mandatory, decimals, minlength, maxlength, formatname, contentlength, id_, categorical) %>% 
  mutate(categorical = if_else(datatype == "double", 0, categorical)) %>% 
  add_row(id = "mrvis", label = "MRI Visit Identifier", datatype = "text", 
          mandatory = "True", decimals = NA, minlength = "1", maxlength = NA, 
          formatname = NA, contentlength = "30", id_ = "MRVIS", categorical = 0) %>% 
  mutate(mandatory = as.logical(mandatory), 
         decimals = as.integer(decimals),
         minlength = as.integer(minlength),
         maxlength = as.integer(maxlength),
         contentlength = as.integer(contentlength)
         )

items <- bind_rows(items, mriitems) %>% 
  distinct(id, .keep_all = TRUE)


items <- items %>%
  left_join(cl, by = "formatname") %>%
  rename_all(tolower) %>% 
  mutate(cols_abb = case_when(
    datatype == "date" ~ "c",
    datatype == "datetime" ~ "c",
    datatype == "double" ~ "d",
    datatype == "integer" ~ "i",
    datatype == "string" ~ "c",
    datatype == "text" ~ "c"
  ))



# Read all datasets

my_read_excel <- function(file){
  my_cols <- read_excel(file, skip = 1, n_max = 1, sheet = "MRI") %>% 
    colnames() %>% 
    as_tibble_col(column_name = "id_")
  
  my_items <- read_excel(file, skip = 1,  sheet = "Items") %>%
    rename_all(tolower) %>%
    mutate(id_ = id,
           id = tolower(id)) %>%
    mutate(categorical = if_else(!is.na(formatname), 2,
                                 if_else(paste0(id, "cd") == lead(id), 1, 0)
    )) %>%
    mutate(formatname = if_else(categorical == 1, lead(formatname), formatname)) %>%
    rename_all(tolower) %>% 
    mutate(cols_abb = case_when(
      datatype == "date" ~ "date",
      datatype == "datetime" ~ "text",
      datatype == "double" ~ "numeric",
      datatype == "integer" ~ "numeric",
      datatype == "string" ~ "text",
      datatype == "text" ~ "text"
    ))
  
  my_col_types <- my_cols %>% 
    left_join(my_items %>% select(id_, cols_abb), by = "id_") %>% 
    mutate(cols_abb = if_else(is.na(cols_abb), "guess", cols_abb)) %>% 
    select(cols_abb) %>% 
    as_vector()

  data <- read_excel(file, col_types = my_col_types, skip = 1, sheet = "MRI") %>% 
    rename_all(tolower)
  return(data)
}


my_read_csv <- function(file){
  my_cols <- read_csv(file, skip = 1, n_max = 1, show_col_types = FALSE) %>% 
    colnames() %>% 
    as_tibble_col(column_name = "id_")
  
  my_col_types <- my_cols %>% 
    left_join(items %>% select(id_, cols_abb), by = "id_") %>% 
    mutate(cols_abb = if_else(is.na(cols_abb), "?", cols_abb)) 
  
  my_col_types <- setNames(as.list(my_col_types$cols_abb), my_col_types$id_)
  
  data <- read_csv(file, col_types = my_col_types, skip = 1, show_col_types = FALSE)
  return(data)
}


raw <- tibble(files = list.files(export_folder)) %>%
  filter(files != "FileData") %>% 
  mutate(
    id = str_remove(files, paste0(prefix, "_")),
    id = str_remove(id, ".csv"),
    id = str_to_lower(id),
    files = file.path(export_folder, files)
  ) %>%
  filter(!(id %in% c("items", "codelists", "readme.txt"))) %>%
  filter(!endsWith(id, ".sas")) %>% 
  mutate(txt = map(files, my_read_csv)) %>%
  mutate(problems = map(txt, problems),
         any_problems = map_chr(problems, \(x) if_else(nrow(x) != 0, "Yes", "No"))) %>% 
  mutate(txt = map(txt, rename_all, tolower)) %>%
  mutate(txt = map(txt, labeliser, codelist = items)) %>%
  mutate(data_lbl = map(txt, factoriser, codelist = items)) %>%
  mutate(data = map(data_lbl, haven::zap_labels)) %>% 
  add_row(files = file.path(export_folder, glue("{export_name}_CodeLists.csv")), 
          id = "codelist", 
          txt = list(cl), 
          data = list(cl), 
          data_lbl = list(cl), 
          any_problems = "No") %>% 
  add_row(files = file.path(export_folder, glue("{export_name}_Items.csv")), 
          id = "items", 
          txt = list(items), 
          data = list(items), 
          data_lbl = list(items),
          any_problems = "No") 

# Read the MR files and add to raw

mridata <- tibble(file = paste0("data/raw/_20250729_133232/FileData/", files)) %>% 
  mutate(raw = map(file, my_read_excel)) %>% 
  unnest(raw) %>%
  mutate(mrvis = eventid) %>% 
  mutate(eventid = str_sub(file, 48, 51)) %>%
  mutate(eventid = str_remove(eventid, "/")) %>% 
  mutate(subjectid = str_extract(subjectid, "[:digit:]+")) %>%
  mutate(subjectid = paste0("ImPR-", subjectid)) %>% 
  select(!c(file, siteseq:subjectseq, eventseq, eventname:mrfuplyncd))

mri <- raw %>% 
  pick("mri") %>% 
  select(siteseq:mrfuplyncd) %>% 
  left_join(mridata, by = c("subjectid", "eventid"))



raw <- raw %>% 
  mutate(txt = if_else(id == "mri", list(mri), data),
         data = if_else(id == "mri", list(mri), data),
         data_lbl = if_else(id == "mri", list(mri), data_lbl
         ))

raw <- raw %>% 
  #mutate(data = map(data,remove_cd)) %>%  # de-comment to remove the "cd" variables
  mutate(data = map(data, labeliser, codelist = items)) %>% 
  mutate(data_lbl = map(data_lbl, remove_fct, codelist = items)) %>% 
  mutate(data_lbl = map(data_lbl, labeliser, codelist = items)) %>% 
  set_variable_labels(
    files = "Path to file",
    id = "ID",
    txt = "Raw import",
    data_lbl = "With value labels for export",
    data = "For analyses (value labels removed)",
    any_problems = "Any problems with the import?",
    problems = "Table of problems"
  ) %>% 
  select(files:id, any_problems, problems, txt, data_lbl:data) %>% 
  mutate(problems = if_else(any_problems == "Yes", problems, as.vector(NA)))

raw_path <- file.path("data", "raw", "raw.rds")


  
write_rds(raw, raw_path)
