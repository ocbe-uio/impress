# install.packages(c("tidyverse","readr","janitor","writexl","stringr"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(writexl)
  library(stringr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# Convert ITEMS + CODELISTS CSVs --> metacore_map.xlsx (P21-style)
# =============================================================================
# Expected content (names are flexible; map via `colmap_*`):
#  - ITEMS CSV: dataset/domain, variable, label, type, length, origin, codelist,
#               key flag, structure (optional), value-level fields (optional)
#  - CODELISTS CSV: codelist name, value/code, decode/meaning
#
# Output Excel has sheets: Datasets, Variables, Codelists, ValueLevel
# =============================================================================

items_codelists_to_metacore <- function(
    items_csv  = "metadata/_20250729_133232_Items.csv",
    codelists_csv = "metadata/_20250729_133232_CodeLists.csv",
    out_xlsx   = "metadata/metacore_map.xlsx",
    
    # Map your CSV headers to logical roles (edit as needed)
    colmap_items = list(
      dataset   = c("dataset","domain","adam_dataset","analysis_dataset"),
      variable  = c("variable","var","variable_name","varname"),
      label     = c("label","var_label","variable_label","description"),
      type      = c("type","datatype","data_type"),
      length    = c("length","len","display_length"),
      origin    = c("origin","source","data_origin"),
      codelist  = c("codelist","ct","controlled_terminology","codelist_name"),
      key       = c("key","is_key","key_flag","primary_key"),
      structure = c("structure","adam_structure"),
      # value-level (optional)
      vl_var    = c("value_level_variable","vl_variable","vl_var"),
      vl_where  = c("where","where_clause","where_expression"),
      unit      = c("unit","units","aval_unit","analysis_unit")
    ),
    colmap_ct = list(
      codelist = c("codelist","ct_name","controlled_terminology"),
      value    = c("value","code","term_code"),
      decode   = c("decode","term","meaning","label","label_text")
    ),
    
    # Guess ADaM dataset structure if not given in ITEMS
    dataset_structure_fn = function(ds) {
      ds <- toupper(ds)
      if (ds %in% c("ADSL","ADSUB","ADSUBJ")) "SPECIAL-PURPOSE"
      else if (str_detect(ds, "^ADTTE")) "BDS-TTE"
      else "BDS"
    },
    
    quiet = FALSE
) {
  stopifnot(file.exists(items_csv))
  if (!file.exists(dirname(out_xlsx))) dir.create(dirname(out_xlsx), recursive = TRUE)
  
  pick <- function(role, df, cmap) {
    hits <- cmap[[role]]
    hit  <- hits[hits %in% names(df)]
    if (length(hit)) hit[1] else NA_character_
  }
  
  # ---- Read & normalize ITEMS ----
  items_raw <- readr::read_csv(items_csv, show_col_types = FALSE) %>% clean_names()
  if (!quiet) cat("\nColumns in ITEMS:", paste(names(items_raw), collapse=", "), "\n")
  
  i_dataset  <- pick("dataset",  items_raw, colmap_items)
  i_variable <- pick("variable", items_raw, colmap_items)
  i_label    <- pick("label",    items_raw, colmap_items)
  i_type     <- pick("type",     items_raw, colmap_items)
  i_length   <- pick("length",   items_raw, colmap_items)
  i_origin   <- pick("origin",   items_raw, colmap_items)
  i_codelist <- pick("codelist", items_raw, colmap_items)
  i_key      <- pick("key",      items_raw, colmap_items)
  i_struct   <- pick("structure",items_raw, colmap_items)
  i_vl_var   <- pick("vl_var",   items_raw, colmap_items)
  i_vl_where <- pick("vl_where", items_raw, colmap_items)
  i_unit     <- pick("unit",     items_raw, colmap_items)
  
  need <- c(i_dataset, i_variable, i_label, i_type)
  if (any(is.na(need))) {
    stop("ITEMS CSV is missing essentials (dataset/variable/label/type). ",
         "Adjust `colmap_items` to your headers.")
  }
  
  items <- items_raw %>%
    transmute(
      dataset   = toupper(as.character(.data[[i_dataset]])),
      variable  = toupper(as.character(.data[[i_variable]])),
      label     = as.character(.data[[i_label]]),
      type      = as.character(.data[[i_type]]),
      length    = if (!is.na(i_length)) suppressWarnings(as.integer(.data[[i_length]])) else NA_integer_,
      origin    = if (!is.na(i_origin)) as.character(.data[[i_origin]]) else NA_character_,
      codelist  = if (!is.na(i_codelist)) as.character(.data[[i_codelist]]) else NA_character_,
      key       = if (!is.na(i_key)) .data[[i_key]] else NA,
      structure = if (!is.na(i_struct)) as.character(.data[[i_struct]]) else NA_character_,
      vl_var    = if (!is.na(i_vl_var)) as.character(.data[[i_vl_var]]) else NA_character_,
      vl_where  = if (!is.na(i_vl_where)) as.character(.data[[i_vl_where]]) else NA_character_,
      unit      = if (!is.na(i_unit)) as.character(.data[[i_unit]]) else NA_character_
    )
  
  # ---- Read & normalize CODELISTS ----
  if (file.exists(codelists_csv)) {
    ct_raw <- readr::read_csv(codelists_csv, show_col_types = FALSE) %>% clean_names()
    if (!quiet) cat("Columns in CODELISTS:", paste(names(ct_raw), collapse=", "), "\n")
    
    ct_name  <- pick("codelist", ct_raw, colmap_ct)
    ct_value <- pick("value",    ct_raw, colmap_ct)
    ct_decode<- pick("decode",   ct_raw, colmap_ct)
    
    if (any(is.na(c(ct_name, ct_value)))) {
      warning("Codelists CSV present but missing a required column (codelist/value). ",
              "Only the Variables/Datasets sheets will be written.")
      codelists <- tibble(Codelist=character(), Value=character(), Decode=character())
    } else {
      codelists <- ct_raw %>%
        transmute(
          Codelist = as.character(.data[[ct_name]]),
          Value    = as.character(.data[[ct_value]]),
          Decode   = as.character((.data[[ct_decode]] %||% .data[[ct_value]]))
        ) %>%
        filter(!is.na(Codelist), !is.na(Value), Codelist != "", Value != "") %>%
        arrange(Codelist, Value) %>%
        distinct()
    }
  } else {
    if (!quiet) message("No CODELISTS CSV found; Codelists sheet will be empty.")
    codelists <- tibble(Codelist=character(), Value=character(), Decode=character())
  }
  
  # ---- Type normalization (CDISC-ish) ----
  normalize_type <- function(x) {
    x <- str_to_lower(as.character(x))
    case_when(
      str_detect(x, "char|str|text|factor") ~ "char",
      str_detect(x, "date")                 ~ "date",
      str_detect(x, "time|datetime|posix")  ~ "datetime",
      str_detect(x, "int|num|dbl|float|dec")~ "num",
      TRUE                                  ~ "char"
    )
  }
  items <- items %>% mutate(type = normalize_type(type))
  
  # ---- Datasets sheet ----
  datasets <- items %>%
    distinct(dataset) %>%
    mutate(
      Structure = ifelse(!is.na(items$structure[match(dataset, items$dataset)]),
                         items$structure[match(dataset, items$dataset)],
                         map_chr(dataset, dataset_structure_fn)),
      Keys  = NA_character_,
      Label = paste0(dataset, " Dataset")
    ) %>%
    transmute(Dataset = dataset, Structure, Keys, Label)
  
  # Use key flags (if present) to populate Keys per dataset
  if (!all(is.na(items$key))) {
    key_tbl <- items %>%
      filter(!is.na(key) & (key %in% c(TRUE,1,"1") | str_to_lower(as.character(key)) %in% c("y","yes","key"))) %>%
      group_by(dataset) %>%
      summarise(Keys = paste(variable, collapse = ", "), .groups = "drop")
    datasets <- datasets %>%
      left_join(key_tbl, by = c("Dataset" = "dataset")) %>%
      mutate(Keys = coalesce(Keys.y, Keys.x)) %>%
      select(-Keys.x, -Keys.y)
  }
  
  # ---- Variables sheet ----
  variables <- items %>%
    transmute(
      Dataset  = dataset,
      Variable = variable,
      Type     = type,
      Length   = length,
      Label    = label,
      Origin   = origin,
      Codelist = na_if(codelist, "")
    ) %>%
    arrange(Dataset, Variable)
  
  # ---- ValueLevel sheet (optional)
  # We take explicit value-level hints if provided (vl_var / vl_where).
  # You can add rows like: Dataset=ADLB, Variable=AVALU, Where="PARAMCD == 'CRP'", Label="CRP Unit"
  valuelevel <- items %>%
    filter(!is.na(vl_var) | !is.na(vl_where)) %>%
    transmute(
      Dataset = dataset,
      Variable = ifelse(is.na(vl_var) | vl_var == "", NA_character_, vl_var),
      Where   = ifelse(is.na(vl_where) | vl_where == "", NA_character_, vl_where),
      Value   = NA_character_,
      Label   = ifelse(!is.na(unit) & unit != "", paste0("Unit: ", unit), NA_character_),
      Type    = NA_character_,
      Length  = NA_integer_
    ) %>%
    distinct()
  
  # ---- Write Excel ----
  writexl::write_xlsx(
    list(
      Datasets  = datasets,
      Variables = variables,
      Codelists = codelists,
      ValueLevel= valuelevel
    ),
    path = out_xlsx
  )
  
  if (!quiet) {
    cat("\nWrote:", out_xlsx, "\n")
    cat("\nDatasets:\n"); print(datasets, n=Inf)
    cat("\nVariables (head):\n"); print(head(variables, 20))
    if (nrow(codelists)) { cat("\nCodelists (head):\n"); print(head(codelists, 20)) }
    if (nrow(valuelevel)) { cat("\nValueLevel (head):\n"); print(head(valuelevel, 10)) }
  }
  
  invisible(list(Datasets=datasets, Variables=variables, Codelists=codelists, ValueLevel=valuelevel))
}

# =============================================================================
# EXAMPLE RUN
# =============================================================================
# items_codelists_to_metacore(
#   items_csv  = "metadata/_20250729_133232_Items.csv",
#   codelists_csv = "metadata/_20250729_133232_CodeLists.csv",
#   out_xlsx   = "metadata/metacore_map.xlsx",
#   quiet      = FALSE
# )
