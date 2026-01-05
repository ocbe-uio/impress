
make_rdbl <- function(adbl, cfg) {
  #' Make analysis dataset for baseline demographics and disease characteristics
  #'
  #' @param adbl ADLB dataset
  #' @param cfg Configuration list
  #'
  #' @return Analysis dataset for baseline demographics and disease characteristics
  #' @import dplyr
  #' @import tidyr
  #' @importFrom labelled set_variable_labels get_variable_labels
  #' @importFrom rtables basic_table analyze_vars build_table
  #' @export


  adbl_base <- adbl |>
    filter(cohortcd == 2) # Newly diagnosed glioblastoma cohort

  num_long <- adbl_base %>% filter(!is.na(aval))
  cat_long <- adbl_base %>% filter(!is.na(avalc))

  # --- 4) Pivot WIDE
  num_wide <- num_long %>%
    select(usubjid, arm, paramcd, aval) %>%
    pivot_wider(names_from = paramcd, values_from = aval)

  cat_wide <- cat_long %>%
    select(usubjid, arm, paramcd, avalc) %>%
    pivot_wider(names_from = paramcd, values_from = avalc) |>
    # convert into factors
    mutate(across(-c(usubjid, arm), ~ factor(.x)))

  # --- 5) Merge (keep only cols that exist)
  wide <- num_wide %>%
    left_join(cat_wide, by = c("usubjid", "arm"), suffix = c("cd", "")) |>
    select(usubjid, arm, tidyselect::any_of(cfg$reporting$baseline_order))

  # --- 6) Attach variable labels from PARAM (to both num and cat columns)
  get_labels <- function(df_long) {
    if (!nrow(df_long)) {
      return(character())
    }
    unique(df_long[c("paramcd", "param")]) # named vector: names=paramcd, values=param
  }
  var_labels <- get_labels(num_long) |>
    add_row(get_labels(cat_long)) |>
    group_by(paramcd) |>
    mutate(
      paramcd = if_else((row_number() > 1), paste0(paramcd, "cd"), paramcd),
      param = if_else((row_number() > 1), paste0(param, " (numeric)"), param)
    ) |>
    ungroup()

  var_labels <- stats::setNames(as.list(var_labels$param), var_labels$paramcd)

  wide <- wide %>%
    set_variable_labels(.labels = var_labels)

  var_labels <- get_variable_labels(wide, unlist = TRUE)
  var_labels <- var_labels[var_labels != ""]




  # Optional: make categoricals factors with desired level order
  # wide[cat_vars] <- lapply(wide[cat_vars], \(x) factor(x, levels = c("Female","Male"))) # example

  # --- 8) One layout, two analyze_colvars() calls (cont then cat)
  lyt <- basic_table() %>%
    # split_cols_by("arm") %>%
    analyze_vars(names(var_labels), var_labels = var_labels) # n (% ) per level

  tbl <- build_table(lyt, wide)
  main_title(tbl) <- "Demographics and Baseline characteristics"
  return(tbl)
}
