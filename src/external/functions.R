

factoriser <- function(data, codelist = items, delabel = TRUE) {
  x <- names(data)
  y <- codelist %>%
    filter(id %in% x) %>%
    filter(categorical == 1) %>%
    select(id,value_labels)
  
  if(length(y$id) == 0) {
    return(data)
  }
  
  for (i in 1:length(y$id)){
    ct <- y %>%
      slice(i) %>%
      unnest(cols = c(value_labels)) %>% 
      unique()
    
    name1 <- y$id[[i]] 
    name2 <- paste0(y$id[[i]],"cd")
    
    
    if(all(ct[["datatype"]] == "integer")){
      
      
      labs <- as.numeric(ct[["codevalue"]])
      names(labs) <- ct[["codetext"]]

      data <- data %>%
      mutate_at(name2, as.numeric) %>%
      mutate_at(name2, haven::labelled, labels = labs) %>%
      mutate(!!name1 := as_factor(!!sym(name2), ordered = FALSE)) 
    }
    
    if(all(ct[["datatype"]] == "string")){
      
      labs <- ct[["codevalue"]]
      names(labs) <- ct[["codetext"]]
      

      data <- data %>%
        mutate_at(name2, as.character) %>%
        mutate_at(name2, haven::labelled, labels = labs) %>%
        mutate(!!name1 := as_factor(!!sym(name2), ordered = FALSE)) 
    }
      
      #if (delabel == TRUE && all(ct$datatype == "integer") ) data[[name2]] <- as.integer(data[[name2]])
      
  }
  
  return(data)
}



labeliser <- function(data, codelist = items){
  x <- names(data)
  labels <- codelist %>%
    filter(id %in% x) %>%
    select(id,label) %>%
    spread(id, label) %>%
    as.list()
  
  labelled::var_label(data) <- labels
  
  return(data)
}

pick <- function(db, name) {
  db %>% dplyr::filter(id == name) %>% purrr::pluck("data",1)
}

remove_fct <- function(data, codelist = items) {
  #Function to remove the factor and only retain the value labelled 
  #variable for export. 
  x <- names(data)
  cd1 <- codelist %>%
    filter(id %in% x) %>%
    filter(categorical == 1) %>%
    select(id) %>% 
    deframe()
  cd2 <- codelist %>%
    filter(id %in% x) %>%
    filter(categorical == 2) %>%
    select(id) %>% 
    deframe()
  
  if(length(cd1) == 0) {
    return(data)
  } else {
  data %>% 
    select(-all_of(cd1)) %>% 
    rename_with(~ str_sub(.x, end = -3), .cols = all_of(cd2))  %>% 
    return()
  }
}

remove_cd <- function(data, codelist = items) {
  #Function to remove the "cd" variables 
  x <- names(data)
  
  cd2 <- codelist %>%
    filter(id %in% x) %>%
    filter(categorical == 2) %>%
    select(id) %>% 
    deframe()
  
  if(length(cd2) == 0) {
    return(data)
  } else {
    data %>% 
      select(-all_of(cd2)) %>% 
      return()
  }
}

###################
# Functions for tables
##################
mean_sd <- function(data, var, group, digits = 1) {
  var <- ensym(var)
  group <- ensym(group)
  data %>% 
    group_by(!!group) %>% 
    summarise(mean = mean(!!var, na.rm = TRUE), 
              sd = sd(!!var, na.rm = TRUE), 
              missing = sum(is.na(!!var))
              , .groups = "drop_last") %>% 
    mutate_at(vars(mean, sd), ~round(., digits = digits)) %>% 
    mutate(txt = paste0(mean, " (", sd, ")")) %>% 
    select(group, txt) %>% 
    deframe
}

median_iqr <- function(data, var, group, digits = 1) {
  var <- ensym(var)
  group <- ensym(group)
  data %>% 
    group_by(!!group) %>% 
    summarise(median = median(!!var, na.rm = TRUE), 
              q1 = quantile(!!var, probs = 0.25, na.rm = TRUE), 
              q3 = quantile(!!var, probs = 0.75, na.rm = TRUE), 
              missing = sum(is.na(!!var))
              , .groups = "drop_last") %>% 
    mutate(across(c(median, q1, q3), ~round(.x, digits = digits))) %>% 
    mutate(txt = paste0(median, " (", q1, " - ", q3, ")")) %>% 
    select(group, txt) %>% 
    deframe
}

n_pct <-  function(data, var, group, level = 1) {
  var <- ensym(var)
  group <- ensym(group)
  data %>% 
    group_by(!!group, !!var, .drop = FALSE) %>% 
    summarise(n = n(),
              tot = n(), 
              .groups = "drop_last") %>% 
    filter(!is.na(!!var)) %>% 
    group_by(!!group, .drop = TRUE) %>% 
    mutate(tot = sum(tot),
           pct = round(n/tot*100, digits = 1)) %>% 
    mutate(txt = paste0(n, " (", pct, "%)")) %>% 
    filter(!!var == !!level & tot>0) %>%
    select(group, txt) %>%
    deframe
}

empty <- function(data, var, group, ...){
  group <- ensym(group)
  data %>% 
    group_by(!!group) %>% 
    summarise(n = n(), 
              .groups = "drop_last") %>% 
    mutate(txt = "") %>% 
    select(group,txt) %>% 
    deframe
}

missing_f <-  function(data, var, group, ...) {
  var <- ensym(var)
  group <- ensym(group)
  data %>% 
    group_by(!!group) %>% 
    summarise(tot = n(),
              non_miss = sum(!is.na(!!var)),
              miss = sum(is.na(!!var)),
              .groups = "drop_last") %>% 
    group_by(!!group) %>%
    mutate(pct = round(miss/tot*100,digits = 1)) %>%
    mutate(txt = paste0(miss, " (", pct, "%)")) %>%
    ungroup %>%
    select(group, txt) %>%
    deframe
}

ae_n_pct <-  function(data, var, group, level = 1) {
  var <- ensym(var)
  group <- ensym(group)
  
  data %>%
    group_by(subjectid, !!group, !!var) %>%
    summarise(n = sum(!!var),
              .groups = "drop_last") %>%
    group_by(!!group, !!var) %>%
    summarise(n_ae = sum(n),
              n_pat = n(),
              .groups = "drop_last") %>%
    group_by(!!group) %>%
    mutate(N_pat = sum(n_pat),
           pct = round(n_pat/N_pat*100,digits = 1),
           txt = paste0(n_pat, " (", pct, "%)")) %>%
    filter(!!var %in% !!level) %>%
    ungroup %>%
    select(!!group, txt) %>%
    deframe
}

ae_N_n_pct <-  function(data, var, group, level = 1) {
  var <- ensym(var)
  group <- ensym(group)
  
  data %>%
    group_by(subjectid, !!group) %>%
    summarise(n = sum(!!var),
              .groups = "drop_last") %>%
    mutate(!!var := if_else(n==0, 0, 1)) %>%
    group_by(!!group, !!var) %>%
    summarise(n_ae = sum(n),
              n_pat = n(),
              .groups = "drop_last") %>%
    group_by(!!group) %>%
    mutate(N_pat = sum(n_pat),
           pct = round(n_pat/N_pat*100,digits = 1),
           txt = paste0("[", n_ae,"] ", n_pat, " (", pct, "%)")) %>%
    mutate(txt = if_else(n_ae == 0, "[0] 0 (0%)", txt)) %>%
    filter(!!var %in% !!level) %>%
    ungroup %>%
    select(!!group, txt) %>%
    deframe
}

stats_exec <- function(f, data, var, group, ...){
  rlang::exec(f, data, var, group, !!!(...))
}


plot_cont_margins <- function(data, ytitle = "Value") {
  data %>%
    gf_line(
      margin ~ studyday,
      color = ~ rantrt,
      group = ~ rantrt,
      position = position_dodge(0.4),
      size = 1
    ) %>%
    gf_point(position = position_dodge(0.4)) %>%
    gf_errorbar(
      ci_lb + ci_ub ~ studyday,
      color = ~ rantrt,
      width = .8,
      position = position_dodge(0.4)
    ) %>%
    gf_labs(x = "Study day",
            y = str2expression(ytitle),
            color = "Treatment") %>%
    gf_theme(theme_classic())
  
}

