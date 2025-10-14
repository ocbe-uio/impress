tmp <- function(tbl) {
  tbl %>%
    mutate(data = list(adae),
           group = "rantrt",
           param = list(level = 1)) %>%
    mutate(res = pmap(list(f, data, var, group, param), rlang::exec)) %>%
    mutate(id = map(res,names)) %>%
    unnest(c(res, id)) %>%
    mutate(id = paste0("txt", id)) %>%
    pivot_wider(values_from = res, names_from = id) %>%
    select(text, starts_with("txt")) %>%
    kable(col.names = c("Parameter", header_ae),
          caption = "Summary of Adverse Events",
          booktabs = TRUE)
}

test <- tmp(ae_summary_table)
write_rds(test, "5.3 Analyser/Datasets/work/test.rds")