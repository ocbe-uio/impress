library(targets)
tar_load(starts_with("ad"))
tar_load(primary_mri_sections)
tar_load(raw)
tar_load(cfg)

deaths <- adtte |>
  filter(cohortcd == 2) |> 
  filter(paramcd == "OS")


deaths |> 
  group_by(trt02p,evntstat) |> 
  summarise(n=n())


 adex_2 <- adex |> 
   filter(cohortcd == 2)


adex_2 |> filter(step != exranst) # None

adex_2 |> filter(exrandoscd != dose01p) #None
adex_2 |> filter(exfudoscd != dose02p) #None

tmp <- adex_2 |> filter(eximpdoscd != dose01p) #None


