########################################
# adbl dataset creattion
# the adbl datatset gathers all demographcs and baseline variables
# from the different souces in the raw data
# According to the protoco, the following information is gathered:
#
# 6.6.1 Before Randomization
# At first visit and/or when relevant, the information collected will include the patients’ personal identification number,
# general patient demographics, clinical status with relevant medical history and current medical conditions, relevant
# concomitant medications, diagnosis and extent of tumor, baseline tumor mutation status (MGMT and IDH for Study A),
# and details of prior anti-neoplastic treatments. In addition, all females of childbearing potential will have a serum pregnancy
# test during screening and before randomization. If positive, the patient is not eligible for the study.
# Clinical status and relevant medical history
# The patients’ clinical status is determined by assessing the following clinical information:
# • Weight with date of the exam
# • Clinical measures on neurologic status (KPS, ECOG, NANO) with date of the exam
# • Clinical information on disease progression with date (Study A – recurrent glioblastoma)
# • Medical history concerning or relevant to the treatment of the brain tumor
# • Any non-research MR exams related to the diagnosis of the brain tumor (section 6.6.2)
# • Blood pressure tests (section 6.6.2)
# • Blood test relevant for the treatment of the brain tumor (section 6.6.2)
# Physical examination
# A complete physical examination will be conducted by the study PI or designee at screening/baseline and will include
# examination of skin, lungs, heart, abdomen, lymph nodes, extremities, and neurological status. If indicated based on
# medical history and/or symptoms, further clinical tests will be performed. Significant findings that were present prior to the
# signing of informed consent must be included in the patient file and the eCRF. Height is measured in centimeters (cm) and
# body weight in kilogram [kg] (to the nearest 0.1 kg) in indoor clothing, but without shoes. Height information will be collected
# at the screening visit only.
# Neurologic performance status
# Assessment of KPS and ECOG performance status (APPENDIX D) and/or the Neurologic Assessment in Neuro-Oncology
# (NANO) scale98 (APPENDIX E) will be performed at screening.
# ImPRESS-losartan study Version no. 7.0 – 17-DEC-2024 revised Page 41 of 87
# Concomitant medication
# All treatment that the study PI or designee considers necessary for a participants’ welfare may be administered at the
# discretion of the investigators in keeping with the community standards of medical care. All concomitant medication will be
# recorded in the eCRF by the responsible research nurse. If changes occur during the study period, documentation of drug
# dosage, frequency, route, and date must also be included in the eCRF.
# Radiographic tumor evaluation
# Participants are evaluated for intracranial disease before treatment onset by the study PI and with the aid the study
# neuroradiologist at will. Patients may have measurable (radiographic) and non-measurable disease (clinical). Traditional
# radiographic response and progression will be evaluated in this study using the modified Response Assessment in Neuro-
# Oncology (RANO) criteria99,100. A trained neuroradiologist blinded to clinical and study data will outline the contrast agent
# enhancing lesion on the MR images as well as the surrounding area of FLAIR hyperintensity for assessment of intracranial
# edema. Lesion volumes are estimated using a volumetric approach101 that summarize pathological image voxels. The
# radiographic tumor evaluation will be recorded in the eCRF by the study PI or the responsible research nurse.


# Specifications of the data including the variable names, labels, types, derivations, and source are given in the metadata/items.csv file.

make_adbl <- function(raw, cfg, adsl, adeff) {

  resolve_visits <- function(cfg) {
    has_map <- function(x) !is.null(x$visits) && !is.null(x$visits$map)
    if (has_map(cfg)) {
      return(visits_from_yaml(cfg))
    }
    fallback_path <- "config/cfg.yml"
    if (fs::file_exists(fallback_path)) {
      fallback <- yaml::read_yaml(fallback_path)
      if (has_map(fallback)) {
        return(visits_from_yaml(fallback))
      }
    }
    stop("No 'visits: map:' section found in configuration.")
  }

  visits <- resolve_visits(cfg)

  subject_keys <- adsl %>%
    select(subjid, usubjid, studyid, cohort, cohortcd, arm, randdt)

  empty_adbl <- tibble::tibble(
    subjid = character(),
    avisit = character(),
    adt = as.Date(character()),
    parcat1 = character(),
    param = character(),
    paramcd = character(),
    aval = numeric(),
    avalc = character(),
    avalu = character()
  )

  label_for <- function(data, var) {
    lab <- attr(data[[var]], "label", exact = TRUE)
    if (is.null(lab) || !nzchar(lab)) var else as.character(lab)
  }

  safe_label <- function(data, var, default) {
    if (is.null(data)) {
      default
    } else {
      lab <- label_for(data, var)
      if (is.null(lab) || !nzchar(lab)) default else lab
    }
  }

  demographics <- dplyr::bind_rows(
    adsl %>%
      transmute(
        subjid,
        avisit = NA_character_,
        adt = as.Date(NA),
        parcat1 = "Demographics",
        param = "Age (years)",
        paramcd = "age",
        aval = as.numeric(age),
        avalc = NA_character_,
        avalu = "years"
      ),
    adsl %>%
      transmute(
        subjid,
        avisit = NA_character_,
        adt = as.Date(NA),
        parcat1 = "Demographics",
        param = "Sex",
        paramcd = "sex",
        aval = NA_real_,
        avalc = as.character(sex),
        avalu = NA_character_
      ),
    adsl %>%
      transmute(
        subjid,
        avisit = NA_character_,
        adt = as.Date(NA),
        parcat1 = "Demographics",
        param = "Cohort",
        paramcd = "cohort",
        aval = NA_real_,
        avalc = as.character(cohort),
        avalu = NA_character_
      ),
    adsl %>%
      transmute(
        subjid,
        avisit = NA_character_,
        adt = as.Date(NA),
        parcat1 = "Demographics",
        param = "Country",
        paramcd = "country",
        aval = NA_real_,
        avalc = as.character(country),
        avalu = NA_character_
      )
  )

  dm <- raw |> get_raw("dm")
  dm_cbp <- if (!is.null(dm)) {
    dm %>%
      transmute(
        subjid = subjectid,
        avisit = NA_character_,
        adt = as.Date(NA),
        parcat1 = "Demographics",
        param = safe_label(dm, "dmcbp", "Childbearing Potential"),
        paramcd = "dmcbp",
        aval = NA_real_,
        avalc = as.character(dmcbp),
        avalu = NA_character_
      ) %>%
      filter(!is.na(subjid))
  } else {
    empty_adbl
  }

  vs <- raw |> get_raw("vs")
  vs_records <- if (!is.null(vs)) {
    vs_vars <- c("vsweight", "vsheight", "vsbmi", "vssys", "vsdia", "vspulse")
    vs_map <- tibble::tibble(
      source = vs_vars,
      paramcd = c("weight", "height", "bmi", "sbp", "dbp", "pulse"),
      avalu = c("kg", "cm", "kg/m2", "mmHg", "mmHg", "bpm"),
      param = purrr::map_chr(vs_vars, ~ label_for(vs, .x))
    ) %>%
      mutate(param = dplyr::if_else(is.na(param) | !nzchar(param),
                                    stringr::str_to_title(source),
                                    param))

    vs %>%
      mutate(
        subjid = subjectid,
        adt = as.Date(vsdat),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      select(subjid, avisit, adt, dplyr::all_of(vs_vars)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(vs_vars),
        names_to = "source",
        values_to = "aval"
      ) %>%
      mutate(
        aval = suppressWarnings(as.numeric(aval))
      ) %>%
      filter(!is.na(aval), !is.na(adt)) %>%
      arrange(subjid, source, adt) %>%
      group_by(subjid, source) %>%
      slice(1) %>%
      ungroup() %>%
      left_join(vs_map, by = "source") %>%
      mutate(
        parcat1 = "Vital Signs",
        avalc = NA_character_
      ) %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1,
        param,
        paramcd,
        aval,
        avalc,
        avalu
      )
  } else {
    empty_adbl
  }

  ecog <- raw |> get_raw("ecog")
  ecog_records <- if (!is.null(ecog)) {
    ecog %>%
      mutate(
        subjid = subjectid,
        adt = as.Date(eventdate),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      filter(!is.na(ecogscd)) %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1 = "Neurologic Status",
        param = safe_label(ecog, "ecogs", "ECOG Score"),
        paramcd = "ecog",
        aval = as.numeric(ecogscd),
        avalc = paste0("ECOG Score", as.character(ecogscd)),
        avalu = NA_character_
      )
  } else {
    empty_adbl
  }

  kps <- raw |> get_raw("kps")
  kps_records <- if (!is.null(kps)) {
    kps %>%
      mutate(
        subjid = subjectid,
        adt = as.Date(eventdate),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      filter(!is.na(kpsscd)) %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1 = "Neurologic Status",
        param = safe_label(kps, "kpss", "Karnofsky Performance Scale"),
        paramcd = "kpss",
        aval = as.numeric(kpsscd),
        avalc = paste0("KPS ", as.character(kpsscd)),
        avalu = NA_character_
      )
  } else {
    empty_adbl
  }

  nano <- raw |> get_raw("nano")
  nano_records <- if (!is.null(nano)) {
    nano %>%
      mutate(
        subjid = subjectid,
        adt = as.Date(eventdate),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      filter(!is.na(nanotot)) %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1 = "Neurologic Status",
        param = safe_label(nano, "nanotot", "Total NANO score"),
        paramcd = "nano_tot",
        aval = as.numeric(nanotot),
        avalc = NA_character_,
        avalu = NA_character_
      )
  } else {
    empty_adbl
  }

  ds <- raw |> get_raw("ds")
  ds_records <- if (!is.null(ds)) {
    base_ds <- ds %>%
      mutate(
        subjid = subjectid,
        adt_record = as.Date(eventdate),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      arrange(subjid, adt_record) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup()

    base_ds <- base_ds %>%
      left_join(subject_keys %>% select(subjid, randdt), by = "subjid")

    dplyr::bind_rows(
      base_ds %>%
        filter(!is.na(dsdia), nzchar(as.character(dsdia))) %>%
        transmute(
          subjid,
          avisit,
          adt = adt_record,
          parcat1 = "Disease History",
          param = safe_label(ds, "dsdia", "Diagnosis"),
          paramcd = "diagnosis",
          aval = NA_real_,
          avalc = as.character(dsdia),
          avalu = NA_character_
        ),
      base_ds %>%
        mutate(
          ds_date = suppressWarnings(as.Date(dsardat)),
          randdt = as.Date(randdt),
          aval = as.numeric(randdt - ds_date)
        ) %>%
        filter(!is.na(ds_date), !is.na(randdt)) %>%
        transmute(
          subjid,
          avisit,
          adt = randdt,
          parcat1 = "Disease History",
          param = "Time from disease progression to randomisation (days)",
          paramcd = "disprog_rand",
          aval = aval,
          avalc = NA_character_,
          avalu = "days"
        ),
      base_ds %>%
        mutate(adt = suppressWarnings(as.Date(dsurgdat))) %>%
        mutate(
          randdt = as.Date(randdt),
          aval = as.numeric(randdt - adt)
        ) %>%
        filter(!is.na(adt), !is.na(randdt)) %>%
        transmute(
          subjid,
          avisit,
          adt = randdt,
          parcat1 = "Disease History",
          param = "Time from surgery to randomisation (days)",
          paramcd = "surg_rand",
          aval = aval,
          avalc = NA_character_,
          avalu = "days"
        ),
      base_ds %>%
        filter(!is.na(dsidh), nzchar(as.character(dsidh))) %>%
        transmute(
          subjid,
          avisit,
          adt = adt_record,
          parcat1 = "Disease History",
          param = safe_label(ds, "dsidh", "IDH1/2 mutation"),
          paramcd = "idh",
          aval = NA_real_,
          avalc = as.character(dsidh),
          avalu = NA_character_
        ),
      base_ds %>%
        filter(!is.na(dsmgmt), nzchar(as.character(dsmgmt))) %>%
        transmute(
          subjid,
          avisit,
          adt = adt_record,
          parcat1 = "Disease History",
          param = safe_label(ds, "dsmgmt", "MGMT status"),
          paramcd = "mgmt",
          aval = NA_real_,
          avalc = as.character(dsmgmt),
          avalu = NA_character_
        ),
      base_ds %>%
        filter(!is.na(dspctrea), nzchar(as.character(dspctrea))) %>%
        transmute(
          subjid,
          avisit,
          adt = adt_record,
          parcat1 = "Prior Treatment",
          param = safe_label(ds, "dspctrea", "Prior chemotherapy"),
          paramcd = "prchemo",
          aval = NA_real_,
          avalc = as.character(dspctrea),
          avalu = NA_character_
        ),
      base_ds %>%
        filter(!is.na(dspitre), nzchar(as.character(dspitre))) %>%
        transmute(
          subjid,
          avisit,
          adt = adt_record,
          parcat1 = "Prior Treatment",
          param = safe_label(ds, "dspitre", "Prior immunotherapy"),
          paramcd = "primmuno",
          aval = NA_real_,
          avalc = as.character(dspitre),
          avalu = NA_character_
        ),
      base_ds %>%
        filter(!is.na(dspotre), nzchar(as.character(dspotre))) %>%
        transmute(
          subjid,
          avisit,
          adt = adt_record,
          parcat1 = "Prior Treatment",
          param = safe_label(ds, "dspotre", "Other cancer therapy"),
          paramcd = "prother",
          aval = NA_real_,
          avalc = as.character(dspotre),
          avalu = NA_character_
        )
    )
  } else {
    empty_adbl
  }

  pt <- raw |> get_raw("pt")
  pt_records <- if (!is.null(pt)) {
    base_pt <- pt %>%
      mutate(
        subjid = subjectid,
        adt = as.Date(eventdate),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup()

    dplyr::bind_rows(
      base_pt %>%
        filter(!is.na(ptchemyn), nzchar(as.character(ptchemyn))) %>%
        transmute(
          subjid,
          avisit,
          adt,
          parcat1 = "Prior Treatment",
          param = safe_label(pt, "ptchemyn", "Received chemotherapy"),
          paramcd = "pt_chemo",
          aval = NA_real_,
          avalc = as.character(ptchemyn),
          avalu = NA_character_
        ),
      base_pt %>%
        filter(!is.na(ptimmyn), nzchar(as.character(ptimmyn))) %>%
        transmute(
          subjid,
          avisit,
          adt,
          parcat1 = "Prior Treatment",
          param = safe_label(pt, "ptimmyn", "Received immunotherapy"),
          paramcd = "pt_immuno",
          aval = NA_real_,
          avalc = as.character(ptimmyn),
          avalu = NA_character_
        ),
      base_pt %>%
        filter(!is.na(ptradyn), nzchar(as.character(ptradyn))) %>%
        transmute(
          subjid,
          avisit,
          adt,
          parcat1 = "Prior Treatment",
          param = safe_label(pt, "ptradyn", "Received radiation therapy to the brain"),
          paramcd = "pt_radiation",
          aval = NA_real_,
          avalc = as.character(ptradyn),
          avalu = NA_character_
        ),
      base_pt %>%
        filter(!is.na(ptcnsyn), nzchar(as.character(ptcnsyn))) %>%
        transmute(
          subjid,
          avisit,
          adt,
          parcat1 = "Prior Treatment",
          param = safe_label(pt, "ptcnsyn", "Received CNS surgery"),
          paramcd = "pt_cnsurgery",
          aval = NA_real_,
          avalc = as.character(ptcnsyn),
          avalu = NA_character_
        )
    )
  } else {
    empty_adbl
  }

  pre <- raw |> get_raw("pre")
  pre_records <- if (!is.null(pre)) {
    pre %>%
      mutate(
        subjid = subjectid,
        adt = as.Date(predat),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      mutate(pregtst_chr = toupper(as.character(pregtst))) %>%
      filter(pregtst_chr == "YES", !is.na(preres), nzchar(as.character(preres))) %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1 = "Pregnancy",
        param = safe_label(pre, "preres", "Pregnancy test result"),
        paramcd = "pregres",
        aval = NA_real_,
        avalc = as.character(preres),
        avalu = NA_character_
      )
  } else {
    empty_adbl
  }

  adeff_records <- adeff %>%
    filter(ablfl == "Y") %>%
    transmute(
      subjid,
      avisit,
      adt,
      parcat1 = "Imaging",
      param,
      paramcd,
      aval = as.numeric(aval),
      avalc = NA_character_,
      avalu = dplyr::case_when(
        stringr::str_to_lower(paramcd) %in% c("brain_volume", "tumor_volume") ~ "mm3",
        TRUE ~ NA_character_
      ),
      paramcd = stringr::str_to_lower(paramcd)
    )

  cs <- raw |> get_raw("cs")
  cs_records <- if (!is.null(cs)) {
    cs %>%
      mutate(
        subjid = subjectid,
        adt = as.Date(eventdate),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE)
      ) %>%
      filter(!is.na(avisit), avisit == "Randomization And Baseline MRI") %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1 = "Concomitant Therapy",
        param = safe_label(cs, "csmriyn", "Steroids at time of baseline MRI"),
        paramcd = "cs_steroids_mri",
        aval = NA_real_,
        avalc = as.character(csmriyn),
        avalu = NA_character_
      )
  } else {
    empty_adbl
  }

  qlq <- raw |> get_raw("qlq")
  qlq_records <- if (!is.null(qlq)) {
    qlq %>%
      mutate(
        subjid = subjectid,
        adt = suppressWarnings(as.Date(qlqc3dat)),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE),
        adt = dplyr::coalesce(adt, as.Date(eventdate))
      ) %>%
      filter(!is.na(qlqc3sc)) %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1 = "Quality of Life",
        param = safe_label(qlq, "qlqc3sc", "EORTC QLQ-C30 score"),
        paramcd = "qlq_c30",
        aval = as.numeric(qlqc3sc),
        avalc = NA_character_,
        avalu = NA_character_
      )
  } else {
    empty_adbl
  }

  qlqbn <- raw |> get_raw("qlqbn")
  qlqbn_records <- if (!is.null(qlqbn)) {
    qlqbn %>%
      mutate(
        subjid = subjectid,
        adt = suppressWarnings(as.Date(qlqbndat)),
        avisit = lookup_avisit(eventname, visits$map, fallback_title = TRUE),
        adt = dplyr::coalesce(adt, as.Date(eventdate))
      ) %>%
      filter(!is.na(qlqbnsc)) %>%
      arrange(subjid, adt) %>%
      group_by(subjid) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(
        subjid,
        avisit,
        adt,
        parcat1 = "Quality of Life",
        param = safe_label(qlqbn, "qlqbnsc", "EORTC QLQ-BN20 score"),
        paramcd = "qlq_bn20",
        aval = as.numeric(qlqbnsc),
        avalc = NA_character_,
        avalu = NA_character_
      )
  } else {
    empty_adbl
  }

  adbl <- dplyr::bind_rows(
    demographics,
    dm_cbp,
    vs_records,
    ecog_records,
    kps_records,
    nano_records,
    ds_records,
    pt_records,
    cs_records,
    qlq_records,
    qlqbn_records,
    pre_records,
    adeff_records
  ) %>%
    mutate(
      avalc = dplyr::na_if(avalc, "")
    ) %>%
    filter(!is.na(aval) | !is.na(avalc)) %>%
    left_join(subject_keys, by = "subjid") %>%
    mutate(
      studyid = dplyr::coalesce(studyid, cfg$studyid),
      domain = "ADBL"
    ) %>%
    select(
      studyid,
      domain,
      usubjid,
      subjid,
      cohort,
      cohortcd,
      randdt,
      arm,
      avisit,
      adt,
      parcat1,
      param,
      paramcd,
      aval,
      avalc,
      avalu
    ) %>%
    mutate(
      paramcd = stringr::str_to_lower(paramcd)
    ) %>%
    arrange(usubjid, parcat1, paramcd, adt)

  adbl
}
