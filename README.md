# Impress

This repository is for code connected with the Impress trial.

## Running the pipeline after cloning

### 1. Add the raw data

The pipeline reads raw data from a single timestamped export folder from Viedoc,
placed under `data/raw/`. Either:

- copy the existing `_20260626_083406` export folder into `data/raw/`, so you end
  up with `data/raw/_20260626_083406/`, or
- download a fresh export from Viedoc and place it under `data/raw/` the same way.

If you use a new export, update the `export:` value in
[config/cfg.yml](config/cfg.yml) to match the new folder's name (e.g.
`export: "_20260626_083406"`), since that value tells the pipeline which
subfolder of `data/raw/` to read.

The biomarker Excel files used by the `adlb` target also need to be present
under `data/raw/biomarkers/`.

### 2. Install dependencies with renv

This project pins package versions with [renv](https://rstudio.github.io/renv/)
(R 4.6.0). From the project root in R:

```r
renv::restore()
```

This installs every package listed in `renv.lock`, including `targets` and
`tarchetypes`, which run the pipeline itself.

### 3. Run the pipeline with targets

The pipeline is defined in [_targets.R](_targets.R). To run it end-to-end:

```r
targets::tar_make()
```

This builds every target in dependency order — importing and cleaning the raw
data, building the ADaM/RDaM datasets, and rendering the final reports:
`reports/impress_statistical_analysis.docx`, `reports/impress_statistical_analysis.pdf`,
and the model-fit reports under `reports/model fit/`.

To build only a specific report (and its dependencies), pass its target name,
e.g. `targets::tar_make(report_pdf)`. To inspect the pipeline before running
it, use `targets::tar_visnetwork()`.
