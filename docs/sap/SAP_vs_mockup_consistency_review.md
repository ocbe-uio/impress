# ImPRESS — Consistency review: SAP v0.2 vs. mock-up report

**Inputs**
- SAP: `docs/sap/ImPRESS AN SAP v0.2.docx`
- Mock-up report: `reports/impress_statistical_analysis.pdf` / `.Rmd`
- Code: `R/make_ad/*.R`, `R/make_rd/*.R`, `_targets.R`, `config/cfg.yml`

**Status:** All issues from the v0.1 and v0.2 reviews are resolved or withdrawn. No open inconsistencies remain. Two sham-data limitations are tracked for database lock; one AE design choice is now documented in SAP §7.2.

---

## Resolved

| Ref | Issue | Resolution |
|-----|-------|------------|
| A1 | Report named "logistic" as 4th candidate model | Report now "linear, piecewise linear, Emax, sigmoid Emax" — matches SAP + code |
| A2 | Co-primary hierarchical gatekeeping not in code | SAP §4.1 clarified: all hypothesis tests are performed; the hierarchical rule is applied **at interpretation** only. Now matches the pipeline (which computes all tests) |
| A3 | "Mod step selects best model (by AIC)" overstated | Report reworded: "the Mod step presents the model fit by AIC and reports ED values" — matches code (all models shown, sorted by AIC) and SAP ("presented descriptively") |
| A4 | Biomarkers out of SAP scope | Added as **Secondary Objective 3** (§1) + **Secondary Biomarker Endpoints** (§5); log scale pre-specified |
| A5 | Safety analyses absent from report | Fully implemented: `adae`/`advs`/`adex`/`adlbsaf` + report §Safety (exposure, AEs, vital-signs LMM, labs with reference-range LOW/NORMAL/HIGH, deaths). Matches SAP §7 |
| A6 | RANO: report 2 endpoints, SAP 1 | SAP §1 now specifies ordinal RANO endpoints "one for each timepoint and one best overall response" — maps to coded `trsdisea` (Days 141/225) + `trorresp` (Day 225) |
| B1 | Protocol date discrepancy | SAP now "17 December 2024" |
| C1 | "key secondary" MRI grouping label | Report now "primary and sensitivity (ROI2)" |
| C2 | KM/log-rank not in SAP | SAP §5.5 now lists Kaplan–Meier plots and overall log-rank as supportive displays |
| C3 | PH diagnostics not shown | SAP references scaled Schoenfeld residual tests/plots and log-log plots in a separate document (`impress_tte_ph_check.pdf`) |
| C4 | Day 29 neuro listed but not collected | SAP §5.3.1.1 now uses Days 15 and 43 only |
| C5 | AN mislabelled "recurrent" | SAP intro now "newly diagnosed glioblastoma (AN)" |
| C7 | Report intro framed as baseline-only review; subtitle typo | Contradicting Introduction/Randomisation-Status framing removed; subtitle typo corrected |

## Withdrawn (not an inconsistency)

| Ref | Note |
|-----|------|
| B2 / C6 | Visit-day "off by one" was a misread. Baseline = Day 1, so "154 days after initiation" = nominal Day 155 (likewise 140→141, 224→225, 238→239). Elapsed-days vs nominal-study-day are two consistent conventions. |

---

## Tracked for database lock (sham-data limitations, not inconsistencies)
- **Safety Set:** `_EX.csv` is empty in the current sham export, so SAF is approximated by the FAS and safety outputs are grouped by **planned** dose. Analysis by **treatment actually received** (per SAP §7) to be wired in at DBL/unblinding.
- **Conventional labs:** only creatinine is populated in the current export; SAP §7 envisages the full conventional panel. Confirm the panel populates at DBL.

## Documented design choice
- AE early/late phasing and dose-at-onset attribution are now described in **SAP §7.2 General Principles** (early window onset ≤ Day 44 by stepped-wedge dose at onset with non-mutually-exclusive 0 mg group; late window onset ≥ Day 45 by long-term dose). Matches `R/make_rd/make_rdae.R`.

---

*No remaining action items. Document reflects SAP v0.2, report, and code as of this review.*
