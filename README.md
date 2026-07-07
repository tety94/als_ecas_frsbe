# ECAS + FrSBe cognitive-behavioural phenotyping in ALS — Reproducibility package

Manuscript: **"The cognitive–behavioural architecture of amyotrophic lateral sclerosis"**
(submitted to *Brain*). PARALS registry, University of Turin.

This package contains everything needed to reproduce, from the data, every
model-based and descriptive number reported in the manuscript's main tables,
plus a set of reviewer-requested sensitivity/robustness/survival analyses and
a script that assembles the manuscript tables as a Word document.

---

## Contents

| File | What it is |
|------|-----------|
| `ALS_ECAS_FrSBe_dataset_reduced.sav` / `.csv` | Analysis dataset. Cohort only (N = 777), **31 variables** — just those used by the scripts. Either format works; `als_common.load_cohort()` tries `.sav` first, then falls back to `.csv`. |
| `DATA_DICTIONARY.csv` | One row per variable: name, non-missing count, plain-language description. |
| `als_common.py` | **Shared module** — single source of truth for the 7 clustering-indicator columns and FrSBe sign-reversal rule, the SAV/CSV cohort loader, the standardised clustering-matrix builder, Hungarian label remapping, phenotype names/order/colors, and p-value formatting. Imported by every script below; nothing in this list re-implements these pieces independently. |
| `reproduce_analyses.py` | **Script 1** — the partition, the BCH weights, and the C9orf72 odds ratios. |
| `reproduce_behaviural_tables.py` | **Script 2** — Table 2 (ECAS-CI) and Table 3 (FrSBe & BBI, premorbid → post-illness). *(Filename as shipped; note the spelling.)* |
| `sensitivity_analyses.py` | **Script 3** — reviewer-requested robustness checks: diagonal vs full-covariance model comparison, StepMix vs k-means, missingness on the clustering indicators. |
| `survival_analysis.py` | **Script 4** — reviewer-requested Kaplan-Meier survival / tracheostomy-free survival by phenotype. **Not reported in the manuscript.** |
| `build_tables.py` | **Script 5** — reads the `results/*.csv` files produced by Scripts 1–2 and writes `results/tables/tables_2_3_4.docx`, a Word document with Table 2, Table 3, and Table 4 styled to match Table 1 of the manuscript. |
| `official_posteriors.npy` / `official_mask.npy` | Posterior responsibilities of the final StepMix model (777×4) and the aligning mask, kept in the package for reference. **Not currently loaded by any script** — `reproduce_analyses.py` and `sensitivity_analyses.py` both refit StepMix directly on the standardised matrix (`random_state=42` makes this deterministic) rather than reading these files back in. |

> The full working database has 555 columns; this reduced file keeps only the
> 31 variables the analyses actually use, to keep the verification unambiguous.
> No direct identifiers are included.

---

## Data source (single source of truth)

`ALS_ECAS_FrSBe_dataset_reduced.sav` or `.csv` — cohort of **N = 777** patients
(rows with a non-missing official partition `ECASFRSBE_CLASS4`). The official
4-class partition is stored in `ECASFRSBE_CLASS4` (1 = Preserved, 2 = Cognitive-only,
3 = Behavioural-only, 4 = Severe/FTD). All loading goes through
`als_common.load_cohort()`.

---

## How to run

```bash
pip install stepmix==3.0.0 scikit-learn statsmodels scipy numpy pandas pyreadstat \
            lifelines matplotlib python-docx

python reproduce_analyses.py            # Script 1: partition + BCH + C9orf72 OR
python reproduce_behaviural_tables.py   # Script 2: Table 2 & Table 3 values
python sensitivity_analyses.py          # Script 3: robustness checks
python survival_analysis.py             # Script 4 (optional): survival by phenotype
python build_tables.py                  # Script 5: assemble tables_2_3_4.docx
                                         #   (run AFTER Scripts 1 and 2 — it reads
                                         #    their results/*.csv output)
```

All scripts read whichever data file is present in the same folder (`.sav`
preferred, `.csv` as fallback — `pyreadstat` is only needed if you use the
`.sav`). Fixed seed `random_state = 42` throughout. Python 3.12.

`als_common.py` is a library, not a script — it is imported by all five
scripts above and is not run directly.

---

## Script 1 — `reproduce_analyses.py`

**Step 1 — the 4-class StepMix partition (N = 777).** Seven clustering indicators
(defined once in `als_common.COGZ` / `als_common.FRS`): four ECAS cognitive
domains as Poletti-2016 z-scores (`ECAS_LIN_Z_POLETTI`, `ECAS_FL_Z_POLETTI`,
`ECAS_EXEC_Z_POLETTI`, `ECAS_MEM_Z_POLETTI`; ECAS Visuospatial excluded a
priori) and three FrSBe post-illness subscales, sign-reversed as `-(T-50)/10`
(`FRSBE_DOPO_APATIA`, `FRSBE_DOPO_DISINIBIZIONE`, `FRSBE_DOPO_DISESECUTIVO`).
All seven z-standardised over the cohort (`als_common.build_clustering_matrix`).

```
StepMix(n_components=4, measurement='continuous', n_init=20, random_state=42, max_iter=200)
```

Reproduces `ECASFRSBE_CLASS4` with **adjusted Rand index = 1.000** and class sizes
258 / 229 / 191 / 99.

> **Two easy-to-miss details are decisive:** (i) `measurement` must be `'continuous'`
> (`gaussian_diag` drifts to ARI ≈ 0.82); (ii) the ECAS inputs must be the
> `*_Z_POLETTI` columns, not `*_PUNTZ`. Any other choice fails to reproduce the
> partition and changes every downstream estimate.

**Step 2 — BCH weights.** Standard Bolck–Croon–Hagenaars correction for
classification uncertainty; the weight matrix is `D⁻¹` (negative entries are
intrinsic to BCH, not an error). Classification quality (diagonal of D):
0.94 / 0.89 / 0.89 / 0.88.

**Step 3 — BCH-weighted logistic regression (custom MLE).** Because statsmodels
rejects negative frequency weights, the odds ratio is obtained by maximising the
BCH-weighted binomial log-likelihood directly (BFGS), with Wald 95% CI / p-values.

**Headline results reproduced:**
- Partition 258 / 229 / 191 / 99, ARI 1.000.
- C9orf72, Severe/FTD vs Preserved (BCH): **OR 2.56** (95% CI 1.05–6.21, p = 0.038).
- C9orf72, combined-impaired vs Preserved (BCH): OR 1.41.
- Genotyped N = 710, 42 C9orf72 carriers.

**Saves:** `results/c9orf72_odds_ratios.csv`, `results/c9orf72_prevalence_by_class.csv`,
`results/partition_class_sizes.csv`.

---

## Script 2 — `reproduce_behaviural_tables.py`

Reproduces the descriptive behavioural tables directly from the dataset.

**Table 2 — ECAS-CI measures (psychosis screen excluded).** The ECAS Carer
Interview here uses the five core frontotemporal-dementia domains (10 binary items;
the 3 psychosis-screen items are excluded, so the total ranges 0–10, not Poletti's
0–13 NoS):
- **ECAS-CI total score** = sum of the 10 items → 1.37 / 1.36 / 2.98 / 3.60.
- **ECAS-CI domains altered** = mean number of altered domains, a domain being
  altered when ≥ 1 item is endorsed, except Disinhibition (3 items) which requires
  ≥ 2 → 0.87 / 0.88 / 1.77 / 2.18.

> Table 2 as currently implemented covers only these two ECAS-CI measures.
> FBI total and the Montuschi criterion percentage are not computed by this
> script; if the manuscript's Table 2 reports them, they need a separate
> pipeline (or an addition to this script) before they can be checked here.

**Table 3 — premorbid → post-illness (FrSBe total T-score; BBI).**
- FrSBe total change: +4.3 / +4.7 / +10.8 / +20.2 T-points.
- BBI change: +4.4 / +5.2 / +10.0 / +18.4.
- All Kruskal–Wallis p < 0.001 across the four phenotypes (exact p printed by the script).

Descriptives are reported as **median [IQR]**, matching the Kruskal-Wallis test
used for the between-class comparison.

**Saves:** `results/table2_table3_values.csv`.

---

## Script 3 — `sensitivity_analyses.py`

Reviewer-requested robustness checks, run on the same 7 clustering indicators
used to build the official partition.

**(A) Diagonal-covariance StepMix vs full-covariance GMM.**
Fits both models on the same standardised matrix and reports BIC, log-likelihood,
and ARI of each against the official partition, plus the direct ARI between the
two solutions. Addresses the reviewer question of whether the diagonal-covariance
assumption (local independence within class) is doing real work, or whether a
full-covariance model gives a comparable or better fit.

**(B) StepMix vs k-means.**
Reports the ARI of a plain k-means (k=4) solution against the official partition —
this is the number behind the manuscript's claim that StepMix was preferred over
k-means "with comparable reproducibility."

**(C) Missingness on the clustering indicators.**
Reports % missing per clustering variable within the N=777 extract, and (only if
any missingness is found) compares age/education between complete- and
incomplete-case patients. Because this reduced extract is, by construction,
already restricted to the 777 patients with a non-missing official partition,
this check is expected to show zero missingness here — in which case the script
explicitly flags that complete-case bias relative to the full eligible PARALS
population **cannot** be assessed from this file alone, and that Table S3 in the
manuscript (included vs eligible) is only a partial proxy for it.

**Saves:** `results/covariance_comparison_summary.csv`,
`results/covariance_comparison_labels.csv`, `results/kmeans_comparison_summary.csv`,
`results/missingness_clustering_vars.csv`.

---

## Script 4 — `survival_analysis.py` (reviewer-requested, not in the manuscript)

Kaplan-Meier overall survival and a tracheostomy-free composite endpoint
(death OR tracheostomy) by cognitive-behavioural phenotype, plus a
multivariate log-rank test across the four groups.

Uses `SURVIVAL`, `STATUS`, `TRACHEO_NEW`, columns kept in the reduced dataset
specifically for this reviewer-facing check. **Limitation flagged by the
script itself:** the composite endpoint reuses the single `SURVIVAL` time
column for both event types, because only one time-to-event column is
available in this reduced extract; if death and tracheostomy do not occur at
the same follow-up time in the source data, this is an approximation.

**Saves:** `results/figures/km_survival_by_phenotype.png`,
`results/figures/km_tracheofree_by_phenotype.png`,
`results/survival_logrank_summary.csv`, `results/survival_median_by_class.csv`.

---

## Script 5 — `build_tables.py`

Assembles Tables 2, 3, and 4 as a single Word document
(`results/tables/tables_2_3_4.docx`), styled to match Table 1 of the
manuscript (Times New Roman, light-blue header shading, thin grey borders).
Must be run **after** Scripts 1 and 2, since it reads their CSV outputs
(`table2_table3_values.csv`, `c9orf72_odds_ratios.csv`,
`c9orf72_prevalence_by_class.csv`) rather than recomputing anything.

Table 1 (cohort characteristics) and Table 5 (detailed genetics, including
SOD1/TARDBP) are **not** produced by this script: the underlying per-patient
values needed for them are not present in the reduced `results/*.csv` files
it reads.

---

## Environment

Python 3.12; `stepmix==3.0.0`, scikit-learn, statsmodels, scipy, numpy, pandas,
`pyreadstat` (optional, only needed for the `.sav` path), `lifelines` and
`matplotlib` (for `survival_analysis.py`), `python-docx` (for `build_tables.py`).
Fixed seed `random_state = 42`.

---

## Notes for the repository

- The BBI, FBI and premorbid FrSBe fields have some missingness by design
  (not every caregiver completed every instrument); the scripts handle this
  pairwise and report the N used for each value.
- `GENETICA` stores the gene symbol; only C9orf72, SOD1 and TARDBP are used in the
  manuscript, with a blank entry meaning non-carrier / negative screening.
- Survival columns (`SURVIVAL`, `STATUS`, `TRACHEO_NEW`) are included only for
  the reviewer-requested `survival_analysis.py`; survival is not reported in
  the manuscript itself.
- `als_common.py` is the single source of truth for the 7 clustering
  indicators, the cohort loader, phenotype naming/order/colors, and p-value
  formatting — if any of these ever need to change, change them there, not
  in the individual scripts.