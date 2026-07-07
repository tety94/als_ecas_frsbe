# ECAS + FrSBe cognitive-behavioural phenotyping in ALS — Reproducibility package

Manuscript: **"The cognitive–behavioural architecture of amyotrophic lateral sclerosis"**
(submitted to *Brain*). PARALS registry, University of Turin.

This package contains everything needed to reproduce, from the data, every
model-based and descriptive number reported in the manuscript's main tables,
plus a set of reviewer-requested sensitivity/robustness analyses.

---

## Contents

| File | What it is |
|------|-----------|
| `ALS_ECAS_FrSBe_dataset_reduced.sav` / `.csv` | Analysis dataset. Cohort only (N = 777), **31 variables** — just those used by the scripts. Either format works; scripts try `.sav` first, then fall back to `.csv`. |
| `DATA_DICTIONARY.csv` | One row per variable: name, non-missing count, plain-language description. |
| `reproduce_analyses.py` | **Script 1** — the partition, the BCH weights, and the C9orf72 odds ratios. |
| `reproduce_behavioural_tables.py` | **Script 2** — Table 2 (ECAS-CI, FBI, BBI) and Table 3 (FrSBe & BBI premorbid → post-illness). |
| `sensitivity_analyses.py` | **Script 3 (new)** — reviewer-requested robustness checks: diagonal vs full-covariance model comparison, StepMix vs k-means, missingness on the clustering indicators. |
| `official_posteriors.npy` | Posterior responsibilities of the final StepMix model (777×4), so the BCH step can be reproduced without refitting. |
| `official_mask.npy` | Boolean mask aligning the posteriors to the dataset rows. |

> The full working database has 555 columns; this reduced file keeps only the
> 31 variables the analyses actually use, to keep the verification unambiguous.
> No direct identifiers are included.

---

## Data source (single source of truth)

`ALS_ECAS_FrSBe_dataset_reduced.sav` or `.csv` — cohort of **N = 777** patients
(rows with a non-missing official partition `ECASFRSBE_CLASS4`). The official
4-class partition is stored in `ECASFRSBE_CLASS4` (1 = Preserved, 2 = Cognitive-only,
3 = Behavioural-only, 4 = Severe/FTD).

---

## How to run

```bash
pip install stepmix==3.0.0 scikit-learn statsmodels scipy numpy pandas pyreadstat
python reproduce_analyses.py            # Script 1: partition + BCH + C9orf72 OR
python reproduce_behavioural_tables.py  # Script 2: Table 2 & Table 3 values
python sensitivity_analyses.py          # Script 3: robustness checks (new)
```

All three scripts read whichever data file is present in the same folder
(`.sav` preferred, `.csv` as fallback — `pyreadstat` is only needed if you use
the `.sav`). Fixed seed `random_state = 42` throughout. Python 3.12.

---

## Script 1 — `reproduce_analyses.py`

**Step 1 — the 4-class StepMix partition (N = 777).** Seven clustering indicators:
four ECAS cognitive domains as Poletti-2016 z-scores
(`ECAS_LIN_Z_POLETTI`, `ECAS_FL_Z_POLETTI`, `ECAS_EXEC_Z_POLETTI`, `ECAS_MEM_Z_POLETTI`;
ECAS Visuospatial excluded a priori) and three FrSBe post-illness subscales,
sign-reversed as `-(T-50)/10` (`FRSBE_DOPO_APATIA`, `FRSBE_DOPO_DISINIBIZIONE`,
`FRSBE_DOPO_DISESECUTIVO`). All seven z-standardised over the cohort.

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

---

## Script 2 — `reproduce_behavioural_tables.py`

Reproduces the descriptive behavioural tables directly from the dataset.

**Table 2 — ECAS-CI measures (psychosis screen excluded).** The ECAS Carer
Interview here uses the five core frontotemporal-dementia domains (10 binary items;
the 3 psychosis-screen items are excluded, so the total ranges 0–10, not Poletti's
0–13 NoS):
- **ECAS-CI total score** = sum of the 10 items → 1.37 / 1.36 / 2.98 / 3.60.
- **ECAS-CI domains altered** = mean number of altered domains, a domain being
  altered when ≥ 1 item is endorsed, except Disinhibition (3 items) which requires
  ≥ 2 → 0.87 / 0.88 / 1.77 / 2.18.

**Table 3 — premorbid → post-illness.**
- FrSBe total change: +4.3 / +4.7 / +10.8 / +20.2 T-points.
- BBI change: +4.4 / +5.2 / +10.0 / +18.4.
- All Kruskal–Wallis p < 0.001 across the four phenotypes (exact p printed by the script).

---

## Script 3 — `sensitivity_analyses.py` (new)

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
k-means "with comparable reproducibility," which was not previously quoted.

**(C) Missingness on the clustering indicators.**
Reports % missing per clustering variable within the N=777 extract, and (only if
any missingness is found) compares age/education between complete- and
incomplete-case patients. Because this reduced extract is, by construction,
already restricted to the 777 patients with a non-missing official partition,
this check is expected to show zero missingness here — in which case the script
explicitly flags that complete-case bias relative to the full eligible PARALS
population **cannot** be assessed from this file alone, and that Table S3 in the
manuscript (included vs eligible) is only a partial proxy for it.

---

## Environment

Python 3.12; `stepmix==3.0.0`, scikit-learn, statsmodels, scipy, numpy, pandas,
`pyreadstat` (optional, only needed for the `.sav` path). Fixed seed
`random_state = 42`.

---

## Notes for the repository

- The BBI, FBI and premorbid FrSBe fields have some missingness by design
  (not every caregiver completed every instrument); the scripts handle this
  pairwise and report the N used for each value.
- `GENETICA` stores the gene symbol; only C9orf72, SOD1 and TARDBP are used in the
  manuscript, with a blank entry meaning non-carrier / negative screening.
- Survival columns (`SURVIVAL`, `STATUS`, `TRACHEO_NEW`) are included only for
  reviewer-requested survival analyses; survival is not reported in the manuscript.