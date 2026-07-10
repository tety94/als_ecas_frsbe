# ECAS + FrSBe cognitive-behavioural phenotyping in ALS — Reproducibility package

Manuscript: **"The cognitive–behavioural architecture of amyotrophic lateral sclerosis"**. PARALS registry, University of Turin.

---

## Contents

| File | What it is |
|------|-----------|
| `als_common.py` | **Shared module** — single source of truth for the 7 clustering-indicator columns and FrSBe sign-reversal rule, the SAV/CSV cohort loader, the standardised clustering-matrix builder, Hungarian label remapping, phenotype names/order/colors, and p-value formatting. Imported by every script below; nothing in this list re-implements these pieces independently. |
| `reproduce_analyses.py` | **Script 1** — the partition, the BCH weights, and the C9orf72 odds ratios. |
| `reproduce_behaviural_tables.py` | **Script 2** — Table 2 (ECAS-CI) and Table 3 (FrSBe & BBI, premorbid → post-illness). *(Filename as shipped; note the spelling.)* |
| `survival_analysis.py` | **Script 4** — reviewer-requested Kaplan-Meier survival / tracheostomy-free survival by phenotype. **Not reported in the manuscript.** |
| `build_tables.py` | **Script 5** — reads the `results/*.csv` files produced by Scripts 1–2 and writes `results/tables/tables_2_3_4.docx`, a Word document with Table 2, Table 3, and Table 4 styled to match Table 1 of the manuscript. |

---

## How to run

```bash
pip install stepmix==3.0.0 scikit-learn statsmodels scipy numpy pandas pyreadstat \
            lifelines matplotlib python-docx

python reproduce_analyses.py            # Script 1: partition + BCH + C9orf72 OR
python reproduce_behaviural_tables.py   # Script 2: Table 2 & Table 3 values
python sensitivity_analyses.py          # Script 3: robustness checks
python build_tables.py                  # Script 4: assemble tables_2_3_4.docx
                                         #   (run AFTER Scripts 1 and 2 — it reads
                                         #    their results/*.csv output)
```

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
- `als_common.py` is the single source of truth for the 7 clustering
  indicators, the cohort loader, phenotype naming/order/colors, and p-value
  formatting — if any of these ever need to change, change them there, not
  in the individual scripts.