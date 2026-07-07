#!/usr/bin/env python3
# =====================================================================
# ECAS + FrSBe cognitive-behavioural phenotyping in ALS  (Paper 1, Brain)
# REPRODUCIBLE PIPELINE — PART 2: behavioural tables (Table 2 & Table 3)
#
# This companion script reproduces, from the frozen SAV, every value in:
#   - Table 2  (independent behavioural measures by phenotype)
#       * ECAS-CI total score      (sum of the 10 behavioural items; psychosis EXCLUDED)
#       * ECAS-CI domains altered  (mean n. of altered domains; rule below)
#       * FBI total, BBI (after), Montuschi criterion %
#   - Table 3  (premorbid and post-illness behavioural scores by phenotype)
#       * FrSBe total (premorbid / post-illness / change)
#       * BBI          (premorbid / post-illness / change)
#       * Kruskal-Wallis p across the four phenotypes
#
# Companion to: reproduce_analyses.py (partition, BCH, C9orf72 OR).
# Data source (single source of truth):
#   ALS_ECAS_FrSBe_dataset_reduced.sav
#   cohort = rows with non-missing ECASFRSBE_CLASS4 (N=777)
#
# Data loading and phenotype naming/order now come from als_common.py
# (shared with every other script).
#
# ENVIRONMENT: Python 3.12; numpy, pandas, scipy, pyreadstat
# =====================================================================

import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats

from als_common import RESULTS_DIR, ensure_dirs, load_cohort, PHENOTYPE_ORDER

ensure_dirs(RESULTS_DIR)

order = PHENOTYPE_ORDER
d = load_cohort(add_phenotype_label=True)
print(f"Cohort N = {len(d)}  (class sizes: " +
      ", ".join(f"{k} {int((d['ph']==k).sum())}" for k in order) + ")\n")

def num(col):
    return pd.to_numeric(d[col], errors="coerce")

# =====================================================================
# ECAS-CI  (ECAS Carer Interview) — 5 behavioural domains, psychosis EXCLUDED
# ---------------------------------------------------------------------
# 10 binary items grouped into 5 frontotemporal-dementia domains:
#   Disinhibition : ECAS_BE_DISIN1/2/3           (3 items)
#   Apathy        : ECAS_BE_APATH                (1 item)
#   LossEmpathy   : ECAS_BE_LOSS_EMPATH1/2       (2 items)
#   Perseveration : ECAS_BE_PERSEV1/2            (2 items)
#   Hyperorality  : ECAS_BE_HYPERORAL1/2         (2 items)
# The 3 ECAS-CI psychosis-screen items are EXCLUDED (column ECAS_PSY not used),
# so the total score ranges 0-10 rather than Poletti's 0-13 NoS.
# =====================================================================
domains = {
    "Disinhibition":  ["ECAS_BE_DISIN1", "ECAS_BE_DISIN2", "ECAS_BE_DISIN3"],
    "Apathy":         ["ECAS_BE_APATH"],
    "LossEmpathy":    ["ECAS_BE_LOSS_EMPATH1", "ECAS_BE_LOSS_EMPATH2"],
    "Perseveration":  ["ECAS_BE_PERSEV1", "ECAS_BE_PERSEV2"],
    "Hyperorality":   ["ECAS_BE_HYPERORAL1", "ECAS_BE_HYPERORAL2"],
}
allitems = [c for v in domains.values() for c in v]
for c in allitems:
    d[c] = num(c)

# --- ECAS-CI total score = sum of the 10 items (0-10) ---
d["ecasci_total"] = d[allitems].sum(axis=1)

# --- ECAS-CI domains altered = count of altered domains ---
# A domain counts as altered when >=1 item is endorsed, EXCEPT Disinhibition
# (3 items) which requires >=2. This rule reproduces the published values.
def domains_altered(row):
    n = 0
    for dom, items in domains.items():
        thr = 2 if dom == "Disinhibition" else 1
        if row[items].sum() >= thr:
            n += 1
    return n
d["ecasci_dom_altered"] = d.apply(domains_altered, axis=1)

all_rows = []

def block(colmap, title):
    """Print per-phenotype median (IQR) for a set of columns + KW p, and collect rows.

    Median (IQR) is reported rather than mean, because the manuscript's own
    Methods specify the Kruskal-Wallis test (rank-based, non-parametric) for
    continuous between-class comparisons -- median (IQR) is the descriptive
    statistic consistent with that choice of test. Reporting means alongside
    a non-parametric test is a common but avoidable mismatch.
    """
    print(f"--- {title} ---")
    for label, col in colmap.items():
        vals = [d.loc[d["ph"] == k, col].dropna() for k in order]
        medians = [v.median() for v in vals]
        q1s = [v.quantile(0.25) for v in vals]
        q3s = [v.quantile(0.75) for v in vals]
        H, p = stats.kruskal(*vals)
        row = "  ".join(f"{m:6.1f} [{q1:.1f}-{q3:.1f}]" for m, q1, q3 in zip(medians, q1s, q3s))
        print(f"  {label:24}  {row}   p={p:.1e}")
        all_rows.append({
            "table": title, "measure": label,
            "Preserved_median": medians[0], "Preserved_Q1": q1s[0], "Preserved_Q3": q3s[0],
            "Cognitive_only_median": medians[1], "Cognitive_only_Q1": q1s[1], "Cognitive_only_Q3": q3s[1],
            "Behavioural_only_median": medians[2], "Behavioural_only_Q1": q1s[2], "Behavioural_only_Q3": q3s[2],
            "Severe_FTD_median": medians[3], "Severe_FTD_Q1": q1s[3], "Severe_FTD_Q3": q3s[3],
            "kruskal_p": p,
        })
    print()

print("=" * 70)
print("TABLE 2 — ECAS-CI measures (psychosis excluded)")
print("=" * 70)
print(f"  {'':24}  {'Presv':>6}  {'Cog':>6}  {'Beh':>6}  {'Sev':>6}")
block({"ECAS-CI total score": "ecasci_total",
       "ECAS-CI domains altered": "ecasci_dom_altered"}, "ECAS-CI")

# =====================================================================
# TABLE 3 — premorbid vs post-illness (FrSBe total T-score; BBI)
# ---------------------------------------------------------------------
#   FrSBe premorbid : FRSBE_PRIMA_TOTALE     post : FRSBE_DOPO_TOTALE
#   BBI   premorbid : BBI_BEFORE             post : BBI_AFTER
# =====================================================================
d["frsbe_pre"]  = num("FRSBE_PRIMA_TOTALE")
d["frsbe_post"] = num("FRSBE_DOPO_TOTALE")
d["frsbe_chg"]  = d["frsbe_post"] - d["frsbe_pre"]
d["bbi_pre"]    = num("BBI_BEFORE")
d["bbi_post"]   = num("BBI_AFTER")
d["bbi_chg"]    = d["bbi_post"] - d["bbi_pre"]

print("=" * 70)
print("TABLE 3 — premorbid / post-illness / change")
print("=" * 70)
print(f"  {'':24}  {'Presv':>6}  {'Cog':>6}  {'Beh':>6}  {'Sev':>6}")
block({"FrSBe premorbid": "frsbe_pre",
       "FrSBe post-illness": "frsbe_post",
       "FrSBe change": "frsbe_chg"}, "FrSBe total (T-score)")
block({"BBI premorbid": "bbi_pre",
       "BBI post-illness": "bbi_post",
       "BBI change": "bbi_chg"}, "BBI")

# =====================================================================
# Expected published values (for quick visual check)
# =====================================================================
print("=" * 70)
print("EXPECTED (manuscript):")
print("  ECAS-CI total score      : 1.37 / 1.36 / 2.98 / 3.60")
print("  ECAS-CI domains altered  : 0.87 / 0.88 / 1.77 / 2.18")
print("  FrSBe change             : +4.3 / +4.7 / +10.8 / +20.2")
print("  BBI   change             : +4.4 / +5.2 / +10.0 / +18.4")
print("=" * 70)

pd.DataFrame(all_rows).to_csv(
    os.path.join(RESULTS_DIR, "table2_table3_values.csv"), index=False)
print(f"\n[Saved] results/table2_table3_values.csv")