#!/usr/bin/env python3
# =====================================================================
# ECAS + FrSBe cognitive-behavioural phenotyping in ALS
# SHARED HELPERS — single source of truth for constants and routines
# duplicated (with drift risk) across the individual analysis scripts:
#   reproduce_analyses.py, reproduce_behavioural_tables.py,
#   sensitivity_analyses.py, survival_analysis.py, build_tables.py,
#   make_tables_figures.py
#
# Centralising these means:
#   - the 7 clustering-indicator columns and the sign-reversal rule for
#     FrSBe are defined ONCE, so the official partition (reproduce_analyses.py)
#     and its sensitivity checks (sensitivity_analyses.py) can never drift
#     apart on what "the clustering matrix" means;
#   - the SAV-with-CSV-fallback loader is used by every script,
#     including sensitivity_analyses.py, which previously read the CSV
#     only and would have silently broken (or silently used stale data)
#     if ever re-run against the frozen .sav;
#   - phenotype names/order/colors and p-value formatting are consistent
#     across every table, figure, and Word document produced.
#
# ENVIRONMENT: Python 3.12; pandas, numpy, scipy, pyreadstat (optional,
# only needed if the .sav file is present)
# =====================================================================

import os
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

# ---------------------------------------------------------------------
# Data file locations (single source of truth)
# ---------------------------------------------------------------------
SAV = "data/ALS_ECAS_FrSBe_dataset_reduced.sav"
CSV = "data/ALS_ECAS_FrSBe_dataset_reduced.csv"
RESULTS_DIR = "results"
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
TAB_DIR = os.path.join(RESULTS_DIR, "tables")

# ---------------------------------------------------------------------
# The 7 clustering indicators (official partition definition)
# ---------------------------------------------------------------------
COGZ = ["ECAS_LIN_Z_POLETTI", "ECAS_FL_Z_POLETTI",
        "ECAS_EXEC_Z_POLETTI", "ECAS_MEM_Z_POLETTI"]
FRS = ["FRSBE_DOPO_APATIA", "FRSBE_DOPO_DISINIBIZIONE",
       "FRSBE_DOPO_DISESECUTIVO"]

# ---------------------------------------------------------------------
# Phenotype naming / ordering / colors (used everywhere a class 1..4
# needs to become a human-readable label, in a fixed display order)
# ---------------------------------------------------------------------
PHENOTYPE_NAMES = {1: "Preserved", 2: "Cognitive-only",
                    3: "Behavioural-only", 4: "Severe/FTD"}
PHENOTYPE_ORDER = ["Preserved", "Cognitive-only", "Behavioural-only", "Severe/FTD"]
PHENOTYPE_COLORS = {"Preserved": "#4c72b0", "Cognitive-only": "#dd8452",
                     "Behavioural-only": "#55a868", "Severe/FTD": "#c44e52"}

# snake_case variants, needed where column names can't contain spaces/slashes
# (e.g. table2_table3_values.csv, built_tables.py's classes list)
PHENOTYPE_SNAKE = ["Preserved", "Cognitive_only", "Behavioural_only", "Severe_FTD"]


def ensure_dirs(*dirs):
    """Create the given directories (default: results/, figures/, tables/)."""
    if not dirs:
        dirs = (RESULTS_DIR,)
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def load_cohort(require_class4=True, add_phenotype_label=False):
    """Load the cohort from the frozen .sav if present, else fall back to
    the .csv extract. This is the ONE place that decides how the source
    data file is found and read; every script should call this rather
    than re-implementing the SAV/CSV fallback.

    Parameters
    ----------
    require_class4 : restrict to rows with a non-missing official
        4-class partition (ECASFRSBE_CLASS4), i.e. the N=777 analysis
        cohort. Set False only for scripts that need the full extract.
    add_phenotype_label : add a "ph" column mapping ECASFRSBE_CLASS4
        (1..4) to its display name via PHENOTYPE_NAMES.
    """
    if os.path.exists(SAV):
        import pyreadstat
        d, meta = pyreadstat.read_sav(SAV)
    elif os.path.exists(CSV):
        d = pd.read_csv(CSV)
    else:
        raise FileNotFoundError(f"Neither {SAV} nor {CSV} found in working directory.")

    if require_class4:
        d = d[d["ECASFRSBE_CLASS4"].notna()].copy()

    if add_phenotype_label:
        d["ph"] = d["ECASFRSBE_CLASS4"].astype(int).map(PHENOTYPE_NAMES)

    return d


def build_clustering_matrix(d):
    """Build the standardised 7-indicator clustering matrix (4 ECAS
    cognitive domains + 3 sign-reversed FrSBe post-illness subscales).

    Mutates `d` in place to coerce the 7 source columns to numeric
    (matches prior script behaviour). Returns (M, mask, Xraw) where:
      M    : standardised matrix, complete cases only (numpy array)
      mask : boolean Series aligned to d.index, True = complete case
      Xraw : the pre-standardisation DataFrame (all rows, for inspection)
    """
    for c in COGZ + FRS:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    X = pd.DataFrame(index=d.index)
    for c in COGZ:
        X[c] = d[c]
    for c in FRS:
        X[c] = -(d[c] - 50) / 10.0          # sign-reverse FrSBe
    Xz = (X - X.mean()) / X.std()           # standardise
    mask = Xz.notna().all(axis=1)
    return Xz[mask].values, mask, X


def remap_labels_to_official(lab, official, k=4):
    """Hungarian-match arbitrary cluster labels (0..k-1) onto the
    official 1..k class numbering, so ARI/BIC comparisons and any
    downstream inspection use the same class identities."""
    cost = np.zeros((k, k))
    for a in range(k):
        for b in range(k):
            cost[a, b] = -np.sum((lab == a) & (official == b + 1))
    ri, ci = linear_sum_assignment(cost)
    remap = {ri[i]: ci[i] + 1 for i in range(k)}
    return np.array([remap[x] for x in lab])


def fmt_p(p):
    """Consistent p-value formatting used across all tables/figures/docx:
    scientific notation below 0.001, 3 decimals otherwise."""
    return f"{p:.1e}" if p < 0.001 else f"{p:.3f}"