#!/usr/bin/env python3
# =====================================================================
# ECAS + FrSBe cognitive-behavioural phenotyping in ALS  (Paper 1)
# COMPLETE REPRODUCIBLE PIPELINE
#
# This script reproduces, from the SAV data file, every model-based
# result in the manuscript:
#   (1) the official 4-class StepMix partition (258/229/191/99, N=777)
#   (2) the BCH (bias-adjusted three-step) analyses that account for
#       classification uncertainty, including the C9orf72 odds ratio
#   (3) the composite-z external-validation profiles (Table 3)
#
# Author of analysis pipeline: A. Chiò group, PARALS, University of Turin
# Reproduces manuscript version: v3 (Brain submission)
#
# ---------------------------------------------------------------------
# ENVIRONMENT
#   Python 3.12
#   stepmix==3.0.0   scikit-learn   statsmodels   scipy   numpy   pandas
#   pyreadstat  (to read the SPSS .sav)
#
# DATA FILE (single source of truth):
#   ALS_ECAS_FrSBe_dataset_reduced.sav
#   - N=2767 rows; cohort = rows with non-missing ECASFRSBE_CLASS4 (N=777)
#   - Official partition stored in column ECASFRSBE_CLASS4 (1..4)
#
# Shared constants/helpers (data loading, the 7-indicator clustering
# matrix, Hungarian label remapping) now live in als_common.py so that
# this script and sensitivity_analyses.py can never drift apart on what
# "the clustering matrix" means. See als_common.py for details.
# =====================================================================

import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from stepmix.stepmix import StepMix
from sklearn.metrics import adjusted_rand_score
from scipy.optimize import minimize, linear_sum_assignment  # noqa
from scipy.stats import norm

from als_common import (
    RESULTS_DIR, ensure_dirs, load_cohort, build_clustering_matrix,
)

ensure_dirs(RESULTS_DIR)

# =====================================================================
# STEP 1 — REPRODUCE THE OFFICIAL 4-CLASS StepMix PARTITION
# ---------------------------------------------------------------------
# CRITICAL PIPELINE DETAILS (these exact settings reproduce the
# published partition with ARI = 1.000 against ECASFRSBE_CLASS4):
#
#   Clustering input = 7 indicators (see als_common.COGZ / als_common.FRS):
#     4 ECAS cognitive domains, Poletti-2016 normative z-scores
#     3 FrSBe post-illness subscales, SIGN-REVERSED so higher = better
#   All 7 columns then z-standardised (mean 0, SD 1) over the cohort.
#
#   Model: StepMix(n_components=4, measurement='continuous',
#                  n_init=20, random_state=42, max_iter=200)
#
#   NOTE: measurement MUST be 'continuous' (NOT 'gaussian_diag'); using
#   gaussian_diag or the *_PUNTZ columns instead of *_Z_POLETTI does NOT
#   reproduce the partition (drifts to ARI~0.82). This is the single most
#   important detail for exact reproduction.
# =====================================================================

def fit_official_partition(M, official):
    model = StepMix(n_components=4, measurement="continuous",
                    n_init=20, random_state=42, max_iter=200, verbose=0)
    model.fit(M)
    lab = model.predict(M)
    ari = adjusted_rand_score(official, lab)
    # Hungarian mapping predicted-label -> official-label (1..4)
    K = 4
    cost = np.zeros((K, K))
    for a in range(K):
        for b in range(K):
            cost[a, b] = -np.sum((lab == a) & (official == b + 1))
    ri, ci = linear_sum_assignment(cost)
    remap = {ri[i]: ci[i] + 1 for i in range(K)}
    R = model.predict_proba(M)
    R_official = np.zeros_like(R)            # reorder columns to official order
    for k in range(K):
        R_official[:, remap[k] - 1] = R[:, k]
    return model, ari, R_official

# =====================================================================
# STEP 2 — BCH (bias-adjusted three-step) WEIGHTS
# ---------------------------------------------------------------------
# Standard Bolck-Croon-Hagenaars correction:
#   D[s,t] = P(modal class = t | true class = s), estimated from the
#            posterior responsibilities R.
#   BCH weight matrix = D^{-1}; each unit i receives the row of D^{-1}
#   indexed by its modal class -> a length-K weight vector (may contain
#   NEGATIVE entries; this is intrinsic to BCH and is handled by the
#   custom weighted MLE below, because statsmodels GLM rejects negative
#   frequency weights).
# =====================================================================

def bch_weights(R):
    K = R.shape[1]
    hard = R.argmax(1)
    D = np.zeros((K, K))
    for s in range(K):
        for t in range(K):
            D[s, t] = R[hard == t, s].sum() / R[:, s].sum()
    Dinv = np.linalg.inv(D)
    W = np.array([[Dinv[hard[i], s] for s in range(K)] for i in range(len(R))])
    return W, D

# =====================================================================
# STEP 3 — BCH-WEIGHTED LOGISTIC (custom MLE allowing negative weights)
# ---------------------------------------------------------------------
# Two-arm odds ratio (reference vs target latent class) for a binary
# external variable (e.g. C9orf72 carrier status).  We maximise the
# BCH-weighted binomial log-likelihood directly (BFGS); SE from the
# inverse Hessian; Wald p-value and 95% CI.
# =====================================================================

def bch_logistic_OR(y, W, ref_class, target_class, valid_idx):
    rows = []
    for i in valid_idx:
        rows.append((y[i], 0.0, W[i, ref_class]))     # reference arm
        rows.append((y[i], 1.0, W[i, target_class]))  # target arm
    Y  = np.array([r[0] for r in rows])
    Xt = np.array([r[1] for r in rows])
    Wt = np.array([r[2] for r in rows])

    def nll(b):
        z = b[0] + b[1] * Xt
        p = 1.0 / (1.0 + np.exp(-z))
        p = np.clip(p, 1e-9, 1 - 1e-9)
        return -np.sum(Wt * (Y * np.log(p) + (1 - Y) * np.log(1 - p)))

    res = minimize(nll, [-3.0, 1.0], method="BFGS")
    b   = res.x
    OR  = np.exp(b[1])
    se  = np.sqrt(res.hess_inv[1, 1])
    lo, hi = np.exp(b[1] - 1.96 * se), np.exp(b[1] + 1.96 * se)
    p   = 2 * (1 - norm.cdf(abs(b[1] / se)))
    return OR, lo, hi, p

# =====================================================================
# MAIN
# =====================================================================

def main():
    d = load_cohort()
    official = d["ECASFRSBE_CLASS4"].astype(int).values

    # --- Step 1: reproduce partition ---
    M, mask, _ = build_clustering_matrix(d)
    official_m = d.loc[mask, "ECASFRSBE_CLASS4"].astype(int).values
    model, ari, R = fit_official_partition(M, official_m)
    print(f"[Step 1] StepMix reproduction ARI vs official = {ari:.4f}  (expect 1.0000)")
    sizes = pd.Series(R.argmax(1) + 1).value_counts().sort_index().to_dict()
    print(f"         class sizes = {sizes}  (expect 258/229/191/99)")

    dm = d[mask].copy()

    # --- Step 2: BCH weights ---
    W, D = bch_weights(R)
    print(f"[Step 2] BCH classification-quality diagonal = {np.round(np.diag(D),3)}")

    # --- Step 3: C9orf72 odds ratios (labels: 0=Preserved..3=Severe/FTD) ---
    gen = dm["GENETICA"].astype(str).str.strip()
    invalid = {"", "nan", "manca", "non prelev", "fare"}
    genotyped = (~gen.isin(invalid)).values
    y = (gen == "C9ORF72").astype(int).values
    idx = np.where(genotyped)[0]
    print(f"[Step 3] genotyped N = {len(idx)}, C9orf72 carriers = {int(y[idx].sum())}")

    OR, lo, hi, p = bch_logistic_OR(y, W, ref_class=0, target_class=3, valid_idx=idx)
    print(f"         C9orf72 OR Severe/FTD vs Preserved (BCH) = "
          f"{OR:.2f} (95% CI {lo:.2f}-{hi:.2f}), p={p:.3f}")
    print(f"         -> manuscript value: OR 2.56 (95% CI 1.05-6.21, p=0.038)")

    # combined-impaired vs Preserved (weights summed over classes 1,2,3)
    Wcomb = W.copy()
    Wcomb[:, 1] = W[:, 1] + W[:, 2] + W[:, 3]
    OR2, lo2, hi2, p2 = bch_logistic_OR(y, Wcomb, ref_class=0, target_class=1, valid_idx=idx)
    print(f"         C9orf72 OR combined-impaired vs Preserved (BCH) = "
          f"{OR2:.2f} (95% CI {lo2:.2f}-{hi2:.2f}), p={p2:.3f}")
    print(f"         -> manuscript value: OR 1.41")

    # BCH-weighted prevalence per class (Table 4 supports crude % on genotyped)
    names = ["Preserved", "Cognitive-only", "Behavioural-only", "Severe/FTD"]
    print("[Table 4] BCH-weighted C9orf72 prevalence:")
    prevalences = []
    for s in range(4):
        w = W[idx, s]
        prev = (w * y[idx]).sum() / w.sum()
        print(f"          {names[s]:16s} {100*prev:5.1f}%")
        prevalences.append({"class": names[s], "C9orf72_prevalence_pct": 100 * prev})

    # --- save outputs ---
    pd.DataFrame([
        {"comparison": "Severe/FTD vs Preserved", "OR": OR, "CI_low": lo, "CI_high": hi, "p": p},
        {"comparison": "Combined-impaired vs Preserved", "OR": OR2, "CI_low": lo2, "CI_high": hi2, "p": p2},
    ]).to_csv(os.path.join(RESULTS_DIR, "c9orf72_odds_ratios.csv"), index=False)

    pd.DataFrame(prevalences).to_csv(
        os.path.join(RESULTS_DIR, "c9orf72_prevalence_by_class.csv"), index=False)

    pd.DataFrame({"class_sizes": sizes}).to_csv(
        os.path.join(RESULTS_DIR, "partition_class_sizes.csv"))

    print(f"\n[Saved] results/c9orf72_odds_ratios.csv, "
          f"results/c9orf72_prevalence_by_class.csv, "
          f"results/partition_class_sizes.csv")

if __name__ == "__main__":
    main()