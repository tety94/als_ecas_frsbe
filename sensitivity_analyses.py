#!/usr/bin/env python3
# =====================================================================
# ECAS + FrSBe cognitive-behavioural phenotyping in ALS
# SENSITIVITY / ROBUSTNESS ANALYSES — companion to reproduce_analyses.py
#
# Adds three checks requested in peer review:
#   (A) diagonal-covariance StepMix vs full-covariance GMM
#       -> BIC, log-likelihood, ARI vs official partition
#   (B) StepMix vs plain k-means
#       -> ARI vs official partition
#   (C) missingness check on the 7 clustering indicators
#       -> % missing per variable, and whether patients with ANY missing
#          clustering indicator differ in severity from complete cases
#          (using variables available in the reduced 31-column extract)
#
# Data loading and the 7-indicator clustering matrix now come from
# als_common.py (shared with reproduce_analyses.py), which also means
# this script picks up the frozen .sav automatically when present
# instead of being CSV-only, as it silently was before.
#
# Run locally, on the real data file, which is NOT uploaded here.
#
# ENVIRONMENT: Python 3.12; stepmix==3.0.0, scikit-learn, scipy, numpy, pandas
# =====================================================================

import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from stepmix.stepmix import StepMix
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from scipy import stats

from als_common import (
    COGZ, FRS, RESULTS_DIR, ensure_dirs, load_cohort,
    build_clustering_matrix, remap_labels_to_official,
)

ensure_dirs(RESULTS_DIR)


# =====================================================================
# (A) Diagonal-covariance StepMix vs full-covariance GMM
# =====================================================================
def check_covariance_structure(M, official_m):
    print("=" * 70)
    print("(A) DIAGONAL vs FULL-COVARIANCE — BIC / log-likelihood / ARI")
    print("=" * 70)

    # -- official model: StepMix, continuous (diagonal) --
    model_diag = StepMix(n_components=4, measurement="continuous",
                          n_init=20, random_state=42, max_iter=200, verbose=0)
    model_diag.fit(M)
    lab_diag = model_diag.predict(M)
    bic_diag = model_diag.bic(M)
    ll_diag  = model_diag.score(M) * len(M)
    ari_diag = adjusted_rand_score(official_m, lab_diag)

    # -- full-covariance GMM, same K, same n_init/seed --
    gmm_full = GaussianMixture(n_components=4, covariance_type="full",
                                n_init=20, random_state=42, max_iter=200)
    gmm_full.fit(M)
    lab_full = gmm_full.predict(M)
    bic_full = gmm_full.bic(M)
    ll_full  = gmm_full.score(M) * len(M)
    ari_full = adjusted_rand_score(official_m, lab_full)

    # -- ARI between the two solutions themselves (do they agree at all?) --
    ari_between = adjusted_rand_score(lab_diag, lab_full)

    print(f"{'Model':22s} {'BIC':>12s} {'LogLik':>12s} {'ARI vs official':>18s}")
    print(f"{'Diagonal (official)':22s} {bic_diag:12.1f} {ll_diag:12.1f} {ari_diag:18.4f}")
    print(f"{'Full covariance':22s} {bic_full:12.1f} {ll_full:12.1f} {ari_full:18.4f}")
    print(f"\nARI between diagonal and full-cov solutions directly: {ari_between:.4f}")
    print("(low value here would indicate the 4-class structure is sensitive")
    print(" to the local-independence assumption, not just to labeling noise)\n")

    if bic_full < bic_diag:
        print(f"NOTE: full-covariance BIC is LOWER by {bic_diag - bic_full:.1f} "
              f"-> full-cov model preferred by BIC over the diagonal model used "
              f"in the manuscript. Report this explicitly if it occurs.")
    else:
        print(f"Diagonal BIC is lower/comparable (full-cov worse by "
              f"{bic_full - bic_diag:.1f}) -> supports the diagonal choice, "
              f"but this number should be quoted in the manuscript, not just implied.")
    print()

    # --- save summary table ---
    summary = pd.DataFrame([
        {"model": "Diagonal (official)", "BIC": bic_diag, "LogLik": ll_diag, "ARI_vs_official": ari_diag},
        {"model": "Full covariance",     "BIC": bic_full, "LogLik": ll_full, "ARI_vs_official": ari_full},
    ])
    summary.to_csv(os.path.join(RESULTS_DIR, "covariance_comparison_summary.csv"), index=False)

    # --- save per-patient labels for downstream plotting (contingency heatmap etc.) ---
    labels_df = pd.DataFrame({
        "official_class": official_m,
        "diagonal_class": remap_labels_to_official(lab_diag, official_m),
        "fullcov_class":  remap_labels_to_official(lab_full, official_m),
    })
    labels_df.to_csv(os.path.join(RESULTS_DIR, "covariance_comparison_labels.csv"), index=False)

    return lab_diag, lab_full


# =====================================================================
# (B) StepMix vs k-means
# =====================================================================
def check_vs_kmeans(M, official_m):
    print("=" * 70)
    print("(B) STEPMIX vs K-MEANS")
    print("=" * 70)
    km = KMeans(n_clusters=4, n_init=20, random_state=42).fit(M)
    ari_km = adjusted_rand_score(official_m, km.labels_)
    print(f"K-means ARI vs official partition = {ari_km:.4f}")
    print("(manuscript states StepMix was preferred over k-means with")
    print(" 'comparable reproducibility' -- this is the number that claim rests on)\n")

    pd.DataFrame([{"model": "K-means", "ARI_vs_official": ari_km}]).to_csv(
        os.path.join(RESULTS_DIR, "kmeans_comparison_summary.csv"), index=False)

    return km.labels_


# =====================================================================
# (C) Missingness on the 7 clustering indicators
# =====================================================================
def check_missingness(d):
    print("=" * 70)
    print("(C) MISSINGNESS ON CLUSTERING INDICATORS")
    print("=" * 70)

    raw = d[COGZ + FRS].apply(pd.to_numeric, errors="coerce")
    miss_pct = raw.isna().mean() * 100
    print("Missing % per clustering variable (within the N=777 reduced extract):")
    for c in COGZ + FRS:
        print(f"  {c:28s} {miss_pct[c]:5.2f}%")

    miss_pct.rename("pct_missing").to_frame().to_csv(
        os.path.join(RESULTS_DIR, "missingness_clustering_vars.csv"))

    any_missing = raw.isna().any(axis=1)
    n_incomplete = int(any_missing.sum())
    print(f"\nRows with >=1 missing clustering indicator: {n_incomplete} / {len(d)}")

    if n_incomplete == 0:
        print("No missingness on the clustering indicators within this reduced")
        print("extract (expected, since README states this file already keeps")
        print("only the N=777 with a non-missing official partition).")
        print("\n*** IMPORTANT ***")
        print("This check CANNOT assess complete-case bias relative to the full")
        print("eligible PARALS population, because that comparison requires the")
        print("original ~555-column database (eligible-but-excluded patients),")
        print("which is not part of this reduced extract. Table S3 in the")
        print("manuscript (included vs eligible) is the only available proxy,")
        print("and it compares aggregate characteristics, not the missingness")
        print("mechanism (MAR vs MNAR) directly. Flag this as an unresolved")
        print("limitation rather than a checked-and-cleared item.")
        return

    # If there IS some missingness (e.g. re-running on the fuller dataset),
    # compare severity between complete vs incomplete cases on whichever
    # covariates are available (age, education, site of onset).
    d["_any_missing_cluster_var"] = any_missing
    covars = ["ETATEST", "Scol"]
    for c in covars:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
            g0 = d.loc[~any_missing, c].dropna()
            g1 = d.loc[any_missing, c].dropna()
            if len(g0) > 1 and len(g1) > 1:
                stat, p = stats.mannwhitneyu(g0, g1, alternative="two-sided")
                print(f"  {c}: complete-case mean={g0.mean():.2f}, "
                      f"incomplete-case mean={g1.mean():.2f}, "
                      f"Mann-Whitney p={p:.3f}")
    print()


# =====================================================================
# MAIN
# =====================================================================
def main():
    d = load_cohort()
    M, mask, _ = build_clustering_matrix(d)
    official_m = d.loc[mask, "ECASFRSBE_CLASS4"].astype(int).values
    print(f"Cohort N = {len(d)}, complete on 7 clustering indicators = {mask.sum()}\n")

    check_missingness(d)
    check_covariance_structure(M, official_m)
    check_vs_kmeans(M, official_m)


if __name__ == "__main__":
    main()