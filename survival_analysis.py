#!/usr/bin/env python3
# =====================================================================
# ECAS + FrSBe cognitive-behavioural phenotyping in ALS
# SURVIVAL ANALYSIS BY PHENOTYPE (reviewer-requested, not in main manuscript)
#
# Uses the SURVIVAL / STATUS / TRACHEO_NEW columns explicitly kept in the
# reduced dataset "for reviewer-requested survival analyses only" (see
# README / DATA_DICTIONARY.csv). These are NOT reported in the manuscript
# itself -- this script exists in case a reviewer asks whether the four
# cognitive-behavioural phenotypes differ in survival / tracheostomy-free
# survival, which is a clinically obvious question for an ALS cohort that
# the current manuscript does not address.
#
# Column definitions (from DATA_DICTIONARY.csv):
#   SURVIVAL     : survival time in years (N=554 non-missing / 777)
#   STATUS       : vital status, 1 = deceased, 2 = alive
#   TRACHEO_NEW  : tracheostomy indicator (0/1), for a composite endpoint
#
# Produces:
#   results/figures/km_survival_by_phenotype.png
#   results/figures/km_tracheofree_by_phenotype.png
#   results/survival_logrank_summary.csv
#   results/survival_median_by_class.csv
#
# Phenotype names/order/colors now come from als_common.py (shared
# with every other script that displays the 4 phenotypes).
#
# ENVIRONMENT: lifelines, pandas, matplotlib
#   pip install lifelines --break-system-packages
# =====================================================================

import os
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from als_common import (
    RESULTS_DIR, FIG_DIR, ensure_dirs, load_cohort,
    PHENOTYPE_ORDER, PHENOTYPE_COLORS,
)

ensure_dirs(RESULTS_DIR, FIG_DIR)


def prep_survival_columns(d):
    d = d.copy()
    d["SURVIVAL"] = pd.to_numeric(d["SURVIVAL"], errors="coerce")
    d["STATUS"] = pd.to_numeric(d["STATUS"], errors="coerce")
    d["TRACHEO_NEW"] = pd.to_numeric(d["TRACHEO_NEW"], errors="coerce")

    # event = 1 (deceased) -> 1 ; event = 2 (alive) -> 0 (censored)
    d["event_death"] = (d["STATUS"] == 1).astype(float)
    d.loc[d["STATUS"].isna(), "event_death"] = pd.NA

    # composite endpoint: death OR tracheostomy (both recorded, same duration
    # column -- this is a simplification given only one time-to-event column
    # is available; flagged as a limitation in the printed output below)
    d["event_composite"] = ((d["STATUS"] == 1) | (d["TRACHEO_NEW"] == 1)).astype(float)
    d.loc[d["STATUS"].isna() & d["TRACHEO_NEW"].isna(), "event_composite"] = pd.NA

    return d


def km_plot_and_test(d, event_col, duration_col, title, out_png):
    valid = d.dropna(subset=[duration_col, event_col, "ph"])
    n_total = len(d)
    n_valid = len(valid)
    print(f"\n--- {title} ---")
    print(f"N with non-missing {duration_col}/{event_col}: {n_valid} / {n_total}")

    fig, ax = plt.subplots(figsize=(7, 5))
    kmf = KaplanMeierFitter()
    medians = []
    for cls in PHENOTYPE_ORDER:
        sub = valid[valid["ph"] == cls]
        if len(sub) == 0:
            continue
        kmf.fit(sub[duration_col], event_observed=sub[event_col], label=f"{cls} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=PHENOTYPE_COLORS[cls])
        med = kmf.median_survival_time_
        medians.append({"phenotype": cls, "n": len(sub),
                         "n_events": int(sub[event_col].sum()),
                         "median_survival_years": med})
        print(f"  {cls:18s} n={len(sub):4d}  events={int(sub[event_col].sum()):4d}  "
              f"median={med if med == med else 'not reached'}")

    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[saved] {out_png}")

    # log-rank test across all 4 groups simultaneously
    result = multivariate_logrank_test(
        valid[duration_col], valid["ph"], valid[event_col]
    )
    print(f"  Multivariate log-rank test across 4 phenotypes: "
          f"chi2={result.test_statistic:.2f}, p={result.p_value:.4g}")

    return medians, result.test_statistic, result.p_value


def main():
    d = load_cohort(add_phenotype_label=True)
    d = prep_survival_columns(d)

    all_summary = []

    # --- overall survival ---
    medians1, chi2_1, p1 = km_plot_and_test(
        d, "event_death", "SURVIVAL",
        "Kaplan-Meier overall survival by cognitive-behavioural phenotype",
        os.path.join(FIG_DIR, "km_survival_by_phenotype.png"),
    )
    for m in medians1:
        m["endpoint"] = "overall_survival"
    all_summary += medians1

    # --- tracheostomy-free survival (composite endpoint) ---
    print("\nNOTE: the composite endpoint (death OR tracheostomy) uses the same "
          "SURVIVAL time column for both event types because only one "
          "time-to-event column is available in this reduced dataset. If "
          "tracheostomy and death do not occur at the same follow-up time in "
          "the source data, this is an approximation -- ideally a dedicated "
          "time-to-tracheostomy-or-death variable should be used instead.")
    medians2, chi2_2, p2 = km_plot_and_test(
        d, "event_composite", "SURVIVAL",
        "Kaplan-Meier tracheostomy-free survival by phenotype",
        os.path.join(FIG_DIR, "km_tracheofree_by_phenotype.png"),
    )
    for m in medians2:
        m["endpoint"] = "tracheostomy_free_survival"
    all_summary += medians2

    # --- save CSV summaries ---
    pd.DataFrame(all_summary).to_csv(
        os.path.join(RESULTS_DIR, "survival_median_by_class.csv"), index=False)

    pd.DataFrame([
        {"endpoint": "overall_survival", "logrank_chi2": chi2_1, "logrank_p": p1},
        {"endpoint": "tracheostomy_free_survival", "logrank_chi2": chi2_2, "logrank_p": p2},
    ]).to_csv(os.path.join(RESULTS_DIR, "survival_logrank_summary.csv"), index=False)

    print(f"\n[saved] results/survival_median_by_class.csv")
    print(f"[saved] results/survival_logrank_summary.csv")


if __name__ == "__main__":
    main()