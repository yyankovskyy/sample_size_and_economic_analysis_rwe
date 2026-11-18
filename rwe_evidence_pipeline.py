#!/usr/bin/env python3
"""Sample Size and Economic Analysis for RWE Observational Studies.

Includes exploratory propensity-score weighted outcome analysis and a separate
pilot-based power/sample-size planning command. The original `analyze` command
by itself does NOT estimate required sample size.

Runs each authorized clinical data source independently. Never joins Cosmos and
NARCOMS patient records. It expects a documented, one-row-per-person target-trial
cohort constructed under an approved study protocol. It DOES NOT connect to Epic
or NARCOMS, assign diagnosis codes, link data, or validate clinical endpoints.

Dependencies: pandas, numpy, scipy, scikit-learn, matplotlib.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

NUMERIC_BASELINE = [
    "age", "baseline_severity", "comorbidity_index", "prior_hospitalizations_365d"
]
CATEGORICAL_BASELINE = ["sex"]
FIXED_OUTCOMES = {
    "clinical_response_365d": ("binary", "Clinical response by day 365"),
    "hospitalizations_365d": ("count", "Hospital admissions in 365 days"),
    "ed_visits_365d": ("count", "Emergency department visits in 365 days"),
    "hospital_days_365d": ("continuous", "Hospital days in 365 days"),
}
REQUIRED = [
    "patient_id", "index_date", "treated", "followup_days",
    *NUMERIC_BASELINE, *CATEGORICAL_BASELINE,
]


def weighted_mean(x, weight):
    x, weight = np.asarray(x, dtype=float), np.asarray(weight, dtype=float)
    good = np.isfinite(x) & np.isfinite(weight) & (weight >= 0)
    if good.sum() == 0 or weight[good].sum() <= 0:
        return float("nan")
    return float(np.sum(x[good] * weight[good]) / np.sum(weight[good]))


def effective_n(weights):
    w = np.asarray(weights, dtype=float)
    s = float(w.sum())
    return s * s / float((w * w).sum()) if (w * w).sum() > 0 else 0.0


def make_ps(df, estimand, seed=2026):
    numeric = [c for c in NUMERIC_BASELINE if c in df]
    categorical = [c for c in CATEGORICAL_BASELINE if c in df]
    preprocess = ColumnTransformer([
        ("num", Pipeline([
            ("impute", SimpleImputer(strategy="median", add_indicator=True)),
            ("scale", StandardScaler()),
        ]), numeric),
        ("cat", Pipeline([
            ("impute", SimpleImputer(strategy="most_frequent")),
            ("encode", OneHotEncoder(handle_unknown="ignore")),
        ]), categorical),
    ])
    model = Pipeline([
        ("prep", preprocess),
        ("logit", LogisticRegression(max_iter=1500, C=1.0, random_state=seed)),
    ])
    model.fit(df[numeric + categorical], df["treated"].astype(int))
    raw = model.predict_proba(df[numeric + categorical])[:, 1]
    ps = np.clip(raw, 1e-5, 1 - 1e-5)
    z = df["treated"].to_numpy(dtype=int)
    if estimand == "ATO":
        weight = np.where(z == 1, 1 - ps, ps)
    elif estimand == "ATT":
        weight = np.where(z == 1, 1.0, ps / (1 - ps))
    elif estimand == "ATE":
        weight = np.where(z == 1, 1 / ps, 1 / (1 - ps))
    else:
        raise ValueError("estimand must be ATO, ATT or ATE")
    return ps, weight, int(((raw <= 1e-5) | (raw >= 1 - 1e-5)).sum())


def bhattacharyya(ps, treatment):
    r = float(np.mean(treatment))
    return float(np.mean(np.sqrt(ps * (1 - ps))) / np.sqrt(r * (1 - r)))


def w_smd(x, treated, weight):
    z = np.asarray(treated, dtype=bool)
    x = np.asarray(x, dtype=float)
    weight = np.asarray(weight, dtype=float)
    m = []
    v = []
    for group in (False, True):
        a = x[z == group]
        w = weight[z == group]
        good = np.isfinite(a) & np.isfinite(w)
        if not good.any() or w[good].sum() <= 0:
            return float("nan")
        a, w = a[good], w[good]
        mean = weighted_mean(a, w)
        m.append(mean)
        v.append(weighted_mean((a - mean) ** 2, w))
    pooled = math.sqrt(max((v[0] + v[1]) / 2, 0))
    return float((m[1] - m[0]) / pooled) if pooled > 1e-12 else 0.0


def balance_table(df, weight):
    tr = df["treated"].to_numpy(dtype=int)
    results = []
    for c in NUMERIC_BASELINE:
        vals = pd.to_numeric(df[c], errors="coerce").to_numpy()
        results.append({"covariate": c,
                        "smd_before": w_smd(vals, tr, np.ones(len(df))),
                        "smd_after": w_smd(vals, tr, weight)})
    for c in CATEGORICAL_BASELINE:
        levels = sorted(df[c].dropna().astype(str).unique())
        for level in levels:
            vals = (df[c].astype("string") == level).fillna(False).astype(float).to_numpy()
            results.append({"covariate": f"{c}={level}",
                            "smd_before": w_smd(vals, tr, np.ones(len(df))),
                            "smd_after": w_smd(vals, tr, weight)})
    return pd.DataFrame(results)


def km_risk_at(df, weight, arm, horizon):
    """Weighted KM cumulative incidence; assumes non-informative censoring.

    Competing events are NOT addressed. Events at a censoring time are processed
    before censoring. A cause-specific event variable alone is not enough when
    death is a meaningful competing risk: use Aalen-Johansen in that setting.
    """
    group = df["treated"].to_numpy(dtype=int) == arm
    time = df.loc[group, "safety_time_days"].to_numpy(dtype=float)
    event = df.loc[group, "safety_event"].to_numpy(dtype=int)
    w = np.asarray(weight)[group]
    good = np.isfinite(time) & np.isfinite(w) & (time > 0)
    time, event, w = time[good], event[good], w[good]
    if len(time) == 0:
        return np.nan
    surv = 1.0
    for t in np.unique(time[(event == 1) & (time <= horizon)]):
        at_risk = w[time >= t].sum()
        occurred = w[(time == t) & (event == 1)].sum()
        if at_risk <= 0:
            break
        surv *= max(0.0, 1.0 - occurred / at_risk)
    return float(np.clip(1.0 - surv, 0.0, 1.0))


def subset_complete(df, field, horizon):
    return df["followup_days"].ge(horizon) & df[field].notna()


def build_metrics(df, weight, horizon, cost_basis, thresholds, responder_wtp=None):
    """Return aggregate rows and dictionary of scalar effects for bootstrap."""
    scalar = {}
    rows = []
    benefit = None
    z = df["treated"].to_numpy(dtype=int)

    if {"safety_event", "safety_time_days"}.issubset(df.columns):
        valid = df["safety_event"].notna() & df["safety_time_days"].notna()
        sub, ww = df.loc[valid], np.asarray(weight)[valid.to_numpy()]
        if (sub["treated"] == 0).sum() >= 10 and (sub["treated"] == 1).sum() >= 10:
            risks = [km_risk_at(sub, ww, arm, horizon) for arm in (0, 1)]
            rd = risks[1] - risks[0]
            scalar["safety_risk_difference"] = rd
            rows.append({"domain": "safety", "outcome": f"Adverse-event risk by day {horizon}",
                         "measure": "weighted_KM_risk_difference", "control": risks[0],
                         "treated": risks[1], "effect": rd, "n_control": int((sub["treated"] == 0).sum()),
                         "n_treated": int((sub["treated"] == 1).sum()),
                         "metric_key": "safety_risk_difference"})
            if risks[0] > 1e-9:
                scalar["safety_risk_ratio"] = risks[1] / risks[0]
                rows.append({"domain": "safety", "outcome": f"Adverse-event risk by day {horizon}",
                             "measure": "weighted_KM_risk_ratio", "control": risks[0],
                             "treated": risks[1], "effect": risks[1] / risks[0],
                             "n_control": int((sub["treated"] == 0).sum()),
                             "n_treated": int((sub["treated"] == 1).sum()),
                             "metric_key": "safety_risk_ratio"})

    for col, (kind, label) in FIXED_OUTCOMES.items():
        if col not in df:
            continue
        valid = subset_complete(df, col, horizon).to_numpy()
        sub, ww = df.loc[valid], np.asarray(weight)[valid]
        if (sub["treated"] == 0).sum() < 10 or (sub["treated"] == 1).sum() < 10:
            continue
        zz = sub["treated"].to_numpy(dtype=int)
        vals = sub[col].to_numpy(dtype=float)
        means = [weighted_mean(vals[zz == arm], ww[zz == arm]) for arm in (0, 1)]
        difference = means[1] - means[0]
        key = col + "_difference"
        scalar[key] = difference
        domain = "effectiveness" if kind == "binary" else "HCRU"
        rows.append({"domain": domain, "outcome": label,
                     "measure": "risk_difference" if kind == "binary" else "weighted_mean_difference",
                     "control": means[0], "treated": means[1], "effect": difference,
                     "n_control": int((zz == 0).sum()), "n_treated": int((zz == 1).sum()),
                     "metric_key": key})
        if kind == "binary" and means[0] > 1e-9:
            key = col + "_risk_ratio"
            scalar[key] = means[1] / means[0]
            rows.append({"domain": domain, "outcome": label, "measure": "risk_ratio",
                         "control": means[0], "treated": means[1], "effect": means[1] / means[0],
                         "n_control": int((zz == 0).sum()), "n_treated": int((zz == 1).sum()),
                         "metric_key": key})

    if cost_basis == "validated_costs" and "cost_365d_usd" in df:
        col = "cost_365d_usd"
        valid = subset_complete(df, col, horizon).to_numpy()
        sub, ww = df.loc[valid], np.asarray(weight)[valid]
        if (sub["treated"] == 0).sum() >= 10 and (sub["treated"] == 1).sum() >= 10:
            zz = sub["treated"].to_numpy(dtype=int)
            vals = sub[col].to_numpy(dtype=float)
            means = [weighted_mean(vals[zz == a], ww[zz == a]) for a in (0, 1)]
            scalar["cost_365d_difference"] = means[1] - means[0]
            rows.append({"domain": "economics", "outcome": "Healthcare costs over 365 days (USD)",
                         "measure": "weighted_mean_difference", "control": means[0],
                         "treated": means[1], "effect": means[1] - means[0],
                         "n_control": int((zz == 0).sum()), "n_treated": int((zz == 1).sum()),
                         "metric_key": "cost_365d_difference"})

        # Joint complete-case economic cohort: never combine mismatched samples.
        benefit = next((name for name in ("qaly_365d", "clinical_response_365d")
                        if name in df and (subset_complete(df, name, horizon) &
                                           df[col].notna()).sum() >= 20), None)
        if benefit:
            mask = subset_complete(df, benefit, horizon) & df[col].notna()
            sub, ww = df.loc[mask], np.asarray(weight)[mask.to_numpy()]
            if (sub["treated"] == 0).sum() >= 10 and (sub["treated"] == 1).sum() >= 10:
                zz = sub["treated"].to_numpy(dtype=int)
                benefit_vals = sub[benefit].to_numpy(dtype=float)
                cost_vals = sub[col].to_numpy(dtype=float)
                b = [weighted_mean(benefit_vals[zz == a], ww[zz == a]) for a in (0, 1)]
                c = [weighted_mean(cost_vals[zz == a], ww[zz == a]) for a in (0, 1)]
                de, dc = b[1] - b[0], c[1] - c[0]
                scalar["economic_joint_delta_benefit"] = de
                scalar["economic_joint_delta_cost"] = dc
                for key, label, val, group in [
                    ("economic_joint_delta_benefit", f"Joint-cohort benefit ({benefit})", de, b),
                    ("economic_joint_delta_cost", "Joint-cohort costs (USD)", dc, c),
                ]:
                    rows.append({"domain": "economics", "outcome": label,
                                 "measure": "weighted_mean_difference", "control": group[0],
                                 "treated": group[1], "effect": val,
                                 "n_control": int((zz == 0).sum()),
                                 "n_treated": int((zz == 1).sum()), "metric_key": key})
                if abs(de) > 1e-12:
                    scalar["economic_icer"] = dc / de
                    rows.append({"domain": "economics", "outcome":
                                 ("USD per QALY gained" if benefit == "qaly_365d"
                                  else "USD per additional responder"),
                                 "measure": "ICER_unstable_ratio_no_CI", "control": np.nan,
                                 "treated": np.nan, "effect": dc / de,
                                 "n_control": int((zz == 0).sum()),
                                 "n_treated": int((zz == 1).sum()),
                                 "metric_key": "economic_icer"})
                # $/QALY and $/responder are not interchangeable.
                # $/responder thresholds must be explicitly chosen by the study.
                chosen_thresholds = (thresholds if benefit == "qaly_365d"
                                     else (responder_wtp or []))
                for threshold in chosen_thresholds:
                    key = f"inmb_{int(threshold)}"
                    scalar[key] = threshold * de - dc
                    nmb_treated = threshold * b[1] - c[1]
                    nmb_control = threshold * b[0] - c[0]
                    rows.append({"domain": "economics", "outcome":
                                 f"Incremental net monetary benefit, WTP=${threshold:,.0f}",
                                 "measure": "INMB_USD", "control": nmb_control,
                                 "treated": nmb_treated, "effect": nmb_treated - nmb_control,
                                 "n_control": int((zz == 0).sum()),
                                 "n_treated": int((zz == 1).sum()),
                                 "metric_key": key})
    # Keep a distinct, semantically complete measure for EACH result row.
    # The `measure_code` is an analysis-specific unique identifier; the
    # `statistic` retains the estimator family for reproducibility.
    canonical_codes = {
        "safety_risk_difference": ("adverse_event_365d_weighted_km_risk_difference", "weighted_KM_risk_difference"),
        "safety_risk_ratio": ("adverse_event_365d_weighted_km_risk_ratio", "weighted_KM_risk_ratio"),
        "clinical_response_365d_difference": ("clinical_response_365d_weighted_risk_difference", "risk_difference"),
        "clinical_response_365d_risk_ratio": ("clinical_response_365d_weighted_risk_ratio", "risk_ratio"),
        "hospitalizations_365d_difference": ("hospital_admissions_365d_weighted_mean_difference", "weighted_mean_difference"),
        "ed_visits_365d_difference": ("emergency_visits_365d_weighted_mean_difference", "weighted_mean_difference"),
        "hospital_days_365d_difference": ("hospital_days_365d_weighted_mean_difference", "weighted_mean_difference"),
        "cost_365d_difference": ("health_cost_365d_usd_weighted_mean_difference", "weighted_mean_difference"),
        "economic_joint_delta_benefit": ("joint_cohort_qaly_weighted_mean_difference" if benefit == "qaly_365d" else
                                         "joint_cohort_response_weighted_risk_difference", "weighted_mean_difference"),
        "economic_joint_delta_cost": ("joint_cohort_cost_usd_weighted_mean_difference", "weighted_mean_difference"),
        "economic_icer": ("joint_cohort_incremental_cost_effectiveness_ratio", "ICER_unstable_ratio_no_CI"),
    }
    for row in rows:
        key = row["metric_key"]
        code, statistic = canonical_codes.get(key, (
            f"joint_cohort_incremental_net_monetary_benefit_{key}", "INMB_USD"))
        row["measure"] = code
        row["statistic"] = statistic
        row["metric_code"] = key
    return rows, scalar


def validate(df, horizon, cost_basis):
    issues = []
    absent = [c for c in REQUIRED if c not in df]
    if absent:
        raise ValueError("Missing required normalized columns: " + ", ".join(absent))
    if df["patient_id"].isna().any() or df["patient_id"].duplicated().any():
        raise ValueError("patient_id must be nonmissing and unique within EACH source")
    z = pd.to_numeric(df["treated"], errors="coerce")
    if not z.isin([0, 1]).all() or z.nunique() < 2:
        raise ValueError("treated must contain 0 and 1, with no missing values")
    df["treated"] = z.astype(int)
    idx = pd.to_datetime(df["index_date"], errors="coerce")
    if idx.isna().any():
        raise ValueError("index_date has missing or invalid dates")
    for c in [*NUMERIC_BASELINE, "followup_days", *FIXED_OUTCOMES,
              "safety_event", "safety_time_days", "cost_365d_usd", "qaly_365d"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if df["followup_days"].isna().any() or (df["followup_days"] <= 0).any():
        raise ValueError("followup_days must be positive and nonmissing")
    if df["age"].isna().all():
        raise ValueError("age is entirely missing")
    if df["baseline_severity"].isna().all():
        raise ValueError("baseline_severity is entirely missing")
    if (df["age"] < 0).any() or (df["age"] > 120).any():
        issues.append("Age outside 0–120: inspect extraction and age masking")
    for c in ("clinical_response_365d", "safety_event"):
        if c in df and not df[c].dropna().isin([0, 1]).all():
            raise ValueError(f"{c} must contain 0/1 or NA")
    for c in ("hospitalizations_365d", "ed_visits_365d", "hospital_days_365d",
              "cost_365d_usd", "qaly_365d"):
        if c in df and (df[c] < 0).any():
            raise ValueError(f"{c} must not be negative")
    if "qaly_365d" in df and (df["qaly_365d"] > 1.0).any():
        issues.append("qaly_365d exceeds 1 in a 365-day horizon; check integration/time horizon")
    if "safety_event" in df or "safety_time_days" in df:
        if not {"safety_event", "safety_time_days"}.issubset(df.columns):
            raise ValueError("Safety analyses require both safety_event and safety_time_days")
        paired = df["safety_event"].notna() == df["safety_time_days"].notna()
        if not paired.all():
            raise ValueError("safety_event and safety_time_days must have matching missingness")
        valid = df["safety_time_days"].notna()
        if (df.loc[valid, "safety_time_days"] <= 0).any() or (
            df.loc[valid, "safety_time_days"] > df.loc[valid, "followup_days"] + 1e-7
        ).any():
            raise ValueError("Safety time must be >0 and <= available followup_days")
    fixed_year_cols = list(FIXED_OUTCOMES) + ["cost_365d_usd", "qaly_365d"]
    for c in fixed_year_cols:
        if c in df:
            invalid = df["followup_days"].lt(horizon) & df[c].notna()
            if invalid.any():
                # Explicitly prevent partial-follow-up outcomes from being analyzed
                # as complete one-year measures.
                df.loc[invalid, c] = np.nan
                issues.append(f"{int(invalid.sum())} {c} values masked: follow-up <{horizon} days")
    if cost_basis == "validated_costs" and "cost_365d_usd" not in df:
        issues.append("Validated cost requested but cost_365d_usd absent: economic analysis skipped")
    if cost_basis != "validated_costs" and "cost_365d_usd" in df:
        issues.append("Cost column present but cost basis not validated: economics suppressed")
    issues.append("All results are exploratory until cohort/exposure/outcome algorithms are validated")
    issues.append("Index dates are syntax-checked only; confirm every exposure, confounder and outcome time window upstream")
    issues.append("Patient bootstrap is not site-cluster robust; use cluster-aware bootstrap for multisite inference")
    issues.append("Complete-case 365-day outcomes assume follow-up availability is ignorable; consider IPCW")
    issues.append("KM safety risk assumes independent censoring and no competing-risk complication")
    return df, issues


def diagnostics(df, ps, weight, balance, estimand, horizon):
    z = df["treated"].to_numpy(dtype=int)
    r = float(z.mean())
    out = {
        "sample_n": int(len(df)), "treated_n": int(z.sum()),
        "control_n": int((z == 0).sum()), "treatment_fraction": r,
        "estimand": estimand, "horizon_days": horizon,
        "phi_bhattacharyya": bhattacharyya(ps, z),
        "ps_min": float(ps.min()), "ps_p01": float(np.quantile(ps, 0.01)),
        "ps_median": float(np.median(ps)), "ps_p99": float(np.quantile(ps, 0.99)),
        "ps_max": float(ps.max()),
        "ps_near_extreme_n": int(((ps < 0.01) | (ps > 0.99)).sum()),
        "weighted_ess_treated": effective_n(weight[z == 1]),
        "weighted_ess_control": effective_n(weight[z == 0]),
        "max_abs_smd_before": float(balance["smd_before"].abs().max()),
        "max_abs_smd_after": float(balance["smd_after"].abs().max()),
        "n_abs_smd_gt_0_1_after": int((balance["smd_after"].abs() > 0.10).sum()),
        "complete_365d_n": int(df["followup_days"].ge(horizon).sum()),
    }
    return out


def demo(n, seed, source):
    """Entirely SYNTHETIC demonstration; no correspondence with real patients."""
    rng = np.random.default_rng(seed)
    age = np.clip(rng.normal(52 if source == "NARCOMS" else 59, 14, n), 18, 90)
    severity = np.clip(rng.normal(0, 1, n), -2.5, 2.5)
    comorbidity = rng.poisson(np.exp(-0.1 + 0.4 * severity), n)
    pri_hosp = rng.poisson(np.exp(-1 + 0.25 * severity), n)
    sex = rng.choice(["F", "M"], n, p=[0.72, 0.28] if source == "NARCOMS" else [0.50, 0.50])
    p = expit(-0.15 + 0.035 * (age - 55) + 0.5 * severity + 0.12 * comorbidity)
    treated = rng.binomial(1, p)
    followup = rng.choice([280, 325, 365, 400, 480], n,
                          p=[0.07, 0.08, 0.48, 0.25, 0.12]).astype(float)
    hazard = 0.10 / 365 * np.exp(0.35 * severity + 0.015 * (age - 55) - 0.25 * treated)
    event_t = rng.exponential(1 / hazard)
    end_t = np.minimum(followup, 365)
    ae = (event_t <= end_t).astype(int)
    surv_t = np.minimum(event_t, end_t)
    response = rng.binomial(1, expit(-0.05 - 0.55 * severity + 0.45 * treated))
    hosp = rng.poisson(np.exp(-0.50 + 0.35 * severity + 0.15 * comorbidity - 0.13 * treated))
    ed = rng.poisson(np.exp(-0.20 + 0.28 * severity - 0.08 * treated))
    days = hosp * np.maximum(1, rng.poisson(2.5, n))
    qaly = np.clip(0.72 - 0.07 * severity + 0.04 * treated + rng.normal(0, 0.10, n), 0.1, 1)
    costs = np.maximum(0, 5400 + 4200 * hosp + 1200 * ed +
                       3100 * np.maximum(severity, 0) + rng.normal(0, 1500, n))
    df = pd.DataFrame({
        "patient_id": [f"SYNTHETIC_{source}_{i:06d}" for i in range(n)],
        "index_date": (pd.Timestamp("2024-01-01") +
                       pd.to_timedelta(rng.integers(0, 300, n), unit="D")).strftime("%Y-%m-%d"),
        "treated": treated, "age": np.round(age, 1), "sex": sex,
        "baseline_severity": np.round(severity, 3),
        "comorbidity_index": comorbidity,
        "prior_hospitalizations_365d": pri_hosp,
        "followup_days": followup,
        "safety_event": ae, "safety_time_days": np.round(surv_t, 3),
        "clinical_response_365d": response.astype(float),
        "hospitalizations_365d": hosp.astype(float),
        "ed_visits_365d": ed.astype(float),
        "hospital_days_365d": days.astype(float),
    })
    if source == "EpicCosmos":
        df["cost_365d_usd"] = costs.round(2)
        df["qaly_365d"] = np.round(qaly, 4)
    for c in list(FIXED_OUTCOMES) + ["cost_365d_usd", "qaly_365d"]:
        if c in df:
            df.loc[df["followup_days"] < 365, c] = np.nan
    return df


def save_charts(df, ps, balance, dest):
    z = df["treated"].to_numpy(dtype=int)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(ps[z == 0], bins=np.linspace(0, 1, 26), alpha=0.65, density=True, label="Control")
    ax.hist(ps[z == 1], bins=np.linspace(0, 1, 26), alpha=0.60, density=True, label="Treated")
    ax.set(xlabel="Estimated propensity score", ylabel="Density", title="Propensity score overlap")
    ax.legend()
    fig.tight_layout()
    fig.savefig(dest / "propensity_overlap.png", dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.50 * len(balance))))
    y = np.arange(len(balance))
    ax.scatter(balance["smd_before"].abs(), y, label="Unadjusted")
    ax.scatter(balance["smd_after"].abs(), y, label="Weighted")
    ax.axvline(0.10, color="gray", linestyle="--", label="|SMD| = 0.10")
    ax.set(yticks=y, yticklabels=balance["covariate"],
           xlabel="Absolute standardized mean difference", title="Covariate balance")
    ax.legend()
    fig.tight_layout()
    fig.savefig(dest / "covariate_balance.png", dpi=150)
    plt.close(fig)


def baseline_table_1(df, weight):
    """Descriptive Table 1 with raw and propensity-weighted baseline distributions."""
    rows = []
    arm = df["treated"].to_numpy(dtype=int)
    for col in NUMERIC_BASELINE:
        values = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        vals = [(col, values, "numeric")]
        for name, v, kind in vals:
            row = {"variable": name, "type": kind,
                   "missing_n": int(np.isnan(v).sum()),
                   "smd_unadjusted": w_smd(v, arm, np.ones(len(df))),
                   "smd_weighted": w_smd(v, arm, weight)}
            for a, label in [(0, "control"), (1, "treated")]:
                mask = arm == a
                row[f"{label}_mean_unadjusted"] = weighted_mean(v[mask], np.ones(mask.sum()))
                row[f"{label}_mean_weighted"] = weighted_mean(v[mask], weight[mask])
                row[f"{label}_n_observed"] = int(np.isfinite(v[mask]).sum())
            rows.append(row)
    for col in CATEGORICAL_BASELINE:
        for category in sorted(df[col].dropna().astype(str).unique()):
            v = (df[col].astype("string") == category).fillna(False).astype(float).to_numpy()
            row = {"variable": f"{col}={category}", "type": "proportion",
                   "missing_n": int(df[col].isna().sum()),
                   "smd_unadjusted": w_smd(v, arm, np.ones(len(df))),
                   "smd_weighted": w_smd(v, arm, weight)}
            for a, label in [(0, "control"), (1, "treated")]:
                mask = arm == a
                row[f"{label}_mean_unadjusted"] = weighted_mean(v[mask], np.ones(mask.sum()))
                row[f"{label}_mean_weighted"] = weighted_mean(v[mask], weight[mask])
                row[f"{label}_n_observed"] = int(df.loc[mask, col].notna().sum())
            rows.append(row)
    return pd.DataFrame(rows)


def aggregate_cohort(df, horizon, min_cell):
    rows = []
    for arm, label in [(0, "control"), (1, "treated")]:
        g = df.loc[df["treated"] == arm]
        n_events = int(g["safety_event"].fillna(0).sum()) if "safety_event" in df else None
        rows.append({
            "group": label, "patients": len(g),
            "full_365d_followup": int(g["followup_days"].ge(horizon).sum()),
            "baseline_mean_age": float(g["age"].mean()),
            "adverse_events": (f"<{min_cell}" if n_events is not None and 0 < n_events < min_cell
                               else n_events),
            "median_followup_days": float(g["followup_days"].median()),
        })
    return pd.DataFrame(rows)


def run_analysis(df, args):
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    # New cohort output invalidates old prospective planning scenarios.
    for stale in ("sample_size_plan.csv", "sample_size_plan.md"):
        (output / stale).unlink(missing_ok=True)
    df, notes = validate(df.copy(), args.horizon, args.cost_basis)
    if len(df) < 30:
        raise ValueError("At least 30 participants are required for this exploratory demonstration")
    ps, weight, n_clipped = make_ps(df, args.estimand, args.seed)
    balance = balance_table(df, weight)
    diag = diagnostics(df, ps, weight, balance, args.estimand, args.horizon)
    diag["source"] = args.source
    diag["dataset_filename"] = Path(args.input).name
    diag["data_origin"] = ("synthetic_demo" if
                           df["patient_id"].astype(str).str.startswith("SYNTHETIC_").all()
                           or Path(args.input).name.upper().startswith("SYNTHETIC_")
                           else "analyst_supplied_unverified")
    diag["bootstrap_replicates_requested"] = args.bootstrap
    diag["ps_predictions_clipped_n"] = n_clipped
    diag["missingness_fraction"] = {col: float(df[col].isna().mean())
                                   for col in list(FIXED_OUTCOMES) + [
                                       "cost_365d_usd", "qaly_365d"] if col in df}
    if diag["max_abs_smd_after"] > 0.10:
        notes.append("WARNING: Residual |SMD| >0.10; confounding balance may be inadequate")
    if diag["weighted_ess_treated"] < 0.2 * diag["treated_n"] or (
        diag["weighted_ess_control"] < 0.2 * diag["control_n"]
    ):
        notes.append("WARNING: Very low effective sample size relative to raw sample")
    if n_clipped:
        notes.append(f"WARNING: {n_clipped} propensity predictions were clipped")
    if args.bootstrap < 500:
        notes.append("For a formal report, use >=500-2000 bootstrap resamples and verify stability")

    rows, original = build_metrics(df, weight, args.horizon, args.cost_basis,
                                  args.thresholds, args.responder_wtp)
    boot = defaultdict(list)
    rng = np.random.default_rng(args.seed)
    i0 = np.flatnonzero(df["treated"].to_numpy(dtype=int) == 0)
    i1 = np.flatnonzero(df["treated"].to_numpy(dtype=int) == 1)
    failures = 0
    for _ in range(args.bootstrap):
        ids = np.r_[rng.choice(i0, len(i0), replace=True),
                    rng.choice(i1, len(i1), replace=True)]
        sample = df.iloc[ids].reset_index(drop=True)
        try:
            _, sw, _ = make_ps(sample, args.estimand, args.seed)
            _, estimates = build_metrics(sample, sw, args.horizon,
                                         args.cost_basis, args.thresholds, args.responder_wtp)
            for k, val in estimates.items():
                if np.isfinite(val):
                    boot[k].append(float(val))
        except (ValueError, FloatingPointError) as exc:
            failures += 1
            if failures == 1:
                notes.append(f"Bootstrap estimation failure (first): {str(exc)[:120]}")
    diag["bootstrap_fit_failures"] = failures
    for row in rows:
        key = row["metric_key"]
        samples = np.asarray(boot.get(key, []), dtype=float)
        row["bootstrap_se"] = (float(samples.std(ddof=1)) if len(samples) >= 2
                               and key != "economic_icer" else np.nan)
        if key == "economic_icer":
            row.update({"ci_low": np.nan, "ci_high": np.nan,
                        "n_bootstrap": len(samples), "bootstrap_probability_positive": np.nan})
        elif len(samples) >= max(20, 0.70 * args.bootstrap):
            row.update({"ci_low": float(np.quantile(samples, 0.025)),
                        "ci_high": float(np.quantile(samples, 0.975)),
                        "n_bootstrap": len(samples),
                        "bootstrap_probability_positive":
                        float(np.mean(samples > 0)) if key.startswith("inmb_") else np.nan})
        else:
            row.update({"ci_low": np.nan, "ci_high": np.nan,
                        "n_bootstrap": len(samples), "bootstrap_probability_positive": np.nan})
            notes.append(f"No stable bootstrap interval for {key} ({len(samples)} replicates)")
        row.pop("metric_key", None)
    results = pd.DataFrame(rows)
    if results.empty:
        raise ValueError("No analyzable outcomes: check 365-day follow-up and outcome fields")

    balance.to_csv(output / "covariate_balance.csv", index=False)
    baseline_table_1(df, weight).to_csv(output / "table1_baseline.csv", index=False)
    results.to_csv(output / "effect_estimates.csv", index=False)
    cohort = aggregate_cohort(df, args.horizon, args.min_cell)
    cohort.to_csv(output / "cohort_summary.csv", index=False)
    pd.DataFrame([{"variable": c, "missing_n": int(df[c].isna().sum()),
                   "missing_fraction": float(df[c].isna().mean())}
                  for c in df.columns if c not in ("patient_id", "index_date")]).to_csv(
        output / "missingness.csv", index=False)
    save_charts(df, ps, balance, output)
    diag["warnings_and_limitations"] = notes
    diag["dependencies"] = {"python": sys.version.split()[0],
                           "pandas": pd.__version__, "numpy": np.__version__}
    (output / "diagnostics.json").write_text(json.dumps(diag, indent=2), encoding="utf-8")
    filename = str(Path(args.input).name) if args.input else "generated synthetic demo"
    txt = [
        f"# Exploratory RWE evidence summary — {args.source}", "",
        "**Not a validated regulatory submission or a peer-reviewed analysis.**", "",
        f"Dataset: `{filename}` | Estimand: **{args.estimand}** | ",
        f"Horizon: **{args.horizon} days** | Bootstraps: **{args.bootstrap}**", "",
        "## Study population and diagnostics", "",
        f"- Participants: {diag['sample_n']:,}; treated: {diag['treated_n']:,}; "
        f"control: {diag['control_n']:,}.",
        f"- Bhattacharyya propensity score overlap: {diag['phi_bhattacharyya']:.3f}.",
        f"- Weighted ESS: treated {diag['weighted_ess_treated']:.3f}, "
        f"control {diag['weighted_ess_control']:.3f}.",
        f"- Maximum weighted absolute SMD: {diag['max_abs_smd_after']:.3f}.", "",
        "## Adjusted estimates (treated versus active comparator)", "",
        "| Domain | Outcome | Measure | Treated | Control | Effect | 95% bootstrap CI |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for r in results.to_dict("records"):
        fmt = lambda v: f"{v:,.3f}" if pd.notna(v) else "NA"
        ci = f"[{fmt(r['ci_low'])}, {fmt(r['ci_high'])}]" if pd.notna(r["ci_low"]) else "not estimated"
        txt.append(f"| {r['domain']} | {r['outcome']} | {r['measure']} | "
                   f"{fmt(r['treated'])} | {fmt(r['control'])} | "
                   f"{fmt(r['effect'])} | {ci} |")
    txt += ["", "## Interpretation and constraints", "",
            "- Positive differences mean higher values in the treated cohort; for adverse "
            "events, higher values generally indicate greater risk.",
            "- The safety estimate is weighted Kaplan-Meier cumulative event risk; "
            "independent censoring is assumed. For competing risks, use weighted Aalen-Johansen.",
            "- The causal interpretation requires consistency, exchangeability given measured "
            "pretreatment covariates, positivity and valid outcome ascertainment.",
            "- Complete-case follow-up is NOT automatically missing at random. Missing outcome "
            "and censoring mechanisms must be justified or handled with IPCW/other sensitivity analyses.",
            "- The outcome columns must be based on validated endpoint definitions and a common "
            "clinical index date in a new-user active-comparator design.",
            "- The ICER is only descriptive; use a cost-effectiveness plane/CEAC, prespecified "
            "willingness-to-pay thresholds, uncertainty and perspective in a formal economic study.",
            "- A real comparative-device study requires confirmed device exposure in BOTH groups. "
            "NARCOMS DMT comparisons are not device evidence without captured device exposure.",
            "- Outputs are aggregate, but all sharing must be reviewed under the source-specific "
            "DUA, IRB, and suppression/disclosure requirements.", "",
            "## Warnings", "",
        ]
    txt += [f"- {w}" for w in notes]
    txt += ["", "## Relevant reporting guidance", "",
            "- FDA, *Use of Real-World Evidence to Support Regulatory Decision-Making "
            "for Medical Devices*, December 2025.",
            "- TARGET statement (2025), target-trial emulation reporting.",
            "- RECORD statement (2015), routinely collected healthcare data reporting.",
            ""]
    (output / "evidence_summary.md").write_text("\n".join(txt), encoding="utf-8")
    if not args.no_word_report:
        word_path = generate_word_report(output, source=args.source)
        print(f"Word report: {word_path.resolve()}")
    print(f"Analyzed {args.source}: N={len(df):,}, results={len(results)} rows")
    print(f"Output directory: {output.resolve()}")
    display = results[["domain", "outcome", "measure", "effect", "ci_low", "ci_high"]].copy()
    print(display.to_string(index=False, formatters={
        key: (lambda x: f"{x:,.3f}" if pd.notna(x) else "NA")
        for key in ("effect", "ci_low", "ci_high")
    }))
    return results





def plan_sample_size_from_pilot(output_dir, metric_code, target_effect,
                                *, alpha=0.05, power=0.80,
                                additional_loss=0.0, design_multiplier=1.0,
                                sidedness="two-sided", refresh_word=True):
    """Exploratory pilot-bootstrap, normal-approximation sample-size scenario.

    This is NOT the Liu/Yang/Li PSpower method: its pilot bootstrap SE
    incorporates the fitted PS weights and observed follow-up pattern, but it
    cannot separately identify unmeasured confounding or extrapolate overlap
    to a different future population. All target effects must be chosen BEFORE
    study outcomes are inspected and justified clinically.

    With var(estimate) ~ V/N, V_hat = N_pilot * (bootstrap_SE)^2.
    Required N = ceil(V_hat * (z_crit + z_power)^2 / target_effect^2
                      * design_multiplier / (1 - additional_loss)).
    This N is enrolled participants, assuming the pilot's usual censoring and
    missingness persist. `additional_loss` means EXTRA losses beyond the pilot.
    """
    output = Path(output_dir)
    path = output / "effect_estimates.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}; first run analyze with this version.")
    if not 0 < alpha < 1 or not .5 < power < 1:
        raise ValueError("alpha must be between 0 and 1; power between 0.5 and 1")
    if not 0 <= additional_loss < 1:
        raise ValueError("additional_loss must be in [0, 1)")
    if not np.isfinite(target_effect) or abs(target_effect) <= 0:
        raise ValueError("target_effect must be a nonzero clinically meaningful difference")
    if not np.isfinite(design_multiplier) or design_multiplier < 1:
        raise ValueError("design_multiplier must be >=1")
    if sidedness not in ("one-sided", "two-sided"):
        raise ValueError("sidedness must be one-sided or two-sided")

    df = pd.read_csv(path)
    if "metric_code" not in df or "bootstrap_se" not in df:
        raise ValueError("This output predates pilot planning. Rerun the analyze command to compute bootstrap_se.")
    selected = df.loc[df["metric_code"] == metric_code]
    if len(selected) != 1:
        choices = df.loc[df["statistic"].isin(["risk_difference", "weighted_KM_risk_difference", "weighted_mean_difference"]),
                         "metric_code"].tolist()
        raise ValueError(f"Unknown/ambiguous metric_code '{metric_code}'. Candidate differences: {choices}")
    r = selected.iloc[0]
    allowed = ("risk_difference", "weighted_KM_risk_difference", "weighted_mean_difference")
    if r["statistic"] not in allowed or r["domain"] not in ("safety", "effectiveness", "HCRU"):
        raise ValueError("Only additive safety, effectiveness or HCRU differences are supported for pilot planning. "
                         "Risk ratios, cost/QALY ratios and INMB require their own scale-specific design.")
    if pd.isna(r["bootstrap_se"]) or float(r["bootstrap_se"]) <= 0 or r["n_bootstrap"] < 100:
        raise ValueError("At least 100 successful bootstrap draws and positive bootstrap_se required")
    diag_path = output / "diagnostics.json"
    if not diag_path.exists():
        raise FileNotFoundError(f"Missing {diag_path}")
    meta = json.loads(diag_path.read_text())
    pilot_n = int(meta["sample_n"])
    se = float(r["bootstrap_se"])
    variance_constant = pilot_n * se**2
    zcrit = float(norm.ppf(1 - alpha / (2 if sidedness == "two-sided" else 1)))
    zpower = float(norm.ppf(power))
    analysis_n = math.ceil(variance_constant * (zcrit + zpower)**2 /
                           target_effect**2 * design_multiplier)
    required_enrolled = math.ceil(analysis_n / (1-additional_loss))
    row = {
        "method": "exploratory_empirical_pilot_bootstrap_normal_approximation",
        "metric_code": str(metric_code), "outcome": str(r["outcome"]),
        "statistic": str(r["statistic"]), "source": str(meta.get("source", "Unknown")),
        "estimand": str(meta.get("estimand", "Unknown")),
        "pilot_n": pilot_n,
        "bootstrap_draws_used": int(r["n_bootstrap"]),
        "pilot_bootstrap_standard_error": se,
        "asymptotic_variance_coefficient": variance_constant,
        "bhattacharyya_phi_in_pilot": float(meta["phi_bhattacharyya"]),
        "target_effect_difference": float(abs(target_effect)),
        "alpha": alpha, "target_power": power, "test": sidedness,
        "design_multiplier": design_multiplier,
        "extra_loss_fraction_beyond_pilot": additional_loss,
        "required_n_before_additional_loss": analysis_n,
        "required_n_total_enrolled": required_enrolled,
    }
    pd.DataFrame([row]).to_csv(output / "sample_size_plan.csv", index=False)
    guide = [
        "# Prospective sample-size planning (pilot-based approximation)", "",
        "This is an exploratory empirical pilot-bootstrap normal approximation, NOT the",
        "Liu/Yang/Li PSpower analytical method. It does not estimate the confounding",
        "coefficient rho-squared or protect against unmeasured confounding.", "",
        f"- Outcome: {r['outcome']} ({metric_code})",
        f"- Estimand: {row['estimand']}; alpha = {alpha:.3f}; power = {power:.3f}.",
        f"- Clinically prespecified effect difference: {abs(target_effect):.3f} in endpoint units.",
        f"- Pilot N = {pilot_n}; bootstrap SE = {se:.3f}; ",
        f"  variance coefficient V = N * SE^2 = {variance_constant:.3f}.",
        f"- Pilot Bhattacharyya overlap = {row['bhattacharyya_phi_in_pilot']:.3f}.",
        f"- Required N before extra loss = {analysis_n:,}; total enrolled N = {required_enrolled:,}.",
        f"- Design multiplier = {design_multiplier:.3f}; extra loss beyond pilot = {additional_loss:.3f}.",
        "", "## Interpretation and limitations", "",
        "The pilot bootstrap variance reflects this cohort's fitted propensity-score weighting,",
        "measured confounders, overlap and existing missingness/follow-up pattern, but is",
        "not transportable without sensitivity analyses. Extra loss is beyond the pilot's loss",
        "already reflected in the SE. A small target effect is not automatically clinically",
        "meaningful. For a design-stage treatment-confounding sensitivity framework",
        "using overlap and a separate rho-squared parameter, consult the validated",
        "Duke PSpower calculator: https://www2.stat.duke.edu/~fl35/pspower.html.",
    ]
    (output / "sample_size_plan.md").write_text("\n".join(guide) + "\n", encoding="utf-8")
    if refresh_word:
        generate_word_report(output, source=str(meta.get("source", "Other")))
    print(f"Pilot planning scenario ({metric_code}): required enrolled N = {required_enrolled:,}")
    print(f"Saved: {(output/'sample_size_plan.csv').resolve()}")
    return row

# -----------------------------------------------------------------------------
# Publication-style Word reporting (exploratory; never implies source validation)
# -----------------------------------------------------------------------------

def _word_number(value, digits=3):
    try:
        number = float(value)
        return f"{number:,.{digits}f}" if np.isfinite(number) else "Not available"
    except (ValueError, TypeError):
        return str(value) if value is not None else "Not available"


def _word_int(value):
    try:
        return f"{int(value):,}" if pd.notna(value) else "Not available"
    except (TypeError, ValueError):
        return str(value)


def _word_extract_label(summary_text):
    import re
    match = re.search(r"Dataset:\s*`([^`]+)`", summary_text)
    return match.group(1) if match else "Unknown (see source audit trail)"


def _word_summary_items(summary_text, heading):
    """Preserve the relevant narrative caveats from the Markdown report."""
    import re
    match = re.search(r"^## " + re.escape(heading) + r"\s*$\n(.*?)(?=^## |\Z)",
                      summary_text, re.MULTILINE | re.DOTALL)
    if not match:
        return []
    return [m.group(1).strip() for m in re.finditer(r"^\s*-\s+(.+)$", match.group(1), re.MULTILINE)]


def _word_expand_terms(value):
    """Expand acronyms in externally generated narrative; preserve the original code."""
    import re
    terms = {
        "ICER": "incremental cost-effectiveness ratio (ICER)",
        "QALYs": "quality-adjusted life years (QALYs)",
        "QALY": "quality-adjusted life year (QALY)",
        "CEAC": "cost-effectiveness acceptability curve (CEAC)",
        "IPCW": "inverse probability of censoring weighting (IPCW)",
        "DUA": "data use agreement (DUA)",
        "IRB": "institutional review board (IRB)",
        "DMT": "disease-modifying therapy (DMT)",
        "HCRU": "healthcare resource utilization (HCRU)",
        "KM": "Kaplan-Meier (KM)",
        "SMD": "standardized mean difference (SMD)",
        "PS": "propensity score (PS)",
        "ESS": "effective sample size (ESS)",
        "FDA": "U.S. Food and Drug Administration (FDA)",
        "RWE": "real-world evidence (RWE)",
    }
    s = str(value)
    # Longest match first avoids partial replacement of QALYs with QALY.
    for key in sorted(terms, key=len, reverse=True):
        s = re.sub(r"(?<![\w(])" + re.escape(key) + r"\b", terms[key], s)
    return s


_WORD_INDICATORS = {
    "sample_n": "Total eligible patients",
    "treated_n": "Patients receiving the studied treatment",
    "control_n": "Patients receiving the active comparator",
    "treatment_fraction": "Proportion receiving the studied treatment",
    "estimand": "Target causal population",
    "horizon_days": "Follow-up horizon (days)",
    "phi_bhattacharyya": "Bhattacharyya propensity-score overlap coefficient",
    "ps_min": "Minimum estimated propensity score",
    "ps_p01": "First percentile of estimated propensity scores",
    "ps_median": "Median estimated propensity score",
    "ps_p99": "Ninety-ninth percentile of estimated propensity scores",
    "ps_max": "Maximum estimated propensity score",
    "ps_near_extreme_n": "Patients with near-extreme estimated propensity scores",
    "weighted_ess_treated": "Weighted effective sample size: treated",
    "weighted_ess_control": "Weighted effective sample size: comparator",
    "max_abs_smd_before": "Largest absolute standardized mean difference before weighting",
    "max_abs_smd_after": "Largest absolute standardized mean difference after weighting",
    "n_abs_smd_gt_0_1_after": "Covariates with absolute standardized mean difference over 0.10 after weighting",
    "complete_365d_n": "Patients with completed 365-day observation",
    "source": "Declared dataset source",
    "dataset_filename": "Input filename",
    "data_origin": "Input data provenance",
    "bootstrap_replicates_requested": "Number of bootstrap resamples requested",
    "ps_predictions_clipped_n": "Propensity scores truncated to numeric bounds",
    "bootstrap_fit_failures": "Bootstrap model-fitting failures",
}


def _word_field_label(name):
    labels = {
        "age": "Age (years)",
        "sex": "Recorded sex",
        "sex=F": "Female sex",
        "sex=M": "Male sex",
        "baseline_severity": "Baseline clinical severity",
        "comorbidity_index": "Comorbidity score",
        "prior_hospitalizations_365d": "Hospital admissions during 365 days before treatment",
        "followup_days": "Observed follow-up (days)",
        "treated": "Treatment assignment",
        "safety_event": "Recorded adverse event",
        "safety_time_days": "Days to adverse event or censoring",
        "clinical_response_365d": "Clinical response at 365 days",
        "hospitalizations_365d": "Hospital admissions within 365 days",
        "ed_visits_365d": "Emergency department visits within 365 days",
        "hospital_days_365d": "Hospital days within 365 days",
        "cost_365d_usd": "Health costs within 365 days (U.S. dollars)",
        "qaly_365d": "Quality-adjusted life years within 365 days",
    }
    return labels.get(str(name), str(name).replace("_", " ").capitalize())


def _word_outcome_label(value):
    s = str(value)
    s = s.replace("Joint-cohort benefit (qaly_365d)",
                  "Joint-cohort benefit in quality-adjusted life years (QALYs)")
    s = s.replace("USD per QALY gained",
                  "U.S. dollars per quality-adjusted life year (USD/QALY)")
    s = s.replace("(USD)", "(U.S. dollars)")
    s = s.replace("WTP=$", "willingness to pay of $")
    return s


def _word_measure_label(value):
    return {
        "weighted_KM_risk_difference": "Weighted Kaplan-Meier risk difference",
        "weighted_KM_risk_ratio": "Weighted Kaplan-Meier risk ratio",
        "risk_difference": "Risk difference",
        "risk_ratio": "Risk ratio",
        "weighted_mean_difference": "Weighted mean difference",
        "ICER_unstable_ratio_no_CI": "Incremental cost-effectiveness ratio (ICER)",
        "INMB_USD": "Incremental net monetary benefit (INMB)",
    }.get(str(value), (
        "Incremental net monetary benefit (INMB)" if "incremental_net_monetary_benefit" in str(value)
        else "Incremental cost-effectiveness ratio (ICER)" if "incremental_cost_effectiveness_ratio" in str(value)
        else "Weighted Kaplan-Meier risk difference" if str(value).endswith("weighted_km_risk_difference")
        else "Weighted Kaplan-Meier risk ratio" if str(value).endswith("weighted_km_risk_ratio")
        else "Weighted risk difference" if str(value).endswith("weighted_risk_difference")
        else "Weighted risk ratio" if str(value).endswith("weighted_risk_ratio")
        else "Weighted mean difference" if str(value).endswith("weighted_mean_difference")
        else str(value).replace("_", " ")))


def _word_interpret_effect(row):
    """Describe an output without claiming statistical or clinical causality."""
    domain = str(row["domain"])
    name = str(row["outcome"])
    measure = str(row["measure"])
    effect = row["effect"]
    if pd.isna(effect):
        return "The point estimate is not available."
    lo, hi = row.get("ci_low", np.nan), row.get("ci_high", np.nan)
    null = 1 if "risk_ratio" in measure else 0
    if pd.notna(lo) and pd.notna(hi):
        contains = float(lo) <= null <= float(hi)
        uncertainty = ("The 95% confidence interval includes the null value." if contains
                       else "The 95% confidence interval excludes the null value.")
    else:
        uncertainty = "A confidence interval was not estimated."
    if measure == "ICER_unstable_ratio_no_CI":
        return ("The descriptive incremental cost-effectiveness ratio (ICER) is "
                f"{_word_number(effect)} U.S. dollars per quality-adjusted life year (QALY). "
                "A negative ratio must be interpreted from the cost and benefit differences together, "
                "not as a standalone cost-per-benefit price. No ratio confidence interval is available.")
    if measure == "INMB_USD":
        return (f"{_word_outcome_label(name)}: the incremental net monetary benefit (INMB) "
                f"is {_word_number(effect)} U.S. dollars. {uncertainty} "
                "This alone does not establish reimbursement value.")
    if measure.endswith("risk_difference") or measure == "risk_difference":
        what = "adverse-event risk" if domain == "safety" else "clinical response probability"
        return (f"The treated-minus-comparator difference in {what} is "
                f"{float(effect) * 100:+.2f} percentage points. {uncertainty}")
    if measure.endswith("risk_ratio") or measure == "risk_ratio":
        return (f"The treated-to-comparator risk ratio is {_word_number(effect)}. {uncertainty}")
    if domain == "HCRU":
        return (f"{name}: an adjusted treated-minus-comparator difference of "
                f"{_word_number(effect)} per patient over 365 days. {uncertainty}")
    if domain == "economics" and "cost" in name.lower():
        return (f"{name}: estimated difference of {_word_number(effect)} U.S. dollars per patient. "
                f"{uncertainty}")
    return (f"{_word_outcome_label(name)}: adjusted treated-minus-comparator difference "
            f"{_word_number(effect)}. {uncertainty}")


def generate_word_report(output_dir, source=None):
    """Build a source-specific report from existing result tables and Markdown.

    The source CSV and Markdown files remain unchanged. The report includes every
    adjusted effect and design indicator, but omits runtime dependency versions.
    Sections 1-4 are tables/metrics only (no Interpretation blocks); Section 6
    retains its separate missing-data interpretation and Sections 8-9 the cautions.
    """
    from docx import Document
    from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Inches, Pt, RGBColor

    output_dir = Path(output_dir)
    required_files = ["effect_estimates.csv", "table1_baseline.csv", "cohort_summary.csv",
                      "covariate_balance.csv", "missingness.csv", "diagnostics.json",
                      "evidence_summary.md"]
    absent = [name for name in required_files if not (output_dir / name).is_file()]
    if absent:
        raise FileNotFoundError(
            f"Cannot create Word report; missing files in {output_dir.resolve()}: "
            f"{', '.join(absent)}. First run the 'analyze' command.")
    effects = pd.read_csv(output_dir / "effect_estimates.csv")
    baseline = pd.read_csv(output_dir / "table1_baseline.csv")
    cohort = pd.read_csv(output_dir / "cohort_summary.csv")
    balance = pd.read_csv(output_dir / "covariate_balance.csv")
    missingness = pd.read_csv(output_dir / "missingness.csv")
    diagnostics_data = json.loads((output_dir / "diagnostics.json").read_text(encoding="utf-8"))
    summary_text = (output_dir / "evidence_summary.md").read_text(encoding="utf-8")
    source = source or diagnostics_data.get("source", "Other")
    label = diagnostics_data.get("dataset_filename") or _word_extract_label(summary_text)
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Inches(8.5), Inches(11)
    sec.left_margin, sec.right_margin = Inches(0.60), Inches(0.60)
    sec.top_margin, sec.bottom_margin = (Inches(0.61), Inches(0.47)) if source == "NARCOMS" else (Inches(0.66), Inches(0.62))
    sec.header_distance, sec.footer_distance = Inches(0.3), Inches(0.33)
    normal = doc.styles["Normal"]
    normal.font.name = "Aptos"
    normal.font.size = Pt(9)
    normal.font.color.rgb = RGBColor(34, 51, 65)
    normal.paragraph_format.space_after = Pt(5)
    normal.paragraph_format.line_spacing = 1.10
    doc.styles["List Bullet"].paragraph_format.space_after = Pt(2)
    for name, size, before, after in (("Title", 17, 0, 9),
                                      ("Heading 1", 12, 13, 5),
                                      ("Heading 2", 10.1, 8, 3)):
        style = doc.styles[name]
        style.font.name = "Aptos"
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(27, 65, 82)
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.keep_with_next = True
    title_ppr = doc.styles["Title"].element.get_or_add_pPr()
    title_border = title_ppr.find(qn("w:pBdr"))
    if title_border is not None:
        title_ppr.remove(title_border)

    def paragraph(text, bold=False):
        p = doc.add_paragraph()
        run = p.add_run(str(text))
        run.bold = bold
        return p

    def heading(text, level=1):
        return doc.add_heading(text, level=level)

    def interpretation(items):
        if not items:
            return
        p = doc.add_paragraph()
        p.paragraph_format.keep_with_next = True
        p.paragraph_format.space_before = Pt(4)
        p.add_run("Interpretation").bold = True
        for item in items:
            doc.add_paragraph(str(item), style="List Bullet")

    def table(headers, rows, widths, *, body_pt=7.5, header_pt=7.3):
        assert len(headers) == len(widths)
        tbl = doc.add_table(rows=1, cols=len(headers))
        tbl.autofit = False
        tbl.style = "Table Grid"
        for j, width in enumerate(widths):
            tbl.columns[j].width = Inches(width)
        for j, head in enumerate(headers):
            cell = tbl.rows[0].cells[j]
            cell.width = Inches(widths[j])
            cell.text = str(head)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            shade = OxmlElement("w:shd")
            shade.set(qn("w:fill"), "E8F0F5")
            cell._tc.get_or_add_tcPr().append(shade)
            for p in cell.paragraphs:
                p.paragraph_format.space_before = Pt(2)
                p.paragraph_format.space_after = Pt(2)
                for run in p.runs:
                    run.font.size = Pt(header_pt)
                    run.font.bold = True
                    run.font.color.rgb = RGBColor(22, 51, 67)
        tr = tbl.rows[0]._tr
        tr_pr = tr.get_or_add_trPr()
        h = OxmlElement("w:tblHeader")
        h.set(qn("w:val"), "true")
        tr_pr.append(h)
        # Body cells remain unshaded, giving the tables a restrained appearance.
        for record in rows:
            if len(record) != len(headers):
                raise ValueError("Report table row length differs from column count")
            cells = tbl.add_row().cells
            for j, item in enumerate(record):
                cells[j].width = Inches(widths[j])
                cells[j].text = str(item)
                cells[j].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
                for p in cells[j].paragraphs:
                    p.paragraph_format.space_before = Pt(1)
                    p.paragraph_format.space_after = Pt(1)
                    for run in p.runs:
                        run.font.size = Pt(body_pt)
        for row_index, row in enumerate(tbl.rows):
            tr_pr = row._tr.get_or_add_trPr()
            no_split = OxmlElement("w:cantSplit")
            tr_pr.append(no_split)
            if len(rows) <= 6 and row_index < len(rows):
                # Keep compact tables (including their header) on one page.
                for cell in row.cells:
                    for cell_para in cell.paragraphs:
                        cell_para.paragraph_format.keep_with_next = True
        spacer = doc.add_paragraph()
        spacer.paragraph_format.space_after = Pt(1)
        return tbl

    def add_bullets(items, compact=False):
        for item in items:
            entry = doc.add_paragraph(_word_expand_terms(item), style="List Bullet")
            if compact:
                # Keep long safety cautions legible without orphaning glossary rows.
                entry.paragraph_format.space_after = Pt(0)
                entry.paragraph_format.line_spacing = 1.02
                for run in entry.runs:
                    run.font.size = Pt(8.5)

    title = doc.add_paragraph(style="Title")
    title.add_run(f"Sample Size and Economic Analysis for RWE Observational Studies - {source}")
    paragraph("Module: exploratory propensity-score weighted comparative analysis; sample-size planning uses a separate, explicitly parameterized pilot-based module.")
    paragraph(f"Dataset: {Path(label).name}  |  Estimand: average treatment effect in the overlap population "
              f"(ATO)" if str(diagnostics_data.get("estimand")) == "ATO" else
              f"Dataset: {Path(label).name}  |  Estimand: {diagnostics_data.get('estimand', 'Unknown')}")
    paragraph(f"Follow-up horizon: {_word_int(diagnostics_data.get('horizon_days'))} days. "
              f"Requested bootstrap resamples (B): {_word_int(diagnostics_data.get('bootstrap_replicates_requested'))}.")
    paragraph("Results are exploratory, not validated clinical conclusions or regulatory evidence. "
              "Causal interpretation requires a prespecified target trial, validated endpoints, "
              "adequate confounding adjustment, positivity and defensible censoring assumptions.")

    heading("1. Study population and design diagnostics")
    stats = [
        ("Eligible patients", _word_int(diagnostics_data.get("sample_n"))),
        ("Treated patients", _word_int(diagnostics_data.get("treated_n"))),
        ("Active-comparator patients", _word_int(diagnostics_data.get("control_n"))),
        ("Bhattacharyya propensity-score overlap coefficient", _word_number(diagnostics_data.get("phi_bhattacharyya"), 3)),
        ("Effective sample size (ESS), treated after weighting", _word_number(diagnostics_data.get("weighted_ess_treated"), 3)),
        ("Effective sample size (ESS), comparator after weighting", _word_number(diagnostics_data.get("weighted_ess_control"), 3)),
        ("Largest absolute standardized mean difference (SMD) after weighting", _word_number(diagnostics_data.get("max_abs_smd_after"), 3)),
        ("Patients with complete 365-day follow-up", _word_int(diagnostics_data.get("complete_365d_n"))),
    ]
    table(["Design indicator", "Result"], stats, [5.6, 1.7], body_pt=8.1, header_pt=8)
    paragraph("The Bhattacharyya propensity-score overlap coefficient (phi) measures the similarity "
              "between the treated and comparator propensity-score distributions. Its theoretical "
              "range is 0 (no overlap) to 1 (identical distributions / complete overlap). "
              "Higher phi indicates better common support, not proof of no confounding. "
              "The estimated coefficient uses baseline propensity scores and the treated allocation "
              "fraction; interpretation requires inspection of score tails and covariate balance. "
              "See Reference [1].")
    heading("2. Cohort characteristics and event counts (cohort_summary.csv)")
    table(["Treatment arm", "Patients", "Full-year patients", "Mean age (years)",
           "Raw adverse events", "Median follow-up (days)"], [
        ["Active comparator" if r["group"] == "control" else "Treated",
         _word_int(r["patients"]), _word_int(r["full_365d_followup"]),
         _word_number(r["baseline_mean_age"], 3), str(r["adverse_events"]),
         _word_number(r["median_followup_days"], 3)]
        for _, r in cohort.iterrows()], [1.33, .69, 1.12, 1.17, 1.15, 1.84], body_pt=7.8)
    heading("3. Adjusted safety, clinical effectiveness, utilization and economic results")
    p_intro = paragraph("For differences, effect = treated minus active comparator; for ratios, effect = treated / comparator. "
                        "The 95% confidence interval (CI) is a percentile interval derived from patient-level bootstrap resampling "
                        "with re-estimation of propensity scores; it does not account for multisite clustering. "
                        "The treated and control sample sizes are the counts contributing to each outcome.")
    p_intro.paragraph_format.keep_with_next = True
    domain_labels = {"safety": "Safety", "effectiveness": "Clinical effectiveness",
                     "HCRU": "Healthcare resource utilization (HCRU)", "economics": "Economic value"}
    for domain in ("safety", "effectiveness", "HCRU", "economics"):
        dfx = effects[effects["domain"] == domain]
        if dfx.empty:
            continue
        heading(domain_labels[domain], 2)
        rows = []
        for _, r in dfx.iterrows():
            ci = (f"[{_word_number(r['ci_low'])}, {_word_number(r['ci_high'])}]" if
                  pd.notna(r.get("ci_low")) and pd.notna(r.get("ci_high")) else "Not estimated")
            rows.append([_word_outcome_label(r["outcome"]), _word_measure_label(r["measure"]),
                         _word_number(r["treated"]), _word_number(r["control"]),
                         _word_number(r["effect"]), ci,
                         _word_int(r["n_treated"]), _word_int(r["n_control"])])
        table(["Outcome", "Effect measure", "Treated value", "Control value", "Effect",
               "95% CI", "Treated N", "Control N"], rows,
              [1.47, 1.21, .72, .72, .72, 1.28, .59, .59], body_pt=7.15, header_pt=7.15)
        counts = sorted({_word_int(v) for v in dfx["n_bootstrap"].dropna().values})
        if len(counts) == 1:
            paragraph(f"Bootstrap resamples (B) for this table: {counts[0]}.")
        elif counts:
            paragraph("Bootstrap resamples (B), by row: " + "; ".join(
                f"{_word_outcome_label(r['outcome'])}: {_word_int(r['n_bootstrap'])}" for _, r in dfx.iterrows()))
        else:
            paragraph("Bootstrap resamples (B): not available in this output.")
        # Preserve the additional probability-positive indicator outside the CI cells.
        ppos = dfx[dfx["bootstrap_probability_positive"].notna()]
        if len(ppos):
            paragraph("Bootstrap P(INMB > 0) by willingness-to-pay: " + "; ".join(
                f"{_word_outcome_label(r['outcome']).split('of $')[-1]}: "
                f"{_word_number(r['bootstrap_probability_positive'], 3)}"
                for _, r in ppos.iterrows()) + ".")
        if domain == "economics":
            paragraph("Economic definitions and reconciliation: displayed values have three "
                      "decimals; calculations retain full precision. "
                      "ICER = (treated cost - control cost) / (treated benefit - control benefit), "
                      "so separate arm-specific ICER values are undefined. "
                      "Net monetary benefit (NMB) for each arm = willingness-to-pay x benefit "
                      "minus cost; INMB = treated NMB - control NMB. A negative ICER alongside "
                      "lower costs and higher benefit denotes dominance, not a negative price "
                      "per QALY. Real-data analyses require validated costs and benefits.")
    if not (effects["domain"] == "economics").any():
        paragraph("Economic outcomes were not estimated. Missing costs or quality-adjusted life years (QALYs) "
                  "must not be interpreted as zero financial impact.")

    heading("4. Baseline characteristics (table1_baseline.csv)")
    p_baseline = paragraph("Unweighted group means (control raw [C raw], treated raw [T raw]) and propensity-score-weighted "
                           "group means (control weighted [C wtd], treated weighted [T wtd]) are provided. "
                           "The standardized mean difference (SMD) is calculated before weighting (SMD raw) "
                           "and after weighting (SMD wtd). For a categorical indicator, the mean is the group proportion.")
    p_baseline.paragraph_format.keep_with_next = True
    table(["Baseline variable", "Type", "Control raw", "Treated raw", "Control weighted",
           "Treated weighted", "SMD before", "SMD after"], [
        [_word_field_label(r["variable"]),
         ("Share" if str(r["type"]).lower() == "proportion" else str(r["type"])),
         _word_number(r["control_mean_unadjusted"], 3),
         _word_number(r["treated_mean_unadjusted"], 3),
         _word_number(r["control_mean_weighted"], 3),
         _word_number(r["treated_mean_weighted"], 3),
         _word_number(r["smd_unadjusted"], 3),
         _word_number(r["smd_weighted"], 3)] for _, r in baseline.iterrows()],
          [1.76, .59, .80, .80, .81, .81, .87, .86], body_pt=7.25, header_pt=7.15)
    table(["Baseline variable", "Missing patients (N)", "Control observed (N)",
           "Treated observed (N)"], [
        [_word_field_label(r["variable"]), _word_int(r["missing_n"]),
         _word_int(r["control_n_observed"]), _word_int(r["treated_n_observed"])]
        for _, r in baseline.iterrows()], [3.0, 1.5, 1.4, 1.4], body_pt=8)
    heading("5. Baseline covariate balance (covariate_balance.csv)")
    table(["Baseline covariate", "SMD before weighting", "SMD after weighting",
           "Absolute SMD after weighting"], [
        [_word_field_label(r["covariate"]), _word_number(r["smd_before"], 3),
         _word_number(r["smd_after"], 3), _word_number(abs(r["smd_after"]), 3)]
        for _, r in balance.iterrows()], [2.7, 1.55, 1.50, 1.55], body_pt=7.9, header_pt=7.8)
    paragraph("Propensity-score overlap: phi = 0 indicates no common support; phi = 1 "
              "indicates identical treated and control score distributions. The 'Unadjusted' "
              "points on the standardized mean difference (SMD) chart compare raw treated "
              "and control means/proportions; 'Weighted' points compare means/proportions "
              "after propensity-score weighting. Each SMD is the treated-minus-control "
              "difference divided by the pooled standard deviation. For the ATO estimand, "
              "weights are 1 - e(X) for treated patients and e(X) for comparators, where "
              "e(X) is the estimated probability of treatment given pretreatment covariates; "
              "ATE and ATT use their own weighting formulas. The dashed |SMD| = 0.100 "
              "line is a conventional screening threshold: values below 0.100 indicate "
              "reasonably balanced measured covariates, NOT absence of unmeasured confounding. "
              "See References [2] and [3].")
    for image_name, caption in (("propensity_overlap.png", "Propensity score (PS) distributions by treatment"),
                                ("covariate_balance.png", "Standardized mean differences (SMDs) before and after weighting")):
        image_path = output_dir / image_name
        if image_path.is_file():
            cap = doc.add_paragraph()
            cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            cap.paragraph_format.keep_with_next = True
            cap.add_run(caption).bold = True
            try:
                doc.add_picture(str(image_path), width=Inches(6.0))
            except Exception as exc:
                paragraph(f"Could not embed diagnostic figure {image_name}: {exc}")

    heading("6. Data completeness and missingness (missingness.csv)")
    table(["Source variable", "Missing patients (N)", "Missing fraction"], [
        [_word_field_label(r["variable"]), _word_int(r["missing_n"]),
         _word_number(r["missing_fraction"], 3)] for _, r in missingness.iterrows()],
          [3.55, 1.65, 2.10], body_pt=7.8)
    interpretation(["Missing follow-up is not automatically missing at random. The current analysis "
                    "uses complete cases for 365-day clinical response, utilization and cost outcomes; "
                    "inverse probability of censoring weighting (IPCW) or other sensitivity analyses "
                    "may be necessary for informative loss to follow-up."])

    heading("7. Full diagnostic metadata (diagnostics.json)")
    paragraph("Each indicator is shown using its full name followed by the exact machine-readable field "
              "code in parentheses.")
    diag_rows = []
    for key, val in diagnostics_data.items():
        if key in ("warnings_and_limitations", "dependencies"):
            continue
        if key == "missingness_fraction":
            if isinstance(val, dict):
                for field, fraction in val.items():
                    diag_rows.append([f"Missingness fraction for {_word_field_label(field)} "
                                      f"(missingness_fraction.{field})", _word_number(fraction, 3)])
            continue
        if isinstance(val, float):
            val = _word_number(val, 3)
        elif isinstance(val, int):
            val = _word_int(val)
        diag_rows.append([f"{_WORD_INDICATORS.get(key, key.replace('_', ' ').capitalize())} ({key})", str(val)])
    table(["Diagnostic indicator and source field", "Value"], diag_rows,
          [5.55, 1.75], body_pt=7.6, header_pt=8.0)

    caution_heading = heading("8. Assumptions, cautions and disclosure")
    if source == "NARCOMS":
        # Keep the source-specific cautions and glossary together on a new page.
        caution_heading.paragraph_format.page_break_before = True
    if source == "EpicCosmos":
        paragraph("Epic Cosmos-specific caution: multisite clustering and patient migration between "
                  "health systems may affect outcome capture and precision. Access to, and export "
                  "of, individual records is subject to Epic and institutional agreements.")
    elif source == "NARCOMS":
        paragraph("North American Research Committee on Multiple Sclerosis (NARCOMS) Registry-specific "
                  "caution: survey-wave timing, self-report, recall, treatment changes and participant "
                  "dropout may complicate inference. Exact adverse-event dates, billing costs and "
                  "device exposure may not be available.")
    cautions = _word_summary_items(summary_text, "Interpretation and constraints")
    cautions += _word_summary_items(summary_text, "Warnings")
    cautions += diagnostics_data.get("warnings_and_limitations", [])
    cautions += [
        "Code lists and treatment, device or outcome algorithms must be validated by clinicians and data owners.",
        "NARCOMS disease-modifying therapy (DMT) comparisons cannot be characterized as device evidence "
        "unless actual device exposure is measured and verified.",
        "All tabulated results, figures and reports must pass data use agreement (DUA), institutional "
        "review board (IRB) and small-cell disclosure review before external release.",
    ]
    # Consolidate substantively duplicated cautions already explained earlier
    # in the report. Their contents remain covered by the more specific bullets.
    duplicate_leads = (
        "all results are exploratory until",  # documented in the title and status
        "complete-case 365-day outcomes assume",  # covered by longer missingness caveat
        "km safety risk assumes",  # covered by the time-to-event caveat above
        "narcoms disease-modifying therapy (dmt) comparisons cannot",  # same device restriction above
        "all tabulated results, figures and reports must pass",  # same source DUA/IRB rule above
    )
    unique = []
    seen = set()
    for item in cautions:
        canonical = str(item).strip().lower().rstrip(".")
        if canonical.startswith(duplicate_leads):
            continue
        if canonical not in seen:
            seen.add(canonical)
            unique.append(item)
    add_bullets(unique, compact=True)

    heading("9. Abbreviations")
    abbreviation_rows = [
        ("ATE", "Average treatment effect"),
        ("ATT", "Average treatment effect among treated patients"),
        ("ATO", "Average treatment effect in the overlap population"),
        ("ESS", "Effective sample size"),
        ("SMD", "Standardized mean difference"),
        ("PS", "Propensity score"),
        ("CI", "Confidence interval"),
        ("B", "Number of bootstrap resamples"),
        ("KM", "Kaplan-Meier survival estimator"),
        ("HCRU", "Healthcare resource utilization"),
        ("QALY", "Quality-adjusted life year"),
        ("ICER", "Incremental cost-effectiveness ratio"),
        ("INMB", "Incremental net monetary benefit"),
        ("WTP", "Willingness to pay"),
        ("DUA", "Data use agreement"),
        ("IRB", "Institutional review board"),
        ("IPCW", "Inverse probability of censoring weighting"),
        ("DMT", "Disease-modifying therapy"),
        ("C raw / T raw", "Control / treated unweighted group means"),
        ("C wtd / T wtd", "Control / treated propensity-score-weighted group means"),
        ("SMD before / after", "Standardized mean difference before / after weighting"),
        ("N", "Number of patients included in a specified analysis"),
        ("USD", "United States dollars"),
    ]
    midpoint = (len(abbreviation_rows) + 1) // 2
    left_terms = abbreviation_rows[:midpoint]
    right_terms = abbreviation_rows[midpoint:]
    glossary = []
    for i, (left_abbrev, left_expanded) in enumerate(left_terms):
        if i < len(right_terms):
            right_abbrev, right_expanded = right_terms[i]
        else:
            right_abbrev, right_expanded = "", ""
        glossary.append([left_abbrev, left_expanded, right_abbrev, right_expanded])
    table(["Shorthand", "Full term", "Shorthand", "Full term"], glossary,
          [.95, 2.65, 1.05, 2.65], body_pt=7.4, header_pt=7.6)

    plan_file = output_dir / "sample_size_plan.csv"
    if plan_file.is_file():
        plan = pd.read_csv(plan_file)
        if len(plan) == 1:
            r = plan.iloc[0]
            planning_heading = heading("10. Prospective sample-size planning")
            planning_heading.paragraph_format.page_break_before = True
            paragraph("This optional module estimates required enrollment using the pilot "
                      "weighted estimator's empirical bootstrap standard error and a "
                      "prespecified clinically meaningful effect (0.100 is illustrative "
                      "in the synthetic example). Consider using the Duke PSpower analytical "
                      "model (https://www2.stat.duke.edu/~fl35/pspower.html) to verify the "
                      "sample size calculations. Because the two methods have different "
                      "variance and outcome-model assumptions, compare their results only "
                      "after aligning the estimand, overlap, effect size and test specifications. "
                      "See Reference [4].")
            planning_rows = [
                ("Target outcome", r["outcome"]),
                ("Measure", _word_measure_label(r["metric_code"])),
                ("Causal target (estimand)", r["estimand"]),
                ("Pilot eligible participant count", _word_int(r["pilot_n"])),
                ("Pilot bootstrap standard error", _word_number(r["pilot_bootstrap_standard_error"])),
                ("Pilot asymptotic variance coefficient V", _word_number(r["asymptotic_variance_coefficient"])),
                ("Prespecified target-effect difference", _word_number(r["target_effect_difference"])),
                ("Significance level (alpha)", _word_number(r["alpha"])),
                ("Target statistical power", _word_number(r["target_power"])),
                ("Test sidedness", r["test"]),
                ("Additional design variance multiplier", _word_number(r["design_multiplier"])),
                ("Extra loss beyond pilot cohort", _word_number(r["extra_loss_fraction_beyond_pilot"])),
                ("Required enrolled patients", _word_int(r["required_n_total_enrolled"])),
            ]
            table(["Planning input or output", "Estimate"], planning_rows,
                  [5.5, 1.8], body_pt=8.0)
            paragraph("Prospective planning assumption: Var(weighted effect estimate) approximately "
                      "equals V/N, where V is calibrated from the pilot's refitted-propensity "
                      "bootstrap. Sample size increases with the square of the inverse of the "
                      "prespecified effect. The existing pilot's incomplete follow-up is already "
                      "reflected in V; apply 'extra loss' only to losses expected beyond the pilot.")

    # Centralized methodological and reporting sources. Section narratives use
    # short reference numbers; complete citations appear here, after any planning.
    heading("11. References" if plan_file.is_file() else "10. References")
    reference_items = [
        "[1] Liu B, Yang C, Li F (2026). Sample size and power calculations "
        "for causal inference in observational studies. Annals of Statistics, "
        "54(4), 2191-2208. https://doi.org/10.1214/26-AOS2649",
        "[2] Austin PC (2009). Balance diagnostics for comparing the distribution "
        "of baseline covariates between treatment groups in propensity-score "
        "matched samples. Statistics in Medicine, 28, 3083-3107. "
        "https://doi.org/10.1002/sim.3697",
        "[3] Li F, Morgan KL, Zaslavsky AM (2018). Balancing covariates "
        "via propensity score weighting. Journal of the American Statistical "
        "Association, 113(521), 390-400. "
        "https://doi.org/10.1080/01621459.2016.1260466",
        "[4] Duke University, PSpower interactive online calculator: Sample size "
        "and power for causal inference in observational studies. "
        "https://www2.stat.duke.edu/~fl35/pspower.html",
    ]
    reporting_sources = _word_summary_items(summary_text, "Relevant reporting guidance")
    reference_items += [f"[{i}] {item}" for i, item in enumerate(reporting_sources, start=5)]
    for citation in reference_items:
        p_ref = doc.add_paragraph()
        p_ref.paragraph_format.left_indent = Inches(0.27)
        p_ref.paragraph_format.first_line_indent = Inches(-0.27)
        p_ref.paragraph_format.space_after = Pt(5)
        p_ref.add_run(citation)

    footer = sec.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    footer.add_run("Observational RWE analysis  -  " + source + "  -  Page ")
    field = OxmlElement("w:fldSimple")
    field.set(qn("w:instr"), "PAGE")
    footer._p.append(field)
    path = output_dir / f"{source}_RWE_Evidence_Report.docx"
    doc.save(path)
    return path



def is_verified_legacy_demo_duplicate(legacy_path: Path, synthetic_path: Path) -> bool:
    """Remove only byte-identical synthetic copies left by older make-demo.

    A user's real/study-specific canonical file must never be deleted, even
    when make-demo is rerun with a different seed or n. If verification is
    impossible, retain the file and tell the user why there are >2 CSVs.
    """
    import filecmp

    if (legacy_path.is_symlink() or synthetic_path.is_symlink()
            or not legacy_path.is_file() or not synthetic_path.is_file()):
        return False
    if not filecmp.cmp(legacy_path, synthetic_path, shallow=False):
        return False
    try:
        patient_ids = pd.read_csv(legacy_path, usecols=["patient_id"],
                                  dtype={"patient_id": "string"})["patient_id"]
        return bool(len(patient_ids) and patient_ids.notna().all()
                    and patient_ids.str.startswith("SYNTHETIC_").all())
    except (ValueError, OSError, pd.errors.ParserError):
        return False


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    gen = sub.add_parser("make-demo", help="Generate clearly labeled synthetic data")
    gen.add_argument("--output-dir", default="demo_data")
    gen.add_argument("--n", type=int, default=900)
    gen.add_argument("--seed", type=int, default=2026)
    an = sub.add_parser("analyze", help="Analyze ONE normalized authorized source")
    an.add_argument("--input", required=True, help="A single-source patient-level CSV")
    an.add_argument("--source", required=True, choices=["EpicCosmos", "NARCOMS", "Other"])
    an.add_argument("--output-dir", required=True)
    an.add_argument("--estimand", default="ATO", choices=["ATO", "ATE", "ATT"])
    an.add_argument("--horizon", default=365, type=int,
                    help="Must match outcome column definitions (defaults to 365)")
    an.add_argument("--bootstrap", default=300, type=int,
                    help=">=500 recommended for stable publication intervals")
    an.add_argument("--seed", default=2026, type=int)
    an.add_argument("--cost-basis", default="none", choices=["none", "validated_costs"],
                    help="Never use charges as costs without explicit transformation/validation")
    an.add_argument("--thresholds", type=float, nargs="+", default=[50000, 100000, 150000],
                    help="USD/QALY thresholds, only when validated QALY exists")
    an.add_argument("--responder-wtp", type=float, nargs="+", default=[],
                    help="Optional, explicitly chosen USD/additional responder thresholds")
    an.add_argument("--min-cell", type=int, default=11,
                    help="Minimal suppression demonstration; actual DUA may require more")
    an.add_argument("--no-word-report", action="store_true",
                    help="Skip automatic DOCX generation (default: generate it)")
    plan = sub.add_parser("plan-sample-size", help="Pilot-bootstrap normal approximation for required N")
    plan.add_argument("--output-dir", required=True, help="Directory containing analyzed pilot outputs")
    plan.add_argument("--metric-code", required=True, help="Unique metric_code in effect_estimates.csv")
    plan.add_argument("--target-effect", required=True, type=float,
                      help="Clinically prespecified absolute difference in outcome units")
    plan.add_argument("--alpha", type=float, default=0.05)
    plan.add_argument("--power", type=float, default=0.80)
    plan.add_argument("--test", choices=["two-sided", "one-sided"], default="two-sided")
    plan.add_argument("--extra-loss", type=float, default=0.0,
                      help="Additional loss beyond the pilot's observed losses")
    plan.add_argument("--design-multiplier", type=float, default=1.0,
                      help="Additional variance inflation factor >=1 for sensitivity planning")
    plan.add_argument("--no-word-refresh", action="store_true")
    rep = sub.add_parser("report", help="Build Word report from existing output tables, without recomputation")
    rep.add_argument("--output-dir", required=True,
                     help="Directory containing an earlier analysis' output CSV/JSON/MD files")
    rep.add_argument("--source", choices=["EpicCosmos", "NARCOMS", "Other"],
                     help="Optional source label; normally read from diagnostics.json")
    args = p.parse_args()
    if args.command == "plan-sample-size":
        plan_sample_size_from_pilot(args.output_dir, args.metric_code, args.target_effect,
                                    alpha=args.alpha, power=args.power,
                                    sidedness=args.test, additional_loss=args.extra_loss,
                                    design_multiplier=args.design_multiplier,
                                    refresh_word=not args.no_word_refresh)
        return
    if args.command == "report":
        word_path = generate_word_report(args.output_dir, source=args.source)
        print(f"Word report created: {word_path.resolve()}")
        return
    if args.command == "make-demo":
        dest = Path(args.output_dir)
        dest.mkdir(parents=True, exist_ok=True)
        for i, source in enumerate(("EpicCosmos", "NARCOMS")):
            synthetic_path = dest / f"SYNTHETIC_{source}_cohort.csv"
            legacy_path = dest / f"{source}_cohort.csv"

            # Previous versions created TWO identical CSVs per source. Remove
            # only a proven synthetic duplicate, BEFORE regenerating the demo
            # (which may be requested with a different seed). Never remove a
            # canonical file containing different data or a symbolic link.
            if legacy_path.exists() or legacy_path.is_symlink():
                if is_verified_legacy_demo_duplicate(legacy_path, synthetic_path):
                    legacy_path.unlink()
                    print(f"Removed verified duplicate synthetic placeholder: {legacy_path.resolve()}")
                else:
                    print(f"Preserved distinct input (not modified): {legacy_path.resolve()}")

            sample = demo(args.n, args.seed + i, source)
            sample.to_csv(synthetic_path, index=False)
            print(f"Generated synthetic {source} cohort: {synthetic_path.resolve()}")
        print("Demo generation complete: one generated CSV per source (two total).")
        return
    if args.bootstrap < 20:
        raise SystemExit("At least 20 bootstrap replicates are required (>=500 preferred)")
    if args.horizon != 365:
        raise SystemExit("This template uses 365-day outcome columns; implement matching "
                         "redefinitions before changing --horizon")
    data = pd.read_csv(args.input, dtype={"patient_id": "string"})
    run_analysis(data, args)


if __name__ == "__main__":
    main()
