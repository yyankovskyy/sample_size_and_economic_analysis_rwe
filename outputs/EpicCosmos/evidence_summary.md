# Exploratory RWE evidence summary — EpicCosmos

**Not a validated regulatory submission or a peer-reviewed analysis.**

Dataset: `SYNTHETIC_EpicCosmos_cohort.csv` | Estimand: **ATO** | 
Horizon: **365 days** | Bootstraps: **1000**

## Study population and diagnostics

- Participants: 450; treated: 246; control: 204.
- Bhattacharyya propensity score overlap: 0.933.
- Weighted ESS: treated 210.503, control 182.659.
- Maximum weighted absolute SMD: 0.005.

## Adjusted estimates (treated versus active comparator)

| Domain | Outcome | Measure | Treated | Control | Effect | 95% bootstrap CI |
|---|---|---|---:|---:|---:|---:|
| safety | Adverse-event risk by day 365 | adverse_event_365d_weighted_km_risk_difference | 0.082 | 0.134 | -0.052 | [-0.120, 0.005] |
| safety | Adverse-event risk by day 365 | adverse_event_365d_weighted_km_risk_ratio | 0.082 | 0.134 | 0.611 | [0.319, 1.057] |
| effectiveness | Clinical response by day 365 | clinical_response_365d_weighted_risk_difference | 0.579 | 0.481 | 0.099 | [-0.006, 0.207] |
| effectiveness | Clinical response by day 365 | clinical_response_365d_weighted_risk_ratio | 0.579 | 0.481 | 1.205 | [0.989, 1.494] |
| HCRU | Hospital admissions in 365 days | hospital_admissions_365d_weighted_mean_difference | 0.700 | 0.939 | -0.238 | [-0.457, -0.026] |
| HCRU | Emergency department visits in 365 days | emergency_visits_365d_weighted_mean_difference | 0.766 | 1.020 | -0.254 | [-0.462, -0.032] |
| HCRU | Hospital days in 365 days | hospital_days_365d_weighted_mean_difference | 1.790 | 2.288 | -0.498 | [-1.126, 0.103] |
| economics | Healthcare costs over 365 days (USD) | health_cost_365d_usd_weighted_mean_difference | 10,520.157 | 11,827.421 | -1,307.265 | [-2,292.618, -245.059] |
| economics | Joint-cohort benefit (qaly_365d) | joint_cohort_qaly_weighted_mean_difference | 0.756 | 0.724 | 0.032 | [0.012, 0.054] |
| economics | Joint-cohort costs (USD) | joint_cohort_cost_usd_weighted_mean_difference | 10,520.157 | 11,827.421 | -1,307.265 | [-2,292.618, -245.059] |
| economics | USD per QALY gained | joint_cohort_incremental_cost_effectiveness_ratio | NA | NA | -40,497.725 | not estimated |
| economics | Incremental net monetary benefit, WTP=$50,000 | joint_cohort_incremental_net_monetary_benefit_inmb_50000 | 27,302.245 | 24,380.983 | 2,921.262 | [1,354.302, 4,488.307] |
| economics | Incremental net monetary benefit, WTP=$100,000 | joint_cohort_incremental_net_monetary_benefit_inmb_100000 | 65,124.647 | 60,589.387 | 4,535.260 | [2,161.910, 6,931.057] |
| economics | Incremental net monetary benefit, WTP=$150,000 | joint_cohort_incremental_net_monetary_benefit_inmb_150000 | 102,947.049 | 96,797.791 | 6,149.258 | [2,858.381, 9,565.688] |

## Interpretation and constraints

- Positive differences mean higher values in the treated cohort; for adverse events, higher values generally indicate greater risk.
- The safety estimate is weighted Kaplan-Meier cumulative event risk; independent censoring is assumed. For competing risks, use weighted Aalen-Johansen.
- The causal interpretation requires consistency, exchangeability given measured pretreatment covariates, positivity and valid outcome ascertainment.
- Complete-case follow-up is NOT automatically missing at random. Missing outcome and censoring mechanisms must be justified or handled with IPCW/other sensitivity analyses.
- The outcome columns must be based on validated endpoint definitions and a common clinical index date in a new-user active-comparator design.
- The ICER is only descriptive; use a cost-effectiveness plane/CEAC, prespecified willingness-to-pay thresholds, uncertainty and perspective in a formal economic study.
- A real comparative-device study requires confirmed device exposure in BOTH groups. NARCOMS DMT comparisons are not device evidence without captured device exposure.
- Outputs are aggregate, but all sharing must be reviewed under the source-specific DUA, IRB, and suppression/disclosure requirements.

## Warnings

- All results are exploratory until cohort/exposure/outcome algorithms are validated
- Index dates are syntax-checked only; confirm every exposure, confounder and outcome time window upstream
- Patient bootstrap is not site-cluster robust; use cluster-aware bootstrap for multisite inference
- Complete-case 365-day outcomes assume follow-up availability is ignorable; consider IPCW
- KM safety risk assumes independent censoring and no competing-risk complication

## Relevant reporting guidance

- FDA, *Use of Real-World Evidence to Support Regulatory Decision-Making for Medical Devices*, December 2025.
- TARGET statement (2025), target-trial emulation reporting.
- RECORD statement (2015), routinely collected healthcare data reporting.
