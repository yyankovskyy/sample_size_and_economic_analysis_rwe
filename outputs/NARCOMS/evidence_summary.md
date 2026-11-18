# Exploratory RWE evidence summary — NARCOMS

**Not a validated regulatory submission or a peer-reviewed analysis.**

Dataset: `SYNTHETIC_NARCOMS_cohort.csv` | Estimand: **ATO** | 
Horizon: **365 days** | Bootstraps: **1000**

## Study population and diagnostics

- Participants: 450; treated: 197; control: 253.
- Bhattacharyya propensity score overlap: 0.910.
- Weighted ESS: treated 171.035, control 204.392.
- Maximum weighted absolute SMD: 0.007.

## Adjusted estimates (treated versus active comparator)

| Domain | Outcome | Measure | Treated | Control | Effect | 95% bootstrap CI |
|---|---|---|---:|---:|---:|---:|
| safety | Adverse-event risk by day 365 | adverse_event_365d_weighted_km_risk_difference | 0.091 | 0.081 | 0.010 | [-0.050, 0.069] |
| safety | Adverse-event risk by day 365 | adverse_event_365d_weighted_km_risk_ratio | 0.091 | 0.081 | 1.121 | [0.534, 2.391] |
| effectiveness | Clinical response by day 365 | clinical_response_365d_weighted_risk_difference | 0.573 | 0.533 | 0.040 | [-0.066, 0.150] |
| effectiveness | Clinical response by day 365 | clinical_response_365d_weighted_risk_ratio | 0.573 | 0.533 | 1.076 | [0.891, 1.313] |
| HCRU | Hospital admissions in 365 days | hospital_admissions_365d_weighted_mean_difference | 0.650 | 0.746 | -0.096 | [-0.284, 0.092] |
| HCRU | Emergency department visits in 365 days | emergency_visits_365d_weighted_mean_difference | 0.818 | 0.935 | -0.117 | [-0.306, 0.079] |
| HCRU | Hospital days in 365 days | hospital_days_365d_weighted_mean_difference | 1.696 | 1.823 | -0.127 | [-0.647, 0.357] |

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
