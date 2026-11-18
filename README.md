# Sample Size and Economic Analysis for RWE Observational Studies (v4.1)

## 1. Purpose and expected source-file inputs

This Python application analyzes **one standardized, patient-level research cohort at a time** to produce exploratory real-world evidence (RWE) for treatment/device **safety**, **comparative clinical effectiveness**, **healthcare resource utilization (HCRU)**, and — **only when cost/benefit input data support it** — economic evaluation. The default estimand is the **average treatment effect in the overlap population (ATO)**. ATE and ATT are also supported.

### Input-file expectations

- **`demo_data/EpicCosmos_cohort.csv`** is supplied as a **synthetic placeholder** in the standardized input format required by this Python application. It can be substituted with a study-specific, authorized **Epic Cosmos-derived analytical cohort** prepared by the user, **within an environment where those data are permitted**.
- **`demo_data/NARCOMS_cohort.csv`** is supplied as a **synthetic placeholder** in the same standardized format. It can be substituted with a study-specific, authorized **NARCOMS-derived analytical cohort** prepared by the user, subject to the registry's permissions and secure processing rules.

**Expected contents:** one row per eligible patient; a unique study identifier; a common clinically valid treatment/index date; treatment-versus-comparator indicator; pretreatment confounders; follow-up time; and any *actually available and validated* safety, response, resource-use or economic outcomes. **Section 6 gives the exact expected column names, human-readable labels, coding and required/optional status.** Do not invent missing clinical variables or set unobserved events to zero.

**Important:** These two placeholders model **this application's analysis-ready schema**, **not the native Epic Cosmos SQL/data model or the native NARCOMS questionnaire/survey-wave schema**. Neither synthetic file has been calibrated or validated against actual proprietary source extracts. Genuine longitudinal extracts must be transformed and their clinical definitions verified before they can replace the placeholders. Replacing the CSV does **not** guarantee study validity or one-step compatibility.

The `demo_data/SYNTHETIC_EpicCosmos_cohort.csv` and `demo_data/SYNTHETIC_NARCOMS_cohort.csv` files are *separately labeled reproducible demonstration copies*. The `make-demo` command regenerates those `SYNTHETIC_*` files and creates the shorter-name placeholder files **only if they do not exist**. It deliberately does **not** overwrite existing `EpicCosmos_cohort.csv` or `NARCOMS_cohort.csv` files, which could contain the user's replacement analytical inputs.

The program **does not connect to Epic, NARCOMS, Clarity, Caboodle, a FHIR server, or a registry download API**, and it does not automatically build source-specific clinical cohorts. All included placeholder patient records, clinical outcomes and economic values are simulated. Treat reports produced from them as software demonstrations, not evidence about real patients or healthcare effectiveness. The application does not merge, pool, link or match individual records across the two sources.

Every `analyze` run writes Markdown, a source-specific **Word (`.docx`)** report, CSV tables and plots to the chosen output directory. Displayed continuous results in the Word and Markdown reports and terminal are rounded to **three decimal places**; raw CSV/JSON retain full floating-point precision to support reproducibility. A separate `report` command produces a Word report from *already calculated outputs* without re-running propensity-score models or bootstrap iterations.

### Two distinct methodological modules (important)

1. **`analyze` — observational RWE outcome analysis.** Fits a propensity-score (PS) model, applies ATE, ATT or ATO weights, evaluates confounder balance/overlap, and estimates safety, comparative effectiveness, resource use and optional health-economic results. **The output of `analyze` alone is NOT a prospective sample-size calculation.**
2. **`plan-sample-size` — exploratory sample-size planning from an eligible pilot cohort.** Uses a bootstrap standard error obtained from `analyze`, a researcher-specified **clinically meaningful target difference**, test alpha, power, optional additional attrition beyond that already seen in the pilot, and an explicit sensitivity variance multiplier. It produces `sample_size_plan.csv`/`.md` and adds a section to the Word report. It assumes the same study design, propensity-score overlap, measured covariates and follow-up conditions as the pilot. **This is not the Liu–Yang–Li / PSpower analytical method.** For design-stage variation in the overlap and separate confounding coefficient (rho-squared), use the validated R package [PSpower](https://cran.r-project.org/web/packages/PSpower/index.html) and/or a prespecified clinical simulation study.


## 2. Installation and location of files

Requirements: Python 3.11 or newer recommended, with `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, **`python-docx`**; `pytest` is optional for tests.

From the application directory on macOS / Linux:

```bash
# Use an existing conda environment, or a new venv:
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

# Confirm where you are and where the demonstration data actually are:
pwd
ls -lh ./demo_data/
python rwe_evidence_pipeline.py --help
```

The expected project layout is:

```text
Sample_Size_and_Economic_Analysis_for_RWE_Observational_Studies/
  rwe_evidence_pipeline.py
  README.md
  requirements.txt
  test_pipeline.py
  demo_data/
    EpicCosmos_cohort.csv            # synthetic placeholder; replace only if authorized
    NARCOMS_cohort.csv               # synthetic placeholder; replace only if authorized
    SYNTHETIC_EpicCosmos_cohort.csv  # clearly labeled reproducible demo
    SYNTHETIC_NARCOMS_cohort.csv     # clearly labeled reproducible demo
  outputs/
    EpicCosmos/ ...
    NARCOMS/ ...
```

**Paths matter:** From the application directory, `./demo_data/file.csv` is a relative path. `/demo_data/file.csv` points to an unrelated directory at the filesystem root and will usually raise `FileNotFoundError`.

## 3. Expected cohort files: synthetic placeholders, replacement, and scope

**The expected application inputs are `demo_data/EpicCosmos_cohort.csv` and `demo_data/NARCOMS_cohort.csv`.** Both files shipped in this project are populated with synthetic placeholder rows so the code can run immediately. They share the analysis-ready column contract defined in **Section 6**. They are not raw database exports, faithful simulations of the native models, or validated representations of clinical event incidence.

| Placeholder path | Initial demonstration content | How to substitute a real study cohort |
|---|---|---|
| `demo_data/EpicCosmos_cohort.csv` | 450 invented patient records, treatment assignment, baseline predictors, simulated safety/response/utilization, and **artificial** costs/QALYs | Prepare a validated, one-row-per-patient Cosmos analytical cohort **inside your permitted Cosmos research environment**, using the Section 6 columns. Do not export restricted line-level Cosmos data to a local `demo_data` folder. |
| `demo_data/NARCOMS_cohort.csv` | 450 invented patient records, treatment assignment, generic baseline predictors, simulated safety/response/utilization; **no cost/QALY columns** | Obtain an approved NARCOMS extract, derive an indexed patient cohort from survey waves/DMT histories, and populate only endpoints genuinely present and valid in the approved data. |

Each canonical placeholder is initially an exact copy of its corresponding `SYNTHETIC_*` demonstration file. The first column (`patient_id`) starts with `SYNTHETIC_`, which lets the pipeline label the data origin as `synthetic_demo` even when the filename is not prefixed `SYNTHETIC_`.

### 3.1 Creating or refreshing demonstration data

```bash
python rwe_evidence_pipeline.py make-demo \
  --output-dir ./demo_data \
  --n 450 \
  --seed 2026
```

`--n 450` generates **450 synthetic patients in each data source scenario**. The generated `SYNTHETIC_*` files can be regenerated; the unprefixed `EpicCosmos_cohort.csv` and `NARCOMS_cohort.csv` files are generated **only if absent** so the command never silently overwrites substituted study cohorts. To intentionally reset a canonical placeholder to newly generated synthetic data, first back up any real extract and explicitly copy the desired `SYNTHETIC_*` file into the corresponding canonical name.

### 3.2 Scope of the illustrative content

The generator simulates age, sex, a numeric severity variable, a generic comorbidity score, preceding hospitalizations, treatment assignment via a logistic propensity model, time to a synthetic safety event, follow-up/censoring, a binary response, and fixed-year resource use. The Epic-labelled example additionally simulates costs and one-year quality-adjusted life years (QALYs). The NARCOMS-labelled example omits those economic variables and uses different **arbitrary simulation settings** for age and sex.

**What this does not demonstrate:** Cosmos SQL table structures, site-level encounter sequences, terminology mapping, actual drug/device identifiers, health-system utilization completeness; nor NARCOMS survey instruments, wave-level records, actual Patient Determined Disease Steps (PDDS) scale, recorded disease-modifying therapies (DMTs), or validated outcomes. No synthetic rate or effect estimate is a real-world estimate from either source.

### 3.3 Replacement rules: when your real extract is ready

1. Work **inside the data source's approved research environment** and retrieve the actual data dictionary, eligibility rules and recorded fields. Both source systems normally contain multiple records per patient, so define cohort construction and index dates upstream.
2. Transform approved records to **one row per eligible patient** in Section 6's schema. Use a source-specific mapping and validated endpoint logic, not simple renaming when the original data are longitudinal.
3. Validate the mandatory columns (including age, severity, comorbidity and prior utilization) and both comparison groups. If a required variable cannot be validly derived, **adapt the model/protocol and code** rather than inventing or filling unsupported clinical data.
4. Supply each available optional outcome with documented ascertainment and time windows; omit unavailable endpoints. Run `--cost-basis none` unless genuine, validated patient-level costs are available.
5. Substitute the normalized approved cohort at a **permitted local path within the source's secure environment** and pass that path to `--input`. The canonical placeholder filenames are only convenient examples; the program accepts any existing permitted CSV path. Do not transmit real Cosmos or NARCOMS records merely to reproduce these examples.
6. Run `analyze` separately for each source and review data completeness, balance and privacy disclosure before interpreting or exporting aggregate findings.

## 4. How to obtain, prepare, and use ACTUAL authorized Epic Cosmos data

**Access and permitted working environment must be established first.** Epic states that Cosmos direct access is for approved users affiliated with participating organizations, and it cannot be purchased: <https://cosmos.epic.com/request-access/>. Institutions such as [Stanford Medicine](https://med.stanford.edu/starr-tools/other-resources/cosmos/faq.html) and [UC Davis](https://health.ucdavis.edu/data/epic-cosmos.html) explicitly state **line-level Cosmos data must not be exported** from Cosmos; only reviewed aggregate results may be transferred out. Confirm the rules for *your particular affiliation and project*. Having access to query records is **not** permission to copy those records onto a personal laptop.

**Operational model for Cosmos research:** move the code into the approved **Cosmos Data Science environment** when institutional policy requires; query and construct the patient-level analytical cohort **inside that environment**, execute `analyze` there if its dependencies are permitted, and export only the aggregated, reviewed results if allowed. If your institution separately approves data use from its **local Epic Clarity/Caboodle** environment, keep that source labeled accurately; local Epic data are **not automatically Cosmos data**.

### 4.1 Step-by-step Cosmos cohort preparation

1. **Obtain approval**: institutional Cosmos access, required certification/data-scientist training, project use approval, IRB/DUA review as applicable, and an approved compute/storage location.
2. **Specify the emulated target trial**: disease/population, eligibility criteria, clinically meaningful treatment/device exposure, a new-user washout period, an active comparator, index date, follow-up, outcomes, confounders, treatment switching, censoring and target estimand (ATO/ATE/ATT).
3. **Identify source objects from the authorized Cosmos data model**: use the actual institutional dictionary, SQL catalog and terminology tables. There is **no universal public SQL schema or export script** for this application. Do not assume generic names in this README are Cosmos table names.
4. **Establish one record per eligible new user**: deterministically identify the first qualifying treatment or device exposure, choose one index date per patient, screen pre-index eligibility, and resolve duplicates across institutions according to the approved research model.
5. **Derive pretreatment features**: age/sex at index; baseline disease severity; comorbidity index from clinically validated diagnosis codes and lookback period; prior hospital admissions in the 365 days **before** index. Handle site and calendar-time differences as needed.
6. **Define exposure rigorously**: determine whether the intended specific **medical device** is actually identifiable (implant record, procedure/model/UDI when available). Procedure code alone may not distinguish device models. For treatments, verify the source reflects therapy initiation versus an order or billing claim only.
7. **Define outcomes and time windows**: adverse events from validated code lists plus clinical ascertainment; `safety_time_days` from index to event/censor; completed 365-day response, utilization and hospital days where observable. Do **not** set unobserved events to zero. Explicitly handle death, competing risks, loss of coverage or site migration.
8. **Validate economics**: Cosmos EHR encounters/charges do not automatically contain **validated patient-level total costs** or QALYs. Link to an institutionally approved cost accounting or utility source when permitted; otherwise omit the economic columns and use `--cost-basis none`.
9. **Perform ETL QA**: study flow diagram; counts before/after every filter; treatment frequencies; code-version mapping; event-date consistency; exposure uniqueness; feature missingness; distributional plausibility; and peer review by the clinical/domain lead.
10. **Create the normalized file INSIDE the permitted environment** with the schema in Section 6. `patient_id` must be de-identified/pseudonymized as permitted and unique within the source. Then run the application **in place** and follow the review/aggregation/suppression procedure before sharing any output.

### 4.2 Typical Cosmos logical source domains (not literal table names)

| Normalized field | Clinical source and derivation |
|---|---|
| `patient_id`, `index_date`, `treated` | Approved cross-encounter person identifier; validated first treatment/device initiation event; new-user active-comparator flags |
| `age`, `sex` | Demographics as available/approved at treatment index, respecting data masking |
| `baseline_severity`, `comorbidity_index` | Pre-index assessments, diagnoses/laboratory/clinical measures, validated scoring |
| `prior_hospitalizations_365d` | Pre-index admissions in prior 365 days, with observation completeness assessed |
| `followup_days` | Days from index to last reliable observation, study end, censoring or competing event, as specified |
| `safety_event`, `safety_time_days` | Validated time-to-event phenotype; event flag and event/censor duration |
| `clinical_response_365d` | Prespecified binary response with sufficient 365-day ascertainment |
| `hospitalizations_365d`, `ed_visits_365d`, `hospital_days_365d` | Fixed-horizon resource use with encounter deduplication and complete follow-up |
| `cost_365d_usd`, `qaly_365d` | **Only** if separately validated and permitted; include price year, currency, perspective and utility source |

**Important:** A commercial market-access or comparative-institution study may face additional Cosmos usage restrictions. Check the Cosmos Rules of the Road and the institution's rules before planning such work.

## 5. How to obtain, prepare, and use ACTUAL NARCOMS research data

NARCOMS is a longitudinal, primarily patient-reported MS registry, not a comprehensive EHR/claims database. The official [Researcher & Provider page](https://www.narcoms.org/researchers-providers) describes approved de-identified dataset access or registry-led analysis, proposal review, IRB documentation, agreements, secure delivery and manuscript review requirements. Contact **msregistry@narcoms.org** for feasibility, the data dictionary and an approved protocol-specific dataset.

### 5.1 Step-by-step NARCOMS cohort preparation

1. **Obtain approvals and permitted extract**: complete the NARCOMS Researcher Information Form and proposal review, provide IRB documentation as requested, sign the required agreements and receive the approved dataset through the provided secure channel. Verify re-use, linkage and publication restrictions.
2. **Review the approved codebook and survey versions**: identify which waves contain age, sex, MS disease duration, **PDDS** or other disability measures, DMT history, relapses, functional scores, comorbidities, and relevant healthcare utilization/outcomes, if any. Fields vary by study/wave: **do not assume they exist**.
3. **Choose the clinical comparison**: e.g., new users of DMT A versus new users of DMT B. If the study instead addresses a medical device, verify that actual **device exposure** and appropriate outcome measurements were captured; otherwise do **not** claim device safety or device effectiveness.
4. **Build longitudinal timelines** from survey waves and reported DMT start/stop dates; derive a defensible treatment index. Distinguish prospective reporting from retrospective recall and establish sufficient pre-index observation to define new-user status.
5. **Select baseline confounders using only pre-index information**: age, sex, MS duration and severity/PDDS, prior DMTs, relapses, disability, and comorbidities when available. The demonstration's numeric `baseline_severity` is *not* an official PDDS field; ordinal PDDS modeling may need a different propensity-score specification.
6. **Establish outcome ascertainment windows**: registry survey dates are not necessarily precise adverse-event dates. An outcome only known to have occurred **between two surveys** is interval-censored; do **not** fabricate a date just to populate `safety_time_days`. Omit survival columns until a defensible exact-time phenotype or an appropriate interval-censoring method exists.
7. **Assess observation and missingness**: account for intermittent responses, dropout, missed biannual surveys, retrospective recall and treatment switching. A complete-case 365-day analysis may be selected; formally test whether its bias assumptions are credible.
8. **Do not assume hospital/ED counts, costs or utility-based QALYs are captured**. Populate these columns only with valid source or authorized linked data; otherwise **omit them**, rather than inserting zeros or guessed values. Running with `--cost-basis none` prevents economic outputs.
9. **Derive a one-row-per-person cohort** with unique registry research ID, index date, pretreatment variables, treatment group, observed follow-up, and only clinically validated endpoints. Produce cohort-flow counts and check survey-wave consistency against the approved data dictionary.
10. **Analyze within authorized workspace; review outputs** against the DUA, IRB and NARCOMS manuscript-review requirements.

### 5.2 Typical NARCOMS logical source domains (not literal field names)

| Normalized field | Potential derivation and special limitation |
|---|---|
| `patient_id` | Stable approved, de-identified registry participant key |
| `index_date`, `treated` | Validated initiation/switch of a captured DMT or **verified device exposure**; index may be approximate |
| `age`, `sex` | Authorized enrollment/update demographic fields |
| `baseline_severity` | Pretreatment PDDS or other disability measurement **with explicitly documented transformation**; ordinal/nonlinear modeling may be needed |
| `comorbidity_index` | Documented score from reliably measured pretreatment comorbidities; don't claim Elixhauser if ICD diagnoses unavailable |
| `prior_hospitalizations_365d` | Only when prior-year hospitalizations are explicitly captured; otherwise protocol/method changes required |
| `followup_days` | Time from treatment index to last verified survey/observation under defined censoring scheme |
| `safety_event`, `safety_time_days` | Only with sufficiently precise, validated event status and onset/censor dates; **not guaranteed** |
| `clinical_response_365d` | Prespecified patient-reported improvement, relapse-free status, or disability progression status, with appropriate longitudinal assessment |
| HCRU columns | Only if the relevant resource-use instruments and time windows support them |
| Cost/QALY columns | Only via validated registry measurements or authorized economic/utility linkage |

**Never concatenate** the Cosmos and NARCOMS patient-level extracts unless linkage rights, harmonization, comparable targets, and scientific/statistical methods have been separately established. Cross-source comparison is not automatically a pooled causal analysis.

## 6. Required standardized patient-level input schema and column labels

**`demo_data/EpicCosmos_cohort.csv` and `demo_data/NARCOMS_cohort.csv` must follow the table below.** These columns are **application-defined research variables**, not verified native column names in either upstream source. Each row represents one eligible patient at a clinically defined, common time zero. Several columns require joining multiple raw encounters or survey waves and applying validated medical definitions.

| CSV column | Human-readable label | Expected data type / coding | Requirement | Meaning / required derivation |
|---|---|---|---|---|
| `patient_id` | Anonymous patient ID | Unique nonempty text | **Required** | Stable research ID within one source; no cross-source linkability assumed; no direct identifiers. |
| `index_date` | Treatment initiation (index) date | ISO date `YYYY-MM-DD` | **Required** | Clinically validated new-user treatment start for both arms. Not a general encounter/survey date. |
| `treated` | Assigned observational treatment cohort | Integer: `1` = intervention, `0` = active comparator | **Required** | Derived from verified treatment/device exposure; both arms must contain eligible patients. Not randomized assignment. |
| `age` | Age at index | Numeric years | **Required** | Pretreatment age at index; obey masking/top-coding in approved source. |
| `sex` | Recorded sex | Categorical text such as `F`, `M`, or valid source code | **Required column** | Source-documented baseline variable; harmonize coding and record missingness. |
| `baseline_severity` | Pretreatment disease-severity measure | Numeric, with declared scale/direction | **Required** | Standardized or specified clinical measure; **do not assume a native Cosmos field or equate a NARCOMS PDDS score to a continuous biomarker without justification**. |
| `comorbidity_index` | Pretreatment comorbidity-burden score | Numeric | **Required** | Validated pre-index score available in this source; document components and lookback period. |
| `prior_hospitalizations_365d` | Hospital admissions in prior year | Nonnegative integer count | **Required** | Pre-index admission count with sufficient observation history; do not replace missing registry utilization with zero. |
| `followup_days` | Observable follow-up duration | Positive numeric days | **Required** | Time since index to last reliable encounter/survey/censoring or study end; must account for source-specific observation limits. |
| `safety_event` | Prespecified safety outcome observed | `1` = event, `0` = observed no event, empty = unknown | **Optional paired** | Clinical event flag, using validated outcome definition; use together with `safety_time_days`. |
| `safety_time_days` | Days to safety event or censoring | Numeric days, >0 and <= `followup_days` | **Optional paired** | Index-to-event/censor time; requires dependable timing and censoring data. Do not assume NARCOMS surveys contain adjudicated event dates. |
| `clinical_response_365d` | Binary response at day 365 | `1` = responder, `0` = nonresponder, empty = unobserved | **Optional** | Define a clinically meaningful, validated response endpoint; available only with evaluable one-year follow-up. |
| `hospitalizations_365d` | Hospital admissions after index | Nonnegative count or empty | **Optional** | Number of admissions in complete post-index year; not equal to pre-index admissions. |
| `ed_visits_365d` | Emergency department visits after index | Nonnegative count or empty | **Optional** | Number of ED visits over complete 365-day post-index observation. |
| `hospital_days_365d` | Hospital inpatient days after index | Nonnegative count/days or empty | **Optional** | Total inpatient days over complete 365-day follow-up. |
| `cost_365d_usd` | Validated one-year healthcare cost (USD) | Nonnegative numeric USD or empty | **Optional** | Genuine patient-level costs with cost-year/perspective documentation; charges alone do not qualify. Only used with `--cost-basis validated_costs`. |
| `qaly_365d` | Quality-adjusted life-years during 365 days | Numeric QALYs (typically 0 to 1) or empty | **Optional** | Derived from valid utility measurements over time; neither source is assumed to provide QALYs directly. |

**Essential versus optional:** `patient_id`, `index_date`, `treated`, `age`, `sex`, `baseline_severity`, `comorbidity_index`, `prior_hospitalizations_365d`, and `followup_days` must all be present. `safety_event` and `safety_time_days` must be included **together or omitted together**. Other outcomes are included only if actually measured. The current code cannot analyze a different baseline covariate set merely by dropping a required column; it would need a study-specific revision.

**Source-specific readiness:** Epic Cosmos clinical data generally require joins across diagnoses, procedures, medications, encounters and observations. NARCOMS generally requires linking enrollment, repeated questionnaires, treatment histories and disability assessments. Neither is guaranteed to contain complete all-cause utilization, validated safety event dates, device identifiers, costs or QALYs. **A synthetic placeholder demonstrates file compatibility, not scientific or source-schema equivalence.**

**Follow-up rule:** The existing code uses a 365-day fixed horizon. Missing responses, utilization and cost outcomes must remain missing when less than 365 days are observed. Do not use apparent zeros to represent unobserved outcomes. Safety uses right-censoring, but the analyst must validate its assumptions. Document clinical definitions, code lists, follow-up ascertainment, missingness, site effects, source terminology and date conventions in the study protocol.

### 6.1 Example field mapping (illustrative placeholders, not source SQL names)

Once **the actual approved data dictionary** is available, edit a source-specific adapter. For example:

```python
import pandas as pd
from pathlib import Path

# Run only inside the workspace permitted for this source's records.
# The capitalized names below are PLACEHOLDERS, not Epic/NARCOMS fields.
raw = pd.read_csv("authorized_source_specific_staging.csv")
column_map = {
    "APPROVED_RESEARCH_PERSON_KEY": "patient_id",
    "VERIFIED_TREATMENT_START": "index_date",
    "NEW_TREATMENT_INDICATOR": "treated",
    "AGE_AT_START": "age",
    "REPORTED_SEX": "sex",
    "PRETREATMENT_SEVERITY": "baseline_severity",
    "PRETREATMENT_COMORBIDITY_SCORE": "comorbidity_index",
    "PRETREATMENT_ADMISSIONS_365D": "prior_hospitalizations_365d",
    "OBSERVED_FOLLOWUP_DAYS": "followup_days",
}
missing = sorted(set(column_map) - set(raw.columns))
if missing:
    raise ValueError(f"Approved mapping incomplete: {missing}")
cohort = raw.rename(columns=column_map)
# WARNING: independently validate source-specific exposure, event,
# response, HCRU and economic derivations BEFORE populating those columns.
assert cohort["patient_id"].is_unique
Path("approved_analysis_dir").mkdir(exist_ok=True)
cohort.to_csv("approved_analysis_dir/normalized_cohort.csv", index=False)
```

This adapter is **not** an extraction method and not sufficient for a causal study without the upstream cohort-building logic, source validation, authorized location, endpoint definitions and medical review. In particular, preserve source lineage: mapping rules, SQL/query revision, approved source version, clinical coding vocabulary/version and access controls.

## 7. Exactly what statistical analysis the code performs

### 7.1 Cohort construction and data checks

The main script **does not construct the real clinical cohort**. It accepts a prespecified normalized cohort and checks required fields, unique patient IDs, valid dates and groups, positive follow-up, binary event indicators, plausible ages/nonnegative outcomes, and consistency of event/censor times. It masks fixed-year outcomes for insufficient follow-up. The actual population definition, coding and clinical validity must be verified upstream.

### 7.2 Confounding adjustment and the chosen estimand

- Fits regularized logistic regression for `treated` using age, baseline severity, comorbidity index, prior hospitalizations and sex. Numeric covariates are median-imputed/standardized; sex is most-frequent-imputed and one-hot encoded. This is a **demonstration PS specification**, not an approved confounder set for any real clinical question.
- Estimates `e(X) = P(treated=1 | pre-index X)`, clipping predictions at `1e-5` and `1-1e-5` to avoid literal division by zero; clipping is reported.
- Applies normalized weighting to estimate treated-versus-control contrasts:

| Estimand | Target | Treated weight | Comparator weight |
|---|---|---|---|
| `ATE` | All eligible patients | `1/e(X)` | `1/[1-e(X)]` |
| `ATT` | Patients actually treated | `1` | `e(X)/[1-e(X)]` |
| `ATO` (default) | Clinical equipoise/overlap population | `1-e(X)` | `e(X)` |

- Calculates the **Bhattacharyya** propensity-score overlap coefficient `phi` (closer to 1 indicates greater overlap), weighted Kish effective sample sizes (ESS) separately in treated/comparator arms, PS minima/percentiles, and signed/absolute standardized mean differences (SMD) before and after weighting. Flags residual maximum `|SMD| > 0.10` and very low ESS.
- **ATC is not implemented in this version**, even though it is defined in some theoretical material. Add it deliberately with validation if required.

The weighting is not a guarantee of unbiased causal inference; it requires consistency, no relevant unmeasured confounding conditional on measured baseline covariates, positivity, valid time zero and appropriate treatment/outcome models.

### 7.3 Safety: weighted time-to-event risk

If **both** `safety_event` and `safety_time_days` exist and each arm has at least 10 analyzable observations, the code estimates **weighted Kaplan–Meier cumulative event risk** at day 365 in each arm, its **risk difference** (treated minus control) and **risk ratio** (treated / control, when comparator risk >0). These are **NOT hazard ratios** and not a competing-risk cumulative incidence function. KM assumes appropriately independent censoring; competing death may require Aalen–Johansen, and informative dropout may require censoring adjustment.

### 7.4 Clinical effectiveness: binary response

For records with complete 365-day follow-up and nonmissing `clinical_response_365d`, the application estimates weighted response probability in each arm, **absolute risk difference** and **risk ratio**. Positive risk difference means a higher response probability with treatment. Response definitions must be prespecified; in NARCOMS, a change in PDDS may require a different model (e.g. ordinal or longitudinal) rather than naive binary coding.

### 7.5 HCRU: admissions, ED visits and hospital days

For each available 365-day resource-use column, estimates weighted mean per person in each arm and the **treated-minus-control mean difference**. A value of `-0.24` hospital admissions means **0.24 fewer admissions per patient over the full one-year horizon** in the weighted comparison (not an incidence-rate ratio, not 0.24 per 100 patient-years). The present code uses mean contrasts, not count-regression, rate offsets, cost-of-illness models or repeated-measures methods.

### 7.6 Economic value: conditional module

Executed **only** when `--cost-basis validated_costs` and a `cost_365d_usd` column are supplied. In real data, documentation must confirm cost validity, perspective, inflation adjustment, period and completeness; synthetic values have no real validity.

- Estimates weighted **mean cost difference** for complete 365-day cost observations.
- For joint cost-effectiveness analysis, uses the **same complete-case patient subset** for both cost and benefit. Prefers `qaly_365d` if available, otherwise falls back to `clinical_response_365d` (which yields cost per additional responder, **not cost per QALY**).
- Computes incremental cost-effectiveness ratio `ICER = (mean cost_T - mean cost_C)/(mean benefit_T - mean benefit_C)` when the benefit difference is nonzero. No ICER confidence interval is reported because ratios can be unstable. **If costs decrease and benefits increase, the treated strategy is descriptively dominant; a negative ICER should not be reported as a meaningful negative cost per QALY.**
- Computes incremental net monetary benefit `INMB = WTP * delta_benefit - delta_cost` at `50,000 / 100,000 / 150,000 USD/QALY` by default **only for a QALY analysis**. For responder-based benefits use explicit `--responder-wtp 10000 20000` with a justified WTP, and never label such values USD/QALY.
- Formal economic decision analysis still requires a defined perspective, uncertainty plane, cost-effectiveness acceptability curve, and sensitivity analyses.

### 7.7 Uncertainty intervals and limitations

The application performs a **stratified patient-level bootstrap** by treatment arm, **re-estimating propensity scores** and weights within every resample; 95% percentile intervals are the 2.5th and 97.5th bootstrap percentiles where a sufficient fraction of resamples succeed. No Wald p-values or false claims of statistical significance are generated. **Intervals crossing 0** (differences) or **1** (ratios) are not strong evidence of a difference at the 95% interval level, but they **do not prove equivalence**.

The code currently lacks cluster-aware bootstrap for multisite Cosmos, inverse probability of censoring weights for incomplete fixed-horizon outcomes, competing-risk models, non-proportional-hazards diagnostics, negative controls, time-varying DMT analyses and comprehensive small-cell suppression. These extensions may be essential before submitting results to regulators, payers or peer-reviewed journals.

## 8. Output files: meaning and Word-report interpretation

One `analyze` run produces the following files **together in the selected output directory**:

| Output | What it contains | Interpretation / caution |
|---|---|---|
| `effect_estimates.csv` | Adjusted treated/control means or risks; differences/ratios; 95% bootstrap intervals; per-outcome arm Ns; bootstrap counts; INMB probability positive | Contrast is for specified **estimand** and eligible analysis subset; not automatically causal |
| `table1_baseline.csv` | Raw and propensity-weighted means/proportions; observed counts; missing counts; signed SMDs | Check baseline comparability before/after weighting |
| `covariate_balance.csv` | SMD before/after for each included baseline variable | `|SMD| < 0.10` is a diagnostic convention, not proof of no confounding |
| `cohort_summary.csv` | Raw arm counts, complete follow-up counts, baseline mean age, selected raw adverse-event totals, median follow-up | Counts are descriptive, not adjusted; small-cell suppression requires independent review |
| `missingness.csv` | Number/fraction missing for every analysis variable | Large or differential missingness can bias complete-case results |
| `diagnostics.json` | Design metadata, overlap, ESS, score range, balance, PS clipping, warnings and software versions | Source provenance is analyst supplied unless explicitly labeled synthetic; no clinical validation is performed |
| `propensity_overlap.png` | Treated versus comparator propensity score distributions | Look for common support and extreme treatment probabilities |
| `covariate_balance.png` | Absolute SMD before and after weights | Look for residual baseline imbalance |
| `evidence_summary.md` | Compact results table and core methodological cautions | Existing report remains unchanged in purpose |
| **`EpicCosmos_RWE_Evidence_Report.docx`** or **`NARCOMS_RWE_Evidence_Report.docx`** | **All output tables, design diagnostics, figures, source-specific limitation statements, and automatically generated descriptive interpretations** | **Automatically generated in the same directory; exploratory, not a validated medical manuscript** |

The Word document includes **all adjusted effect rows** reported in `evidence_summary.md`, along with the study population, overlap coefficient, weighted effective sample sizes (ESS), balance diagnostic, baseline table, group counts, observed sample sizes, missingness, full diagnostics, figures and narrative caveats. Its tables use **human-readable technical terms**, while `diagnostics.json` and the CSV outputs retain the exact machine-readable codes.

**Word report table conventions:**

- `95% CI` means **95% confidence interval (CI)**. The number of bootstrap resamples (`B`) is displayed **below each results table**, not combined into the confidence interval column. The `n_bootstrap` value is available in `effect_estimates.csv` for individual indicators.
- `Treated N` and `Control N` are **separate** sample-count columns; they can differ across outcomes because of analysis-specific missing data and follow-up.
- In the baseline table, raw control/treated means (`C raw`, `T raw`) are before weighting, and weighted means (`C wtd`, `T wtd`) are after propensity-score adjustment. The standardized mean difference (`SMD`) columns distinguish **before** and **after** weighting. Expanded terms appear in the Word report and its abbreviation glossary.
- All diagnostic field identifiers are accompanied by descriptive full names, and nested missingness indicators are included. Software runtime/dependency versions remain in `diagnostics.json` for audit purposes but **are not printed in the Word report**.
- The incremental cost-effectiveness ratio (`ICER`), quality-adjusted life year (`QALY`), incremental net monetary benefit (`INMB`) and willingness-to-pay threshold (`WTP`) are expanded in the report. **QALY**, not `CALY`, is the correct term.
- The prominent synthetic demonstration status banner at the top of Word reports has been removed for cleaner formatting. **Data origin remains disclosed** through the `SYNTHETIC_` input filename and the `Input data provenance (data_origin)` entry in the full diagnostic metadata. The examples are synthetic; removing the banner does not make them clinical evidence.
- **Sections 1–4 of the Word report contain only their tables, metrics, definitions and necessary method notes; their former `Interpretation` blocks have been removed.** This applies to the safety, effectiveness, healthcare utilization and economic result tables within Section 3. The missing-data `Interpretation` block in Section 6 and the later methodological cautions remain. Section headings use `1.`, `2.`, etc., rather than pipe characters.
- For economic outcomes, the bootstrap proportion of resamples with a positive incremental effect is reported in an additional narrative note, **not mixed into the confidence interval column**.

**Consistency validation:** `test_pipeline.py` contains a report-reconciliation test for both synthetic sources, asserting that every row in `evidence_summary.md` is reproduced in the Word document with the same 4-decimal estimates and confidence intervals. This checks output transcription, **not the scientific validity of the underlying analysis**.

### 8.1 Interpretation of the supplied 1,000-bootstrap SYNTHETIC reports

These examples are from the supplied Markdown demonstration outputs; they are **not effects observed in Cosmos or NARCOMS**.

- **Epic Cosmos-style:** 450 simulated participants; ATO propensity-score overlap `phi=0.933`, weighted ESS treated `210.5`, control `182.7`, largest weighted `|SMD|=0.005`. The adverse-event 365-day risk difference was `-0.0524` (95% bootstrap CI `[-0.1200, +0.0054]`): point estimate suggests 5.24 percentage points fewer events, **but interval includes no difference**. Hospital admissions difference `-0.2382` (CI `[-0.4566, -0.0261]`). Simulated cost difference `-$1,307` and simulated QALY difference `+0.0323` represent **descriptive dominance in the synthetic cost/utility model**, not demonstrated healthcare value.
- **NARCOMS-style:** 450 simulated participants; ATO `phi=0.910`, weighted ESS treated `171.0`, control `204.4`, largest weighted `|SMD|=0.007`. Adverse-event 365-day risk difference `+0.0098` (CI `[-0.0501, +0.0686]`); clinical response difference `+0.0405` (CI `[-0.0661, +0.1502]`). **Both intervals include no difference**. No economic module ran because this synthetic file lacks validated cost/QALY inputs.

Results may change when data are regenerated, covariates/outcomes are redefined or bootstrap settings change. Never cite these demonstration quantities as published clinical effects.

## 9. Commands: demonstration, current analysis and report-only

### 9.1 Regenerate data and run the demonstration

Run from the folder containing `rwe_evidence_pipeline.py`:

```bash
# Synthetic data: 450 records per source
python rwe_evidence_pipeline.py make-demo --output-dir ./demo_data --n 450 --seed 2026

# Epic-like synthetic analysis: economics enabled ONLY to test simulated data
python rwe_evidence_pipeline.py analyze \
  --input ./demo_data/EpicCosmos_cohort.csv \
  --source EpicCosmos \
  --output-dir ./results/cosmos \
  --estimand ATO \
  --cost-basis validated_costs \
  --bootstrap 1000

# NARCOMS-like synthetic analysis: no economics
python rwe_evidence_pipeline.py analyze \
  --input ./demo_data/NARCOMS_cohort.csv \
  --source NARCOMS \
  --output-dir ./results/narcoms \
  --estimand ATO \
  --cost-basis none \
  --bootstrap 1000
```

In each output directory you will find **both** `evidence_summary.md` **and** the new `.docx` report. Reduce bootstraps for a quick smoke test (`--bootstrap 100`); **1,000** is a practical exploratory setting, not a validation guarantee.

### 9.2 Build Word documents from analysis outputs already on disk

**No rerun or patient-level input needed**:

```bash
python rwe_evidence_pipeline.py report \
  --output-dir ./results/cosmos \
  --source EpicCosmos

python rwe_evidence_pipeline.py report \
  --output-dir ./results/narcoms \
  --source NARCOMS
```

The `report` command reads the required `.csv`, `diagnostics.json`, `evidence_summary.md` and available `.png` files from that directory. It **does not recalculate effects** or substitute different bootstrap estimates. It overwrites the existing Word report of the same name; the original CSV/JSON/Markdown outputs remain unchanged. If a required file is missing, it lists that file and asks you to run `analyze` first.

### 9.3 Run an approved REAL cohort (only in permitted environment)

```bash
python rwe_evidence_pipeline.py analyze \
  --input ./approved_analysis_dir/normalized_cohort.csv \
  --source EpicCosmos \
  --output-dir ./approved_results/cosmos \
  --estimand ATO \
  --cost-basis none \
  --bootstrap 1000
```

Use `--source NARCOMS` and an independently prepared NARCOMS cohort for registry analysis. For actual Cosmos data, **do not copy restricted line-level data to local or cloud storage** unless your institution and Epic explicitly permit it; use the approved virtual environment when required.

### 9.4 Additional CLI options and tests

- `--estimand ATE` or `--estimand ATT` (ATO by default); each estimand targets a different population.
- `--no-word-report` on `analyze` skips Word generation (useful if `python-docx` is intentionally unavailable).
- `--cost-basis none` prohibits economic outputs even when a cost column is present; for real cost analysis use `--cost-basis validated_costs` only after validation.
- `--seed 2026` fixes stochastic bootstrap initialization. Data generation has a separate `--seed` on `make-demo`.
- `--min-cell 11` is only the **minimal demonstration rule** for some raw event counts; NOT a comprehensive disclosure-control system.

Tests:

```bash
python -m pip install pytest
python -m pytest -q test_pipeline.py
```

## 10. Required validation before scientific, clinical or regulatory use

1. **Study design**: prespecified protocol, target trial, eligibility, intervention/comparator, estimand, time zero, follow-up, censoring, treatment switching, negative controls and sensitivity analyses.
2. **Clinical code lists**: ICD-9/10-CM, CPT/HCPCS, device model/UDI, RxNorm where relevant; documented vocabulary versions, data lineage, source drift and local validation samples.
3. **Positivity and balance**: propensity score overlap, extreme weights, ESS, nonlinearity, interactions, appropriate site and calendar-year covariates; consider diagnostics for continuous covariates beyond simple linear/logit models.
4. **Outcome definitions and missing data**: precise event dates, self-report limitations, competing risks, informative censoring/IPCW, interval censoring, hospital/network leakage and survey response dropout.
5. **Statistical precision**: pre-planned sample-size/power, adequate effective events, source/site-clustered inference where needed, bootstrap diagnostics and simulation reflecting the actual target estimand.
6. **Economic validity**: actual resource unit costs and verified cost perspective/price year; clinically valid benefits/QALYs; uncertainty, subgroup and cost-effectiveness acceptability analyses.
7. **Source restrictions**: Epic Cosmos rules of the road, institutional compute and export rules, NARCOMS research agreements/manuscript review, IRB/HIPAA and source-specific reporting/suppression requirements.
8. **Reproducibility**: pinned dependencies, source-specific extraction queries and dictionary, review log, cohort attrition diagram, stable seed, methods/limitations in a statistical analysis plan.

**Disclosure**: The program generates aggregate tables but **does not comprehensively suppress small cells, complementary cells or other identifying combinations**. Every output, Word report and figure must undergo the source's required disclosure review before being shared outside the approved research environment.

## 11. Primary references and official access guidance

- [Epic Cosmos: Request Access](https://cosmos.epic.com/request-access/)
- [Stanford Medicine: Cosmos Frequently Asked Questions (line-level export restrictions)](https://med.stanford.edu/starr-tools/other-resources/cosmos/faq.html)
- [UC Davis Health: Epic Cosmos data-science access and sharing rules](https://health.ucdavis.edu/data/epic-cosmos.html)
- [NARCOMS: Researchers and Providers / Data Access and Collection](https://www.narcoms.org/researchers-providers)
- [NARCOMS: Contact Information](https://www.narcoms.org/contact-us)
- [FDA: Use of Real-World Evidence to Support Regulatory Decision-Making for Medical Devices](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/use-real-world-evidence-support-regulatory-decision-making-medical-devices)
- [TARGET reporting guidance (BMJ, 2025)](https://www.bmj.com/content/390/bmj-2025-087179)
- [RECORD reporting statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC4595218/)
- Hernán MA, Robins JM (2016). *Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available.* American Journal of Epidemiology. DOI: 10.1093/aje/kwv254.


## 12. Consolidated version history (formerly CHANGES_v2.md)

### Version 3 — 2026-11-17

1. Reorganized the Word reporting component of `rwe_evidence_pipeline.py` without changing the statistical estimators or the existing result CSV/Markdown formats.
2. Verified that **all 14 adjusted effect indicators** in the Epic Cosmos-style summary and **all 7** in the NARCOMS-style summary are included in the respective Word reports, with exact rounded estimates and 95% confidence intervals. Reconciliation also covers study-population diagnostics, overlap, weighted sample size, and balance.
3. Split `95% CI / B` into `95% CI` in the effects tables; bootstrap resample counts (`B`) now appear **outside** the tables. Split the combined `N T/C` into distinct `Treated N` and `Control N` columns.
4. Expanded clinical, causal-inference, methodological and economic abbreviations. Replaced cryptic diagnostic field labels with long names plus the original machine-readable keys in parentheses. Added the previously nested outcome-specific diagnostic missingness fractions.
5. Removed the sentence listing Python/pandas/NumPy runtime versions from Word files; those versions remain available in the JSON diagnostics metadata.
6. Changed numbered Word section headings to use a **period** rather than a vertical bar. Replaced repetitive `Interpretation -` prefixes with a single heading followed by concise findings for each table.
7. Refined economically dense tables to separate model results from bootstrap uncertainty summaries, without removing any output indicators. Preserved the strict **synthetic demonstration** disclaimer.
8. Expanded report assertions in `test_pipeline.py` to reconcile every reported adjusted effect and confidence interval from both Markdown reports to the generated Word documents.
9. Consolidated these notes and the complete version 2 history into this **single README.md**; the separate `CHANGES_v2.md` file is no longer distributed.

### Version 2 — 2026-11-17 (merged from the former change log)

1. Added automatic source-specific Word (`.docx`) reports to `analyze`. Reports are assembled from the resulting CSV tables, diagnostics JSON, Markdown summary and existing plots; saved as `<Source>_RWE_Evidence_Report.docx` in the chosen output directory.
2. Added `report` subcommand, allowing Word-generation from existing aggregated results **without loading patient-level data or repeating bootstrap calculations**.
3. Added `--no-word-report` flag to skip Word generation when desired; added `python-docx` to `requirements.txt`.
4. Extended Word reporting with every adjusted effect, baseline and cohort table, covariate balance, missingness, diagnostic metadata, two propensity-score figures, source-specific methodological cautions, bootstrap interpretation and controlled provenance labeling.
5. Expanded approved real-source cohort construction guidance for Epic Cosmos and NARCOMS. Explicitly documented that **both `SYNTHETIC_*` CSVs are generated by the application** and were not extracted from those restricted sources; never invent actual source table names or permitted export rights.
6. Reproduced the two reference 450-person synthetic cohorts with 1,000 stratified patient-bootstrap replicates, matching the supplied Markdown summaries at their displayed precision.
7. Introduced automated tests covering report-only generation and explicit failure on missing prerequisite result files.

### Important scope limitations (unchanged)

- Neither version extracts or downloads patient-level Epic Cosmos or NARCOMS records; approved source-specific extraction and data dictionary mapping are prerequisites for applying it to genuine data.
- Statistical methods and clinical/cost inputs must undergo study-specific validation before peer-reviewed, regulatory or payer-facing claims can be supported. The demo outputs are synthetic and are not evidence of actual safety, effectiveness or cost savings.


## v3.1 — Word report revisions (2026-11-17)

Removed the `Interpretation` blocks from sections 1–4 of the Word reports for both sources, while preserving all tables and later cautions.

## v3.2 — Report title cleanup (2026-11-17)

Removed the standalone, all-caps synthetic demonstration status line from the top of generated Word reports (and the corresponding analyst-supplied unverified status line). The dataset filename and full diagnostic metadata still disclose data origin; methodological limitations remain. No statistical estimators, input schemas, results, or other report sections changed. To refresh an existing report without rerunning analysis:

```bash
python rwe_evidence_pipeline.py report --output-dir ./results/cosmos --source EpicCosmos
python rwe_evidence_pipeline.py report --output-dir ./results/narcoms --source NARCOMS
```

## v3.3 — Standardized input placeholders and clearer source-schema expectations (2026-11-17)

1. Added replaceable `demo_data/EpicCosmos_cohort.csv` and `demo_data/NARCOMS_cohort.csv` example inputs; both are initially **synthetic placeholders**, not source extracts or native-schema reproductions.
2. Updated `make-demo` to refresh the clearly labeled `SYNTHETIC_*` files, but **never overwrite existing** canonical study input paths. Canonical placeholders are created automatically only if absent.
3. Rewrote Sections 1, 3, and 6 with expected inputs, required/optional columns, human-readable labels, coding, source-specific mapping expectations and explicit caveats for unavailable clinical variables.
4. The computational estimators, Word report generator, and previously calculated example results remain unchanged.


## Appendix — Pilot-based sample-size planning (new in v4.0)

**Example: a hypothetical prespecified 0.100 absolute difference in 365-day clinical response**, 5% two-sided type-I error, 80% target power, no additional loss beyond pilot. The value 0.100 is *illustrative only*, not an established clinically meaningful difference. Run `analyze` first to generate the bootstrap standard error:

```bash
python rwe_evidence_pipeline.py analyze \
  --input ./demo_data/SYNTHETIC_EpicCosmos_cohort.csv \
  --source EpicCosmos --output-dir ./outputs/EpicCosmos \
  --estimand ATO --cost-basis validated_costs --bootstrap 1000

python rwe_evidence_pipeline.py plan-sample-size \
  --output-dir ./outputs/EpicCosmos \
  --metric-code clinical_response_365d_difference \
  --target-effect 0.10 --alpha 0.05 --power 0.80 \
  --test two-sided --extra-loss 0.00 --design-multiplier 1.0
```

Repeat with NARCOMS, using `--source NARCOMS` and `--cost-basis none`.

The planner saves `sample_size_plan.csv`, `sample_size_plan.md` and updates the Word report with a clearly labeled **illustrative planning** section. It uses the refitted-propensity-score bootstrap **standard error**, `SE`, to estimate an asymptotic variance coefficient `V = N_pilot * SE^2`. The approximate required enrolled sample is:

`N ≈ ceil[V * (z_(1-alpha/2) + z_(power))^2 / target_effect^2 * design_multiplier / (1 - extra_loss)]`.

The pilot's **existing outcome missingness / incomplete follow-up are already reflected in the empirical standard error**. The `--extra-loss` flag must represent only losses **beyond those seen in the pilot**, or sample size will be inflated twice. Do not use a treatment effect estimated from the same pilot to define the clinically important effect. Avoid post-hoc/achieved power claims after observing outcomes. Ratios and economic endpoints are intentionally ineligible for this additive-difference planning shortcut.

**What is not modeled:** changes in baseline case mix, PS overlap, unmeasured confounding, future sites, informative censoring or cluster correlation are not estimated automatically. Use scenario multipliers only when justified and cross-check with [CRAN PSpower v2.0](https://cran.r-project.org/web/packages/PSpower/index.html) for its specific estimators and design assumptions. An actual RWE protocol needs prespecified estimand, endpoint, target effect, allocation, follow-up, source validation and sensitivity analysis before sample size can be finalized.

### Re-generating the Word report after a plan

```bash
python rwe_evidence_pipeline.py report --output-dir ./outputs/EpicCosmos --source EpicCosmos
python rwe_evidence_pipeline.py report --output-dir ./outputs/NARCOMS --source NARCOMS
```

### Changes from v3.3 to v4.0

- Updated scientific name, with distinct cohort-analysis and pilot-sample-size modules.
- Unique outcome-specific measure codes and unambiguous terminal columns.
- Economic INMB treated/control arm-level NMBs; ICER clarified as ratio of increments.
- Three-decimal display in human-facing outputs; full precision preserved in raw data tables.
- Explanations for overlap coefficient and SMD weights/0.100 chart threshold in Word reports.
- New optional `plan-sample-size` command; regression tests and synthetic worked examples.

### Publishing the application to GitHub

The included `.gitignore` deliberately excludes short-name canonical input CSVs (`EpicCosmos_cohort.csv`, `NARCOMS_cohort.csv`), generated `outputs/`, and arbitrary CSV files. Only explicitly named `SYNTHETIC_*` demonstration CSVs are whitelisted. The ZIP archive still includes the generated Word examples so you can inspect them locally; add them to a public repository only after verifying they are synthetic and disclosure-safe. `.gitignore` cannot protect files that were already tracked in Git. Do not publish identifiable, restricted or DUA-covered source data.


### Changes from v4.0 to v4.1

- Renamed the application **Sample Size and Economic Analysis for RWE Observational Studies**.
- Changed Word-report section 10 to **Prospective sample-size planning**. The estimator remains a *pilot-bootstrap normal approximation*, and the report recommends independently checking assumptions and sample-size results with the [Duke PSpower calculator](https://www2.stat.duke.edu/~fl35/pspower.html). A match is not guaranteed because the two methods are not identical.
- Consolidated full citations and reporting guidance into a numbered **References** section at the end of each Word report, leaving only short cross-references in the body. Section 9 is now **Abbreviations**. The references appear as section 11 when the sample-size planning section is present, or section 10 otherwise.
- No changes to the underlying causal estimators, economic calculations, bootstrap outputs, or sample-size formula.
