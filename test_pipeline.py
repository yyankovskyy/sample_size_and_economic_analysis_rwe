"""Unit tests for the analytical template (synthetic only)."""
import numpy as np
import pytest
import rwe_evidence_pipeline as rw


def test_bhattacharyya_equals_one_when_ps_equals_allocation_probability():
    z = np.array([0] * 60 + [1] * 40)
    ps = np.full(len(z), 0.4)
    assert np.isclose(rw.bhattacharyya(ps, z), 1.0)


def test_weighted_mean_and_ess():
    assert np.isclose(rw.weighted_mean([0, 10], [3, 1]), 2.5)
    assert np.isclose(rw.effective_n(np.ones(10)), 10.0)


def test_km_risk_without_censoring_matches_weighted_event_fraction():
    df = rw.demo(n=120, seed=21, source="EpicCosmos")
    df["followup_days"] = 365
    df["safety_time_days"] = np.where(df["safety_event"] == 1, 111., 365.)
    ps, w, _ = rw.make_ps(df, "ATO")
    treated = df["treated"].eq(1).to_numpy()
    expected = rw.weighted_mean(df.loc[treated, "safety_event"], w[treated])
    assert np.isclose(rw.km_risk_at(df, w, arm=1, horizon=365), expected)


def test_economics_requires_explicit_validated_cost_basis():
    df = rw.demo(n=160, seed=29, source="EpicCosmos")
    _, w, _ = rw.make_ps(df, "ATO")
    _, without = rw.build_metrics(df, w, 365, "none", [100000])
    _, with_cost = rw.build_metrics(df, w, 365, "validated_costs", [100000])
    assert "cost_365d_difference" not in without
    assert "economic_joint_delta_cost" in with_cost


def test_duplicate_patient_id_is_rejected():
    df = rw.demo(n=110, seed=44, source="NARCOMS")
    df.loc[1, "patient_id"] = df.loc[0, "patient_id"]
    with pytest.raises(ValueError, match="unique"):
        rw.validate(df, horizon=365, cost_basis="none")


def test_partial_year_costs_are_excluded():
    df = rw.demo(n=120, seed=14, source="EpicCosmos")
    mask = df["followup_days"].lt(365)
    df.loc[mask, "cost_365d_usd"] = 1234.0
    validated, issues = rw.validate(df, horizon=365, cost_basis="validated_costs")
    assert validated.loc[mask, "cost_365d_usd"].isna().all()
    assert any("masked" in str(w) for w in issues)


def test_word_report_from_complete_existing_outputs(tmp_path):
    """Report-only mode must not need the underlying patient-level CSV."""
    from pathlib import Path
    import shutil
    from docx import Document
    folder = Path(__file__).parent / "outputs" / "NARCOMS"
    for fname in ["effect_estimates.csv", "table1_baseline.csv", "cohort_summary.csv",
                  "covariate_balance.csv", "missingness.csv", "diagnostics.json",
                  "evidence_summary.md", "propensity_overlap.png", "covariate_balance.png"]:
        shutil.copy2(folder / fname, tmp_path / fname)
    report = rw.generate_word_report(tmp_path, source="NARCOMS")
    assert report.is_file()
    loaded = Document(report)
    combined = " ".join(p.text for p in loaded.paragraphs)
    assert "SYNTHETIC DEMONSTRATION - NOT EVIDENCE FROM ACTUAL PATIENT RECORDS" not in combined
    assert "ANALYST-SUPPLIED SOURCE - PATIENT-LEVEL PROVENANCE NOT VERIFIED BY THIS SOFTWARE" not in combined
    # Dataset provenance is still explicit in the filename and diagnostic metadata.
    assert "SYNTHETIC_NARCOMS_cohort.csv" in combined
    assert "synthetic_demo" in " ".join(
        c.text for tbl in loaded.tables for row in tbl.rows for c in row.cells
    )
    assert "Baseline covariate balance" in combined
    assert "missingness" in combined.lower()
    assert len(loaded.tables) >= 6


def test_word_report_fails_loudly_for_missing_results(tmp_path):
    with pytest.raises(FileNotFoundError, match="effect_estimates.csv"):
        rw.generate_word_report(tmp_path, source="EpicCosmos")


def _extract_md_effect_rows(path):
    """Read all machine-readable estimates in the human Markdown table."""
    rows = []
    for line in path.read_text().splitlines():
        cells = [c.strip() for c in line.split("|")[1:-1]]
        if len(cells) == 7 and cells[0] in {"safety", "effectiveness", "HCRU", "economics"}:
            rows.append(dict(zip(["domain", "outcome", "measure", "treated", "control", "effect", "ci"], cells)))
    return rows


def _all_effect_table_rows(doc):
    rows = []
    for table in doc.tables:
        headings = [cell.text.strip() for cell in table.rows[0].cells]
        if headings[:2] == ["Outcome", "Effect measure"]:
            assert headings[5:] == ["95% CI", "Treated N", "Control N"]
            for row in table.rows[1:]:
                values = [cell.text.strip() for cell in row.cells]
                rows.append(values)
    return rows


@pytest.mark.parametrize("source,expected", [("EpicCosmos", 14), ("NARCOMS", 7)])
def test_report_preserves_all_markdown_effect_estimates_and_ci(source, expected):
    """100% reconciliation to both Markdown summaries, not just spot-checks."""
    from pathlib import Path
    from docx import Document
    root = Path(__file__).parent / "outputs" / source
    doc = Document(root / f"{source}_RWE_Evidence_Report.docx")
    md_rows = _extract_md_effect_rows(root / "evidence_summary.md")
    word_rows = _all_effect_table_rows(doc)
    assert len(md_rows) == len(word_rows) == expected

    for reported, produced in zip(md_rows, word_rows):
        assert rw._word_outcome_label(reported["outcome"]) == produced[0]
        assert rw._word_measure_label(reported["measure"]) == produced[1]
        for src_value, word_value in zip([reported["treated"], reported["control"], reported["effect"]],
                                         produced[2:5]):
            if src_value == "NA":
                assert word_value == "Not available"
            else:
                assert np.isclose(float(src_value.replace(",", "")), float(word_value.replace(",", "")), atol=0.001)
        assert reported["ci"].replace("not estimated", "Not estimated") == produced[5]
        assert "B=" not in produced[5]  # Bootstrap replication count is outside the CI table.

    # All summary population and diagnostic indicators appear in the Word report.
    import json
    diag = json.loads((root / "diagnostics.json").read_text())
    all_cells = " ".join(c.text for t in doc.tables for row in t.rows for c in row.cells)
    assert rw._word_int(diag["sample_n"]) in all_cells
    assert rw._word_number(diag["phi_bhattacharyya"], 3) in all_cells
    assert rw._word_number(diag["weighted_ess_treated"], 3) in all_cells
    assert rw._word_number(diag["weighted_ess_control"], 3) in all_cells
    assert rw._word_number(diag["max_abs_smd_after"], 3) in all_cells
    for key in diag:
        if key not in {"dependencies", "warnings_and_limitations", "missingness_fraction"}:
            assert "(" + key + ")" in all_cells
    for key in diag["missingness_fraction"]:
        assert "(missingness_fraction." + key + ")" in all_cells

    paras = " ".join(p.text for p in doc.paragraphs)
    assert "SYNTHETIC DEMONSTRATION - NOT EVIDENCE FROM ACTUAL PATIENT RECORDS" not in paras
    assert "ANALYST-SUPPLIED SOURCE - PATIENT-LEVEL PROVENANCE NOT VERIFIED BY THIS SOFTWARE" not in paras
    assert "Runtime versions" not in paras
    assert "python=3.13.5" not in all_cells
    assert "Sample Size and Economic Analysis for RWE Observational Studies" in paras
    assert "Bhattacharyya propensity-score overlap coefficient" in paras
    assert "0 (no overlap)" in paras
    assert "Unadjusted" in paras and "Weighted" in paras
    assert "0.100" in paras
    assert "1. Study population" in paras
    assert "8. Assumptions" in paras
    assert "1 | Study population" not in paras
    assert "Interpretation -" not in paras
    assert "Bootstrap resamples (B) for this table: 1,000" in paras
    if source == "EpicCosmos":
        assert "Incremental cost-effectiveness ratio (ICER)" in all_cells
    else:
        assert "Economic outcomes were not estimated" in paras


@pytest.mark.parametrize("source", ["EpicCosmos", "NARCOMS"])
def test_word_sections_one_to_four_have_no_interpretation_blocks(source):
    """Section-specific check: remove narrative blocks without removing later cautions."""
    from pathlib import Path
    from docx import Document
    root = Path(__file__).parent / "outputs" / source
    doc = Document(root / f"{source}_RWE_Evidence_Report.docx")
    in_first_four = False
    first_four_paragraphs = []
    later_paragraphs = []
    for par in doc.paragraphs:
        text = par.text.strip()
        if text.startswith("1. Study population and design diagnostics"):
            in_first_four = True
        elif text.startswith("5. Baseline covariate balance"):
            in_first_four = False
        if in_first_four:
            first_four_paragraphs.append(text)
        else:
            later_paragraphs.append(text)
    assert first_four_paragraphs
    assert not any(p == "Interpretation" or p.startswith("Interpretation -")
                   for p in first_four_paragraphs)
    assert "6. Data completeness and missingness (missingness.csv)" in later_paragraphs
    assert "Interpretation" in later_paragraphs
    assert "4. Baseline characteristics (table1_baseline.csv)" in first_four_paragraphs


def test_distributed_placeholders_match_documented_schema():
    """Both canonical files are populated synthetic compatible examples, not native exports."""
    from pathlib import Path
    import pandas as pd
    root = Path(__file__).parent / "demo_data"
    for src in ("EpicCosmos", "NARCOMS"):
        canonical = pd.read_csv(root / f"{src}_cohort.csv", dtype={"patient_id": "string"})
        marked = pd.read_csv(root / f"SYNTHETIC_{src}_cohort.csv", dtype={"patient_id": "string"})
        pd.testing.assert_frame_equal(canonical, marked)
        assert set(rw.REQUIRED).issubset(canonical.columns)
        assert canonical["patient_id"].str.startswith("SYNTHETIC_").all()
        assert canonical["patient_id"].is_unique


def test_make_demo_does_not_overwrite_replaced_canonical_paths(tmp_path, monkeypatch):
    """Explicit regression test for avoiding overwrite of user substituted files."""
    import subprocess
    import sys
    from pathlib import Path
    import pandas as pd
    project = Path(__file__).resolve().parent
    canonical = tmp_path / "EpicCosmos_cohort.csv"
    canonical.write_text("patient_id,important_data\nuser_data,keep_me\n")
    args = [sys.executable, str(project / "rwe_evidence_pipeline.py"), "make-demo",
            "--output-dir", str(tmp_path), "--n", "50", "--seed", "2026"]
    subprocess.run(args, check=True, capture_output=True, text=True)
    assert canonical.read_text() == "patient_id,important_data\nuser_data,keep_me\n"
    assert (tmp_path / "SYNTHETIC_EpicCosmos_cohort.csv").is_file()
    assert (tmp_path / "SYNTHETIC_NARCOMS_cohort.csv").is_file()
    assert (tmp_path / "NARCOMS_cohort.csv").is_file()
    assert (tmp_path / "NARCOMS_cohort.csv").read_text().startswith("patient_id,index_date,treated")


@pytest.mark.parametrize('source,expected', [('EpicCosmos', 14), ('NARCOMS', 7)])
def test_unique_measure_codes_and_explanations(source, expected):
    from pathlib import Path
    from docx import Document
    import pandas as pd
    root = Path(__file__).resolve().parent / 'outputs' / source
    effects = pd.read_csv(root / 'effect_estimates.csv')
    assert len(effects) == expected
    assert effects['measure'].is_unique
    assert effects['metric_code'].is_unique
    assert effects['statistic'].notna().all()
    assert 'bootstrap_se' in effects.columns
    assert source != 'NARCOMS' or (effects['domain'] != 'economics').all()
    doc = Document(root / f'{source}_RWE_Evidence_Report.docx')
    paragraphs = ' '.join(p.text for p in doc.paragraphs)
    assert '0 (no overlap)' in paragraphs
    assert '0.100' in paragraphs
    assert '10. Prospective sample-size planning' in paragraphs


def test_economic_icer_and_inmb_decomposition():
    from pathlib import Path
    import pandas as pd
    root = Path(__file__).resolve().parent / 'outputs' / 'EpicCosmos'
    estimates = pd.read_csv(root / 'effect_estimates.csv')
    def metric(code):
        return estimates.loc[estimates.metric_code.eq(code)].iloc[0]
    dc = metric('economic_joint_delta_cost')
    db = metric('economic_joint_delta_benefit')
    assert np.isclose(dc.effect, dc.treated - dc.control)
    assert np.isclose(db.effect, db.treated - db.control)
    icer = metric('economic_icer')
    assert np.isnan(icer.treated) and np.isnan(icer.control)
    assert np.isclose(icer.effect, dc.effect / db.effect, rtol=1e-10)
    assert np.isnan(icer.ci_low) and np.isnan(icer.ci_high)
    for wtp in (50000, 100000, 150000):
        nm = metric(f'inmb_{wtp}')
        assert np.isclose(nm.treated, wtp * db.treated - dc.treated)
        assert np.isclose(nm.control, wtp * db.control - dc.control)
        assert np.isclose(nm.effect, nm.treated - nm.control)
        assert np.isclose(nm.effect, wtp * db.effect - dc.effect)


def test_pilot_planning_scaling_and_invalid_metric(tmp_path):
    from pathlib import Path
    import shutil
    root = Path(__file__).resolve().parent / 'outputs' / 'EpicCosmos'
    for name in ('effect_estimates.csv', 'diagnostics.json'):
        shutil.copy2(root / name, tmp_path / name)
    n10 = rw.plan_sample_size_from_pilot(
        tmp_path, 'clinical_response_365d_difference', 0.10, refresh_word=False)
    n20 = rw.plan_sample_size_from_pilot(
        tmp_path, 'clinical_response_365d_difference', 0.20, refresh_word=False)
    assert n20['required_n_total_enrolled'] < n10['required_n_total_enrolled']
    assert n20['required_n_total_enrolled'] <= n10['required_n_total_enrolled'] / 4 + 2
    extra = rw.plan_sample_size_from_pilot(
        tmp_path, 'clinical_response_365d_difference', 0.10,
        additional_loss=0.15, refresh_word=False)
    assert extra['required_n_total_enrolled'] > n10['required_n_total_enrolled']
    with pytest.raises(ValueError, match='Only additive'):
        rw.plan_sample_size_from_pilot(tmp_path, 'clinical_response_365d_risk_ratio',
                                       0.10, refresh_word=False)


@pytest.mark.parametrize("source", ["EpicCosmos", "NARCOMS"])
def test_v41_word_sources_centralized(source):
    from pathlib import Path
    from docx import Document
    path = Path(__file__).parent / "outputs" / source / f"{source}_RWE_Evidence_Report.docx"
    doc = Document(path)
    ps = [p.text for p in doc.paragraphs]
    text = "\n".join(ps)
    assert "Sample Size and Economic Analysis for RWE Observational Studies" in text
    assert "9. Abbreviations" in ps
    assert "10. Prospective sample-size planning" in ps
    assert "11. References" in ps
    assert "Consider using the Duke PSpower analytical model" in text
    assert "https://www2.stat.duke.edu/~fl35/pspower.html" in text
    assert "It is not a replication of the Duke PSpower analytical model" not in text
    idx = ps.index("11. References")
    refs = "\n".join(ps[idx+1:])
    body = "\n".join(ps[:idx])
    for source_id in ("10.1214/26-AOS2649", "10.1002/sim.3697",
                      "10.1080/01621459.2016.1260466"):
        assert source_id in refs
        assert source_id not in body
    assert "Relevant reporting guidance" not in body
