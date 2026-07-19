import pandas as pd

from oips_repro.reports import (
    build_analysis_report,
    build_numeric_crosscheck,
    build_rebuild_report,
    build_validation_report,
)
from oips_repro.validation import CheckResult, ValidationReport


def _manifest():
    return {
        "schema_version": 1,
        "configuration": {"path": "config/manuscript.yaml", "sha256": "a" * 64},
        "inputs": [{"path": "data/features.csv", "sha256": "b" * 64, "bytes": 10}],
        "outputs": [
            {"path": "analysis/final_candidate_prioritization_metrics.csv", "sha256": "c" * 64, "bytes": 20},
            {"path": "analysis/target_level_candidate_prioritization.csv", "sha256": "d" * 64, "bytes": 30},
        ],
        "managed_paths": ["analysis/metrics.csv", "reports/rebuild_report.md"],
        "summary": {
            "records": 16, "mapped": 15, "excluded": 1, "same_tool_units": 15, "clusters": 9,
            "targets": 3, "boundary_sensitive": 0, "maximum_diameter_A": 0.2,
            "maximum_formal_votes_per_cluster_tool": 1,
        },
    }


def _tables(value: float = 0.75):
    return {
        "final_candidate_prioritization_metrics.csv": pd.DataFrame(
            [{"analysis_level": "overall", "group": "all_targets", "Top1": value}]
        ),
        "target_level_candidate_prioritization.csv": pd.DataFrame(
            [{"pdb_id": "5J89", "Static_top1_reference_success": True}]
        ),
    }


def test_rebuild_report_cites_relative_inputs_outputs_and_hashes():
    text = build_rebuild_report(_manifest())

    assert "config/manuscript.yaml" in text
    assert "data/features.csv" in text and "b" * 64 in text
    assert "analysis/final_candidate_prioritization_metrics.csv" in text and "c" * 64 in text
    assert "records" in text and "16" in text and "maximum_diameter_A" in text
    assert not any(token in text for token in ("C:\\", "E:\\", "file://"))


def test_analysis_report_is_derived_from_current_tables():
    first = build_analysis_report(_tables(0.75), _manifest())
    second = build_analysis_report(_tables(0.50), _manifest())

    assert first != second
    assert "0.75" in first and "0.5" in second
    assert "5J89" in first
    assert "c" * 64 in first and "d" * 64 in first


def test_analysis_report_changes_when_a_nonpreview_row_changes():
    tables = _tables()
    tables["target_level_candidate_prioritization.csv"] = pd.DataFrame([
        {"pdb_id": f"T{index:03d}", "Static_top1_reference_success": True}
        for index in range(12)
    ])
    first = build_analysis_report(tables, _manifest())
    tables["target_level_candidate_prioritization.csv"].loc[10, "Static_top1_reference_success"] = False

    assert build_analysis_report(tables, _manifest()) != first


def test_numeric_crosscheck_uses_named_absolute_tolerances():
    expected = {
        "analysis_metrics": {
            "top1_reference_success": {"value": 0.8, "abs_tolerance": 0.05},
            "first_supported_mrr": {"value": 0.7, "abs_tolerance": 0.01},
        }
    }
    text = build_numeric_crosscheck(
        {"top1_reference_success": 0.78, "first_supported_mrr": 0.6}, expected
    )

    assert "top1_reference_success" in text and "PASS" in text
    assert "first_supported_mrr" in text and "FAIL" in text
    assert "abs_tolerance" in text and "0.02" in text


def test_validation_report_lists_every_check_and_count():
    report = ValidationReport((
        CheckResult("schema", "pass", "valid", {"files": 2}),
        CheckResult("optional_md", "warning", "not supplied", {}),
        CheckResult("checksum", "fail", "mismatch", {"files": 1}),
    ))
    text = build_validation_report(report)

    assert all(value in text for value in ("schema", "optional_md", "checksum"))
    assert "Pass: 1" in text and "Fail: 1" in text and "Warning: 1" in text
