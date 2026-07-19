from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

from oips_repro.statistics import (
    AblationResult,
    SingleToolResult,
    WeightSensitivityResult,
    bootstrap_intervals,
    build_target_summary,
    leave_one_family_out,
    orel_ablation,
    single_tool_prioritization,
    summarize_by_group,
    summarize_metrics,
    weight_sensitivity,
)


ROOT = Path(__file__).parents[2]
MODULE_WEIGHTS = {
    "C_cons": 0.22,
    "G_geo": 0.18,
    "P_lig": 0.24,
    "O_rel_formal": 0.24,
    "Q_evidence": 0.12,
}


def _small_static_inputs():
    rankings = pd.DataFrame(
        [
            ["AAAA", "AAAA_V2C001", 1.0, 90.0, 80.0, np.nan, 70.0, 60.0],
            ["AAAA", "AAAA_V2C002", 2.0, 80.0, 70.0, 50.0, 10.0, 40.0],
            ["BBBB", "BBBB_V2C001", 1.0, 60.0, 60.0, 60.0, 60.0, 60.0],
            ["CCCC", "CCCC_V2C001", 1.0, 50.0, 50.0, 50.0, 50.0, 50.0],
        ],
        columns=[
            "pdb_id", "cluster_v2_id", "Within_PDB_rank", "C_cons", "G_geo",
            "P_lig", "O_rel_formal", "Q_evidence",
        ],
    )
    labels = pd.DataFrame(
        [
            ["AAAA", "AAAA_V2C001", 1.0, "A_auto"],
            ["AAAA", "AAAA_V2C002", 2.0, "U_auto"],
            ["BBBB", "BBBB_V2C001", 1.0, "U_auto"],
            ["CCCC", "CCCC_V2C001", 1.0, "R_auto"],
        ],
        columns=["pdb_id", "cluster_v2_id", "Within_PDB_rank", "automated_evidence_label"],
    )
    systems = pd.DataFrame(
        [
            ["AAAA", "orthosteric", "family-a"],
            ["BBBB", "allosteric", "family-b"],
            ["CCCC", "orthosteric", "family-c"],
        ],
        columns=["pdb_id", "pocket_category", "protein_family"],
    )
    evaluable = pd.Series({"AAAA": True, "BBBB": False, "CCCC": True})
    return rankings, labels, systems, evaluable


def test_target_summary_uses_explicit_reference_denominator_and_zero_for_no_hit():
    rankings, labels, systems, evaluable = _small_static_inputs()

    target = build_target_summary(rankings, labels, evaluable, systems)
    by_id = target.set_index("pdb_id")

    assert bool(by_id.loc["AAAA", "reference_evaluable"])
    assert np.isnan(by_id.loc["AAAA", "reference_first_rank"])
    assert not bool(by_id.loc["AAAA", "reference_top1"])
    assert by_id.loc["AAAA", "reference_reciprocal_rank"] == 0.0
    assert not bool(by_id.loc["BBBB", "reference_evaluable"])
    assert np.isnan(by_id.loc["BBBB", "reference_top1"])
    assert np.isnan(by_id.loc["BBBB", "reference_reciprocal_rank"])
    assert by_id.loc["CCCC", "reference_rank_percentile"] == 1.0

    metrics = summarize_metrics(target)
    assert metrics["reference_evaluable_N"] == 2.0
    assert metrics["Reference_Top1"] == 0.5
    assert metrics["Reference_MRR"] == 0.5
    assert metrics["Reference_median_rank"] == 1.0


def test_metric_boolean_parser_treats_text_false_as_false_and_empty_denominators_as_na():
    frame = pd.DataFrame(
        {
            "reference_evaluable": ["False", "false"],
            "reference_top1": ["True", "False"],
            "reference_top3": ["True", "False"],
            "reference_top5": ["True", "False"],
            "reference_reciprocal_rank": [1.0, 0.0],
            "reference_first_rank": [1.0, np.nan],
            "reference_rank_percentile": [1.0, np.nan],
            "first_supported_evaluable": ["False", "False"],
            "first_supported_top1": ["True", "False"],
            "first_supported_top3": ["True", "False"],
            "first_supported_reciprocal_rank": [1.0, 0.0],
        }
    )

    metrics = summarize_metrics(frame)

    assert metrics["reference_evaluable_N"] == 0.0
    assert np.isnan(metrics["Reference_Top1"])
    assert metrics["first_supported_evaluable_N"] == 0.0
    assert np.isnan(metrics["First_supported_MRR"])


def test_bootstrap_is_repeatable_shuffle_invariant_and_honors_nondefault_config():
    rankings, labels, systems, evaluable = _small_static_inputs()
    target = build_target_summary(rankings, labels, evaluable, systems)

    first = bootstrap_intervals(target, iterations=31, seed=17)
    second = bootstrap_intervals(target.sample(frac=1, random_state=9), iterations=31, seed=17)

    pd.testing.assert_frame_equal(first, second)
    assert first.shape == (14, 7)
    assert first["iterations"].eq(31).all()
    assert first["random_seed"].eq(17).all()
    assert first["bootstrap_method"].drop_duplicates().tolist() == [
        "target_resampling", "family_clustered",
    ]
    assert leave_one_family_out(target)["excluded_family"].tolist() == [
        "family-a", "family-b", "family-c",
    ]
    grouped = summarize_by_group(target, "Pocket_category")
    assert grouped["group"].tolist() == ["allosteric", "orthosteric"]


def test_orel_ablation_is_missing_aware_and_uses_average_tie_ranks():
    rankings, labels, systems, evaluable = _small_static_inputs()
    tied = pd.DataFrame(
        [
            ["DDDD", "DDDD_V2C002", 1.0, 100.0, np.nan, np.nan, 0.0, np.nan],
            ["DDDD", "DDDD_V2C001", 2.0, 100.0, np.nan, np.nan, 100.0, np.nan],
        ],
        columns=rankings.columns,
    )
    rankings = pd.concat([rankings, tied], ignore_index=True)
    labels = pd.concat(
        [
            labels,
            pd.DataFrame(
                [
                    ["DDDD", "DDDD_V2C002", 1.0, "R_auto"],
                    ["DDDD", "DDDD_V2C001", 2.0, "U_auto"],
                ],
                columns=labels.columns,
            ),
        ],
        ignore_index=True,
    )
    systems = pd.concat(
        [systems, pd.DataFrame([["DDDD", "allosteric", "family-d"]], columns=systems.columns)],
        ignore_index=True,
    )
    evaluable.loc["DDDD"] = True

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = orel_ablation(
            rankings, labels, evaluable, systems, MODULE_WEIGHTS, tie_method="average"
        )

    assert isinstance(result, AblationResult)
    dddd = result.cluster_rankings.loc[result.cluster_rankings["pdb_id"].eq("DDDD")]
    assert dddd["Without_O_rel_score"].tolist() == pytest.approx([100.0, 100.0])
    assert dddd["Without_O_rel_rank"].tolist() == [1.5, 1.5]
    target = result.targets.set_index("pdb_id")
    assert target.loc["DDDD", "Without_O_rel_top_cluster_v2_id"] == "DDDD_V2C001"
    assert target.loc["DDDD", "Without_O_rel_reference_rank"] == 1.5
    assert result.targets.shape[1] == 17
    assert result.categories.shape[1] == 11
    assert result.cluster_rankings.shape[1] == 5
    allosteric = result.categories.set_index("Pocket_category").loc["allosteric"]
    assert allosteric["Without_O_rel_reference_MRR"] == pytest.approx(2 / 3)


def test_weight_sensitivity_builds_ten_scenarios_without_mutating_inputs():
    rankings, labels, _, _ = _small_static_inputs()
    before_rankings = rankings.copy(deep=True)
    before_labels = labels.copy(deep=True)

    result = weight_sensitivity(rankings, labels, MODULE_WEIGHTS)

    assert isinstance(result, WeightSensitivityResult)
    assert result.scenarios.shape == (10, 12)
    assert set(result.scenarios["perturbed_module"]) == set(MODULE_WEIGHTS)
    assert set(result.scenarios["direction"]) == {"decrease", "increase"}
    assert result.scenarios["target_N"].eq(3).all()
    assert result.targets["top1_retention_count"].between(0, 10).all()
    pd.testing.assert_frame_equal(rankings, before_rankings)
    pd.testing.assert_frame_equal(labels, before_labels)


def test_single_tool_comparison_uses_native_order_and_common_complete_cases():
    tools = ["CavityPlus", "DoGSiteScorer", "DoGSite3", "CASTpFold", "SiteMap"]
    feature_rows = []
    mapping_rows = []
    row_id = 1
    for pdb_id in ("AAAA", "BBBB"):
        for tool in tools:
            cluster = f"{pdb_id}_V2C002" if pdb_id == "AAAA" else f"{pdb_id}_V2C001"
            feature_rows.append([row_id, pdb_id, tool, 1, np.nan])
            mapping_rows.append([
                row_id, pdb_id, tool, f"{pdb_id}_{tool}_U001", "true", cluster,
                "mapped_to_cluster_v2",
            ])
            row_id += 1
    features = pd.DataFrame(
        feature_rows, columns=["row_id", "pdb_id", "tool", "display_order", "sitemap_rank"],
    )
    mapping = pd.DataFrame(
        mapping_rows,
        columns=[
            "row_id", "pdb_id", "tool", "same_tool_unit_id",
            "representative_for_tool_unit", "cluster_v2_id", "mapping_status",
        ],
    )
    rankings = pd.DataFrame(
        [
            ["AAAA", "AAAA_V2C001", 1.0], ["AAAA", "AAAA_V2C002", 2.0],
            ["BBBB", "BBBB_V2C001", 1.0], ["BBBB", "BBBB_V2C002", 2.0],
        ], columns=["pdb_id", "cluster_v2_id", "Within_PDB_rank"],
    )
    labels = pd.DataFrame(
        [
            ["AAAA", "AAAA_V2C001", "U_auto"], ["AAAA", "AAAA_V2C002", "R_auto"],
            ["BBBB", "BBBB_V2C001", "R_auto"], ["BBBB", "BBBB_V2C002", "U_auto"],
        ], columns=["pdb_id", "cluster_v2_id", "automated_evidence_label"],
    )

    result = single_tool_prioritization(features, mapping, rankings, labels)

    assert isinstance(result, SingleToolResult)
    assert result.target_ranks.shape == (12, 7)
    assert result.target_ranks["complete_case"].all()
    metrics = result.complete_case_metrics.set_index("method")
    assert metrics["N"].eq(2).all()
    assert metrics.loc["CavityPlus", "Top1"] == 1.0
    assert metrics.loc["OIPS-P", "Top1"] == 0.5
    assert metrics.loc["OIPS-P", "MRR"] == pytest.approx(0.75)


def test_importing_statistics_creates_no_files(tmp_path: Path):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [sys.executable, "-c", "import oips_repro.statistics"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert list(tmp_path.iterdir()) == []
