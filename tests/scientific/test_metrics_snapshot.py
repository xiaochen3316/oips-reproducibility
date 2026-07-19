from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
from PIL import Image

from oips_repro import io, statistics
from oips_repro.cli import main
from oips_repro.figures import SOURCE_COLUMNS


ROOT = Path(__file__).parents[2]
EXPECTED = json.loads(
    (Path(__file__).parent / "data" / "expected_summary.json").read_text(encoding="utf-8")
)
ANALYSIS_FILES = {
    "final_top3_automated_QC.csv",
    "final_reference_mapping.csv",
    "final_md_cluster_v2_mapping.csv",
    "final_automated_evidence_labels.csv",
    "final_redocking_cluster_v2_mapping.csv",
    "unresolved_cases.csv",
    "final_cluster_v2_master_table.csv",
    "target_level_candidate_prioritization.csv",
    "final_candidate_prioritization_metrics.csv",
    "final_category_metrics.csv",
    "final_bootstrap_intervals.csv",
    "final_family_sensitivity.csv",
    "final_orel_ablation_targets.csv",
    "final_orel_ablation_categories.csv",
    "orel_ablation_cluster_rankings.csv",
    "representative_case_results.csv",
    "weight_sensitivity_scenarios.csv",
    "weight_sensitivity_targets.csv",
    "single_tool_target_ranks.csv",
    "single_tool_complete_case_metrics.csv",
}
STEMS = {
    "repository_summary_figure_1_candidate_landscape",
    "repository_summary_figure_2_qc_and_orel_ablation",
    "repository_summary_figure_3_posthoc_evidence",
    "repository_summary_figure_4_representative_cases",
}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_metrics_analysis_and_figure_bundle(tmp_path: Path):
    config = ROOT / "config" / "manuscript.yaml"
    cluster_dir = tmp_path / "cluster"
    static_dir = tmp_path / "static"
    analysis_dir = tmp_path / "analysis"
    figure_dir = tmp_path / "figures"
    assert main(["cluster", "--config", str(config), "--output", str(cluster_dir)]) == 0
    assert main([
        "score", "--config", str(config), "--cluster-dir", str(cluster_dir),
        "--output", str(static_dir),
    ]) == 0
    analyze_args = [
        "analyze", "--config", str(config), "--cluster-dir", str(cluster_dir),
        "--static-dir", str(static_dir),
        "--posthoc-data", str(ROOT / "data" / "posthoc"), "--output", str(analysis_dir),
    ]
    assert main(analyze_args) == 0
    assert {path.name for path in analysis_dir.iterdir()} == ANALYSIS_FILES

    metrics = pd.read_csv(analysis_dir / "final_candidate_prioritization_metrics.csv").iloc[0]
    for name, spec in EXPECTED["analysis_metrics"].items():
        assert metrics[name] == pytest.approx(spec["value"], abs=spec["abs_tolerance"])

    bootstrap = pd.read_csv(analysis_dir / "final_bootstrap_intervals.csv")
    loaded_config = io.load_manuscript_config(config)
    assert bootstrap.shape == (2 * len(statistics.BOOTSTRAP_METRICS), len(statistics.BOOTSTRAP_COLUMNS))
    assert bootstrap["iterations"].eq(loaded_config.data["statistics"]["bootstrap_iterations"]).all()
    family = bootstrap.loc[bootstrap["bootstrap_method"].eq("family_clustered")].set_index("metric")
    for key, spec in EXPECTED["bootstrap_family_intervals"].items():
        metric, column = key.split("__", 1)
        column = column.replace("CI_2_5_percent", "CI_2.5_percent").replace(
            "CI_97_5_percent", "CI_97.5_percent",
        )
        assert family.loc[metric, column] == pytest.approx(
            spec["value"], abs=spec["abs_tolerance"],
        )
    analysis_counts = EXPECTED["analysis_counts"]
    assert len(pd.read_csv(analysis_dir / "final_family_sensitivity.csv")) == analysis_counts["family_sensitivity_rows"]

    ablation = pd.read_csv(analysis_dir / "final_orel_ablation_targets.csv")
    categories = pd.read_csv(analysis_dir / "final_orel_ablation_categories.csv")
    cluster_ablation = pd.read_csv(analysis_dir / "orel_ablation_cluster_rankings.csv")
    assert (len(ablation), len(categories), len(cluster_ablation)) == (
        analysis_counts["ablation_target_rows"],
        analysis_counts["ablation_category_rows"],
        analysis_counts["ablation_cluster_rows"],
    )
    ablation_actual = {
        "Without_O_rel_Reference_Top1": float(ablation["Without_O_rel_reference_rank"].le(1).mean()),
        "Without_O_rel_Reference_Top3": float(ablation["Without_O_rel_reference_rank"].le(3).mean()),
        "Without_O_rel_Reference_MRR": float((1 / ablation["Without_O_rel_reference_rank"]).mean()),
        "mean_rank_Spearman_rho": float(ablation["rank_Spearman_rho"].mean()),
        "top_cluster_identity_change_fraction": float(ablation["top_cluster_identity_changed"].mean()),
    }
    for name, spec in EXPECTED["orel_ablation"].items():
        assert ablation_actual[name] == pytest.approx(
            spec["value"], abs=spec["abs_tolerance"],
        )
    assert pd.read_csv(analysis_dir / "representative_case_results.csv")["pdb_id"].tolist() == EXPECTED["representative_case_order"]

    weight_scenarios = pd.read_csv(analysis_dir / "weight_sensitivity_scenarios.csv")
    weight_targets = pd.read_csv(analysis_dir / "weight_sensitivity_targets.csv")
    assert len(weight_scenarios) == 10
    assert (weight_scenarios["baseline_top1_retained_N"].min(), weight_scenarios["baseline_top1_retained_N"].max()) == (20, 21)
    assert (weight_scenarios["Reference_Top1_N"].min(), weight_scenarios["Reference_Top1_N"].max()) == (11, 12)
    assert (weight_scenarios["Reference_Top3_N"].min(), weight_scenarios["Reference_Top3_N"].max()) == (17, 19)
    assert weight_scenarios["First_supported_Top1_N"].eq(20).all()
    assert weight_scenarios["First_supported_Top3_N"].eq(21).all()
    changed = weight_targets.loc[weight_targets["top1_retention_count"].lt(10)]
    assert changed[["pdb_id", "top1_retention_count", "top1_changed_scenarios"]].to_dict("records") == [{
        "pdb_id": "7O2I", "top1_retention_count": 9,
        "top1_changed_scenarios": "P_lig_decrease_20pct",
    }]

    single_targets = pd.read_csv(analysis_dir / "single_tool_target_ranks.csv")
    complete_ids = sorted(single_targets.loc[
        single_targets["method"].eq("OIPS-P") & single_targets["complete_case"], "pdb_id"
    ].tolist())
    assert complete_ids == EXPECTED["single_tool_complete_case_targets"]
    single_metrics = pd.read_csv(analysis_dir / "single_tool_complete_case_metrics.csv").set_index("method")
    assert single_metrics.loc["DoGSite3", ["Top1_N", "Top3_N", "Top5_N"]].tolist() == [12, 12, 13]
    assert single_metrics.loc["OIPS-P", ["Top1_N", "Top3_N", "Top5_N"]].tolist() == [7, 12, 13]
    assert single_metrics.loc["OIPS-P", "MRR"] == pytest.approx(0.679761904761905, abs=1e-12)

    hashes = {name: _sha(analysis_dir / name) for name in ANALYSIS_FILES}
    assert main(analyze_args) == 0
    assert hashes == {name: _sha(analysis_dir / name) for name in ANALYSIS_FILES}

    assert main([
        "figures", "--config", str(config), "--analysis", str(analysis_dir),
        "--output", str(figure_dir),
    ]) == 0
    expected_figures = {
        filename
        for stem in STEMS
        for filename in (f"{stem}.svg", f"{stem}.png", f"{stem}_source_data.csv")
    }
    assert {path.name for path in figure_dir.iterdir()} == expected_figures
    for stem in STEMS:
        assert "<text" in (figure_dir / f"{stem}.svg").read_text(encoding="utf-8")
        source = pd.read_csv(figure_dir / f"{stem}_source_data.csv")
        assert source.columns.tolist() == list(SOURCE_COLUMNS)
        with Image.open(figure_dir / f"{stem}.png") as image:
            assert abs(image.info["dpi"][0] - 600) < 1
    target_path = analysis_dir / "target_level_candidate_prioritization.csv"
    original_target = target_path.read_bytes()
    target = pd.read_csv(target_path, dtype=str, keep_default_na=False)
    target.loc[0, "reference_top1"] = "not-a-boolean"
    target.to_csv(target_path, index=False, lineterminator="\n")
    assert main([
        "figures", "--config", str(config), "--analysis", str(analysis_dir),
        "--output", str(figure_dir),
    ]) == 1
    target_path.write_bytes(original_target)
    (figure_dir / "unexpected.txt").write_text("x", encoding="utf-8")
    assert main([
        "figures", "--config", str(config), "--analysis", str(analysis_dir),
        "--output", str(figure_dir),
    ]) == 1


def test_analyze_and_figures_reject_unexpected_managed_directory_entries(tmp_path: Path):
    config = ROOT / "config" / "manuscript.yaml"
    bad = tmp_path / "analysis"
    bad.mkdir()
    (bad / "unexpected.txt").write_text("x", encoding="utf-8")
    assert main([
        "analyze", "--config", str(config), "--cluster-dir", str(tmp_path / "missing-cluster"),
        "--static-dir", str(tmp_path / "missing"),
        "--posthoc-data", str(ROOT / "data" / "posthoc"), "--output", str(bad),
    ]) == 1
