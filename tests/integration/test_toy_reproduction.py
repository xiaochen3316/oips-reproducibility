from pathlib import Path
import json
import socket
import subprocess
import sys

import pandas as pd
import pytest

from oips_repro.cli import main


ROOT = Path(__file__).parents[2]
TOY = ROOT / "examples" / "toy"
CONFIG = TOY / "config" / "manuscript.yaml"
SNAPSHOT = TOY / "expected_summary.json"

CLUSTER_FILES = {
    "cluster_v2_candidates.csv", "cluster_v2_membership.csv",
    "tool_record_to_cluster_v2_mapping.csv", "excluded_unmappable_records.csv",
    "cluster_v2_boundary_audit.csv",
}
STATIC_FILES = {"cluster_v2_master_table.csv", "cluster_v2_static_rankings.csv"}
ANALYSIS_FILES = {
    "final_top3_automated_QC.csv", "final_reference_mapping.csv",
    "final_md_cluster_v2_mapping.csv", "final_automated_evidence_labels.csv",
    "final_redocking_cluster_v2_mapping.csv", "unresolved_cases.csv",
    "final_cluster_v2_master_table.csv", "target_level_candidate_prioritization.csv",
    "final_candidate_prioritization_metrics.csv", "final_category_metrics.csv",
    "final_bootstrap_intervals.csv", "final_family_sensitivity.csv",
    "final_orel_ablation_targets.csv", "final_orel_ablation_categories.csv",
    "orel_ablation_cluster_rankings.csv", "representative_case_results.csv",
    "weight_sensitivity_scenarios.csv", "weight_sensitivity_targets.csv",
    "single_tool_target_ranks.csv", "single_tool_complete_case_metrics.csv",
}
FIGURE_FILES = {
    f"repository_summary_figure_{number}_{stem}{suffix}"
    for number, stem in (
        (1, "candidate_landscape"), (2, "qc_and_orel_ablation"),
        (3, "posthoc_evidence"), (4, "representative_cases"),
    )
    for suffix in (".svg", ".png", "_source_data.csv")
}
REPORT_FILES = {
    "rebuild_report.md", "analysis_report.md", "numeric_crosscheck.md",
    "validation_report.md",
}


def _names(path: Path) -> set[str]:
    return {entry.name for entry in path.iterdir()}


def _bytes(path: Path) -> dict[str, bytes]:
    return {
        item.relative_to(path).as_posix(): item.read_bytes()
        for item in path.rglob("*") if item.is_file() and not item.is_symlink()
    }


def test_genuine_toy_runs_all_seven_commands_without_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys
):
    def blocked(*args, **kwargs):
        raise AssertionError("network access is forbidden in toy reproduction")

    monkeypatch.setattr(socket, "create_connection", blocked)
    monkeypatch.setattr(socket.socket, "connect", blocked)
    monkeypatch.setattr(socket.socket, "connect_ex", blocked)
    monkeypatch.setattr(socket, "getaddrinfo", blocked)
    monkeypatch.delitem(__import__("os").environ, "SCHRODINGER", raising=False)
    stage_bundle = tmp_path / "stage-bundle"
    cluster_dir, static_dir = stage_bundle / "clustering", stage_bundle / "static"
    analysis_dir, figure_dir = stage_bundle / "analysis", stage_bundle / "figures"

    assert main(["validate-data", "--config", str(CONFIG)]) == 0
    assert main(["cluster", "--config", str(CONFIG), "--output", str(cluster_dir)]) == 0
    assert main([
        "score", "--config", str(CONFIG), "--cluster-dir", str(cluster_dir),
        "--output", str(static_dir),
    ]) == 0
    assert main([
        "analyze", "--config", str(CONFIG), "--cluster-dir", str(cluster_dir),
        "--static-dir", str(static_dir),
        "--posthoc-data", str(TOY / "data" / "posthoc"), "--output", str(analysis_dir),
    ]) == 0
    assert main([
        "figures", "--config", str(CONFIG), "--analysis", str(analysis_dir),
        "--output", str(figure_dir),
    ]) == 0
    assert (_names(cluster_dir), _names(static_dir), _names(analysis_dir), _names(figure_dir)) == (
        CLUSTER_FILES, STATIC_FILES, ANALYSIS_FILES, FIGURE_FILES,
    )
    assert main([
        "verify", "--config", str(CONFIG), "--bundle", str(stage_bundle),
        "--snapshot", str(SNAPSHOT),
    ]) == 0
    assert _names(stage_bundle / "reports") == REPORT_FILES
    assert not (stage_bundle / "run_manifest.json").exists()

    bundle = tmp_path / "bundle"
    assert main([
        "reproduce", "--config", str(CONFIG), "--output", str(bundle),
    ]) == 0
    reports_before = _bytes(bundle / "reports")
    assert main([
        "verify", "--config", str(CONFIG), "--bundle", str(bundle),
        "--snapshot", str(SNAPSHOT),
    ]) == 0
    assert _bytes(bundle / "reports") == reports_before

    assert _names(bundle) == {
        "clustering", "static", "analysis", "figures", "reports", "run_manifest.json"
    }
    assert _names(bundle / "clustering") == CLUSTER_FILES
    assert _names(bundle / "static") == STATIC_FILES
    assert _names(bundle / "analysis") == ANALYSIS_FILES
    assert _names(bundle / "figures") == FIGURE_FILES
    assert _names(bundle / "reports") == REPORT_FILES
    manifest = json.loads((bundle / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 1
    assert manifest["configuration"]["path"] == "config/manuscript.yaml"
    assert all(".." not in entry["path"] and not Path(entry["path"]).is_absolute()
               for key in ("inputs", "outputs") for entry in manifest[key])
    assert "run_manifest.json" not in {entry["path"] for entry in manifest["outputs"]}
    output = capsys.readouterr().out
    assert "status" in output


def test_reproduce_preflight_is_nonmutating_and_force_is_manifest_bounded(tmp_path: Path):
    output = tmp_path / "existing"
    output.mkdir()
    (output / "sentinel.txt").write_bytes(b"preserve-me")
    before = _bytes(output)

    assert main(["reproduce", "--config", str(CONFIG), "--output", str(output)]) == 1
    assert _bytes(output) == before

    (output / "run_manifest.json").write_text(json.dumps({
        "schema_version": 1,
        "managed_paths": ["../victim.txt"],
        "outputs": [],
    }), encoding="utf-8")
    before = _bytes(output)
    assert main([
        "reproduce", "--config", str(CONFIG), "--output", str(output), "--force",
    ]) == 1
    assert _bytes(output) == before


def test_force_rejects_undeclared_file_and_reference_output_without_mutation(tmp_path: Path):
    bundle = tmp_path / "bundle"
    assert main(["reproduce", "--config", str(CONFIG), "--output", str(bundle)]) == 0
    assert main([
        "reproduce", "--config", str(CONFIG), "--output", str(bundle), "--force",
    ]) == 0
    manifest_path = bundle / "run_manifest.json"
    original_manifest = manifest_path.read_bytes()
    report_bytes = _bytes(bundle / "reports")
    manifest = json.loads(original_manifest)
    manifest["configuration"]["sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")
    assert main([
        "verify", "--config", str(CONFIG), "--bundle", str(bundle),
        "--snapshot", str(SNAPSHOT),
    ]) == 1
    assert _bytes(bundle / "reports") == report_bytes
    manifest_path.write_bytes(original_manifest)
    manifest = json.loads(original_manifest)
    manifest["summary"]["targets"] += 1
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")
    assert main([
        "verify", "--config", str(CONFIG), "--bundle", str(bundle),
        "--snapshot", str(SNAPSHOT),
    ]) == 1
    assert _bytes(bundle / "reports") == report_bytes
    manifest_path.write_bytes(original_manifest)
    (bundle / "undeclared.txt").write_bytes(b"unknown")
    before = _bytes(bundle)

    assert main([
        "reproduce", "--config", str(CONFIG), "--output", str(bundle), "--force",
    ]) == 1
    assert _bytes(bundle) == before

    reference = TOY / "results" / "reference"
    assert not reference.exists()
    assert main([
        "reproduce", "--config", str(CONFIG), "--output", str(reference), "--force",
    ]) == 1
    assert not reference.exists()


def test_verify_aggregates_manifest_context_stage_and_snapshot_failures(
    tmp_path: Path, capsys,
):
    bundle = tmp_path / "bundle"
    assert main(["reproduce", "--config", str(CONFIG), "--output", str(bundle)]) == 0
    reports_before = _bytes(bundle / "reports")
    manifest_path = bundle / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["configuration"]["sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")
    rankings_path = bundle / "static" / "cluster_v2_static_rankings.csv"
    rankings = pd.read_csv(rankings_path)
    rankings.loc[0, "Within_PDB_rank"] = 99
    rankings.to_csv(rankings_path, index=False, lineterminator="\n")
    snapshot = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    snapshot["maximum_diameter_A"] += 1
    snapshot_path = tmp_path / "bad-snapshot.json"
    snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
    capsys.readouterr()

    assert main([
        "verify", "--config", str(CONFIG), "--bundle", str(bundle),
        "--snapshot", str(snapshot_path),
    ]) == 1
    payload = json.loads(capsys.readouterr().out.strip())
    assert {"manifest_integrity", "manifest_context", "bundle_contract", "snapshot_cluster_summary"}.issubset(
        set(payload["failures"])
    )
    assert _bytes(bundle / "reports") == reports_before


def test_release_scanner_accepts_declared_public_payload():
    completed = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "release_check.py")],
        cwd=ROOT, capture_output=True, text=True, check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["status"] == "pass"
