from dataclasses import FrozenInstanceError
import importlib.util
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from oips_repro import io, posthoc, scoring
from oips_repro.config import load_manuscript_config
from oips_repro.provenance import (
    build_run_manifest,
    prepare_output_bundle,
    scan_sensitive_content,
    validate_run_manifest,
    write_run_manifest,
)
from oips_repro.validation import (
    CheckResult, ValidationReport, _safe_repo_path, validate_figure_source_files,
    validate_bundle_relationships, validate_public_inputs, validate_static_science,
)


ROOT = Path(__file__).parents[2]
CONFIG = ROOT / "config" / "manuscript.yaml"


def _release_scanner_module():
    path = ROOT / "scripts" / "release_check.py"
    spec = importlib.util.spec_from_file_location("oips_release_check", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_config(tmp_path: Path, mutate) -> Path:
    data = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    mutate(data)
    project = tmp_path / "project"
    path = project / "config" / "manuscript.yaml"
    path.parent.mkdir(parents=True)
    (project / "config" / "figure_contract.yaml").write_bytes(
        (ROOT / "config" / "figure_contract.yaml").read_bytes()
    )
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


def test_configuration_is_strict_and_repository_relative():
    config = load_manuscript_config(CONFIG)

    assert config.repository_root == ROOT.resolve()
    assert config.resolve_configured_path("asset_rights") == (
        ROOT / "data" / "metadata" / "asset-rights.tsv"
    ).resolve()
    assert tuple(config.data) == (
        "schema_version", "analysis", "paths", "tool_weights", "clustering",
        "scoring", "ranking", "posthoc", "statistics",
    )


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda data: data.update({"surprise": 1}), "unknown.*surprise"),
        (lambda data: data["scoring"].update({"surprise": 1}), "unknown.*surprise"),
        (lambda data: data["statistics"].__setitem__("bootstrap_iterations", True), "bootstrap_iterations"),
        (lambda data: data["statistics"].__setitem__("random_seed", True), "random_seed"),
        (lambda data: data["tool_weights"].__setitem__("SiteMap", float("nan")), "finite"),
        (lambda data: data["scoring"]["module_weights"].__setitem__("C_cons", float("inf")), "finite"),
        (lambda data: data["clustering"]["cross_tool"].__setitem__("close_center_max_A", 11.0), "close.*conditional"),
        (lambda data: data["posthoc"]["reference"].__setitem__("near_dcc_max_A", 13.0), "near.*middle.*far"),
        (lambda data: data["posthoc"]["md"].__setitem__("concordant_center_max_A", 9.0), "center"),
        (lambda data: data["paths"].__setitem__("schema", "../schema.json"), r"relative|\.\."),
        (lambda data: data["paths"].__setitem__("schema", "C:/private/schema.json"), "relative"),
        (lambda data: data["paths"].__setitem__("schema", "config\\schema.json"), "POSIX|backslash"),
    ],
)
def test_configuration_rejects_adversarial_values(tmp_path: Path, mutate, message):
    with pytest.raises(ValueError, match=message):
        load_manuscript_config(_write_config(tmp_path, mutate))


def test_configuration_rejects_missing_nested_key(tmp_path: Path):
    path = _write_config(
        tmp_path,
        lambda data: data["posthoc"]["reference"].pop("near_recall_min"),
    )

    with pytest.raises(ValueError, match="missing.*near_recall_min"):
        load_manuscript_config(path)


def test_configuration_requires_a_regular_figure_contract(tmp_path: Path):
    path = _write_config(tmp_path, lambda data: None)
    (path.parent / "figure_contract.yaml").unlink()

    with pytest.raises((FileNotFoundError, ValueError), match="figure contract"):
        load_manuscript_config(path)


def test_validation_types_are_frozen_and_aggregate_all_failures():
    checks = (
        CheckResult("manifest", "fail", "bad manifest", {"count": 2}),
        CheckResult("checksums", "fail", "bad checksum", {"count": 1}),
        CheckResult("optional_md", "warning", "MD absent", {}),
        CheckResult("schema", "pass", "valid", {}),
    )
    report = ValidationReport(checks)

    assert (report.pass_count, report.failure_count, report.warning_count) == (1, 2, 1)
    with pytest.raises(ValueError, match="manifest.*checksums"):
        report.raise_for_failures()
    with pytest.raises(FrozenInstanceError):
        checks[0].status = "pass"
    with pytest.raises(ValueError, match="status"):
        CheckResult("bad", "unknown", "bad", {})


def test_provenance_manifest_has_stable_safe_relative_inventory(tmp_path: Path):
    config = load_manuscript_config(CONFIG)
    bundle = tmp_path / "bundle"
    output = bundle / "analysis" / "result.csv"
    output.parent.mkdir(parents=True)
    output.write_text("value\n1\n", encoding="utf-8", newline="\n")

    manifest = build_run_manifest(
        config,
        bundle,
        consumed_inputs=[config.resolve_configured_path("feature_table")],
        output_paths=[output],
        managed_paths=["analysis/result.csv", "run_manifest.json"],
        summary={
            "records": 2, "mapped": 1, "excluded": 1, "same_tool_units": 1,
            "clusters": 1, "targets": 1, "boundary_sensitive": 0,
            "maximum_diameter_A": 0.0,
            "maximum_formal_votes_per_cluster_tool": 1,
        },
        generated_at_utc="2026-07-17T00:00:00Z",
    )
    written = write_run_manifest(bundle, manifest)

    assert written == bundle / "run_manifest.json"
    assert manifest["schema_version"] == 1
    assert manifest["generated_at_utc"] == "2026-07-17T00:00:00Z"
    assert manifest["configuration"]["path"] == "config/manuscript.yaml"
    assert list(manifest["dependencies"]) == sorted(manifest["dependencies"])
    assert manifest["outputs"][0]["path"] == "analysis/result.csv"
    assert "run_manifest.json" not in {entry["path"] for entry in manifest["outputs"]}
    assert "run_manifest.json" in manifest["managed_paths"]
    assert manifest["summary"]["mapped"] + manifest["summary"]["excluded"] == manifest["summary"]["records"]
    assert validate_run_manifest(bundle) == manifest
    assert written.read_bytes().endswith(b"\n")


@pytest.mark.parametrize(
    "mutation",
    [
        lambda manifest: manifest["managed_paths"].append("../victim.txt"),
        lambda manifest: manifest["outputs"][0].__setitem__("sha256", "bad"),
        lambda manifest: manifest["outputs"][0].__setitem__("path", "C:/private/result.csv"),
        lambda manifest: manifest["outputs"].append({"path": "missing.csv", "sha256": "0" * 64, "bytes": 1}),
        lambda manifest: manifest["package"].update({"unexpected": "value"}),
    ],
)
def test_manifest_validation_rejects_traversal_malformed_and_missing_files(
    tmp_path: Path, mutation
):
    config = load_manuscript_config(CONFIG)
    bundle = tmp_path / "bundle"
    output = bundle / "analysis" / "result.csv"
    output.parent.mkdir(parents=True)
    output.write_text("value\n1\n", encoding="utf-8")
    manifest = build_run_manifest(
        config, bundle,
        consumed_inputs=[config.resolve_configured_path("feature_table")],
        output_paths=[output], managed_paths=["analysis/result.csv", "run_manifest.json"],
        summary={
            "records": 1, "mapped": 1, "excluded": 0, "same_tool_units": 1,
            "clusters": 1, "targets": 1, "boundary_sensitive": 0,
            "maximum_diameter_A": 0.0,
            "maximum_formal_votes_per_cluster_tool": 1,
        },
        generated_at_utc="2026-07-17T00:00:00Z",
    )
    mutation(manifest)
    write_run_manifest(bundle, manifest)

    with pytest.raises(ValueError):
        validate_run_manifest(bundle)


def test_sensitive_scan_collects_multiple_failure_classes():
    result = scan_sensitive_content({
        "reports/a.md": "TOKEN=secret-value\nfile:///private/input\n",
        "reports/b.md": "C:\\Users\\person\\data\nC:/private/data\nhttps://host/session/abc/result/123\nhttps://host/token/abc\n",
    })

    assert result.status == "fail"
    assert result.evidence["match_count"] >= 6
    assert {"credential", "file_uri", "drive_path", "session_url"}.issubset(
        set(result.evidence["classes"])
    )


def test_output_root_symlink_is_rejected_before_resolution(tmp_path: Path):
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "linked-output"
    try:
        link.symlink_to(real, target_is_directory=True)
    except OSError:
        pytest.skip("directory symlinks are unavailable on this Windows host")

    with pytest.raises(ValueError, match="symlink|junction"):
        prepare_output_bundle(link, tmp_path / "reference")
    with pytest.raises(ValueError, match="symlink|junction"):
        validate_run_manifest(link)


def test_output_parent_symlink_is_rejected_before_resolution(tmp_path: Path):
    real = tmp_path / "real-parent"
    real.mkdir()
    link = tmp_path / "linked-parent"
    try:
        link.symlink_to(real, target_is_directory=True)
    except OSError:
        pytest.skip("directory symlinks are unavailable on this Windows host")

    output = link / "new-bundle"
    with pytest.raises(ValueError, match="symlink|junction"):
        prepare_output_bundle(output, tmp_path / "reference")
    assert not (real / "new-bundle").exists()


def test_gitkeep_symlink_is_not_treated_as_an_empty_output(tmp_path: Path):
    output = tmp_path / "output"
    output.mkdir()
    target = tmp_path / "marker"
    target.write_text("external", encoding="utf-8")
    try:
        (output / ".gitkeep").symlink_to(target)
    except OSError:
        pytest.skip("file symlinks are unavailable on this Windows host")

    with pytest.raises(ValueError, match="symlink|junction"):
        prepare_output_bundle(output, tmp_path / "reference")
    assert target.read_text(encoding="utf-8") == "external"


def test_broken_output_symlink_is_rejected_without_creating_target(tmp_path: Path):
    output = tmp_path / "broken-output"
    target = tmp_path / "missing-target"
    try:
        output.symlink_to(target, target_is_directory=True)
    except OSError:
        pytest.skip("directory symlinks are unavailable on this Windows host")

    with pytest.raises(ValueError, match="symlink|junction"):
        prepare_output_bundle(output, tmp_path / "reference")
    assert not target.exists()


def test_public_release_validation_includes_manual_decisions_schema():
    report = validate_public_inputs(load_manuscript_config(CONFIG))
    tables = next(check for check in report.checks if check.check_id == "release_tables")

    assert tables.status == "pass"
    assert tables.evidence["tables"] == 9


def test_repository_path_rejects_symlink_before_resolving(tmp_path: Path):
    target = tmp_path / "target"
    target.mkdir()
    (target / "payload.txt").write_text("x", encoding="utf-8")
    link = tmp_path / "linked"
    try:
        link.symlink_to(target, target_is_directory=True)
    except OSError:
        pytest.skip("directory symlinks are unavailable on this Windows host")

    with pytest.raises(ValueError, match="symlink|junction"):
        _safe_repo_path("linked/payload.txt", tmp_path, "payload")


@pytest.mark.parametrize(
    ("content", "failure_class"),
    [
        ("API_KEY=not-a-public-value\n", "credential_assignment"),
        ("C:/Users/person/private/result.csv\n", "drive_path"),
        (r"\\server\share\private\result.csv", "unc_path"),
        ("/tmp/private/result.csv\n", "private_posix_path"),
        ("https://example.test/api?job_id=private-run\n", "session_query"),
    ],
)
def test_release_scanner_rejects_sensitive_attack_payloads(
    tmp_path: Path, content: str, failure_class: str,
):
    scanner = _release_scanner_module()
    payload = tmp_path / "payload.txt"
    payload.write_text(content, encoding="utf-8")

    failures = scanner._scan_tracked_file(tmp_path, "payload.txt")

    assert {failure["check"] for failure in failures} >= {failure_class}


@pytest.mark.parametrize(
    "relative",
    ["scripts/release_check.py", "src/oips_repro/provenance.py"],
)
def test_release_scanner_does_not_match_its_own_security_code(relative: str):
    scanner = _release_scanner_module()

    failures = scanner._scan_tracked_file(ROOT, relative)

    assert failures == []


def test_release_scanner_rejects_malformed_release_manifest(tmp_path: Path):
    scanner = _release_scanner_module()
    release = tmp_path / "release"
    release.mkdir()
    (release / "manifest.json").write_text('{"schema_version":1}\n', encoding="utf-8")

    failures = scanner._validate_release_manifest(
        tmp_path, "release/manifest.json", {"release/manifest.json"}, required=True,
    )

    assert failures == [{
        "check": "release_manifest_root_schema", "path": "release/manifest.json",
    }]


def test_release_scanner_exact_root_schema_fails_closed_without_crashing(tmp_path: Path):
    scanner = _release_scanner_module()
    release = tmp_path / "release"
    release.mkdir()
    document = {key: None for key in scanner.RELEASE_KEYS}
    (release / "manifest.json").write_text(json.dumps(document), encoding="utf-8")

    failures = scanner._validate_release_manifest(
        tmp_path, "release/manifest.json", {"release/manifest.json"}, required=True,
    )

    assert failures
    assert {failure["check"] for failure in failures} >= {
        "release_manifest_schema_version", "release_manifest_inventory",
    }


def _release_entry(path: Path, root: Path) -> dict[str, object]:
    content = path.read_bytes()
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": hashlib.sha256(content).hexdigest(),
        "bytes": len(content),
    }


def test_release_scanner_accepts_manuscript_tag_and_enforces_version_consistency(
    tmp_path: Path,
):
    scanner = _release_scanner_module()
    paths = [
        tmp_path / "config" / "manuscript.yaml",
        tmp_path / "environment" / "constraints.txt",
        tmp_path / "data" / "input.csv",
        tmp_path / "results" / "reference" / "result.csv",
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"public release fixture: {path.name}\n", encoding="utf-8")
    release = tmp_path / "release" / "manifest.json"
    release.parent.mkdir()
    document = {
        "schema_version": 1,
        "generated_at_utc": "2026-07-17T00:00:00Z",
        "version": "1.0.0",
        "release_status": "pre_publication_incomplete",
        "git": {"commit": None, "dirty": True},
        "release_tag": {"tag": "v1.0.0-manuscript", "status": "not_created"},
        "identifiers": {
            "repository_url": None, "code_doi": None,
            "data_doi": None, "manuscript_doi": None,
        },
        "configuration": _release_entry(paths[0], tmp_path),
        "inputs": [_release_entry(paths[2], tmp_path)],
        "environment": _release_entry(paths[1], tmp_path),
        "key_results": [_release_entry(paths[3], tmp_path)],
        "publication_payload": [_release_entry(paths[3], tmp_path)],
    }
    release.write_text(json.dumps(document), encoding="utf-8")
    tracked = {path.relative_to(tmp_path).as_posix() for path in paths} | {
        "release/manifest.json"
    }

    assert scanner._validate_release_manifest(
        tmp_path, "release/manifest.json", tracked, required=True,
        expected_inputs={"data/input.csv"},
    ) == []

    document["release_tag"]["tag"] = "v0.9.0-manuscript"
    release.write_text(json.dumps(document), encoding="utf-8")
    failures = scanner._validate_release_manifest(
        tmp_path, "release/manifest.json", tracked, required=True,
        expected_inputs={"data/input.csv"},
    )
    assert {failure["check"] for failure in failures} >= {
        "release_manifest_tag_version"
    }


def _synthetic_static_tables():
    rows = []
    for index, module_value in enumerate((20.0, 10.0), start=1):
        row = {column: 0.0 for column in io.MASTER_COLUMNS}
        row.update({
            "pdb_id": "TEST", "cluster_v2_id": f"TEST_V2C00{index}",
            "C_cons": module_value, "G_geo": module_value,
            "P_lig": module_value, "O_rel_formal": module_value,
            "Q_evidence": module_value, "OIPS-P_static": module_value,
        })
        rows.append(row)
    master = pd.DataFrame(rows, columns=io.MASTER_COLUMNS)
    config = load_manuscript_config(CONFIG)
    weights = dict(config.data["scoring"]["module_weights"])
    master.attrs.update(module_weights=weights, ranking_direction="descending")
    rankings = scoring.rank_within_target(master, rank_method="average")
    return config, master, rankings


def test_static_science_validation_recomputes_formula_rank_and_shared_fields():
    config, master, rankings = _synthetic_static_tables()
    validate_static_science(master, rankings, config)

    bad_rank = rankings.copy()
    bad_rank.loc[0, "Within_PDB_rank"] = 99
    with pytest.raises(ValueError, match="ranking"):
        validate_static_science(master, bad_rank, config)

    bad_master = master.copy()
    bad_master.loc[0, "C_cons"] += 1
    with pytest.raises(ValueError, match="master|formula"):
        validate_static_science(bad_master, rankings, config)


def test_bundle_relationships_enforce_diameter_and_cross_table_identity():
    config, _master, rankings = _synthetic_static_tables()
    keys = rankings.loc[:, ["pdb_id", "cluster_v2_id"]]
    labels = keys.copy()
    references = keys.copy()
    for column in posthoc.CONVENIENCE_COLUMNS[:4]:
        labels[column] = False if column == "unresolved_flag" else "label"
    for column in posthoc.CONVENIENCE_COLUMNS[4:]:
        references[column] = False if column == "R_auto_rule_pass" else 1.0
    final = rankings.copy()
    for column in posthoc.CONVENIENCE_COLUMNS[:4]:
        final[column] = labels[column]
    for column in posthoc.CONVENIENCE_COLUMNS[4:]:
        final[column] = references[column]
    frames = {
        "final_cluster_v2_master_table.csv": final,
        "final_automated_evidence_labels.csv": labels,
        "final_reference_mapping.csv": references,
    }
    candidates = rankings.loc[:, list(io.CANDIDATE_COLUMNS)].copy()
    candidates["cluster_diameter_A"] = [1.0, 2.0]
    candidates.loc[0, "cluster_chain_entropy"] = float("nan")
    candidates.loc[1, "cluster_chain_entropy"] = 0.75
    rankings = rankings.copy()
    rankings["cluster_diameter_A"] = candidates["cluster_diameter_A"]
    rankings.loc[1, "cluster_chain_entropy"] = 0.75
    rankings = rankings.iloc[::-1].reset_index(drop=True)
    frames["final_cluster_v2_master_table.csv"]["cluster_diameter_A"] = candidates["cluster_diameter_A"]
    frames["final_cluster_v2_master_table.csv"].loc[1, "cluster_chain_entropy"] = 0.75
    validate_bundle_relationships(candidates, rankings, frames, config)

    too_wide = candidates.copy()
    too_wide.loc[0, "cluster_diameter_A"] = 13.0
    with pytest.raises(ValueError, match="diameter"):
        validate_bundle_relationships(too_wide, rankings, frames, config)
    inconsistent_candidate = candidates.copy()
    inconsistent_candidate.loc[0, "tool_support_count"] = 99
    with pytest.raises(ValueError, match="clustering.*static"):
        validate_bundle_relationships(inconsistent_candidate, rankings, frames, config)
    bad_entropy = rankings.copy()
    bad_entropy.loc[0, "cluster_chain_entropy"] = 0.25
    with pytest.raises(ValueError, match="entropy"):
        validate_bundle_relationships(candidates, bad_entropy, frames, config)
    bad_final = final.copy()
    bad_final.loc[0, "OIPS-P_static"] += 1
    with pytest.raises(ValueError, match="static ranking"):
        validate_bundle_relationships(candidates, rankings, {**frames, "final_cluster_v2_master_table.csv": bad_final}, config)
    bad_label = labels.copy()
    bad_label.loc[0, "automated_evidence_label"] = "different"
    with pytest.raises(ValueError, match="evidence labels"):
        validate_bundle_relationships(candidates, rankings, {**frames, "final_automated_evidence_labels.csv": bad_label}, config)


def test_figure_source_validation_rejects_value_tampering(tmp_path: Path):
    source = pd.DataFrame({
        "panel": ["a"], "record_id": ["row"], "series": ["metric"],
        "group": ["group"], "x": [1.25], "y": [2.0], "x_label": ["x"],
        "y_label": ["y"], "lower": [1.0], "upper": [1.5],
        "annotation": [""],
    })
    io.write_stable_csv(
        source, tmp_path / "figure_source_data.csv",
        columns=tuple(source.columns), sort_by=("panel", "record_id", "series"),
    )
    validate_figure_source_files(tmp_path, {"figure": source})

    tampered = source.copy()
    tampered.loc[0, "x"] = 9.0
    tampered.to_csv(tmp_path / "figure_source_data.csv", index=False)
    with pytest.raises(ValueError, match="source data"):
        validate_figure_source_files(tmp_path, {"figure": source})
