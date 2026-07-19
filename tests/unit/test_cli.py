from pathlib import Path
import json

import pytest
import yaml

from oips_repro.cli import _safe_error_message, main
from oips_repro.io import path_is_equal_or_within


STABLE_COMMANDS = (
    "validate-data",
    "cluster",
    "score",
    "analyze",
    "figures",
    "verify",
    "reproduce",
)

EXPECTED_PATHS = {
    "schema": "config/schema.json",
    "data_manifest": "data/manifest.tsv",
    "data_checksums": "data/SHA256SUMS",
    "asset_rights": "data/metadata/asset-rights.tsv",
    "external_archive_manifest": "data/external_archive_manifest.tsv",
    "expected_summary": "tests/scientific/data/expected_summary.json",
    "feature_table": "data/static/tool_pocket_features.csv",
    "systems": "data/metadata/systems.tsv",
    "structure_manifest": "data/manifest.tsv",
    "structures_dir": "data/structures",
    "reference_annotations": "data/posthoc/reference_annotations.csv",
    "md_evidence": "data/posthoc/md_evidence.csv",
    "redocking_evidence": "data/posthoc/redocking_evidence.csv",
    "legacy_crosswalk": "legacy/working-cluster-to-cluster-v2.csv",
    "figure_contract": "config/figure_contract.yaml",
    "reference_results": "results/reference",
    "default_output": "results/reproduced",
}

EXPECTED_MODULE_WEIGHTS = {
    "C_cons": 0.22,
    "G_geo": 0.18,
    "P_lig": 0.24,
    "O_rel_formal": 0.24,
    "Q_evidence": 0.12,
}


def test_version_command(capsys):
    assert main(["--version"]) == 0
    assert capsys.readouterr().out.strip().startswith("oips-repro ")


def test_parser_exposes_approved_stable_commands(capsys):
    assert main(["--help"]) == 0
    help_text = capsys.readouterr().out

    for command in STABLE_COMMANDS:
        assert command in help_text


def test_manuscript_configuration_freezes_reviewed_contract():
    config_path = Path(__file__).parents[2] / "config" / "manuscript.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    assert config["schema_version"] == 1
    assert config["paths"] == EXPECTED_PATHS
    assert all(not Path(path).is_absolute() for path in config["paths"].values())
    assert config["scoring"]["module_weights"] == EXPECTED_MODULE_WEIGHTS
    assert sum(config["scoring"]["module_weights"].values()) == pytest.approx(1.0)
    assert (
        config["clustering"]["cross_tool"]["cluster_center_diameter_max_A"]
        == 12.0
    )
    assert config["statistics"] == {
        "bootstrap_iterations": 10000,
        "random_seed": 20260710,
    }


def test_validate_data_prints_compact_aggregated_json(capsys):
    config_path = Path(__file__).parents[2] / "config" / "manuscript.yaml"

    assert main(["validate-data", "--config", str(config_path)]) == 0
    text = capsys.readouterr().out.strip()
    payload = json.loads(text)

    assert set(payload) == {"status", "checks", "failures"}
    assert payload["status"] == "pass" and payload["failures"] == []
    assert payload["checks"] >= 6
    assert "\n" not in text


def test_missing_config_error_does_not_leak_absolute_path_or_traceback(tmp_path: Path, capsys):
    missing = tmp_path / "private-user" / "TOKEN=secret" / "missing.yaml"

    assert main(["validate-data", "--config", str(missing)]) == 1
    error = capsys.readouterr().err

    assert str(tmp_path) not in error
    assert "private-user" not in error and "secret" not in error
    assert "Traceback" not in error


def test_error_sanitizer_redacts_posix_home_paths():
    message = _safe_error_message(FileNotFoundError("/home/private-user/token/input.csv"))

    assert "/home/" not in message and "private-user" not in message


def test_reference_path_classification_includes_nested_bundles(tmp_path: Path):
    reference = tmp_path / "reference"

    assert path_is_equal_or_within(reference, reference)
    assert path_is_equal_or_within(reference / "nested", reference)
    assert not path_is_equal_or_within(tmp_path / "reference-sibling", reference)
