import csv
import hashlib
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path
from urllib.parse import urlparse

import jsonschema
import pandas as pd
import pytest
import yaml


ROOT = Path(__file__).parents[2]
DATA = ROOT / "data"

FEATURE_COLUMNS = [
    "row_id",
    "pdb_id",
    "tool",
    "pocket_id",
    "display_order",
    "sitemap_rank",
    "center_x",
    "center_y",
    "center_z",
    "center_method",
    "residue_count",
    "residue_set_json",
    "pocket_geometry_score",
    "pocket_ligandability_score",
]
SYSTEM_COLUMNS = [
    "pdb_id",
    "protein_system",
    "oligomeric_state",
    "reference_ligand_annotation",
    "dataset_role",
    "pocket_category",
    "category_mapping_status",
    "protein_family",
    "structure_source",
    "structure_source_version_or_date",
    "biological_assembly_id",
    "chain_mapping",
    "retained_hetero_groups",
    "deleted_hetero_groups",
    "missing_residue_atom_handling",
    "protonation_tool",
    "protonation_ph",
    "force_field",
    "water_model",
    "md_engine",
    "md_engine_version",
    "md_timestep_fs",
    "md_ensemble",
    "md_random_seed",
    "trajectory_alignment",
    "trajectory_sampling",
    "pocket_tools_and_versions",
    "id_mapping_provenance",
    "metadata_completeness",
    "metadata_status",
]
REFERENCE_COLUMNS = [
    "pdb_id",
    "selected_ligand",
    "ligand_status",
    "reference_ligand_annotation",
    "override_applied",
    "decision_id",
]
MD_COLUMNS = [
    "pdb_id",
    "MD_run",
    "Simulation_context",
    "Persistent_MD_contact_residues",
    "MD_contact_center",
    "D_dyn_run_score",
]
REDOCKING_COLUMNS = [
    "pdb_id",
    "Docking_software",
    "GlideScore_kcal_per_mol",
    "Raw_ligand_RMSD_A",
    "RMSD_threshold_call",
    "Failure_or_warning_reason",
]
LEGACY_COLUMNS = [
    "pdb_id",
    "working_cluster_id",
    "cluster_v2_id",
    "mapped_record_count",
    "working_cluster_member_count",
    "mapping_fraction",
    "unmappable_record_count",
    "mapping_type",
]
RIGHTS_COLUMNS = [
    "asset_class",
    "scientific_role",
    "source_provider",
    "content_scope",
    "repository_action",
    "review_status",
    "license_scope",
    "terms_or_evidence_url",
    "approval_basis",
    "review_date",
    "notes",
]
DECISION_COLUMNS = [
    "decision_id",
    "pdb_id_or_scope",
    "decision_type",
    "field_name",
    "previous_value",
    "released_value",
    "rationale",
    "evidence_source",
    "approval_status",
]
MANIFEST_COLUMNS = [
    "asset_id",
    "scientific_role",
    "local_path_or_url",
    "repository",
    "persistent_id",
    "sha256",
    "bytes",
    "media_type",
    "license_or_terms",
    "access_status",
    "source_version_or_date",
]
EXTERNAL_COLUMNS = [
    "asset_id",
    "pdb_id_or_scope",
    "scientific_role",
    "archive_filename_or_label",
    "local_availability",
    "repository",
    "persistent_id",
    "bytes",
    "sha256",
    "access_status",
    "completeness",
    "license_status",
    "source_version_or_date",
    "notes",
]

TABLE_CONTRACTS = {
    DATA / "static" / "tool_pocket_features.csv": (
        FEATURE_COLUMNS,
        ",",
        "tool_pocket_features_row",
    ),
    DATA / "metadata" / "systems.tsv": (SYSTEM_COLUMNS, "\t", "systems_row"),
    DATA / "metadata" / "asset-rights.tsv": (
        RIGHTS_COLUMNS,
        "\t",
        "asset_rights_row",
    ),
    DATA / "metadata" / "manual_decisions.tsv": (
        DECISION_COLUMNS,
        "\t",
        "manual_decisions_row",
    ),
    DATA / "posthoc" / "reference_annotations.csv": (
        REFERENCE_COLUMNS,
        ",",
        "reference_annotations_row",
    ),
    DATA / "posthoc" / "md_evidence.csv": (MD_COLUMNS, ",", "md_evidence_row"),
    DATA / "posthoc" / "redocking_evidence.csv": (
        REDOCKING_COLUMNS,
        ",",
        "redocking_evidence_row",
    ),
    DATA / "external_archive_manifest.tsv": (
        EXTERNAL_COLUMNS,
        "\t",
        "external_archive_manifest_row",
    ),
    DATA / "manifest.tsv": (MANIFEST_COLUMNS, "\t", "manifest_row"),
    ROOT / "legacy" / "working-cluster-to-cluster-v2.csv": (
        LEGACY_COLUMNS,
        ",",
        "legacy_crosswalk_row",
    ),
}

ANALYSIS_TABLES = [
    DATA / "static" / "tool_pocket_features.csv",
    DATA / "metadata" / "systems.tsv",
    DATA / "posthoc" / "reference_annotations.csv",
    DATA / "posthoc" / "md_evidence.csv",
    DATA / "posthoc" / "redocking_evidence.csv",
    ROOT / "legacy" / "working-cluster-to-cluster-v2.csv",
]


def read_records(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), list(reader)


def as_schema_value(value: str, definition: dict):
    if "oneOf" in definition:
        for option in definition["oneOf"]:
            if option.get("const") == value:
                return value
        for option in definition["oneOf"]:
            if option.get("type") == "integer" and value != "":
                return int(value)
    kinds = definition.get("type", "string")
    kinds = [kinds] if isinstance(kinds, str) else kinds
    if value == "" and "null" in kinds:
        return None
    if "integer" in kinds:
        return int(value)
    if "number" in kinds:
        return float(value)
    if "boolean" in kinds:
        assert value in {"true", "false"}
        return value == "true"
    return value


def test_release_files_exist_and_headers_match_exact_contracts():
    assert (DATA / "README.md").is_file()
    assert (ROOT / "config" / "schema.json").is_file()
    assert (ROOT / "scripts" / "prepare_public_data.py").is_file()
    for path, (columns, delimiter, _) in TABLE_CONTRACTS.items():
        assert path.is_file(), path
        header, _ = read_records(path, delimiter)
        assert header == columns


def test_feature_snapshot_has_expected_rows_targets_types_and_sorting():
    path = DATA / "static" / "tool_pocket_features.csv"
    frame = pd.read_csv(path, keep_default_na=False)

    assert len(frame) == 1742
    assert frame["pdb_id"].nunique() == 21
    assert frame["row_id"].is_unique
    assert frame["row_id"].tolist() == sorted(frame["row_id"].tolist()) or list(
        frame[["pdb_id", "tool", "row_id"]].itertuples(index=False, name=None)
    ) == sorted(frame[["pdb_id", "tool", "row_id"]].itertuples(index=False, name=None))
    assert list(frame[["pdb_id", "tool", "row_id"]].itertuples(index=False, name=None)) == sorted(
        frame[["pdb_id", "tool", "row_id"]].itertuples(index=False, name=None)
    )
    assert frame["row_id"].map(type).eq(int).all()
    assert frame["residue_count"].map(type).eq(int).all()
    for raw in frame["residue_set_json"]:
        assert isinstance(json.loads(raw), list)


def test_target_tables_and_reviewed_5lge_reference_are_complete():
    systems = pd.read_csv(DATA / "metadata" / "systems.tsv", sep="\t", dtype=str)
    refs = pd.read_csv(
        DATA / "posthoc" / "reference_annotations.csv",
        dtype=str,
        keep_default_na=False,
    )

    assert len(systems) == systems["pdb_id"].nunique() == 21
    assert len(refs) == refs["pdb_id"].nunique() == 21
    assert systems["pdb_id"].tolist() == sorted(systems["pdb_id"])
    assert refs["pdb_id"].tolist() == sorted(refs["pdb_id"])
    unavailable_fields = SYSTEM_COLUMNS[9:28]
    assert (systems[unavailable_fields] == "not_available_in_source_materials").all().all()
    assert (systems["metadata_status"] == "author_confirmation_required").all()
    reviewed = refs.set_index("pdb_id").loc["5LGE"]
    assert reviewed["reference_ligand_annotation"] == "6VN:D:503:"
    assert reviewed["override_applied"] == "true"
    assert reviewed["decision_id"] != ""


def test_posthoc_and_legacy_tables_are_minimal_unique_and_stably_sorted():
    md = pd.read_csv(DATA / "posthoc" / "md_evidence.csv", keep_default_na=False)
    redocking = pd.read_csv(
        DATA / "posthoc" / "redocking_evidence.csv", keep_default_na=False
    )
    legacy = pd.read_csv(
        ROOT / "legacy" / "working-cluster-to-cluster-v2.csv",
        keep_default_na=False,
    )

    assert len(md) > 0
    assert not md.duplicated(["pdb_id", "Simulation_context", "MD_run"]).any()
    assert list(md[["pdb_id", "Simulation_context", "MD_run"]].itertuples(index=False, name=None)) == sorted(
        md[["pdb_id", "Simulation_context", "MD_run"]].itertuples(index=False, name=None)
    )
    assert len(redocking) == redocking["pdb_id"].nunique() == 21
    assert redocking["pdb_id"].tolist() == sorted(redocking["pdb_id"])
    assert "Raw_data_directory" not in md.columns
    assert "Source" not in redocking.columns
    assert len(legacy) == 911
    assert legacy["working_cluster_id"].nunique() == 369
    assert (legacy["cluster_v2_id"] == "").sum() == 178
    assert legacy.loc[legacy["cluster_v2_id"] != "", "cluster_v2_id"].nunique() == 733
    assert not legacy.duplicated(["pdb_id", "working_cluster_id", "cluster_v2_id"]).any()
    keys = list(
        legacy[["pdb_id", "working_cluster_id", "cluster_v2_id"]].itertuples(
            index=False, name=None
        )
    )
    assert keys == sorted(keys, key=lambda row: (row[0], row[1], row[2] == "", row[2]))
    assert not (
        {"member_rows", "member_pockets", "tools", "OIPS_S2_current", "OIPS_V1_score"}
        & set(legacy.columns)
    )


def test_json_schema_is_valid_and_validates_representative_rows():
    schema = json.loads((ROOT / "config" / "schema.json").read_text(encoding="utf-8"))
    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    assert schema["schema_version"] == 1
    jsonschema.Draft202012Validator.check_schema(schema)

    for path, (columns, delimiter, definition_name) in TABLE_CONTRACTS.items():
        _, rows = read_records(path, delimiter)
        definition = schema["$defs"][definition_name]
        assert definition["x-columns"] == columns
        validator = jsonschema.Draft202012Validator(
            definition, format_checker=jsonschema.FormatChecker()
        )
        for row in range(len(rows)):
            typed = {
                key: as_schema_value(value, definition["properties"][key])
                for key, value in rows[row].items()
            }
            validator.validate(typed)


def test_analysis_tables_do_not_leak_internal_headers_paths_or_urls():
    forbidden_header = re.compile(r"path|file|dir|url|sess|job|run_id", re.I)
    forbidden_value = re.compile(r"[A-Z]:\\|file://|https?://", re.I)
    missing_token = re.compile(r"^(?:na|n/a|nan|none|null)$", re.I)

    for path in ANALYSIS_TABLES:
        delimiter = "\t" if path.suffix == ".tsv" else ","
        header, rows = read_records(path, delimiter)
        assert not any(forbidden_header.search(column) for column in header)
        for row in rows:
            assert not any(forbidden_value.search(value) for value in row.values())
            assert not any(missing_token.fullmatch(value.strip()) for value in row.values())

    _, decisions = read_records(DATA / "metadata" / "manual_decisions.tsv", "\t")
    for row in decisions:
        assert not any(forbidden_value.search(value) for value in row.values())


def test_manual_decision_audit_covers_all_required_reviewed_choices():
    _, rows = read_records(DATA / "metadata" / "manual_decisions.tsv", "\t")
    assert len(rows) == 26
    assert len({row["decision_id"] for row in rows}) == len(rows)
    assert rows == sorted(rows, key=lambda row: row["decision_id"])
    assert all(all(row[column] for column in DECISION_COLUMNS) for row in rows)
    assert all(row["approval_status"] == "approved_for_repository" for row in rows)

    decisions_by_id = {row["decision_id"]: row for row in rows}
    ref_row = next(row for row in rows if row["field_name"] == "reference_ligand_annotation")
    assert ref_row["pdb_id_or_scope"] == "5LGE"
    assert ref_row["previous_value"] == "NAP:C:501:"
    assert ref_row["released_value"] == "6VN:D:503:"
    _, reference_rows = read_records(
        DATA / "posthoc" / "reference_annotations.csv", ","
    )
    released_5lge = next(row for row in reference_rows if row["pdb_id"] == "5LGE")
    assert released_5lge["decision_id"] == ref_row["decision_id"]
    assert decisions_by_id[released_5lge["decision_id"]] == ref_row

    categories = [row for row in rows if row["decision_type"] == "category_annotation"]
    assert len(categories) == 21
    assert len({row["pdb_id_or_scope"] for row in categories}) == 21
    parameters = [row for row in rows if row["field_name"] == "cluster_v2_parameter_set"]
    assert len(parameters) == 1
    cases = [row for row in rows if row["decision_type"] == "representative_case_selection"]
    assert {row["pdb_id_or_scope"] for row in cases} == {"5J89", "5TBM", "4W9H"}


def test_rights_gate_approves_every_local_payload_and_keeps_restricted_assets_manifest_only():
    _, rights_rows = read_records(DATA / "metadata" / "asset-rights.tsv", "\t")
    _, manifest_rows = read_records(DATA / "manifest.tsv", "\t")
    approved_by_role = {
        row["scientific_role"]: row
        for row in rights_rows
        if row["review_status"] == "approved_for_repository"
    }

    for asset in manifest_rows:
        right = approved_by_role[asset["scientific_role"]]
        assert right["license_scope"]
        assert right["terms_or_evidence_url"]
        assert right["repository_action"] == "commit_local_payload"

    restricted = [
        row
        for row in rights_rows
        if row["review_status"] in {"manifest_only", "excluded"}
    ]
    assert restricted
    local_roles = {row["scientific_role"] for row in manifest_rows}
    assert all(row["scientific_role"] not in local_roles for row in restricted)
    assert all(row["repository_action"] != "commit_local_payload" for row in restricted)

    _, external_rows = read_records(DATA / "external_archive_manifest.tsv", "\t")
    rights_by_role = {row["scientific_role"]: row for row in rights_rows}
    for asset in external_rows:
        right = rights_by_role[asset["scientific_role"]]
        assert right["review_status"] in {"manifest_only", "excluded"}
        assert right["repository_action"] != "commit_local_payload"
        assert right["license_scope"]
        assert right["terms_or_evidence_url"]


def test_manifest_exactly_matches_local_payload_and_all_structure_hashes():
    _, rows = read_records(DATA / "manifest.tsv", "\t")
    assert len({row["asset_id"] for row in rows}) == len(rows)
    listed = {row["local_path_or_url"] for row in rows}
    intended = {
        path.relative_to(ROOT).as_posix()
        for path in DATA.rglob("*")
        if path.is_file() and path.name not in {"manifest.tsv", "SHA256SUMS"}
    }
    assert listed == intended

    for row in rows:
        relative = Path(row["local_path_or_url"])
        assert not relative.is_absolute()
        path = ROOT / relative
        content = path.read_bytes()
        assert int(row["bytes"]) == len(content)
        assert row["sha256"] == hashlib.sha256(content).hexdigest()
        for value in row.values():
            for url in re.findall(r"https?://[^\s;,]+", value):
                parsed = urlparse(url)
                assert parsed.scheme == "https"
                assert not parsed.query
                assert not re.search(r"sess|job|result|token", url, re.I)

    structures = [row for row in rows if row["scientific_role"] == "prepared paired structure"]
    assert len(structures) == 21
    assert all(row["local_path_or_url"].startswith("data/structures/") for row in structures)


def test_external_inventory_covers_required_non_git_asset_classes_without_internal_paths():
    _, rows = read_records(DATA / "external_archive_manifest.tsv", "\t")
    assert len({row["asset_id"] for row in rows}) == len(rows)
    assert rows == sorted(rows, key=lambda row: row["asset_id"])
    roles = {row["scientific_role"] for row in rows}
    assert {
        "full MD inputs",
        "full MD trajectories",
        "restart or checkpoint files",
        "representative MD snapshots",
        "large raw pocket outputs",
    } <= roles
    forbidden = re.compile(r"[A-Z]:\\|file://|https?://.*(?:sess|job|result|token)", re.I)
    assert not any(forbidden.search(value) for row in rows for value in row.values())
    for row in rows:
        if row["local_availability"] != "available_locally":
            assert row["bytes"] == "not_available_in_source_materials"
            assert row["sha256"] == "not_available_in_source_materials"


def test_sha256sums_covers_every_data_file_except_itself():
    checksum_path = DATA / "SHA256SUMS"
    parsed = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        parsed[relative] = digest

    intended = {
        path.relative_to(ROOT).as_posix()
        for path in DATA.rglob("*")
        if path.is_file() and path != checksum_path
    }
    assert set(parsed) == intended
    for relative, digest in parsed.items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest


def test_release_files_use_lf_and_stay_below_git_size_limit():
    for path in DATA.rglob("*"):
        if not path.is_file():
            continue
        assert path.stat().st_size <= 10 * 1024 * 1024
        if path.suffix.lower() in {".csv", ".tsv", ".md"} or path.name == "SHA256SUMS":
            assert b"\r\n" not in path.read_bytes()


def test_publication_metadata_version_license_citation_and_constraints_are_release_ready():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert project["project"]["version"] == "1.0.0"
    assert "cffconvert>=2,<3" in project["project"]["optional-dependencies"]["release"]
    assert (ROOT / "src" / "oips_repro" / "__init__.py").read_text(
        encoding="utf-8"
    ).count('__version__ = "1.0.0"') == 1

    constraints = (ROOT / "environment" / "constraints.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    pins = [line for line in constraints if line and not line.startswith("#")]
    package_name = lambda pin: re.sub(r"[-_.]+", "-", pin.split("==", 1)[0]).casefold()
    assert pins == sorted(pins, key=package_name)
    assert all(re.fullmatch(r"[A-Za-z0-9_.-]+==[^=<>!~]+", pin) for pin in pins)
    for required in ("numpy", "pandas", "scipy", "matplotlib", "PyYAML", "jsonschema", "pytest", "cffconvert"):
        assert any(pin.casefold().startswith(required.casefold() + "==") for pin in pins)

    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert license_text.startswith("BSD 3-Clause License\n")
    assert "Copyright (c) 2026" in license_text
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    assert citation["cff-version"] == "1.2.0"
    assert citation["version"] == "1.0.0"
    assert citation["title"] == "OIPS Reproducibility Package"
    assert [author["family-names"] for author in citation["authors"]] == [
        "Chen", "Zhu", "Zhao", "Zhang", "Sun", "Chen", "Xue",
    ]
    assert citation["authors"][-1]["email"] == "xuexin@njucm.edu.cn"


def test_release_documentation_equations_and_data_dictionary_cover_public_outputs():
    required_docs = {
        "methods.md", "oips-formula.md", "data-dictionary.md",
        "provenance.md", "limitations.md", "journal-statements.md",
    }
    assert required_docs <= {path.name for path in (ROOT / "docs").glob("*.md")}
    formula = (ROOT / "docs" / "oips-formula.md").read_text(encoding="utf-8")
    for token in (
        "OIPS-P_static", "C_cons", "G_geo", "P_lig", "Q_evidence",
        "O_rel_formal", "DCC", "IoU", "Top-k", "MRR", "bootstrap",
    ):
        assert token in formula
    assert "\\sum" in formula and "\\operatorname" in formula

    dictionary = (ROOT / "docs" / "data-dictionary.md").read_text(encoding="utf-8")
    public_csvs = {
        path.relative_to(ROOT).as_posix()
        for path in (ROOT / "results" / "reference").rglob("*.csv")
    }
    assert public_csvs
    assert all(f"`{relative}`" in dictionary for relative in public_csvs)
    assert "empty numeric cell" in dictionary.casefold()


def test_ci_runs_constrained_matrix_and_all_release_gates():
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    for version in ('"3.11"', '"3.12"'):
        assert version in workflow
    assert "environment/constraints.txt" in workflow
    for command in (
        "pytest tests/unit", "pytest tests/integration", "pytest tests/scientific",
        "oips-repro validate-data", "oips-repro verify",
        "cffconvert --validate", "scripts/release_check.py",
    ):
        assert command in workflow


def test_release_manifest_exactly_hashes_reference_payload_and_declares_version():
    manifest = json.loads((ROOT / "release" / "manifest.json").read_text(encoding="utf-8"))
    assert set(manifest) == {
        "schema_version", "generated_at_utc", "version", "release_status", "git",
        "release_tag", "identifiers", "configuration", "inputs", "environment",
        "key_results", "publication_payload",
    }
    assert manifest["version"] == "1.0.0"
    assert manifest["release_tag"]["tag"] == "v1.0.0-manuscript"
    declared = {entry["path"]: entry for entry in manifest["publication_payload"]}
    actual = {
        path.relative_to(ROOT).as_posix(): path
        for prefix in (ROOT / "results" / "reference", ROOT / "figures" / "manuscript")
        for path in prefix.rglob("*")
        if path.is_file()
    }
    assert set(declared) == set(actual)
    for relative, path in actual.items():
        content = path.read_bytes()
        assert declared[relative] == {
            "path": relative,
            "sha256": hashlib.sha256(content).hexdigest(),
            "bytes": len(content),
        }


def test_curator_cli_exposes_audit_force_and_explicit_roots():
    completed = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "prepare_public_data.py"), "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0
    assert "--source-root" in completed.stdout
    assert "--destination-root" in completed.stdout
    assert "--rights-policy" in completed.stdout
    assert "--audit-only" in completed.stdout
    assert "--force" in completed.stdout
