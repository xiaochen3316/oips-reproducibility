"""Pure table builders for the curated OIPS public snapshot."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


UNAVAILABLE = "not_available_in_source_materials"
MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null"}

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
SYSTEM_SOURCE_COLUMNS = [
    "pdb_id",
    "Protein_system",
    "Oligomeric_state",
    "Reference_ligand_annotation",
    "Dataset_role",
    "Pocket_category",
    "Category_mapping_status",
    "Protein_family",
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
MD_SOURCE_COLUMNS = [
    "PDB",
    "MD_run",
    "Simulation_context",
    "Persistent_MD_contact_residues",
    "MD_contact_center",
    "D_dyn_run_score",
]
MD_COLUMNS = [
    "pdb_id",
    "MD_run",
    "Simulation_context",
    "Persistent_MD_contact_residues",
    "MD_contact_center",
    "D_dyn_run_score",
]
REDOCKING_SOURCE_COLUMNS = [
    "PDB",
    "Docking_software",
    "GlideScore_kcal_per_mol",
    "Raw_ligand_RMSD_A",
    "RMSD_threshold_call",
    "Failure_or_warning_reason",
]
REDOCKING_COLUMNS = ["pdb_id", *REDOCKING_SOURCE_COLUMNS[1:]]
LEGACY_ALLOWLIST = ["pdb_id", "cluster_id", "member_rows", "member_pockets"]
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


def read_source_csv(path: Path) -> pd.DataFrame:
    """Read source CSV while preserving exact source strings and blank fields."""

    if not path.is_file():
        raise FileNotFoundError(f"required source file is missing: {path}")
    return pd.read_csv(path, keep_default_na=False, low_memory=False, encoding="utf-8")


def require_columns(frame: pd.DataFrame, columns: Iterable[str], source: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{source} is missing allow-listed columns: {missing}")


def normalize_missing_value(value: object) -> object:
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.lower() in MISSING_TOKENS:
            return pd.NA
        return stripped
    return value


def normalize_missing(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.map(normalize_missing_value)


def numeric_column(
    frame: pd.DataFrame, column: str, *, integer: bool = False, nullable: bool = True
) -> None:
    original = frame[column]
    converted = pd.to_numeric(original, errors="coerce")
    invalid = original.notna() & converted.isna()
    if invalid.any():
        examples = original.loc[invalid].astype(str).head(3).tolist()
        raise ValueError(f"non-numeric values in {column}: {examples}")
    if not nullable and converted.isna().any():
        raise ValueError(f"required numeric column {column} contains missing values")
    if integer:
        nonmissing = converted.dropna()
        if not (nonmissing % 1 == 0).all():
            raise ValueError(f"non-integer values in {column}")
        frame[column] = converted.astype("Int64")
    else:
        frame[column] = converted.astype(float)


def fill_required_text(frame: pd.DataFrame, columns: Iterable[str]) -> None:
    for column in columns:
        frame[column] = frame[column].fillna(UNAVAILABLE).astype(str)


def build_features(source: pd.DataFrame) -> pd.DataFrame:
    require_columns(source, FEATURE_COLUMNS, "feature source")
    output = normalize_missing(source.loc[:, FEATURE_COLUMNS].copy())
    output["pdb_id"] = output["pdb_id"].astype(str).str.upper()
    for column in ["row_id", "residue_count"]:
        numeric_column(output, column, integer=True, nullable=False)
    for column in ["display_order", "sitemap_rank"]:
        numeric_column(output, column, integer=True)
    for column in [
        "center_x",
        "center_y",
        "center_z",
        "pocket_geometry_score",
        "pocket_ligandability_score",
    ]:
        numeric_column(output, column)
    fill_required_text(
        output,
        ["pdb_id", "tool", "pocket_id", "center_method", "residue_set_json"],
    )
    if len(output) != 1742 or output["pdb_id"].nunique() != 21:
        raise ValueError("feature source must contain exactly 1,742 rows and 21 targets")
    if not output["row_id"].is_unique:
        raise ValueError("feature row_id values must be unique")
    for value in output["residue_set_json"]:
        parsed = json.loads(value)
        if not isinstance(parsed, list):
            raise ValueError("residue_set_json must encode a JSON list")
    return output


def build_structure_index(source: pd.DataFrame, *, manifest_path: Path) -> pd.DataFrame:
    require_columns(source, ["pdb_id", "paired_structure"], "structure manifest")
    output = normalize_missing(source.loc[:, ["pdb_id", "paired_structure"]].copy())
    fill_required_text(output, ["pdb_id", "paired_structure"])
    output["pdb_id"] = output["pdb_id"].str.upper()
    if len(output) != 21 or not output["pdb_id"].is_unique:
        raise ValueError("structure manifest must contain 21 unique targets")
    for row in output.itertuples(index=False):
        relative = Path(row.paired_structure)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("paired_structure must be a source-root-relative path")
        if not (manifest_path.parent / relative).is_file():
            raise FileNotFoundError(f"prepared structure is missing for {row.pdb_id}")
    return output


def build_systems(category_source: pd.DataFrame, structures: pd.DataFrame) -> pd.DataFrame:
    require_columns(category_source, SYSTEM_SOURCE_COLUMNS, "category source")
    output = normalize_missing(category_source.loc[:, SYSTEM_SOURCE_COLUMNS].copy())
    fill_required_text(output, SYSTEM_SOURCE_COLUMNS)
    output = output.rename(
        columns={
            "Protein_system": "protein_system",
            "Oligomeric_state": "oligomeric_state",
            "Reference_ligand_annotation": "reference_ligand_annotation",
            "Dataset_role": "dataset_role",
            "Pocket_category": "pocket_category",
            "Category_mapping_status": "category_mapping_status",
            "Protein_family": "protein_family",
        }
    )
    output["pdb_id"] = output["pdb_id"].str.upper()
    if set(output["pdb_id"]) != set(structures["pdb_id"]):
        raise ValueError("category and structure manifests must contain the same targets")
    output["structure_source"] = "RCSB PDB-derived prepared structure"
    for column in SYSTEM_COLUMNS[9:28]:
        output[column] = UNAVAILABLE
    output["metadata_completeness"] = "approved_category_fields_and_structure_provenance_only"
    output["metadata_status"] = "author_confirmation_required"
    return output.loc[:, SYSTEM_COLUMNS]


def build_reference_annotations(source: pd.DataFrame) -> pd.DataFrame:
    required = ["pdb_id", "selected_ligand", "ligand_status"]
    require_columns(source, required, "target summary")
    output = normalize_missing(source.loc[:, required].copy())
    fill_required_text(output, required)
    output["pdb_id"] = output["pdb_id"].str.upper()
    output["reference_ligand_annotation"] = output["selected_ligand"]
    output["override_applied"] = False
    output["decision_id"] = "not_applicable"
    reviewed = output["pdb_id"] == "5LGE"
    if reviewed.sum() != 1:
        raise ValueError("target summary must contain one 5LGE row")
    output.loc[reviewed, "reference_ligand_annotation"] = "6VN:D:503:"
    output.loc[reviewed, "override_applied"] = True
    output.loc[reviewed, "decision_id"] = "DEC-REF-5LGE-001"
    if len(output) != 21 or not output["pdb_id"].is_unique:
        raise ValueError("reference source must contain 21 unique targets")
    return output.loc[:, REFERENCE_COLUMNS]


def build_md_evidence(source: pd.DataFrame) -> pd.DataFrame:
    require_columns(source, MD_SOURCE_COLUMNS, "MD source")
    output = normalize_missing(source.loc[:, MD_SOURCE_COLUMNS].copy())
    output = output.rename(columns={"PDB": "pdb_id"})
    fill_required_text(output, ["pdb_id", "MD_run", "Simulation_context"])
    output["pdb_id"] = output["pdb_id"].str.upper()
    numeric_column(output, "D_dyn_run_score")
    if output.duplicated(["pdb_id", "Simulation_context", "MD_run"]).any():
        raise ValueError("MD evidence composite keys must be unique")
    return output.loc[:, MD_COLUMNS]


def build_redocking_evidence(source: pd.DataFrame) -> pd.DataFrame:
    require_columns(source, REDOCKING_SOURCE_COLUMNS, "redocking source")
    output = normalize_missing(source.loc[:, REDOCKING_SOURCE_COLUMNS].copy())
    output = output.rename(columns={"PDB": "pdb_id"})
    output["pdb_id"] = output["pdb_id"].astype(str).str.upper()
    for column in ["GlideScore_kcal_per_mol", "Raw_ligand_RMSD_A"]:
        numeric_column(output, column)
    fill_required_text(
        output,
        ["pdb_id", "Docking_software", "RMSD_threshold_call", "Failure_or_warning_reason"],
    )
    if len(output) != 21 or not output["pdb_id"].is_unique:
        raise ValueError("redocking source must contain 21 unique targets")
    return output.loc[:, REDOCKING_COLUMNS]


def build_legacy_crosswalk(
    legacy_source: pd.DataFrame, mapping_source: pd.DataFrame
) -> pd.DataFrame:
    require_columns(legacy_source, LEGACY_ALLOWLIST, "legacy score source")
    require_columns(mapping_source, LEGACY_COLUMNS, "formal crosswalk")
    legacy = normalize_missing(legacy_source.loc[:, LEGACY_ALLOWLIST].copy())
    mapping = normalize_missing(mapping_source.loc[:, LEGACY_COLUMNS].copy())
    legacy_keys = set(zip(legacy["pdb_id"], legacy["cluster_id"]))
    mapping_keys = set(zip(mapping["pdb_id"], mapping["working_cluster_id"]))
    if mapping_keys != legacy_keys:
        raise ValueError("formal crosswalk IDs do not exactly match the legacy allowlist")
    for column in [
        "mapped_record_count",
        "working_cluster_member_count",
        "unmappable_record_count",
    ]:
        numeric_column(mapping, column, integer=True, nullable=False)
    numeric_column(mapping, "mapping_fraction", nullable=False)
    mapping["pdb_id"] = mapping["pdb_id"].astype(str).str.upper()
    fill_required_text(mapping, ["pdb_id", "working_cluster_id", "mapping_type"])
    if len(mapping) != 911:
        raise ValueError("formal crosswalk must contain 911 rows")
    return mapping.loc[:, LEGACY_COLUMNS]


def build_manual_decisions(
    category_source: pd.DataFrame,
    destination_root: Path,
    representative_source: pd.DataFrame,
) -> pd.DataFrame:
    require_columns(
        category_source,
        [
            "pdb_id",
            "Initial_pocket_annotation",
            "Pocket_category",
            "Category_mapping_status",
        ],
        "category source",
    )
    require_columns(representative_source, ["pdb_id"], "representative cases")
    representative_ids = set(representative_source["pdb_id"].astype(str).str.upper())
    if representative_ids != {"5J89", "5TBM", "4W9H"}:
        raise ValueError("representative-case source must contain exactly the approved cases")

    rows: list[dict[str, str]] = []
    for source_row in category_source.sort_values("pdb_id").to_dict("records"):
        pdb_id = str(source_row["pdb_id"]).upper()
        rows.append(
            {
                "decision_id": f"DEC-CATEGORY-{pdb_id}",
                "pdb_id_or_scope": pdb_id,
                "decision_type": "category_annotation",
                "field_name": "pocket_category",
                "previous_value": str(source_row["Initial_pocket_annotation"]),
                "released_value": str(source_row["Pocket_category"]),
                "rationale": (
                    "Reviewed project category retained with mapping status "
                    f"{source_row['Category_mapping_status']}."
                ),
                "evidence_source": (
                    "OIPS_FINAL_ANALYSIS/00_configuration/"
                    f"pocket_category_mapping.csv#pdb_id={pdb_id}"
                ),
                "approval_status": "approved_for_repository",
            }
        )
    rows.append(
        {
            "decision_id": "DEC-REF-5LGE-001",
            "pdb_id_or_scope": "5LGE",
            "decision_type": "reference_key_correction",
            "field_name": "reference_ligand_annotation",
            "previous_value": "NAP:C:501:",
            "released_value": "6VN:D:503:",
            "rationale": (
                "Reviewed correction selects the redocked 6VN ligand at chain D residue 503; "
                "the prepared structure contains this exact residue key."
            ),
            "evidence_source": (
                "zhu-pack/pack/pdb-resource/5LGE/5LGE-p.pdb#6VN:D:503:; "
                "OIPS_FINAL_ANALYSIS/07_redocking/redocking_complete_summary.csv#PDB=5LGE"
            ),
            "approval_status": "approved_for_repository",
        }
    )
    config_path = destination_root / "config" / "manuscript.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    parameter_set = json.dumps(
        config["clustering"], ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    rows.append(
        {
            "decision_id": "DEC-PARAM-CLUSTER-V2-001",
            "pdb_id_or_scope": "all_targets",
            "decision_type": "parameter_override",
            "field_name": "cluster_v2_parameter_set",
            "previous_value": "not_formalized_as_release_contract",
            "released_value": parameter_set,
            "rationale": "Freeze the reviewed same-tool and cross-tool cluster-v2 parameters for deterministic reproduction.",
            "evidence_source": "config/manuscript.yaml#clustering",
            "approval_status": "approved_for_repository",
        }
    )
    cases = {
        "5J89": (
            "supported_interface_site_with_boundary_uncertainty",
            "V2C003 is static/reference rank 1 with R_auto, four-tool support, interface overlap 0.369, and redocking RMSD 0.368 A; it supports site prioritization while retaining exact-boundary uncertainty.",
        ),
        "5TBM": (
            "competing_assembly_and_internal_allosteric_hypotheses",
            "V2C001 is static rank 1 A_auto while established reference V2C002 is rank 2 and becomes rank 1 without O_rel, with RMSD 0.819 A; V2C001 is not claimed as a new allosteric site.",
        ),
        "4W9H": (
            "unresolved_cross_layer_conflict",
            "Assembly-associated V2C003 is static rank 1 while reference V2C001 is rank 10, MD maps locally, and RMSD is 3.544 A; discordance is preserved without a majority-vote biological conclusion.",
        ),
    }
    for pdb_id, (label, rationale) in cases.items():
        rows.append(
            {
                "decision_id": f"DEC-CASE-{pdb_id}",
                "pdb_id_or_scope": pdb_id,
                "decision_type": "representative_case_selection",
                "field_name": "representative_case_label",
                "previous_value": "not_selected",
                "released_value": label,
                "rationale": rationale,
                "evidence_source": (
                    "OIPS_FINAL_ANALYSIS_0711/representative_case_results.csv#"
                    f"pdb_id={pdb_id}; OIPS_FINAL_ANALYSIS_0711/claim_evidence_matrix.csv"
                ),
                "approval_status": "approved_for_repository",
            }
        )
    output = pd.DataFrame(rows, columns=DECISION_COLUMNS)
    if len(output) != 26 or not output["decision_id"].is_unique:
        raise ValueError("manual decision audit must contain 26 stable decisions")
    return output
