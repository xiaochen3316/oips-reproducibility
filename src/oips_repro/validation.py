"""Aggregating validation primitives for OIPS input and result bundles."""
from __future__ import annotations
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
import math
from pathlib import Path, PurePosixPath
import json
import re
import pandas as pd
from jsonschema import Draft202012Validator, FormatChecker
from . import io, posthoc, scoring, statistics
VALID_STATUSES = frozenset({"pass", "fail", "warning"})
@dataclass(frozen=True)
class CheckResult:
    check_id: str
    status: str
    summary: str
    evidence: Mapping[str, object]
    def __post_init__(self) -> None:
        if not isinstance(self.check_id, str) or not self.check_id:
            raise ValueError("check_id must be a nonempty string")
        if self.status not in VALID_STATUSES:
            raise ValueError("status must be exactly pass, fail, or warning")
        if not isinstance(self.summary, str) or not self.summary:
            raise ValueError("summary must be a nonempty string")
        if not isinstance(self.evidence, Mapping):
            raise ValueError("evidence must be a mapping")
        object.__setattr__(self, "evidence", MappingProxyType(dict(self.evidence)))
@dataclass(frozen=True)
class ValidationReport:
    checks: Sequence[CheckResult]
    def __post_init__(self) -> None:
        values = tuple(self.checks)
        if any(not isinstance(value, CheckResult) for value in values):
            raise ValueError("checks must contain CheckResult values")
        identifiers = [value.check_id for value in values]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("validation check IDs must be unique")
        object.__setattr__(self, "checks", values)
    @property
    def pass_count(self) -> int:
        return sum(check.status == "pass" for check in self.checks)
    @property
    def failure_count(self) -> int:
        return sum(check.status == "fail" for check in self.checks)
    @property
    def warning_count(self) -> int:
        return sum(check.status == "warning" for check in self.checks)
    @property
    def status(self) -> str:
        return "fail" if self.failure_count else "pass"
    def raise_for_failures(self) -> None:
        failed = [check.check_id for check in self.checks if check.status == "fail"]
        if failed:
            raise ValueError("validation failed checks: " + ", ".join(failed))
POSTHOC_INPUT_FILES = {
    "reference_annotations.csv", "md_evidence.csv", "redocking_evidence.csv",
}
ANALYSIS_SCHEMAS = {
    "final_top3_automated_QC.csv": posthoc.TOPK_QC_COLUMNS,
    "final_reference_mapping.csv": posthoc.REFERENCE_COLUMNS,
    "final_md_cluster_v2_mapping.csv": posthoc.MD_COLUMNS,
    "final_automated_evidence_labels.csv": posthoc.LABEL_COLUMNS,
    "final_redocking_cluster_v2_mapping.csv": posthoc.REDOCKING_COLUMNS,
    "unresolved_cases.csv": posthoc.UNRESOLVED_COLUMNS,
    "final_cluster_v2_master_table.csv": (*io.RANKING_COLUMNS, *posthoc.CONVENIENCE_COLUMNS),
    "target_level_candidate_prioritization.csv": statistics.TARGET_COLUMNS,
    "final_candidate_prioritization_metrics.csv": statistics.POINT_COLUMNS,
    "final_category_metrics.csv": statistics.POINT_COLUMNS,
    "final_bootstrap_intervals.csv": statistics.BOOTSTRAP_COLUMNS,
    "final_family_sensitivity.csv": statistics.FAMILY_COLUMNS,
    "final_orel_ablation_targets.csv": statistics.ABLATION_TARGET_COLUMNS,
    "final_orel_ablation_categories.csv": statistics.ABLATION_CATEGORY_COLUMNS,
    "orel_ablation_cluster_rankings.csv": statistics.ABLATION_CLUSTER_COLUMNS,
    "representative_case_results.csv": statistics.REPRESENTATIVE_COLUMNS,
    "weight_sensitivity_scenarios.csv": statistics.WEIGHT_SCENARIO_COLUMNS,
    "weight_sensitivity_targets.csv": statistics.WEIGHT_TARGET_COLUMNS,
    "single_tool_target_ranks.csv": statistics.SINGLE_TOOL_TARGET_COLUMNS,
    "single_tool_complete_case_metrics.csv": statistics.SINGLE_TOOL_METRIC_COLUMNS,
}
ANALYSIS_SORTS = {
    "final_top3_automated_QC.csv": ("pdb_id", "Within_PDB_rank", "cluster_v2_id"),
    "final_reference_mapping.csv": ("pdb_id", "Within_PDB_rank", "cluster_v2_id"),
    "final_md_cluster_v2_mapping.csv": ("pdb_id", "Simulation_context", "MD_run"),
    "final_automated_evidence_labels.csv": ("pdb_id", "Within_PDB_rank", "cluster_v2_id"),
    "final_redocking_cluster_v2_mapping.csv": ("pdb_id",),
    "unresolved_cases.csv": ("issue_id",),
    "final_cluster_v2_master_table.csv": ("pdb_id", "Within_PDB_rank", "cluster_v2_id"),
    "target_level_candidate_prioritization.csv": ("pdb_id",),
    "final_candidate_prioritization_metrics.csv": ("analysis_level", "group"),
    "final_category_metrics.csv": ("group",),
    "final_family_sensitivity.csv": ("excluded_family",),
    "final_orel_ablation_targets.csv": ("pdb_id",),
    "final_orel_ablation_categories.csv": ("Pocket_category",),
    "orel_ablation_cluster_rankings.csv": ("pdb_id", "Without_O_rel_rank", "cluster_v2_id"),
    "weight_sensitivity_scenarios.csv": ("perturbed_module", "direction"),
    "weight_sensitivity_targets.csv": ("pdb_id",),
    "single_tool_target_ranks.csv": ("pdb_id", "method"),
    "single_tool_complete_case_metrics.csv": ("method",),
}
ANALYSIS_BOOLEAN_COLUMNS = {
    "final_top3_automated_QC.csv": (
        "center_available", "residue_set_available", "spatial_continuity",
        "possible_over_merging", "clear_exclusion_flag", "possible_split_subpocket",
        "possible_surface_noise_cluster", "missing_residue_proximity_flag", "boundary_sensitive",
    ),
    "final_reference_mapping.csv": ("R_auto_rule_pass", "reference_selection_unresolved"),
    "final_automated_evidence_labels.csv": ("unresolved_flag",),
    "final_redocking_cluster_v2_mapping.csv": ("Reference_pose_recovered",),
    "final_cluster_v2_master_table.csv": (
        "spatial_continuity", "boundary_sensitive", "tie_flag", "unresolved_flag", "R_auto_rule_pass",
    ),
    "target_level_candidate_prioritization.csv": (
        "reference_evaluable", "reference_top1", "reference_top3", "reference_top5",
        "first_supported_evaluable", "first_supported_top1", "first_supported_top3",
    ),
    "final_orel_ablation_targets.csv": (
        "Full_reference_Top3", "Without_O_rel_reference_Top3", "top_cluster_identity_changed",
    ),
    "single_tool_target_ranks.csv": ("complete_case",),
}
ANALYSIS_NUMERIC_COLUMNS = {
    "final_top3_automated_QC.csv": (
        "Within_PDB_rank", "OIPS-P_static", "cluster_diameter_A", "tool_support_count",
        "center_dispersion_A", "core_envelope_ratio", "contributing_chain_count",
        "dominant_chain_fraction", "pocket_interface_overlap", "distance_to_interface_A",
        "nearest_missing_residue_distance_sequence_positions",
    ),
    "final_reference_mapping.csv": (
        "Within_PDB_rank", "DCC_A", "reference_contact_residue_count",
        "reference_contact_overlap_count", "contact_precision", "contact_recall", "residue_IoU",
        "functional_residue_overlap", "interface_overlap", "distance_to_interface_A",
        "chain_contribution_count", "reference_ligand_atom_count",
    ),
    "final_md_cluster_v2_mapping.csv": (
        "static_top_rank", "persistent_MD_contact_residue_count", "best_MD_cluster_Jaccard",
        "best_MD_cluster_precision", "best_MD_cluster_recall", "best_MD_cluster_center_distance_A",
        "static_top_MD_overlap_count", "static_top_MD_cluster_coverage",
        "static_top_MD_contact_coverage", "static_top_MD_Jaccard",
        "static_top_MD_center_distance_A", "D_dyn_run_score",
    ),
    "final_automated_evidence_labels.csv": (
        "Within_PDB_rank", "R_auto_DCC_A", "R_auto_contact_recall", "R_auto_residue_IoU",
    ),
    "final_redocking_cluster_v2_mapping.csv": (
        "reference_cluster_static_rank", "Raw_ligand_RMSD_A", "GlideScore_kcal_per_mol",
    ),
    "final_cluster_v2_master_table.csv": tuple(
        column for column in (*io.RANKING_COLUMNS, *posthoc.CONVENIENCE_COLUMNS)
        if column in io._RESULT_FLOAT_COLUMNS or column in io._RESULT_INTEGER_COLUMNS
        or column in {"DCC_A", "contact_recall", "contact_precision", "residue_IoU"}
    ),
    "target_level_candidate_prioritization.csv": (
        "formal_cluster_v2_count", "R_auto_cluster_count", "reference_first_rank",
        "reference_rank_percentile", "reference_reciprocal_rank", "first_supported_rank",
        "first_supported_reciprocal_rank",
    ),
    "final_candidate_prioritization_metrics.csv": statistics.METRIC_FIELDS,
    "final_category_metrics.csv": statistics.METRIC_FIELDS,
    "final_bootstrap_intervals.csv": (
        "point_estimate", "CI_2.5_percent", "CI_97.5_percent", "iterations", "random_seed",
    ),
    "final_family_sensitivity.csv": ("excluded_target_N", "remaining_target_N", *statistics.METRIC_FIELDS),
    "final_orel_ablation_targets.csv": (
        "Full_reference_rank", "Without_O_rel_reference_rank",
        "reference_rank_change_without_minus_full", "Full_first_supported_rank",
        "Without_O_rel_first_supported_rank", "first_supported_rank_change_without_minus_full",
        "reference_Top3_inclusion_change", "rank_Spearman_rho",
    ),
    "final_orel_ablation_categories.csv": (
        "target_N", "Full_reference_Top1", "Without_O_rel_reference_Top1",
        "Full_reference_Top3", "Without_O_rel_reference_Top3", "Full_reference_MRR",
        "Without_O_rel_reference_MRR", "mean_rank_Spearman_rho",
        "top_cluster_identity_change_fraction",
    ),
    "orel_ablation_cluster_rankings.csv": ("Within_PDB_rank", "Without_O_rel_score", "Without_O_rel_rank"),
    "representative_case_results.csv": (
        "static_top_score", "static_top_tool_support", "static_top_interface_overlap",
        "reference_rank", "reference_DCC_A", "without_O_rel_reference_rank", "redocking_RMSD_A",
    ),
    "weight_sensitivity_scenarios.csv": (
        "multiplier", "target_N", "baseline_top1_retained_N", "mean_top3_jaccard",
        "median_spearman_rho", "Reference_Top1_N", "Reference_Top3_N",
        "First_supported_Top1_N", "First_supported_Top3_N",
    ),
    "weight_sensitivity_targets.csv": (
        "top1_retention_count", "baseline_reference_rank", "minimum_reference_rank",
        "maximum_reference_rank", "baseline_first_supported_rank",
        "minimum_first_supported_rank", "maximum_first_supported_rank",
    ),
    "single_tool_target_ranks.csv": (
        "mappable_region_count", "first_reference_associated_rank",
    ),
    "single_tool_complete_case_metrics.csv": (
        "N", "Top1_N", "Top1", "Top3_N", "Top3", "Top5_N", "Top5", "MRR",
    ),
}
ANALYSIS_ENUM_COLUMNS = {
    "final_top3_automated_QC.csv": {
        "QC_status": {"QC_pass", "QC_boundary_sensitive", "QC_possible_split", "QC_possible_overmerge", "QC_insufficient_evidence", "QC_unmappable"},
        "mappability": set(posthoc.VALID_MAPPABILITY),
    },
    "final_reference_mapping.csv": {"reference_selection_status": {"exact_project_reference_key", "unique_resname_fallback", "ambiguous_resname_largest_group_fallback", "invalid_project_reference_key", "reference_ligand_not_found"}},
    "final_md_cluster_v2_mapping.csv": {"Concordance_call": {"concordant", "partially_concordant", "boundary_shift", "static_dynamic_conflict", "apo_only_context", "insufficient_MD_evidence", "MD_not_available"}},
    "final_automated_evidence_labels.csv": {"automated_evidence_label": {"R_auto", "A_auto", "U_auto", "X_auto"}},
    "final_redocking_cluster_v2_mapping.csv": {"RMSD_threshold_call": {"RMSD <= 2 A", "2 A < RMSD <= 3 A", "RMSD > 3 A"}},
    "final_cluster_v2_master_table.csv": {
        "mappability": set(posthoc.VALID_MAPPABILITY),
        "automated_evidence_label": {"R_auto", "A_auto", "U_auto", "X_auto"},
    },
    "target_level_candidate_prioritization.csv": {"static_top_evidence_label": {"R_auto", "A_auto", "U_auto", "X_auto"}},
    "weight_sensitivity_scenarios.csv": {
        "perturbed_module": set(statistics.STATIC_MODULES),
        "direction": {"decrease", "increase"},
    },
    "single_tool_target_ranks.csv": {
        "method": {"CavityPlus", "DoGSiteScorer", "DoGSite3", "CASTpFold", "SiteMap", "OIPS-P"},
        "output_status": {"unavailable", "unmappable", "no_hit", "reference_hit"},
    },
    "single_tool_complete_case_metrics.csv": {
        "method": {"CavityPlus", "DoGSiteScorer", "DoGSite3", "CASTpFold", "SiteMap", "OIPS-P"},
    },
    "final_bootstrap_intervals.csv": {
        "bootstrap_method": {"target_resampling", "family_clustered"},
        "metric": set(statistics.BOOTSTRAP_METRICS),
    },
}
ANALYSIS_KEYS = {
    "final_top3_automated_QC.csv": ("pdb_id", "cluster_v2_id"),
    "final_reference_mapping.csv": ("pdb_id", "cluster_v2_id"),
    "final_md_cluster_v2_mapping.csv": ("pdb_id", "Simulation_context", "MD_run"),
    "final_automated_evidence_labels.csv": ("pdb_id", "cluster_v2_id"),
    "final_redocking_cluster_v2_mapping.csv": ("pdb_id",),
    "unresolved_cases.csv": ("issue_id",),
    "final_cluster_v2_master_table.csv": ("pdb_id", "cluster_v2_id"),
    "target_level_candidate_prioritization.csv": ("pdb_id",),
    "final_candidate_prioritization_metrics.csv": ("analysis_level", "group"),
    "final_category_metrics.csv": ("analysis_level", "group"),
    "final_bootstrap_intervals.csv": ("bootstrap_method", "metric"),
    "final_family_sensitivity.csv": ("excluded_family",),
    "final_orel_ablation_targets.csv": ("pdb_id",),
    "final_orel_ablation_categories.csv": ("Pocket_category",),
    "orel_ablation_cluster_rankings.csv": ("pdb_id", "cluster_v2_id"),
    "representative_case_results.csv": ("pdb_id",),
    "weight_sensitivity_scenarios.csv": ("scenario_id",),
    "weight_sensitivity_targets.csv": ("pdb_id",),
    "single_tool_target_ranks.csv": ("pdb_id", "method"),
    "single_tool_complete_case_metrics.csv": ("method",),
}
def validate_analysis_csv(filename: str, raw: pd.DataFrame) -> pd.DataFrame:
    columns = ANALYSIS_SCHEMAS[filename]
    if raw.columns.tolist() != list(columns):
        raise ValueError(f"{filename} has incorrect schema: {raw.columns.tolist()}")
    text, parsed = raw.astype(str).copy(), raw.astype(str).copy()
    for column in ANALYSIS_BOOLEAN_COLUMNS.get(filename, ()):
        invalid = ~text[column].isin(("", io.BOOL_TRUE, io.BOOL_FALSE))
        if invalid.any():
            raise ValueError(f"{filename} boolean column {column} must use lowercase true/false")
        parsed[column] = text[column].map({io.BOOL_TRUE: True, io.BOOL_FALSE: False, "": None})
    for column in ANALYSIS_NUMERIC_COLUMNS.get(filename, ()):
        missing = text[column].eq("")
        numeric = pd.to_numeric(text[column].mask(missing), errors="coerce")
        invalid = ~missing & (numeric.isna() | ~numeric.map(lambda value: bool(pd.isna(value)) or math.isfinite(float(value))))
        if invalid.any():
            raise ValueError(f"{filename} numeric column {column} must be finite or missing")
        parsed[column] = numeric
    for column, allowed in ANALYSIS_ENUM_COLUMNS.get(filename, {}).items():
        if (~text[column].isin(allowed)).any():
            raise ValueError(f"{filename} enum column {column} contains an invalid value")
    keys = ANALYSIS_KEYS[filename]
    if parsed.duplicated(list(keys)).any():
        raise ValueError(f"{filename} contains duplicate keys: {list(keys)}")
    if filename in ANALYSIS_SORTS:
        expected = parsed.sort_values(list(ANALYSIS_SORTS[filename]), kind="mergesort", na_position="last").reset_index(drop=True)
        if not parsed.reset_index(drop=True).equals(expected):
            raise ValueError(f"{filename} is not sorted by {list(ANALYSIS_SORTS[filename])}")
    return parsed
def _safe_repo_path(value: str, root: Path, label: str) -> Path:
    if not value or "\\" in value:
        raise ValueError(f"{label} must be a relative POSIX path")
    posix = PurePosixPath(value)
    if posix.is_absolute() or ".." in posix.parts or re.match(r"^[A-Za-z]:", value):
        raise ValueError(f"{label} must be a contained relative POSIX path")
    candidate = root.resolve()
    for part in posix.parts:
        candidate = candidate / part
        if candidate.exists() and _link_or_reparse(candidate):
            raise ValueError(f"{label} may not traverse a symlink or junction")
    path = (root / Path(*posix.parts)).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError(f"{label} escapes the repository") from error
    return path
def _link_or_reparse(path: Path) -> bool:
    try:
        stat = path.lstat()
    except OSError:
        return True
    return path.is_symlink() or bool(getattr(stat, "st_file_attributes", 0) & 0x400)


def _checksums(config) -> dict[str, tuple[str, int]]:
    from .provenance import sha256_file

    source = config.resolve_configured_path("data_checksums")
    rows: dict[str, tuple[str, int]] = {}
    previous = ""
    for number, line in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (\S+)", line)
        if not match:
            raise ValueError(f"checksum line {number} is malformed")
        digest, relative = match.groups()
        if relative <= previous or relative in rows:
            raise ValueError("checksum entries must be unique and sorted")
        previous = relative
        path = _safe_repo_path(relative, config.repository_root, "checksum path")
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"checksum path is not a regular file: {relative}")
        if sha256_file(path) != digest:
            raise ValueError(f"checksum mismatch: {relative}")
        rows[relative] = (digest, path.stat().st_size)
    data_root = source.parent
    entries = list(data_root.rglob("*"))
    if any(_link_or_reparse(path) for path in entries):
        raise ValueError("data payload may not contain symlinks or junctions")
    actual = {path.relative_to(config.repository_root).as_posix() for path in entries if path.is_file() and path != source}
    if set(rows) != actual:
        raise ValueError("checksum inventory must exactly cover the data payload")
    return rows


def _coerce_schema(value: str, schema: Mapping[str, object], label: str) -> object:
    types = schema.get("type")
    allowed = set(types if isinstance(types, list) else [types])
    if value == "" and "null" in allowed:
        return None
    if "boolean" in allowed:
        if value not in {"true", "false"}:
            raise ValueError(f"{label} must use lowercase true/false")
        return value == "true"
    if "integer" in allowed:
        if not re.fullmatch(r"-?\d+", value):
            raise ValueError(f"{label} must be an integer")
        return int(value)
    if "number" in allowed:
        try:
            result = float(value)
        except ValueError as error:
            raise ValueError(f"{label} must be numeric") from error
        if not math.isfinite(result):
            raise ValueError(f"{label} must be finite")
        return result
    one_of = schema.get("oneOf")
    if isinstance(one_of, list):
        if re.fullmatch(r"-?\d+", value) and any(item.get("type") == "integer" for item in one_of if isinstance(item, Mapping)):
            return int(value)
    return value


def _validate_schema_table(
    path: Path, definition: Mapping[str, object], *, delimiter: str,
) -> pd.DataFrame:
    columns = definition.get("x-columns")
    properties = definition.get("properties")
    if not isinstance(columns, list) or not isinstance(properties, Mapping):
        raise ValueError("table schema definition is incomplete")
    raw = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False)
    if raw.columns.tolist() != columns:
        raise ValueError(f"{path.name} has incorrect exact schema")
    records: list[dict[str, object]] = []
    for row_number, record in enumerate(raw.to_dict("records"), 2):
        records.append({
            column: _coerce_schema(str(record[column]), properties[column], f"{path.name}:{row_number}:{column}")
            for column in columns
        })
    validator = Draft202012Validator(definition, format_checker=FormatChecker())
    errors = [error for record in records for error in validator.iter_errors(record)]
    if errors:
        raise ValueError(f"{path.name} schema validation failed: {errors[0].message}")
    keys = definition.get("x-primary-key", ())
    if keys and raw.duplicated(list(keys)).any():
        raise ValueError(f"{path.name} contains duplicate primary keys")
    sort_keys = definition.get("x-sort-key", ())
    if sort_keys:
        parsed = pd.DataFrame(records, columns=columns)
        expected = parsed.sort_values(list(sort_keys), kind="mergesort").reset_index(drop=True)
        if not parsed.reset_index(drop=True).equals(expected):
            raise ValueError(f"{path.name} is not sorted by its declared key")
    return raw


def _validate_release_tables(config) -> dict[str, pd.DataFrame]:
    schema = json.loads(config.resolve_configured_path("schema").read_text(encoding="utf-8"))
    definitions = schema.get("$defs")
    if not isinstance(definitions, Mapping):
        raise ValueError("configured JSON schema has no $defs mapping")
    specs = {
        "feature_table": ("tool_pocket_features_row", ","),
        "systems": ("systems_row", "\t"),
        "reference_annotations": ("reference_annotations_row", ","),
        "md_evidence": ("md_evidence_row", ","),
        "redocking_evidence": ("redocking_evidence_row", ","),
        "data_manifest": ("manifest_row", "\t"),
        "asset_rights": ("asset_rights_row", "\t"),
        "external_archive_manifest": ("external_archive_manifest_row", "\t"),
    }
    tables: dict[str, pd.DataFrame] = {}
    for path_key, (definition_key, delimiter) in specs.items():
        definition = definitions.get(definition_key)
        if not isinstance(definition, Mapping):
            raise ValueError(f"schema definition is missing: {definition_key}")
        tables[path_key] = _validate_schema_table(
            config.resolve_configured_path(path_key), definition, delimiter=delimiter,
        )
    decisions = config.repository_root / "data" / "metadata" / "manual_decisions.tsv"
    if decisions.is_file():
        definition = definitions.get("manual_decisions_row")
        if not isinstance(definition, Mapping):
            raise ValueError("schema definition is missing: manual_decisions_row")
        tables["manual_decisions"] = _validate_schema_table(
            decisions, definition, delimiter="\t",
        )
    return tables


def _validate_manifest(config, checksums: Mapping[str, tuple[str, int]], table: pd.DataFrame) -> None:
    from .provenance import sha256_file

    paths = set(table["local_path_or_url"])
    manifest_relative = config.resolve_configured_path("data_manifest").relative_to(config.repository_root).as_posix()
    checksum_relative = config.resolve_configured_path("data_checksums").relative_to(config.repository_root).as_posix()
    expected = set(checksums) - {manifest_relative}
    expected.discard(checksum_relative)
    if paths != expected:
        raise ValueError("data manifest must exactly declare payload files except itself and checksums")
    for row in table.to_dict("records"):
        path = _safe_repo_path(row["local_path_or_url"], config.repository_root, "manifest path")
        if sha256_file(path) != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise ValueError(f"manifest hash or byte count mismatch: {row['local_path_or_url']}")


def _coverage(config, tables: Mapping[str, pd.DataFrame]) -> None:
    features = io.load_feature_table(config.resolve_configured_path("feature_table"))
    targets = set(features["pdb_id"])
    for key in ("systems", "reference_annotations", "redocking_evidence"):
        values = tables[key]["pdb_id"].astype(str)
        if len(values) != len(targets) or values.nunique() != len(targets) or set(values) != targets:
            raise ValueError(f"{key} must contain one row per feature target")
    md = tables["md_evidence"]
    if not set(md["pdb_id"].astype(str)).issubset(targets):
        raise ValueError("MD evidence targets must be a subset of feature targets")
    io.load_structure_paths(config, target_ids=targets)


def _result(check_id: str, action) -> CheckResult:
    try:
        evidence = action()
        return CheckResult(check_id, "pass", "validated", evidence or {})
    except Exception as error:  # aggregate independent validation failures
        return CheckResult(check_id, "fail", type(error).__name__, {"reason": _safe_reason(error)})


def _safe_reason(error: BaseException) -> str:
    message = str(error).replace("\r", " ").replace("\n", " ")
    message = re.sub(r"(?i)[A-Z]:[\\/].*", "<redacted-path>", message)
    message = re.sub(r"(?i)(?:file|https?)://\S+", "<redacted-url>", message)
    message = re.sub(r"(?i)(?<![:\w])/(?:home|users|private|tmp|var|opt|mnt)/\S*", "<redacted-path>", message)
    return re.sub(r"(?i)\b(?:token|password|secret|api[_-]?key)\s*[:=]\s*\S+", "<redacted-credential>", message)


def validate_public_inputs(config) -> ValidationReport:
    """Collect independent general input failures without stopping at the first."""
    cache: dict[str, object] = {}

    def checksums_action():
        cache["checksums"] = _checksums(config)
        return {"files": len(cache["checksums"])}

    def tables_action():
        cache["tables"] = _validate_release_tables(config)
        return {"tables": len(cache["tables"])}

    def manifest_action():
        checksums = cache.get("checksums") or _checksums(config)
        tables = cache.get("tables") or _validate_release_tables(config)
        _validate_manifest(config, checksums, tables["data_manifest"])
        return {"rows": len(tables["data_manifest"])}

    def coverage_action():
        tables = cache.get("tables") or _validate_release_tables(config)
        _coverage(config, tables)
        features = io.load_feature_table(config.resolve_configured_path("feature_table"))
        return {"records": len(features), "targets": features["pdb_id"].nunique()}

    def snapshot_action():
        value = json.loads(config.resolve_configured_path("expected_summary").read_text(encoding="utf-8"))
        if not isinstance(value, dict):
            raise ValueError("expected summary must be a JSON object")
        return {"keys": len(value)}

    def sensitive_action():
        from .provenance import scan_sensitive_content

        named: dict[str, str] = {}
        for relative in sorted((_checksums(config))):
            path = _safe_repo_path(relative, config.repository_root, "scan path")
            try:
                named[relative] = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
        result = scan_sensitive_content(named)
        if result.status == "fail":
            raise ValueError("sensitive content classes: " + ", ".join(result.evidence["classes"]))
        return result.evidence

    return ValidationReport((
        _result("checksums", checksums_action),
        _result("release_tables", tables_action),
        _result("data_manifest", manifest_action),
        _result("scientific_coverage", coverage_action),
        _result("expected_summary", snapshot_action),
        _result("sensitive_content", sensitive_action),
    ))


def consumed_input_paths(config) -> list[Path]:
    """Return only files actually read by validation and scientific stages."""
    paths = {
        config.resolve_configured_path(key)
        for key in (
            "schema", "data_manifest", "data_checksums", "asset_rights",
            "external_archive_manifest", "expected_summary", "feature_table", "systems",
            "reference_annotations", "md_evidence", "redocking_evidence", "figure_contract",
        )
    }
    for relative in _checksums(config):
        paths.add(_safe_repo_path(relative, config.repository_root, "consumed input"))
    return sorted(paths, key=lambda path: path.relative_to(config.repository_root).as_posix())


def _assert_frame_match(
    left: pd.DataFrame, right: pd.DataFrame, label: str, *, tolerance: float = 1e-12,
) -> None:
    try:
        pd.testing.assert_frame_equal(
            left.reset_index(drop=True), right.reset_index(drop=True),
            check_dtype=False, check_exact=False, rtol=0.0, atol=tolerance,
        )
    except AssertionError as error:
        raise ValueError(f"{label} values differ") from error


def validate_static_science(
    master: pd.DataFrame, rankings: pd.DataFrame, config,
) -> None:
    """Recompute the configured score, ranks and ties from the static master."""
    if master.columns.tolist() != list(io.MASTER_COLUMNS):
        raise ValueError("static master schema differs")
    if rankings.columns.tolist() != list(io.RANKING_COLUMNS):
        raise ValueError("static ranking schema differs")
    master_sorted = master.sort_values(
        ["pdb_id", "cluster_v2_id"], kind="mergesort",
    ).reset_index(drop=True)
    ranked_shared = rankings.loc[:, list(io.MASTER_COLUMNS)].sort_values(
        ["pdb_id", "cluster_v2_id"], kind="mergesort",
    ).reset_index(drop=True)
    _assert_frame_match(master_sorted, ranked_shared, "static master and ranking shared fields")
    scoring_data = config.data["scoring"]
    ranking_data = config.data["ranking"]
    weights = dict(scoring_data["module_weights"])
    reconstructed = master_sorted.copy()
    reconstructed.attrs.update(
        module_weights=weights, ranking_direction=str(ranking_data["direction"]),
    )
    try:
        expected = scoring.rank_within_target(
            reconstructed, rank_method=str(ranking_data["tie_method"]),
        )
    except (AssertionError, KeyError, TypeError, ValueError) as error:
        raise ValueError("static master violates the configured score formula") from error
    observed = rankings.sort_values(
        ["pdb_id", "Within_PDB_rank", "cluster_v2_id"], kind="mergesort",
    ).reset_index(drop=True)
    _assert_frame_match(expected, observed, "static ranking, tie, or recomputation")
def validate_bundle_relationships(
    candidates: pd.DataFrame, rankings: pd.DataFrame,
    frames: Mapping[str, pd.DataFrame], config,
) -> None:
    """Validate configured diameter and exact identities across stage handoffs."""
    cap = float(config.data["clustering"]["cross_tool"]["cluster_center_diameter_max_A"])
    if float(pd.to_numeric(candidates["cluster_diameter_A"], errors="raise").max()) > cap:
        raise ValueError("cluster diameter exceeds configured cap")
    keys = ["pdb_id", "cluster_v2_id"]
    shared = [column for column in io.CANDIDATE_COLUMNS if column != "cluster_chain_entropy"]
    expected_candidates = candidates.loc[:, shared].sort_values(keys).reset_index(drop=True)
    observed_candidates = rankings.loc[:, shared].sort_values(keys).reset_index(drop=True)
    _assert_frame_match(expected_candidates, observed_candidates, "clustering to static candidate fields")
    expected_entropy = pd.to_numeric(candidates.sort_values(keys)["cluster_chain_entropy"], errors="raise").fillna(0.0)
    observed_entropy = pd.to_numeric(rankings.sort_values(keys)["cluster_chain_entropy"], errors="raise")
    _assert_frame_match(expected_entropy.to_frame(), observed_entropy.to_frame(), "clustering to static entropy")
    final = frames["final_cluster_v2_master_table.csv"]
    expected_static = rankings.loc[:, list(io.RANKING_COLUMNS)].sort_values(keys).reset_index(drop=True)
    observed_static = final.loc[:, list(io.RANKING_COLUMNS)].sort_values(keys).reset_index(drop=True)
    _assert_frame_match(expected_static, observed_static, "static ranking to final master")
    for label, filename, columns in (
        ("evidence labels", "final_automated_evidence_labels.csv", posthoc.CONVENIENCE_COLUMNS[:4]),
        ("reference mapping", "final_reference_mapping.csv", posthoc.CONVENIENCE_COLUMNS[4:]),
    ):
        expected = frames[filename].loc[:, [*keys, *columns]].sort_values(keys).reset_index(drop=True)
        observed = final.loc[:, [*keys, *columns]].sort_values(keys).reset_index(drop=True)
        _assert_frame_match(expected, observed, f"{label} to convenience master")
    supplementary = {
        "weight_sensitivity_scenarios.csv", "weight_sensitivity_targets.csv",
        "single_tool_target_ranks.csv", "single_tool_complete_case_metrics.csv",
    }
    if supplementary.issubset(frames):
        scenario = frames["weight_sensitivity_scenarios.csv"]
        expected_pairs = {
            (module, direction)
            for module in statistics.STATIC_MODULES
            for direction in ("decrease", "increase")
        }
        observed_pairs = set(zip(scenario["perturbed_module"], scenario["direction"]))
        if observed_pairs != expected_pairs or len(scenario) != len(expected_pairs):
            raise ValueError("weight sensitivity must contain the five-by-two scenario set")
        target_ids = set(rankings["pdb_id"].astype(str))
        weight_targets = frames["weight_sensitivity_targets.csv"]
        if set(weight_targets["pdb_id"].astype(str)) != target_ids or len(weight_targets) != len(target_ids):
            raise ValueError("weight-sensitivity targets must equal the static target set")
        scenario_n = len(expected_pairs)
        retained = pd.to_numeric(weight_targets["top1_retention_count"], errors="raise")
        if ((retained < 0) | (retained > scenario_n)).any():
            raise ValueError("weight-sensitivity Top-1 retention count is out of range")

        single = frames["single_tool_target_ranks.csv"]
        methods = {"CavityPlus", "DoGSiteScorer", "DoGSite3", "CASTpFold", "SiteMap", "OIPS-P"}
        if set(single["method"].astype(str)) != methods:
            raise ValueError("single-tool target table has an unexpected method set")
        method_counts = single.groupby("pdb_id")["method"].nunique()
        if set(single["pdb_id"].astype(str)) != target_ids or not method_counts.eq(len(methods)).all():
            raise ValueError("single-tool target table must contain six methods per target")
        complete_text = single["complete_case"].map(lambda value: str(value).strip().casefold())
        complete_targets = set(single.loc[complete_text.eq("true"), "pdb_id"].astype(str))
        metrics = frames["single_tool_complete_case_metrics.csv"]
        if set(metrics["method"].astype(str)) != methods or len(metrics) != len(methods):
            raise ValueError("single-tool metric table must contain each method exactly once")
        metric_n = pd.to_numeric(metrics["N"], errors="raise")
        if not metric_n.eq(len(complete_targets)).all():
            raise ValueError("single-tool metric denominators differ from complete-case targets")
def validate_figure_source_files(
    directory: str | Path, expected_sources: Mapping[str, pd.DataFrame],
) -> None:
    """Require source CSV bytes to equal deterministic tables rebuilt in memory."""
    root = Path(directory)
    for stem, frame in sorted(expected_sources.items()):
        if not isinstance(stem, str) or not stem or frame.empty:
            raise ValueError("figure source expectations must be named nonempty tables")
        output = frame.sort_values(
            ["panel", "record_id", "series"], kind="mergesort", na_position="last",
        ).reset_index(drop=True).copy()
        output = output.map(
            lambda value: "true" if value is True else "false" if value is False else value
        )
        expected = output.to_csv(
            index=False, lineterminator=io.LINE_TERMINATOR,
            float_format=io.FLOAT_FORMAT, na_rep=io.NA_REP,
        )
        path = root / f"{stem}_source_data.csv"
        try:
            observed = path.read_text(encoding=io.ENCODING)
        except OSError as error:
            raise ValueError("figure source data is missing") from error
        if observed != expected:
            raise ValueError(f"figure source data differs from rebuilt values: {stem}")


def _snapshot_numeric_check(
    check_id: str, actual: Mapping[str, object], expected: object,
) -> CheckResult:
    try:
        if not isinstance(expected, Mapping) or set(expected) != set(actual):
            raise ValueError("expected and actual metric names differ")
        failures: list[str] = []
        for name, value in actual.items():
            spec = expected[name]
            if not isinstance(spec, Mapping) or set(spec) != {"value", "abs_tolerance"}:
                raise ValueError(f"{name} must define value and abs_tolerance")
            target, tolerance = spec["value"], spec["abs_tolerance"]
            if any(isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(float(item)) for item in (value, target, tolerance)) or float(tolerance) < 0:
                raise ValueError(f"{name} snapshot values must be finite with nonnegative tolerance")
            if abs(float(value) - float(target)) > float(tolerance):
                failures.append(name)
        if failures:
            raise ValueError("out of tolerance: " + ", ".join(failures))
        return CheckResult(check_id, "pass", "snapshot values match", {"metrics": len(actual)})
    except Exception as error:
        return CheckResult(check_id, "fail", type(error).__name__, {"reason": str(error)})


def _snapshot_exact_check(check_id: str, actual: object, expected: object) -> CheckResult:
    matches = actual == expected
    return CheckResult(
        check_id, "pass" if matches else "fail",
        "snapshot values match" if matches else "snapshot values differ",
        {} if matches else {"actual": actual},
    )


def validate_snapshot(bundle: str | Path, snapshot: str | Path) -> ValidationReport:
    try:
        expected = json.loads(Path(snapshot).read_text(encoding="utf-8"))
        if not isinstance(expected, Mapping):
            raise ValueError("snapshot root must be an object")
    except Exception as error:
        return ValidationReport((CheckResult("snapshot_file", "fail", type(error).__name__, {"reason": _safe_reason(error)}),))
    from .provenance import cluster_summary, snapshot_actual
    summary = cluster_summary(bundle)
    actual = snapshot_actual(bundle, expected)
    count_failures: list[str] = []
    for key, value in summary.items():
        if key == "maximum_diameter_A":
            continue
        if key not in expected:
            count_failures.append(key); continue
        if isinstance(expected[key], bool) or not isinstance(expected[key], (int, float)) or float(value) != float(expected[key]):
            count_failures.append(key)
    diameter = expected.get("cluster_numeric", {}).get("maximum_diameter_A", {})
    flat_diameter = expected.get("maximum_diameter_A")
    if (
        not isinstance(diameter, Mapping)
        or set(diameter) != {"value", "abs_tolerance"}
        or isinstance(flat_diameter, bool)
        or not isinstance(flat_diameter, (int, float))
        or flat_diameter != diameter.get("value")
        or not isinstance(diameter.get("abs_tolerance"), (int, float))
        or abs(float(summary["maximum_diameter_A"]) - float(diameter.get("value", math.nan)))
        > float(diameter.get("abs_tolerance", -1))
    ):
        count_failures.append("maximum_diameter_A")
    counts = CheckResult(
        "snapshot_cluster_summary", "fail" if count_failures else "pass",
        "cluster snapshot mismatch" if count_failures else "cluster snapshot matches",
        {"mismatches": tuple(count_failures)},
    )
    cases_expected = expected.get("representative_case_order")
    cases_ok = isinstance(cases_expected, list) and cases_expected == actual["representative_case_order"]
    supplementary_checks: tuple[CheckResult, ...] = ()
    if "supplementary_analysis_metrics" in expected or "single_tool_complete_case_targets" in expected:
        supplementary_checks = (
            _snapshot_numeric_check(
                "snapshot_supplementary_analyses", actual["supplementary_analysis_metrics"],
                expected.get("supplementary_analysis_metrics"),
            ),
            _snapshot_exact_check(
                "snapshot_single_tool_complete_cases", actual["single_tool_complete_case_targets"],
                expected.get("single_tool_complete_case_targets"),
            ),
        )
    return ValidationReport((
        counts,
        _snapshot_numeric_check("snapshot_cluster_numeric", actual["cluster_numeric"], expected.get("cluster_numeric")),
        _snapshot_exact_check("snapshot_cluster_distributions", actual["cluster_distributions"], expected.get("cluster_distributions")),
        _snapshot_numeric_check("snapshot_exemplar_static", actual.get("exemplar_static_metrics", {}), expected.get("exemplar_static_metrics")),
        _snapshot_exact_check("snapshot_exemplar_state", actual.get("exemplar_static_state"), expected.get("exemplar_static_state")),
        _snapshot_numeric_check("snapshot_analysis_metrics", actual["analysis_metrics"], expected.get("analysis_metrics")),
        _snapshot_numeric_check("snapshot_bootstrap_family", actual["bootstrap_family_intervals"], expected.get("bootstrap_family_intervals")),
        _snapshot_numeric_check("snapshot_orel_ablation", actual["orel_ablation"], expected.get("orel_ablation")),
        *supplementary_checks,
        _snapshot_exact_check("snapshot_posthoc_counts", actual["posthoc_counts"], expected.get("posthoc_counts")),
        _snapshot_exact_check("snapshot_posthoc_distributions", actual["posthoc_distributions"], expected.get("posthoc_distributions")),
        _snapshot_exact_check("snapshot_posthoc_reference_case", actual.get("posthoc_reference_case"), expected.get("posthoc_reference_case")),
        _snapshot_exact_check("snapshot_analysis_counts", actual["analysis_counts"], expected.get("analysis_counts")),
        CheckResult(
            "snapshot_representative_order", "pass" if cases_ok else "fail",
            "representative case order matches" if cases_ok else "representative case order mismatch",
            {"actual": tuple(actual["representative_case_order"])},
        ),
    ))


def numeric_snapshot_sections(snapshot: Mapping[str, object]) -> dict[str, object]:
    return {
        key: snapshot[key]
        for key in (
            "cluster_numeric", "exemplar_static_metrics", "analysis_metrics",
            "bootstrap_family_intervals", "orel_ablation",
        )
    }
