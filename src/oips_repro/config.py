"""Strict, side-effect-free loading of the public OIPS configuration."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from pathlib import Path, PurePosixPath
import re
from typing import Any

import yaml


ROOT_KEYS = (
    "schema_version", "analysis", "paths", "tool_weights", "clustering",
    "scoring", "ranking", "posthoc", "statistics",
)
PATH_KEYS = (
    "feature_table", "systems", "structure_manifest", "structures_dir",
    "reference_annotations", "md_evidence", "redocking_evidence",
    "legacy_crosswalk", "figure_contract", "reference_results",
    "default_output", "schema", "data_manifest", "data_checksums",
    "asset_rights", "external_archive_manifest", "expected_summary",
)
MODULE_KEYS = ("C_cons", "G_geo", "P_lig", "O_rel_formal", "Q_evidence")

SAME_TOOL_KEYS = (
    "residue_iou_strong", "residue_containment_nested", "nested_center_max_A",
    "near_center_max_A", "near_center_min_iou", "center_only_duplicate_max_A",
    "dogsite_hierarchy_center_max_A", "dogsite_hierarchy_min_iou",
    "dogsite_hierarchy_containment_min", "group_center_diameter_max_A",
    "representative_missing_center_penalty",
    "representative_missing_residue_penalty",
)
CROSS_TOOL_KEYS = (
    "close_center_max_A", "conditional_center_max_A",
    "conditional_residue_iou_min", "residue_iou_min",
    "cluster_center_diameter_max_A",
)
Q_EVIDENCE_KEYS = (
    "base", "per_representative", "center_fraction", "residue_fraction",
    "sitemap_bonus",
)
INTERFACE_WEIGHT_KEYS = (
    "interface_fraction", "chain_context", "distance", "interface_extent",
)
DISTANCE_SCORE_KEYS = (
    "near_max_A", "near_score", "intermediate_max_A", "intermediate_start",
    "intermediate_loss_per_A", "far_max_A", "far_start",
    "far_loss_per_A", "tail_start", "tail_loss_per_A", "tail_floor",
)
INTERFACE_KEYS = (
    "contact_cutoff_A", "weights", "distance_score", "interface_fraction_base",
    "interface_fraction_multiplier", "multi_chain_base",
    "chain_entropy_multiplier", "single_chain_interface_fraction_min",
    "single_chain_supported_score", "single_chain_other_score", "context_base",
    "context_per_extra_chain", "context_extra_chain_cap",
    "context_interface_fraction_cap", "context_interface_fraction_multiplier",
    "no_interface_monomer_score", "no_interface_multichain_score",
)
BOUNDARY_KEYS = (
    "diameter_gt_A", "dispersion_gt_A", "core_envelope_ratio_lt",
    "median_residue_iou_lt",
)
SCORING_KEYS = (
    "module_weights", "geometry_top_n", "ligandability_top_n", "q_evidence",
    "interface", "boundary",
)
TOP3_QC_KEYS = (
    "top_k", "overmerge_diameter_gt_A", "overmerge_core_ratio_lt",
    "overmerge_secondary_units_min", "overmerge_tool_support_max",
    "overmerge_median_iou_lt", "exclusion_diameter_gt_A",
    "exclusion_core_ratio_lt", "exclusion_tool_support_min",
    "split_tool_support_min", "split_dispersion_gt_A", "split_median_iou_lt",
    "split_core_ratio_lt", "surface_noise_tool_support",
    "surface_noise_geometry_lt", "surface_noise_ligandability_lt",
    "surface_noise_interface_fraction_lt",
)
REFERENCE_KEYS = (
    "contact_cutoff_A", "near_dcc_max_A", "near_recall_min",
    "middle_dcc_max_A", "middle_recall_min", "middle_iou_min",
    "far_dcc_max_A", "far_recall_min", "far_iou_min",
)
MD_KEYS = (
    "concordant_iou_min", "concordant_precision_min",
    "concordant_center_max_A", "partial_iou_min", "alternative_iou_min",
    "alternative_center_max_A", "boundary_shift_center_max_A",
)


def _mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    if any(not isinstance(key, str) for key in value):
        raise ValueError(f"{label} keys must be strings")
    return dict(value)


def _exact(mapping: Mapping[str, object], keys: tuple[str, ...], label: str) -> None:
    actual, expected = set(mapping), set(keys)
    unknown, missing = sorted(actual - expected), sorted(expected - actual)
    if unknown:
        raise ValueError(f"unknown {label} keys: {unknown}")
    if missing:
        raise ValueError(f"missing {label} keys: {missing}")


def _number(value: object, label: str, *, minimum: float | None = None,
            maximum: float | None = None, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    if positive and result <= 0:
        raise ValueError(f"{label} must be positive")
    if minimum is not None and result < minimum:
        raise ValueError(f"{label} must be >= {minimum}")
    if maximum is not None and result > maximum:
        raise ValueError(f"{label} must be <= {maximum}")
    return result


def _integer(value: object, label: str, *, positive: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer (not a boolean)")
    if positive and value <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return value


def _fraction(value: object, label: str) -> float:
    return _number(value, label, minimum=0.0, maximum=1.0)


def _numbers(section: Mapping[str, object], label: str) -> None:
    for key, value in section.items():
        _number(value, f"{label}.{key}", minimum=0.0)


def _strict_path(value: object, label: str, root: Path) -> None:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a nonempty repository-relative POSIX path")
    if "\\" in value:
        raise ValueError(f"{label} must use POSIX separators; backslash is forbidden")
    path = PurePosixPath(value)
    if path.is_absolute() or re.match(r"^[A-Za-z]:", value):
        raise ValueError(f"{label} must be repository-relative")
    if ".." in path.parts or "." in path.parts:
        raise ValueError(f"{label} must not contain . or .. segments")
    resolved = (root / Path(*path.parts)).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} escapes the repository root") from error


def _validate_config(data: dict[str, Any], root: Path) -> None:
    _exact(data, ROOT_KEYS, "root")
    if data["schema_version"] != 1 or isinstance(data["schema_version"], bool):
        raise ValueError("schema_version must be integer 1")

    analysis = _mapping(data["analysis"], "analysis")
    _exact(analysis, ("id", "representative_case_ids"), "analysis")
    if not isinstance(analysis["id"], str) or not analysis["id"].strip():
        raise ValueError("analysis.id must be a nonempty string")
    cases = analysis["representative_case_ids"]
    if not isinstance(cases, list) or not cases:
        raise ValueError("analysis.representative_case_ids must be a nonempty list")
    if any(not isinstance(value, str) or not re.fullmatch(r"[0-9][A-Z0-9]{3}", value) for value in cases):
        raise ValueError("representative case IDs must be unique uppercase PDB IDs")
    if len(set(cases)) != len(cases):
        raise ValueError("representative case IDs must be unique uppercase PDB IDs")

    paths = _mapping(data["paths"], "paths")
    _exact(paths, PATH_KEYS, "paths")
    for key, value in paths.items():
        _strict_path(value, f"paths.{key}", root)

    tools = _mapping(data["tool_weights"], "tool_weights")
    if not tools or any(not isinstance(key, str) or not key for key in tools):
        raise ValueError("tool_weights must contain named tools")
    for key, value in tools.items():
        _number(value, f"tool_weights.{key}", positive=True)

    clustering = _mapping(data["clustering"], "clustering")
    _exact(clustering, ("same_tool", "cross_tool"), "clustering")
    same = _mapping(clustering["same_tool"], "clustering.same_tool")
    cross = _mapping(clustering["cross_tool"], "clustering.cross_tool")
    _exact(same, SAME_TOOL_KEYS, "clustering.same_tool")
    _exact(cross, CROSS_TOOL_KEYS, "clustering.cross_tool")
    for key in SAME_TOOL_KEYS:
        (_fraction if any(word in key for word in ("iou", "containment", "penalty")) else _number)(
            same[key], f"clustering.same_tool.{key}", **({} if any(word in key for word in ("iou", "containment", "penalty")) else {"minimum": 0.0})
        )
    for key in CROSS_TOOL_KEYS:
        (_fraction if "iou" in key else _number)(
            cross[key], f"clustering.cross_tool.{key}", **({} if "iou" in key else {"minimum": 0.0})
        )
    if not (cross["close_center_max_A"] <= cross["conditional_center_max_A"] <= cross["cluster_center_diameter_max_A"]):
        raise ValueError("clustering close center must be <= conditional center <= diameter cap")

    scoring = _mapping(data["scoring"], "scoring")
    _exact(scoring, SCORING_KEYS, "scoring")
    modules = _mapping(scoring["module_weights"], "scoring.module_weights")
    _exact(modules, MODULE_KEYS, "scoring.module_weights")
    for key, value in modules.items():
        _number(value, f"scoring.module_weights.{key}", positive=True)
    if not math.isclose(sum(map(float, modules.values())), 1.0, abs_tol=1e-12):
        raise ValueError("scoring.module_weights must sum to one")
    _integer(scoring["geometry_top_n"], "scoring.geometry_top_n", positive=True)
    _integer(scoring["ligandability_top_n"], "scoring.ligandability_top_n", positive=True)
    quality = _mapping(scoring["q_evidence"], "scoring.q_evidence")
    _exact(quality, Q_EVIDENCE_KEYS, "scoring.q_evidence"); _numbers(quality, "scoring.q_evidence")
    boundary = _mapping(scoring["boundary"], "scoring.boundary")
    _exact(boundary, BOUNDARY_KEYS, "scoring.boundary")
    for key, value in boundary.items():
        (_fraction if key.endswith("_lt") else _number)(
            value, f"scoring.boundary.{key}", **({} if key.endswith("_lt") else {"minimum": 0.0})
        )
    interface = _mapping(scoring["interface"], "scoring.interface")
    _exact(interface, INTERFACE_KEYS, "scoring.interface")
    weights = _mapping(interface["weights"], "scoring.interface.weights")
    _exact(weights, INTERFACE_WEIGHT_KEYS, "scoring.interface.weights")
    for key, value in weights.items(): _number(value, f"scoring.interface.weights.{key}", positive=True)
    if not math.isclose(sum(map(float, weights.values())), 1.0, abs_tol=1e-12):
        raise ValueError("scoring.interface.weights must sum to one")
    distance = _mapping(interface["distance_score"], "scoring.interface.distance_score")
    _exact(distance, DISTANCE_SCORE_KEYS, "scoring.interface.distance_score"); _numbers(distance, "scoring.interface.distance_score")
    if not (distance["near_max_A"] <= distance["intermediate_max_A"] <= distance["far_max_A"]):
        raise ValueError("interface distance thresholds must be ordered")
    for key in set(INTERFACE_KEYS) - {"weights", "distance_score", "context_extra_chain_cap"}:
        if key in {"single_chain_interface_fraction_min", "context_interface_fraction_cap"}:
            _fraction(interface[key], f"scoring.interface.{key}")
        else:
            _number(interface[key], f"scoring.interface.{key}", minimum=0.0)
    _integer(interface["context_extra_chain_cap"], "scoring.interface.context_extra_chain_cap", positive=True)

    ranking = _mapping(data["ranking"], "ranking")
    _exact(ranking, ("direction", "tie_method"), "ranking")
    if ranking != {"direction": "descending", "tie_method": "average"}:
        raise ValueError("ranking must be exactly descending/average")

    posthoc = _mapping(data["posthoc"], "posthoc")
    _exact(posthoc, ("top3_qc", "reference", "md"), "posthoc")
    qc = _mapping(posthoc["top3_qc"], "posthoc.top3_qc")
    ref = _mapping(posthoc["reference"], "posthoc.reference")
    md = _mapping(posthoc["md"], "posthoc.md")
    _exact(qc, TOP3_QC_KEYS, "posthoc.top3_qc")
    _exact(ref, REFERENCE_KEYS, "posthoc.reference")
    _exact(md, MD_KEYS, "posthoc.md")
    integer_qc = {"top_k", "overmerge_secondary_units_min", "overmerge_tool_support_max", "exclusion_tool_support_min", "split_tool_support_min", "surface_noise_tool_support"}
    fraction_qc = {key for key in TOP3_QC_KEYS if any(word in key for word in ("ratio", "iou", "fraction"))}
    for key, value in qc.items():
        if key in integer_qc: _integer(value, f"posthoc.top3_qc.{key}", positive=True)
        elif key in fraction_qc: _fraction(value, f"posthoc.top3_qc.{key}")
        else: _number(value, f"posthoc.top3_qc.{key}", minimum=0.0)
    for key, value in ref.items():
        (_fraction if key.endswith(("recall_min", "iou_min")) else _number)(
            value, f"posthoc.reference.{key}", **({} if key.endswith(("recall_min", "iou_min")) else {"minimum": 0.0})
        )
    if not (ref["near_dcc_max_A"] <= ref["middle_dcc_max_A"] <= ref["far_dcc_max_A"]):
        raise ValueError("reference near, middle and far distance thresholds must be ordered")
    for key, value in md.items():
        (_fraction if "iou" in key or "precision" in key else _number)(
            value, f"posthoc.md.{key}", **({} if "iou" in key or "precision" in key else {"minimum": 0.0})
        )
    if not (md["concordant_center_max_A"] <= md["alternative_center_max_A"] <= md["boundary_shift_center_max_A"]):
        raise ValueError("MD center thresholds must be ordered")

    statistics = _mapping(data["statistics"], "statistics")
    _exact(statistics, ("bootstrap_iterations", "random_seed"), "statistics")
    _integer(statistics["bootstrap_iterations"], "statistics.bootstrap_iterations", positive=True)
    _integer(statistics["random_seed"], "statistics.random_seed")

    contract_path = (root / Path(*PurePosixPath(paths["figure_contract"]).parts)).resolve()
    if not contract_path.is_file() or contract_path.is_symlink():
        raise FileNotFoundError("figure contract must be a regular file")
    contract = yaml.safe_load(contract_path.read_text(encoding="utf-8"))
    titles = contract.get("case_titles") if isinstance(contract, Mapping) else None
    if not isinstance(titles, Mapping) or set(titles) != set(cases):
        raise ValueError("figure contract case titles must cover representative cases exactly")


@dataclass(frozen=True)
class ManuscriptConfig:
    path: Path
    repository_root: Path
    data: Mapping[str, Any]

    def resolve_configured_path(self, key: str) -> Path:
        paths = _mapping(self.data.get("paths"), "paths")
        if key not in paths:
            raise KeyError(f"configuration path is missing: {key}")
        posix = PurePosixPath(str(paths[key]))
        resolved = (self.repository_root / Path(*posix.parts)).resolve()
        try:
            resolved.relative_to(self.repository_root)
        except ValueError as error:
            raise ValueError(f"configured path escapes repository root: {key}") from error
        return resolved


def load_manuscript_config(path: str | Path) -> ManuscriptConfig:
    config_path = Path(path).resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    data = _mapping(raw, "configuration root")
    repository_root = config_path.parent.parent.resolve()
    _validate_config(data, repository_root)
    return ManuscriptConfig(config_path, repository_root, data)


def dump_manuscript_config(data: Mapping[str, object]) -> str:
    """Serialize configuration safely for fixture construction and review."""
    return yaml.safe_dump(dict(data), sort_keys=False, allow_unicode=True)
