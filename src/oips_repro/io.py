"""Validated public inputs and deterministic cluster-v2 serialization."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import json
import os
from pathlib import Path, PurePosixPath
import re

import numpy as np
import pandas as pd

from .config import ManuscriptConfig, load_manuscript_config


ENCODING = "utf-8"
LINE_TERMINATOR = "\n"
FLOAT_FORMAT = "%.15g"
NA_REP = ""
BOOL_TRUE = "true"
BOOL_FALSE = "false"

FEATURE_COLUMNS = (
    "row_id", "pdb_id", "tool", "pocket_id", "display_order",
    "sitemap_rank", "center_x", "center_y", "center_z", "center_method",
    "residue_count", "residue_set_json", "pocket_geometry_score",
    "pocket_ligandability_score",
)
CANDIDATE_COLUMNS = (
    "pdb_id", "cluster_v2_id", "medoid_unit_id", "medoid_tool",
    "medoid_pocket_id", "medoid_center_x", "medoid_center_y",
    "medoid_center_z", "cluster_diameter_A", "center_dispersion_A",
    "tool_support_count", "supporting_tools", "representative_pockets_per_tool",
    "same_tool_units", "same_tool_secondary_unit_count", "raw_record_count",
    "core_residue_count", "envelope_residue_count", "core_envelope_ratio",
    "core_residues", "envelope_residues", "contributing_chains",
    "contributing_chain_count", "dominant_chain_fraction",
    "cluster_chain_entropy", "mappability", "center_available_representatives",
    "residue_available_representatives", "pairwise_residue_iou_median",
    "pairwise_residue_iou_min", "spatial_continuity", "boundary_sensitive",
)
MEMBERSHIP_COLUMNS = (
    "pdb_id", "cluster_v2_id", "tool", "same_tool_unit_id", "raw_row_id",
    "raw_pocket_id", "formal_tool_representative", "representative_pocket_id",
    "formal_vote_count",
)
MAPPING_COLUMNS = (
    "row_id", "pdb_id", "tool", "pocket_id", "center_x", "center_y",
    "center_z", "residue_count", "same_tool_unit_id", "same_tool_group_size",
    "representative_for_tool_unit", "representative_pocket_id", "cluster_v2_id",
    "mapping_status", "exclusion_reason",
)
EXCLUDED_COLUMNS = (
    "row_id", "pdb_id", "tool", "pocket_id", "center_method",
    "residue_count", "exclusion_reason", "retained_in_audit",
)
BOUNDARY_COLUMNS = (
    "pdb_id", "cluster_v2_id", "tool_support_count",
    "same_tool_secondary_unit_count", "cluster_diameter_A",
    "center_dispersion_A", "pairwise_residue_iou_median",
    "pairwise_residue_iou_min", "core_envelope_ratio", "spatial_continuity",
    "boundary_sensitive", "boundary_reason",
)
MASTER_COLUMNS = CANDIDATE_COLUMNS[:-1] + (
    "C_cons", "G_geo", "P_lig", "Q_evidence", "O_rel_formal",
    "interface_fraction", "interface_recall", "cluster_interface_residue_count",
    "cluster_chain_count", "distance_to_interface_A", "interface_distance_score",
    "boundary_sensitive", "OIPS-P_static",
)
RANKING_COLUMNS = MASTER_COLUMNS + (
    "OIPS-P_static_recomputed", "Within_PDB_rank", "tie_flag", "tie_size",
)

CLUSTER_FILE_SCHEMAS = {
    "cluster_v2_candidates.csv": (CANDIDATE_COLUMNS, ("pdb_id", "cluster_v2_id")),
    "cluster_v2_membership.csv": (MEMBERSHIP_COLUMNS, ("pdb_id", "cluster_v2_id", "tool", "raw_row_id")),
    "tool_record_to_cluster_v2_mapping.csv": (MAPPING_COLUMNS, ("pdb_id", "tool", "row_id")),
    "excluded_unmappable_records.csv": (EXCLUDED_COLUMNS, ("pdb_id", "tool", "row_id")),
    "cluster_v2_boundary_audit.csv": (BOUNDARY_COLUMNS, ("pdb_id", "cluster_v2_id")),
}
SCORE_FILE_SCHEMAS = {
    "cluster_v2_master_table.csv": (MASTER_COLUMNS, ("pdb_id", "cluster_v2_id")),
    "cluster_v2_static_rankings.csv": (RANKING_COLUMNS, ("pdb_id", "Within_PDB_rank", "cluster_v2_id")),
}
_NULLABLE_NUMERIC_FEATURES = (
    "display_order", "sitemap_rank", "center_x", "center_y", "center_z",
    "pocket_geometry_score", "pocket_ligandability_score",
)
_FORBIDDEN_PREFIXES = ("reference", "dcc", "redocking", "literature")


def _normalized_header(header: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", header.casefold()).strip("_")


def _forbidden_header(header: str) -> bool:
    normalized = _normalized_header(header)
    tokens = normalized.split("_") if normalized else []
    if "md" in tokens:
        return True
    if any(
        left == "ligand" and right == "contact"
        for left, right in zip(tokens, tokens[1:])
    ):
        return True
    return any(prefix in tokens for prefix in _FORBIDDEN_PREFIXES)


def _read_header(path: Path) -> list[str]:
    try:
        return pd.read_csv(path, nrows=0, keep_default_na=False).columns.tolist()
    except (OSError, pd.errors.ParserError) as exc:
        raise ValueError(f"unable to read CSV header: {path}") from exc


def _numeric_from_text(series: pd.Series, column: str, *, nullable: bool) -> pd.Series:
    text = series.astype(str)
    missing = text.eq("")
    if not nullable and missing.any():
        raise ValueError(f"{column} may not be missing")
    parsed = pd.to_numeric(text.mask(missing), errors="coerce")
    invalid = ~missing & parsed.isna()
    if invalid.any():
        raise ValueError(f"numeric field {column} contains invalid text")
    nonfinite = ~missing & ~np.isfinite(parsed.astype(float))
    if nonfinite.any():
        raise ValueError(f"numeric field {column} must be finite when present")
    return parsed


def load_feature_table(path: str | Path) -> pd.DataFrame:
    """Load a positive, sorted feature table without post-hoc leakage."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"feature table does not exist: {source}")
    headers = _read_header(source)
    forbidden = [header for header in headers if _forbidden_header(header)]
    if forbidden:
        raise ValueError(f"forbidden post-hoc evidence header(s): {forbidden}")
    if headers != list(FEATURE_COLUMNS):
        raise ValueError(
            f"feature table must have the exact 14-column allowlist in order; got {headers}"
        )
    frame = pd.read_csv(source, dtype=str, keep_default_na=False, encoding=ENCODING)
    if frame.empty:
        raise ValueError("feature table must contain at least one record")
    row_ids = _numeric_from_text(frame["row_id"], "row_id", nullable=False)
    if (row_ids % 1 != 0).any() or row_ids.duplicated().any() or row_ids.nunique() != len(frame):
        raise ValueError("row_id must contain unique integers")
    frame["row_id"] = row_ids.astype("int64")
    for column in ("pdb_id", "tool", "pocket_id"):
        if frame[column].eq("").any():
            raise ValueError(f"{column} must contain nonempty strings")
        frame[column] = frame[column].astype(str)
    if not frame["pdb_id"].str.fullmatch(r"[0-9][A-Z0-9]{3}").all():
        raise ValueError("pdb_id values must be uppercase four-character PDB IDs")
    for column in _NULLABLE_NUMERIC_FEATURES:
        frame[column] = _numeric_from_text(frame[column], column, nullable=True)
    for column in ("display_order", "sitemap_rank"):
        present = frame[column].dropna()
        if (present % 1 != 0).any():
            raise ValueError(f"{column} must contain integers when present")
    residue_counts = _numeric_from_text(frame["residue_count"], "residue_count", nullable=False)
    if (residue_counts % 1 != 0).any() or (residue_counts < 0).any():
        raise ValueError("residue_count must be a nonnegative integer")
    frame["residue_count"] = residue_counts.astype("int64")
    parsed_residues: list[list[object]] = []
    for row_number, value in zip(frame["row_id"], frame["residue_set_json"]):
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError as exc:
            raise ValueError(f"residue_set_json must contain a JSON array (row_id={row_number})") from exc
        if not isinstance(parsed, list):
            raise ValueError(f"residue_set_json must contain a JSON array (row_id={row_number})")
        if any(not isinstance(item, str) or not item for item in parsed):
            raise ValueError(f"residue_set_json entries must be nonempty strings (row_id={row_number})")
        parsed_residues.append(parsed)
    center_counts = frame[["center_x", "center_y", "center_z"]].notna().sum(axis=1)
    if (~center_counts.isin([0, 3])).any():
        raise ValueError("center coordinates must be all present or all missing")
    lengths = pd.Series([len(values) for values in parsed_residues], index=frame.index)
    if not lengths.eq(frame["residue_count"]).all():
        raise ValueError("residue_count must equal len(residue_set_json)")
    sort_keys = ["pdb_id", "tool", "row_id"]
    sorted_frame = frame.sort_values(sort_keys, kind="mergesort").reset_index(drop=True)
    if not frame.reset_index(drop=True)[sort_keys].equals(sorted_frame[sort_keys]):
        raise ValueError("feature table must be stably sorted by pdb_id, tool, row_id")
    return frame.loc[:, list(FEATURE_COLUMNS)].copy()


def _lowercase_boolean(value: object) -> object:
    if isinstance(value, (bool, np.bool_)):
        return BOOL_TRUE if bool(value) else BOOL_FALSE
    return value


def write_stable_csv(
    frame: pd.DataFrame, path: str | Path, *, columns: Sequence[str],
    sort_by: Sequence[str], delimiter: str = ",",
) -> None:
    missing = [column for column in [*columns, *sort_by] if column not in frame]
    if missing:
        raise KeyError(f"missing columns required for serialization: {sorted(set(missing))}")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output = frame.loc[:, list(columns)].copy()
    output = output.sort_values(list(sort_by), kind="mergesort", na_position="last")
    output = output.map(_lowercase_boolean)
    with destination.open("w", encoding=ENCODING, newline="") as handle:
        output.to_csv(
            handle, index=False, sep=delimiter, encoding=ENCODING,
            lineterminator=LINE_TERMINATOR, float_format=FLOAT_FORMAT,
            na_rep=NA_REP,
        )


def _validate_exact_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    if frame.columns.tolist() != list(columns):
        raise ValueError(f"{label} has incorrect schema: {frame.columns.tolist()}")


def _validate_clustering_relationships(result: object) -> None:
    candidates = result.candidates
    membership = result.membership
    mapping = result.mapping
    excluded = result.excluded
    boundary = result.boundary
    if candidates["cluster_v2_id"].duplicated().any():
        raise ValueError("candidate cluster IDs must be unique")
    if boundary["cluster_v2_id"].duplicated().any():
        raise ValueError("boundary cluster IDs must be unique")
    valid_mappability = {"center_and_residue_mappable", "center_only_mappable", "residue_only_mappable"}
    if not set(candidates["mappability"]).issubset(valid_mappability):
        raise ValueError("candidate mappability contains an invalid value")
    valid_status = {"mapped_to_cluster_v2", "excluded_unmappable"}
    if not set(mapping["mapping_status"]).issubset(valid_status):
        raise ValueError("mapping_status contains an invalid value")
    if mapping["row_id"].duplicated().any():
        raise ValueError("mapping row_id values must be unique")
    if membership["raw_row_id"].duplicated().any():
        raise ValueError("membership raw_row_id values must be unique")
    candidate_ids = set(candidates["cluster_v2_id"])
    if set(membership["cluster_v2_id"]) != candidate_ids:
        raise ValueError("membership cluster IDs must equal candidate cluster IDs")
    if set(boundary["cluster_v2_id"]) != candidate_ids:
        raise ValueError("boundary cluster IDs must equal candidate cluster IDs")
    mapped = mapping.loc[mapping["mapping_status"].eq("mapped_to_cluster_v2")]
    if set(mapped["row_id"]) != set(membership["raw_row_id"]):
        raise ValueError("mapped row IDs must equal membership raw row IDs")
    mapped_clusters = dict(zip(mapped["row_id"], mapped["cluster_v2_id"]))
    member_clusters = dict(zip(membership["raw_row_id"], membership["cluster_v2_id"]))
    if mapped_clusters != member_clusters:
        raise ValueError("mapping and membership must assign each row to the same cluster")
    if not set(mapped["cluster_v2_id"]).issubset(candidate_ids):
        raise ValueError("mapping contains an unknown candidate cluster ID")
    if set(excluded["row_id"]) != set(mapping.loc[mapping["mapping_status"].eq("excluded_unmappable"), "row_id"]):
        raise ValueError("excluded row IDs must equal excluded mapping row IDs")
    votes = membership.groupby(["cluster_v2_id", "tool"])["formal_vote_count"].sum()
    if not votes.empty and (votes > 1).any():
        raise ValueError("formal votes may not exceed one per cluster and tool")
    formal_flags = membership["formal_tool_representative"].eq(True)
    if not membership["formal_vote_count"].eq(formal_flags.astype(int)).all():
        raise ValueError("formal vote counts must match representative flags")


def _numeric_identity(value: object) -> str:
    if value is None or pd.isna(value):
        return NA_REP
    return FLOAT_FORMAT % float(value)


def _validate_cluster_handoff(result: object, features: pd.DataFrame) -> None:
    required = {"row_id", "pdb_id", "tool", "pocket_id", "center_x", "center_y", "center_z", "residue_count"}
    if not required.issubset(features.columns):
        raise ValueError("supplied feature table lacks mapping identity fields")
    if features["row_id"].duplicated().any():
        raise ValueError("supplied feature table row_id values must be unique")

    observed = result.mapping.set_index("row_id").sort_index()
    expected = features.set_index("row_id").sort_index()
    if not observed.index.equals(expected.index):
        raise ValueError("cluster mapping row IDs do not match the supplied feature table")
    for column in ("pdb_id", "tool", "pocket_id"):
        left = observed[column].map(str)
        right = expected[column].map(str)
        mismatch = ~left.eq(right)
        if mismatch.any():
            row_id = int(mismatch.index[mismatch.argmax()])
            raise ValueError(
                f"cluster mapping identity mismatch for {column} at row_id {row_id}"
            )
    residue_mismatch = ~pd.to_numeric(observed["residue_count"]).eq(
        pd.to_numeric(expected["residue_count"])
    )
    if residue_mismatch.any():
        row_id = int(residue_mismatch.index[residue_mismatch.argmax()])
        raise ValueError(
            f"cluster mapping identity mismatch for residue_count at row_id {row_id}"
        )
    for column in ("center_x", "center_y", "center_z"):
        left = observed[column].map(_numeric_identity)
        right = expected[column].map(_numeric_identity)
        mismatch = ~left.eq(right)
        if mismatch.any():
            row_id = int(mismatch.index[mismatch.argmax()])
            raise ValueError(
                f"cluster mapping identity mismatch for {column} at row_id {row_id}"
            )


def _prepare_output_directory(path: Path, expected_files: set[str]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    extras = {entry.name for entry in path.iterdir()} - expected_files
    if extras:
        raise ValueError(f"output directory contains unexpected entries: {sorted(extras)}")


def write_clustering_result(result: object, output_dir: str | Path) -> None:
    frames = {
        "cluster_v2_candidates.csv": result.candidates,
        "cluster_v2_membership.csv": result.membership,
        "tool_record_to_cluster_v2_mapping.csv": result.mapping,
        "excluded_unmappable_records.csv": result.excluded,
        "cluster_v2_boundary_audit.csv": result.boundary,
    }
    for filename, frame in frames.items():
        _validate_exact_columns(frame, CLUSTER_FILE_SCHEMAS[filename][0], filename)
    _validate_clustering_relationships(result)
    destination = Path(output_dir)
    _prepare_output_directory(destination, set(CLUSTER_FILE_SCHEMAS))
    for filename, frame in frames.items():
        columns, sort_by = CLUSTER_FILE_SCHEMAS[filename]
        write_stable_csv(frame, destination / filename, columns=columns, sort_by=sort_by)


_RESULT_FLOAT_COLUMNS = {
    "medoid_center_x", "medoid_center_y", "medoid_center_z", "cluster_diameter_A",
    "center_dispersion_A", "core_envelope_ratio", "dominant_chain_fraction",
    "cluster_chain_entropy", "pairwise_residue_iou_median", "pairwise_residue_iou_min",
    "center_x", "center_y", "center_z", "C_cons", "G_geo", "P_lig",
    "Q_evidence", "O_rel_formal", "interface_fraction", "interface_recall",
    "distance_to_interface_A", "interface_distance_score", "OIPS-P_static",
    "OIPS-P_static_recomputed", "Within_PDB_rank",
}
_RESULT_INTEGER_COLUMNS = {
    "tool_support_count", "same_tool_secondary_unit_count", "raw_record_count",
    "core_residue_count", "envelope_residue_count", "contributing_chain_count",
    "center_available_representatives", "residue_available_representatives",
    "raw_row_id", "formal_vote_count", "row_id", "residue_count",
    "same_tool_group_size", "cluster_interface_residue_count", "cluster_chain_count",
    "tie_size",
}
_RESULT_BOOLEAN_COLUMNS = {
    "spatial_continuity", "boundary_sensitive", "formal_tool_representative",
    "representative_for_tool_unit", "retained_in_audit", "tie_flag",
}


def _parse_boolean_series(series: pd.Series, column: str) -> pd.Series:
    values = series.astype(str)
    invalid = ~values.isin(["", BOOL_TRUE, BOOL_FALSE])
    if invalid.any():
        raise ValueError(f"boolean column {column} must use lowercase true/false")
    parsed = values.map({BOOL_TRUE: True, BOOL_FALSE: False, "": None})
    return parsed.astype(object)


def _read_result_csv(path: Path, columns: Sequence[str], sort_by: Sequence[str]) -> pd.DataFrame:
    headers = _read_header(path)
    if headers != list(columns):
        raise ValueError(f"{path.name} has incorrect schema: {headers}")
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, encoding=ENCODING)
    for column in frame.columns:
        if column in _RESULT_FLOAT_COLUMNS:
            frame[column] = _numeric_from_text(frame[column], column, nullable=True)
        elif column in _RESULT_INTEGER_COLUMNS:
            numeric = _numeric_from_text(frame[column], column, nullable=True)
            if ((numeric.dropna() % 1) != 0).any():
                raise ValueError(f"integer column {column} contains fractional values")
            frame[column] = numeric.astype("Int64") if numeric.isna().any() else numeric.astype("int64")
        elif column in _RESULT_BOOLEAN_COLUMNS:
            frame[column] = _parse_boolean_series(frame[column], column)
    expected = frame.sort_values(list(sort_by), kind="mergesort", na_position="last").reset_index(drop=True)
    if not frame.reset_index(drop=True).equals(expected):
        raise ValueError(f"{path.name} is not sorted by {list(sort_by)}")
    return frame


def load_clustering_result(
    cluster_dir: str | Path, *, features: pd.DataFrame
):
    from .clustering import ClusteringResult

    source = Path(cluster_dir)
    if not source.is_dir():
        raise FileNotFoundError(f"cluster result directory does not exist: {source}")
    actual = {entry.name for entry in source.iterdir() if entry.is_file()}
    if actual != set(CLUSTER_FILE_SCHEMAS):
        raise ValueError(f"cluster directory file set is invalid: {sorted(actual)}")
    loaded = {
        filename: _read_result_csv(source / filename, *schema)
        for filename, schema in CLUSTER_FILE_SCHEMAS.items()
    }
    result = ClusteringResult(
        candidates=loaded["cluster_v2_candidates.csv"],
        membership=loaded["cluster_v2_membership.csv"],
        mapping=loaded["tool_record_to_cluster_v2_mapping.csv"],
        excluded=loaded["excluded_unmappable_records.csv"],
        boundary=loaded["cluster_v2_boundary_audit.csv"],
        features=features.loc[:, list(FEATURE_COLUMNS)].copy(),
    )
    _validate_clustering_relationships(result)
    _validate_cluster_handoff(result, features)
    return result


def write_scoring_result(result: object, output_dir: str | Path) -> None:
    frames = {
        "cluster_v2_master_table.csv": result.master,
        "cluster_v2_static_rankings.csv": result.rankings,
    }
    for filename, frame in frames.items():
        _validate_exact_columns(frame, SCORE_FILE_SCHEMAS[filename][0], filename)
    if set(result.master["cluster_v2_id"]) != set(result.rankings["cluster_v2_id"]):
        raise ValueError("master and ranking cluster IDs must match")
    destination = Path(output_dir)
    _prepare_output_directory(destination, set(SCORE_FILE_SCHEMAS))
    for filename, frame in frames.items():
        columns, sort_by = SCORE_FILE_SCHEMAS[filename]
        write_stable_csv(frame, destination / filename, columns=columns, sort_by=sort_by)


def path_is_equal_or_within(path: str | Path, root: str | Path) -> bool:
    # Treat repository safety boundaries case-insensitively on every platform.
    # A checkout may be prepared on Linux and later used on a case-insensitive
    # Windows or macOS filesystem, so relying on the host-specific behaviour of
    # ``normcase`` alone can approve an output path that aliases the frozen
    # reference directory after transfer.
    candidate_key = os.path.normcase(str(Path(path).resolve())).casefold()
    root_key = os.path.normcase(str(Path(root).resolve())).casefold()
    try:
        return os.path.commonpath([candidate_key, root_key]) == root_key
    except ValueError:
        return False


def ensure_safe_output_path(output: str | Path, reference: str | Path) -> Path:
    destination = Path(output).resolve()
    if path_is_equal_or_within(destination, reference):
        raise ValueError("refusing to write equal to or within configured reference results")
    return destination


def load_structure_paths(
    config: ManuscriptConfig, *, target_ids: Sequence[str] | None = None,
) -> dict[str, Path]:
    manifest_path = config.resolve_configured_path("structure_manifest")
    if not manifest_path.is_file():
        raise FileNotFoundError(f"structure manifest does not exist: {manifest_path}")
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"asset_id", "scientific_role", "local_path_or_url"}
    if not required.issubset(manifest.columns):
        raise ValueError("structure manifest is missing required columns")
    selected = manifest.loc[manifest["scientific_role"].eq("prepared paired structure")]
    if selected.empty:
        raise ValueError("structure manifest must identify prepared paired structures")
    if selected["asset_id"].duplicated().any() or selected["local_path_or_url"].duplicated().any():
        raise ValueError("prepared structure manifest IDs and paths must be unique")
    structures_root = config.resolve_configured_path("structures_dir")
    structures: dict[str, Path] = {}
    for row in selected.sort_values("asset_id", kind="mergesort").itertuples(index=False):
        match = re.fullmatch(r"structure-([a-z0-9]{4})", str(row.asset_id))
        if not match:
            raise ValueError(f"invalid prepared structure asset ID: {row.asset_id}")
        relative = PurePosixPath(str(row.local_path_or_url))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"prepared structure path must be repository-relative: {relative}")
        local_path = (config.repository_root / Path(*relative.parts)).resolve()
        try:
            local_path.relative_to(structures_root)
        except ValueError as error:
            raise ValueError("prepared structure path must be within structures_dir") from error
        if not local_path.is_file() or local_path.is_symlink():
            raise FileNotFoundError(f"prepared structure file does not exist: {local_path}")
        pdb_id = match.group(1).upper()
        if pdb_id in structures:
            raise ValueError(f"duplicate prepared structure asset ID for {pdb_id}")
        structures[pdb_id] = local_path
    if target_ids is not None:
        expected = {str(value) for value in target_ids}
        if structures.keys() != expected:
            missing = sorted(expected - structures.keys())
            extra = sorted(structures.keys() - expected)
            raise ValueError(f"prepared structures must exactly cover targets; missing={missing}, extra={extra}")
    return structures
