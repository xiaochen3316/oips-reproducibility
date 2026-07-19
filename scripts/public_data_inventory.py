"""Build a path-scrubbed inventory of large assets excluded from Git."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from public_data_tables import UNAVAILABLE, normalize_missing_value, require_columns


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
POCKET_PATH_FIELDS = [
    "result_table_file",
    "raw_download_file",
    "measure_file",
    "status_file",
    "sitemap_result_dir",
    "sitemap_smap_file",
    "sitemap_log_file",
    "sitemap_out_maegz",
]
MD_ARCHIVE_CLASSES = {
    "full MD inputs": {
        ".cms",
        ".cfg",
        ".msj",
        ".mae",
        ".maegz",
        ".top",
        ".itp",
        ".prm",
        ".tpr",
    },
    "full MD trajectories": {".trj", ".xtc", ".trr", ".dcd", ".nc", ".traj"},
    "restart or checkpoint files": {
        ".cpt",
        ".chk",
        ".rst",
        ".restart",
        ".checkpoint",
    },
    "representative MD snapshots": {".pdb", ".cif", ".mmcif", ".gro"},
}


def slug(value: object) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "-", str(value).lower()).strip("-")
    return normalized or "unlabeled"


def _under_source_root(path: Path, source_root: Path) -> bool:
    try:
        path.relative_to(source_root)
    except ValueError:
        return False
    return True


def resolve_internal_reference(raw: object, source_root: Path) -> Path | None:
    cleaned = normalize_missing_value(raw)
    if cleaned is pd.NA or pd.isna(cleaned):
        return None
    candidate = Path(str(cleaned))
    if not candidate.is_absolute():
        candidate = source_root / candidate
    candidate = candidate.resolve(strict=False)
    source_root = source_root.resolve(strict=True)
    if not _under_source_root(candidate, source_root):
        return None
    return candidate


def expand_files(paths: Iterable[Path]) -> list[Path]:
    files: set[Path] = set()
    for path in paths:
        if path.is_file():
            files.add(path.resolve())
        elif path.is_dir():
            files.update(item.resolve() for item in path.rglob("*") if item.is_file())
    return sorted(files, key=lambda item: str(item).casefold())


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_collection(files: Iterable[Path], source_root: Path) -> tuple[int, str]:
    unique = sorted(set(files), key=lambda item: str(item).casefold())
    total = 0
    collection_digest = hashlib.sha256()
    for path in unique:
        if not _under_source_root(path.resolve(), source_root.resolve()):
            raise ValueError("external inventory attempted to fingerprint outside source root")
        total += path.stat().st_size
        relative = path.resolve().relative_to(source_root.resolve()).as_posix()
        content_digest = file_sha256(path)
        collection_digest.update(relative.encode("utf-8"))
        collection_digest.update(b"\0")
        collection_digest.update(content_digest.encode("ascii"))
        collection_digest.update(b"\n")
    return total, collection_digest.hexdigest()


def availability_fields(
    files: list[Path], source_root: Path
) -> tuple[str, int | str, str, str, str]:
    if not files:
        return (
            UNAVAILABLE,
            UNAVAILABLE,
            UNAVAILABLE,
            UNAVAILABLE,
            "requested artifact class was not present in the discoverable source tree",
        )
    byte_count, digest = fingerprint_collection(files, source_root)
    return (
        "available_locally",
        byte_count,
        digest,
        "not_archived",
        "local collection fingerprint verified; persistent deposit not assigned",
    )


def build_md_archive_rows(md_source: pd.DataFrame, source_root: Path) -> list[dict[str, object]]:
    required = ["PDB", "MD_run", "Simulation_context", "Raw_data_directory"]
    require_columns(md_source, required, "MD external-inventory source")
    rows: list[dict[str, object]] = []
    keys = md_source.loc[:, required].drop_duplicates()
    for source_row in keys.to_dict("records"):
        pdb_id = str(source_row["PDB"]).upper()
        run = str(source_row["MD_run"])
        context = str(source_row["Simulation_context"])
        internal = resolve_internal_reference(source_row["Raw_data_directory"], source_root)
        run_files: list[Path] = []
        if internal is not None and internal.exists():
            run_root = internal.parent if internal.name.lower() == "raw-data" else internal
            run_files = expand_files([run_root])
        for role, suffixes in MD_ARCHIVE_CLASSES.items():
            matching = [path for path in run_files if path.suffix.lower() in suffixes]
            local, byte_count, digest, access, completeness = availability_fields(
                matching, source_root
            )
            role_slug = slug(role)
            stable_base = f"{pdb_id.lower()}-{slug(context)}-{slug(run)}"
            rows.append(
                {
                    "asset_id": f"md-{stable_base}-{role_slug}",
                    "pdb_id_or_scope": pdb_id,
                    "scientific_role": role,
                    "archive_filename_or_label": (
                        f"md/{pdb_id.lower()}/{slug(context)}-{slug(run)}/{role_slug}"
                    ),
                    "local_availability": local,
                    "repository": "not_archived",
                    "persistent_id": "",
                    "bytes": byte_count,
                    "sha256": digest,
                    "access_status": access,
                    "completeness": completeness,
                    "license_status": (
                        "manifest_only; team publication authorization confirmed; "
                        "licensed software artifacts excluded"
                    ),
                    "source_version_or_date": UNAVAILABLE,
                    "notes": (
                        "Stable inventory label only; analysis images and summary exports "
                        "do not establish availability of this artifact class."
                    ),
                }
            )
    return rows


def pocket_license_status(tool: str) -> str:
    if tool == "SiteMap":
        return "manifest_only under Schrodinger terms; binary formats and logs excluded"
    if tool in {"DoGSiteScorer", "DoGSite3"}:
        return "manifest_only; ProteinsPlus attribution retained; no raw-package relicensing"
    if tool == "CASTpFold":
        return "manifest_only; CASTpFold citation terms retained; no raw-package relicensing"
    if tool == "CavityPlus":
        return "manifest_only; CavityPlus citation terms retained; no raw-package relicensing"
    return "manifest_only; provider terms retained; no raw-package relicensing"


def build_raw_pocket_rows(
    feature_source: pd.DataFrame, source_root: Path
) -> list[dict[str, object]]:
    required = ["pdb_id", "tool", *POCKET_PATH_FIELDS]
    require_columns(feature_source, required, "pocket external-inventory source")
    rows: list[dict[str, object]] = []
    for (pdb_id, tool), group in feature_source.groupby(["pdb_id", "tool"], sort=True):
        discovered: list[Path] = []
        for field in POCKET_PATH_FIELDS:
            for raw in group[field].drop_duplicates().tolist():
                candidate = resolve_internal_reference(raw, source_root)
                if candidate is not None:
                    discovered.append(candidate)
        files = expand_files(discovered)
        local, byte_count, digest, access, completeness = availability_fields(
            files, source_root
        )
        pdb_id = str(pdb_id).upper()
        tool = str(tool)
        rows.append(
            {
                "asset_id": f"raw-pocket-{pdb_id.lower()}-{slug(tool)}",
                "pdb_id_or_scope": pdb_id,
                "scientific_role": "large raw pocket outputs",
                "archive_filename_or_label": (
                    f"raw-pocket-output/{pdb_id.lower()}/{slug(tool)}"
                ),
                "local_availability": local,
                "repository": "not_archived",
                "persistent_id": "",
                "bytes": byte_count,
                "sha256": digest,
                "access_status": access,
                "completeness": completeness,
                "license_status": pocket_license_status(tool),
                "source_version_or_date": UNAVAILABLE,
                "notes": (
                    "Stable collection label and aggregate fingerprint only; original paths, "
                    "filenames, service identifiers, and raw packages are not released."
                ),
            }
        )
    return rows


def build_external_inventory(
    feature_source: pd.DataFrame, md_source: pd.DataFrame, source_root: Path
) -> pd.DataFrame:
    rows = build_md_archive_rows(md_source, source_root)
    rows.extend(build_raw_pocket_rows(feature_source, source_root))
    output = pd.DataFrame(rows, columns=EXTERNAL_COLUMNS)
    if output.empty or not output["asset_id"].is_unique:
        raise ValueError("external archive inventory must have unique stable asset IDs")
    return output
