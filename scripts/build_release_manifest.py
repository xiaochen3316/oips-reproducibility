"""Build the deterministic pre-publication release manifest.

The manifest deliberately records no repository URL, DOI, or Git commit until
the corresponding release artifacts exist.  Inventories are derived from the
committed data checksum contract and the on-disk publication payload.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess


VERSION = "1.0.0"
PLANNED_TAG = "v1.0.0-manuscript"
INPUT_AUXILIARIES = (
    "config/figure_contract.yaml",
    "config/schema.json",
    "tests/scientific/data/expected_summary.json",
)
KEY_RESULTS = (
    "figures/manuscript/repository_summary_figure_1_candidate_landscape.svg",
    "figures/manuscript/repository_summary_figure_2_qc_and_orel_ablation.svg",
    "figures/manuscript/repository_summary_figure_3_posthoc_evidence.svg",
    "figures/manuscript/repository_summary_figure_4_representative_cases.svg",
    "results/reference/analysis/final_bootstrap_intervals.csv",
    "results/reference/analysis/final_candidate_prioritization_metrics.csv",
    "results/reference/analysis/final_family_sensitivity.csv",
    "results/reference/analysis/final_reference_mapping.csv",
    "results/reference/analysis/representative_case_results.csv",
    "results/reference/static/cluster_v2_static_rankings.csv",
)
TIMESTAMP_RE = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ")
CHECKSUM_RE = re.compile(r"([0-9a-f]{64})  (\S+)")


def _relative_path(value: str) -> str:
    if not value or "\\" in value:
        raise ValueError(f"invalid repository-relative path: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"invalid repository-relative path: {value!r}")
    if re.match(r"^[A-Za-z]:", value):
        raise ValueError(f"invalid repository-relative path: {value!r}")
    return path.as_posix()


def _repository_path(root: Path, relative: str) -> Path:
    relative = _relative_path(relative)
    path = root.joinpath(*PurePosixPath(relative).parts)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"release entry is not a regular file: {relative}")
    resolved_root = root.resolve()
    try:
        path.resolve(strict=True).relative_to(resolved_root)
    except ValueError as error:
        raise ValueError(f"release entry escapes repository: {relative}") from error
    return path


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _entry(root: Path, relative: str) -> dict[str, object]:
    relative = _relative_path(relative)
    path = _repository_path(root, relative)
    return {
        "path": relative,
        "sha256": _sha256(path),
        "bytes": path.stat().st_size,
    }


def _declared_inputs(root: Path) -> list[str]:
    checksum_relative = "data/SHA256SUMS"
    checksum_path = _repository_path(root, checksum_relative)
    declared: list[str] = []
    previous = ""
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="utf-8").splitlines(), start=1,
    ):
        match = CHECKSUM_RE.fullmatch(line)
        if match is None:
            raise ValueError(f"malformed data/SHA256SUMS line {line_number}")
        expected_digest, relative = match.groups()
        relative = _relative_path(relative)
        if relative <= previous:
            raise ValueError("data/SHA256SUMS entries must be unique and sorted")
        path = _repository_path(root, relative)
        if _sha256(path) != expected_digest:
            raise ValueError(f"data checksum mismatch: {relative}")
        declared.append(relative)
        previous = relative
    return sorted({*declared, checksum_relative, *INPUT_AUXILIARIES})


def _payload_paths(root: Path) -> list[str]:
    paths: list[str] = []
    for relative_directory in ("figures/manuscript", "results/reference"):
        directory = root.joinpath(*PurePosixPath(relative_directory).parts)
        if directory.is_symlink() or not directory.is_dir():
            raise ValueError(f"publication directory is missing: {relative_directory}")
        for path in directory.rglob("*"):
            if path.is_symlink():
                raise ValueError(f"publication payload contains a link: {path.name}")
            if path.is_file():
                relative = path.relative_to(root).as_posix()
                _repository_path(root, relative)
                paths.append(relative)
    if not paths:
        raise ValueError("publication payload is empty")
    if len(paths) != len(set(paths)):
        raise ValueError("publication payload paths are not unique")
    return sorted(paths)


def _timestamp_from_existing(path: Path) -> str | None:
    if not path.is_file() or path.is_symlink():
        return None
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    value = document.get("generated_at_utc") if isinstance(document, dict) else None
    return value if isinstance(value, str) and TIMESTAMP_RE.fullmatch(value) else None


def _format_epoch(value: str) -> str:
    try:
        epoch = int(value)
    except ValueError as error:
        raise ValueError("SOURCE_DATE_EPOCH must be an integer") from error
    return datetime.fromtimestamp(epoch, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _deterministic_timestamp(root: Path, output: Path, explicit: str | None) -> str:
    if explicit is not None:
        if not TIMESTAMP_RE.fullmatch(explicit):
            raise ValueError("--generated-at-utc must use YYYY-MM-DDTHH:MM:SSZ")
        datetime.strptime(explicit, "%Y-%m-%dT%H:%M:%SZ")
        return explicit
    source_epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if source_epoch is not None:
        return _format_epoch(source_epoch)
    existing = _timestamp_from_existing(output)
    if existing is not None:
        return existing
    completed = subprocess.run(
        ["git", "log", "-1", "--format=%ct"], cwd=root, check=True,
        capture_output=True, text=True,
    )
    return _format_epoch(completed.stdout.strip())


def build_manifest(root: Path, output: Path, generated_at_utc: str | None) -> dict[str, object]:
    root = root.resolve(strict=True)
    inputs = _declared_inputs(root)
    payload = _payload_paths(root)
    missing_key_results = sorted(set(KEY_RESULTS) - set(payload))
    if missing_key_results:
        raise ValueError(f"missing key result: {missing_key_results[0]}")
    return {
        "schema_version": 1,
        "generated_at_utc": _deterministic_timestamp(root, output, generated_at_utc),
        "version": VERSION,
        "release_status": "pre_publication_incomplete",
        "git": {"commit": None, "dirty": True},
        "release_tag": {"tag": PLANNED_TAG, "status": "not_created"},
        "identifiers": {
            "repository_url": None,
            "code_doi": None,
            "data_doi": None,
            "manuscript_doi": None,
        },
        "configuration": _entry(root, "config/manuscript.yaml"),
        "inputs": [_entry(root, path) for path in inputs],
        "environment": _entry(root, "environment/constraints.txt"),
        "key_results": [_entry(root, path) for path in KEY_RESULTS],
        "publication_payload": [_entry(root, path) for path in payload],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository-root", type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--generated-at-utc")
    arguments = parser.parse_args()
    root = arguments.repository_root.resolve(strict=True)
    output = arguments.output or root / "release" / "manifest.json"
    if not output.is_absolute():
        output = root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest(root, output, arguments.generated_at_utc)
    rendered = json.dumps(manifest, ensure_ascii=True, indent=2) + "\n"
    output.write_text(rendered, encoding="utf-8", newline="\n")
    print(output.relative_to(root).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
