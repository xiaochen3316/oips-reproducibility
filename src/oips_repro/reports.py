"""Deterministic Markdown report builders for validated OIPS results.

This module deliberately contains no filesystem I/O.  Callers validate the
objects first, pass them here, and decide where (or whether) to persist the
returned Markdown.
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from hashlib import sha256
from html import escape
import json
import math
from numbers import Real
from pathlib import PurePosixPath
import re
from typing import Any


_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_DRIVE_RE = re.compile(r"[A-Za-z]:[\\/]")
_URI_RE = re.compile(r"[A-Za-z][A-Za-z0-9+.-]*://")
_TABLE_PREVIEW_ROWS = 10


def _markdown_cell(value: object) -> str:
    """Render a scalar safely inside a Markdown table cell."""
    text = _format_scalar(value)
    text = escape(text, quote=False)
    return text.replace("|", "\\|").replace("\r", " ").replace("\n", "<br>")


def _format_scalar(value: object) -> str:
    if value is None or type(value).__name__ in {"NAType", "NaTType"}:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, Real):
        number = float(value)
        if math.isnan(number):
            return "NA"
        if math.isinf(number):
            return "Inf" if number > 0 else "-Inf"
        return format(number, ".12g")
    return str(value)


def _stable_value(value: object) -> object:
    """Convert evidence values into a deterministic JSON-compatible form."""
    if isinstance(value, Mapping):
        return {
            str(key): _stable_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_stable_value(item) for item in value]
    if isinstance(value, (set, frozenset)):
        converted = [_stable_value(item) for item in value]
        return sorted(converted, key=lambda item: json.dumps(item, sort_keys=True))
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Real):
        number = float(value)
        return number if math.isfinite(number) else _format_scalar(number)
    item_method = getattr(value, "item", None)
    if callable(item_method):
        try:
            return _stable_value(item_method())
        except (TypeError, ValueError):
            pass
    return str(value)


def _stable_json(value: object) -> str:
    return json.dumps(
        _stable_value(value), ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )


def _safe_relative_path(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a nonempty relative POSIX path")
    if "\\" in value or _DRIVE_RE.search(value) or _URI_RE.search(value):
        raise ValueError(f"{label} must be a relative POSIX path")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or "." in path.parts:
        raise ValueError(f"{label} must be a contained relative POSIX path")
    return path.as_posix()


def _require_hash(value: object, label: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def _manifest_entries(manifest: Mapping[str, object], key: str) -> list[Mapping[str, object]]:
    raw = manifest.get(key, ())
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray)):
        raise ValueError(f"manifest {key} must be a sequence")
    entries: list[Mapping[str, object]] = []
    for index, entry in enumerate(raw):
        if not isinstance(entry, Mapping):
            raise ValueError(f"manifest {key}[{index}] must be a mapping")
        entries.append(entry)
    return sorted(entries, key=lambda entry: str(entry.get("path", "")))


def _file_inventory(entries: Sequence[Mapping[str, object]], label: str) -> list[str]:
    lines = ["| Path | Bytes | SHA-256 |", "|---|---:|---|"]
    if not entries:
        lines.append("| _none_ | 0 | _none_ |")
        return lines
    for index, entry in enumerate(entries):
        path = _safe_relative_path(entry.get("path"), f"{label}[{index}].path")
        digest = _require_hash(entry.get("sha256"), f"{label}[{index}].sha256")
        size = entry.get("bytes")
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise ValueError(f"{label}[{index}].bytes must be a nonnegative integer")
        lines.append(
            f"| {_markdown_cell(path)} | {size} | `{_markdown_cell(digest)}` |"
        )
    return lines


def build_rebuild_report(manifest: Mapping[str, object]) -> str:
    """Summarize a validated run manifest as deterministic Markdown."""
    if not isinstance(manifest, Mapping):
        raise ValueError("manifest must be a mapping")
    config = manifest.get("configuration")
    if not isinstance(config, Mapping):
        raise ValueError("manifest configuration must be a mapping")
    config_path = _safe_relative_path(config.get("path"), "manifest configuration path")
    config_hash = _require_hash(config.get("sha256"), "manifest configuration hash")
    inputs = _manifest_entries(manifest, "inputs")
    outputs = _manifest_entries(manifest, "outputs")

    schema_version = manifest.get("schema_version", "NA")
    seed = manifest.get("seed", manifest.get("random_seed", "NA"))
    lines = [
        "# Rebuild report",
        "",
        "This report is derived from the validated run manifest.",
        "",
        "## Run configuration",
        "",
        "| Field | Value |",
        "|---|---|",
        f"| Manifest schema | {_markdown_cell(schema_version)} |",
        f"| Configuration | `{_markdown_cell(config_path)}` |",
        f"| Configuration SHA-256 | `{_markdown_cell(config_hash)}` |",
        f"| Random seed | {_markdown_cell(seed)} |",
        "",
        "## Consumed inputs",
        "",
        *_file_inventory(inputs, "inputs"),
        "",
        "## Produced outputs",
        "",
        *_file_inventory(outputs, "outputs"),
        "",
        "## Managed paths",
        "",
    ]
    managed = manifest.get("managed_paths", ())
    if not isinstance(managed, Sequence) or isinstance(managed, (str, bytes, bytearray)):
        raise ValueError("manifest managed_paths must be a sequence")
    managed_paths = sorted(
        _safe_relative_path(value, "managed path") for value in managed
    )
    lines.extend(
        [f"- `{_markdown_cell(path)}`" for path in managed_paths]
        or ["- _none_"]
    )

    summary = manifest.get("summary")
    if summary is not None:
        if not isinstance(summary, Mapping):
            raise ValueError("manifest summary must be a mapping")
        lines.extend(
            [
                "",
                "## Rebuild summary",
                "",
                "| Field | Value |",
                "|---|---|",
                *[
                    f"| {_markdown_cell(key)} | {_markdown_cell(value)} |"
                    for key, value in sorted(summary.items(), key=lambda pair: str(pair[0]))
                ],
            ]
        )
    return "\n".join(lines) + "\n"


def _table_source(name: str) -> str:
    raw = name if "/" in name else f"analysis/{name}"
    return _safe_relative_path(raw, "analysis table source")


def _table_csv(frame: object) -> str:
    to_csv = getattr(frame, "to_csv", None)
    if not callable(to_csv):
        raise ValueError("analysis tables must be pandas DataFrames")
    try:
        return to_csv(index=False, lineterminator="\n")
    except TypeError as exc:
        raise ValueError("analysis tables must be pandas DataFrames") from exc


def _table_preview(frame: object) -> list[str]:
    columns = [str(value) for value in getattr(frame, "columns", ())]
    if not columns:
        return ["_No columns._"]
    head = getattr(frame, "head", None)
    rows = getattr(head(_TABLE_PREVIEW_ROWS), "itertuples", None) if callable(head) else None
    if not callable(rows):
        raise ValueError("analysis tables must be pandas DataFrames")
    rendered = [
        "| " + " | ".join(_markdown_cell(column) for column in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows(index=False, name=None):
        rendered.append("| " + " | ".join(_markdown_cell(value) for value in row) + " |")
    if len(frame) == 0:
        rendered.append("| " + " | ".join("_none_" for _ in columns) + " |")
    elif len(frame) > _TABLE_PREVIEW_ROWS:
        rendered.append("")
        rendered.append(f"_Preview shows {_TABLE_PREVIEW_ROWS} of {len(frame)} rows._")
    return rendered


def _analysis_group(name: str) -> str:
    lowered = name.lower()
    if any(token in lowered for token in ("bootstrap", "sensitivity", "ablation")):
        return "Uncertainty and sensitivity"
    if any(token in lowered for token in ("mapping", "evidence", "unresolved", "top3")):
        return "Evidence mappings"
    if "representative" in lowered:
        return "Representative cases"
    if any(token in lowered for token in ("metrics", "target_level")):
        return "Headline metrics"
    return "Other validated analysis tables"


def build_analysis_report(
    tables: Mapping[str, object], manifest: Mapping[str, object] | None = None,
) -> str:
    """Render a content-addressed summary of current validated analysis tables."""
    if not isinstance(tables, Mapping) or not tables:
        raise ValueError("tables must be a nonempty mapping")
    output_hashes: dict[str, str] = {}
    if manifest is not None:
        if not isinstance(manifest, Mapping):
            raise ValueError("manifest must be a mapping")
        output_hashes = {
            _safe_relative_path(entry.get("path"), "manifest output path"):
            _require_hash(entry.get("sha256"), "manifest output hash")
            for entry in _manifest_entries(manifest, "outputs")
        }
    records: list[tuple[str, str, object, str, str]] = []
    for raw_name, frame in sorted(tables.items(), key=lambda pair: str(pair[0])):
        if not isinstance(raw_name, str) or not raw_name:
            raise ValueError("analysis table names must be nonempty strings")
        source = _table_source(raw_name)
        content_digest = sha256(_table_csv(frame).encode("utf-8")).hexdigest()
        if manifest is None:
            digest = content_digest
        else:
            if source not in output_hashes:
                raise ValueError(f"analysis table is absent from manifest outputs: {source}")
            digest = output_hashes[source]
        records.append((raw_name, source, frame, digest, content_digest))

    lines = [
        "# Analysis report",
        "",
        "All values are derived from the current validated tables; source hashes are "
        "taken from the validated output manifest when supplied.",
        "",
        "## Table inventory",
        "",
        "| Source | Rows | Columns | Source SHA-256 | Validated-table fingerprint |",
        "|---|---:|---:|---|---|",
    ]
    for _name, source, frame, digest, content_digest in records:
        lines.append(
            f"| `{_markdown_cell(source)}` | {len(frame)} | {len(frame.columns)} | `{digest}` | `{content_digest}` |"
        )

    groups = (
        "Headline metrics",
        "Evidence mappings",
        "Uncertainty and sensitivity",
        "Representative cases",
        "Other validated analysis tables",
    )
    for group in groups:
        selected = [record for record in records if _analysis_group(record[0]) == group]
        lines.extend(["", f"## {group}", ""])
        if not selected:
            lines.append("_No table supplied for this section._")
            continue
        for name, source, frame, digest, content_digest in selected:
            lines.extend(
                [
                    f"### {_markdown_cell(name)}",
                    "",
                    f"Source: `{_markdown_cell(source)}`  ",
                    f"SHA-256: `{digest}`",
                    f"Validated-table fingerprint: `{content_digest}`",
                    "",
                    *_table_preview(frame),
                    "",
                ]
            )
        while lines and lines[-1] == "":
            lines.pop()

    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "This report summarizes validated computational tables; it does not add "
            "evidence beyond the cited bundle sources.",
        ]
    )
    return "\n".join(lines) + "\n"


def _expected_specs(
    value: Mapping[str, object], prefix: tuple[str, ...] = ()
) -> list[tuple[tuple[str, ...], Mapping[str, object]]]:
    found: list[tuple[tuple[str, ...], Mapping[str, object]]] = []
    for key, item in sorted(value.items(), key=lambda pair: str(pair[0])):
        path = prefix + (str(key),)
        if isinstance(item, Mapping) and {"value", "abs_tolerance"}.issubset(item):
            found.append((path, item))
        elif isinstance(item, Mapping):
            found.extend(_expected_specs(item, path))
    return found


def _flatten_actual(value: Mapping[str, object], prefix: tuple[str, ...] = ()) -> dict[str, object]:
    flattened: dict[str, object] = {}
    for key, item in value.items():
        path = prefix + (str(key),)
        dotted = ".".join(path)
        if isinstance(item, Mapping):
            flattened.update(_flatten_actual(item, path))
        else:
            flattened[dotted] = item
    return flattened


def _finite_number(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, Real):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def build_numeric_crosscheck(
    actual: Mapping[str, object], expected: Mapping[str, object]
) -> str:
    """Compare current metrics with named absolute tolerances from a snapshot."""
    if not isinstance(actual, Mapping) or not isinstance(expected, Mapping):
        raise ValueError("actual and expected metrics must be mappings")
    specs = _expected_specs(expected)
    if not specs:
        raise ValueError("expected metrics contain no named absolute tolerances")
    actual_values = _flatten_actual(actual)
    leaf_counts = Counter(path[-1] for path, _spec in specs)
    lines = [
        "# Numeric crosscheck",
        "",
        "| Metric | Actual | Expected | Delta (actual - expected) | Tolerance | Status |",
        "|---|---:|---:|---:|---|---|",
    ]
    for path, spec in specs:
        full_name = ".".join(path)
        metric_name = path[-1] if leaf_counts[path[-1]] == 1 else full_name
        candidates = (full_name, ".".join(path[1:]), path[-1])
        actual_raw = next((actual_values[key] for key in candidates if key in actual_values), None)
        actual_number = _finite_number(actual_raw)
        expected_number = _finite_number(spec.get("value"))
        tolerance = _finite_number(spec.get("abs_tolerance"))
        if expected_number is None:
            raise ValueError(f"expected value for {full_name} must be finite numeric")
        if tolerance is None or tolerance < 0:
            raise ValueError(f"abs_tolerance for {full_name} must be finite and nonnegative")
        if actual_number is None:
            actual_text, delta_text, status = "NA", "NA", "FAIL"
        else:
            delta = actual_number - expected_number
            actual_text = _format_scalar(actual_number)
            delta_text = _format_scalar(delta)
            status = "PASS" if abs(delta) <= tolerance else "FAIL"
        lines.append(
            "| "
            + " | ".join(
                (
                    _markdown_cell(metric_name),
                    actual_text,
                    _format_scalar(expected_number),
                    delta_text,
                    f"abs_tolerance={_format_scalar(tolerance)}",
                    status,
                )
            )
            + " |"
        )
    return "\n".join(lines) + "\n"


def build_validation_report(report: object) -> str:
    """Render every validation check, including stable evidence details."""
    checks = getattr(report, "checks", None)
    if not isinstance(checks, Sequence) or isinstance(checks, (str, bytes, bytearray)):
        raise ValueError("report must expose a sequence of checks")
    statuses = [getattr(check, "status", None) for check in checks]
    if any(status not in {"pass", "fail", "warning"} for status in statuses):
        raise ValueError("validation checks have invalid statuses")
    pass_count = sum(status == "pass" for status in statuses)
    fail_count = sum(status == "fail" for status in statuses)
    warning_count = sum(status == "warning" for status in statuses)
    overall = "FAIL" if fail_count else "PASS"
    lines = [
        "# Validation report",
        "",
        f"Overall status: **{overall}**",
        "",
        f"Pass: {pass_count}  ",
        f"Fail: {fail_count}  ",
        f"Warning: {warning_count}",
        "",
        "## Checks",
        "",
        "| Check ID | Status | Summary | Evidence |",
        "|---|---|---|---|",
    ]
    for check in checks:
        check_id = getattr(check, "check_id", None)
        summary = getattr(check, "summary", None)
        evidence = getattr(check, "evidence", None)
        if not isinstance(check_id, str) or not isinstance(summary, str):
            raise ValueError("validation checks must expose check_id and summary strings")
        if not isinstance(evidence, Mapping):
            raise ValueError("validation check evidence must be a mapping")
        lines.append(
            "| "
            + " | ".join(
                (
                    _markdown_cell(check_id),
                    str(getattr(check, "status")).upper(),
                    _markdown_cell(summary),
                    _markdown_cell(_stable_json(evidence)),
                )
            )
            + " |"
        )
    return "\n".join(lines) + "\n"
