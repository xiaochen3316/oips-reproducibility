"""Fail closed when the tracked publication payload is unsafe or undeclared."""
from __future__ import annotations

import argparse
import ast
import csv
from hashlib import sha256
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tokenize
from typing import Mapping


MAX_FILE_BYTES = 10 * 1024 * 1024
MAX_REPOSITORY_BYTES = 25 * 1024 * 1024
HASH_RE = re.compile(r"[0-9a-f]{64}")
PROHIBITED_SUFFIXES = frozenset({
    ".7z", ".bin", ".chk", ".ckpt", ".cms", ".dcd", ".dll", ".dtr",
    ".docx", ".dylib", ".env", ".exe", ".gz", ".joblib", ".key",
    ".log", ".mae", ".maegz", ".ndjson", ".npz", ".p12", ".pdf",
    ".pem", ".pfx", ".pickle", ".pkl", ".pyc", ".rar", ".smap",
    ".so", ".tar", ".trr", ".vis", ".xlsx", ".xtc", ".zip",
})
PROHIBITED_NAMES = frozenset({".env", ".ds_store", "id_dsa", "id_ed25519", "id_rsa"})
SENSITIVE_PATTERNS: Mapping[str, re.Pattern[str]] = {
    "credential_assignment": re.compile(
        r"(?i)\b(?:api[_-]?key|access[_-]?token|auth[_-]?token|token|password|"
        r"passwd|client[_-]?secret|secret)\b\s*[:=]\s*[\"']?[^\s\"']+"
    ),
    "drive_path": re.compile(r"(?i)(?<![A-Za-z0-9_\\])[A-Z]:[\\/][^\s\"']+"),
    "file_uri": re.compile(r"(?i)\bfile://[^\s\"']+"),
    "private_posix_path": re.compile(
        r"(?i)(?<![:A-Za-z0-9_])/(?:home|Users|private|tmp|var|opt|mnt|srv|root)/[^\s\"']+"
    ),
    "unc_path": re.compile(
        r"(?i)(?<![\\])\\\\[A-Za-z0-9_$-][A-Za-z0-9._$-]+\\[^\s\"']+"
    ),
    "session_url": re.compile(
        r"(?i)https?://[^\s<>\"']*/(?:session|job|result|token)"
        r"(?:[/\?#]|$)[^\s<>\"']*"
    ),
    "session_query": re.compile(
        r"(?i)https?://[^\s<>\"']*[?&](?:job|job_id|session|session_id|"
        r"result|result_id|token)=[^&\s<>\"']+"
    ),
    "private_identifier": re.compile(
        r"(?i)\b(?:job|session|result)_?id\s*[:=]\s*[\"']?[A-Za-z0-9._-]+"
    ),
}
CONTENT_EXEMPT_PREFIXES = ("docs/superpowers/",)
RULE_DEFINITION_FILES = frozenset({
    "scripts/release_check.py",
    "src/oips_repro/provenance.py",
})
ATTACK_TEST_NAME_PARTS = (
    "attack", "sensitive", "leak", "internal_path", "redact", "private",
    "traversal", "symlink", "junction", "scanner", "adversarial",
    "relative_inputs",
)
MANIFEST_COLUMNS = (
    "asset_id", "scientific_role", "local_path_or_url", "repository",
    "persistent_id", "sha256", "bytes", "media_type", "license_or_terms",
    "access_status", "source_version_or_date",
)
RELEASE_KEYS = frozenset({
    "schema_version", "generated_at_utc", "version", "release_status", "git",
    "release_tag", "identifiers", "configuration", "inputs", "environment",
    "key_results", "publication_payload",
})
RELEASE_INVENTORIES = ("inputs", "key_results", "publication_payload")
RELEASE_PAYLOAD_PREFIXES = ("results/reference/", "figures/manuscript/")
RELEASE_AUXILIARY_INPUTS = frozenset({
    "config/schema.json", "config/figure_contract.yaml",
    "tests/scientific/data/expected_summary.json",
})


def _failure(check: str, path: str) -> dict[str, str]:
    """Return a non-secret failure record suitable for public CI logs."""
    return {"check": check, "path": path}


def _safe_relative(value: str) -> str:
    if not value or "\\" in value:
        raise ValueError("path is not relative POSIX")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or "." in path.parts:
        raise ValueError("path is not contained")
    if re.match(r"^[A-Za-z]:", value):
        raise ValueError("drive path is forbidden")
    return path.as_posix()


def _is_link_or_reparse(path: Path) -> bool:
    try:
        status = path.lstat()
    except OSError:
        return True
    return path.is_symlink() or bool(getattr(status, "st_file_attributes", 0) & 0x400)


def _reject_link_chain(root: Path, path: Path) -> None:
    try:
        relative = path.absolute().relative_to(root.absolute())
    except ValueError as error:
        raise ValueError("tracked path escapes repository") from error
    current = root.absolute()
    for part in ((), *relative.parts):
        if part:
            current = current / part
        if os.path.lexists(current) and _is_link_or_reparse(current):
            raise ValueError("tracked path traverses a symlink or junction")


def _digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def _tracked_paths(root: Path) -> list[str]:
    completed = subprocess.run(
        ["git", "ls-files", "-z", "--cached"], cwd=root, check=True,
        capture_output=True,
    )
    values = completed.stdout.decode("utf-8").split("\0")
    paths = sorted(_safe_relative(value) for value in values if value)
    if len(set(paths)) != len(paths):
        raise ValueError("tracked paths are not unique")
    return paths


def _walk_files(directory: Path) -> tuple[set[Path], set[Path]]:
    """Walk without following links, returning regular files and unsafe entries."""
    files: set[Path] = set()
    unsafe: set[Path] = set()
    pending = [directory]
    while pending:
        current = pending.pop()
        try:
            entries = list(os.scandir(current))
        except OSError:
            unsafe.add(current)
            continue
        for entry in entries:
            path = Path(entry.path)
            if entry.is_symlink() or _is_link_or_reparse(path):
                unsafe.add(path)
            elif entry.is_dir(follow_symlinks=False):
                pending.append(path)
            elif entry.is_file(follow_symlinks=False):
                files.add(path)
            else:
                unsafe.add(path)
    return files, unsafe


def _has_sensitive_text(value: str) -> bool:
    return (
        any(pattern.search(value) for pattern in SENSITIVE_PATTERNS.values())
        or bool(re.fullmatch(r"(?i)[A-Z]:[\\/]", value))
        or value.lower() == "file://"
    )


def _literal_exemption_lines(relative: str, content: str) -> set[int]:
    try:
        tree = ast.parse(content)
    except SyntaxError:
        return set()
    lines: set[int] = set()
    if relative.startswith("tests/"):
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and any(
                part in node.name.lower() for part in ATTACK_TEST_NAME_PARTS
            ):
                start, end = node.lineno, getattr(node, "end_lineno", node.lineno)
                lines.update(range(start, end + 1))
                for decorator in node.decorator_list:
                    lines.update(range(
                        decorator.lineno,
                        getattr(decorator, "end_lineno", decorator.lineno) + 1,
                    ))
    elif relative in RULE_DEFINITION_FILES:
        for node in ast.walk(tree):
            names: set[str] = set()
            if isinstance(node, ast.Assign):
                names = {
                    target.id for target in node.targets if isinstance(target, ast.Name)
                }
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                names = {node.target.id}
            if names & {"SENSITIVE_PATTERNS", "_SENSITIVE_PATTERNS"}:
                lines.update(range(
                    node.lineno, getattr(node, "end_lineno", node.lineno) + 1,
                ))
    return lines


def _mask_attack_literals(relative: str, content: str) -> str:
    """Mask sensitive literals only inside explicit attack/rule AST regions."""
    exempt_lines = _literal_exemption_lines(relative, content)
    if not exempt_lines:
        return content
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(content).readline))
    except (IndentationError, tokenize.TokenError):
        return content
    masked: list[tokenize.TokenInfo] = []
    for token_info in tokens:
        if token_info.type == tokenize.STRING and token_info.start[0] in exempt_lines:
            try:
                value = ast.literal_eval(token_info.string)
            except (SyntaxError, ValueError):
                value = None
            if isinstance(value, str) and _has_sensitive_text(value):
                token_info = tokenize.TokenInfo(
                    token_info.type, "''", token_info.start, token_info.end,
                    token_info.line,
                )
        masked.append(token_info)
    try:
        return tokenize.untokenize(masked)
    except (IndentationError, tokenize.TokenError):
        return content


def _content_for_scan(relative: str, content: str) -> str | None:
    if relative.startswith(CONTENT_EXEMPT_PREFIXES):
        return None
    if relative in RULE_DEFINITION_FILES or relative.startswith("tests/"):
        return _mask_attack_literals(relative, content)
    return content


def _scan_tracked_file(root: Path, relative: str) -> list[dict[str, str]]:
    failures: list[dict[str, str]] = []
    path = root / Path(*PurePosixPath(relative).parts)
    try:
        _reject_link_chain(root, path)
    except ValueError:
        return [_failure("unsafe_path_chain", relative)]
    if not path.is_file() or _is_link_or_reparse(path):
        return [_failure("not_regular_file", relative)]
    size = path.stat().st_size
    if size > MAX_FILE_BYTES:
        failures.append(_failure("oversized_file", relative))
    name = path.name.lower()
    suffixes = {suffix.lower() for suffix in path.suffixes}
    if name in PROHIBITED_NAMES or suffixes & PROHIBITED_SUFFIXES:
        failures.append(_failure("prohibited_extension", relative))
    try:
        content = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return failures
    content = _content_for_scan(relative, content)
    if content is None:
        return failures
    for label, pattern in SENSITIVE_PATTERNS.items():
        if pattern.search(content):
            failures.append(_failure(label, relative))
    return failures


def _parse_checksums(path: Path, base: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    declared: dict[str, str] = {}
    failures: list[dict[str, str]] = []
    label = path.relative_to(base).as_posix() if path.is_relative_to(base) else path.name
    previous = ""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return declared, [_failure("missing_checksum_manifest", label)]
    for line in lines:
        match = re.fullmatch(r"([0-9a-f]{64})  (\S+)", line)
        if not match:
            failures.append(_failure("malformed_checksum_manifest", label))
            continue
        digest, relative = match.groups()
        try:
            relative = _safe_relative(relative)
        except ValueError:
            failures.append(_failure("unsafe_checksum_path", label))
            continue
        if relative <= previous or relative in declared:
            failures.append(_failure("unsorted_checksum_manifest", label))
        previous = relative
        declared[relative] = digest
    return declared, failures


def _validate_tabular_manifest(
    base: Path, manifest: Path, checksum_paths: set[str],
) -> list[dict[str, str]]:
    relative_manifest = manifest.relative_to(base).as_posix()
    failures: list[dict[str, str]] = []
    try:
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
                return [_failure("manifest_schema", relative_manifest)]
            rows = list(reader)
    except (OSError, csv.Error):
        return [_failure("malformed_data_manifest", relative_manifest)]
    row_paths: set[str] = set()
    asset_ids: set[str] = set()
    for row in rows:
        try:
            relative = _safe_relative(row["local_path_or_url"])
        except (KeyError, ValueError):
            failures.append(_failure("unsafe_manifest_path", relative_manifest))
            continue
        if relative in row_paths or row.get("asset_id", "") in asset_ids:
            failures.append(_failure("duplicate_manifest_key", relative_manifest))
        row_paths.add(relative)
        asset_ids.add(row.get("asset_id", ""))
        payload = base / Path(*PurePosixPath(relative).parts)
        try:
            size = int(row["bytes"])
        except (KeyError, TypeError, ValueError):
            failures.append(_failure("manifest_byte_count", relative_manifest))
            continue
        digest = row.get("sha256", "")
        if (
            not HASH_RE.fullmatch(digest) or not payload.is_file()
            or _is_link_or_reparse(payload) or payload.stat().st_size != size
            or _digest(payload) != digest
        ):
            failures.append(_failure("manifest_payload_mismatch", relative))
    expected = checksum_paths - {relative_manifest}
    if row_paths != expected:
        failures.append(_failure("manifest_payload_set", relative_manifest))
    return failures


def _validate_data_tree(
    root: Path, checksum_relative: str, tracked: set[str],
) -> list[dict[str, str]]:
    failures: list[dict[str, str]] = []
    checksum_path = root / Path(*PurePosixPath(checksum_relative).parts)
    data_dir = checksum_path.parent
    base = data_dir.parent
    declared, parse_failures = _parse_checksums(checksum_path, base)
    failures.extend(parse_failures)
    files, unsafe = _walk_files(data_dir)
    for path in sorted(unsafe):
        failures.append(_failure("unsafe_payload_entry", path.relative_to(root).as_posix()))
    checksum_from_base = checksum_path.relative_to(base).as_posix()
    actual = {
        path.relative_to(base).as_posix() for path in files if path != checksum_path
    }
    if set(declared) != actual:
        failures.append(_failure("undeclared_payload", data_dir.relative_to(root).as_posix()))
    expected_tracked = {
        (base.relative_to(root) / Path(*PurePosixPath(relative).parts)).as_posix()
        if base != root else relative
        for relative in declared
    } | {checksum_relative}
    tracked_under_data = {
        relative for relative in tracked
        if relative == data_dir.relative_to(root).as_posix()
        or relative.startswith(data_dir.relative_to(root).as_posix() + "/")
    }
    if tracked_under_data != expected_tracked:
        failures.append(_failure("untracked_or_undeclared_payload", data_dir.relative_to(root).as_posix()))
    for relative, digest in declared.items():
        payload = base / Path(*PurePosixPath(relative).parts)
        if (
            not payload.is_file() or _is_link_or_reparse(payload)
            or _digest(payload) != digest
        ):
            failures.append(_failure("checksum_mismatch", relative))
    manifest = data_dir / "manifest.tsv"
    if manifest.is_file():
        failures.extend(_validate_tabular_manifest(base, manifest, set(declared)))
    else:
        failures.append(_failure("missing_data_manifest", checksum_from_base))
    return failures


def _release_entry(
    root: Path, entry: object, tracked: set[str], label: str,
) -> tuple[str | None, list[dict[str, str]]]:
    failures: list[dict[str, str]] = []
    if not isinstance(entry, Mapping) or set(entry) != {"path", "sha256", "bytes"}:
        return None, [_failure("release_manifest_entry_schema", label)]
    try:
        relative = _safe_relative(entry["path"])
        size, digest = entry["bytes"], entry["sha256"]
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise ValueError
        if not isinstance(digest, str) or not HASH_RE.fullmatch(digest):
            raise ValueError
    except (TypeError, ValueError):
        return None, [_failure("release_manifest_entry_schema", label)]
    payload = root / Path(*PurePosixPath(relative).parts)
    try:
        _reject_link_chain(root, payload)
    except ValueError:
        failures.append(_failure("release_manifest_path_chain", relative))
    if relative not in tracked:
        failures.append(_failure("release_manifest_untracked_path", relative))
    if (
        not payload.is_file() or _is_link_or_reparse(payload)
        or payload.stat().st_size != size or _digest(payload) != digest
    ):
        failures.append(_failure("release_manifest_payload_mismatch", relative))
    return relative, failures


def _release_inventory(
    root: Path, value: object, tracked: set[str], label: str,
) -> tuple[list[str], list[dict[str, str]]]:
    if not isinstance(value, list) or not value:
        return [], [_failure("release_manifest_inventory", label)]
    paths: list[str] = []
    failures: list[dict[str, str]] = []
    for index, entry in enumerate(value):
        relative, entry_failures = _release_entry(root, entry, tracked, f"{label}[{index}]")
        failures.extend(entry_failures)
        if relative is not None:
            paths.append(relative)
    if paths != sorted(paths) or len(set(paths)) != len(paths):
        failures.append(_failure("release_manifest_inventory_order", label))
    return paths, failures


def _validate_release_manifest(
    root: Path, relative: str, tracked: set[str], *, required: bool,
    expected_inputs: set[str] | None = None,
) -> list[dict[str, str]]:
    path = root / Path(*PurePosixPath(relative).parts)
    if not path.exists():
        return [_failure("missing_release_manifest", relative)] if required else []
    if relative not in tracked:
        return [_failure("untracked_release_manifest", relative)]
    try:
        _reject_link_chain(root, path)
    except ValueError:
        return [_failure("release_manifest_path_chain", relative)]
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return [_failure("malformed_release_manifest", relative)]
    if not isinstance(document, Mapping) or set(document) != RELEASE_KEYS:
        return [_failure("release_manifest_root_schema", relative)]
    failures: list[dict[str, str]] = []
    if document["schema_version"] != 1 or isinstance(document["schema_version"], bool):
        failures.append(_failure("release_manifest_schema_version", relative))
    if not isinstance(document["version"], str) or not re.fullmatch(r"\d+\.\d+\.\d+", document["version"]):
        failures.append(_failure("release_manifest_version", relative))
    if document["release_status"] not in {"pre_publication_incomplete", "ready"}:
        failures.append(_failure("release_manifest_status", relative))
    if not isinstance(document["generated_at_utc"], str) or not re.fullmatch(
        r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", document["generated_at_utc"],
    ):
        failures.append(_failure("release_manifest_timestamp", relative))
    git = document["git"]
    if not isinstance(git, Mapping) or set(git) != {"commit", "dirty"}:
        failures.append(_failure("release_manifest_git", relative))
    elif (
        git["commit"] is not None
        and (not isinstance(git["commit"], str) or not re.fullmatch(r"[0-9a-f]{40}", git["commit"]))
    ) or not isinstance(git["dirty"], bool):
        failures.append(_failure("release_manifest_git", relative))
    tag = document["release_tag"]
    if not isinstance(tag, Mapping) or set(tag) != {"tag", "status"}:
        failures.append(_failure("release_manifest_tag", relative))
    elif (
        not isinstance(tag["tag"], str)
        or not re.fullmatch(r"v\d+\.\d+\.\d+-manuscript", tag["tag"])
        or tag["status"] not in {"not_created", "draft", "released"}
    ):
        failures.append(_failure("release_manifest_tag", relative))
    elif tag["tag"] != f"v{document['version']}-manuscript":
        failures.append(_failure("release_manifest_tag_version", relative))
    identifiers = document["identifiers"]
    identifier_keys = {"repository_url", "code_doi", "data_doi", "manuscript_doi"}
    if not isinstance(identifiers, Mapping) or set(identifiers) != identifier_keys:
        failures.append(_failure("release_manifest_identifiers", relative))
    elif any(value is not None and (not isinstance(value, str) or not value) for value in identifiers.values()):
        failures.append(_failure("release_manifest_identifiers", relative))
    if document["release_status"] == "ready" and (
        not isinstance(identifiers, Mapping) or any(value is None for value in identifiers.values())
    ):
        failures.append(_failure("release_manifest_ready_identifiers", relative))
    single_paths: dict[str, str | None] = {}
    for label in ("configuration", "environment"):
        value, entry_failures = _release_entry(root, document[label], tracked, label)
        single_paths[label] = value
        failures.extend(entry_failures)
    inventories: dict[str, list[str]] = {}
    for label in RELEASE_INVENTORIES:
        values, inventory_failures = _release_inventory(root, document[label], tracked, label)
        inventories[label] = values
        failures.extend(inventory_failures)
    if expected_inputs is not None and set(inventories.get("inputs", ())) != expected_inputs:
        failures.append(_failure("release_manifest_input_set", relative))
    actual_publication = {
        item for item in tracked if item.startswith(RELEASE_PAYLOAD_PREFIXES)
    }
    if set(inventories.get("publication_payload", ())) != actual_publication:
        failures.append(_failure("release_manifest_payload_set", relative))
    if not set(inventories.get("key_results", ())).issubset(actual_publication):
        failures.append(_failure("release_manifest_key_results", relative))
    if single_paths.get("configuration") != "config/manuscript.yaml":
        failures.append(_failure("release_manifest_configuration", relative))
    if single_paths.get("environment") != "environment/constraints.txt":
        failures.append(_failure("release_manifest_environment", relative))
    return failures


def scan_repository(
    root: str | Path, *, data_manifest: str = "data/manifest.tsv",
    release_manifest: str | None = None,
) -> dict[str, object]:
    original = Path(root).absolute()
    if os.path.lexists(original) and _is_link_or_reparse(original):
        return {"status": "fail", "checked": 0, "failures": [_failure("repository_path_chain", ".")]}
    repository = original.resolve()
    failures: list[dict[str, str]] = []
    try:
        data_manifest = _safe_relative(data_manifest)
        selected_release = _safe_relative(release_manifest) if release_manifest else "release/manifest.json"
        tracked = _tracked_paths(repository)
    except (OSError, subprocess.SubprocessError, UnicodeError, ValueError):
        return {"status": "fail", "checked": 0, "failures": [_failure("git_inventory", ".")]}
    tracked_set = set(tracked)
    total_bytes = 0
    for relative in tracked:
        failures.extend(_scan_tracked_file(repository, relative))
        path = repository / Path(*PurePosixPath(relative).parts)
        if path.is_file() and not _is_link_or_reparse(path):
            total_bytes += path.stat().st_size
    if total_bytes > MAX_REPOSITORY_BYTES:
        failures.append(_failure("repository_size", "."))
    checksum_manifests = [
        relative for relative in tracked
        if relative == "data/SHA256SUMS" or relative.endswith("/data/SHA256SUMS")
    ]
    if "data/SHA256SUMS" not in checksum_manifests:
        failures.append(_failure("missing_checksum_manifest", "data/SHA256SUMS"))
    if data_manifest not in tracked_set:
        failures.append(_failure("missing_data_manifest", data_manifest))
    if PurePosixPath(data_manifest).name != "manifest.tsv":
        failures.append(_failure("data_manifest_contract", data_manifest))
    if PurePosixPath(data_manifest).parent / "SHA256SUMS" not in {
        PurePosixPath(value) for value in checksum_manifests
    }:
        failures.append(_failure("data_manifest_without_checksums", data_manifest))
    for relative in checksum_manifests:
        failures.extend(_validate_data_tree(repository, relative, tracked_set))
    expected_release_inputs: set[str] = set(RELEASE_AUXILIARY_INPUTS)
    selected_checksum = (PurePosixPath(data_manifest).parent / "SHA256SUMS").as_posix()
    if selected_checksum in tracked_set:
        checksum_path = repository / Path(*PurePosixPath(selected_checksum).parts)
        base = checksum_path.parent.parent
        declared, _ = _parse_checksums(checksum_path, base)
        base_relative = base.relative_to(repository)
        expected_release_inputs.update(
            (base_relative / Path(*PurePosixPath(value).parts)).as_posix()
            if base != repository else value
            for value in declared
        )
        expected_release_inputs.add(selected_checksum)
    failures.extend(_validate_release_manifest(
        repository, selected_release, tracked_set, required=release_manifest is not None,
        expected_inputs=expected_release_inputs,
    ))
    unique = sorted(
        {(item["check"], item["path"]) for item in failures},
        key=lambda item: (item[1], item[0]),
    )
    result_failures = [_failure(check, path) for check, path in unique]
    return {
        "status": "fail" if result_failures else "pass",
        "checked": len(tracked),
        "failures": result_failures,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository-root", type=Path, default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--data-manifest", default="data/manifest.tsv")
    parser.add_argument("--release-manifest")
    arguments = parser.parse_args(argv)
    result = scan_repository(
        arguments.repository_root, data_manifest=arguments.data_manifest,
        release_manifest=arguments.release_manifest,
    )
    print(json.dumps(result, ensure_ascii=True, sort_keys=True, separators=(",", ":")))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())
