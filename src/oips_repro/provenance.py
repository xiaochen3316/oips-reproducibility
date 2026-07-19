"""Stable provenance manifests and manifest-bounded overwrite safety."""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import subprocess

import pandas as pd

from oips_repro import __version__
from .config import ManuscriptConfig
from .validation import CheckResult
from . import statistics, validation


HASH_RE = re.compile(r"[0-9a-f]{64}")
DEPENDENCIES = ("PyYAML", "jsonschema", "matplotlib", "numpy", "pandas", "scipy")
MANIFEST_KEYS = (
    "schema_version", "generated_at_utc", "package", "python", "dependencies",
    "git", "configuration", "random_seed", "inputs", "outputs", "managed_paths",
    "summary",
)
SUMMARY_KEYS = (
    "records", "mapped", "excluded", "same_tool_units", "clusters", "targets",
    "boundary_sensitive", "maximum_diameter_A",
    "maximum_formal_votes_per_cluster_tool",
)


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _relative(path: Path, root: Path, label: str) -> str:
    try:
        value = path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError as error:
        raise ValueError(f"{label} must be contained by its declared root") from error
    return _safe_relative(value, label)


def _safe_relative(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"{label} must be a nonempty relative POSIX path")
    path = PurePosixPath(value)
    if path.is_absolute() or re.match(r"^[A-Za-z]:", value) or ".." in path.parts or "." in path.parts:
        raise ValueError(f"{label} must be a contained relative POSIX path")
    return path.as_posix()


def _entry(path: Path, root: Path, label: str) -> dict[str, object]:
    _reject_link_chain(path, root, label)
    if not path.is_file() or _is_link_or_reparse(path):
        raise ValueError(f"{label} must be a regular file")
    return {
        "path": _relative(path, root, label),
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
    }


def _git_state(root: Path) -> dict[str, object]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True,
            capture_output=True, text=True, timeout=5,
        ).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "status", "--porcelain"], cwd=root, check=True,
            capture_output=True, text=True, timeout=5,
        ).stdout.strip())
        if not re.fullmatch(r"[0-9a-f]{40}", commit):
            commit = None
        return {"commit": commit, "dirty": dirty}
    except (OSError, subprocess.SubprocessError):
        return {"commit": None, "dirty": None}


def _dependency_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in sorted(DEPENDENCIES):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def build_run_manifest(
    config: ManuscriptConfig,
    bundle: str | Path,
    *,
    consumed_inputs: Sequence[str | Path],
    output_paths: Sequence[str | Path],
    managed_paths: Sequence[str],
    summary: Mapping[str, object],
    generated_at_utc: str | None = None,
) -> dict[str, object]:
    """Build a stable manifest mapping without writing it."""
    root, output_root = config.repository_root.resolve(), Path(bundle).resolve()
    generated = generated_at_utc or _utc_now()
    if not re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", generated):
        raise ValueError("generated_at_utc must be RFC3339 UTC ending in Z")
    inputs = sorted(
        (_entry(Path(path), root, "input") for path in consumed_inputs),
        key=lambda value: str(value["path"]),
    )
    outputs = sorted(
        (_entry(Path(path), output_root, "output") for path in output_paths),
        key=lambda value: str(value["path"]),
    )
    if any(entry["path"] == "run_manifest.json" for entry in outputs):
        raise ValueError("run_manifest.json must not hash itself")
    if len({entry["path"] for entry in inputs}) != len(inputs):
        raise ValueError("consumed input paths must be unique")
    if len({entry["path"] for entry in outputs}) != len(outputs):
        raise ValueError("output paths must be unique")
    managed = sorted(_safe_relative(value, "managed path") for value in managed_paths)
    if len(set(managed)) != len(managed):
        raise ValueError("managed paths must be unique")
    expected_managed = sorted([*(str(entry["path"]) for entry in outputs), "run_manifest.json"])
    if managed != expected_managed:
        raise ValueError("managed paths must exactly equal outputs plus run_manifest.json")
    _validate_summary(summary)
    statistics = config.data.get("statistics")
    if not isinstance(statistics, Mapping):
        raise ValueError("configuration statistics section is missing")
    manifest: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": generated,
        "package": {"name": "oips-repro", "version": __version__},
        "python": {"version": platform.python_version()},
        "dependencies": _dependency_versions(),
        "git": _git_state(root),
        "configuration": {
            "path": _relative(config.path, root, "configuration"),
            "sha256": sha256_file(config.path),
        },
        "random_seed": statistics["random_seed"],
        "inputs": inputs,
        "outputs": outputs,
        "managed_paths": managed,
        "summary": dict(summary),
    }
    return manifest


def write_run_manifest(bundle: str | Path, manifest: Mapping[str, object]) -> Path:
    destination = Path(bundle) / "run_manifest.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(dict(manifest), ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    destination.write_text(text, encoding="utf-8", newline="\n")
    return destination


def _is_link_or_reparse(path: Path) -> bool:
    try:
        stat = path.lstat()
    except OSError:
        return True
    return path.is_symlink() or bool(getattr(stat, "st_file_attributes", 0) & 0x400)


def _reject_link_chain(path: Path, root: Path, label: str) -> None:
    base, candidate = Path(root).absolute(), Path(path).absolute()
    try:
        relative = candidate.relative_to(base)
    except ValueError as error:
        raise ValueError(f"{label} must be contained by its declared root") from error
    current = base
    for part in ((), *relative.parts):
        if part:
            current = current / part
        if os.path.lexists(current) and _is_link_or_reparse(current):
            raise ValueError(f"{label} may not traverse a symlink or junction")


def _validate_entry(entry: object, label: str) -> tuple[str, str, int]:
    if not isinstance(entry, Mapping) or set(entry) != {"path", "sha256", "bytes"}:
        raise ValueError(f"{label} must contain exactly path, sha256 and bytes")
    path = _safe_relative(entry["path"], f"{label}.path")
    digest = entry["sha256"]
    size = entry["bytes"]
    if not isinstance(digest, str) or not HASH_RE.fullmatch(digest):
        raise ValueError(f"{label}.sha256 must be lowercase 64-hex")
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise ValueError(f"{label}.bytes must be a nonnegative integer")
    return path, digest, size


def _validate_summary(summary: object) -> None:
    if not isinstance(summary, Mapping) or set(summary) != set(SUMMARY_KEYS):
        raise ValueError("run summary has incorrect keys")
    for key in set(SUMMARY_KEYS) - {"maximum_diameter_A"}:
        value = summary[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"run summary {key} must be a nonnegative integer")
    diameter = summary["maximum_diameter_A"]
    if isinstance(diameter, bool) or not isinstance(diameter, (int, float)) or not __import__("math").isfinite(float(diameter)) or diameter < 0:
        raise ValueError("run summary maximum_diameter_A must be finite and nonnegative")
    if summary["mapped"] + summary["excluded"] != summary["records"]:
        raise ValueError("run summary mapped plus excluded must equal records")


def _walk_bundle(root: Path) -> tuple[set[str], set[str]]:
    files: set[str] = set()
    directories: set[str] = set()
    pending = [root]
    while pending:
        current = pending.pop()
        with os.scandir(current) as entries:
            for entry in entries:
                path = Path(entry.path)
                relative = path.relative_to(root).as_posix()
                if entry.is_symlink() or _is_link_or_reparse(path):
                    raise ValueError(f"bundle contains a symlink or junction: {relative}")
                if entry.is_dir(follow_symlinks=False):
                    directories.add(relative); pending.append(path)
                elif entry.is_file(follow_symlinks=False):
                    files.add(relative)
                else:
                    raise ValueError(f"bundle contains a non-regular entry: {relative}")
    return files, directories


def validate_run_manifest(bundle: str | Path) -> dict[str, object]:
    """Validate a complete old bundle before any forced replacement."""
    original = Path(bundle).absolute()
    if os.path.lexists(original) and _is_link_or_reparse(original):
        raise ValueError("bundle root may not be a symlink or junction")
    root = original.resolve()
    manifest_path = root / "run_manifest.json"
    if not root.is_dir() or not manifest_path.is_file() or _is_link_or_reparse(manifest_path):
        raise ValueError("bundle must contain a regular run_manifest.json")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("run manifest is malformed") from error
    if not isinstance(manifest, dict) or set(manifest) != set(MANIFEST_KEYS):
        raise ValueError("run manifest has incorrect root keys")
    if manifest["schema_version"] != 1 or isinstance(manifest["schema_version"], bool):
        raise ValueError("run manifest schema_version must be 1")
    if not isinstance(manifest["generated_at_utc"], str) or not re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", manifest["generated_at_utc"]):
        raise ValueError("run manifest generated_at_utc is invalid")
    configuration = manifest["configuration"]
    if not isinstance(configuration, Mapping) or set(configuration) != {"path", "sha256"}:
        raise ValueError("run manifest configuration is invalid")
    _safe_relative(configuration["path"], "configuration.path")
    if not isinstance(configuration["sha256"], str) or not HASH_RE.fullmatch(configuration["sha256"]):
        raise ValueError("configuration.sha256 is invalid")
    package, python_data, dependencies, git = (
        manifest["package"], manifest["python"], manifest["dependencies"], manifest["git"]
    )
    if not isinstance(package, Mapping) or set(package) != {"name", "version"} or package["name"] != "oips-repro" or not isinstance(package["version"], str):
        raise ValueError("run manifest package schema is invalid")
    if not isinstance(python_data, Mapping) or set(python_data) != {"version"} or not isinstance(python_data["version"], str):
        raise ValueError("run manifest python schema is invalid")
    if not isinstance(dependencies, Mapping) or list(dependencies) != sorted(DEPENDENCIES) or any(not isinstance(value, str) or not value for value in dependencies.values()):
        raise ValueError("run manifest dependencies schema is invalid")
    if not isinstance(git, Mapping) or set(git) != {"commit", "dirty"}:
        raise ValueError("run manifest git schema is invalid")
    if git["commit"] is not None and (not isinstance(git["commit"], str) or not re.fullmatch(r"[0-9a-f]{40}", git["commit"])):
        raise ValueError("run manifest git commit is invalid")
    if git["dirty"] is not None and not isinstance(git["dirty"], bool):
        raise ValueError("run manifest git dirty flag is invalid")
    if isinstance(manifest["random_seed"], bool) or not isinstance(manifest["random_seed"], int):
        raise ValueError("run manifest random_seed must be an integer")
    inputs_raw, outputs_raw = manifest["inputs"], manifest["outputs"]
    managed_raw = manifest["managed_paths"]
    if not isinstance(inputs_raw, list) or not isinstance(outputs_raw, list) or not isinstance(managed_raw, list):
        raise ValueError("run manifest inventories must be lists")
    inputs = [_validate_entry(value, f"inputs[{index}]") for index, value in enumerate(inputs_raw)]
    outputs = [_validate_entry(value, f"outputs[{index}]") for index, value in enumerate(outputs_raw)]
    if inputs != sorted(inputs) or outputs != sorted(outputs):
        raise ValueError("run manifest inventories must be sorted")
    managed = [_safe_relative(value, "managed path") for value in managed_raw]
    if managed != sorted(managed) or len(set(managed)) != len(managed):
        raise ValueError("run manifest managed paths must be sorted and unique")
    output_names = [value[0] for value in outputs]
    if managed != sorted([*output_names, "run_manifest.json"]) or len(set(output_names)) != len(output_names):
        raise ValueError("managed paths must exactly equal unique outputs plus run_manifest.json")
    _validate_summary(manifest["summary"])
    files, directories = _walk_bundle(root)
    allowed_files = set(output_names) | {"run_manifest.json", ".gitkeep"}
    if files != allowed_files and files != allowed_files - {".gitkeep"}:
        raise ValueError("bundle contains undeclared or missing output files")
    allowed_dirs = {
        PurePosixPath(path).parents[index].as_posix()
        for path in output_names
        for index in range(len(PurePosixPath(path).parents) - 1)
    }
    allowed_dirs.discard(".")
    if directories != allowed_dirs:
        raise ValueError("bundle contains undeclared or missing directories")
    for relative, digest, size in outputs:
        path = root / Path(*PurePosixPath(relative).parts)
        if not path.is_file() or _is_link_or_reparse(path):
            raise ValueError(f"manifest output is not a regular file: {relative}")
        if path.stat().st_size != size or sha256_file(path) != digest:
            raise ValueError(f"manifest output hash or byte count mismatch: {relative}")
    return manifest


def prepare_output_bundle(
    output: str | Path, reference: str | Path, *, force: bool = False,
) -> Path:
    """Preflight and, only after full validation, clear a managed old bundle."""
    from .io import ensure_safe_output_path

    original = Path(output).absolute()
    _reject_link_chain(original, Path(original.anchor), "output path")
    if os.path.lexists(original) and _is_link_or_reparse(original):
        raise ValueError("output root may not be a symlink or junction")
    destination = ensure_safe_output_path(original, reference)
    if not destination.exists():
        destination.mkdir(parents=True)
        return destination
    if not destination.is_dir() or _is_link_or_reparse(destination):
        raise ValueError("output must be a regular directory")
    entries = list(destination.iterdir())
    if any(_is_link_or_reparse(entry) for entry in entries):
        raise ValueError("output entries may not be symlinks or junctions")
    if not entries or all(entry.name == ".gitkeep" and entry.is_file() for entry in entries):
        return destination
    if not force:
        raise ValueError("output bundle already exists; use --force after validating its manifest")
    manifest = validate_run_manifest(destination)
    managed = [str(value) for value in manifest["managed_paths"]]
    for relative in managed:
        path = destination / Path(*PurePosixPath(relative).parts)
        path.unlink()
    directories = sorted({
        destination / Path(*parent.parts)
        for relative in managed for parent in PurePosixPath(relative).parents
        if parent.as_posix() != "."
    }, key=lambda value: len(value.parts), reverse=True)
    for directory in directories:
        if directory.exists() and not any(directory.iterdir()):
            directory.rmdir()
    return destination


def validate_manifest_context(
    config: ManuscriptConfig, manifest: Mapping[str, object], consumed_inputs: Sequence[str | Path],
) -> None:
    expected_configuration = {
        "path": _relative(config.path, config.repository_root, "configuration"),
        "sha256": sha256_file(config.path),
    }
    if manifest.get("configuration") != expected_configuration:
        raise ValueError("manifest configuration does not match the current configuration")
    statistics = config.data.get("statistics")
    if not isinstance(statistics, Mapping) or manifest.get("random_seed") != statistics.get("random_seed"):
        raise ValueError("manifest random seed does not match the current configuration")
    expected_inputs = sorted(
        (_entry(Path(path), config.repository_root, "input") for path in consumed_inputs),
        key=lambda value: str(value["path"]),
    )
    if manifest.get("inputs") != expected_inputs:
        raise ValueError("manifest inputs do not match current consumed inputs")


def collect_verification(
    config: ManuscriptConfig, bundle: Path, snapshot_path: Path,
    stage_validator: Callable[[], Mapping[str, pd.DataFrame]], *, immutable: bool,
) -> tuple[object, object, validation.ValidationReport, validation.ValidationReport, object, bool]:
    """Aggregate independent manifest, stage, input and snapshot validation."""
    state: dict[str, object] = {}
    manifest_path, present = bundle / "run_manifest.json", (bundle / "run_manifest.json").is_file()

    def manifest_action():
        manifest = validate_run_manifest(bundle)
        if manifest["summary"] != cluster_summary(bundle):
            raise ValueError("manifest summary does not match clustering outputs")
        state["manifest"] = manifest
        return {}

    def context_action():
        manifest = state.get("manifest")
        if manifest is None:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        validate_manifest_context(config, manifest, validation.consumed_input_paths(config))
        return {}

    if present:
        manifest_checks = (
            validation._result("manifest_integrity", manifest_action),
            validation._result("manifest_context", context_action),
        )
    else:
        status = "fail" if immutable else "warning"
        summary = "immutable reference manifest is missing" if immutable else "stage-only bundle has no manifest"
        manifest_checks = tuple(
            validation.CheckResult(name, status, summary, {})
            for name in ("manifest_integrity", "manifest_context")
        )

    def stage_action():
        state["frames"] = stage_validator()
        return {}

    bundle_check = validation._result("bundle_contract", stage_action)
    if bundle_check.status == "pass":
        bundle_check = validation.CheckResult("bundle_contract", "pass", "stage contracts validated", {})
    try:
        snapshot_report = validation.validate_snapshot(bundle, snapshot_path)
    except Exception as error:
        snapshot_report = validation.ValidationReport((validation.CheckResult(
            "snapshot_validation", "fail", type(error).__name__,
            {"reason": validation._safe_reason(error)},
        ),))
    input_checks = validation.validate_public_inputs(config).checks
    content = validation.ValidationReport((*input_checks, bundle_check, *snapshot_report.checks))
    aggregate = validation.ValidationReport((
        *input_checks, *manifest_checks, bundle_check, *snapshot_report.checks,
    ))
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8")) if not aggregate.failure_count else None
    return state.get("manifest"), state.get("frames"), aggregate, content, snapshot, present


def _counts(series: pd.Series) -> dict[str, int]:
    return {str(key): int(value) for key, value in series.value_counts().to_dict().items()}


def cluster_summary(bundle: str | Path) -> dict[str, object]:
    root = Path(bundle)
    candidates = pd.read_csv(root / "clustering" / "cluster_v2_candidates.csv")
    membership = pd.read_csv(root / "clustering" / "cluster_v2_membership.csv")
    mapping = pd.read_csv(root / "clustering" / "tool_record_to_cluster_v2_mapping.csv")
    excluded = pd.read_csv(root / "clustering" / "excluded_unmappable_records.csv")
    votes = membership.groupby(["cluster_v2_id", "tool"])["formal_vote_count"].sum()
    return {
        "records": len(mapping), "mapped": int(mapping["cluster_v2_id"].notna().sum()),
        "excluded": len(excluded), "same_tool_units": int(membership["same_tool_unit_id"].nunique()),
        "clusters": len(candidates), "targets": int(candidates["pdb_id"].nunique()),
        "boundary_sensitive": int(candidates["boundary_sensitive"].sum()),
        "maximum_diameter_A": float(candidates["cluster_diameter_A"].max()),
        "maximum_formal_votes_per_cluster_tool": int(votes.max()),
    }


def snapshot_actual(
    bundle: str | Path, expected: Mapping[str, object] | None = None,
) -> dict[str, object]:
    root, analysis = Path(bundle), Path(bundle) / "analysis"
    candidates = pd.read_csv(root / "clustering" / "cluster_v2_candidates.csv")
    master = pd.read_csv(root / "static" / "cluster_v2_master_table.csv")
    rankings = pd.read_csv(root / "static" / "cluster_v2_static_rankings.csv")
    metrics = pd.read_csv(analysis / "final_candidate_prioritization_metrics.csv").iloc[0]
    bootstrap = pd.read_csv(analysis / "final_bootstrap_intervals.csv")
    family = bootstrap.loc[bootstrap["bootstrap_method"].eq("family_clustered")]
    qc = pd.read_csv(analysis / "final_top3_automated_QC.csv")
    labels = pd.read_csv(analysis / "final_automated_evidence_labels.csv")
    reference = pd.read_csv(analysis / "final_reference_mapping.csv")
    md = pd.read_csv(analysis / "final_md_cluster_v2_mapping.csv", keep_default_na=False)
    redocking = pd.read_csv(analysis / "final_redocking_cluster_v2_mapping.csv")
    unresolved = pd.read_csv(analysis / "unresolved_cases.csv")
    ablation = pd.read_csv(analysis / "final_orel_ablation_targets.csv")
    weight_scenarios = pd.read_csv(analysis / "weight_sensitivity_scenarios.csv")
    weight_targets = pd.read_csv(analysis / "weight_sensitivity_targets.csv")
    single_targets = pd.read_csv(analysis / "single_tool_target_ranks.csv")
    single_metrics = pd.read_csv(analysis / "single_tool_complete_case_metrics.csv")
    maximum = candidates.loc[candidates["cluster_diameter_A"].idxmax()]
    result: dict[str, object] = {
        "cluster_numeric": {
            "maximum_diameter_A": float(maximum["cluster_diameter_A"]),
            "maximum_recomputation_difference": float((rankings["OIPS-P_static"] - rankings["OIPS-P_static_recomputed"]).abs().max()),
        },
        "cluster_distributions": {
            "mappability": _counts(candidates["mappability"]),
            "tool_support_count": _counts(candidates["tool_support_count"]),
            "center_only_clusters": int(candidates["mappability"].eq("center_only_mappable").sum()),
            "missing_ligandability": int(master["P_lig"].isna().sum()),
            "tie_count": int(rankings["tie_flag"].sum()),
            "spatial_continuity_true": int(candidates["spatial_continuity"].sum()),
            "maximum_diameter_cluster_id": str(maximum["cluster_v2_id"]),
        },
        "analysis_metrics": {name: float(metrics[name]) for name in statistics.METRIC_FIELDS},
        "bootstrap_family_intervals": {
            **{f"{row['metric']}__CI_2_5_percent": float(row["CI_2.5_percent"]) for _, row in family.iterrows()},
            **{f"{row['metric']}__CI_97_5_percent": float(row["CI_97.5_percent"]) for _, row in family.iterrows()},
        },
        "orel_ablation": {
            "Without_O_rel_Reference_Top1": float(ablation["Without_O_rel_reference_rank"].le(1).mean()),
            "Without_O_rel_Reference_Top3": float(ablation["Without_O_rel_reference_rank"].le(3).mean()),
            "Without_O_rel_Reference_MRR": float((1.0 / ablation["Without_O_rel_reference_rank"]).mean()),
            "mean_rank_Spearman_rho": float(ablation["rank_Spearman_rho"].mean()),
            "top_cluster_identity_change_fraction": float(ablation["top_cluster_identity_changed"].mean()),
        },
        "supplementary_analysis_metrics": {
            "weight_scenario_rows": float(len(weight_scenarios)),
            "weight_min_top1_retained_N": float(weight_scenarios["baseline_top1_retained_N"].min()),
            "weight_max_top1_retained_N": float(weight_scenarios["baseline_top1_retained_N"].max()),
            "weight_min_mean_top3_jaccard": float(weight_scenarios["mean_top3_jaccard"].min()),
            "weight_max_mean_top3_jaccard": float(weight_scenarios["mean_top3_jaccard"].max()),
            "weight_min_median_spearman_rho": float(weight_scenarios["median_spearman_rho"].min()),
            "weight_max_median_spearman_rho": float(weight_scenarios["median_spearman_rho"].max()),
            "weight_min_reference_top1_N": float(weight_scenarios["Reference_Top1_N"].min()),
            "weight_max_reference_top1_N": float(weight_scenarios["Reference_Top1_N"].max()),
            "weight_min_reference_top3_N": float(weight_scenarios["Reference_Top3_N"].min()),
            "weight_max_reference_top3_N": float(weight_scenarios["Reference_Top3_N"].max()),
            "weight_min_first_supported_top1_N": float(weight_scenarios["First_supported_Top1_N"].min()),
            "weight_max_first_supported_top1_N": float(weight_scenarios["First_supported_Top1_N"].max()),
            "weight_min_first_supported_top3_N": float(weight_scenarios["First_supported_Top3_N"].min()),
            "weight_max_first_supported_top3_N": float(weight_scenarios["First_supported_Top3_N"].max()),
            "weight_min_target_top1_retention": float(weight_targets["top1_retention_count"].min()),
            **{
                f"single_tool_{row['method']}_{field}": float(row[field])
                for _, row in single_metrics.iterrows()
                for field in ("N", "Top1_N", "Top3_N", "Top5_N", "MRR")
            },
        },
        "single_tool_complete_case_targets": sorted(single_targets.loc[
            single_targets["method"].eq("OIPS-P") & single_targets["complete_case"].map(
                lambda value: str(value).strip().lower() == "true"
            ), "pdb_id"
        ].astype(str).tolist()),
        "posthoc_counts": {
            "topk_qc_rows": len(qc), "reference_rows": len(reference), "label_rows": len(labels),
            "md_rows": len(md), "md_real_runs": int(md["MD_run"].astype(str).ne("").sum()),
            "md_targets_with_runs": int(md.loc[md["MD_run"].astype(str).ne(""), "pdb_id"].nunique()),
            "redocking_rows": len(redocking), "unresolved_rows": len(unresolved),
            "reference_evaluable": int(reference.loc[
                reference["reference_selection_unresolved"].map(
                    lambda value: str(value).strip().lower() == "false"
                ), "pdb_id"
            ].nunique()),
        },
        "posthoc_distributions": {
            "QC_status": _counts(qc["QC_status"]),
            "automated_evidence_label": _counts(labels["automated_evidence_label"]),
            "top3_evidence_label": _counts(labels.loc[labels["Within_PDB_rank"].le(3), "automated_evidence_label"]),
            "Concordance_call": _counts(md["Concordance_call"]),
            "reference_selection_status": _counts(reference["reference_selection_status"]),
            "redocking_RMSD_bins": {
                "le_2A": int(redocking["Raw_ligand_RMSD_A"].le(2).sum()),
                "gt_2_le_3A": int((redocking["Raw_ligand_RMSD_A"].gt(2) & redocking["Raw_ligand_RMSD_A"].le(3)).sum()),
                "gt_3A": int(redocking["Raw_ligand_RMSD_A"].gt(3).sum()),
            },
            "unresolved_issue_type": _counts(unresolved["issue_type"]),
        },
        "analysis_counts": {
            "family_sensitivity_rows": len(pd.read_csv(analysis / "final_family_sensitivity.csv")),
            "ablation_target_rows": len(ablation),
            "ablation_category_rows": len(pd.read_csv(analysis / "final_orel_ablation_categories.csv")),
            "ablation_cluster_rows": len(pd.read_csv(analysis / "orel_ablation_cluster_rankings.csv")),
        },
        "representative_case_order": pd.read_csv(analysis / "representative_case_results.csv")["pdb_id"].astype(str).tolist(),
    }
    if expected:
        exemplar_id = expected.get("exemplar_cluster_id")
        if isinstance(exemplar_id, str):
            row = rankings.set_index("cluster_v2_id").loc[exemplar_id]
            names = expected.get("exemplar_static_metrics", {})
            result["exemplar_static_metrics"] = {name: float(row[name]) for name in names}
            result["exemplar_static_state"] = {
                "cluster_v2_id": exemplar_id, "tie_flag": bool(row["tie_flag"]),
            }
        case = expected.get("posthoc_reference_case")
        if isinstance(case, Mapping) and isinstance(case.get("pdb_id"), str):
            rows = reference.loc[reference["pdb_id"].eq(case["pdb_id"])]
            result["posthoc_reference_case"] = {
                "pdb_id": case["pdb_id"],
                "selection_status": str(rows["reference_selection_status"].iloc[0]),
                "selected_ligand_key": str(rows["selected_reference_ligand_key"].iloc[0]),
                "ligand_atom_count": int(rows["reference_ligand_atom_count"].iloc[0]),
            }
    return result


_SENSITIVE_PATTERNS = {
    "credential": re.compile(r"(?i)\b(?:api[_-]?key|token|password|secret)\s*[:=]\s*\S+"),
    "file_uri": re.compile(r"(?i)\bfile://\S+"),
    "drive_path": re.compile(r"(?i)(?:^|[\s\"'])\b[A-Z]:[\\/][^\s\"']+"),
    "session_url": re.compile(r"(?i)https?://\S+/(?:session|job|result|token)(?:/|\?|$)\S*"),
}


def scan_sensitive_content(named_text: Mapping[str, str]) -> CheckResult:
    classes: set[str] = set()
    files: set[str] = set()
    count = 0
    for name, content in named_text.items():
        if not isinstance(name, str) or not isinstance(content, str):
            raise ValueError("sensitive scan inputs must map names to text")
        for category, pattern in _SENSITIVE_PATTERNS.items():
            matches = list(pattern.finditer(content))
            if matches:
                classes.add(category); files.add(name); count += len(matches)
    status = "fail" if count else "pass"
    summary = "sensitive content detected" if count else "no sensitive content detected"
    return CheckResult(
        "sensitive_content", status, summary,
        {"match_count": count, "classes": tuple(sorted(classes)), "files": tuple(sorted(files))},
    )
