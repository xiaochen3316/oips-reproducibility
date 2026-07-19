"""Command-line interface for the OIPS reproducibility workflow."""
from __future__ import annotations
import argparse
from collections.abc import Mapping, Sequence
from contextlib import contextmanager
import json
from pathlib import Path
import re
import sys
import pandas as pd
import yaml
from oips_repro import __version__
from . import clustering, io, posthoc, scoring, statistics, structure, validation
from . import provenance, release, reports
STABLE_COMMANDS = (
    "validate-data", "cluster", "score", "analyze", "figures", "verify",
    "reproduce",
)
def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="oips-repro")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command")
    for command in STABLE_COMMANDS:
        subparser = subparsers.add_parser(command)
        if command == "validate-data":
            subparser.add_argument("--config", type=Path, required=True)
        elif command == "cluster":
            subparser.add_argument("--config", type=Path, required=True)
            subparser.add_argument("--output", type=Path, required=True)
        elif command == "score":
            subparser.add_argument("--config", type=Path, required=True)
            subparser.add_argument("--cluster-dir", type=Path, required=True)
            subparser.add_argument("--output", type=Path, required=True)
        elif command == "analyze":
            subparser.add_argument("--config", type=Path, required=True)
            subparser.add_argument("--cluster-dir", type=Path, required=True)
            subparser.add_argument("--static-dir", type=Path, required=True)
            subparser.add_argument("--posthoc-data", type=Path, required=True)
            subparser.add_argument("--output", type=Path, required=True)
        elif command == "figures":
            subparser.add_argument("--config", type=Path, required=True)
            subparser.add_argument("--analysis", type=Path, required=True)
            subparser.add_argument("--output", type=Path, required=True)
        elif command == "verify":
            subparser.add_argument("--config", type=Path, required=True)
            subparser.add_argument("--bundle", type=Path, required=True)
            subparser.add_argument("--snapshot", type=Path, required=True)
        elif command == "reproduce":
            subparser.add_argument("--config", type=Path, required=True)
            subparser.add_argument("--output", type=Path)
            subparser.add_argument("--force", action="store_true")
    return parser
def _compact(payload: Mapping[str, object]) -> None:
    print(json.dumps(payload, separators=(",", ":"), sort_keys=True))
def _safe_error_message(error: BaseException) -> str:
    message = str(error).replace("\r", " ").replace("\n", " ")
    message = re.sub(r"(?i)[A-Z]:[\\/].*", "<redacted-path>", message)
    message = re.sub(r"(?i)(?:file|https?)://\S+", "<redacted-url>", message)
    message = re.sub(r"(?i)(?<![:\w])/(?:home|users|private|tmp|var|opt|mnt)/\S*", "<redacted-path>", message)
    message = re.sub(r"(?i)\b(?:token|password|secret|api[_-]?key)\s*[:=]\s*\S+", "<redacted-credential>", message)
    return f"{type(error).__name__}: {message}"
def _run_validate_data(arguments: argparse.Namespace) -> int:
    config = io.load_manuscript_config(arguments.config)
    report = validation.validate_public_inputs(config)
    failures = [check.check_id for check in report.checks if check.status == "fail"]
    _compact({"status": report.status, "checks": len(report.checks), "failures": failures})
    return 1 if failures else 0
def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"configuration section must be a mapping: {label}")
    return value
def _same_tool_config(config: io.ManuscriptConfig) -> clustering.SameToolConfig:
    cluster_data = _mapping(config.data.get("clustering"), "clustering")
    values = _mapping(cluster_data.get("same_tool"), "clustering.same_tool")
    return clustering.SameToolConfig(**{key: value for key, value in values.items()})
def _cross_tool_config(config: io.ManuscriptConfig) -> clustering.CrossToolConfig:
    cluster_data = _mapping(config.data.get("clustering"), "clustering")
    values = _mapping(cluster_data.get("cross_tool"), "clustering.cross_tool")
    return clustering.CrossToolConfig(**{key: value for key, value in values.items()})
def _tool_weights(config: io.ManuscriptConfig) -> dict[str, float]:
    values = _mapping(config.data.get("tool_weights"), "tool_weights")
    return {str(key): float(value) for key, value in values.items()}
def _scoring_config(config: io.ManuscriptConfig) -> scoring.ScoringConfig:
    values = _mapping(config.data.get("scoring"), "scoring")
    modules = _mapping(values.get("module_weights"), "scoring.module_weights")
    quality = _mapping(values.get("q_evidence"), "scoring.q_evidence")
    interface = _mapping(values.get("interface"), "scoring.interface")
    weights = _mapping(interface.get("weights"), "scoring.interface.weights")
    distance = _mapping(interface.get("distance_score"), "scoring.interface.distance_score")
    ranking = _mapping(config.data.get("ranking"), "ranking")
    return scoring.ScoringConfig(
        module_weights={str(key): float(value) for key, value in modules.items()},
        tool_weights=_tool_weights(config),
        geometry_top_n=int(values["geometry_top_n"]),
        ligandability_top_n=int(values["ligandability_top_n"]),
        q_base=float(quality["base"]),
        q_per_representative=float(quality["per_representative"]),
        q_center_fraction=float(quality["center_fraction"]),
        q_residue_fraction=float(quality["residue_fraction"]),
        q_sitemap_bonus=float(quality["sitemap_bonus"]),
        ranking_direction=str(ranking["direction"]),
        ranking_tie_method=str(ranking["tie_method"]),
        interface_fraction_weight=float(weights["interface_fraction"]),
        chain_context_weight=float(weights["chain_context"]),
        distance_weight=float(weights["distance"]),
        interface_extent_weight=float(weights["interface_extent"]),
        distance_near_max_A=float(distance["near_max_A"]),
        distance_near_score=float(distance["near_score"]),
        distance_intermediate_max_A=float(distance["intermediate_max_A"]),
        distance_intermediate_start=float(distance["intermediate_start"]),
        distance_intermediate_loss_per_A=float(distance["intermediate_loss_per_A"]),
        distance_far_max_A=float(distance["far_max_A"]),
        distance_far_start=float(distance["far_start"]),
        distance_far_loss_per_A=float(distance["far_loss_per_A"]),
        distance_tail_start=float(distance["tail_start"]),
        distance_tail_loss_per_A=float(distance["tail_loss_per_A"]),
        distance_tail_floor=float(distance["tail_floor"]),
        interface_fraction_base=float(interface["interface_fraction_base"]),
        interface_fraction_multiplier=float(interface["interface_fraction_multiplier"]),
        multi_chain_base=float(interface["multi_chain_base"]),
        chain_entropy_multiplier=float(interface["chain_entropy_multiplier"]),
        single_chain_interface_fraction_min=float(interface["single_chain_interface_fraction_min"]),
        single_chain_supported_score=float(interface["single_chain_supported_score"]),
        single_chain_other_score=float(interface["single_chain_other_score"]),
        context_base=float(interface["context_base"]),
        context_per_extra_chain=float(interface["context_per_extra_chain"]),
        context_extra_chain_cap=int(interface["context_extra_chain_cap"]),
        context_interface_fraction_cap=float(interface["context_interface_fraction_cap"]),
        context_interface_fraction_multiplier=float(interface["context_interface_fraction_multiplier"]),
        no_interface_monomer_score=float(interface["no_interface_monomer_score"]),
        no_interface_multichain_score=float(interface["no_interface_multichain_score"]),
    )
def _safe_output(config: io.ManuscriptConfig, output: Path) -> Path:
    return io.ensure_safe_output_path(
        output, config.resolve_configured_path("reference_results")
    )
def _run_cluster(arguments: argparse.Namespace) -> int:
    config = io.load_manuscript_config(arguments.config)
    output = _safe_output(config, arguments.output)
    features = io.load_feature_table(config.resolve_configured_path("feature_table"))
    result = clustering.build_clusters(
        features,
        same_tool=_same_tool_config(config),
        cross_tool=_cross_tool_config(config),
        tool_weights=_tool_weights(config),
    )
    score_data = _mapping(config.data.get("scoring"), "scoring")
    boundary = _mapping(score_data.get("boundary"), "scoring.boundary")
    result = clustering.annotate_boundaries(
        result,
        diameter_gt_A=float(boundary["diameter_gt_A"]),
        dispersion_gt_A=float(boundary["dispersion_gt_A"]),
        core_envelope_ratio_lt=float(boundary["core_envelope_ratio_lt"]),
        median_residue_iou_lt=float(boundary["median_residue_iou_lt"]),
    )
    io.write_clustering_result(result, output)
    summary = {
        "records": len(result.mapping),
        "mapped": int(result.mapping["cluster_v2_id"].notna().sum()),
        "excluded": len(result.excluded),
        "same_tool_units": len(result.units),
        "clusters": len(result.candidates),
        "targets": int(result.candidates["pdb_id"].nunique()),
        "boundary_sensitive": int(result.candidates["boundary_sensitive"].sum()),
    }
    print(json.dumps(summary, indent=2))
    return 0
def _run_score(arguments: argparse.Namespace) -> int:
    config = io.load_manuscript_config(arguments.config)
    output = _safe_output(config, arguments.output)
    features = io.load_feature_table(config.resolve_configured_path("feature_table"))
    clustered = io.load_clustering_result(arguments.cluster_dir.resolve(), features=features)
    score_config = _scoring_config(config)
    score_data = _mapping(config.data.get("scoring"), "scoring")
    interface_data = _mapping(score_data.get("interface"), "scoring.interface")
    target_ids = set(clustered.candidates["pdb_id"].astype(str))
    structure_paths = io.load_structure_paths(config, target_ids=target_ids)
    interfaces = {
        pdb_id: structure.build_interface_profile(
            structure.read_structure_atoms(path),
            contact_cutoff_A=float(interface_data["contact_cutoff_A"]),
        )
        for pdb_id, path in sorted(structure_paths.items())
    }
    result = scoring.score_clusters(clustered, interfaces=interfaces, config=score_config)
    io.write_scoring_result(result, output)
    summary = {
        "clusters": len(result.master),
        "targets": int(result.master["pdb_id"].nunique()),
        "rankings": len(result.rankings),
        "ties": int(result.rankings["tie_flag"].sum()),
        "maximum_recomputation_difference": result.maximum_recomputation_difference,
    }
    print(json.dumps(summary, indent=2))
    return 0
POSTHOC_INPUT_FILES = validation.POSTHOC_INPUT_FILES
ANALYSIS_SCHEMAS = validation.ANALYSIS_SCHEMAS
ANALYSIS_SORTS = validation.ANALYSIS_SORTS
_validate_analysis_csv = validation.validate_analysis_csv
def _directory_files(path: Path, label: str) -> set[str]:
    if not path.is_dir():
        raise FileNotFoundError(f"{label} directory does not exist: {path}")
    if any(not entry.is_file() for entry in path.iterdir()):
        raise ValueError(f"{label} directory may contain files only")
    return {entry.name for entry in path.iterdir()}
def _prepare_managed_directory(path: Path, expected: set[str], label: str) -> None:
    if path.exists() and not path.is_dir():
        raise ValueError(f"{label} output must be a directory")
    path.mkdir(parents=True, exist_ok=True)
    actual = {entry.name for entry in path.iterdir()}
    unexpected = actual - expected
    invalid_managed = [entry.name for entry in path.iterdir() if entry.name in expected and not entry.is_file()]
    if unexpected or invalid_managed:
        raise ValueError(
            f"{label} output contains unexpected entries: {sorted(unexpected | set(invalid_managed))}"
        )
def _load_static_tables(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    actual = _directory_files(path, "static result")
    if actual != set(io.SCORE_FILE_SCHEMAS):
        raise ValueError(f"static result directory file set is invalid: {sorted(actual)}")
    loaded = {
        filename: io._read_result_csv(path / filename, *schema)
        for filename, schema in io.SCORE_FILE_SCHEMAS.items()
    }
    rankings = loaded["cluster_v2_static_rankings.csv"]
    master = loaded["cluster_v2_master_table.csv"]
    if rankings.empty or rankings["pdb_id"].nunique() < 1:
        raise ValueError("analysis requires nonempty static rankings")
    if set(rankings["cluster_v2_id"]) != set(master["cluster_v2_id"]):
        raise ValueError("static master and ranking IDs do not match")
    if rankings.duplicated(["pdb_id", "cluster_v2_id"]).any():
        raise ValueError("static rankings contain duplicate cluster keys")
    return master, rankings
def _load_static_rankings(path: Path) -> pd.DataFrame:
    return _load_static_tables(path)[1]
def _load_posthoc_inputs(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    actual = _directory_files(path, "post-hoc data")
    if actual != POSTHOC_INPUT_FILES:
        raise ValueError(f"post-hoc data directory file set is invalid: {sorted(actual)}")
    return (
        pd.read_csv(path / "reference_annotations.csv", keep_default_na=False),
        pd.read_csv(path / "md_evidence.csv", keep_default_na=False),
        pd.read_csv(path / "redocking_evidence.csv", keep_default_na=False),
    )
def _statistics_config(config: io.ManuscriptConfig) -> tuple[int, int]:
    values = _mapping(config.data.get("statistics"), "statistics")
    iterations = values.get("bootstrap_iterations")
    seed = values.get("random_seed")
    if not isinstance(iterations, int) or iterations < 1:
        raise ValueError("statistics.bootstrap_iterations must be a positive integer")
    if not isinstance(seed, int):
        raise ValueError("statistics.random_seed must be an integer")
    return iterations, seed
def _module_weights(config: io.ManuscriptConfig) -> dict[str, float]:
    scoring_data = _mapping(config.data.get("scoring"), "scoring")
    weights = _mapping(scoring_data.get("module_weights"), "scoring.module_weights")
    return {str(key): float(value) for key, value in weights.items()}
def _case_ids(config: io.ManuscriptConfig) -> tuple[str, ...]:
    return posthoc.representative_case_ids(config)
def _write_in_order(frame: pd.DataFrame, path: Path, columns: Sequence[str]) -> None:
    if frame.columns.tolist() != list(columns):
        raise ValueError(f"{path.name} has incorrect schema: {frame.columns.tolist()}")
    output = frame.copy().map(io._lowercase_boolean)
    with path.open("w", encoding=io.ENCODING, newline="") as handle:
        output.to_csv(
            handle, index=False, lineterminator=io.LINE_TERMINATOR,
            float_format=io.FLOAT_FORMAT, na_rep=io.NA_REP,
        )
def _write_analysis_bundle(frames: Mapping[str, pd.DataFrame], output: Path) -> None:
    if set(frames) != set(ANALYSIS_SCHEMAS):
        raise ValueError("analysis bundle is incomplete")
    _prepare_managed_directory(output, set(ANALYSIS_SCHEMAS), "analysis")
    for filename, columns in ANALYSIS_SCHEMAS.items():
        frame = frames[filename]
        if frame.columns.tolist() != list(columns):
            raise ValueError(f"{filename} has incorrect schema: {frame.columns.tolist()}")
        if filename in {"final_bootstrap_intervals.csv", "representative_case_results.csv"}:
            _write_in_order(frame, output / filename, columns)
        else:
            io.write_stable_csv(
                frame, output / filename, columns=columns,
                sort_by=ANALYSIS_SORTS[filename],
            )
def _run_analyze(arguments: argparse.Namespace) -> int:
    config = io.load_manuscript_config(arguments.config)
    output = _safe_output(config, arguments.output)
    _prepare_managed_directory(output, set(ANALYSIS_SCHEMAS), "analysis")
    static_dir = arguments.static_dir.resolve()
    rankings = _load_static_rankings(static_dir)
    feature_table = io.load_feature_table(config.resolve_configured_path("feature_table"))
    record_mapping_path = arguments.cluster_dir.resolve() / "tool_record_to_cluster_v2_mapping.csv"
    if not record_mapping_path.is_file():
        raise FileNotFoundError(
            "single-tool analysis requires the sibling clustering/"
            "tool_record_to_cluster_v2_mapping.csv output"
        )
    record_mapping = pd.read_csv(record_mapping_path, keep_default_na=True)
    references, md_runs, redocking = _load_posthoc_inputs(arguments.posthoc_data.resolve())
    systems = pd.read_csv(
        config.resolve_configured_path("systems"), sep="\t", keep_default_na=False,
    )
    mapped = posthoc.run_posthoc(
        posthoc.PosthocInputs(
            rankings=rankings,
            structures=io.load_structure_paths(
                config, target_ids=set(rankings["pdb_id"].astype(str)),
            ),
            references=references,
            md_runs=md_runs,
            redocking=redocking,
            systems=systems,
        ),
        config,
    )
    target = statistics.build_target_summary(
        rankings, mapped.evidence_labels, mapped.reference_evaluable, systems,
    )
    overall = pd.DataFrame([
        {"analysis_level": "overall", "group": "all_targets", **statistics.summarize_metrics(target)}
    ], columns=statistics.POINT_COLUMNS)
    categories = statistics.summarize_by_group(target, "Pocket_category")
    iterations, seed = _statistics_config(config)
    bootstrap = statistics.bootstrap_intervals(target, iterations=iterations, seed=seed)
    family = statistics.leave_one_family_out(target)
    ranking_data = _mapping(config.data.get("ranking"), "ranking")
    ablation = statistics.orel_ablation(
        rankings, mapped.evidence_labels, mapped.reference_evaluable, systems,
        _module_weights(config), tie_method=str(ranking_data["tie_method"]),
    )
    weight = statistics.weight_sensitivity(
        rankings, mapped.evidence_labels, _module_weights(config),
        perturbation_fraction=0.20, tie_method=str(ranking_data["tie_method"]),
    )
    single_tool = statistics.single_tool_prioritization(
        feature_table, record_mapping, rankings, mapped.evidence_labels,
    )
    cases = statistics.build_representative_cases(
        rankings, mapped.evidence_labels, mapped.reference_mapping,
        mapped.md_mapping, mapped.redocking_mapping, ablation.targets, _case_ids(config),
    )
    frames = {
        "final_top3_automated_QC.csv": mapped.topk_qc,
        "final_reference_mapping.csv": mapped.reference_mapping,
        "final_md_cluster_v2_mapping.csv": mapped.md_mapping,
        "final_automated_evidence_labels.csv": mapped.evidence_labels,
        "final_redocking_cluster_v2_mapping.csv": mapped.redocking_mapping,
        "unresolved_cases.csv": mapped.unresolved_cases,
        "final_cluster_v2_master_table.csv": mapped.convenience_master,
        "target_level_candidate_prioritization.csv": target,
        "final_candidate_prioritization_metrics.csv": overall,
        "final_category_metrics.csv": categories,
        "final_bootstrap_intervals.csv": bootstrap,
        "final_family_sensitivity.csv": family,
        "final_orel_ablation_targets.csv": ablation.targets,
        "final_orel_ablation_categories.csv": ablation.categories,
        "orel_ablation_cluster_rankings.csv": ablation.cluster_rankings,
        "representative_case_results.csv": cases,
        "weight_sensitivity_scenarios.csv": weight.scenarios,
        "weight_sensitivity_targets.csv": weight.targets,
        "single_tool_target_ranks.csv": single_tool.target_ranks,
        "single_tool_complete_case_metrics.csv": single_tool.complete_case_metrics,
    }
    _write_analysis_bundle(frames, output)
    print(json.dumps({
        "targets": len(target), "clusters": len(rankings), "bootstrap_rows": len(bootstrap),
        "analysis_files": len(frames), "weight_scenarios": len(weight.scenarios),
        "single_tool_complete_cases": int(single_tool.target_ranks.loc[
            single_tool.target_ranks["method"].eq("OIPS-P"), "complete_case"
        ].sum()),
    }, indent=2))
    return 0
def _validate_analysis_relationships(
    frames: Mapping[str, pd.DataFrame], *, top_k: int = 3,
    case_ids: Sequence[str] | None = None, iterations: int | None = None,
    seed: int | None = None,
) -> None:
    target = frames["target_level_candidate_prioritization.csv"]; targets = set(target["pdb_id"])
    categories, families = set(target["Pocket_category"]), set(target["Protein_family"])
    checks = (
        ("final_md_cluster_v2_mapping.csv", "pdb_id", targets, "MD mapping target set"), ("final_redocking_cluster_v2_mapping.csv", "pdb_id", targets, "redocking target set"),
        ("final_orel_ablation_targets.csv", "pdb_id", targets, "ablation target set"), ("final_category_metrics.csv", "group", categories, "category metric set"),
        ("final_orel_ablation_categories.csv", "Pocket_category", categories, "ablation category set"), ("final_family_sensitivity.csv", "excluded_family", families, "family sensitivity set"),
    )
    for filename, column, expected, label in checks:
        if set(frames[filename][column]) != expected:
            raise ValueError(f"{label} does not match the target summary")
    master = frames["final_cluster_v2_master_table.csv"]; columns = ["pdb_id", "cluster_v2_id", "Within_PDB_rank"]
    expected_qc = master.loc[master["Within_PDB_rank"].le(top_k), columns]
    qc = frames["final_top3_automated_QC.csv"]
    expected_rows = set(expected_qc.itertuples(index=False, name=None))
    observed_rows = set(qc[columns].itertuples(index=False, name=None))
    if observed_rows != expected_rows or len(qc) != len(expected_qc):
        raise ValueError("QC keys and ranks do not exactly match configured master Top-k")
    master_keys = set(zip(master["pdb_id"], master["cluster_v2_id"]))
    for filename in (
        "final_reference_mapping.csv", "final_automated_evidence_labels.csv",
        "orel_ablation_cluster_rankings.csv",
    ):
        keys = set(zip(frames[filename]["pdb_id"], frames[filename]["cluster_v2_id"]))
        if keys != master_keys or len(frames[filename]) != len(master):
            raise ValueError(f"{filename} cluster keys do not match the final master")
    overall = frames["final_candidate_prioritization_metrics.csv"]
    if len(overall) != 1 or overall.loc[0, "analysis_level"] != "overall":
        raise ValueError("overall analysis metrics must contain exactly one overall row")
    expected_bootstrap = [
        (method, metric) for method in ("target_resampling", "family_clustered")
        for metric in statistics.BOOTSTRAP_METRICS
    ]
    bootstrap = frames["final_bootstrap_intervals.csv"]
    if list(zip(bootstrap["bootstrap_method"], bootstrap["metric"])) != expected_bootstrap:
        raise ValueError("bootstrap rows must be exactly two methods by seven metrics")
    if iterations is not None and not bootstrap["iterations"].eq(iterations).all():
        raise ValueError("bootstrap iteration metadata does not match configuration")
    if seed is not None and not bootstrap["random_seed"].eq(seed).all():
        raise ValueError("bootstrap seed metadata does not match configuration")
    if case_ids is not None and tuple(frames["representative_case_results.csv"]["pdb_id"]) != tuple(case_ids):
        raise ValueError("representative cases are not in configured order")
def _load_analysis_bundle(
    path: Path, config: io.ManuscriptConfig | None = None,
) -> dict[str, pd.DataFrame]:
    actual = _directory_files(path, "analysis")
    if actual != set(ANALYSIS_SCHEMAS):
        raise ValueError(f"analysis directory file set is invalid: {sorted(actual)}")
    frames: dict[str, pd.DataFrame] = {}
    for filename, columns in ANALYSIS_SCHEMAS.items():
        raw = pd.read_csv(path / filename, dtype=str, keep_default_na=False)
        frames[filename] = _validate_analysis_csv(filename, raw)
    targets = set(frames["target_level_candidate_prioritization.csv"]["pdb_id"])
    if not targets or any(not __import__("re").fullmatch(r"[0-9][A-Z0-9]{3}", str(value)) for value in targets):
        raise ValueError("analysis bundle must contain uppercase PDB target IDs")
    for filename, frame in frames.items():
        if "pdb_id" in frame and not set(frame["pdb_id"]).issubset(targets):
            raise ValueError(f"{filename} contains a target outside the target summary")
    if config is None:
        _validate_analysis_relationships(frames)
    else:
        posthoc_data = _mapping(config.data["posthoc"], "posthoc")
        qc_data = _mapping(posthoc_data["top3_qc"], "posthoc.top3_qc")
        iterations, seed = _statistics_config(config)
        _validate_analysis_relationships(
            frames, top_k=int(qc_data["top_k"]), case_ids=_case_ids(config),
            iterations=iterations, seed=seed,
        )
    return frames
@contextmanager
def _figure_runtime():
    import importlib, os, tempfile
    cache = tempfile.TemporaryDirectory(prefix="oips-matplotlib-")
    previous = os.environ.get("MPLCONFIGDIR")
    os.environ["MPLCONFIGDIR"] = cache.name
    try:
        yield importlib.import_module(".figures", __package__)
    finally:
        os.environ.pop("MPLCONFIGDIR", None) if previous is None else os.environ.__setitem__("MPLCONFIGDIR", previous)
        cache.cleanup()
def _run_figures(arguments: argparse.Namespace) -> int:
    with _figure_runtime() as figures:
        return _render_figures(arguments, figures)
def _render_figures(arguments: argparse.Namespace, figures: object) -> int:
    config = io.load_manuscript_config(arguments.config)
    output = _safe_output(config, arguments.output)
    contract = figures.load_figure_contract(config.resolve_configured_path("figure_contract"))
    expected = {
        filename
        for spec in contract.figures.values()
        for filename in (
            f"{spec.stem}.svg", f"{spec.stem}.png", f"{spec.stem}_source_data.csv",
        )
    }
    _prepare_managed_directory(output, expected, "figure")
    frames = _load_analysis_bundle(arguments.analysis.resolve(), config)
    configured_cases = _case_ids(config)
    observed_cases = tuple(frames["representative_case_results.csv"]["pdb_id"])
    if observed_cases != configured_cases:
        raise ValueError("representative cases are not in the configured order")
    tables = _figure_tables(frames)
    written = figures.render_all_figures(tables, contract, configured_cases, output)
    if {path.name for path in written} != expected:
        raise ValueError("figure renderer did not produce the managed file set")
    print(json.dumps({"figures": 4, "files": len(written), "png_dpi": contract.png_dpi}, indent=2))
    return 0
def _figure_tables(frames: Mapping[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    return {
        "target_summary": frames["target_level_candidate_prioritization.csv"],
        "bootstrap": frames["final_bootstrap_intervals.csv"],
        "topk_qc": frames["final_top3_automated_QC.csv"],
        "orel_targets": frames["final_orel_ablation_targets.csv"],
        "evidence_labels": frames["final_automated_evidence_labels.csv"],
        "md_mapping": frames["final_md_cluster_v2_mapping.csv"],
        "redocking_mapping": frames["final_redocking_cluster_v2_mapping.csv"],
        "category_metrics": frames["final_category_metrics.csv"],
        "master": frames["final_cluster_v2_master_table.csv"],
        "representative_cases": frames["representative_case_results.csv"],
    }
REPORT_NAMES = release.REPORT_NAMES

def _bundle_stage_paths(bundle: Path, *, reference_layout: bool = False) -> list[Path]:
    paths: list[Path] = []
    expected = {
        "clustering": set(io.CLUSTER_FILE_SCHEMAS), "static": set(io.SCORE_FILE_SCHEMAS),
        "analysis": set(ANALYSIS_SCHEMAS),
    }
    for directory, names in expected.items():
        if _directory_files(bundle / directory, directory) != names:
            raise ValueError(f"{directory} bundle file set is invalid")
        paths.extend(bundle / directory / name for name in sorted(names))
    paths.extend(release.figure_stage_paths(bundle, reference_layout=reference_layout))
    return paths


def _validate_bundle_stages(
    config: io.ManuscriptConfig, bundle: Path, *, reference_layout: bool = False,
) -> dict[str, pd.DataFrame]:
    features = io.load_feature_table(config.resolve_configured_path("feature_table"))
    clustered = io.load_clustering_result(bundle / "clustering", features=features)
    master, rankings = _load_static_tables(bundle / "static")
    if set(rankings["cluster_v2_id"]) != set(clustered.candidates["cluster_v2_id"]):
        raise ValueError("static rankings do not match clustering candidates")
    validation.validate_static_science(master, rankings, config)
    frames = _load_analysis_bundle(bundle / "analysis", config)
    validation.validate_bundle_relationships(clustered.candidates, rankings, frames, config)
    with _figure_runtime() as figure_module:
        contract = figure_module.load_figure_contract(config.resolve_configured_path("figure_contract"))
        release.validate_figure_artifacts(
            bundle / ("figure_source_data" if reference_layout else "figures"),
            _figure_tables(frames), contract, _case_ids(config), figure_module,
            source_only=reference_layout,
        )
    return frames


def _report_manifest_view(manifest: Mapping[str, object]) -> dict[str, object]:
    view = dict(manifest)
    view["outputs"] = [entry for entry in manifest["outputs"] if not str(entry["path"]).startswith("reports/")]
    view["managed_paths"] = [path for path in manifest["managed_paths"] if not str(path).startswith("reports/")]
    return view


def _report_texts(
    bundle: Path, manifest_view: Mapping[str, object], frames: Mapping[str, pd.DataFrame],
    check_report: validation.ValidationReport, snapshot: Mapping[str, object],
) -> dict[str, str]:
    actual = provenance.snapshot_actual(bundle, snapshot)
    initial = {
        "rebuild_report.md": reports.build_rebuild_report(manifest_view),
        "analysis_report.md": reports.build_analysis_report(frames, manifest_view),
        "numeric_crosscheck.md": reports.build_numeric_crosscheck(
            actual, validation.numeric_snapshot_sections(snapshot),
        ),
    }
    scanned = provenance.scan_sensitive_content(initial)
    sensitive = validation.CheckResult(
        "report_sensitive_content", scanned.status, scanned.summary, scanned.evidence,
    )
    combined = validation.ValidationReport((*check_report.checks, sensitive))
    combined.raise_for_failures()
    final = {**initial, "validation_report.md": reports.build_validation_report(combined)}
    final_scan = provenance.scan_sensitive_content(final)
    if final_scan.status == "fail":
        raise ValueError("sensitive content detected in final reports")
    return final


def _write_reports(bundle: Path, texts: Mapping[str, str]) -> list[Path]:
    directory = bundle / "reports"
    _prepare_managed_directory(directory, REPORT_NAMES, "reports")
    paths: list[Path] = []
    for name in sorted(REPORT_NAMES):
        path = directory / name
        path.write_text(texts[name], encoding="utf-8", newline="\n")
        paths.append(path)
    return paths


def _pre_manifest_checks(
    config: io.ManuscriptConfig, bundle: Path, snapshot_path: Path,
) -> tuple[dict[str, pd.DataFrame], validation.ValidationReport, dict[str, object]]:
    inputs = validation.validate_public_inputs(config)
    frames = _validate_bundle_stages(config, bundle)
    snapshot_report = validation.validate_snapshot(bundle, snapshot_path)
    report = validation.ValidationReport((
        *inputs.checks,
        validation.CheckResult("bundle_contract", "pass", "stage contracts validated", {}),
        *snapshot_report.checks,
    ))
    report.raise_for_failures()
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    return frames, report, snapshot


def _run_reproduce(arguments: argparse.Namespace) -> int:
    config = io.load_manuscript_config(arguments.config)
    validation.validate_public_inputs(config).raise_for_failures()
    output = arguments.output or config.resolve_configured_path("default_output")
    bundle = provenance.prepare_output_bundle(
        output, config.resolve_configured_path("reference_results"), force=arguments.force,
    )
    _run_cluster(argparse.Namespace(config=arguments.config, output=bundle / "clustering"))
    _run_score(argparse.Namespace(config=arguments.config, cluster_dir=bundle / "clustering", output=bundle / "static"))
    _run_analyze(argparse.Namespace(
        config=arguments.config, cluster_dir=bundle / "clustering", static_dir=bundle / "static",
        posthoc_data=config.resolve_configured_path("reference_annotations").parent,
        output=bundle / "analysis",
    ))
    _run_figures(argparse.Namespace(config=arguments.config, analysis=bundle / "analysis", output=bundle / "figures"))
    snapshot_path = config.resolve_configured_path("expected_summary")
    frames, check_report, snapshot = _pre_manifest_checks(config, bundle, snapshot_path)
    stage_paths = _bundle_stage_paths(bundle)
    summary = provenance.cluster_summary(bundle)
    draft = provenance.build_run_manifest(
        config, bundle, consumed_inputs=validation.consumed_input_paths(config),
        output_paths=stage_paths,
        managed_paths=[*(path.relative_to(bundle).as_posix() for path in stage_paths), "run_manifest.json"],
        summary=summary,
    )
    report_paths = _write_reports(bundle, _report_texts(bundle, draft, frames, check_report, snapshot))
    outputs = [*stage_paths, *report_paths]
    final = provenance.build_run_manifest(
        config, bundle, consumed_inputs=validation.consumed_input_paths(config),
        output_paths=outputs,
        managed_paths=[*(path.relative_to(bundle).as_posix() for path in outputs), "run_manifest.json"],
        summary=summary,
    )
    provenance.write_run_manifest(bundle, final)
    checked = provenance.validate_run_manifest(bundle)
    provenance.validate_manifest_context(config, checked, validation.consumed_input_paths(config))
    _compact({"status": "pass", "files": len(outputs) + 1, "failures": []})
    return 0


def _run_verify(arguments: argparse.Namespace) -> int:
    config = io.load_manuscript_config(arguments.config)
    bundle = arguments.bundle.absolute()
    reference = config.resolve_configured_path("reference_results")
    immutable = io.path_is_equal_or_within(bundle, reference)
    reference_layout = immutable and io.path_is_equal_or_within(reference, bundle)
    manifest, frames, check_report, content_report, snapshot, present = provenance.collect_verification(
        config, bundle, arguments.snapshot,
        lambda: _validate_bundle_stages(
            config, bundle, reference_layout=reference_layout,
        ),
        immutable=immutable and not reference_layout,
    )
    failures = [check.check_id for check in check_report.checks if check.status == "fail"]
    if failures:
        _compact({"status": "fail", "checks": len(check_report.checks), "failures": failures})
        return 1
    if manifest is None:
        stage_paths = _bundle_stage_paths(bundle, reference_layout=reference_layout)
        managed = [path.relative_to(bundle).as_posix() for path in stage_paths]
        manifest = provenance.build_run_manifest(
            config, bundle, consumed_inputs=validation.consumed_input_paths(config),
            output_paths=stage_paths,
            managed_paths=[*managed, "run_manifest.json"],
            summary=provenance.cluster_summary(bundle),
        )
    texts = _report_texts(bundle, _report_manifest_view(manifest), frames, content_report, snapshot)
    if reference_layout:
        release.compare_reference_reports(bundle / "reports", texts)
    elif immutable:
        for name, text in texts.items():
            if (bundle / "reports" / name).read_text(encoding="utf-8") != text:
                raise ValueError(f"immutable reference report mismatch: {name}")
    else:
        _write_reports(bundle, texts)
        if present:
            provenance.validate_run_manifest(bundle)
    _compact({"status": "pass", "checks": len(check_report.checks), "failures": []})
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Parse arguments, execute a stage, and return a process exit code."""
    try:
        arguments = _build_parser().parse_args(argv)
    except SystemExit as error:
        return int(error.code)
    try:
        if arguments.command == "validate-data":
            return _run_validate_data(arguments)
        if arguments.command == "cluster":
            return _run_cluster(arguments)
        if arguments.command == "score":
            return _run_score(arguments)
        if arguments.command == "analyze":
            return _run_analyze(arguments)
        if arguments.command == "figures":
            return _run_figures(arguments)
        if arguments.command == "verify":
            return _run_verify(arguments)
        if arguments.command == "reproduce":
            return _run_reproduce(arguments)
        return 0
    except (AssertionError, KeyError, OSError, TypeError, ValueError,
            pd.errors.ParserError, yaml.YAMLError) as error:
        print(f"error: {_safe_error_message(error)}", file=sys.stderr)
        return 1
