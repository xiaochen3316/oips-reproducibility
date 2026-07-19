"""Pure statistics, uncertainty, and O-rel ablation for OIPS."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math

import numpy as np
import pandas as pd

from .single_tool import (
    SINGLE_TOOL_METRIC_COLUMNS,
    SINGLE_TOOL_TARGET_COLUMNS,
    SingleToolResult,
    single_tool_prioritization,
)


METRIC_FIELDS = (
    "reference_evaluable_N", "Reference_Top1", "Reference_Top3",
    "Reference_Top5", "Reference_MRR", "Reference_median_rank",
    "Reference_median_rank_percentile", "first_supported_evaluable_N",
    "First_supported_Top1", "First_supported_Top3", "First_supported_MRR",
)
TARGET_COLUMNS = (
    "pdb_id", "Pocket_category", "Protein_family", "formal_cluster_v2_count",
    "reference_evaluable", "R_auto_cluster_count", "reference_first_rank",
    "reference_rank_percentile", "reference_top1", "reference_top3",
    "reference_top5", "reference_reciprocal_rank",
    "first_supported_evaluable", "first_supported_rank",
    "first_supported_top1", "first_supported_top3",
    "first_supported_reciprocal_rank", "static_top_cluster_v2_id",
    "static_top_evidence_label",
)
POINT_COLUMNS = ("analysis_level", "group", *METRIC_FIELDS)
BOOTSTRAP_COLUMNS = (
    "bootstrap_method", "metric", "point_estimate", "CI_2.5_percent",
    "CI_97.5_percent", "iterations", "random_seed",
)
FAMILY_COLUMNS = (
    "excluded_family", "excluded_target_N", "remaining_target_N", *METRIC_FIELDS,
)
ABLATION_TARGET_COLUMNS = (
    "analysis_level", "pdb_id", "Pocket_category", "Protein_family",
    "Full_reference_rank", "Without_O_rel_reference_rank",
    "reference_rank_change_without_minus_full", "Full_first_supported_rank",
    "Without_O_rel_first_supported_rank",
    "first_supported_rank_change_without_minus_full", "Full_reference_Top3",
    "Without_O_rel_reference_Top3", "reference_Top3_inclusion_change",
    "rank_Spearman_rho", "Full_top_cluster_v2_id",
    "Without_O_rel_top_cluster_v2_id", "top_cluster_identity_changed",
)
ABLATION_CATEGORY_COLUMNS = (
    "analysis_level", "Pocket_category", "target_N", "Full_reference_Top1",
    "Without_O_rel_reference_Top1", "Full_reference_Top3",
    "Without_O_rel_reference_Top3", "Full_reference_MRR",
    "Without_O_rel_reference_MRR", "mean_rank_Spearman_rho",
    "top_cluster_identity_change_fraction",
)
ABLATION_CLUSTER_COLUMNS = (
    "pdb_id", "cluster_v2_id", "Within_PDB_rank", "Without_O_rel_score",
    "Without_O_rel_rank",
)
REPRESENTATIVE_COLUMNS = (
    "pdb_id", "static_top_cluster_v2_id", "static_top_score",
    "static_top_label", "static_top_tool_support", "static_top_chains",
    "static_top_interface_overlap", "reference_cluster_v2_id",
    "reference_rank", "reference_DCC_A", "without_O_rel_reference_rank",
    "without_O_rel_top_cluster_v2_id", "MD_calls", "redocking_RMSD_A",
    "redocking_status",
)
WEIGHT_SCENARIO_COLUMNS = (
    "scenario_id", "perturbed_module", "direction", "multiplier", "target_N",
    "baseline_top1_retained_N", "mean_top3_jaccard", "median_spearman_rho",
    "Reference_Top1_N", "Reference_Top3_N", "First_supported_Top1_N",
    "First_supported_Top3_N",
)
WEIGHT_TARGET_COLUMNS = (
    "pdb_id", "baseline_top_cluster_v2_id", "top1_retention_count",
    "top1_changed_scenarios", "baseline_reference_rank", "minimum_reference_rank",
    "maximum_reference_rank", "baseline_first_supported_rank",
    "minimum_first_supported_rank", "maximum_first_supported_rank",
)
BOOTSTRAP_METRICS = (
    "Reference_Top1", "Reference_Top3", "Reference_Top5", "Reference_MRR",
    "First_supported_Top1", "First_supported_Top3", "First_supported_MRR",
)
STATIC_MODULES = ("C_cons", "G_geo", "P_lig", "O_rel_formal", "Q_evidence")


@dataclass(frozen=True)
class AblationResult:
    targets: pd.DataFrame
    categories: pd.DataFrame
    cluster_rankings: pd.DataFrame


@dataclass(frozen=True)
class WeightSensitivityResult:
    scenarios: pd.DataFrame
    targets: pd.DataFrame


def _require(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def _strict_bool(value: object, label: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str) and value.strip().casefold() in {"true", "false"}:
        return value.strip().casefold() == "true"
    raise ValueError(f"{label} must be true or false")


def _boolean_series(series: pd.Series, label: str) -> pd.Series:
    return pd.Series([_strict_bool(value, label) for value in series], index=series.index, dtype=bool)


def _metadata(systems: pd.DataFrame) -> pd.DataFrame:
    category = "pocket_category" if "pocket_category" in systems else "Pocket_category"
    family = "protein_family" if "protein_family" in systems else "Protein_family"
    _require(systems, ("pdb_id", category, family), "systems")
    result = systems.loc[:, ["pdb_id", category, family]].copy()
    result.columns = ["pdb_id", "Pocket_category", "Protein_family"]
    result["pdb_id"] = result["pdb_id"].astype(str).str.upper()
    if result["pdb_id"].duplicated().any():
        raise ValueError("systems contains duplicate target IDs")
    return result.set_index("pdb_id")


def _reference_flags(
    reference_evaluable: pd.Series | Mapping[str, object], targets: set[str]
) -> dict[str, bool]:
    series = pd.Series(reference_evaluable, dtype=object)
    series.index = series.index.astype(str).str.upper()
    if series.index.duplicated().any() or set(series.index) != targets:
        raise ValueError("reference_evaluable must cover each ranked target exactly once")
    return {str(key): _strict_bool(value, "reference_evaluable") for key, value in series.items()}


def build_target_summary(
    rankings: pd.DataFrame,
    labels: pd.DataFrame,
    reference_evaluable: pd.Series | Mapping[str, object],
    systems: pd.DataFrame,
) -> pd.DataFrame:
    """Build one row per target using an explicit reference denominator."""
    _require(rankings, ("pdb_id", "cluster_v2_id", "Within_PDB_rank"), "rankings")
    _require(labels, ("pdb_id", "cluster_v2_id", "automated_evidence_label"), "labels")
    ranked = rankings.copy(deep=True)
    marked = labels.copy(deep=True)
    for frame in (ranked, marked):
        frame["pdb_id"] = frame["pdb_id"].astype(str).str.upper()
    keys = ["pdb_id", "cluster_v2_id"]
    if ranked.duplicated(keys).any() or marked.duplicated(keys).any():
        raise ValueError("rankings and labels require unique target/cluster keys")
    ranked_keys = pd.MultiIndex.from_frame(ranked[keys])
    marked_lookup = marked.set_index(keys)["automated_evidence_label"]
    if set(ranked_keys) != set(marked_lookup.index):
        raise ValueError("label coverage does not match static rankings")
    ranked["label"] = marked_lookup.reindex(ranked_keys).to_numpy()
    targets = set(ranked["pdb_id"])
    flags = _reference_flags(reference_evaluable, targets)
    metadata = _metadata(systems)
    if set(metadata.index) != targets:
        raise ValueError("systems target coverage does not match static rankings")

    rows: list[dict[str, object]] = []
    for pdb_id, group in ranked.groupby("pdb_id", sort=True):
        group = group.copy()
        group["Within_PDB_rank"] = pd.to_numeric(group["Within_PDB_rank"], errors="raise")
        ordered = group.sort_values(["Within_PDB_rank", "cluster_v2_id"], kind="mergesort")
        references = ordered.loc[ordered["label"].eq("R_auto")]
        supported = ordered.loc[ordered["label"].isin(["R_auto", "A_auto"])]
        reference_rank = float(references["Within_PDB_rank"].min()) if len(references) else math.nan
        supported_rank = float(supported["Within_PDB_rank"].min()) if len(supported) else math.nan
        evaluable = flags[pdb_id]
        count = len(group)
        percentile = (
            math.nan if math.isnan(reference_rank)
            else 1.0 if count <= 1
            else 1.0 - (reference_rank - 1.0) / (count - 1.0)
        )
        reference_reciprocal = (
            math.nan if not evaluable else 0.0 if math.isnan(reference_rank) else 1.0 / reference_rank
        )
        first_evaluable = not math.isnan(supported_rank)
        rows.append({
            "pdb_id": pdb_id,
            "Pocket_category": metadata.loc[pdb_id, "Pocket_category"],
            "Protein_family": metadata.loc[pdb_id, "Protein_family"],
            "formal_cluster_v2_count": count,
            "reference_evaluable": evaluable,
            "R_auto_cluster_count": len(references),
            "reference_first_rank": reference_rank,
            "reference_rank_percentile": percentile,
            "reference_top1": (reference_rank <= 1) if evaluable else math.nan,
            "reference_top3": (reference_rank <= 3) if evaluable else math.nan,
            "reference_top5": (reference_rank <= 5) if evaluable else math.nan,
            "reference_reciprocal_rank": reference_reciprocal,
            "first_supported_evaluable": first_evaluable,
            "first_supported_rank": supported_rank,
            "first_supported_top1": supported_rank <= 1 if first_evaluable else False,
            "first_supported_top3": supported_rank <= 3 if first_evaluable else False,
            "first_supported_reciprocal_rank": 1.0 / supported_rank if first_evaluable else math.nan,
            "static_top_cluster_v2_id": ordered.iloc[0]["cluster_v2_id"],
            "static_top_evidence_label": ordered.iloc[0]["label"],
        })
    return pd.DataFrame(rows, columns=TARGET_COLUMNS).sort_values("pdb_id", kind="mergesort").reset_index(drop=True)


def _endpoint_mean(frame: pd.DataFrame, column: str) -> float:
    if not len(frame):
        return math.nan
    return float(_boolean_series(frame[column], column).mean())


def summarize_metrics(frame: pd.DataFrame) -> dict[str, float]:
    """Summarize the two non-equivalent target-level endpoints."""
    required = (
        "reference_evaluable", "reference_top1", "reference_top3", "reference_top5",
        "reference_reciprocal_rank", "reference_first_rank", "reference_rank_percentile",
        "first_supported_evaluable", "first_supported_top1", "first_supported_top3",
        "first_supported_reciprocal_rank",
    )
    _require(frame, required, "target summary")
    reference_mask = _boolean_series(frame["reference_evaluable"], "reference_evaluable")
    supported_mask = _boolean_series(frame["first_supported_evaluable"], "first_supported_evaluable")
    reference = frame.loc[reference_mask]
    supported = frame.loc[supported_mask]
    reciprocal = pd.to_numeric(reference["reference_reciprocal_rank"], errors="raise")
    supported_reciprocal = pd.to_numeric(supported["first_supported_reciprocal_rank"], errors="raise")
    ranks = pd.to_numeric(reference["reference_first_rank"], errors="coerce").dropna()
    percentiles = pd.to_numeric(reference["reference_rank_percentile"], errors="coerce").dropna()
    return {
        "reference_evaluable_N": float(len(reference)),
        "Reference_Top1": _endpoint_mean(reference, "reference_top1"),
        "Reference_Top3": _endpoint_mean(reference, "reference_top3"),
        "Reference_Top5": _endpoint_mean(reference, "reference_top5"),
        "Reference_MRR": float(reciprocal.mean()) if len(reference) else math.nan,
        "Reference_median_rank": float(ranks.median()) if len(ranks) else math.nan,
        "Reference_median_rank_percentile": float(percentiles.median()) if len(percentiles) else math.nan,
        "first_supported_evaluable_N": float(len(supported)),
        "First_supported_Top1": _endpoint_mean(supported, "first_supported_top1"),
        "First_supported_Top3": _endpoint_mean(supported, "first_supported_top3"),
        "First_supported_MRR": float(supported_reciprocal.mean()) if len(supported) else math.nan,
    }


def summarize_by_group(frame: pd.DataFrame, group_column: str = "Pocket_category") -> pd.DataFrame:
    _require(frame, (group_column,), "target summary")
    rows = [
        {"analysis_level": "pocket_archetype", "group": name, **summarize_metrics(group)}
        for name, group in frame.groupby(group_column, sort=True)
    ]
    return pd.DataFrame(rows, columns=POINT_COLUMNS).sort_values("group", kind="mergesort").reset_index(drop=True)


def bootstrap_intervals(
    target_summary: pd.DataFrame, *, iterations: int, seed: int
) -> pd.DataFrame:
    if not isinstance(iterations, int) or iterations < 1:
        raise ValueError("bootstrap iterations must be a positive integer")
    if not isinstance(seed, int):
        raise ValueError("bootstrap seed must be an integer")
    _require(target_summary, ("pdb_id", "Protein_family"), "target summary")
    target = target_summary.sort_values("pdb_id", kind="mergesort").reset_index(drop=True)
    point = summarize_metrics(target)
    rows: list[dict[str, object]] = []
    for method in ("target_resampling", "family_clustered"):
        rng = np.random.default_rng(seed)
        values = {metric: [] for metric in BOOTSTRAP_METRICS}
        if method == "target_resampling":
            positions = np.arange(len(target))
            for _ in range(iterations):
                sample = target.iloc[rng.choice(positions, size=len(positions), replace=True)]
                metrics = summarize_metrics(sample)
                for metric in BOOTSTRAP_METRICS:
                    values[metric].append(metrics[metric])
        else:
            families = sorted(target["Protein_family"].astype(str).unique())
            groups = {name: target.loc[target["Protein_family"].astype(str).eq(name)] for name in families}
            for _ in range(iterations):
                chosen = rng.choice(families, size=len(families), replace=True)
                sample = pd.concat([groups[name] for name in chosen], ignore_index=True)
                metrics = summarize_metrics(sample)
                for metric in BOOTSTRAP_METRICS:
                    values[metric].append(metrics[metric])
        for metric in BOOTSTRAP_METRICS:
            vector = np.asarray(values[metric], dtype=float)
            rows.append({
                "bootstrap_method": method,
                "metric": metric,
                "point_estimate": point[metric],
                "CI_2.5_percent": float(np.nanquantile(vector, 0.025)),
                "CI_97.5_percent": float(np.nanquantile(vector, 0.975)),
                "iterations": iterations,
                "random_seed": seed,
            })
    return pd.DataFrame(rows, columns=BOOTSTRAP_COLUMNS)


def leave_one_family_out(target_summary: pd.DataFrame) -> pd.DataFrame:
    _require(target_summary, ("Protein_family",), "target summary")
    rows: list[dict[str, object]] = []
    families = target_summary["Protein_family"].astype(str)
    for family in sorted(families.unique()):
        keep = ~families.eq(family)
        rows.append({
            "excluded_family": family,
            "excluded_target_N": int((~keep).sum()),
            "remaining_target_N": int(keep.sum()),
            **summarize_metrics(target_summary.loc[keep]),
        })
    return pd.DataFrame(rows, columns=FAMILY_COLUMNS)


def _weighted_available(frame: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    if set(weights) != set(STATIC_MODULES):
        raise ValueError(f"module weights must contain exactly {list(STATIC_MODULES)}")
    numeric_weights = {key: float(value) for key, value in weights.items()}
    if any(not math.isfinite(value) or value <= 0 for value in numeric_weights.values()):
        raise ValueError("module weights must be finite and positive")
    numerator = np.zeros(len(frame), dtype=float)
    denominator = np.zeros(len(frame), dtype=float)
    for module in STATIC_MODULES:
        values = pd.to_numeric(frame[module], errors="coerce").to_numpy(float)
        valid = np.isfinite(values)
        numerator[valid] += numeric_weights[module] * values[valid]
        denominator[valid] += numeric_weights[module]
    result = np.full(len(frame), np.nan)
    valid = denominator > 0
    result[valid] = numerator[valid] / denominator[valid]
    return pd.Series(result, index=frame.index)


def weight_sensitivity(
    rankings: pd.DataFrame,
    labels: pd.DataFrame,
    module_weights: Mapping[str, float],
    *,
    perturbation_fraction: float = 0.20,
    tie_method: str = "average",
) -> WeightSensitivityResult:
    """Run ten one-at-a-time ±20% static-module perturbations.

    Candidate definitions, module values, missingness, and evidence labels are
    fixed.  Changing one positive weight and applying the missing-aware weighted
    mean is algebraically equivalent to explicitly normalizing all five weights
    before scoring.
    """
    _require(rankings, ("pdb_id", "cluster_v2_id", "Within_PDB_rank", *STATIC_MODULES), "rankings")
    _require(labels, ("pdb_id", "cluster_v2_id", "automated_evidence_label"), "labels")
    if tie_method != "average":
        raise ValueError("weight sensitivity requires average tie ranking")
    if not math.isfinite(float(perturbation_fraction)) or not 0 < float(perturbation_fraction) < 1:
        raise ValueError("perturbation_fraction must lie strictly between zero and one")
    work = rankings.copy(deep=True)
    marked = labels.copy(deep=True)
    for frame in (work, marked):
        frame["pdb_id"] = frame["pdb_id"].astype(str).str.upper()
    keys = ["pdb_id", "cluster_v2_id"]
    if work.duplicated(keys).any() or marked.duplicated(keys).any():
        raise ValueError("weight sensitivity requires unique target/cluster keys")
    lookup = marked.set_index(keys)["automated_evidence_label"]
    work_keys = pd.MultiIndex.from_frame(work[keys])
    if set(lookup.index) != set(work_keys):
        raise ValueError("weight-sensitivity label coverage does not match rankings")
    work["label"] = lookup.reindex(work_keys).to_numpy()
    work["Within_PDB_rank"] = pd.to_numeric(work["Within_PDB_rank"], errors="raise")
    base_weights = {module: float(module_weights[module]) for module in STATIC_MODULES}
    _weighted_available(work, base_weights)  # validate weights and module values

    target_state: dict[str, dict[str, object]] = {}
    for pdb_id, group in work.groupby("pdb_id", sort=True):
        ordered = group.sort_values(["Within_PDB_rank", "cluster_v2_id"], kind="mergesort")
        references = ordered.loc[ordered["label"].eq("R_auto"), "Within_PDB_rank"]
        supported = ordered.loc[ordered["label"].isin(["R_auto", "A_auto"]), "Within_PDB_rank"]
        target_state[pdb_id] = {
            "baseline_top": str(ordered.iloc[0]["cluster_v2_id"]),
            "baseline_top3": set(ordered.loc[ordered["Within_PDB_rank"].le(3), "cluster_v2_id"]),
            "reference": float(references.min()) if len(references) else math.nan,
            "supported": float(supported.min()) if len(supported) else math.nan,
            "retained": 0,
            "changed": [],
            "reference_ranks": [],
            "supported_ranks": [],
        }

    scenario_rows: list[dict[str, object]] = []
    for module in STATIC_MODULES:
        for direction, multiplier in (("decrease", 1.0 - perturbation_fraction), ("increase", 1.0 + perturbation_fraction)):
            scenario_id = f"{module}_{direction}_20pct"
            weights = dict(base_weights)
            weights[module] *= multiplier
            work["perturbed_score"] = _weighted_available(work, weights)
            work["perturbed_rank"] = work.groupby("pdb_id")["perturbed_score"].rank(
                method=tie_method, ascending=False,
            )
            retained = 0
            jaccards: list[float] = []
            rhos: list[float] = []
            reference_top1 = reference_top3 = supported_top1 = supported_top3 = 0
            for pdb_id, group in work.groupby("pdb_id", sort=True):
                ordered = group.sort_values(["perturbed_rank", "cluster_v2_id"], kind="mergesort")
                state = target_state[pdb_id]
                top = str(ordered.iloc[0]["cluster_v2_id"])
                if top == state["baseline_top"]:
                    retained += 1
                    state["retained"] = int(state["retained"]) + 1
                else:
                    state["changed"].append(scenario_id)
                top3 = set(ordered.loc[ordered["perturbed_rank"].le(3), "cluster_v2_id"])
                union = top3 | state["baseline_top3"]
                jaccards.append(len(top3 & state["baseline_top3"]) / len(union) if union else 1.0)
                baseline = pd.to_numeric(group["Within_PDB_rank"], errors="raise")
                perturbed = pd.to_numeric(group["perturbed_rank"], errors="raise")
                rho = baseline.corr(perturbed, method="spearman")
                rhos.append(float(rho))
                refs = group.loc[group["label"].eq("R_auto"), "perturbed_rank"]
                supp = group.loc[group["label"].isin(["R_auto", "A_auto"]), "perturbed_rank"]
                reference_rank = float(refs.min()) if len(refs) else math.nan
                supported_rank = float(supp.min()) if len(supp) else math.nan
                state["reference_ranks"].append(reference_rank)
                state["supported_ranks"].append(supported_rank)
                reference_top1 += int(not math.isnan(reference_rank) and reference_rank <= 1)
                reference_top3 += int(not math.isnan(reference_rank) and reference_rank <= 3)
                supported_top1 += int(not math.isnan(supported_rank) and supported_rank <= 1)
                supported_top3 += int(not math.isnan(supported_rank) and supported_rank <= 3)
            scenario_rows.append({
                "scenario_id": scenario_id,
                "perturbed_module": module,
                "direction": direction,
                "multiplier": multiplier,
                "target_N": len(target_state),
                "baseline_top1_retained_N": retained,
                "mean_top3_jaccard": float(np.mean(jaccards)),
                "median_spearman_rho": float(np.median(rhos)),
                "Reference_Top1_N": reference_top1,
                "Reference_Top3_N": reference_top3,
                "First_supported_Top1_N": supported_top1,
                "First_supported_Top3_N": supported_top3,
            })

    def finite_min(values: Sequence[float]) -> float:
        finite = [float(value) for value in values if math.isfinite(float(value))]
        return min(finite) if finite else math.nan

    def finite_max(values: Sequence[float]) -> float:
        finite = [float(value) for value in values if math.isfinite(float(value))]
        return max(finite) if finite else math.nan

    target_rows = []
    for pdb_id, state in sorted(target_state.items()):
        target_rows.append({
            "pdb_id": pdb_id,
            "baseline_top_cluster_v2_id": state["baseline_top"],
            "top1_retention_count": state["retained"],
            "top1_changed_scenarios": ";".join(state["changed"]),
            "baseline_reference_rank": state["reference"],
            "minimum_reference_rank": finite_min(state["reference_ranks"]),
            "maximum_reference_rank": finite_max(state["reference_ranks"]),
            "baseline_first_supported_rank": state["supported"],
            "minimum_first_supported_rank": finite_min(state["supported_ranks"]),
            "maximum_first_supported_rank": finite_max(state["supported_ranks"]),
        })
    return WeightSensitivityResult(
        scenarios=pd.DataFrame(scenario_rows, columns=WEIGHT_SCENARIO_COLUMNS),
        targets=pd.DataFrame(target_rows, columns=WEIGHT_TARGET_COLUMNS),
    )


def _weighted_without_orel(frame: pd.DataFrame, weights: Mapping[str, float]) -> pd.Series:
    if set(weights) != set(STATIC_MODULES):
        raise ValueError(f"module weights must contain exactly {list(STATIC_MODULES)}")
    numeric_weights = {key: float(value) for key, value in weights.items()}
    if any(not math.isfinite(value) or value <= 0 for value in numeric_weights.values()):
        raise ValueError("module weights must be finite and positive")
    numerator = np.zeros(len(frame), dtype=float)
    denominator = np.zeros(len(frame), dtype=float)
    for module in STATIC_MODULES:
        if module == "O_rel_formal":
            continue
        values = pd.to_numeric(frame[module], errors="coerce").to_numpy(float)
        valid = np.isfinite(values)
        numerator[valid] += numeric_weights[module] * values[valid]
        denominator[valid] += numeric_weights[module]
    result = np.full(len(frame), np.nan)
    valid = denominator > 0
    result[valid] = numerator[valid] / denominator[valid]
    return pd.Series(result, index=frame.index)


def orel_ablation(
    rankings: pd.DataFrame,
    labels: pd.DataFrame,
    reference_evaluable: pd.Series | Mapping[str, object],
    systems: pd.DataFrame,
    module_weights: Mapping[str, float],
    *,
    tie_method: str = "average",
) -> AblationResult:
    """Remove only O-rel and return physically separate ablation tables."""
    _require(rankings, ("pdb_id", "cluster_v2_id", "Within_PDB_rank", *STATIC_MODULES), "rankings")
    _require(labels, ("pdb_id", "cluster_v2_id", "automated_evidence_label"), "labels")
    if tie_method != "average":
        raise ValueError("O-rel ablation requires average tie ranking")
    work = rankings.copy(deep=True)
    work["pdb_id"] = work["pdb_id"].astype(str).str.upper()
    work["Within_PDB_rank"] = pd.to_numeric(work["Within_PDB_rank"], errors="raise")
    labels_work = labels.copy(deep=True)
    labels_work["pdb_id"] = labels_work["pdb_id"].astype(str).str.upper()
    keys = ["pdb_id", "cluster_v2_id"]
    if work.duplicated(keys).any() or labels_work.duplicated(keys).any():
        raise ValueError("ablation inputs require unique target/cluster keys")
    label_lookup = labels_work.set_index(keys)["automated_evidence_label"]
    work_keys = pd.MultiIndex.from_frame(work[keys])
    if set(label_lookup.index) != set(work_keys):
        raise ValueError("ablation label coverage does not match rankings")
    work["label"] = label_lookup.reindex(work_keys).to_numpy()
    work["Without_O_rel_score"] = _weighted_without_orel(work, module_weights)
    work["Without_O_rel_rank"] = work.groupby("pdb_id")["Without_O_rel_score"].rank(
        method=tie_method, ascending=False,
    )
    metadata = _metadata(systems)
    if set(metadata.index) != set(work["pdb_id"]):
        raise ValueError("systems target coverage does not match rankings")
    reference_flags = _reference_flags(reference_evaluable, set(work["pdb_id"]))

    rows: list[dict[str, object]] = []
    for pdb_id, group in work.groupby("pdb_id", sort=True):
        references = group.loc[group["label"].eq("R_auto")]
        supported = group.loc[group["label"].isin(["R_auto", "A_auto"])]
        full_reference = float(references["Within_PDB_rank"].min()) if len(references) else math.nan
        no_reference = float(references["Without_O_rel_rank"].min()) if len(references) else math.nan
        full_supported = float(supported["Within_PDB_rank"].min()) if len(supported) else math.nan
        no_supported = float(supported["Without_O_rel_rank"].min()) if len(supported) else math.nan
        full_top = group.sort_values(["Within_PDB_rank", "cluster_v2_id"], kind="mergesort").iloc[0]
        no_top = group.sort_values(["Without_O_rel_rank", "cluster_v2_id"], kind="mergesort").iloc[0]
        full_ranks = pd.to_numeric(group["Within_PDB_rank"], errors="raise")
        ablated_ranks = pd.to_numeric(group["Without_O_rel_rank"], errors="raise")
        rho = (
            math.nan
            if full_ranks.nunique(dropna=True) < 2 or ablated_ranks.nunique(dropna=True) < 2
            else full_ranks.corr(ablated_ranks, method="spearman")
        )
        rows.append({
            "analysis_level": "target",
            "pdb_id": pdb_id,
            "Pocket_category": metadata.loc[pdb_id, "Pocket_category"],
            "Protein_family": metadata.loc[pdb_id, "Protein_family"],
            "Full_reference_rank": full_reference,
            "Without_O_rel_reference_rank": no_reference,
            "reference_rank_change_without_minus_full": no_reference - full_reference if len(references) else math.nan,
            "Full_first_supported_rank": full_supported,
            "Without_O_rel_first_supported_rank": no_supported,
            "first_supported_rank_change_without_minus_full": no_supported - full_supported if len(supported) else math.nan,
            "Full_reference_Top3": bool(len(references) and full_reference <= 3),
            "Without_O_rel_reference_Top3": bool(len(references) and no_reference <= 3),
            "reference_Top3_inclusion_change": int(len(references) and no_reference <= 3) - int(len(references) and full_reference <= 3),
            "rank_Spearman_rho": float(rho),
            "Full_top_cluster_v2_id": full_top["cluster_v2_id"],
            "Without_O_rel_top_cluster_v2_id": no_top["cluster_v2_id"],
            "top_cluster_identity_changed": full_top["cluster_v2_id"] != no_top["cluster_v2_id"],
        })
    targets = pd.DataFrame(rows, columns=ABLATION_TARGET_COLUMNS).sort_values("pdb_id", kind="mergesort").reset_index(drop=True)
    category_rows: list[dict[str, object]] = []
    for category, group in targets.groupby("Pocket_category", sort=True):
        evaluable = group.loc[group["pdb_id"].map(reference_flags)]
        full = pd.to_numeric(evaluable["Full_reference_rank"], errors="coerce")
        without = pd.to_numeric(evaluable["Without_O_rel_reference_rank"], errors="coerce")
        category_rows.append({
            "analysis_level": "pocket_archetype",
            "Pocket_category": category,
            "target_N": len(group),
            "Full_reference_Top1": float(full.le(1).mean()) if len(evaluable) else math.nan,
            "Without_O_rel_reference_Top1": float(without.le(1).mean()) if len(evaluable) else math.nan,
            "Full_reference_Top3": float(full.le(3).mean()) if len(evaluable) else math.nan,
            "Without_O_rel_reference_Top3": float(without.le(3).mean()) if len(evaluable) else math.nan,
            "Full_reference_MRR": float((1.0 / full).fillna(0.0).mean()) if len(evaluable) else math.nan,
            "Without_O_rel_reference_MRR": float((1.0 / without).fillna(0.0).mean()) if len(evaluable) else math.nan,
            "mean_rank_Spearman_rho": float(group["rank_Spearman_rho"].mean()),
            "top_cluster_identity_change_fraction": float(group["top_cluster_identity_changed"].mean()),
        })
    categories = pd.DataFrame(category_rows, columns=ABLATION_CATEGORY_COLUMNS)
    clusters = work.loc[:, list(ABLATION_CLUSTER_COLUMNS)].sort_values(
        ["pdb_id", "Without_O_rel_rank", "cluster_v2_id"], kind="mergesort"
    ).reset_index(drop=True)
    return AblationResult(targets, categories, clusters)


def build_representative_cases(
    rankings: pd.DataFrame,
    labels: pd.DataFrame,
    reference_mapping: pd.DataFrame,
    md_mapping: pd.DataFrame,
    redocking_mapping: pd.DataFrame,
    ablation_targets: pd.DataFrame,
    case_ids: Sequence[str],
) -> pd.DataFrame:
    """Materialize the configured representative cases after ablation."""
    ranking_required = (
        "pdb_id", "cluster_v2_id", "Within_PDB_rank", "OIPS-P_static_recomputed",
        "tool_support_count", "contributing_chains", "interface_fraction",
    )
    _require(rankings, ranking_required, "rankings")
    _require(labels, ("pdb_id", "cluster_v2_id", "automated_evidence_label"), "labels")
    _require(reference_mapping, ("pdb_id", "cluster_v2_id", "DCC_A"), "reference mapping")
    _require(md_mapping, ("pdb_id", "Concordance_call"), "MD mapping")
    _require(redocking_mapping, ("pdb_id", "Raw_ligand_RMSD_A", "RMSD_threshold_call"), "redocking mapping")
    _require(ablation_targets, ("pdb_id", "Without_O_rel_reference_rank", "Without_O_rel_top_cluster_v2_id"), "ablation targets")
    requested = tuple(str(value).upper() for value in case_ids)
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("representative case IDs must be nonempty and unique")
    ranked = rankings.copy()
    ranked["pdb_id"] = ranked["pdb_id"].astype(str).str.upper()
    label_lookup = labels.assign(pdb_id=labels["pdb_id"].astype(str).str.upper()).set_index(["pdb_id", "cluster_v2_id"])["automated_evidence_label"]
    ref_lookup = reference_mapping.assign(pdb_id=reference_mapping["pdb_id"].astype(str).str.upper()).set_index(["pdb_id", "cluster_v2_id"])
    abl_lookup = ablation_targets.assign(pdb_id=ablation_targets["pdb_id"].astype(str).str.upper()).set_index("pdb_id")
    red_lookup = redocking_mapping.assign(pdb_id=redocking_mapping["pdb_id"].astype(str).str.upper()).set_index("pdb_id")
    md = md_mapping.assign(pdb_id=md_mapping["pdb_id"].astype(str).str.upper())
    rows: list[dict[str, object]] = []
    for pdb_id in requested:
        group = ranked.loc[ranked["pdb_id"].eq(pdb_id)].sort_values(["Within_PDB_rank", "cluster_v2_id"], kind="mergesort")
        if group.empty or pdb_id not in abl_lookup.index or pdb_id not in red_lookup.index:
            raise ValueError(f"representative case is missing required inputs: {pdb_id}")
        top = group.iloc[0]
        labeled = group.copy()
        labeled["label"] = [label_lookup.loc[(pdb_id, value)] for value in labeled["cluster_v2_id"]]
        references = labeled.loc[labeled["label"].eq("R_auto")].sort_values(["Within_PDB_rank", "cluster_v2_id"], kind="mergesort")
        reference = references.iloc[0] if len(references) else None
        ref_metrics = ref_lookup.loc[(pdb_id, reference["cluster_v2_id"])] if reference is not None else None
        calls = sorted(set(md.loc[md["pdb_id"].eq(pdb_id) & ~md["Concordance_call"].eq("apo_only_context"), "Concordance_call"].astype(str)))
        calls = [value for value in calls if value]
        red = red_lookup.loc[pdb_id]
        abl = abl_lookup.loc[pdb_id]
        rows.append({
            "pdb_id": pdb_id,
            "static_top_cluster_v2_id": top["cluster_v2_id"],
            "static_top_score": top["OIPS-P_static_recomputed"],
            "static_top_label": label_lookup.loc[(pdb_id, top["cluster_v2_id"])],
            "static_top_tool_support": top["tool_support_count"],
            "static_top_chains": top["contributing_chains"],
            "static_top_interface_overlap": top["interface_fraction"],
            "reference_cluster_v2_id": reference["cluster_v2_id"] if reference is not None else "",
            "reference_rank": reference["Within_PDB_rank"] if reference is not None else math.nan,
            "reference_DCC_A": ref_metrics["DCC_A"] if ref_metrics is not None else math.nan,
            "without_O_rel_reference_rank": abl["Without_O_rel_reference_rank"],
            "without_O_rel_top_cluster_v2_id": abl["Without_O_rel_top_cluster_v2_id"],
            "MD_calls": ";".join(calls) if calls else "MD_not_available",
            "redocking_RMSD_A": red["Raw_ligand_RMSD_A"],
            "redocking_status": red["RMSD_threshold_call"],
        })
    return pd.DataFrame(rows, columns=REPRESENTATIVE_COLUMNS)
