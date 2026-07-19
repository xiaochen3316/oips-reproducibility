"""Deterministic, leakage-free reconstruction of OIPS cluster-v2 candidates."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
import json
import math

import numpy as np
import pandas as pd


INPUT_COLUMNS = (
    "row_id", "pdb_id", "tool", "pocket_id", "display_order",
    "sitemap_rank", "center_x", "center_y", "center_z", "center_method",
    "residue_count", "residue_set_json", "pocket_geometry_score",
    "pocket_ligandability_score",
)


@dataclass(frozen=True)
class SameToolConfig:
    residue_iou_strong: float = 0.70
    residue_containment_nested: float = 0.85
    nested_center_max_A: float = 8.0
    near_center_max_A: float = 2.5
    near_center_min_iou: float = 0.25
    center_only_duplicate_max_A: float = 1.5
    dogsite_hierarchy_center_max_A: float = 6.0
    dogsite_hierarchy_min_iou: float = 0.35
    dogsite_hierarchy_containment_min: float = 0.70
    group_center_diameter_max_A: float = 8.0
    representative_missing_center_penalty: float = 0.20
    representative_missing_residue_penalty: float = 0.20


@dataclass(frozen=True)
class CrossToolConfig:
    close_center_max_A: float = 6.0
    conditional_center_max_A: float = 10.0
    conditional_residue_iou_min: float = 0.20
    residue_iou_min: float = 0.35
    cluster_center_diameter_max_A: float = 12.0


@dataclass(frozen=True)
class Unit:
    pdb_id: str
    tool: str
    unit_id: str
    representative_row_id: int
    representative_pocket_id: str
    raw_row_ids: tuple[int, ...]
    raw_pocket_ids: tuple[str, ...]
    center: tuple[float, float, float] | None
    residues: frozenset[str]
    geometry: float
    ligandability: float
    native_rank: float


@dataclass(frozen=True)
class ClusteringResult:
    candidates: pd.DataFrame
    membership: pd.DataFrame
    mapping: pd.DataFrame
    excluded: pd.DataFrame
    boundary: pd.DataFrame
    units: tuple[Unit, ...] = ()
    features: pd.DataFrame | None = None
    tool_weights: Mapping[str, float] | None = None


def _as_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def parse_residues(value: object) -> set[str]:
    """Parse a JSON array or compatibility semicolon list into residue IDs."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return set()
    text = str(value)
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return {item.strip() for item in text.split(";") if item.strip()}
    if not isinstance(parsed, list):
        raise ValueError("residue input must be a JSON array or semicolon list")
    return {str(item) for item in parsed if str(item)}


def simple_residue(value: str) -> str:
    parts = str(value).split(":")
    if len(parts) >= 3:
        insertion = parts[3] if len(parts) > 3 else ""
        return f"{parts[0]}:{parts[2]}:{insertion}"
    return str(value)


def simple_set(values: Iterable[str]) -> set[str]:
    return {simple_residue(value) for value in values if value}


def residue_iou(a: set[str] | frozenset[str], b: set[str] | frozenset[str]) -> float:
    aa, bb = simple_set(a), simple_set(b)
    if not aa or not bb:
        return math.nan
    return len(aa & bb) / len(aa | bb)


def residue_containment(
    a: set[str] | frozenset[str], b: set[str] | frozenset[str]
) -> float:
    aa, bb = simple_set(a), simple_set(b)
    if not aa or not bb:
        return math.nan
    return len(aa & bb) / min(len(aa), len(bb))


def _center_from_row(row: Mapping[str, object] | pd.Series) -> tuple[float, float, float] | None:
    values = [_as_float(row.get(column)) for column in ("center_x", "center_y", "center_z")]
    if any(math.isnan(value) for value in values):
        return None
    return values[0], values[1], values[2]


def distance(
    a: tuple[float, float, float] | None,
    b: tuple[float, float, float] | None,
) -> float:
    if a is None or b is None:
        return math.nan
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def max_center_diameter(records: list[Unit] | list[dict[str, object]]) -> float:
    centers = [record.center if isinstance(record, Unit) else record.get("center") for record in records]
    present = [center for center in centers if center is not None]
    if len(present) <= 1:
        return 0.0
    return max(
        distance(present[i], present[j])
        for i in range(len(present))
        for j in range(i + 1, len(present))
    )


def dogsite_root(pocket_id: str) -> str:
    parts = str(pocket_id).split("_")
    if len(parts) >= 2 and parts[0].upper() == "P" and parts[1].isdigit():
        return f"P_{parts[1]}"
    return str(pocket_id)


def native_rank(row: Mapping[str, object] | pd.Series) -> float:
    tool, pocket = str(row.get("tool", "")), str(row.get("pocket_id", ""))
    if tool == "SiteMap":
        rank = _as_float(row.get("sitemap_rank"))
        if math.isfinite(rank):
            return rank
    if tool in {"DoGSiteScorer", "DoGSite3"}:
        try:
            root_index = int(dogsite_root(pocket).split("_")[1])
            return float(root_index + (1 if tool == "DoGSiteScorer" else 0))
        except (IndexError, ValueError):
            pass
    order = _as_float(row.get("display_order"))
    if math.isfinite(order):
        return order
    digits = "".join(character for character in pocket if character.isdigit())
    return float(digits) if digits else 1e9


def same_tool_pair_is_duplicate(
    a: Mapping[str, object], b: Mapping[str, object], config: SameToolConfig
) -> bool:
    if a["tool"] != b["tool"]:
        return False
    d = distance(a["center"], b["center"])  # type: ignore[arg-type]
    iou = residue_iou(a["residues"], b["residues"])  # type: ignore[arg-type]
    containment = residue_containment(a["residues"], b["residues"])  # type: ignore[arg-type]
    hierarchical = (
        a["tool"] in {"DoGSiteScorer", "DoGSite3"}
        and dogsite_root(str(a["pocket_id"])) == dogsite_root(str(b["pocket_id"]))
    )
    if hierarchical:
        if not math.isnan(d) and d <= config.dogsite_hierarchy_center_max_A:
            return True
        if not math.isnan(iou) and iou >= config.dogsite_hierarchy_min_iou:
            return True
        if not math.isnan(containment) and containment >= config.dogsite_hierarchy_containment_min:
            return True
    if not math.isnan(iou) and iou >= config.residue_iou_strong:
        return True
    if not math.isnan(containment) and containment >= config.residue_containment_nested:
        return math.isnan(d) or d <= config.nested_center_max_A
    if not math.isnan(d) and d <= config.near_center_max_A:
        if not math.isnan(iou) and iou >= config.near_center_min_iou:
            return True
        if not a["residues"] and not b["residues"]:
            return d <= config.center_only_duplicate_max_A
    return False


def _raw_pair_similarity(a: Mapping[str, object], b: Mapping[str, object]) -> float:
    d = distance(a["center"], b["center"])  # type: ignore[arg-type]
    iou = residue_iou(a["residues"], b["residues"])  # type: ignore[arg-type]
    parts: list[float] = []
    if not math.isnan(d):
        parts.append(max(0.0, 1.0 - d / 12.0))
    if not math.isnan(iou):
        parts.append(iou)
    return float(np.mean(parts)) if parts else 0.0


def _group_same_tool(records: list[dict[str, object]], config: SameToolConfig) -> list[list[dict[str, object]]]:
    ordered = sorted(records, key=lambda record: (native_rank(record), int(record["row_id"])))
    groups: list[list[dict[str, object]]] = []
    for record in ordered:
        candidates: list[tuple[float, int]] = []
        for index, group in enumerate(groups):
            if all(same_tool_pair_is_duplicate(record, member, config) for member in group):
                if max_center_diameter(group + [record]) <= config.group_center_diameter_max_A:
                    candidates.append((float(np.mean([_raw_pair_similarity(record, member) for member in group])), index))
        if not candidates:
            groups.append([record])
        else:
            _, best = max(candidates, key=lambda item: (item[0], -item[1]))
            groups[best].append(record)
    return groups


def _representative_cost(
    candidate: Mapping[str, object], group: list[dict[str, object]], config: SameToolConfig
) -> tuple[float, float, int, int]:
    costs: list[float] = []
    for other in group:
        if other["row_id"] == candidate["row_id"]:
            continue
        d = distance(candidate["center"], other["center"])  # type: ignore[arg-type]
        iou = residue_iou(candidate["residues"], other["residues"])  # type: ignore[arg-type]
        terms: list[float] = []
        if not math.isnan(d):
            terms.append(min(2.0, d / 8.0))
        if not math.isnan(iou):
            terms.append(1.0 - iou)
        costs.append(float(np.mean(terms)) if terms else 1.0)
    cost = float(np.mean(costs)) if costs else 0.0
    if candidate["center"] is None:
        cost += config.representative_missing_center_penalty
    if not candidate["residues"]:
        cost += config.representative_missing_residue_penalty
    hierarchy_depth = max(0, str(candidate["pocket_id"]).count("_") - 1)
    return cost, native_rank(candidate), hierarchy_depth, int(candidate["row_id"])


def _build_units(
    features: pd.DataFrame, config: SameToolConfig
) -> tuple[list[Unit], pd.DataFrame, pd.DataFrame]:
    units: list[Unit] = []
    mapping_rows: list[dict[str, object]] = []
    excluded_rows: list[dict[str, object]] = []
    work = features.copy()
    residues_by_index = work["residue_set_json"].map(parse_residues)
    work["_has_center"] = work[["center_x", "center_y", "center_z"]].notna().all(axis=1)
    work["_has_residues"] = residues_by_index.map(bool)
    work["_mappable"] = work["_has_center"] | work["_has_residues"]
    for row in work.loc[~work["_mappable"]].itertuples(index=False):
        excluded_rows.append({
            "row_id": int(row.row_id), "pdb_id": str(row.pdb_id), "tool": str(row.tool),
            "pocket_id": str(row.pocket_id), "center_method": str(row.center_method),
            "residue_count": int(row.residue_count or 0),
            "exclusion_reason": "unmappable_no_center_and_no_residue_set",
            "retained_in_audit": True,
        })
    mapped = work.loc[work["_mappable"]]
    for (pdb_id, tool), subset in mapped.groupby(["pdb_id", "tool"], sort=True):
        records: list[dict[str, object]] = []
        for _, row in subset.sort_values("row_id").iterrows():
            records.append({
                "row_id": int(row["row_id"]), "pdb_id": str(pdb_id), "tool": str(tool),
                "pocket_id": str(row["pocket_id"]), "center": _center_from_row(row),
                "residues": parse_residues(row["residue_set_json"]),
                "geometry": _as_float(row["pocket_geometry_score"]),
                "ligandability": _as_float(row["pocket_ligandability_score"]),
                "display_order": row["display_order"], "sitemap_rank": row["sitemap_rank"],
            })
        groups = sorted(_group_same_tool(records, config), key=lambda group: min(int(item["row_id"]) for item in group))
        for number, group in enumerate(groups, start=1):
            representative = min(group, key=lambda item: _representative_cost(item, group, config))
            unit_id = f"{pdb_id}_{tool}_U{number:03d}"
            unit = Unit(
                pdb_id=str(pdb_id), tool=str(tool), unit_id=unit_id,
                representative_row_id=int(representative["row_id"]),
                representative_pocket_id=str(representative["pocket_id"]),
                raw_row_ids=tuple(sorted(int(item["row_id"]) for item in group)),
                raw_pocket_ids=tuple(str(item["pocket_id"]) for item in sorted(group, key=lambda item: int(item["row_id"]))),
                center=representative["center"], residues=frozenset(representative["residues"]),  # type: ignore[arg-type]
                geometry=float(representative["geometry"]), ligandability=float(representative["ligandability"]),
                native_rank=float(native_rank(representative)),
            )
            units.append(unit)
            for raw in group:
                mapping_rows.append({
                    "row_id": int(raw["row_id"]), "pdb_id": str(pdb_id), "tool": str(tool),
                    "pocket_id": str(raw["pocket_id"]), "same_tool_unit_id": unit_id,
                    "same_tool_group_size": len(group),
                    "representative_for_tool_unit": int(raw["row_id"]) == unit.representative_row_id,
                    "representative_pocket_id": unit.representative_pocket_id,
                })
    unit_mapping = pd.DataFrame(mapping_rows, columns=(
        "row_id", "pdb_id", "tool", "pocket_id", "same_tool_unit_id",
        "same_tool_group_size", "representative_for_tool_unit", "representative_pocket_id",
    ))
    excluded = pd.DataFrame(excluded_rows, columns=(
        "row_id", "pdb_id", "tool", "pocket_id", "center_method", "residue_count",
        "exclusion_reason", "retained_in_audit",
    ))
    return units, unit_mapping, excluded


def site_pair_compatible(a: Unit, b: Unit, config: CrossToolConfig) -> bool:
    d, iou = distance(a.center, b.center), residue_iou(a.residues, b.residues)
    if not math.isnan(d) and d <= config.close_center_max_A:
        return True
    if (not math.isnan(d) and d <= config.conditional_center_max_A
            and not math.isnan(iou) and iou >= config.conditional_residue_iou_min):
        return True
    return not math.isnan(iou) and iou >= config.residue_iou_min


def _unit_similarity(a: Unit, b: Unit, config: CrossToolConfig) -> float:
    d, iou = distance(a.center, b.center), residue_iou(a.residues, b.residues)
    parts: list[float] = []
    if not math.isnan(d):
        parts.append(max(0.0, 1.0 - d / config.cluster_center_diameter_max_A))
    if not math.isnan(iou):
        parts.append(iou)
    return float(np.mean(parts)) if parts else 0.0


def medoid_unit(cluster: list[Unit]) -> Unit:
    candidates = [unit for unit in cluster if unit.center is not None] or cluster
    def cost(candidate: Unit) -> tuple[float, int]:
        pair_costs: list[float] = []
        for other in cluster:
            if other.unit_id == candidate.unit_id:
                continue
            d, iou = distance(candidate.center, other.center), residue_iou(candidate.residues, other.residues)
            terms: list[float] = []
            if not math.isnan(d):
                terms.append(min(2.0, d / 10.0))
            if not math.isnan(iou):
                terms.append(1.0 - iou)
            pair_costs.append(float(np.mean(terms)) if terms else 1.0)
        return (float(np.mean(pair_costs)) if pair_costs else 0.0, candidate.representative_row_id)
    return min(candidates, key=cost)


def medoid_constrained_group_is_valid(cluster: list[Unit], config: CrossToolConfig) -> bool:
    if len(cluster) <= 1:
        return True
    if max_center_diameter(cluster) > config.cluster_center_diameter_max_A:
        return False
    medoid = medoid_unit(cluster)
    return all(unit.unit_id == medoid.unit_id or site_pair_compatible(unit, medoid, config) for unit in cluster)


def _cross_tool_pair_compatible(a: Unit, b: Unit, config: CrossToolConfig) -> bool:
    return a.tool != b.tool and site_pair_compatible(a, b, config)


def _build_cross_tool_clusters(
    units: list[Unit], config: CrossToolConfig, tool_weights: Mapping[str, float]
) -> dict[str, list[list[Unit]]]:
    by_target: dict[str, list[list[Unit]]] = {}
    for pdb_id in sorted({unit.pdb_id for unit in units}):
        target = [unit for unit in units if unit.pdb_id == pdb_id]
        target.sort(key=lambda unit: (
            -(int(unit.center is not None) + int(bool(unit.residues))),
            -float(tool_weights.get(unit.tool, 1.0)), unit.native_rank,
            unit.representative_row_id,
        ))
        clusters: list[list[Unit]] = []
        for unit in target:
            candidates: list[tuple[float, int]] = []
            for index, cluster in enumerate(clusters):
                if (any(_cross_tool_pair_compatible(unit, other, config) for other in cluster)
                        and medoid_constrained_group_is_valid(cluster + [unit], config)):
                    candidates.append((float(np.mean([_unit_similarity(unit, other, config) for other in cluster])), index))
            if not candidates:
                clusters.append([unit])
            else:
                _, best = max(candidates, key=lambda item: (item[0], -item[1]))
                clusters[best].append(unit)
        while True:
            merges: list[tuple[float, int, int]] = []
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    if not any(_cross_tool_pair_compatible(a, b, config) for a in clusters[i] for b in clusters[j]):
                        continue
                    if not medoid_constrained_group_is_valid(clusters[i] + clusters[j], config):
                        continue
                    cross = [_unit_similarity(a, b, config) for a in clusters[i] for b in clusters[j] if a.tool != b.tool]
                    merges.append((float(np.mean(cross)) if cross else 0.0, i, j))
            if not merges:
                break
            _, i, j = max(merges, key=lambda item: (item[0], -item[1], -item[2]))
            clusters[i] = clusters[i] + clusters[j]
            del clusters[j]
        by_target[pdb_id] = clusters
    return by_target


def formal_representatives(
    cluster: list[Unit], config: CrossToolConfig = CrossToolConfig()
) -> list[Unit]:
    representatives: list[Unit] = []
    for tool in sorted({unit.tool for unit in cluster}):
        candidates = [unit for unit in cluster if unit.tool == tool]
        others = [unit for unit in cluster if unit.tool != tool]
        def key(candidate: Unit) -> tuple[float, int, int, float, int]:
            dissimilarity = float(np.mean([1.0 - _unit_similarity(candidate, other, config) for other in others])) if others else 0.0
            return (dissimilarity, int(candidate.center is None) + int(not candidate.residues),
                    -len(candidate.residues), candidate.native_rank, candidate.representative_row_id)
        representatives.append(min(candidates, key=key))
    return representatives


def residue_consensus(cluster: list[Unit]) -> tuple[set[str], set[str], dict[str, int]]:
    canonical: dict[str, str] = {}
    support: Counter[str] = Counter()
    for unit in cluster:
        present = simple_set(unit.residues)
        support.update(present)
        for residue in sorted(unit.residues):
            canonical.setdefault(simple_residue(residue), residue)
    envelope = {canonical[value] for value in support}
    threshold = 1 if len(cluster) == 1 else max(2, math.ceil(len(cluster) / 2))
    core = {canonical[value] for value, count in support.items() if count >= threshold}
    return core, envelope, dict(support)


def _chain_metrics(residues: set[str]) -> tuple[list[str], float, float]:
    counts = Counter(residue.split(":")[0] for residue in residues if ":" in residue)
    if not counts:
        return [], math.nan, math.nan
    total = sum(counts.values())
    dominant = max(counts.values()) / total
    if len(counts) == 1:
        entropy = 0.0
    else:
        proportions = np.asarray(list(counts.values()), dtype=float) / total
        entropy = float(-(proportions * np.log(proportions)).sum() / math.log(len(proportions)))
    return sorted(counts), dominant, entropy


def _cluster_pair_metrics(cluster: list[Unit]) -> tuple[float, float, float]:
    distances: list[float] = []
    ious: list[float] = []
    for i in range(len(cluster)):
        for j in range(i + 1, len(cluster)):
            d = distance(cluster[i].center, cluster[j].center)
            if not math.isnan(d):
                distances.append(d)
            iou = residue_iou(cluster[i].residues, cluster[j].residues)
            if not math.isnan(iou):
                ious.append(iou)
    return (max(distances) if distances else 0.0,
            float(np.median(ious)) if ious else math.nan,
            min(ious) if ious else math.nan)


def _build_candidate_tables(
    target_clusters: dict[str, list[list[Unit]]], config: CrossToolConfig
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, str]]:
    candidate_rows: list[dict[str, object]] = []
    membership_rows: list[dict[str, object]] = []
    boundary_rows: list[dict[str, object]] = []
    raw_to_cluster: dict[int, str] = {}
    for pdb_id, clusters in sorted(target_clusters.items()):
        keyed: list[tuple[list[Unit], tuple[object, ...]]] = []
        for cluster in clusters:
            reps = formal_representatives(cluster, config)
            medoid = medoid_unit(reps)
            center = medoid.center or (math.inf, math.inf, math.inf)
            keyed.append((cluster, (-len(reps), *center, min(unit.representative_row_id for unit in cluster))))
        for number, cluster in enumerate((item[0] for item in sorted(keyed, key=lambda item: item[1])), start=1):
            cluster_id = f"{pdb_id}_V2C{number:03d}"
            reps = formal_representatives(cluster, config)
            formal_ids = {unit.unit_id for unit in reps}
            medoid = medoid_unit(reps)
            core, envelope, _ = residue_consensus(reps)
            chains, dominant, entropy = _chain_metrics(envelope)
            diameter = max_center_diameter(cluster)
            _, median_iou, minimum_iou = _cluster_pair_metrics(reps)
            if medoid.center is None:
                dispersion = math.nan
            else:
                distances = [distance(medoid.center, unit.center) for unit in reps if unit.center is not None]
                dispersion = float(np.sqrt(np.mean(np.square(distances)))) if distances else 0.0
            core_ratio = len(core) / len(envelope) if envelope else math.nan
            center_count = sum(unit.center is not None for unit in reps)
            residue_count = sum(bool(unit.residues) for unit in reps)
            mappability = ("center_and_residue_mappable" if center_count and residue_count else
                           "center_only_mappable" if center_count else "residue_only_mappable")
            row = {
                "pdb_id": pdb_id, "cluster_v2_id": cluster_id,
                "medoid_unit_id": medoid.unit_id, "medoid_tool": medoid.tool,
                "medoid_pocket_id": medoid.representative_pocket_id,
                "medoid_center_x": medoid.center[0] if medoid.center else math.nan,
                "medoid_center_y": medoid.center[1] if medoid.center else math.nan,
                "medoid_center_z": medoid.center[2] if medoid.center else math.nan,
                "cluster_diameter_A": diameter, "center_dispersion_A": dispersion,
                "tool_support_count": len(reps),
                "supporting_tools": ";".join(sorted(unit.tool for unit in reps)),
                "representative_pockets_per_tool": ";".join(f"{unit.tool}:{unit.representative_pocket_id}" for unit in sorted(reps, key=lambda unit: unit.tool)),
                "same_tool_units": ";".join(sorted(unit.unit_id for unit in cluster)),
                "same_tool_secondary_unit_count": len(cluster) - len(reps),
                "raw_record_count": sum(len(unit.raw_row_ids) for unit in cluster),
                "core_residue_count": len(core), "envelope_residue_count": len(envelope),
                "core_envelope_ratio": core_ratio, "core_residues": ";".join(sorted(core)),
                "envelope_residues": ";".join(sorted(envelope)),
                "contributing_chains": ";".join(chains), "contributing_chain_count": len(chains),
                "dominant_chain_fraction": dominant, "cluster_chain_entropy": entropy,
                "mappability": mappability, "center_available_representatives": center_count,
                "residue_available_representatives": residue_count,
                "pairwise_residue_iou_median": median_iou, "pairwise_residue_iou_min": minimum_iou,
                "spatial_continuity": bool(diameter <= config.cluster_center_diameter_max_A),
                "boundary_sensitive": False,
            }
            candidate_rows.append(row)
            boundary_rows.append({
                "pdb_id": pdb_id, "cluster_v2_id": cluster_id, "tool_support_count": len(reps),
                "same_tool_secondary_unit_count": len(cluster) - len(reps),
                "cluster_diameter_A": diameter, "center_dispersion_A": dispersion,
                "pairwise_residue_iou_median": median_iou, "pairwise_residue_iou_min": minimum_iou,
                "core_envelope_ratio": core_ratio,
                "spatial_continuity": bool(diameter <= config.cluster_center_diameter_max_A),
                "boundary_sensitive": False, "boundary_reason": "",
            })
            formal_by_tool = {unit.tool: unit for unit in reps}
            for unit in cluster:
                for raw_id, pocket_id in zip(unit.raw_row_ids, unit.raw_pocket_ids):
                    raw_to_cluster[raw_id] = cluster_id
                    is_formal = unit.unit_id in formal_ids and raw_id == unit.representative_row_id
                    membership_rows.append({
                        "pdb_id": pdb_id, "cluster_v2_id": cluster_id, "tool": unit.tool,
                        "same_tool_unit_id": unit.unit_id, "raw_row_id": raw_id,
                        "raw_pocket_id": pocket_id, "formal_tool_representative": is_formal,
                        "representative_pocket_id": formal_by_tool[unit.tool].representative_pocket_id,
                        "formal_vote_count": 1 if is_formal else 0,
                    })
    return (pd.DataFrame(candidate_rows), pd.DataFrame(membership_rows),
            pd.DataFrame(boundary_rows), raw_to_cluster)


def _build_mapping(
    features: pd.DataFrame, unit_mapping: pd.DataFrame, excluded: pd.DataFrame,
    raw_to_cluster: Mapping[int, str],
) -> pd.DataFrame:
    core = features[["row_id", "pdb_id", "tool", "pocket_id", "center_x", "center_y", "center_z", "residue_count"]].copy()
    core = core.merge(unit_mapping, on=["row_id", "pdb_id", "tool", "pocket_id"], how="left")
    core["cluster_v2_id"] = core["row_id"].map(raw_to_cluster)
    reasons = dict(zip(excluded.get("row_id", pd.Series(dtype=int)), excluded.get("exclusion_reason", pd.Series(dtype=str))))
    core["mapping_status"] = np.where(core["cluster_v2_id"].notna(), "mapped_to_cluster_v2", "excluded_unmappable")
    core["exclusion_reason"] = core["row_id"].map(reasons)
    return core.sort_values(["pdb_id", "tool", "row_id"], kind="mergesort").reset_index(drop=True)


def build_clusters(
    features: pd.DataFrame, *, same_tool: SameToolConfig,
    cross_tool: CrossToolConfig, tool_weights: Mapping[str, float],
) -> ClusteringResult:
    """Build deterministic cluster-v2 structural candidates from allow-listed fields."""
    missing = [column for column in INPUT_COLUMNS if column not in features]
    if missing:
        raise KeyError(f"missing clustering input columns: {missing}")
    selected = features.loc[:, list(INPUT_COLUMNS)].copy()
    selected["pdb_id"] = selected["pdb_id"].map(str)
    selected["tool"] = selected["tool"].map(str)
    selected["pocket_id"] = selected["pocket_id"].map(str)
    units, unit_mapping, excluded = _build_units(selected, same_tool)
    target_clusters = _build_cross_tool_clusters(units, cross_tool, tool_weights)
    candidates, membership, boundary, raw_to_cluster = _build_candidate_tables(target_clusters, cross_tool)
    mapping = _build_mapping(selected, unit_mapping, excluded, raw_to_cluster)
    candidates = candidates.sort_values(["pdb_id", "cluster_v2_id"], kind="mergesort").reset_index(drop=True)
    membership = membership.sort_values(["pdb_id", "cluster_v2_id", "tool", "raw_row_id"], kind="mergesort").reset_index(drop=True)
    excluded = excluded.sort_values(["pdb_id", "tool", "row_id"], kind="mergesort").reset_index(drop=True)
    boundary = boundary.sort_values(["pdb_id", "cluster_v2_id"], kind="mergesort").reset_index(drop=True)
    if mapping["row_id"].duplicated().any() or len(mapping) != len(selected):
        raise AssertionError("every raw record must appear exactly once in the mapping")
    if int(mapping["cluster_v2_id"].notna().sum()) + len(excluded) != len(selected):
        raise AssertionError("mapped plus excluded records must equal source records")
    votes = membership.groupby(["cluster_v2_id", "tool"])["formal_vote_count"].sum()
    if not votes.empty and (votes > 1).any():
        raise AssertionError("a cluster has more than one formal vote from a tool")
    return ClusteringResult(
        candidates=candidates, membership=membership, mapping=mapping,
        excluded=excluded, boundary=boundary, units=tuple(units),
        features=selected, tool_weights=dict(tool_weights),
    )


def boundary_decision(
    *, tool_support_count: int, diameter: float, dispersion: float,
    core_ratio: float, median_iou: float, diameter_gt_A: float,
    dispersion_gt_A: float, core_envelope_ratio_lt: float,
    median_residue_iou_lt: float,
) -> tuple[bool, str]:
    if tool_support_count < 2:
        return False, ""
    reasons = [
        (diameter > diameter_gt_A, "diameter_gt_9A"),
        (not math.isnan(dispersion) and dispersion > dispersion_gt_A, "dispersion_gt_4A"),
        (not math.isnan(core_ratio) and core_ratio < core_envelope_ratio_lt, "low_core_envelope_ratio"),
        (not math.isnan(median_iou) and median_iou < median_residue_iou_lt, "low_median_residue_iou"),
    ]
    labels = [label for flag, label in reasons if flag]
    return bool(labels), ";".join(labels)


def annotate_boundaries(
    result: ClusteringResult, *, diameter_gt_A: float, dispersion_gt_A: float,
    core_envelope_ratio_lt: float, median_residue_iou_lt: float,
) -> ClusteringResult:
    candidates = result.candidates.copy()
    boundary = result.boundary.copy()
    decisions = [boundary_decision(
        tool_support_count=int(row.tool_support_count), diameter=float(row.cluster_diameter_A),
        dispersion=float(row.center_dispersion_A), core_ratio=float(row.core_envelope_ratio),
        median_iou=float(row.pairwise_residue_iou_median), diameter_gt_A=diameter_gt_A,
        dispersion_gt_A=dispersion_gt_A, core_envelope_ratio_lt=core_envelope_ratio_lt,
        median_residue_iou_lt=median_residue_iou_lt,
    ) for row in boundary.itertuples(index=False)]
    boundary["boundary_sensitive"] = [item[0] for item in decisions]
    boundary["boundary_reason"] = [item[1] for item in decisions]
    flags = dict(zip(boundary["cluster_v2_id"], boundary["boundary_sensitive"]))
    candidates["boundary_sensitive"] = candidates["cluster_v2_id"].map(flags).astype(bool)
    return replace(result, candidates=candidates, boundary=boundary)
