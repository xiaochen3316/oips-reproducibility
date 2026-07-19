"""Pure static OIPS scoring primitives.

Chain entropy has one definition here: normalized Shannon entropy over the counts
of unique observed residues per chain. Empty and single-chain residue collections
have entropy ``0.0``; equally populated chains have entropy ``1.0``.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
import math
from types import MappingProxyType

import numpy as np
import pandas as pd

from .structure import InterfaceProfile, euclidean_distance


@dataclass(frozen=True)
class ScoringConfig:
    """Frozen static and interface-scoring constants from manuscript YAML."""

    module_weights: Mapping[str, float] = field(
        default_factory=lambda: {
            "C_cons": 0.22,
            "G_geo": 0.18,
            "P_lig": 0.24,
            "O_rel_formal": 0.24,
            "Q_evidence": 0.12,
        }
    )
    tool_weights: Mapping[str, float] = field(
        default_factory=lambda: {
            "DoGSiteScorer": 1.00,
            "DoGSite3": 1.00,
            "CavityPlus": 1.00,
            "CASTpFold": 0.75,
            "SiteMap": 1.15,
        }
    )
    geometry_top_n: int = 3
    ligandability_top_n: int = 2
    q_base: float = 30.0
    q_per_representative: float = 9.0
    q_center_fraction: float = 18.0
    q_residue_fraction: float = 14.0
    q_sitemap_bonus: float = 4.0
    ranking_direction: str = "descending"
    ranking_tie_method: str = "average"

    interface_fraction_weight: float = 0.42
    chain_context_weight: float = 0.25
    distance_weight: float = 0.20
    interface_extent_weight: float = 0.13

    distance_near_max_A: float = 4.0
    distance_near_score: float = 100.0
    distance_intermediate_max_A: float = 8.0
    distance_intermediate_start: float = 90.0
    distance_intermediate_loss_per_A: float = 5.0
    distance_far_max_A: float = 16.0
    distance_far_start: float = 70.0
    distance_far_loss_per_A: float = 6.0
    distance_tail_start: float = 22.0
    distance_tail_loss_per_A: float = 1.2
    distance_tail_floor: float = 10.0

    interface_fraction_base: float = 25.0
    interface_fraction_multiplier: float = 130.0
    multi_chain_base: float = 75.0
    chain_entropy_multiplier: float = 25.0
    single_chain_interface_fraction_min: float = 0.35
    single_chain_supported_score: float = 65.0
    single_chain_other_score: float = 30.0
    context_base: float = 45.0
    context_per_extra_chain: float = 10.0
    context_extra_chain_cap: int = 4
    context_interface_fraction_cap: float = 0.5
    context_interface_fraction_multiplier: float = 80.0
    no_interface_monomer_score: float = 35.0
    no_interface_multichain_score: float = 45.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "module_weights", MappingProxyType(dict(self.module_weights)))
        object.__setattr__(self, "tool_weights", MappingProxyType(dict(self.tool_weights)))


@dataclass(frozen=True)
class ScoringResult:
    """The exact static master and within-target ranking tables."""

    master: "pd.DataFrame"
    rankings: "pd.DataFrame"
    maximum_recomputation_difference: float


def _canonical_residue_set(residues: Iterable[str]) -> frozenset[str]:
    return frozenset(str(residue) for residue in residues if residue)


def _simple_residue(residue: str) -> str:
    parts = residue.split(":")
    if len(parts) >= 3:
        insertion_code = parts[3] if len(parts) > 3 else ""
        return f"{parts[0]}:{parts[2]}:{insertion_code}"
    return residue


def normalized_chain_entropy(residues: Iterable[str]) -> float:
    """Return normalized Shannon entropy of unique residue counts by chain.

    Empty and single-chain inputs return ``0.0`` rather than NaN. Residue names
    do not affect chain membership, and repeated residue IDs count only once.
    """

    unique_residues = _canonical_residue_set(residues)
    counts = Counter(
        residue.split(":", 1)[0]
        for residue in unique_residues
        if ":" in residue
    )
    if len(counts) <= 1:
        return 0.0
    total = sum(counts.values())
    proportions = (count / total for count in counts.values())
    return -sum(proportion * math.log(proportion) for proportion in proportions) / math.log(len(counts))


def _clamp(value: float, low: float = 0.0, high: float = 100.0) -> float:
    return max(low, min(high, value))


def _distance_to_interface(
    center: Sequence[float] | None,
    interface_coordinates: Sequence[Sequence[float]],
) -> float:
    if center is None or not interface_coordinates:
        return math.nan
    return min(euclidean_distance(center, coordinate) for coordinate in interface_coordinates)


def _distance_score(distance: float, config: ScoringConfig) -> float:
    if math.isnan(distance):
        return math.nan
    if distance <= config.distance_near_max_A:
        return config.distance_near_score
    if distance <= config.distance_intermediate_max_A:
        return config.distance_intermediate_start - config.distance_intermediate_loss_per_A * (
            distance - config.distance_near_max_A
        )
    if distance <= config.distance_far_max_A:
        return config.distance_far_start - config.distance_far_loss_per_A * (
            distance - config.distance_intermediate_max_A
        )
    return max(
        config.distance_tail_floor,
        config.distance_tail_start
        - config.distance_tail_loss_per_A * (distance - config.distance_far_max_A),
    )


def weighted_module_score(
    modules: Mapping[str, float | None],
    weights: Mapping[str, float],
) -> float:
    """Return the positive-weight mean of available modules.

    ``None`` and numeric NaN values are unavailable. The denominator is
    renormalized over available modules; all-missing input returns numeric NaN.
    """

    numerator = 0.0
    denominator = 0.0
    for module, raw_weight in weights.items():
        weight = float(raw_weight)
        if weight <= 0.0 or math.isnan(weight):
            continue
        value = modules.get(module)
        if value is None:
            continue
        numeric_value = float(value)
        if math.isnan(numeric_value):
            continue
        numerator += weight * numeric_value
        denominator += weight
    return numerator / denominator if denominator else math.nan


def score_oligomer_relevance(
    residues: Iterable[str],
    center: Sequence[float] | None,
    interface: InterfaceProfile,
    *,
    config: ScoringConfig = ScoringConfig(),
) -> dict[str, float | int]:
    """Score a residue-defined pocket against a static oligomer interface.

    Empty ``residues`` remain empty: their interface fraction is NaN, recall is
    zero when an interface exists, chain count and entropy are zero, and the
    formal score uses the configured fraction baseline with NA renormalization.
    """

    cluster_residues = _canonical_residue_set(residues)
    simple_cluster = {_simple_residue(residue) for residue in cluster_residues}
    simple_interface = {
        _simple_residue(residue) for residue in interface.interface_residues
    }
    overlap = simple_cluster & simple_interface

    interface_fraction = (
        len(overlap) / len(simple_cluster) if simple_cluster else math.nan
    )
    interface_recall = (
        len(overlap) / len(simple_interface) if simple_interface else math.nan
    )
    cluster_chain_count = len(
        {
            residue.split(":", 1)[0]
            for residue in cluster_residues
            if ":" in residue
        }
    )
    entropy = normalized_chain_entropy(cluster_residues)
    distance = _distance_to_interface(center, interface.interface_atom_coordinates)
    distance_score = _distance_score(distance, config)
    context_fraction = (
        len(simple_interface) / interface.protein_residue_count
        if interface.protein_residue_count
        else 0.0
    )
    context_score = _clamp(
        config.context_base
        + config.context_per_extra_chain
        * min(config.context_extra_chain_cap, interface.chain_count - 1)
        + config.context_interface_fraction_multiplier
        * min(config.context_interface_fraction_cap, context_fraction)
    )

    common: dict[str, float | int] = {
        "interface_fraction": interface_fraction,
        "interface_recall": interface_recall,
        "cluster_interface_residue_count": len(overlap),
        "cluster_chain_count": cluster_chain_count,
        "cluster_chain_entropy": entropy,
        "distance_to_interface_A": distance,
        "interface_distance_score": distance_score,
        "context_score": context_score,
    }

    if not simple_interface or interface.chain_count < 2:
        common["O_rel_formal"] = (
            config.no_interface_monomer_score
            if interface.chain_count < 2
            else config.no_interface_multichain_score
        )
        return {"O_rel_formal": common.pop("O_rel_formal"), **common}

    fraction_score = _clamp(
        config.interface_fraction_base
        + config.interface_fraction_multiplier
        * (0.0 if math.isnan(interface_fraction) else interface_fraction)
    )
    if cluster_chain_count >= 2:
        chain_score = config.multi_chain_base + config.chain_entropy_multiplier * entropy
    elif not math.isnan(interface_fraction) and (
        interface_fraction >= config.single_chain_interface_fraction_min
    ):
        chain_score = config.single_chain_supported_score
    else:
        chain_score = config.single_chain_other_score

    formal = weighted_module_score(
        {
            "interface_fraction": fraction_score,
            "chain_context": chain_score,
            "distance": distance_score,
            "interface_extent": context_score,
        },
        {
            "interface_fraction": config.interface_fraction_weight,
            "chain_context": config.chain_context_weight,
            "distance": config.distance_weight,
            "interface_extent": config.interface_extent_weight,
        },
    )
    return {"O_rel_formal": formal, **common}


def _top_mean(values: Iterable[object], count: int) -> float:
    available: list[float] = []
    for value in values:
        if value is None:
            continue
        numeric = float(value)
        if not math.isnan(numeric):
            available.append(numeric)
    return float(np.mean(sorted(available, reverse=True)[:count])) if available else math.nan


def _evidence_quality(representatives: pd.DataFrame, config: ScoringConfig) -> float:
    count = len(representatives)
    if count == 0:
        return math.nan
    center_fraction = float(
        representatives[["center_x", "center_y", "center_z"]].notna().all(axis=1).mean()
    )
    residue_fraction = float((representatives["residue_count"] > 0).mean())
    return _clamp(
        config.q_base
        + config.q_per_representative * count
        + config.q_center_fraction * center_fraction
        + config.q_residue_fraction * residue_fraction
        + (
            config.q_sitemap_bonus
            if representatives["tool"].eq("SiteMap").any()
            else 0.0
        )
    )


def _consensus_score(representatives: pd.DataFrame, config: ScoringConfig) -> float:
    denominator = sum(float(weight) for weight in config.tool_weights.values())
    if denominator <= 0.0:
        return math.nan
    support = sum(
        float(config.tool_weights.get(str(tool), 1.0))
        for tool in representatives["tool"]
    )
    single_tool_factor = 0.5 if len(representatives) == 1 else 1.0
    return 100.0 * support / denominator * single_tool_factor


def score_clusters(
    clustering,
    *,
    interfaces: Mapping[str, InterfaceProfile],
    config: ScoringConfig,
) -> ScoringResult:
    """Score serialized cluster-stage candidates without rerunning clustering."""
    from .clustering import parse_residues
    from .io import MASTER_COLUMNS

    if clustering.features is None:
        raise ValueError("cluster scoring requires the public feature table")
    features = clustering.features.copy()
    membership = clustering.membership.copy()
    formal = membership.loc[membership["formal_tool_representative"].eq(True)].copy()
    if formal.empty and not clustering.candidates.empty:
        raise ValueError("cluster membership contains no formal representatives")
    representatives = formal.merge(
        features,
        left_on=["pdb_id", "tool", "raw_row_id"],
        right_on=["pdb_id", "tool", "row_id"],
        how="left",
        validate="one_to_one",
    )
    if representatives["row_id"].isna().any():
        raise ValueError("a formal representative could not be rejoined to public features")
    vote_sizes = representatives.groupby(["cluster_v2_id", "tool"]).size()
    if not vote_sizes.empty and (vote_sizes > 1).any():
        raise ValueError("more than one formal representative exists for a cluster and tool")

    master_rows: list[dict[str, object]] = []
    for candidate in clustering.candidates.sort_values(
        ["pdb_id", "cluster_v2_id"], kind="mergesort"
    ).to_dict("records"):
        pdb_id = str(candidate["pdb_id"])
        cluster_id = str(candidate["cluster_v2_id"])
        if pdb_id not in interfaces:
            raise KeyError(f"interface profile is missing for target {pdb_id}")
        reps = representatives.loc[
            representatives["cluster_v2_id"].eq(cluster_id)
        ].sort_values(["tool", "raw_row_id"], kind="mergesort")
        if len(reps) != int(candidate["tool_support_count"]):
            raise ValueError(f"formal representative count mismatch for {cluster_id}")
        center_values = [candidate[key] for key in (
            "medoid_center_x", "medoid_center_y", "medoid_center_z"
        )]
        center = (
            tuple(float(value) for value in center_values)
            if all(pd.notna(value) for value in center_values)
            else None
        )
        envelope = parse_residues(candidate["envelope_residues"])
        interface = score_oligomer_relevance(
            envelope, center, interfaces[pdb_id], config=config
        )
        modules = {
            "C_cons": _consensus_score(reps, config),
            "G_geo": _top_mean(reps["pocket_geometry_score"], config.geometry_top_n),
            "P_lig": _top_mean(
                reps["pocket_ligandability_score"], config.ligandability_top_n
            ),
            "Q_evidence": _evidence_quality(reps, config),
            "O_rel_formal": float(interface["O_rel_formal"]),
        }
        row = {
            key: value
            for key, value in candidate.items()
            if key != "boundary_sensitive"
        }
        row.update(modules)
        row.update(
            {
                "interface_fraction": interface["interface_fraction"],
                "interface_recall": interface["interface_recall"],
                "cluster_interface_residue_count": interface[
                    "cluster_interface_residue_count"
                ],
                "cluster_chain_count": interface["cluster_chain_count"],
                "cluster_chain_entropy": interface["cluster_chain_entropy"],
                "distance_to_interface_A": interface["distance_to_interface_A"],
                "interface_distance_score": interface["interface_distance_score"],
                "boundary_sensitive": bool(candidate["boundary_sensitive"]),
            }
        )
        row["OIPS-P_static"] = weighted_module_score(row, config.module_weights)
        master_rows.append(row)
    master = pd.DataFrame(master_rows, columns=list(MASTER_COLUMNS))
    master = master.sort_values(
        ["pdb_id", "cluster_v2_id"], kind="mergesort"
    ).reset_index(drop=True)
    master.attrs["module_weights"] = dict(config.module_weights)
    master.attrs["ranking_direction"] = config.ranking_direction
    rankings = rank_within_target(
        master,
        score_column="OIPS-P_static",
        rank_method=config.ranking_tie_method,
    )
    difference = (
        rankings["OIPS-P_static_recomputed"] - rankings["OIPS-P_static"]
    ).abs()
    maximum = float(difference.max()) if len(difference) else 0.0
    return ScoringResult(master=master, rankings=rankings, maximum_recomputation_difference=maximum)


def rank_within_target(
    master: pd.DataFrame,
    *,
    score_column: str = "OIPS-P_static",
    rank_method: str = "average",
) -> pd.DataFrame:
    """Recompute static OIPS and apply exact-value average rank within targets."""
    if score_column not in master or not {"pdb_id", "cluster_v2_id"}.issubset(master):
        raise KeyError("ranking requires pdb_id, cluster_v2_id, and the score column")
    output = master.copy()
    module_names = (
        "C_cons", "G_geo", "P_lig", "O_rel_formal", "Q_evidence"
    )
    if set(module_names).issubset(output.columns):
        weights = output.attrs.get("module_weights")
        if not isinstance(weights, Mapping):
            weights = ScoringConfig().module_weights
        recomputed = [
            weighted_module_score(row, weights)
            for row in output.loc[:, list(module_names)].to_dict("records")
        ]
        output["OIPS-P_static_recomputed"] = recomputed
        difference = (
            output["OIPS-P_static_recomputed"] - output[score_column]
        ).abs()
        if len(difference) and float(difference.max()) > 1e-10:
            raise AssertionError(
                "stored OIPS-P_static differs from the configured formula: "
                f"max difference={float(difference.max())}"
            )
    else:
        output["OIPS-P_static_recomputed"] = pd.to_numeric(
            output[score_column], errors="coerce"
        )
    direction = output.attrs.get("ranking_direction", "descending")
    if direction not in {"descending", "ascending"}:
        raise ValueError(f"unsupported ranking direction: {direction}")
    output["Within_PDB_rank"] = output.groupby("pdb_id")[
        "OIPS-P_static_recomputed"
    ].rank(method=rank_method, ascending=direction == "ascending")
    tie_size = output.groupby(
        ["pdb_id", "OIPS-P_static_recomputed"], dropna=False
    )["cluster_v2_id"].transform("size")
    output["tie_flag"] = tie_size > 1
    output["tie_size"] = tie_size.astype(int)
    return output.sort_values(
        ["pdb_id", "Within_PDB_rank", "cluster_v2_id"], kind="mergesort"
    ).reset_index(drop=True)
