"""Native single-tool complete-case comparison."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

import pandas as pd


SINGLE_TOOL_TARGET_COLUMNS = (
    "pdb_id", "method", "output_status", "mappable_region_count",
    "native_top_cluster_v2_id", "first_reference_associated_rank", "complete_case",
)
SINGLE_TOOL_METRIC_COLUMNS = (
    "method", "N", "Top1_N", "Top1", "Top3_N", "Top3", "Top5_N", "Top5", "MRR",
)


@dataclass(frozen=True)
class SingleToolResult:
    target_ranks: pd.DataFrame
    complete_case_metrics: pd.DataFrame


def _require(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def single_tool_prioritization(
    features: pd.DataFrame,
    record_mapping: pd.DataFrame,
    rankings: pd.DataFrame,
    labels: pd.DataFrame,
) -> SingleToolResult:
    """Compare native regional priorities on the five-tool complete cases."""
    tools = ("CavityPlus", "DoGSiteScorer", "DoGSite3", "CASTpFold", "SiteMap")
    _require(features, ("row_id", "pdb_id", "tool", "display_order", "sitemap_rank"), "features")
    _require(record_mapping, (
        "row_id", "pdb_id", "tool", "same_tool_unit_id", "representative_for_tool_unit",
        "cluster_v2_id", "mapping_status",
    ), "record mapping")
    _require(rankings, ("pdb_id", "cluster_v2_id", "Within_PDB_rank"), "rankings")
    _require(labels, ("pdb_id", "cluster_v2_id", "automated_evidence_label"), "labels")
    feature = features.copy(deep=True)
    mapping = record_mapping.copy(deep=True)
    marked = labels.copy(deep=True)
    ranked = rankings.copy(deep=True)
    for frame in (feature, mapping, marked, ranked):
        frame["pdb_id"] = frame["pdb_id"].astype(str).str.upper()
    targets = sorted(ranked["pdb_id"].unique())
    if mapping["row_id"].duplicated().any() or feature["row_id"].duplicated().any():
        raise ValueError("single-tool analysis requires unique row_id values")
    representative_text = mapping["representative_for_tool_unit"].fillna("").astype(str).str.strip().str.casefold()
    invalid_representative = ~representative_text.isin({"", "true", "false"})
    if invalid_representative.any():
        raise ValueError("representative_for_tool_unit must be true, false, or missing")
    representatives = mapping.loc[representative_text.eq("true")].copy()
    representatives = representatives.merge(
        feature[["row_id", "display_order", "sitemap_rank"]], on="row_id", how="left", validate="one_to_one",
    )
    representatives["native_order"] = pd.to_numeric(
        representatives["display_order"], errors="coerce",
    ).fillna(pd.to_numeric(representatives["sitemap_rank"], errors="coerce")).fillna(
        pd.to_numeric(representatives["row_id"], errors="raise")
    )
    mapped = representatives.loc[representatives["cluster_v2_id"].notna()].copy()
    mapped = mapped.sort_values(
        ["pdb_id", "tool", "native_order", "row_id"], kind="mergesort",
    ).drop_duplicates(["pdb_id", "tool", "cluster_v2_id"], keep="first")
    mapped["native_region_rank"] = mapped.groupby(["pdb_id", "tool"]).cumcount() + 1
    label_lookup = marked.set_index(["pdb_id", "cluster_v2_id"])["automated_evidence_label"]
    mapped_keys = pd.MultiIndex.from_frame(mapped[["pdb_id", "cluster_v2_id"]])
    if not set(mapped_keys).issubset(set(label_lookup.index)):
        raise ValueError("single-tool mappings contain unknown cluster-v2 IDs")
    mapped["label"] = label_lookup.reindex(mapped_keys).to_numpy()
    coverage = mapped.groupby(["pdb_id", "tool"]).size().unstack(fill_value=0).reindex(
        index=targets, columns=tools, fill_value=0,
    )
    complete_targets = set(coverage.index[(coverage > 0).all(axis=1)])
    raw_coverage = feature.groupby(["pdb_id", "tool"]).size()

    rows: list[dict[str, object]] = []
    for pdb_id in targets:
        for tool in tools:
            group = mapped.loc[mapped["pdb_id"].eq(pdb_id) & mapped["tool"].eq(tool)]
            raw_count = int(raw_coverage.get((pdb_id, tool), 0))
            refs = group.loc[group["label"].eq("R_auto"), "native_region_rank"]
            if raw_count == 0:
                status = "unavailable"
            elif group.empty:
                status = "unmappable"
            elif refs.empty:
                status = "no_hit"
            else:
                status = "reference_hit"
            ordered = group.sort_values(["native_region_rank", "cluster_v2_id"], kind="mergesort")
            rows.append({
                "pdb_id": pdb_id,
                "method": tool,
                "output_status": status,
                "mappable_region_count": len(group),
                "native_top_cluster_v2_id": "" if group.empty else ordered.iloc[0]["cluster_v2_id"],
                "first_reference_associated_rank": math.nan if refs.empty else float(refs.min()),
                "complete_case": pdb_id in complete_targets,
            })
        oips = ranked.loc[ranked["pdb_id"].eq(pdb_id)].sort_values(
            ["Within_PDB_rank", "cluster_v2_id"], kind="mergesort",
        )
        oips_keys = pd.MultiIndex.from_frame(oips[["pdb_id", "cluster_v2_id"]])
        oips_labels = label_lookup.reindex(oips_keys).to_numpy()
        reference_ranks = pd.to_numeric(
            oips.loc[pd.Series(oips_labels, index=oips.index).eq("R_auto"), "Within_PDB_rank"], errors="raise",
        )
        rows.append({
            "pdb_id": pdb_id,
            "method": "OIPS-P",
            "output_status": "reference_hit" if len(reference_ranks) else "no_hit",
            "mappable_region_count": len(oips),
            "native_top_cluster_v2_id": oips.iloc[0]["cluster_v2_id"],
            "first_reference_associated_rank": math.nan if not len(reference_ranks) else float(reference_ranks.min()),
            "complete_case": pdb_id in complete_targets,
        })
    target_ranks = pd.DataFrame(rows, columns=SINGLE_TOOL_TARGET_COLUMNS)
    metric_rows = []
    subset = target_ranks.loc[target_ranks["complete_case"]]
    for method in (*tools, "OIPS-P"):
        group = subset.loc[subset["method"].eq(method)]
        ranks = pd.to_numeric(group["first_reference_associated_rank"], errors="coerce")
        n = len(group)
        metric_rows.append({
            "method": method,
            "N": n,
            "Top1_N": int(ranks.le(1).sum()),
            "Top1": float(ranks.le(1).sum() / n),
            "Top3_N": int(ranks.le(3).sum()),
            "Top3": float(ranks.le(3).sum() / n),
            "Top5_N": int(ranks.le(5).sum()),
            "Top5": float(ranks.le(5).sum() / n),
            "MRR": float((1.0 / ranks).fillna(0.0).mean()),
        })
    return SingleToolResult(
        target_ranks=target_ranks,
        complete_case_metrics=pd.DataFrame(metric_rows, columns=SINGLE_TOOL_METRIC_COLUMNS),
    )
