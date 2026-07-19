"""Contract-driven manuscript figures built only from validated analysis tables."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import yaml

from . import io

if TYPE_CHECKING:
    from matplotlib.figure import Figure


SOURCE_COLUMNS = (
    "panel", "record_id", "series", "group", "x", "y", "x_label",
    "y_label", "lower", "upper", "annotation",
)
MODULES = ("C_cons", "G_geo", "P_lig", "O_rel_formal", "Q_evidence")
MODULE_LABELS = {
    "C_cons": "Consensus", "G_geo": "Geometry", "P_lig": "Ligandability",
    "O_rel_formal": "Oligomer relevance", "Q_evidence": "Evidence quality",
}


@dataclass(frozen=True)
class PanelSpec:
    panel_id: str
    title: str
    x_label: str
    y_label: str


@dataclass(frozen=True)
class FigureSpec:
    figure_id: str
    stem: str
    size_inches: tuple[float, float]
    panels: Mapping[str, PanelSpec]


@dataclass(frozen=True)
class FigureContract:
    schema_version: int
    formats: tuple[str, ...]
    png_dpi: int
    editable_svg_text: bool
    palette: Mapping[str, str]
    category_order: tuple[str, ...]
    category_colors: Mapping[str, str]
    figures: Mapping[str, FigureSpec]
    case_titles: Mapping[str, str]


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"figure contract section must be a mapping: {label}")
    return value


def load_figure_contract(path: str | Path) -> FigureContract:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"figure contract does not exist: {source}")
    data = yaml.safe_load(source.read_text(encoding=io.ENCODING))
    root = _mapping(data, "root")
    if int(root.get("schema_version", 0)) != 1:
        raise ValueError("unsupported figure contract schema version")
    output = _mapping(root.get("output"), "output")
    formats = tuple(str(value).lower() for value in output.get("formats", ()))
    if formats != ("svg", "png"):
        raise ValueError("figure output formats must be [svg, png]")
    dpi = int(output.get("png_dpi", 0))
    if dpi != 600 or output.get("editable_svg_text") is not True:
        raise ValueError("figure contract requires 600-dpi PNG and editable SVG text")
    figure_data = _mapping(root.get("figures"), "figures")
    figures: dict[str, FigureSpec] = {}
    for figure_id in ("figure1", "figure2", "figure3", "figure4"):
        values = _mapping(figure_data.get(figure_id), f"figures.{figure_id}")
        size = values.get("size_inches")
        if not isinstance(size, Sequence) or isinstance(size, (str, bytes)) or len(size) != 2:
            raise ValueError(f"invalid figure size: {figure_id}")
        panels_data = _mapping(values.get("panels"), f"figures.{figure_id}.panels")
        panels: dict[str, PanelSpec] = {}
        for panel_id, panel_value in panels_data.items():
            panel = _mapping(panel_value, f"{figure_id}.{panel_id}")
            panels[str(panel_id)] = PanelSpec(
                str(panel_id), str(panel["title"]), str(panel["x_label"]), str(panel["y_label"])
            )
        expected = 6 if figure_id == "figure4" else 4
        if len(panels) != expected:
            raise ValueError(f"{figure_id} must declare exactly {expected} panels")
        figures[figure_id] = FigureSpec(
            figure_id, str(values["stem"]), (float(size[0]), float(size[1])), panels
        )
    stems = [spec.stem for spec in figures.values()]
    if len(set(stems)) != 4:
        raise ValueError("figure stems must be unique")
    return FigureContract(
        1, formats, dpi, True,
        {str(key): str(value) for key, value in _mapping(root.get("palette"), "palette").items()},
        tuple(str(value) for value in root.get("category_order", ())),
        {str(key): str(value) for key, value in _mapping(root.get("category_colors"), "category_colors").items()},
        figures,
        {str(key).upper(): str(value) for key, value in _mapping(root.get("case_titles"), "case_titles").items()},
    )


def _style():
    import matplotlib as mpl
    mpl.use("Agg", force=True)
    import matplotlib.pyplot as plt

    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans", "sans-serif"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "legend.frameon": False,
        "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.facecolor": "white",
    })
    return plt


def _source(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=SOURCE_COLUMNS).sort_values(
        ["panel", "record_id", "series"], kind="mergesort"
    ).reset_index(drop=True)


def _row(panel: PanelSpec, record_id: object, series: str, **values: object) -> dict[str, object]:
    return {
        "panel": panel.panel_id, "record_id": str(record_id), "series": series,
        "group": values.get("group", ""), "x": values.get("x", math.nan),
        "y": values.get("y", math.nan), "x_label": panel.x_label,
        "y_label": panel.y_label, "lower": values.get("lower", math.nan),
        "upper": values.get("upper", math.nan), "annotation": values.get("annotation", ""),
    }


def _panel(axis, panel: PanelSpec, label: str) -> None:
    axis.set_title(panel.title, loc="left", fontweight="bold")
    axis.set_xlabel(panel.x_label)
    axis.set_ylabel(panel.y_label)
    axis.text(-0.12, 1.06, label, transform=axis.transAxes, fontsize=9, fontweight="bold", va="top")


def _gid(artist, panel: str, series: str) -> None:
    artist.set_gid(f"source:{panel}:{series}")


def _bar_gids(container, panel: str, series: str) -> None:
    for patch in container.patches:
        _gid(patch, panel, series)


def _strict_bool(value: object, label: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str) and value.strip().casefold() in {"true", "false"}:
        return value.strip().casefold() == "true"
    raise ValueError(f"{label} must be true or false")


def _endpoint_fraction(frame: pd.DataFrame, value: str, evaluable: str) -> float:
    if evaluable not in frame or value not in frame:
        raise ValueError(f"target summary is missing endpoint fields: {evaluable}, {value}")
    mask = pd.Series(
        [_strict_bool(item, evaluable) for item in frame[evaluable]], index=frame.index,
    )
    selected = frame.loc[mask, value]
    if not len(selected):
        return math.nan
    return float(np.mean([_strict_bool(item, value) for item in selected]))


def _rank_threshold_fraction(
    frame: pd.DataFrame, rank: str, evaluable: str, threshold: float
) -> float:
    if evaluable not in frame or rank not in frame:
        raise ValueError(f"target summary is missing endpoint fields: {evaluable}, {rank}")
    mask = pd.Series(
        [_strict_bool(item, evaluable) for item in frame[evaluable]], index=frame.index,
    )
    if not mask.any():
        return math.nan
    return float(pd.to_numeric(frame.loc[mask, rank], errors="raise").le(threshold).mean())


def build_figure1(tables: Mapping[str, pd.DataFrame], contract: FigureContract) -> tuple[Figure, pd.DataFrame]:
    spec = contract.figures["figure1"]
    panels = list(spec.panels.values())
    target = tables["target_summary"].copy()
    bootstrap = tables["bootstrap"].copy()
    rows: list[dict[str, object]] = []
    landscape = target.sort_values(["formal_cluster_v2_count", "pdb_id"], kind="mergesort").reset_index(drop=True)
    for index, item in landscape.iterrows():
        rows.append(_row(panels[0], item["pdb_id"], "formal_cluster_v2_count", group=item["Pocket_category"], x=item["formal_cluster_v2_count"], y=index))
    ranks = target.sort_values(["reference_first_rank", "pdb_id"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    for index, item in ranks.iterrows():
        rows.append(_row(panels[1], item["pdb_id"], "First R_auto", group=item["Pocket_category"], x=item["reference_first_rank"], y=index))
        rows.append(_row(panels[1], item["pdb_id"], "First R_auto/A_auto", group=item["Pocket_category"], x=item["first_supported_rank"], y=index))
    endpoint_data = (
        ("Top-1", "Reference-associated", _endpoint_fraction(target, "reference_top1", "reference_evaluable")),
        ("Top-3", "Reference-associated", _endpoint_fraction(target, "reference_top3", "reference_evaluable")),
        ("Top-5", "Reference-associated", _endpoint_fraction(target, "reference_top5", "reference_evaluable")),
        ("Top-1", "First supported", _endpoint_fraction(target, "first_supported_top1", "first_supported_evaluable")),
        ("Top-3", "First supported", _endpoint_fraction(target, "first_supported_top3", "first_supported_evaluable")),
        ("Top-5", "First supported", _rank_threshold_fraction(target, "first_supported_rank", "first_supported_evaluable", 5)),
    )
    for threshold, series, value in endpoint_data:
        rows.append(_row(panels[2], threshold, series, group=threshold, x=("Top-1", "Top-3", "Top-5").index(threshold), y=value, annotation=f"{value:.2f}"))
    metrics = ("Reference_Top1", "Reference_Top3", "Reference_MRR", "First_supported_Top1", "First_supported_Top3", "First_supported_MRR")
    lookup = bootstrap.set_index(["bootstrap_method", "metric"])
    for index, metric in enumerate(metrics):
        for method in ("target_resampling", "family_clustered"):
            item = lookup.loc[(method, metric)]
            rows.append(_row(panels[3], metric, method, group=metric, x=item["point_estimate"], y=index, lower=item["CI_2.5_percent"], upper=item["CI_97.5_percent"]))
    source = _source(rows)

    plt = _style()
    fig, axes = plt.subplots(2, 2, figsize=spec.size_inches, constrained_layout=True)
    axes = axes.ravel()
    a = source.loc[source["panel"].eq("a")]
    colors = [contract.category_colors.get(value, contract.palette["gray"]) for value in a["group"]]
    bars = axes[0].barh(a["y"], a["x"], color=colors, height=0.72)
    _bar_gids(bars, "a", "formal_cluster_v2_count")
    axes[0].set_yticks(a["y"], a["record_id"])
    _panel(axes[0], panels[0], "a")
    b = source.loc[source["panel"].eq("b")]
    for series, color, marker in (("First R_auto", contract.palette["coral"], "o"), ("First R_auto/A_auto", contract.palette["teal"], "s")):
        view = b.loc[b["series"].eq(series)]
        artist = axes[1].scatter(view["x"], view["y"], s=22, color=color, marker=marker, label=series, zorder=3)
        _gid(artist, "b", series)
    axes[1].set_yticks(sorted(b["y"].unique()), b.drop_duplicates("y").sort_values("y")["record_id"])
    axes[1].legend(fontsize=6.2)
    _panel(axes[1], panels[1], "b")
    c = source.loc[source["panel"].eq("c")]
    width = 0.34
    for series, offset, color in (("Reference-associated", -width / 2, contract.palette["coral"]), ("First supported", width / 2, contract.palette["teal"])):
        view = c.loc[c["series"].eq(series)].sort_values("x")
        bars = axes[2].bar(view["x"] + offset, view["y"], width, color=color, label=series)
        _bar_gids(bars, "c", series)
    axes[2].set_xticks([0, 1, 2], ["Top-1", "Top-3", "Top-5"])
    axes[2].set_ylim(0, 1.12)
    axes[2].legend(fontsize=6.2)
    _panel(axes[2], panels[2], "c")
    d = source.loc[source["panel"].eq("d")]
    for method, offset, color, marker in (("target_resampling", -0.1, contract.palette["blue"], "o"), ("family_clustered", 0.1, contract.palette["navy"], "s")):
        view = d.loc[d["series"].eq(method)].sort_values("y")
        container = axes[3].errorbar(view["x"], view["y"] + offset, xerr=np.vstack([view["x"] - view["lower"], view["upper"] - view["x"]]), fmt=marker, color=color, markersize=3.5, capsize=2, lw=0.8, label=method.replace("_", " "))
        _gid(container.lines[0], "d", method)
    axes[3].set_yticks(range(len(metrics)), metrics)
    axes[3].set_xlim(0, 1.04)
    axes[3].legend(fontsize=6.2)
    _panel(axes[3], panels[3], "d")
    return fig, source


def build_figure2(tables: Mapping[str, pd.DataFrame], contract: FigureContract) -> tuple[Figure, pd.DataFrame]:
    spec = contract.figures["figure2"]
    panels = list(spec.panels.values())
    qc = tables["topk_qc"]
    ablation = tables["orel_targets"].copy()
    rows: list[dict[str, object]] = []
    qc_order = ("QC_pass", "QC_boundary_sensitive", "QC_possible_split", "QC_possible_overmerge", "QC_insufficient_evidence")
    counts = qc["QC_status"].value_counts()
    for index, status in enumerate(qc_order):
        rows.append(_row(panels[0], status, "QC_count", group=status, x=int(counts.get(status, 0)), y=index))
    for item in ablation.to_dict("records"):
        changed = str(item["top_cluster_identity_changed"]).casefold() == "true"
        rows.append(_row(panels[1], item["pdb_id"], "reference_rank", group=str(changed).lower(), x=item["Full_reference_rank"], y=item["Without_O_rel_reference_rank"], annotation=item["pdb_id"] if abs(float(item["reference_rank_change_without_minus_full"])) >= 2 else ""))
        rows.append(_row(panels[2], item["pdb_id"], "first_supported_rank", group=str(changed).lower(), x=item["Full_first_supported_rank"], y=item["Without_O_rel_first_supported_rank"], annotation=item["pdb_id"] if abs(float(item["first_supported_rank_change_without_minus_full"])) >= 2 else ""))
    for index, item in ablation.sort_values("rank_Spearman_rho", kind="mergesort").reset_index(drop=True).iterrows():
        rows.append(_row(panels[3], item["pdb_id"], "rank_Spearman_rho", group=str(item["top_cluster_identity_changed"]).lower(), x=item["rank_Spearman_rho"], y=index))
    source = _source(rows)
    plt = _style()
    fig, axes = plt.subplots(2, 2, figsize=spec.size_inches, constrained_layout=True)
    axes = axes.ravel()
    a = source.loc[source["panel"].eq("a")].sort_values("y")
    colors = [contract.palette[name] for name in ("teal", "blue", "gold", "coral", "gray")]
    bars = axes[0].barh(a["y"], a["x"], color=colors, height=0.7)
    _bar_gids(bars, "a", "QC_count")
    axes[0].set_yticks(a["y"], a["record_id"].str.replace("QC_", "", regex=False))
    _panel(axes[0], panels[0], "a")
    for axis, panel_id, series, panel in ((axes[1], "b", "reference_rank", panels[1]), (axes[2], "c", "first_supported_rank", panels[2])):
        view = source.loc[source["panel"].eq(panel_id)]
        colors = np.where(view["group"].eq("true"), contract.palette["coral"], contract.palette["navy"] if panel_id == "b" else contract.palette["teal"])
        artist = axis.scatter(view["x"], view["y"], c=colors, s=24)
        _gid(artist, panel_id, series)
        finite = pd.concat([pd.to_numeric(view["x"], errors="coerce"), pd.to_numeric(view["y"], errors="coerce")]).dropna()
        limit = max(4.0, float(finite.max()) if len(finite) else 4.0) + 0.5
        axis.plot([0.5, limit], [0.5, limit], color=contract.palette["gray"], lw=0.8, ls="--")
        axis.set_xlim(0.5, limit)
        axis.set_ylim(0.5, limit)
        for item in view.loc[view["annotation"].astype(str).ne("")].itertuples():
            axis.text(item.x + 0.12, item.y + 0.12, item.annotation, fontsize=6)
        _panel(axis, panel, panel_id)
    d = source.loc[source["panel"].eq("d")].sort_values("y")
    artist = axes[3].scatter(d["x"], d["y"], c=np.where(d["group"].eq("true"), contract.palette["coral"], contract.palette["navy"]), s=20)
    _gid(artist, "d", "rank_Spearman_rho")
    axes[3].axvline(pd.to_numeric(d["x"]).mean(), color=contract.palette["gray"], ls="--", lw=0.8)
    axes[3].set_yticks(d["y"], d["record_id"])
    axes[3].set_xlim(0, 1.03)
    _panel(axes[3], panels[3], "d")
    return fig, source


def build_figure3(tables: Mapping[str, pd.DataFrame], contract: FigureContract) -> tuple[Figure, pd.DataFrame]:
    spec = contract.figures["figure3"]
    panels = list(spec.panels.values())
    labels = tables["evidence_labels"].loc[tables["evidence_labels"]["Within_PDB_rank"].le(3)]
    md = tables["md_mapping"]
    redocking = tables["redocking_mapping"].sort_values("Raw_ligand_RMSD_A", kind="mergesort").reset_index(drop=True)
    category = tables["category_metrics"]
    rows: list[dict[str, object]] = []
    label_order = ("R_auto", "A_auto", "U_auto", "X_auto")
    for rank in (1, 2, 3):
        for label in label_order:
            count = int((labels["Within_PDB_rank"].eq(rank) & labels["automated_evidence_label"].eq(label)).sum())
            rows.append(_row(panels[0], f"rank_{rank}", label, group=label, x=rank, y=count))
    call_order = ("concordant", "partially_concordant", "boundary_shift", "static_dynamic_conflict", "insufficient_MD_evidence", "apo_only_context")
    real = md.loc[~md["Concordance_call"].eq("MD_not_available")]
    missing_targets = md.loc[md["Concordance_call"].eq("MD_not_available"), "pdb_id"].astype(str).nunique()
    for index, call in enumerate(call_order):
        rows.append(_row(panels[1], call, "MD_runs", group=call, x=int(real["Concordance_call"].eq(call).sum()), y=index, annotation=f"{missing_targets} targets without MD" if index == 0 else ""))
    for index, item in redocking.iterrows():
        rows.append(_row(panels[2], item["pdb_id"], "Raw_ligand_RMSD_A", group="<=2" if item["Raw_ligand_RMSD_A"] <= 2 else "<=3" if item["Raw_ligand_RMSD_A"] <= 3 else ">3", x=item["Raw_ligand_RMSD_A"], y=index))
    category_lookup = category.set_index("group")
    category_names = [value for value in contract.category_order if value in category_lookup.index]
    category_names.extend(sorted(set(category_lookup.index) - set(category_names)))
    for index, name in enumerate(category_names):
        for series, column in (("Reference Top-1", "Reference_Top1"), ("First-supported Top-1", "First_supported_Top1")):
            rows.append(_row(panels[3], name, series, group=name, x=category_lookup.loc[name, column], y=index))
    source = _source(rows)
    plt = _style()
    fig, axes = plt.subplots(2, 2, figsize=spec.size_inches, constrained_layout=True)
    axes = axes.ravel()
    a = source.loc[source["panel"].eq("a")]
    bottom = np.zeros(3)
    for label, color in zip(label_order, (contract.palette["coral"], contract.palette["teal"], contract.palette["gold"], contract.palette["gray"])):
        view = a.loc[a["series"].eq(label)].sort_values("x")
        bars = axes[0].bar(view["x"], view["y"], bottom=bottom, color=color, width=0.68, label=label)
        _bar_gids(bars, "a", label)
        bottom += view["y"].to_numpy(float)
    axes[0].set_xticks([1, 2, 3], ["Rank 1", "Rank 2", "Rank 3"])
    axes[0].legend(ncol=2, fontsize=6.2)
    _panel(axes[0], panels[0], "a")
    b = source.loc[source["panel"].eq("b")].sort_values("y")
    bars = axes[1].barh(b["y"], b["x"], color=[contract.palette[name] for name in ("teal", "blue", "gold", "coral", "light", "gray")], height=0.7)
    _bar_gids(bars, "b", "MD_runs")
    axes[1].set_yticks(b["y"], b["record_id"])
    annotation = b["annotation"].loc[b["annotation"].astype(str).ne("")]
    if len(annotation):
        axes[1].text(0.98, 0.02, annotation.iloc[0], transform=axes[1].transAxes, ha="right", va="bottom", fontsize=6.2, color=contract.palette["gray"])
    _panel(axes[1], panels[1], "b")
    c = source.loc[source["panel"].eq("c")].sort_values("y")
    colors = np.where(c["group"].eq("<=2"), contract.palette["teal"], np.where(c["group"].eq("<=3"), contract.palette["gold"], contract.palette["coral"]))
    axes[2].hlines(c["y"], 0, c["x"], color="#D5DDE1", lw=0.8)
    artist = axes[2].scatter(c["x"], c["y"], c=colors, s=20)
    _gid(artist, "c", "Raw_ligand_RMSD_A")
    axes[2].axvline(2, color=contract.palette["gray"], ls="--", lw=0.8)
    axes[2].axvline(3, color=contract.palette["gray"], ls=":", lw=0.8)
    axes[2].set_yticks(c["y"], c["record_id"])
    _panel(axes[2], panels[2], "c")
    d = source.loc[source["panel"].eq("d")]
    height = 0.34
    for series, offset, color in (("Reference Top-1", -height / 2, contract.palette["coral"]), ("First-supported Top-1", height / 2, contract.palette["teal"])):
        view = d.loc[d["series"].eq(series)].sort_values("y")
        bars = axes[3].barh(view["y"] + offset, view["x"], height, color=color, label=series)
        _bar_gids(bars, "d", series)
    labels_y = d.drop_duplicates("y").sort_values("y")
    axes[3].set_yticks(labels_y["y"], labels_y["record_id"])
    axes[3].set_xlim(0, 1.05)
    axes[3].legend(fontsize=6.2)
    _panel(axes[3], panels[3], "d")
    return fig, source


def build_figure4(
    tables: Mapping[str, pd.DataFrame], contract: FigureContract, case_ids: Sequence[str]
) -> tuple[Figure, pd.DataFrame]:
    spec = contract.figures["figure4"]
    panels = list(spec.panels.values())
    cases = tables["representative_cases"].copy()
    cases["pdb_id"] = cases["pdb_id"].astype(str).str.upper()
    case_lookup = cases.set_index("pdb_id")
    master = tables["master"].copy()
    master["pdb_id"] = master["pdb_id"].astype(str).str.upper()
    requested = tuple(str(value).upper() for value in case_ids)
    if len(panels) != 2 * len(requested) or not set(requested).issubset(case_lookup.index):
        raise ValueError("configured representative cases do not match the figure contract")
    rows: list[dict[str, object]] = []
    for index, pdb_id in enumerate(requested):
        case = case_lookup.loc[pdb_id]
        module_panel, rank_panel = panels[index * 2:index * 2 + 2]
        top = master.loc[master["pdb_id"].eq(pdb_id) & master["cluster_v2_id"].eq(case["static_top_cluster_v2_id"])]
        if len(top) != 1 or pdb_id not in contract.case_titles:
            raise ValueError(f"missing representative figure inputs: {pdb_id}")
        top_row = top.iloc[0]
        for module_index, module in enumerate(MODULES):
            value = pd.to_numeric(pd.Series([top_row[module]]), errors="coerce").iloc[0]
            rows.append(_row(module_panel, pdb_id, module, group=MODULE_LABELS[module], x=value, y=module_index, annotation="NA" if pd.isna(value) else ""))
        rank_annotation = f"MD: {case['MD_calls']}; Redocking RMSD: {float(case['redocking_RMSD_A']):.3f} A"
        rows.append(_row(rank_panel, f"{pdb_id}_full", "Reference rank", group=pdb_id, x=0, y=case["reference_rank"], annotation=rank_annotation))
        rows.append(_row(rank_panel, f"{pdb_id}_without", "Reference rank", group=pdb_id, x=1, y=case["without_O_rel_reference_rank"]))
        rows.append(_row(rank_panel, f"{pdb_id}_static", "Static top", group=pdb_id, x=0, y=1.0))
    source = _source(rows)
    plt = _style()
    fig, axes = plt.subplots(len(requested), 2, figsize=spec.size_inches, constrained_layout=True, gridspec_kw={"width_ratios": [1.25, 1]})
    axes = np.asarray(axes).reshape(-1)
    for index, pdb_id in enumerate(requested):
        module_panel, rank_panel = panels[index * 2:index * 2 + 2]
        module_axis, rank_axis = axes[index * 2:index * 2 + 2]
        module_data = source.loc[source["panel"].eq(module_panel.panel_id)].sort_values("y")
        colors = [contract.palette[name] for name in ("navy", "blue", "gold", "teal", "gray")]
        bars = module_axis.barh(module_data["y"], module_data["x"], color=colors, height=0.65)
        for patch, module in zip(bars.patches, module_data["series"]):
            _gid(patch, module_panel.panel_id, module)
        module_axis.set_yticks(module_data["y"], module_data["group"])
        module_axis.invert_yaxis()
        module_axis.set_xlim(0, 105)
        for item in module_data.loc[module_data["annotation"].eq("NA")].itertuples():
            module_axis.text(2, item.y, "NA", va="center", fontsize=6.5, color=contract.palette["gray"])
        title = contract.case_titles[pdb_id]
        module_axis.set_title(f"{pdb_id} | {title}", loc="left", fontweight="bold")
        module_axis.set_xlabel(module_panel.x_label)
        module_axis.set_ylabel(module_panel.y_label)
        module_axis.text(-0.12, 1.06, module_panel.panel_id, transform=module_axis.transAxes, fontsize=9, fontweight="bold", va="top")
        rank_data = source.loc[source["panel"].eq(rank_panel.panel_id)]
        reference = rank_data.loc[rank_data["series"].eq("Reference rank")].sort_values("x")
        line = rank_axis.plot(reference["x"], reference["y"], color=contract.palette["coral"], marker="o", lw=1.1)[0]
        _gid(line, rank_panel.panel_id, "Reference rank")
        static = rank_data.loc[rank_data["series"].eq("Static top")]
        artist = rank_axis.scatter(static["x"], static["y"], marker="s", color=contract.palette["teal"], s=28)
        _gid(artist, rank_panel.panel_id, "Static top")
        rank_axis.set_xticks([0, 1], ["Full", "Without O_rel"])
        finite = pd.to_numeric(reference["y"], errors="coerce").dropna()
        rank_axis.set_ylim(max(4.0, float(finite.max()) if len(finite) else 4.0) + 1, 0.5)
        annotation = reference["annotation"].loc[reference["annotation"].astype(str).ne("")]
        if len(annotation):
            rank_axis.text(0.04, 0.05, annotation.iloc[0], transform=rank_axis.transAxes, fontsize=6.2)
        _panel(rank_axis, rank_panel, rank_panel.panel_id)
    return fig, source


def save_figure(
    figure: Figure,
    source_data: pd.DataFrame,
    output_dir: str | Path,
    spec: FigureSpec,
    contract: FigureContract,
) -> tuple[Path, ...]:
    """Save one figure and its complete, stable source-data table."""
    if source_data.columns.tolist() != list(SOURCE_COLUMNS):
        raise ValueError("figure source data has an invalid schema")
    plt = _style()
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for fmt in contract.formats:
        path = destination / f"{spec.stem}.{fmt}"
        if fmt == "png":
            figure.savefig(path, dpi=contract.png_dpi, bbox_inches="tight")
        else:
            figure.savefig(path, bbox_inches="tight")
        written.append(path)
    source_path = destination / f"{spec.stem}_source_data.csv"
    io.write_stable_csv(source_data, source_path, columns=SOURCE_COLUMNS, sort_by=("panel", "record_id", "series"))
    written.append(source_path)
    plt.close(figure)
    return tuple(written)


def render_all_figures(
    tables: Mapping[str, pd.DataFrame],
    contract: FigureContract,
    case_ids: Sequence[str],
    output_dir: str | Path,
) -> tuple[Path, ...]:
    builders = (
        (build_figure1, (tables, contract)),
        (build_figure2, (tables, contract)),
        (build_figure3, (tables, contract)),
        (build_figure4, (tables, contract, case_ids)),
    )
    written: list[Path] = []
    for builder, arguments in builders:
        figure, source = builder(*arguments)
        figure_id = builder.__name__.replace("build_", "")
        written.extend(save_figure(figure, source, output_dir, contract.figures[figure_id], contract))
    return tuple(written)
