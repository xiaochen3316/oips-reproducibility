"""Canonical reference-release inventory and read-only comparison helpers."""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from . import validation


FIGURE_STEMS = (
    "repository_summary_figure_1_candidate_landscape",
    "repository_summary_figure_2_qc_and_orel_ablation",
    "repository_summary_figure_3_posthoc_evidence",
    "repository_summary_figure_4_representative_cases",
)
REPORT_NAMES = {
    "rebuild_report.md", "analysis_report.md", "numeric_crosscheck.md",
    "validation_report.md",
}
REFERENCE_DIRECTORIES = {
    "clustering", "static", "analysis", "figure_source_data", "reports",
}


def _directory_files(path: Path, label: str) -> set[str]:
    if not path.is_dir():
        raise FileNotFoundError(f"{label} directory does not exist")
    if any(not entry.is_file() for entry in path.iterdir()):
        raise ValueError(f"{label} directory may contain files only")
    return {entry.name for entry in path.iterdir()}


def figure_artifact_names(*, source_only: bool) -> set[str]:
    sources = {f"{stem}_source_data.csv" for stem in FIGURE_STEMS}
    if source_only:
        return sources
    return sources | {
        f"{stem}.{suffix}" for stem in FIGURE_STEMS for suffix in ("svg", "png")
    }


def validate_reference_root(bundle: Path) -> None:
    if not bundle.is_dir():
        raise FileNotFoundError("reference root does not exist")
    entries = {entry.name: entry for entry in bundle.iterdir()}
    if set(entries) != REFERENCE_DIRECTORIES or any(
        not entry.is_dir() for entry in entries.values()
    ):
        raise ValueError("reference root directory file set is invalid")


def figure_stage_paths(bundle: Path, *, reference_layout: bool) -> list[Path]:
    if reference_layout:
        validate_reference_root(bundle)
    label = "figure source-data" if reference_layout else "figure"
    directory = bundle / ("figure_source_data" if reference_layout else "figures")
    expected = figure_artifact_names(source_only=reference_layout)
    if _directory_files(directory, label) != expected:
        raise ValueError(f"{label} bundle file set is invalid")
    return [directory / name for name in sorted(expected)]


def compare_reference_reports(
    directory: Path, expected_texts: Mapping[str, str],
) -> None:
    expected = set(expected_texts)
    if expected != REPORT_NAMES or _directory_files(directory, "reference report") != expected:
        raise ValueError("reference report file set is invalid")
    for name in sorted(expected):
        if directory.joinpath(name).read_text(encoding="utf-8") != expected_texts[name]:
            raise ValueError(f"immutable reference report mismatch: {name}")


def validate_figure_artifacts(
    directory: Path,
    tables: Mapping[str, pd.DataFrame],
    contract: object,
    case_ids: Sequence[str],
    figure_module: object,
    *,
    source_only: bool,
) -> None:
    figures = getattr(contract, "figures")
    stems = tuple(spec.stem for spec in figures.values())
    if stems != FIGURE_STEMS:
        raise ValueError("figure contract stems do not match the release contract")
    expected = figure_artifact_names(source_only=source_only)
    label = "figure source-data" if source_only else "figure"
    if _directory_files(directory, label) != expected:
        raise ValueError(f"{label} bundle file set is invalid")
    source_columns = list(getattr(figure_module, "SOURCE_COLUMNS"))
    for name in expected:
        path = directory / name
        if name.endswith("_source_data.csv"):
            if pd.read_csv(path, nrows=0).columns.tolist() != source_columns:
                raise ValueError("figure source-data schema is invalid")
        elif name.endswith(".svg"):
            if "<text" not in path.read_text(encoding="utf-8"):
                raise ValueError("SVG text must remain editable")
        elif name.endswith(".png"):
            from PIL import Image

            with Image.open(path) as image:
                if abs(float(image.info.get("dpi", (0,))[0]) - contract.png_dpi) >= 1:
                    raise ValueError("PNG DPI does not match figure contract")
    builders = (
        ("figure1", figure_module.build_figure1, (tables, contract)),
        ("figure2", figure_module.build_figure2, (tables, contract)),
        ("figure3", figure_module.build_figure3, (tables, contract)),
        ("figure4", figure_module.build_figure4, (tables, contract, tuple(case_ids))),
    )
    rebuilt: dict[str, pd.DataFrame] = {}
    for figure_id, builder, values in builders:
        figure, source = builder(*values)
        rebuilt[figures[figure_id].stem] = source
        figure_module._style().close(figure)
    validation.validate_figure_source_files(directory, rebuilt)
