from importlib import import_module, util
from inspect import signature
from pathlib import Path

import pytest

from oips_repro import io, validation
from oips_repro.cli import _bundle_stage_paths


FIGURE_STEMS = (
    "repository_summary_figure_1_candidate_landscape",
    "repository_summary_figure_2_qc_and_orel_ablation",
    "repository_summary_figure_3_posthoc_evidence",
    "repository_summary_figure_4_representative_cases",
)
REPORT_NAMES = {
    "analysis_report.md",
    "numeric_crosscheck.md",
    "rebuild_report.md",
    "validation_report.md",
}


def _touch_files(directory: Path, names: set[str]) -> None:
    directory.mkdir(parents=True)
    for name in names:
        (directory / name).write_text("placeholder\n", encoding="utf-8")


def _touch_tabular_stages(bundle: Path) -> None:
    _touch_files(bundle / "clustering", set(io.CLUSTER_FILE_SCHEMAS))
    _touch_files(bundle / "static", set(io.SCORE_FILE_SCHEMAS))
    _touch_files(bundle / "analysis", set(validation.ANALYSIS_SCHEMAS))


def test_reference_stage_inventory_uses_only_four_figure_source_csvs(tmp_path: Path):
    assert "reference_layout" in signature(_bundle_stage_paths).parameters
    bundle = tmp_path / "reference"
    _touch_tabular_stages(bundle)
    source_names = {f"{stem}_source_data.csv" for stem in FIGURE_STEMS}
    _touch_files(bundle / "figure_source_data", source_names)
    _touch_files(bundle / "reports", REPORT_NAMES)

    paths = _bundle_stage_paths(bundle, reference_layout=True)
    relative = {path.relative_to(bundle).as_posix() for path in paths}

    assert len(relative) == 31
    assert {path for path in relative if path.startswith("figure_source_data/")} == {
        f"figure_source_data/{name}" for name in source_names
    }
    assert not any(path.startswith("figures/") for path in relative)

    (bundle / "figure_source_data" / "repository_summary_figure_1_candidate_landscape.png").touch()
    with pytest.raises(ValueError, match="figure source-data.*file set"):
        _bundle_stage_paths(bundle, reference_layout=True)
    (bundle / "figure_source_data" / "repository_summary_figure_1_candidate_landscape.png").unlink()
    (bundle / "run_manifest.json").touch()
    with pytest.raises(ValueError, match="reference root.*file set"):
        _bundle_stage_paths(bundle, reference_layout=True)


def test_nonreference_stage_inventory_still_requires_all_twelve_figure_files(
    tmp_path: Path,
):
    assert "reference_layout" in signature(_bundle_stage_paths).parameters
    bundle = tmp_path / "reproduced"
    _touch_tabular_stages(bundle)
    source_names = {f"{stem}_source_data.csv" for stem in FIGURE_STEMS}
    _touch_files(bundle / "figures", source_names)

    with pytest.raises(ValueError, match="figure.*file set"):
        _bundle_stage_paths(bundle, reference_layout=False)

    for stem in FIGURE_STEMS:
        (bundle / "figures" / f"{stem}.svg").touch()
        (bundle / "figures" / f"{stem}.png").touch()
    assert len(_bundle_stage_paths(bundle, reference_layout=False)) == 39


def test_reference_report_comparison_is_exact_and_read_only(tmp_path: Path):
    spec = util.find_spec("oips_repro.release")
    assert spec is not None, "reference-release helper is not implemented"
    release = import_module("oips_repro.release")
    reports = {
        "analysis_report.md": "analysis\n",
        "numeric_crosscheck.md": "numeric\n",
        "rebuild_report.md": "rebuild\n",
        "validation_report.md": "validation\n",
    }
    directory = tmp_path / "reference" / "reports"
    _touch_files(directory, set(reports))
    for name, text in reports.items():
        (directory / name).write_text(text, encoding="utf-8", newline="\n")
    before = {name: (directory / name).read_bytes() for name in reports}

    release.compare_reference_reports(directory, reports)

    assert {name: (directory / name).read_bytes() for name in reports} == before
    (directory / "analysis_report.md").write_text("tampered\n", encoding="utf-8")
    with pytest.raises(ValueError, match="report mismatch"):
        release.compare_reference_reports(directory, reports)
    (directory / "analysis_report.md").write_text(reports["analysis_report.md"], encoding="utf-8")
    (directory / "unexpected.md").touch()
    with pytest.raises(ValueError, match="report file set"):
        release.compare_reference_reports(directory, reports)
