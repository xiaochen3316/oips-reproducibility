from __future__ import annotations

from dataclasses import replace
import os
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import pytest

from oips_repro.cli import (
    ANALYSIS_SCHEMAS,
    _figure_runtime,
    _validate_analysis_csv,
    _validate_analysis_relationships,
)
from oips_repro.figures import (
    SOURCE_COLUMNS,
    build_figure1,
    build_figure2,
    build_figure3,
    build_figure4,
    load_figure_contract,
    save_figure,
)


ROOT = Path(__file__).parents[2]


def _tables() -> dict[str, pd.DataFrame]:
    target = pd.DataFrame(
        {
            "pdb_id": ["AAAA", "BBBB"],
            "Pocket_category": ["orthosteric", "allosteric"],
            "formal_cluster_v2_count": [2, 1],
            "reference_first_rank": [1.0, 2.0],
            "first_supported_rank": [1.0, 1.0],
            "reference_top1": [True, False],
            "reference_top3": [True, True],
            "reference_top5": [True, True],
            "first_supported_top1": [True, True],
            "first_supported_top3": [True, True],
            "reference_evaluable": [True, False],
            "first_supported_evaluable": [True, True],
        }
    )
    metric_names = [
        "Reference_Top1", "Reference_Top3", "Reference_Top5", "Reference_MRR",
        "First_supported_Top1", "First_supported_Top3", "First_supported_MRR",
    ]
    bootstrap = pd.DataFrame(
        [
            [method, metric, 0.75, 0.5, 1.0, 10, 1]
            for method in ["target_resampling", "family_clustered"]
            for metric in metric_names
        ],
        columns=[
            "bootstrap_method", "metric", "point_estimate", "CI_2.5_percent",
            "CI_97.5_percent", "iterations", "random_seed",
        ],
    )
    ablation = pd.DataFrame(
        {
            "pdb_id": ["AAAA", "BBBB"],
            "Full_reference_rank": [1.0, 2.0],
            "Without_O_rel_reference_rank": [2.0, 1.0],
            "reference_rank_change_without_minus_full": [1.0, -1.0],
            "Full_first_supported_rank": [1.0, 1.0],
            "Without_O_rel_first_supported_rank": [1.0, 2.0],
            "first_supported_rank_change_without_minus_full": [0.0, 1.0],
            "rank_Spearman_rho": [0.8, 0.9],
            "top_cluster_identity_changed": [True, False],
        }
    )
    labels = pd.DataFrame(
        {
            "pdb_id": ["AAAA", "AAAA", "BBBB"],
            "cluster_v2_id": ["AAAA_V2C001", "AAAA_V2C002", "BBBB_V2C001"],
            "Within_PDB_rank": [1.0, 2.0, 1.0],
            "automated_evidence_label": ["R_auto", "A_auto", "U_auto"],
        }
    )
    md = pd.DataFrame(
        {
            "pdb_id": ["AAAA", "BBBB"],
            "MD_run": ["run1", ""],
            "Concordance_call": ["concordant", "MD_not_available"],
        }
    )
    redocking = pd.DataFrame(
        {"pdb_id": ["AAAA", "BBBB"], "Raw_ligand_RMSD_A": [1.0, 3.2]}
    )
    category = pd.DataFrame(
        {
            "group": ["allosteric", "orthosteric"],
            "Reference_Top1": [0.0, 1.0],
            "First_supported_Top1": [1.0, 1.0],
        }
    )
    master = pd.DataFrame(
        {
            "pdb_id": ["AAAA", "BBBB", "CCCC"],
            "cluster_v2_id": ["AAAA_V2C001", "BBBB_V2C001", "CCCC_V2C001"],
            "C_cons": [90.0, 80.0, 75.0],
            "G_geo": [80.0, 70.0, 65.0],
            "P_lig": [np.nan, 60.0, 55.0],
            "O_rel_formal": [70.0, 50.0, 45.0],
            "Q_evidence": [60.0, 40.0, 35.0],
        }
    )
    cases = pd.DataFrame(
        {
            "pdb_id": ["AAAA", "BBBB", "CCCC"],
            "static_top_cluster_v2_id": ["AAAA_V2C001", "BBBB_V2C001", "CCCC_V2C001"],
            "static_top_score": [80.0, 70.0, 65.0],
            "static_top_label": ["R_auto", "A_auto", "U_auto"],
            "static_top_tool_support": [2, 2, 1],
            "static_top_chains": ["A;B", "A;B", "A"],
            "static_top_interface_overlap": [0.5, 0.4, 0.2],
            "reference_cluster_v2_id": ["AAAA_V2C001", "BBBB_V2C001", "CCCC_V2C001"],
            "reference_rank": [1.0, 1.0, 2.0],
            "reference_DCC_A": [2.0, 3.0, 4.0],
            "without_O_rel_reference_rank": [2.0, 1.0, 2.0],
            "without_O_rel_top_cluster_v2_id": ["AAAA_V2C001", "BBBB_V2C001", "CCCC_V2C001"],
            "MD_calls": ["concordant", "MD_not_available", "partially_concordant"],
            "redocking_RMSD_A": [1.0, 3.2, 2.5],
            "redocking_status": ["RMSD <= 2 A", "RMSD > 3 A", "2 A < RMSD <= 3 A"],
        }
    )
    qc = pd.DataFrame({"QC_status": ["QC_pass", "QC_boundary_sensitive"]})
    return {
        "target_summary": target,
        "bootstrap": bootstrap,
        "topk_qc": qc,
        "orel_targets": ablation,
        "evidence_labels": labels,
        "md_mapping": md,
        "redocking_mapping": redocking,
        "category_metrics": category,
        "master": master,
        "representative_cases": cases,
    }


def test_figure_contract_builders_return_all_panels_and_common_source_schema():
    contract = load_figure_contract(ROOT / "config" / "figure_contract.yaml")
    synthetic_contract = replace(
        contract, case_titles={"AAAA": "Case A", "BBBB": "Case B", "CCCC": "Case C"}
    )
    tables = _tables()
    builders = [
        (build_figure1, (tables, contract)),
        (build_figure2, (tables, contract)),
        (build_figure3, (tables, contract)),
        (build_figure4, (tables, synthetic_contract, ("AAAA", "BBBB", "CCCC"))),
    ]

    for builder, args in builders:
        fig, source = builder(*args)
        figure_id = builder.__name__.replace("build_", "")
        assert source.columns.tolist() == list(SOURCE_COLUMNS)
        assert source.equals(
            source.sort_values(["panel", "record_id", "series"], kind="mergesort").reset_index(drop=True)
        )
        active_contract = synthetic_contract if figure_id == "figure4" else contract
        assert set(source["panel"]) == set(active_contract.figures[figure_id].panels)
        assert len(fig.axes) == len(active_contract.figures[figure_id].panels)
        for axis, panel in zip(fig.axes, active_contract.figures[figure_id].panels.values()):
            assert axis.get_xlabel() == panel.x_label
            assert axis.get_ylabel() == panel.y_label
        source_groups = {(row.panel, row.series) for row in source.itertuples()}
        artist_groups = {
            tuple(str(artist.get_gid()).removeprefix("source:").split(":", 1))
            for axis in fig.axes
            for artist in axis.get_children()
            if str(artist.get_gid()).startswith("source:")
        }
        assert artist_groups
        assert artist_groups.issubset(source_groups)
        plt.close(fig)

    fig4, source4 = build_figure4(
        tables, synthetic_contract, ("AAAA", "BBBB", "CCCC")
    )
    missing = source4.loc[
        source4["record_id"].eq("AAAA") & source4["series"].eq("P_lig"), "x"
    ]
    assert len(missing) == 1 and missing.isna().all()
    assert any(text.get_text() == "NA" for axis in fig4.axes for text in axis.texts)
    plt.close(fig4)


def test_md_missing_annotation_is_computed_from_mapping():
    contract = load_figure_contract(ROOT / "config" / "figure_contract.yaml")
    tables = _tables()
    fig, source = build_figure3(tables, contract)
    annotations = source.loc[source["panel"].eq("b"), "annotation"].dropna().tolist()
    assert "1 targets without MD" in annotations
    tables["md_mapping"] = pd.concat(
        [tables["md_mapping"], pd.DataFrame([["CCCC", "", "MD_not_available"]], columns=tables["md_mapping"].columns)],
        ignore_index=True,
    )
    _, changed = build_figure3(tables, contract)
    assert "2 targets without MD" in changed.loc[changed["panel"].eq("b"), "annotation"].dropna().tolist()
    assert "insufficient_MD_evidence" in set(changed.loc[changed["panel"].eq("b"), "record_id"])
    assert changed.loc[
        changed["panel"].eq("b") & changed["record_id"].eq("insufficient_MD_evidence"), "x"
    ].item() == 0
    plt.close(fig)


def test_figure1_strictly_parses_textual_false_endpoint_values():
    contract = load_figure_contract(ROOT / "config" / "figure_contract.yaml")
    tables = _tables()
    tables["target_summary"]["reference_top1"] = ["False", "True"]

    fig, source = build_figure1(tables, contract)

    value = source.loc[
        source["panel"].eq("c")
        & source["record_id"].eq("Top-1")
        & source["series"].eq("Reference-associated"),
        "y",
    ].item()
    assert value == 0.0
    plt.close(fig)


def test_save_figure_writes_editable_svg_and_600_dpi_png(tmp_path: Path):
    contract = load_figure_contract(ROOT / "config" / "figure_contract.yaml")
    fig, source = build_figure1(_tables(), contract)
    spec = contract.figures["figure1"]

    save_figure(fig, source, tmp_path, spec, contract)

    assert {path.name for path in tmp_path.iterdir()} == {
        f"{spec.stem}.svg", f"{spec.stem}.png", f"{spec.stem}_source_data.csv",
    }
    assert "<text" in (tmp_path / f"{spec.stem}.svg").read_text(encoding="utf-8")
    with Image.open(tmp_path / f"{spec.stem}.png") as image:
        dpi = image.info["dpi"]
    assert abs(dpi[0] - 600) < 1 and abs(dpi[1] - 600) < 1


def test_source_coordinates_match_figure3_bars_and_figure4_module_axes():
    contract = load_figure_contract(ROOT / "config" / "figure_contract.yaml")
    tables = _tables()
    fig3, source3 = build_figure3(tables, contract)
    source_r = source3.loc[
        source3["panel"].eq("a") & source3["series"].eq("R_auto")
    ].sort_values("x")
    patches = [
        patch for patch in fig3.axes[0].patches
        if patch.get_gid() == "source:a:R_auto"
    ]
    assert [patch.get_x() + patch.get_width() / 2 for patch in patches] == source_r["x"].tolist()
    assert [patch.get_height() for patch in patches] == source_r["y"].tolist()
    plt.close(fig3)

    synthetic_contract = replace(
        contract, case_titles={"AAAA": "Case A", "BBBB": "Case B", "CCCC": "Case C"}
    )
    fig4, source4 = build_figure4(
        tables, synthetic_contract, ("AAAA", "BBBB", "CCCC")
    )
    c_cons = source4.loc[
        source4["panel"].eq("a") & source4["series"].eq("C_cons")
    ].iloc[0]
    patch = next(
        patch for patch in fig4.axes[0].patches
        if patch.get_gid() == "source:a:C_cons"
    )
    assert c_cons["x"] == 90.0
    assert c_cons["y"] == 0
    assert patch.get_width() == c_cons["x"]
    assert patch.get_y() + patch.get_height() / 2 == c_cons["y"]
    plt.close(fig4)


def test_importing_figures_creates_no_files(tmp_path: Path):
    mpl_config = tmp_path / "mpl-config"
    mpl_config.mkdir()
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    environment["MPLCONFIGDIR"] = str(mpl_config)
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import oips_repro.figures; "
            "raise SystemExit(1 if 'matplotlib.pyplot' in sys.modules else 0)",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert list(mpl_config.iterdir()) == []


def test_importing_cli_does_not_eagerly_import_matplotlib(tmp_path: Path):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import oips_repro.cli; "
            "raise SystemExit(1 if 'matplotlib.pyplot' in sys.modules else 0)",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert list(tmp_path.iterdir()) == []


def test_cli_figure_loader_uses_a_writable_ephemeral_matplotlib_cache(tmp_path: Path):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            "import os; from pathlib import Path; from oips_repro.cli import _figure_runtime; "
            "before=os.environ.get('MPLCONFIGDIR'); "
            "ctx=_figure_runtime(); module=ctx.__enter__(); cache=Path(os.environ['MPLCONFIGDIR']); "
            "module._style(); live=cache.is_dir(); ctx.__exit__(None,None,None); "
            "restored=os.environ.get('MPLCONFIGDIR')==before; "
            "raise SystemExit(0 if live and restored and not cache.exists() else 1)",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "font_manager cache" not in completed.stderr
    assert list(tmp_path.iterdir()) == []


def _relationship_frames() -> dict[str, pd.DataFrame]:
    targets = ["AAAA", "BBBB", "CCCC"]
    categories = ["cat-a", "cat-b", "cat-c"]
    families = ["family-a", "family-b", "family-c"]
    target = pd.DataFrame({
        "pdb_id": targets,
        "Pocket_category": categories,
        "Protein_family": families,
    })
    master_rows = []
    for pdb_id in targets:
        for rank in (1.0, 2.0, 3.0, 4.0):
            master_rows.append({
                "pdb_id": pdb_id,
                "cluster_v2_id": f"{pdb_id}_V2C{int(rank):03d}",
                "Within_PDB_rank": rank,
            })
    master = pd.DataFrame(master_rows)
    qc = master.loc[master["Within_PDB_rank"].le(3)].copy()
    return {
        "target_level_candidate_prioritization.csv": target,
        "final_cluster_v2_master_table.csv": master,
        "final_md_cluster_v2_mapping.csv": pd.DataFrame({
            "pdb_id": targets,
            "Simulation_context": ["", "native", "native"],
            "MD_run": ["", "run-b", "run-c"],
        }),
        "final_top3_automated_QC.csv": qc,
        "final_redocking_cluster_v2_mapping.csv": pd.DataFrame({"pdb_id": targets}),
        "final_orel_ablation_targets.csv": pd.DataFrame({"pdb_id": targets}),
        "final_category_metrics.csv": pd.DataFrame({"group": categories}),
        "final_orel_ablation_categories.csv": pd.DataFrame({"Pocket_category": categories}),
        "final_family_sensitivity.csv": pd.DataFrame({"excluded_family": families}),
    }


def _replace_md_target(frames):
    frames["final_md_cluster_v2_mapping.csv"].loc[2, "pdb_id"] = "BBBB"


def _replace_qc_key(frames):
    qc = frames["final_top3_automated_QC.csv"]
    qc.loc[qc["cluster_v2_id"].eq("AAAA_V2C003"), "cluster_v2_id"] = "AAAA_V2C004"


def _replace_values(filename, column, values):
    def mutate(frames):
        frames[filename].loc[:, column] = values
    return mutate


@pytest.mark.parametrize(
    "mutator, message",
    [
        (_replace_md_target, "MD mapping target set"),
        (_replace_qc_key, "QC keys and ranks"),
        (_replace_values("final_redocking_cluster_v2_mapping.csv", "pdb_id", ["AAAA", "AAAA", "BBBB"]), "redocking target set"),
        (_replace_values("final_orel_ablation_targets.csv", "pdb_id", ["AAAA", "AAAA", "BBBB"]), "ablation target set"),
        (_replace_values("final_category_metrics.csv", "group", ["cat-a", "cat-b", "invented"]), "category metric set"),
        (_replace_values("final_orel_ablation_categories.csv", "Pocket_category", ["cat-a", "cat-b", "invented"]), "ablation category set"),
        (_replace_values("final_family_sensitivity.csv", "excluded_family", ["family-a", "family-b", "invented"]), "family sensitivity set"),
    ],
)
def test_analysis_relationship_validation_rejects_cross_table_replacements(mutator, message):
    frames = _relationship_frames()
    mutator(frames)
    if message == "MD mapping target set":
        md = frames["final_md_cluster_v2_mapping.csv"]
        keys = ["pdb_id", "Simulation_context", "MD_run"]
        assert not md.duplicated(keys).any()
        pd.testing.assert_frame_equal(md, md.sort_values(keys).reset_index(drop=True))
    with pytest.raises(ValueError, match=message):
        _validate_analysis_relationships(frames)


def test_analysis_csv_validation_rejects_numeric_enum_key_and_sort_tampering():
    columns = ANALYSIS_SCHEMAS["target_level_candidate_prioritization.csv"]
    rows = []
    for pdb_id in ("AAAA", "BBBB"):
        row = {column: "" for column in columns}
        row.update({
            "pdb_id": pdb_id,
            "Pocket_category": "orthosteric",
            "Protein_family": "family",
            "formal_cluster_v2_count": "2",
            "reference_evaluable": "true",
            "R_auto_cluster_count": "1",
            "reference_first_rank": "1",
            "reference_rank_percentile": "1",
            "reference_top1": "true",
            "reference_top3": "true",
            "reference_top5": "true",
            "reference_reciprocal_rank": "1",
            "first_supported_evaluable": "true",
            "first_supported_rank": "1",
            "first_supported_top1": "true",
            "first_supported_top3": "true",
            "first_supported_reciprocal_rank": "1",
            "static_top_cluster_v2_id": f"{pdb_id}_V2C001",
            "static_top_evidence_label": "R_auto",
        })
        rows.append(row)
    valid = pd.DataFrame(rows, columns=columns)
    parsed = _validate_analysis_csv("target_level_candidate_prioritization.csv", valid)
    assert parsed["formal_cluster_v2_count"].tolist() == [2, 2]

    bad_numeric = valid.copy()
    bad_numeric.loc[0, "reference_first_rank"] = "not-numeric"
    with pytest.raises(ValueError, match="numeric"):
        _validate_analysis_csv("target_level_candidate_prioritization.csv", bad_numeric)
    bad_boolean = valid.copy()
    bad_boolean.loc[0, "reference_top1"] = "yes"
    with pytest.raises(ValueError, match="lowercase true/false"):
        _validate_analysis_csv("target_level_candidate_prioritization.csv", bad_boolean)
    with pytest.raises(ValueError, match="sorted"):
        _validate_analysis_csv(
            "target_level_candidate_prioritization.csv", valid.iloc[::-1].reset_index(drop=True)
        )
    duplicate = pd.concat([valid.iloc[[0]], valid.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        _validate_analysis_csv("target_level_candidate_prioritization.csv", duplicate)

    qc_columns = ANALYSIS_SCHEMAS["final_top3_automated_QC.csv"]
    qc = pd.DataFrame([{column: "" for column in qc_columns}], columns=qc_columns)
    qc.loc[0, ["pdb_id", "cluster_v2_id", "Within_PDB_rank", "QC_status"]] = [
        "AAAA", "AAAA_V2C001", "1", "not-a-qc-state",
    ]
    qc.loc[0, "mappability"] = "center_and_residue_mappable"
    with pytest.raises(ValueError, match="QC_status"):
        _validate_analysis_csv("final_top3_automated_QC.csv", qc)


def test_task6_production_modules_remain_below_700_lines():
    for name in ("cli.py", "statistics.py", "figures.py"):
        lines = (ROOT / "src" / "oips_repro" / name).read_text(encoding="utf-8").splitlines()
        assert len(lines) < 700, f"{name} has {len(lines)} lines"
