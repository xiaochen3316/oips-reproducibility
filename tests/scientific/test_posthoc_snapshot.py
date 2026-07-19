from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

import pandas as pd
import pytest

from oips_repro import io
from oips_repro.cli import main
from oips_repro.posthoc import PosthocInputs, run_posthoc


ROOT = Path(__file__).parents[2]
EXPECTED = json.loads(
    (Path(__file__).parent / "data" / "expected_summary.json").read_text(encoding="utf-8")
)
CONVENIENCE_APPEND = [
    "automated_evidence_label",
    "label_reason",
    "unresolved_flag",
    "QC_status",
    "DCC_A",
    "contact_recall",
    "contact_precision",
    "residue_IoU",
    "R_auto_rule_pass",
]


def _static_rankings(tmp_path: Path) -> tuple[io.ManuscriptConfig, pd.DataFrame]:
    config_path = ROOT / "config" / "manuscript.yaml"
    cluster_dir = tmp_path / "cluster"
    score_dir = tmp_path / "static"
    assert main(["cluster", "--config", str(config_path), "--output", str(cluster_dir)]) == 0
    assert main(
        [
            "score",
            "--config",
            str(config_path),
            "--cluster-dir",
            str(cluster_dir),
            "--output",
            str(score_dir),
        ]
    ) == 0
    return (
        io.load_manuscript_config(config_path),
        pd.read_csv(score_dir / "cluster_v2_static_rankings.csv"),
    )


def test_frozen_posthoc_snapshot_and_static_boundary(tmp_path: Path):
    config, rankings = _static_rankings(tmp_path)
    original = rankings.copy(deep=True)
    inputs = PosthocInputs(
        rankings=rankings,
        structures=io.load_structure_paths(config),
        references=pd.read_csv(config.resolve_configured_path("reference_annotations"), keep_default_na=False),
        md_runs=pd.read_csv(config.resolve_configured_path("md_evidence"), keep_default_na=False),
        redocking=pd.read_csv(config.resolve_configured_path("redocking_evidence"), keep_default_na=False),
        systems=pd.read_csv(config.resolve_configured_path("systems"), sep="\t", keep_default_na=False),
    )

    result = run_posthoc(inputs, config)

    pd.testing.assert_frame_equal(rankings, original)
    assert rankings.columns.tolist() == list(io.RANKING_COLUMNS)

    counts = EXPECTED["posthoc_counts"]
    distributions = EXPECTED["posthoc_distributions"]
    assert len(result.topk_qc) == counts["topk_qc_rows"]
    assert result.topk_qc["QC_status"].value_counts().to_dict() == distributions["QC_status"]

    assert result.evidence_labels["automated_evidence_label"].value_counts().to_dict() == distributions["automated_evidence_label"]
    assert "X_auto" not in set(result.evidence_labels["automated_evidence_label"])
    assert result.evidence_labels.loc[
        result.evidence_labels["Within_PDB_rank"].le(3), "automated_evidence_label"
    ].value_counts().to_dict() == distributions["top3_evidence_label"]

    assert len(result.reference_mapping) == counts["reference_rows"]
    assert result.reference_mapping["reference_selection_status"].value_counts().to_dict() == distributions["reference_selection_status"]
    reference_case = EXPECTED["posthoc_reference_case"]
    selected = result.reference_mapping.loc[
        result.reference_mapping["pdb_id"].eq(reference_case["pdb_id"])
    ]
    assert selected["reference_selection_status"].eq(reference_case["selection_status"]).all()
    assert selected["selected_reference_ligand_key"].eq(reference_case["selected_ligand_key"]).all()
    assert selected["reference_ligand_atom_count"].eq(reference_case["ligand_atom_count"]).all()
    assert result.reference_evaluable.to_dict() == {
        pdb_id: True for pdb_id in sorted(rankings["pdb_id"].unique())
    }

    assert len(result.md_mapping) == counts["md_rows"]
    real_md = result.md_mapping.loc[result.md_mapping["MD_run"].astype(str).ne("")]
    assert len(real_md) == counts["md_real_runs"]
    assert real_md["pdb_id"].nunique() == counts["md_targets_with_runs"]
    assert result.md_mapping["Concordance_call"].value_counts().to_dict() == distributions["Concordance_call"]

    assert len(result.redocking_mapping) == counts["redocking_rows"]
    rmsd = result.redocking_mapping["Raw_ligand_RMSD_A"]
    rmsd_bins = distributions["redocking_RMSD_bins"]
    assert int(rmsd.le(2.0).sum()) == rmsd_bins["le_2A"]
    assert int(rmsd.gt(2.0).mul(rmsd.le(3.0)).sum()) == rmsd_bins["gt_2_le_3A"]
    assert int(rmsd.gt(3.0).sum()) == rmsd_bins["gt_3A"]

    assert len(result.unresolved_cases) == counts["unresolved_rows"]
    assert result.unresolved_cases["issue_type"].value_counts().to_dict() == distributions["unresolved_issue_type"]
    assert result.representative_case_ids == tuple(EXPECTED["representative_case_order"])

    assert len(result.convenience_master) == EXPECTED["clusters"]
    assert result.convenience_master.columns.tolist() == [
        *io.RANKING_COLUMNS,
        *CONVENIENCE_APPEND,
    ]
    pd.testing.assert_frame_equal(
        result.convenience_master.loc[:, list(io.RANKING_COLUMNS)].reset_index(drop=True),
        original.reset_index(drop=True),
        check_dtype=False,
    )

    duplicate_md = inputs.md_runs.copy(deep=True)
    duplicate_md.iloc[1] = duplicate_md.iloc[0]
    duplicate_md.loc[1, "pdb_id"] = str(duplicate_md.loc[0, "pdb_id"]).lower()
    with pytest.raises(ValueError, match="duplicate keys"):
        run_posthoc(replace(inputs, md_runs=duplicate_md), config)
    with pytest.raises(ValueError, match="missing required columns.*pdb_id"):
        run_posthoc(replace(inputs, md_runs=inputs.md_runs.drop(columns="pdb_id")), config)
