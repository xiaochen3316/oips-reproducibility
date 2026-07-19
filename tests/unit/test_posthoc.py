from __future__ import annotations

from dataclasses import replace
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from oips_repro import io
from oips_repro.posthoc import (
    MD_COLUMNS,
    MDRules,
    QCRules,
    REFERENCE_COLUMNS,
    ReferenceRules,
    TOPK_QC_COLUMNS,
    assign_evidence_labels,
    derive_reference_evaluable,
    map_md_evidence,
    map_reference_evidence,
    map_redocking_evidence,
    reference_rule_pass,
    representative_case_ids,
    run_topk_qc,
)


ROOT = Path(__file__).parents[2]


def _pdb_atom(
    *,
    record: str = "ATOM",
    serial: int = 1,
    atom_name: str = "CA",
    residue_name: str = "ALA",
    chain: str = "A",
    residue_number: int = 1,
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    element: str = "C",
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom_name:^4} {residue_name:>3} {chain:1}"
        f"{residue_number:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{1.0:>6.2f}{20.0:>6.2f}          {element:>2}\n"
    )


def _qc_row(cluster_id: str, rank: float, **updates: object) -> dict[str, object]:
    row: dict[str, object] = {
        "pdb_id": "TEST",
        "cluster_v2_id": cluster_id,
        "Within_PDB_rank": rank,
        "OIPS-P_static_recomputed": 80.0 - rank,
        "mappability": "center_and_residue_mappable",
        "cluster_diameter_A": 5.0,
        "spatial_continuity": True,
        "tool_support_count": 2,
        "center_dispersion_A": 1.0,
        "core_envelope_ratio": 0.5,
        "contributing_chain_count": 2,
        "dominant_chain_fraction": 0.5,
        "interface_fraction": 0.4,
        "distance_to_interface_A": 2.0,
        "same_tool_secondary_unit_count": 0,
        "G_geo": 70.0,
        "P_lig": 70.0,
        "center_available_representatives": 2,
        "residue_available_representatives": 2,
        "pairwise_residue_iou_median": 0.5,
        "boundary_sensitive": False,
        "envelope_residues": "A:ALA:1:",
    }
    row.update(updates)
    return row


@pytest.mark.parametrize(
    ("dcc", "recall", "iou", "expected"),
    [
        (6.0, 0.10, 0.0, True),
        (6.0, np.nextafter(0.10, 0.0), 0.0, False),
        (np.nextafter(6.0, np.inf), 0.10, 0.0, False),
        (10.0, 0.20, 0.0, True),
        (10.0, 0.0, 0.15, True),
        (10.0, np.nextafter(0.20, 0.0), np.nextafter(0.15, 0.0), False),
        (np.nextafter(10.0, np.inf), 0.20, 0.15, False),
        (12.0, 0.50, 0.15, True),
        (12.0, np.nextafter(0.50, 0.0), 0.15, False),
        (12.0, 0.50, np.nextafter(0.15, 0.0), False),
        (np.nextafter(12.0, np.inf), 0.50, 0.15, False),
        (np.nan, 1.0, 1.0, False),
        (1.0, np.nan, 1.0, False),
    ],
)
def test_reference_rule_boundaries_are_inclusive_only_at_the_frozen_values(
    dcc: float, recall: float, iou: float, expected: bool
):
    assert reference_rule_pass(dcc, recall, iou, ReferenceRules()) is expected


def test_reference_evaluable_depends_on_available_reference_not_an_r_label():
    mapping = pd.DataFrame(
        {
            "pdb_id": ["TEST", "TEST", "MISS"],
            "selected_reference_ligand_key": ["LIG:A:1:", "LIG:A:1:", ""],
            "reference_ligand_atom_count": [8, 8, 0],
            "DCC_A": [np.nan, np.nan, np.nan],
            "R_auto_rule_pass": [False, False, False],
        }
    )

    evaluable = derive_reference_evaluable(mapping)

    assert evaluable.to_dict() == {"MISS": False, "TEST": True}


def test_reference_override_uses_the_curated_after_key(tmp_path: Path):
    structure = tmp_path / "override.pdb"
    structure.write_text(
        _pdb_atom(record="HETATM", residue_name="NAP", chain="C", residue_number=501, x=20.0)
        + _pdb_atom(record="HETATM", serial=2, residue_name="6VN", chain="D", residue_number=503)
        + _pdb_atom(serial=3, residue_name="ALA", chain="A", residue_number=1, z=4.5),
        encoding="utf-8",
    )
    rankings = pd.DataFrame(
        [{
            "pdb_id": "TEST", "cluster_v2_id": "TEST_V2C001", "Within_PDB_rank": 1.0,
            "medoid_center_x": 0.0, "medoid_center_y": 0.0, "medoid_center_z": 0.0,
            "envelope_residues": "A:ALA:1:", "interface_fraction": 1.0,
            "distance_to_interface_A": 0.0, "contributing_chain_count": 1,
        }]
    )
    references = pd.DataFrame(
        [["TEST", "NAP:C:501:", "ok", "6VN:D:503:", "True", "DEC-001"]],
        columns=[
            "pdb_id", "selected_ligand", "ligand_status", "reference_ligand_annotation",
            "override_applied", "decision_id",
        ],
    )

    mapped = map_reference_evidence(rankings, {"TEST": structure}, references)

    assert mapped.loc[0, "requested_reference_ligand_key"] == "6VN:D:503:"
    assert mapped.loc[0, "selected_reference_ligand_key"] == "6VN:D:503:"


def test_reference_mapping_rejects_a_nonfinite_selected_ligand_centroid(tmp_path: Path):
    structure = tmp_path / "nonfinite-ligand.pdb"
    structure.write_text(
        _pdb_atom(record="HETATM", residue_name="LIG", chain="A", residue_number=1, x=np.nan)
        + _pdb_atom(serial=2, residue_name="ALA", chain="B", residue_number=2),
        encoding="utf-8",
    )
    rankings = pd.DataFrame([{
        "pdb_id": "TEST", "cluster_v2_id": "TEST_V2C001", "Within_PDB_rank": 1.0,
        "medoid_center_x": 0.0, "medoid_center_y": 0.0, "medoid_center_z": 0.0,
        "envelope_residues": "B:ALA:2:", "interface_fraction": 1.0,
        "distance_to_interface_A": 0.0, "contributing_chain_count": 1,
    }])
    references = pd.DataFrame(
        [["TEST", "LIG:A:1:", "ok", "LIG:A:1:", False, "not_applicable"]],
        columns=[
            "pdb_id", "selected_ligand", "ligand_status", "reference_ligand_annotation",
            "override_applied", "decision_id",
        ],
    )

    with pytest.raises(ValueError, match="ligand centroid.*finite"):
        map_reference_evidence(rankings, {"TEST": structure}, references)


def test_qc_parses_text_false_and_applies_precedence_and_strict_thresholds(
    tmp_path: Path,
):
    structure = tmp_path / "test.pdb"
    structure.write_text(_pdb_atom(), encoding="utf-8")
    rows = [
        _qc_row("TEST_V2C006", 6, mappability="invalid", spatial_continuity="False"),
        _qc_row("TEST_V2C002", 2, spatial_continuity="False"),
        _qc_row("TEST_V2C003", 3, center_dispersion_A=np.nextafter(4.0, np.inf)),
        _qc_row("TEST_V2C004", 4, boundary_sensitive=True, center_dispersion_A=4.0),
        _qc_row("TEST_V2C005", 5, mappability="center_only_mappable"),
        _qc_row(
            "TEST_V2C001",
            1,
            cluster_diameter_A=10.5,
            core_envelope_ratio=0.25,
            pairwise_residue_iou_median=0.15,
        ),
    ]

    qc = run_topk_qc(pd.DataFrame(rows), {"TEST": structure}, top_k=6, rules=QCRules())

    statuses = qc.set_index("cluster_v2_id")["QC_status"].to_dict()
    assert statuses == {
        "TEST_V2C001": "QC_pass",
        "TEST_V2C002": "QC_possible_overmerge",
        "TEST_V2C003": "QC_possible_split",
        "TEST_V2C004": "QC_boundary_sensitive",
        "TEST_V2C005": "QC_insufficient_evidence",
        "TEST_V2C006": "QC_unmappable",
    }
    excluded = qc.set_index("cluster_v2_id")["clear_exclusion_flag"]
    assert bool(excluded["TEST_V2C002"])
    assert bool(excluded["TEST_V2C006"])
    assert qc["cluster_v2_id"].tolist() == [f"TEST_V2C00{i}" for i in range(1, 7)]


def test_qc_rejects_malformed_external_boolean(tmp_path: Path):
    structure = tmp_path / "test.pdb"
    structure.write_text(_pdb_atom(), encoding="utf-8")
    ranking = pd.DataFrame([_qc_row("TEST_V2C001", 1, spatial_continuity="not-a-bool")])

    with pytest.raises(ValueError, match="spatial_continuity"):
        run_topk_qc(ranking, {"TEST": structure})


def _md_rankings() -> pd.DataFrame:
    top_residues = ";".join(f"A:ALA:{number}:" for number in range(1, 21))
    alternative_residues = ";".join(f"B:GLY:{number}:" for number in range(101, 121))
    return pd.DataFrame(
        [
            {
                "pdb_id": "TEST",
                "cluster_v2_id": "TEST_V2C002",
                "Within_PDB_rank": 1.0,
                "envelope_residues": alternative_residues,
                "medoid_center_x": 20.0,
                "medoid_center_y": 0.0,
                "medoid_center_z": 0.0,
            },
            {
                "pdb_id": "TEST",
                "cluster_v2_id": "TEST_V2C001",
                "Within_PDB_rank": 1.0,
                "envelope_residues": top_residues,
                "medoid_center_x": 0.0,
                "medoid_center_y": 0.0,
                "medoid_center_z": 0.0,
            },
            {
                "pdb_id": "MISS",
                "cluster_v2_id": "MISS_V2C001",
                "Within_PDB_rank": 1.0,
                "envelope_residues": "M:ALA:1:",
                "medoid_center_x": 0.0,
                "medoid_center_y": 0.0,
                "medoid_center_z": 0.0,
            },
        ]
    )


def _md_runs() -> pd.DataFrame:
    four_top_and_many_other = [f"A:ALA:{number}:" for number in range(1, 5)] + [
        f"C:SER:{number}:" for number in range(201, 297)
    ]
    return pd.DataFrame(
        [
            ["TEST", "concordant-iou", "native", "A:ALA:1:;A:ALA:2:;A:ALA:3:;A:ALA:4:", "0;0;0", 1],
            ["TEST", "concordant-precision", "native", ";".join(four_top_and_many_other), "6;0;0", 2],
            ["TEST", "partial", "native", "A:ALA:1:", "20;0;0", 3],
            ["TEST", "alternative", "native", "B:GLY:101:", "20;0;0", 4],
            ["TEST", "boundary", "native", "C:SER:999:", "8;0;0", 5],
            ["TEST", "fallback", "native", "C:SER:998:", "10;0;0", 6],
            ["TEST", "apo", "apo", "", "", ""],
            ["TEST", "insufficient", "native", "", "", ""],
        ],
        columns=[
            "pdb_id",
            "MD_run",
            "Simulation_context",
            "Persistent_MD_contact_residues",
            "MD_contact_center",
            "D_dyn_run_score",
        ],
    )


def test_md_mapping_covers_every_branch_and_uses_deterministic_tie_breaks():
    rankings = _md_rankings()
    runs = _md_runs()

    first = map_md_evidence(rankings, runs, rules=MDRules())
    second = map_md_evidence(
        rankings.sample(frac=1, random_state=7),
        runs.sample(frac=1, random_state=8),
        rules=MDRules(),
    )

    pd.testing.assert_frame_equal(first, second)
    calls = first.set_index("MD_run")["Concordance_call"].to_dict()
    assert calls == {
        "": "MD_not_available",
        "alternative": "static_dynamic_conflict",
        "apo": "apo_only_context",
        "boundary": "boundary_shift",
        "concordant-iou": "concordant",
        "concordant-precision": "concordant",
        "fallback": "static_dynamic_conflict",
        "insufficient": "insufficient_MD_evidence",
        "partial": "partially_concordant",
    }
    assert first.loc[first["MD_run"].eq("alternative"), "best_MD_mapped_cluster_v2_id"].item() == "TEST_V2C002"
    assert first.loc[first["MD_run"].eq("boundary"), "best_MD_mapped_cluster_v2_id"].item() == "TEST_V2C001"
    assert first.loc[first["MD_run"].eq("fallback"), "best_MD_mapped_cluster_v2_id"].item() == "TEST_V2C001"
    assert first.loc[first["MD_run"].eq("apo"), "source_persistent_contacts"].item() == (
        "previously_chain_reconciled_MD_contact_region"
    )
    assert first.loc[first["Concordance_call"].eq("MD_not_available"), "source_persistent_contacts"].item() == (
        "not_available"
    )


def test_md_mapping_rejects_case_normalized_duplicate_keys_and_missing_target_column():
    duplicate = _md_runs().iloc[[0, 0]].copy().reset_index(drop=True)
    duplicate.loc[1, "pdb_id"] = "test"

    with pytest.raises(ValueError, match="duplicate keys"):
        map_md_evidence(_md_rankings(), duplicate)
    with pytest.raises(ValueError, match="missing required columns.*pdb_id"):
        map_md_evidence(_md_rankings(), duplicate.drop(columns="pdb_id"))


def test_evidence_label_dynamic_support_uses_supplied_md_rules():
    rankings = pd.DataFrame([{
        "pdb_id": "TEST", "cluster_v2_id": "TEST_V2C001", "Within_PDB_rank": 1.0,
        "tool_support_count": 1, "interface_fraction": np.nan,
        "distance_to_interface_A": np.nan, "spatial_continuity": True,
        "mappability": "center_and_residue_mappable",
    }])
    reference_row = {column: np.nan for column in REFERENCE_COLUMNS}
    reference_row.update({
        "pdb_id": "TEST", "cluster_v2_id": "TEST_V2C001", "Within_PDB_rank": 1.0,
        "R_auto_rule_pass": False, "reference_selection_unresolved": False,
    })
    md_row = {column: np.nan for column in MD_COLUMNS}
    md_row.update({
        "pdb_id": "TEST", "MD_run": "run-1", "Simulation_context": "native",
        "persistent_MD_contact_residue_count": 1,
        "best_MD_mapped_cluster_v2_id": "TEST_V2C001",
        "best_MD_cluster_Jaccard": 0.10,
        "best_MD_cluster_center_distance_A": 4.0,
    })
    reference = pd.DataFrame([reference_row], columns=REFERENCE_COLUMNS)
    qc = pd.DataFrame(columns=TOPK_QC_COLUMNS)
    md = pd.DataFrame([md_row], columns=MD_COLUMNS)

    default = assign_evidence_labels(rankings, reference, qc, md)
    strict = assign_evidence_labels(
        rankings,
        reference,
        qc,
        md,
        rules=replace(MDRules(), alternative_iou_min=0.20, alternative_center_max_A=3.0),
    )

    assert default.loc[0, "automated_evidence_label"] == "A_auto"
    assert strict.loc[0, "automated_evidence_label"] == "U_auto"


def test_redocking_recomputes_threshold_categories_at_two_and_three_angstroms():
    labels = pd.DataFrame(
        [
            ["A000", "A000_V2C001", 1.0, "R_auto"],
            ["A001", "A001_V2C001", 1.0, "R_auto"],
            ["A002", "A002_V2C001", 2.0, "R_auto"],
            ["A003", "A003_V2C001", 3.0, "R_auto"],
        ],
        columns=["pdb_id", "cluster_v2_id", "Within_PDB_rank", "automated_evidence_label"],
    )
    redocking = pd.DataFrame(
        [
            ["A000", "tool", 0.0, 0.0, "contradictory", "zero"],
            ["A001", "tool", -1.0, 2.0, "contradictory", "none"],
            ["A002", "tool", -2.0, 3.0, "contradictory", "review"],
            ["A003", "tool", -3.0, 3.0001, "contradictory", "warning"],
        ],
        columns=[
            "pdb_id",
            "Docking_software",
            "GlideScore_kcal_per_mol",
            "Raw_ligand_RMSD_A",
            "RMSD_threshold_call",
            "Failure_or_warning_reason",
        ],
    )

    mapped = map_redocking_evidence(labels, redocking)

    assert mapped["RMSD_threshold_call"].tolist() == [
        "RMSD <= 2 A",
        "RMSD <= 2 A",
        "2 A < RMSD <= 3 A",
        "RMSD > 3 A",
    ]
    assert mapped["Reference_pose_recovered"].tolist() == [True, True, False, False]
    assert mapped["Failure_or_warning_reason"].tolist() == ["zero", "none", "review", "warning"]


@pytest.mark.parametrize(("rmsd", "message"), [(-0.2, "nonnegative"), (np.inf, "finite")])
def test_redocking_rejects_negative_and_nonfinite_rmsd(rmsd: float, message: str):
    labels = pd.DataFrame(
        [["TEST", "TEST_V2C001", 1.0, "R_auto"]],
        columns=["pdb_id", "cluster_v2_id", "Within_PDB_rank", "automated_evidence_label"],
    )
    redocking = pd.DataFrame(
        [["TEST", "tool", -1.0, rmsd, "untrusted", "review"]],
        columns=[
            "pdb_id", "Docking_software", "GlideScore_kcal_per_mol",
            "Raw_ligand_RMSD_A", "RMSD_threshold_call", "Failure_or_warning_reason",
        ],
    )

    with pytest.raises(ValueError, match=message):
        map_redocking_evidence(labels, redocking)


def test_representative_case_ids_are_loaded_from_configuration():
    config = io.load_manuscript_config(ROOT / "config" / "manuscript.yaml")

    assert representative_case_ids(config) == ("5J89", "5TBM", "4W9H")


def test_importing_posthoc_creates_no_files(tmp_path: Path):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(ROOT / "src")

    completed = subprocess.run(
        [sys.executable, "-c", "import oips_repro.posthoc"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert list(tmp_path.iterdir()) == []
