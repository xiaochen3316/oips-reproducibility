from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path

import pandas as pd
import pytest

from oips_repro import io
from oips_repro.cli import main


ROOT = Path(__file__).parents[2]
EXPECTED = json.loads(
    (Path(__file__).parent / "data" / "expected_summary.json").read_text(encoding="utf-8")
)
def _contains_posthoc_semantics(columns: list[str]) -> bool:
    normalized = [column.lower().replace("-", "_") for column in columns]
    return any(
        value == "dcc"
        or value.startswith(("reference_", "dcc_", "ligand_contact_", "redocking_", "md_", "literature_"))
        for value in normalized
    )


def test_frozen_reference_csv_hashes_match_release_manifest():
    manifest = json.loads((ROOT / "release" / "manifest.json").read_text(encoding="utf-8"))
    declared = {
        entry["path"]: entry["sha256"]
        for entry in manifest["publication_payload"]
        if entry["path"].startswith("results/reference/")
        and entry["path"].endswith(".csv")
    }
    actual = {
        path.relative_to(ROOT).as_posix(): __import__("hashlib").sha256(
            path.read_bytes()
        ).hexdigest()
        for path in (ROOT / "results" / "reference").rglob("*.csv")
    }
    assert declared == actual


def test_frozen_cluster_v2_and_static_score_snapshot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    config = ROOT / "config" / "manuscript.yaml"
    cluster_dir = tmp_path / "cluster"
    score_dir = tmp_path / "score"
    monkeypatch.chdir(tmp_path)

    assert main(["cluster", "--config", str(config), "--output", str(cluster_dir)]) == 0
    assert main(
        [
            "score",
            "--config",
            str(config),
            "--cluster-dir",
            str(cluster_dir),
            "--output",
            str(score_dir),
        ]
    ) == 0

    assert {path.name for path in cluster_dir.iterdir()} == set(io.CLUSTER_FILE_SCHEMAS)
    assert {path.name for path in score_dir.iterdir()} == set(io.SCORE_FILE_SCHEMAS)

    candidates = pd.read_csv(cluster_dir / "cluster_v2_candidates.csv")
    membership = pd.read_csv(cluster_dir / "cluster_v2_membership.csv")
    mapping = pd.read_csv(cluster_dir / "tool_record_to_cluster_v2_mapping.csv")
    excluded = pd.read_csv(cluster_dir / "excluded_unmappable_records.csv")
    boundary = pd.read_csv(cluster_dir / "cluster_v2_boundary_audit.csv")
    master = pd.read_csv(score_dir / "cluster_v2_master_table.csv")
    rankings = pd.read_csv(score_dir / "cluster_v2_static_rankings.csv")

    assert len(mapping) == EXPECTED["records"]
    assert mapping["cluster_v2_id"].notna().sum() == EXPECTED["mapped"]
    assert len(excluded) == EXPECTED["excluded"]
    assert membership["same_tool_unit_id"].nunique() == EXPECTED["same_tool_units"]
    assert len(candidates) == EXPECTED["clusters"]
    assert candidates["pdb_id"].nunique() == EXPECTED["targets"]
    assert candidates["boundary_sensitive"].sum() == EXPECTED["boundary_sensitive"]
    assert len(boundary) == EXPECTED["clusters"]

    mappability = candidates["mappability"].value_counts().to_dict()
    cluster_expected = EXPECTED["cluster_distributions"]
    assert mappability == cluster_expected["mappability"]
    tool_support = {
        str(key): value
        for key, value in candidates["tool_support_count"].value_counts().sort_index().to_dict().items()
    }
    assert tool_support == cluster_expected["tool_support_count"]
    assert int(candidates["spatial_continuity"].sum()) == cluster_expected["spatial_continuity_true"]
    maximum = candidates.loc[candidates["cluster_diameter_A"].idxmax()]
    assert maximum["cluster_v2_id"] == cluster_expected["maximum_diameter_cluster_id"]
    diameter = EXPECTED["cluster_numeric"]["maximum_diameter_A"]
    assert maximum["cluster_diameter_A"] == pytest.approx(
        diameter["value"], abs=diameter["abs_tolerance"]
    )

    votes = membership.groupby(["cluster_v2_id", "tool"])["formal_vote_count"].sum()
    assert votes.max() == EXPECTED["maximum_formal_votes_per_cluster_tool"]
    assert membership.loc[membership["formal_tool_representative"]].groupby(
        ["cluster_v2_id", "tool"]
    ).size().max() == 1

    assert len(master) == EXPECTED["clusters"]
    assert master["P_lig"].isna().sum() == cluster_expected["missing_ligandability"]
    assert master["OIPS-P_static"].map(math.isfinite).all()
    center_only_master = master.loc[master["mappability"].eq("center_only_mappable")]
    center_only_rankings = rankings.loc[
        rankings["mappability"].eq("center_only_mappable")
    ]
    assert len(center_only_master) == cluster_expected["center_only_clusters"]
    assert len(center_only_rankings) == cluster_expected["center_only_clusters"]
    assert center_only_master["cluster_chain_entropy"].eq(0.0).all()
    assert center_only_rankings["cluster_chain_entropy"].eq(0.0).all()
    exemplar_id = EXPECTED["exemplar_cluster_id"]
    exemplar = master.set_index("cluster_v2_id").loc[exemplar_id]
    for name, spec in EXPECTED["exemplar_static_metrics"].items():
        if name in master:
            assert exemplar[name] == pytest.approx(spec["value"], abs=spec["abs_tolerance"])

    assert len(rankings) == EXPECTED["clusters"]
    ranked_exemplar = rankings.set_index("cluster_v2_id").loc[exemplar_id]
    rank_spec = EXPECTED["exemplar_static_metrics"]["Within_PDB_rank"]
    assert ranked_exemplar["Within_PDB_rank"] == pytest.approx(
        rank_spec["value"], abs=rank_spec["abs_tolerance"],
    )
    assert bool(ranked_exemplar["tie_flag"]) == EXPECTED["exemplar_static_state"]["tie_flag"]
    assert rankings["tie_flag"].sum() == cluster_expected["tie_count"]
    assert rankings["tie_size"].eq(1).all()
    assert rankings["Within_PDB_rank"].map(lambda value: float(value).is_integer()).all()
    recomputation = EXPECTED["cluster_numeric"]["maximum_recomputation_difference"]
    assert (
        rankings["OIPS-P_static"] - rankings["OIPS-P_static_recomputed"]
    ).abs().max() == pytest.approx(recomputation["value"], abs=recomputation["abs_tolerance"])

    for filename, (columns, sort_by) in io.CLUSTER_FILE_SCHEMAS.items():
        frame = pd.read_csv(cluster_dir / filename)
        assert frame.columns.tolist() == list(columns)
        assert frame.reset_index(drop=True).equals(
            frame.sort_values(list(sort_by), kind="mergesort", na_position="last").reset_index(drop=True)
        )
        assert not _contains_posthoc_semantics(frame.columns.tolist())
    for filename, (columns, sort_by) in io.SCORE_FILE_SCHEMAS.items():
        frame = pd.read_csv(score_dir / filename)
        assert frame.columns.tolist() == list(columns)
        assert frame.reset_index(drop=True).equals(
            frame.sort_values(list(sort_by), kind="mergesort", na_position="last").reset_index(drop=True)
        )
        assert not _contains_posthoc_semantics(frame.columns.tolist())


def _build_cluster_result(tmp_path: Path, capsys):
    config = ROOT / "config" / "manuscript.yaml"
    cluster_dir = tmp_path / "original-cluster"
    assert main(
        ["cluster", "--config", str(config), "--output", str(cluster_dir)]
    ) == 0
    capsys.readouterr()
    features = io.load_feature_table(ROOT / "data" / "static" / "tool_pocket_features.csv")
    result = io.load_clustering_result(cluster_dir, features=features)
    return config, features, result


def test_score_cli_rejects_a_coherently_truncated_cluster_handoff(
    tmp_path: Path, capsys
):
    config, _, result = _build_cluster_result(tmp_path, capsys)
    singleton = result.candidates.loc[result.candidates["raw_record_count"].eq(1)].iloc[0]
    cluster_id = singleton["cluster_v2_id"]
    tampered = replace(
        result,
        candidates=result.candidates.loc[
            ~result.candidates["cluster_v2_id"].eq(cluster_id)
        ].reset_index(drop=True),
        membership=result.membership.loc[
            ~result.membership["cluster_v2_id"].eq(cluster_id)
        ].reset_index(drop=True),
        mapping=result.mapping.loc[
            ~result.mapping["cluster_v2_id"].eq(cluster_id)
        ].reset_index(drop=True),
        boundary=result.boundary.loc[
            ~result.boundary["cluster_v2_id"].eq(cluster_id)
        ].reset_index(drop=True),
    )
    tampered_dir = tmp_path / "truncated-cluster"
    io.write_clustering_result(tampered, tampered_dir)
    assert len(tampered.candidates) == len(result.candidates) - 1
    assert len(tampered.mapping) == len(result.mapping) - 1

    exit_code = main(
        [
            "score", "--config", str(config), "--cluster-dir", str(tampered_dir),
            "--output", str(tmp_path / "truncated-score"),
        ]
    )
    captured = capsys.readouterr()

    assert exit_code != 0
    assert "mapping row IDs do not match" in captured.err
    assert "Traceback" not in captured.err


def test_score_cli_rejects_mapping_identity_tampering(tmp_path: Path, capsys):
    config, _, result = _build_cluster_result(tmp_path, capsys)
    mapping = result.mapping.copy()
    mapping.loc[mapping.index[0], "pocket_id"] = "BOGUS_POCKET"
    tampered = replace(result, mapping=mapping)
    tampered_dir = tmp_path / "identity-tampered-cluster"
    io.write_clustering_result(tampered, tampered_dir)

    exit_code = main(
        [
            "score", "--config", str(config), "--cluster-dir", str(tampered_dir),
            "--output", str(tmp_path / "identity-tampered-score"),
        ]
    )
    captured = capsys.readouterr()

    assert exit_code != 0
    assert "mapping identity mismatch" in captured.err
    assert "BOGUS_POCKET" not in captured.err
    assert "Traceback" not in captured.err
