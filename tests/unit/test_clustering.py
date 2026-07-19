from __future__ import annotations

from dataclasses import FrozenInstanceError
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).parents[2]
PUBLIC_FEATURES = ROOT / "data" / "static" / "tool_pocket_features.csv"


def load_modules():
    try:
        from oips_repro import clustering, io, scoring, structure
    except ImportError as exc:
        pytest.fail(f"public cluster-v2 API is missing: {exc}")
    return clustering, io, scoring, structure


def make_feature(
    row_id: int,
    *,
    pdb_id: str = "TEST",
    tool: str = "CavityPlus",
    pocket_id: str | None = None,
    center: tuple[float, float, float] | None = (0.0, 0.0, 0.0),
    residues: tuple[str, ...] = ("A:ALA:1:",),
    display_order: float | None = 1.0,
    sitemap_rank: float | None = None,
    geometry: float | None = 50.0,
    ligandability: float | None = 50.0,
) -> dict[str, object]:
    x, y, z = center if center is not None else (np.nan, np.nan, np.nan)
    return {
        "row_id": row_id,
        "pdb_id": pdb_id,
        "tool": tool,
        "pocket_id": pocket_id or f"P{row_id}",
        "display_order": display_order,
        "sitemap_rank": sitemap_rank,
        "center_x": x,
        "center_y": y,
        "center_z": z,
        "center_method": "synthetic" if center is not None else "missing",
        "residue_count": len(residues),
        "residue_set_json": json.dumps(list(residues)),
        "pocket_geometry_score": geometry,
        "pocket_ligandability_score": ligandability,
    }


def make_unit(
    row_id: int,
    *,
    tool: str = "CavityPlus",
    center: tuple[float, float, float] | None = (0.0, 0.0, 0.0),
    residues: frozenset[str] = frozenset({"A:ALA:1:"}),
    native_rank: float = 1.0,
):
    clustering, _, _, _ = load_modules()
    return clustering.Unit(
        pdb_id="TEST",
        tool=tool,
        unit_id=f"TEST_{tool}_U{row_id:03d}",
        representative_row_id=row_id,
        representative_pocket_id=f"P{row_id}",
        raw_row_ids=(row_id,),
        raw_pocket_ids=(f"P{row_id}",),
        center=center,
        residues=residues,
        geometry=50.0,
        ligandability=50.0,
        native_rank=native_rank,
    )


def test_configuration_and_result_types_are_immutable():
    clustering, _, scoring, _ = load_modules()

    same = clustering.SameToolConfig()
    cross = clustering.CrossToolConfig()
    score = scoring.ScoringConfig()

    assert same.residue_iou_strong == 0.70
    assert same.dogsite_hierarchy_containment_min == 0.70
    assert same.representative_missing_center_penalty == 0.20
    assert cross.cluster_center_diameter_max_A == 12.0
    assert score.geometry_top_n == 3
    assert score.ligandability_top_n == 2
    with pytest.raises(FrozenInstanceError):
        same.residue_iou_strong = 0.5
    with pytest.raises(FrozenInstanceError):
        cross.close_center_max_A = 5.0


def test_residue_parsing_and_simple_conversion_are_canonical():
    clustering, _, _, _ = load_modules()
    residues = {"A:ALA:10:", "B:GLY:11:A"}

    assert clustering.simple_residue("A:ALA:10:") == "A:10:"
    assert clustering.simple_residue("B:GLY:11:A") == "B:11:A"
    assert clustering.parse_residues(json.dumps(sorted(residues))) == residues
    assert clustering.parse_residues("A:ALA:10:;B:GLY:11:A") == residues
    assert clustering.parse_residues(None) == set()


def test_residue_overlap_metrics_simplify_names_and_keep_empty_nan():
    clustering, _, _, _ = load_modules()
    a = {"A:ALA:1:", "A:SER:2:", "B:GLY:3:"}
    b = {"A:VAL:1:", "B:GLY:3:", "C:TYR:4:"}

    assert clustering.residue_iou(a, b) == pytest.approx(2 / 4)
    assert clustering.residue_containment(a, b) == pytest.approx(2 / 3)
    assert math.isnan(clustering.residue_iou(set(), b))
    assert math.isnan(clustering.residue_containment(a, set()))


def test_same_tool_grouping_is_complete_link_not_single_link():
    clustering, _, _, _ = load_modules()
    features = pd.DataFrame(
        [
            make_feature(1, center=(0.0, 0.0, 0.0), residues=()),
            make_feature(2, center=(1.4, 0.0, 0.0), residues=()),
            make_feature(3, center=(2.8, 0.0, 0.0), residues=()),
        ]
    )

    result = clustering.build_clusters(
        features,
        same_tool=clustering.SameToolConfig(),
        cross_tool=clustering.CrossToolConfig(),
        tool_weights={"CavityPlus": 1.0},
    )

    assert len(result.units) == 2
    assert sorted(unit.raw_row_ids for unit in result.units) == [(1, 2), (3,)]


def _raw_record(
    *,
    tool: str = "CavityPlus",
    pocket_id: str = "P_1",
    center: tuple[float, float, float] | None = (0.0, 0.0, 0.0),
    residues: set[str] | None = None,
) -> dict[str, object]:
    return {
        "row_id": 1,
        "tool": tool,
        "pocket_id": pocket_id,
        "center": center,
        "residues": residues if residues is not None else {"A:ALA:1:"},
        "display_order": 1,
        "sitemap_rank": np.nan,
    }


def _overlap_sets(intersection: int, a_only: int, b_only: int) -> tuple[set[str], set[str]]:
    common = {f"A:ALA:{index}:" for index in range(intersection)}
    a = common | {f"A:SER:{100 + index}:" for index in range(a_only)}
    b = common | {f"A:GLY:{200 + index}:" for index in range(b_only)}
    return a, b


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (
            _raw_record(tool="DoGSite3", pocket_id="P_1", center=(0.0, 0.0, 0.0)),
            _raw_record(tool="DoGSite3", pocket_id="P_1_1", center=(6.0, 0.0, 0.0)),
        ),
        (
            _raw_record(tool="DoGSite3", pocket_id="P_1", center=None, residues=_overlap_sets(7, 6, 7)[0]),
            _raw_record(tool="DoGSite3", pocket_id="P_1_1", center=None, residues=_overlap_sets(7, 6, 7)[1]),
        ),
        (
            _raw_record(tool="DoGSite3", pocket_id="P_1", center=None, residues=_overlap_sets(7, 3, 20)[0]),
            _raw_record(tool="DoGSite3", pocket_id="P_1_1", center=None, residues=_overlap_sets(7, 3, 20)[1]),
        ),
        (
            _raw_record(center=None, residues=_overlap_sets(7, 1, 2)[0]),
            _raw_record(center=None, residues=_overlap_sets(7, 1, 2)[1]),
        ),
        (
            _raw_record(center=(0.0, 0.0, 0.0), residues=_overlap_sets(17, 3, 20)[0]),
            _raw_record(center=(8.0, 0.0, 0.0), residues=_overlap_sets(17, 3, 20)[1]),
        ),
        (
            _raw_record(center=(0.0, 0.0, 0.0), residues=_overlap_sets(2, 3, 3)[0]),
            _raw_record(center=(2.5, 0.0, 0.0), residues=_overlap_sets(2, 3, 3)[1]),
        ),
        (
            _raw_record(center=(0.0, 0.0, 0.0), residues=set()),
            _raw_record(center=(1.5, 0.0, 0.0), residues=set()),
        ),
    ],
)
def test_every_same_tool_threshold_is_inclusive(left, right):
    clustering, _, _, _ = load_modules()

    assert clustering.same_tool_pair_is_duplicate(
        left, right, clustering.SameToolConfig()
    )


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (
            _raw_record(tool="DoGSite3", pocket_id="P_1", center=(0.0, 0.0, 0.0), residues=set()),
            _raw_record(tool="DoGSite3", pocket_id="P_1_1", center=(6.000001, 0.0, 0.0), residues=set()),
        ),
        (
            _raw_record(center=(0.0, 0.0, 0.0), residues=_overlap_sets(17, 3, 20)[0]),
            _raw_record(center=(8.000001, 0.0, 0.0), residues=_overlap_sets(17, 3, 20)[1]),
        ),
        (
            _raw_record(center=(0.0, 0.0, 0.0), residues=_overlap_sets(2, 3, 3)[0]),
            _raw_record(center=(2.500001, 0.0, 0.0), residues=_overlap_sets(2, 3, 3)[1]),
        ),
        (
            _raw_record(center=(0.0, 0.0, 0.0), residues=set()),
            _raw_record(center=(1.500001, 0.0, 0.0), residues=set()),
        ),
    ],
)
def test_same_tool_distance_thresholds_reject_just_beyond(left, right):
    clustering, _, _, _ = load_modules()

    assert not clustering.same_tool_pair_is_duplicate(
        left, right, clustering.SameToolConfig()
    )


def test_cross_tool_compatibility_thresholds_are_inclusive_and_bounded():
    clustering, _, _, _ = load_modules()
    config = clustering.CrossToolConfig()
    empty = frozenset()
    iou_20 = tuple(frozenset(x) for x in _overlap_sets(2, 4, 4))
    iou_35 = tuple(frozenset(x) for x in _overlap_sets(7, 6, 7))

    assert clustering.site_pair_compatible(
        make_unit(1, tool="A", center=(0.0, 0.0, 0.0), residues=empty),
        make_unit(2, tool="B", center=(6.0, 0.0, 0.0), residues=empty),
        config,
    )
    assert not clustering.site_pair_compatible(
        make_unit(1, tool="A", center=(0.0, 0.0, 0.0), residues=empty),
        make_unit(2, tool="B", center=(6.000001, 0.0, 0.0), residues=empty),
        config,
    )
    assert clustering.site_pair_compatible(
        make_unit(1, tool="A", center=(0.0, 0.0, 0.0), residues=iou_20[0]),
        make_unit(2, tool="B", center=(10.0, 0.0, 0.0), residues=iou_20[1]),
        config,
    )
    assert not clustering.site_pair_compatible(
        make_unit(1, tool="A", center=(0.0, 0.0, 0.0), residues=iou_20[0]),
        make_unit(2, tool="B", center=(10.000001, 0.0, 0.0), residues=iou_20[1]),
        config,
    )
    assert clustering.site_pair_compatible(
        make_unit(1, tool="A", center=None, residues=iou_35[0]),
        make_unit(2, tool="B", center=None, residues=iou_35[1]),
        config,
    )


def test_global_diameter_cap_is_inclusive_and_rejects_just_beyond():
    clustering, _, _, _ = load_modules()
    config = clustering.CrossToolConfig()
    at_cap = [
        make_unit(1, tool="A", center=(0.0, 0.0, 0.0), residues=frozenset()),
        make_unit(2, tool="B", center=(6.0, 0.0, 0.0), residues=frozenset()),
        make_unit(3, tool="C", center=(12.0, 0.0, 0.0), residues=frozenset()),
    ]
    beyond = [at_cap[0], at_cap[1], make_unit(3, tool="C", center=(12.000001, 0.0, 0.0), residues=frozenset())]

    assert clustering.medoid_constrained_group_is_valid(at_cap, config)
    assert not clustering.medoid_constrained_group_is_valid(beyond, config)


def test_medoid_and_formal_representative_ties_use_stable_row_ids():
    clustering, _, _, _ = load_modules()
    units = [
        make_unit(3, tool="A", center=(0.0, 0.0, 0.0), residues=frozenset(), native_rank=1.0),
        make_unit(1, tool="A", center=(5.0, 0.0, 0.0), residues=frozenset(), native_rank=1.0),
        make_unit(2, tool="B", center=(2.5, 0.0, 0.0), residues=frozenset(), native_rank=1.0),
    ]

    assert clustering.medoid_unit(units).representative_row_id == 2
    representatives = clustering.formal_representatives(units)
    assert [(unit.tool, unit.representative_row_id) for unit in representatives] == [
        ("A", 1),
        ("B", 2),
    ]


def test_equal_cost_medoid_candidates_choose_lower_representative_row_id():
    clustering, _, _, _ = load_modules()
    units = [
        make_unit(9, tool="A", center=(0.0, 0.0, 0.0), residues=frozenset()),
        make_unit(2, tool="B", center=(10.0, 0.0, 0.0), residues=frozenset()),
    ]

    assert clustering.medoid_unit(units).representative_row_id == 2


def test_cluster_ids_and_formal_votes_are_deterministic():
    clustering, _, _, _ = load_modules()
    features = pd.DataFrame(
        [
            make_feature(3, tool="A", center=(0.0, 0.0, 0.0), residues=(), display_order=1),
            make_feature(1, tool="A", center=(5.0, 0.0, 0.0), residues=(), display_order=1),
            make_feature(2, tool="B", center=(2.5, 0.0, 0.0), residues=(), display_order=1),
        ]
    )

    result = clustering.build_clusters(
        features,
        same_tool=clustering.SameToolConfig(),
        cross_tool=clustering.CrossToolConfig(),
        tool_weights={"A": 1.0, "B": 1.0},
    )

    assert result.candidates["cluster_v2_id"].tolist() == ["TEST_V2C001"]
    votes = result.membership.groupby(["cluster_v2_id", "tool"])["formal_vote_count"].sum()
    assert votes.max() == 1
    formal = result.membership.loc[result.membership["formal_tool_representative"]]
    assert formal.groupby(["cluster_v2_id", "tool"]).size().max() == 1
    assert formal.loc[formal["tool"].eq("A"), "raw_row_id"].item() == 1


def test_residue_consensus_uses_formal_representatives():
    clustering, _, _, _ = load_modules()
    units = [
        make_unit(1, tool="A", residues=frozenset({"A:ALA:1:", "A:SER:2:"})),
        make_unit(2, tool="B", residues=frozenset({"A:VAL:1:", "B:GLY:3:"})),
        make_unit(3, tool="C", residues=frozenset({"A:THR:1:", "B:GLY:3:", "C:TYR:4:"})),
    ]

    core, envelope, support = clustering.residue_consensus(units)

    assert {clustering.simple_residue(value) for value in core} == {"A:1:", "B:3:"}
    assert {clustering.simple_residue(value) for value in envelope} == {
        "A:1:", "A:2:", "B:3:", "C:4:"
    }
    assert support["A:1:"] == 3


def test_boundary_rules_are_strict_and_require_two_formal_tools():
    clustering, _, _, _ = load_modules()
    thresholds = dict(
        diameter_gt_A=9.0,
        dispersion_gt_A=4.0,
        core_envelope_ratio_lt=0.20,
        median_residue_iou_lt=0.15,
    )

    exact = clustering.boundary_decision(
        tool_support_count=2,
        diameter=9.0,
        dispersion=4.0,
        core_ratio=0.20,
        median_iou=0.15,
        **thresholds,
    )
    over = clustering.boundary_decision(
        tool_support_count=2,
        diameter=9.000001,
        dispersion=4.0,
        core_ratio=0.20,
        median_iou=0.15,
        **thresholds,
    )
    single = clustering.boundary_decision(
        tool_support_count=1,
        diameter=99.0,
        dispersion=99.0,
        core_ratio=0.0,
        median_iou=0.0,
        **thresholds,
    )

    assert exact == (False, "")
    assert over == (True, "diameter_gt_9A")
    assert single == (False, "")


def test_missing_aware_cluster_modules_and_oips_are_renormalized():
    clustering, _, scoring, structure = load_modules()
    features = pd.DataFrame(
        [
            make_feature(
                1,
                geometry=90.0,
                ligandability=None,
                residues=("A:ALA:1:",),
            )
        ]
    )
    clustered = clustering.build_clusters(
        features,
        same_tool=clustering.SameToolConfig(),
        cross_tool=clustering.CrossToolConfig(),
        tool_weights={"CavityPlus": 1.0, "Other": 1.0},
    )
    interface = structure.InterfaceProfile(
        chain_count=1,
        protein_residue_count=1,
        interface_residues=(),
        interface_atom_coordinates=(),
        pair_contact_counts={},
    )
    config = scoring.ScoringConfig(
        tool_weights={"CavityPlus": 1.0, "Other": 1.0},
    )

    result = scoring.score_clusters(clustered, interfaces={"TEST": interface}, config=config)
    row = result.master.iloc[0]

    assert row["C_cons"] == 25.0
    assert row["G_geo"] == 90.0
    assert math.isnan(row["P_lig"])
    assert row["Q_evidence"] == 71.0
    assert row["O_rel_formal"] == 35.0
    expected = scoring.weighted_module_score(
        {"C_cons": 25.0, "G_geo": 90.0, "P_lig": np.nan, "O_rel_formal": 35.0, "Q_evidence": 71.0},
        config.module_weights,
    )
    assert row["OIPS-P_static"] == pytest.approx(expected)


def test_average_ranking_preserves_fractional_ties_and_stable_sort():
    _, _, scoring, _ = load_modules()
    master = pd.DataFrame(
        {
            "pdb_id": ["B", "A", "A", "A"],
            "cluster_v2_id": ["B_V2C001", "A_V2C002", "A_V2C001", "A_V2C003"],
            "OIPS-P_static": [7.0, 10.0, 10.0, 5.0],
        }
    )

    ranked = scoring.rank_within_target(master)
    target = ranked.loc[ranked["pdb_id"].eq("A")]

    assert target["cluster_v2_id"].tolist() == ["A_V2C001", "A_V2C002", "A_V2C003"]
    assert target["Within_PDB_rank"].tolist() == [1.5, 1.5, 3.0]
    assert target["tie_flag"].tolist() == [True, True, False]
    assert target["tie_size"].tolist() == [2, 2, 1]


def test_low_level_clustering_ignores_unrelated_in_memory_columns():
    clustering, _, _, _ = load_modules()
    features = pd.DataFrame(
        [
            make_feature(2, tool="A", center=(4.0, 0.0, 0.0)),
            make_feature(1, tool="B", center=(0.0, 0.0, 0.0)),
        ]
    )
    kwargs = dict(
        same_tool=clustering.SameToolConfig(),
        cross_tool=clustering.CrossToolConfig(),
        tool_weights={"A": 1.0, "B": 1.0},
    )

    baseline = clustering.build_clusters(features, **kwargs)
    extended = clustering.build_clusters(features.assign(unrelated_note=["x", "y"]), **kwargs)

    pd.testing.assert_frame_equal(baseline.candidates, extended.candidates)
    pd.testing.assert_frame_equal(baseline.membership, extended.membership)
    pd.testing.assert_frame_equal(baseline.mapping, extended.mapping)


def test_shuffled_input_produces_byte_equivalent_cluster_outputs(tmp_path: Path):
    clustering, io, _, _ = load_modules()
    features = pd.DataFrame(
        [
            make_feature(4, tool="A", center=(20.0, 0.0, 0.0), residues=("A:ALA:4:",)),
            make_feature(2, tool="A", center=(0.0, 0.0, 0.0), residues=("A:ALA:1:",)),
            make_feature(3, tool="B", center=(4.0, 0.0, 0.0), residues=("A:VAL:1:",)),
            make_feature(1, tool="B", center=(20.0, 0.0, 0.0), residues=("A:ALA:4:",)),
        ]
    )
    kwargs = dict(
        same_tool=clustering.SameToolConfig(),
        cross_tool=clustering.CrossToolConfig(),
        tool_weights={"A": 1.0, "B": 1.0},
    )
    first = clustering.build_clusters(features, **kwargs)
    second = clustering.build_clusters(features.sample(frac=1.0, random_state=17), **kwargs)
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"

    io.write_clustering_result(first, first_dir)
    io.write_clustering_result(second, second_dir)

    assert {
        path.name: path.read_bytes() for path in first_dir.iterdir()
    } == {
        path.name: path.read_bytes() for path in second_dir.iterdir()
    }


@pytest.mark.parametrize(
    "forbidden_header",
    [
        "reference_ligand",
        "DCC_A",
        "ligand_contact_recovery",
        "redocking_score",
        "md_support",
        "literature_annotation",
    ],
)
def test_public_loader_rejects_every_posthoc_header_class(
    tmp_path: Path, forbidden_header: str
):
    _, io, _, _ = load_modules()
    frame = pd.DataFrame([make_feature(1)])
    frame[forbidden_header] = "leak"
    path = tmp_path / "forbidden.csv"
    frame.to_csv(path, index=False)

    with pytest.raises(ValueError, match="post-hoc|forbidden"):
        io.load_feature_table(path)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda frame: frame.__setitem__("row_id", [1] * len(frame)), "row_id"),
        (lambda frame: frame.__setitem__("pdb_id", [""] + frame["pdb_id"].tolist()[1:]), "pdb_id"),
        (lambda frame: frame.__setitem__("center_x", ["bad"] + frame["center_x"].tolist()[1:]), "numeric"),
        (lambda frame: frame.__setitem__("residue_count", [-1] + frame["residue_count"].tolist()[1:]), "residue_count"),
        (lambda frame: frame.__setitem__("residue_set_json", ["{}"] + frame["residue_set_json"].tolist()[1:]), "JSON array"),
        (lambda frame: frame.loc.__setitem__((0, "center_x"), ""), "center"),
        (lambda frame: frame.loc.__setitem__((0, "residue_count"), int(frame.loc[0, "residue_count"]) + 1), "residue_count"),
    ],
)
def test_public_loader_validates_the_frozen_data_contract(tmp_path: Path, mutation, message):
    _, io, _, _ = load_modules()
    frame = pd.read_csv(PUBLIC_FEATURES, keep_default_na=False)
    mutation(frame)
    path = tmp_path / "invalid.csv"
    frame.to_csv(path, index=False)

    with pytest.raises(ValueError, match=message):
        io.load_feature_table(path)


def test_public_loader_preserves_schema_directed_missingness():
    _, io, _, _ = load_modules()

    features = io.load_feature_table(PUBLIC_FEATURES)

    assert list(features.columns) == list(io.FEATURE_COLUMNS)
    assert len(features) > 0
    assert features["row_id"].dtype.kind in "iu"
    assert features["pdb_id"].map(type).eq(str).all()
    for column in ("center_x", "sitemap_rank", "pocket_ligandability_score"):
        assert features[column].isna().any()
        assert features[column].notna().any()
    assert not features["pdb_id"].eq("NaN").any()


def test_public_loader_accepts_arbitrary_positive_sorted_row_count(tmp_path: Path):
    _, io, _, _ = load_modules()
    rows = pd.DataFrame([
        make_feature(1, pdb_id="4W9H", tool="CavityPlus", pocket_id="P_1"),
        make_feature(2, pdb_id="5J89", tool="DoGSite3", pocket_id="P_1"),
    ], columns=io.FEATURE_COLUMNS)
    path = tmp_path / "features.csv"
    rows.to_csv(path, index=False)

    loaded = io.load_feature_table(path)

    assert len(loaded) == 2
    assert loaded["row_id"].tolist() == [1, 2]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda frame: frame.iloc[[1, 0]].reset_index(drop=True), "sorted"),
        (lambda frame: frame.assign(center_x=[float("inf"), 0.0]), "finite"),
        (lambda frame: frame.assign(row_id=[1, 1]), "row_id"),
    ],
)
def test_public_loader_rejects_unstable_or_nonfinite_small_input(
    tmp_path: Path, mutation, message
):
    _, io, _, _ = load_modules()
    rows = pd.DataFrame([
        make_feature(1, pdb_id="4W9H", tool="CavityPlus", pocket_id="P_1"),
        make_feature(2, pdb_id="5J89", tool="DoGSite3", pocket_id="P_1"),
    ], columns=io.FEATURE_COLUMNS)
    path = tmp_path / "features.csv"
    mutation(rows).to_csv(path, index=False)

    with pytest.raises(ValueError, match=message):
        io.load_feature_table(path)


def test_cluster_output_schemas_sorts_booleans_and_readback_missingness(tmp_path: Path):
    clustering, io, _, _ = load_modules()
    features = pd.DataFrame(
        [make_feature(1, center=None, residues=("A:ALA:1:",), ligandability=None)]
    )
    result = clustering.build_clusters(
        features,
        same_tool=clustering.SameToolConfig(),
        cross_tool=clustering.CrossToolConfig(),
        tool_weights={"CavityPlus": 1.0},
    )
    result = clustering.annotate_boundaries(
        result,
        diameter_gt_A=9.0,
        dispersion_gt_A=4.0,
        core_envelope_ratio_lt=0.20,
        median_residue_iou_lt=0.15,
    )

    io.write_clustering_result(result, tmp_path)

    expected_files = set(io.CLUSTER_FILE_SCHEMAS)
    assert {path.name for path in tmp_path.iterdir()} == expected_files
    for filename, (columns, sort_by) in io.CLUSTER_FILE_SCHEMAS.items():
        raw = pd.read_csv(tmp_path / filename, keep_default_na=False)
        assert raw.columns.tolist() == list(columns)
        assert raw.reset_index(drop=True).equals(
            raw.sort_values(list(sort_by), kind="mergesort", na_position="last").reset_index(drop=True)
        )
    candidate_bytes = (tmp_path / "cluster_v2_candidates.csv").read_bytes()
    assert b",true," in candidate_bytes or b",false," in candidate_bytes
    assert b"True" not in candidate_bytes and b"False" not in candidate_bytes

    roundtrip = pd.read_csv(tmp_path / "cluster_v2_candidates.csv")
    assert math.isnan(roundtrip.loc[0, "medoid_center_x"])
    assert bool(roundtrip.loc[0, "spatial_continuity"]) is True
    assert bool(roundtrip.loc[0, "boundary_sensitive"]) is False


def test_cluster_loader_rejects_invalid_candidate_enum(tmp_path: Path):
    clustering, io, _, _ = load_modules()
    features = pd.DataFrame([make_feature(1)])
    result = clustering.build_clusters(
        features,
        same_tool=clustering.SameToolConfig(),
        cross_tool=clustering.CrossToolConfig(),
        tool_weights={"CavityPlus": 1.0},
    )
    result = clustering.annotate_boundaries(
        result, diameter_gt_A=9.0, dispersion_gt_A=4.0,
        core_envelope_ratio_lt=0.2, median_residue_iou_lt=0.15,
    )
    io.write_clustering_result(result, tmp_path)
    candidate = pd.read_csv(tmp_path / "cluster_v2_candidates.csv", dtype=str, keep_default_na=False)
    candidate.loc[0, "mappability"] = "invalid_state"
    candidate.to_csv(tmp_path / "cluster_v2_candidates.csv", index=False, lineterminator="\n")

    with pytest.raises(ValueError, match="mappability.*invalid"):
        io.load_clustering_result(tmp_path, features=features)


def test_repository_root_resolution_and_reference_output_guard_are_cwd_independent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    _, io, _, _ = load_modules()
    project = tmp_path / "Project"
    config_dir = project / "config"
    config_dir.mkdir(parents=True)
    config_path = config_dir / "manuscript.yaml"
    source_config = Path(__file__).parents[2] / "config" / "manuscript.yaml"
    config_path.write_text(source_config.read_text(encoding="utf-8"), encoding="utf-8")
    source_contract = Path(__file__).parents[2] / "config" / "figure_contract.yaml"
    (config_dir / "figure_contract.yaml").write_bytes(source_contract.read_bytes())
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    loaded = io.load_manuscript_config(config_path)

    assert loaded.repository_root == project.resolve()
    assert loaded.resolve_configured_path("feature_table") == (
        project / "data" / "static" / "tool_pocket_features.csv"
    ).resolve()
    reference = project / "Results" / "Reference"
    with pytest.raises(ValueError, match="reference"):
        io.ensure_safe_output_path(project / "results" / "reference" / "nested", reference)
    assert io.ensure_safe_output_path(project / "results" / "reproduced", reference).name == "reproduced"


def test_cli_stage_arguments_are_explicit_and_posthoc_inputs_are_rejected(capsys):
    _, _, _, _ = load_modules()
    from oips_repro.cli import main

    assert main(["cluster", "--help"]) == 0
    cluster_help = capsys.readouterr().out
    assert "--config" in cluster_help and "--output" in cluster_help
    assert "--reference" not in cluster_help
    assert "--md" not in cluster_help

    assert main(["score", "--help"]) == 0
    score_help = capsys.readouterr().out
    assert "--config" in score_help
    assert "--cluster-dir" in score_help
    assert "--output" in score_help
    assert "--redocking" not in score_help
    assert "--literature" not in score_help

    assert main([
        "cluster", "--config", "config/manuscript.yaml", "--output", "out",
        "--reference", "secret.csv",
    ]) != 0
    assert "unrecognized arguments" in capsys.readouterr().err
