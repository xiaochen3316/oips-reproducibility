from dataclasses import FrozenInstanceError
import math

import pytest


def load_scoring_modules():
    try:
        from oips_repro import scoring, structure
    except ImportError as exc:
        pytest.fail(f"public scoring primitives are missing: {exc}")
    return scoring, structure


def make_interface(
    *,
    chain_count: int = 2,
    protein_residue_count: int = 100,
    interface_residues: tuple[str, ...] = ("A:ALA:1:", "B:GLY:2:"),
    interface_atom_coordinates: tuple[tuple[float, float, float], ...] = ((0.0, 0.0, 0.0),),
):
    _, structure = load_scoring_modules()
    return structure.InterfaceProfile(
        chain_count=chain_count,
        protein_residue_count=protein_residue_count,
        interface_residues=interface_residues,
        interface_atom_coordinates=interface_atom_coordinates,
        pair_contact_counts={"A-B": 1} if interface_residues else {},
    )


def test_scoring_config_freezes_manuscript_defaults():
    scoring, _ = load_scoring_modules()

    config = scoring.ScoringConfig()

    assert config.interface_fraction_weight == 0.42
    assert config.chain_context_weight == 0.25
    assert config.distance_weight == 0.20
    assert config.interface_extent_weight == 0.13
    assert config.distance_near_max_A == 4.0
    assert config.distance_intermediate_max_A == 8.0
    assert config.distance_far_max_A == 16.0
    assert config.distance_tail_floor == 10.0
    assert config.interface_fraction_base == 25.0
    assert config.interface_fraction_multiplier == 130.0
    assert config.multi_chain_base == 75.0
    assert config.chain_entropy_multiplier == 25.0
    assert config.single_chain_interface_fraction_min == 0.35
    assert config.single_chain_supported_score == 65.0
    assert config.single_chain_other_score == 30.0
    assert config.context_base == 45.0
    assert config.context_per_extra_chain == 10.0
    assert config.context_extra_chain_cap == 4
    assert config.context_interface_fraction_cap == 0.5
    assert config.context_interface_fraction_multiplier == 80.0
    assert config.no_interface_monomer_score == 35.0
    assert config.no_interface_multichain_score == 45.0
    with pytest.raises(FrozenInstanceError):
        config.distance_near_max_A = 5.0


def test_weighted_module_score_renormalizes_available_positive_weights():
    scoring, _ = load_scoring_modules()
    modules = {
        "C_cons": 100.0,
        "G_geo": 50.0,
        "P_lig": None,
        "O_rel_formal": None,
        "Q_evidence": None,
    }
    weights = {
        "C_cons": 0.22,
        "G_geo": 0.18,
        "P_lig": 0.24,
        "O_rel_formal": 0.24,
        "Q_evidence": 0.12,
    }

    assert scoring.weighted_module_score(modules, weights) == pytest.approx(77.5)


def test_weighted_module_score_ignores_none_nan_and_nonpositive_weights():
    scoring, _ = load_scoring_modules()

    assert scoring.weighted_module_score(
        {"valid": 80.0, "none": None, "nan": float("nan"), "zero": 0.0},
        {"valid": 2.0, "none": 3.0, "nan": 4.0, "zero": 0.0},
    ) == 80.0
    assert math.isnan(
        scoring.weighted_module_score(
            {"none": None, "nan": float("nan"), "negative": 90.0},
            {"none": 1.0, "nan": 1.0, "negative": -1.0},
        )
    )


@pytest.mark.parametrize(
    ("distance", "expected"),
    [
        (4.0, 100.0),
        (4.000001, 89.999995),
        (8.0, 70.0),
        (8.000001, 69.999994),
        (16.0, 22.0),
        (16.000001, 21.9999988),
        (26.0, 10.0),
        (100.0, 10.0),
    ],
)
def test_distance_score_breakpoints_are_inclusive_and_tail_is_floored(distance, expected):
    scoring, _ = load_scoring_modules()
    result = scoring.score_oligomer_relevance(
        {"A:ALA:1:"},
        (distance, 0.0, 0.0),
        make_interface(),
    )

    assert result["distance_to_interface_A"] == pytest.approx(distance)
    assert result["interface_distance_score"] == pytest.approx(expected)


def test_no_interface_profiles_use_frozen_monomer_and_multichain_bases():
    scoring, _ = load_scoring_modules()
    monomer = make_interface(
        chain_count=1,
        protein_residue_count=1,
        interface_residues=(),
        interface_atom_coordinates=(),
    )
    multichain = make_interface(
        chain_count=2,
        protein_residue_count=2,
        interface_residues=(),
        interface_atom_coordinates=(),
    )

    monomer_result = scoring.score_oligomer_relevance({"A:ALA:1:"}, None, monomer)
    multichain_result = scoring.score_oligomer_relevance({"A:ALA:1:"}, None, multichain)

    assert monomer_result["O_rel_formal"] == 35.0
    assert multichain_result["O_rel_formal"] == 45.0
    assert math.isnan(monomer_result["interface_recall"])
    assert math.isnan(multichain_result["interface_distance_score"])
    assert "context_score" in monomer_result
    assert "context_score" in multichain_result


def test_exact_synthetic_interface_match_reproduces_formal_score():
    scoring, _ = load_scoring_modules()
    residues = {"A:ALA:1:", "B:GLY:2:"}

    result = scoring.score_oligomer_relevance(
        residues,
        (0.0, 0.0, 0.0),
        make_interface(),
    )

    assert result["interface_fraction"] == 1.0
    assert result["interface_recall"] == 1.0
    assert result["cluster_interface_residue_count"] == 2
    assert result["cluster_chain_count"] == 2
    assert result["cluster_chain_entropy"] == 1.0
    assert result["distance_to_interface_A"] == 0.0
    assert result["interface_distance_score"] == 100.0
    assert result["context_score"] == 56.6
    assert result["O_rel_formal"] == pytest.approx(94.358)


def test_chain_entropy_is_normalized_and_empty_residues_are_explicit():
    scoring, _ = load_scoring_modules()
    interface = make_interface()

    unequal = scoring.score_oligomer_relevance(
        {"A:ALA:1:", "A:SER:2:", "B:GLY:3:"},
        None,
        interface,
    )
    single = scoring.score_oligomer_relevance({"A:ALA:1:"}, None, interface)
    empty = scoring.score_oligomer_relevance((), None, interface)

    expected_unequal = -(
        (2 / 3) * math.log(2 / 3) + (1 / 3) * math.log(1 / 3)
    ) / math.log(2)
    assert unequal["cluster_chain_entropy"] == pytest.approx(expected_unequal)
    assert single["cluster_chain_entropy"] == 0.0
    assert empty["cluster_chain_entropy"] == 0.0
    assert empty["cluster_chain_count"] == 0
    assert empty["cluster_interface_residue_count"] == 0
    assert math.isnan(empty["interface_fraction"])
    assert empty["interface_recall"] == 0.0
    assert math.isnan(empty["distance_to_interface_A"])
    assert empty["O_rel_formal"] == pytest.approx(31.6975)


def test_formal_score_renormalizes_when_interface_distance_is_unavailable():
    scoring, _ = load_scoring_modules()

    result = scoring.score_oligomer_relevance(
        {"A:ALA:1:", "B:GLY:2:"},
        None,
        make_interface(),
    )

    assert math.isnan(result["interface_distance_score"])
    assert result["O_rel_formal"] == pytest.approx(92.9475)
