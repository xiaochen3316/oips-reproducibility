from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[2]
STRUCTURES = ROOT / "data" / "structures"


def load_structure_module():
    try:
        from oips_repro import structure
    except ImportError as exc:
        pytest.fail(f"public structure module is missing: {exc}")
    return structure


def pdb_atom_line(
    *,
    record: str = "ATOM",
    serial: int = 1,
    atom_name: str = "CA",
    residue_name: str = "ALA",
    chain: str = "A",
    residue_number: int = 1,
    insertion_code: str = " ",
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0),
    element: str = "C",
) -> str:
    x, y, z = coordinates
    return (
        f"{record:<6}{serial:>5} {atom_name:^4} {residue_name:>3} {chain:1}"
        f"{residue_number:>4}{insertion_code:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{1.0:>6.2f}{20.0:>6.2f}          {element:>2}\n"
    )


def make_atom(
    *,
    record: str = "ATOM",
    atom_name: str = "CA",
    residue_name: str = "ALA",
    chain: str = "A",
    residue_number: str = "1",
    insertion_code: str = "",
    element: str = "C",
    coordinates: tuple[float, float, float] = (0.0, 0.0, 0.0),
):
    structure = load_structure_module()
    return structure.Atom(
        record=record,
        atom_name=atom_name,
        residue_name=residue_name,
        chain=chain,
        residue_number=residue_number,
        insertion_code=insertion_code,
        element=element,
        coordinates=coordinates,
    )


def test_pdb_parser_preserves_records_and_canonicalizes_residue_ids(tmp_path: Path):
    structure = load_structure_module()
    path = tmp_path / "synthetic.ent"
    path.write_text(
        pdb_atom_line(
            residue_name="gly",
            chain=" ",
            residue_number=7,
            insertion_code="A",
            coordinates=(1.0, 2.0, 3.0),
        )
        + pdb_atom_line(
            record="HETATM",
            serial=2,
            atom_name="O",
            residue_name="HOH",
            chain="B",
            residue_number=8,
            coordinates=(4.0, 5.0, 6.0),
            element="O",
        ),
        encoding="utf-8",
    )

    atoms = structure.read_structure_atoms(path)

    assert len(atoms) == 2
    assert atoms[0].record == "ATOM"
    assert atoms[0].atom_name == "CA"
    assert atoms[0].residue_name == "GLY"
    assert atoms[0].chain == ""
    assert atoms[0].residue_number == "7"
    assert atoms[0].insertion_code == "A"
    assert atoms[0].element == "C"
    assert atoms[0].coordinates == (1.0, 2.0, 3.0)
    assert atoms[0].residue_id == "_:GLY:7:A"
    assert atoms[1].record == "HETATM"
    assert atoms[1].residue_id == "B:HOH:8:"
    with pytest.raises(FrozenInstanceError):
        atoms[0].chain = "X"


def test_mmcif_parser_uses_atom_site_headers_instead_of_fixed_positions(tmp_path: Path):
    structure = load_structure_module()
    path = tmp_path / "reordered.mmcif"
    path.write_text(
        """data_reordered
loop_
_atom_site.Cartn_z
_atom_site.label_comp_id
_atom_site.group_PDB
_atom_site.auth_asym_id
_atom_site.label_atom_id
_atom_site.Cartn_x
_atom_site.pdbx_PDB_ins_code
_atom_site.label_seq_id
_atom_site.type_symbol
_atom_site.label_asym_id
_atom_site.Cartn_y
3.5 ser ATOM ? OG 1.5 A 42 O . 2.5
#
""",
        encoding="utf-8",
    )

    atoms = structure.read_structure_atoms(path)

    assert len(atoms) == 1
    assert atoms[0].record == "ATOM"
    assert atoms[0].atom_name == "OG"
    assert atoms[0].residue_name == "SER"
    assert atoms[0].chain == ""
    assert atoms[0].element == "O"
    assert atoms[0].coordinates == (1.5, 2.5, 3.5)
    assert atoms[0].residue_id == "_:SER:42:A"


@pytest.mark.parametrize("missing_element", [".", "?"])
def test_mmcif_missing_element_is_inferred_for_hydrogen_exclusion(
    tmp_path: Path,
    missing_element: str,
):
    structure = load_structure_module()
    path = tmp_path / f"missing-element-{ord(missing_element)}.cif"
    path.write_text(
        f"""data_missing_element
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM C CA ALA A 1 ? 0.0 0.0 0.0
ATOM C CA GLY B 2 ? 4.0 0.0 0.0
ATOM {missing_element} H1 SER C 3 ? 2.0 0.0 0.0
#
""",
        encoding="utf-8",
    )

    atoms = structure.read_structure_atoms(path)
    profile = structure.build_interface_profile(atoms)

    assert atoms[2].element == "H"
    assert profile.chain_count == 2
    assert profile.protein_residue_count == 2
    assert dict(profile.pair_contact_counts) == {"A-B": 1}


def test_committed_header_driven_cif_parses():
    structure = load_structure_module()

    atoms = structure.read_structure_atoms(STRUCTURES / "5tbm_paired.cif")

    assert atoms
    assert atoms[0].record == "ATOM"
    assert atoms[0].residue_id == "A:LEU:239:"


def test_all_21_committed_structures_parse_to_nonempty_atom_tuples():
    structure = load_structure_module()
    paths = sorted([*STRUCTURES.glob("*.pdb"), *STRUCTURES.glob("*.cif")])

    assert len(paths) == 21
    assert all(structure.read_structure_atoms(path) for path in paths)


def test_structure_reader_rejects_missing_and_unsupported_paths(tmp_path: Path):
    structure = load_structure_module()

    with pytest.raises(FileNotFoundError, match="structure file does not exist"):
        structure.read_structure_atoms(tmp_path / "missing.pdb")

    unsupported = tmp_path / "structure.xyz"
    unsupported.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported structure format"):
        structure.read_structure_atoms(unsupported)


def test_centroid_and_euclidean_distance_include_empty_centroid_behavior():
    structure = load_structure_module()

    assert structure.centroid([]) is None
    assert structure.centroid([(0.0, 0.0, 0.0), (2.0, 4.0, 6.0)]) == (
        1.0,
        2.0,
        3.0,
    )
    assert structure.euclidean_distance((0.0, 0.0, 0.0), (3.0, 4.0, 12.0)) == 13.0


def test_two_chain_contact_includes_the_exact_cutoff_boundary():
    structure = load_structure_module()
    atoms = (
        make_atom(chain="B", residue_name="GLY", residue_number="2", coordinates=(5.0, 0.0, 0.0)),
        make_atom(chain="A", residue_name="ALA", residue_number="1", coordinates=(0.0, 0.0, 0.0)),
    )

    profile = structure.build_interface_profile(atoms, contact_cutoff_A=5.0)

    assert profile.chain_count == 2
    assert profile.protein_residue_count == 2
    assert profile.interface_residues == ("A:ALA:1:", "B:GLY:2:")
    assert profile.interface_atom_coordinates == ((0.0, 0.0, 0.0), (5.0, 0.0, 0.0))
    assert dict(profile.pair_contact_counts) == {"A-B": 1}
    assert profile.interface_pair_count == 1
    assert profile.interface_atom_contact_count == 1
    with pytest.raises(TypeError):
        profile.pair_contact_counts["A-B"] = 2


def test_interface_construction_excludes_hetero_hydrogen_and_nonstandard_atoms():
    structure = load_structure_module()
    atoms = (
        make_atom(chain="A", residue_name="ALA", residue_number="1", coordinates=(0.0, 0.0, 0.0)),
        make_atom(chain="B", residue_name="GLY", residue_number="2", coordinates=(4.0, 0.0, 0.0)),
        make_atom(record="HETATM", chain="C", residue_name="ALA", residue_number="3", coordinates=(1.0, 0.0, 0.0)),
        make_atom(atom_name="H", chain="D", residue_name="SER", residue_number="4", element="H", coordinates=(2.0, 0.0, 0.0)),
        make_atom(chain="E", residue_name="MSE", residue_number="5", coordinates=(3.0, 0.0, 0.0)),
    )

    profile = structure.build_interface_profile(atoms)

    assert profile.chain_count == 2
    assert profile.protein_residue_count == 2
    assert profile.interface_residues == ("A:ALA:1:", "B:GLY:2:")
    assert profile.interface_atom_contact_count == 1


def test_interface_output_is_input_order_invariant():
    structure = load_structure_module()
    atoms = (
        make_atom(chain="C", residue_name="SER", residue_number="3", coordinates=(9.0, 0.0, 0.0)),
        make_atom(chain="A", residue_name="ALA", residue_number="1", coordinates=(0.0, 0.0, 0.0)),
        make_atom(chain="B", residue_name="GLY", residue_number="2", coordinates=(4.0, 0.0, 0.0)),
    )

    forward = structure.build_interface_profile(atoms)
    reverse = structure.build_interface_profile(tuple(reversed(atoms)))

    assert forward == reverse
    assert tuple(forward.pair_contact_counts) == ("A-B", "B-C")
    assert forward.interface_atom_coordinates == (
        (0.0, 0.0, 0.0),
        (4.0, 0.0, 0.0),
        (4.0, 0.0, 0.0),
        (9.0, 0.0, 0.0),
    )
