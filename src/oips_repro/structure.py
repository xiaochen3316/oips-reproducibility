"""Pure structure parsing and protein-interface construction primitives."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
import math
from pathlib import Path
import shlex
from types import MappingProxyType

from scipy.spatial import cKDTree


Coordinate = tuple[float, float, float]

STANDARD_AMINO_ACIDS = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)


def _optional_text(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in {".", "?"} else text


def _residue_id(chain: str, residue_name: str, residue_number: str, insertion_code: str) -> str:
    number = _optional_text(residue_number)
    if not number:
        return ""
    canonical_chain = _optional_text(chain) or "_"
    name = _optional_text(residue_name).upper()
    insertion = _optional_text(insertion_code)
    return f"{canonical_chain}:{name}:{number}:{insertion}"


@dataclass(frozen=True)
class Atom:
    """One immutable coordinate record parsed from PDB or mmCIF."""

    record: str
    atom_name: str
    residue_name: str
    chain: str
    residue_number: str
    insertion_code: str
    element: str
    coordinates: Coordinate

    def __post_init__(self) -> None:
        object.__setattr__(self, "record", self.record.strip().upper())
        object.__setattr__(self, "atom_name", self.atom_name.strip())
        object.__setattr__(self, "residue_name", self.residue_name.strip().upper())
        object.__setattr__(self, "chain", _optional_text(self.chain))
        object.__setattr__(self, "residue_number", _optional_text(self.residue_number))
        object.__setattr__(self, "insertion_code", _optional_text(self.insertion_code))
        object.__setattr__(self, "element", self.element.strip().upper())
        object.__setattr__(
            self,
            "coordinates",
            tuple(float(value) for value in self.coordinates),
        )

    @property
    def residue_id(self) -> str:
        return _residue_id(
            self.chain,
            self.residue_name,
            self.residue_number,
            self.insertion_code,
        )

    # Read-only aliases preserve the naming used by the private reference scripts.
    @property
    def atom(self) -> str:
        return self.atom_name

    @property
    def resname(self) -> str:
        return self.residue_name

    @property
    def resnum(self) -> str:
        return self.residue_number

    @property
    def icode(self) -> str:
        return self.insertion_code

    @property
    def elem(self) -> str:
        return self.element

    @property
    def xyz(self) -> Coordinate:
        return self.coordinates

    @property
    def residue(self) -> str:
        return self.residue_id


@dataclass(frozen=True)
class InterfaceProfile:
    """Immutable summary of inter-chain heavy-atom protein contacts."""

    chain_count: int
    protein_residue_count: int
    interface_residues: tuple[str, ...] = ()
    interface_atom_coordinates: tuple[Coordinate, ...] = ()
    pair_contact_counts: Mapping[str, int] = field(default_factory=dict)
    interface_pair_count: int = field(init=False)
    interface_atom_contact_count: int = field(init=False)

    def __post_init__(self) -> None:
        residues = tuple(sorted(set(self.interface_residues)))
        coordinates = tuple(
            sorted(tuple(float(value) for value in coordinate) for coordinate in self.interface_atom_coordinates)
        )
        pair_counts = {
            str(pair): int(count)
            for pair, count in sorted(self.pair_contact_counts.items())
        }
        object.__setattr__(self, "chain_count", int(self.chain_count))
        object.__setattr__(self, "protein_residue_count", int(self.protein_residue_count))
        object.__setattr__(self, "interface_residues", residues)
        object.__setattr__(self, "interface_atom_coordinates", coordinates)
        object.__setattr__(self, "pair_contact_counts", MappingProxyType(pair_counts))
        object.__setattr__(self, "interface_pair_count", len(pair_counts))
        object.__setattr__(self, "interface_atom_contact_count", sum(pair_counts.values()))

    @property
    def interface_atom_coords(self) -> tuple[Coordinate, ...]:
        """Reference-script-compatible alias for interface atom coordinates."""

        return self.interface_atom_coordinates


def _parse_pdb(path: Path) -> tuple[Atom, ...]:
    atoms: list[Atom] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[:6].strip().upper()
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                coordinates = (
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            except ValueError:
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper()
            if not element:
                element = next((character.upper() for character in atom_name if character.isalpha()), "")
            atoms.append(
                Atom(
                    record=record,
                    atom_name=atom_name,
                    residue_name=line[17:20],
                    chain=line[21:22],
                    residue_number=line[22:26],
                    insertion_code=line[26:27],
                    element=element,
                    coordinates=coordinates,
                )
            )
    return tuple(atoms)


def _tokenize_cif_line(line: str) -> list[str]:
    try:
        return shlex.split(line, comments=True, posix=True)
    except ValueError:
        return line.split()


def _cif_value(row: Mapping[str, str], *names: str) -> str:
    for name in names:
        value = _optional_text(row.get(name.lower(), ""))
        if value:
            return value
    return ""


def _atoms_from_cif_rows(headers: Sequence[str], tokens: Sequence[str]) -> list[Atom]:
    width = len(headers)
    if not width:
        return []
    lowered_headers = tuple(header.lower() for header in headers)
    atoms: list[Atom] = []
    for start in range(0, len(tokens) - width + 1, width):
        row = dict(zip(lowered_headers, tokens[start : start + width]))
        record = _cif_value(row, "_atom_site.group_PDB").upper()
        if record not in {"ATOM", "HETATM"}:
            continue
        try:
            coordinates = (
                float(_cif_value(row, "_atom_site.Cartn_x")),
                float(_cif_value(row, "_atom_site.Cartn_y")),
                float(_cif_value(row, "_atom_site.Cartn_z")),
            )
        except ValueError:
            continue
        atom_name = _cif_value(row, "_atom_site.auth_atom_id", "_atom_site.label_atom_id")
        element = _cif_value(row, "_atom_site.type_symbol")
        if not element:
            element = next((character.upper() for character in atom_name if character.isalpha()), "")
        atoms.append(
            Atom(
                record=record,
                atom_name=atom_name,
                residue_name=_cif_value(row, "_atom_site.auth_comp_id", "_atom_site.label_comp_id"),
                chain=_cif_value(row, "_atom_site.auth_asym_id", "_atom_site.label_asym_id"),
                residue_number=_cif_value(row, "_atom_site.auth_seq_id", "_atom_site.label_seq_id"),
                insertion_code=_cif_value(row, "_atom_site.pdbx_PDB_ins_code"),
                element=element,
                coordinates=coordinates,
            )
        )
    return atoms


def _parse_mmcif(path: Path) -> tuple[Atom, ...]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    atoms: list[Atom] = []
    index = 0
    while index < len(lines):
        if lines[index].strip().lower() != "loop_":
            index += 1
            continue

        index += 1
        headers: list[str] = []
        while index < len(lines) and lines[index].lstrip().startswith("_"):
            header_tokens = _tokenize_cif_line(lines[index])
            if header_tokens:
                headers.append(header_tokens[0])
            index += 1

        tokens: list[str] = []
        while index < len(lines):
            stripped = lines[index].strip()
            if not stripped:
                index += 1
                continue
            if stripped.startswith("#"):
                index += 1
                break
            lowered = stripped.lower()
            if lowered == "loop_" or stripped.startswith("_") or lowered.startswith(("data_", "save_")):
                break
            tokens.extend(_tokenize_cif_line(stripped))
            index += 1

        if headers and all(header.lower().startswith("_atom_site.") for header in headers):
            atoms.extend(_atoms_from_cif_rows(headers, tokens))

    return tuple(atoms)


def read_structure_atoms(path: str | Path) -> tuple[Atom, ...]:
    """Read all ATOM/HETATM records from a supported coordinate file.

    PDB ``.pdb``/``.ent`` and mmCIF ``.cif``/``.mmcif`` files are supported.
    """

    structure_path = Path(path)
    if not structure_path.is_file():
        raise FileNotFoundError(f"structure file does not exist: {structure_path}")
    suffix = structure_path.suffix.lower()
    if suffix in {".pdb", ".ent"}:
        return _parse_pdb(structure_path)
    if suffix in {".cif", ".mmcif"}:
        return _parse_mmcif(structure_path)
    raise ValueError(f"unsupported structure format {suffix!r}: {structure_path}")


def centroid(coordinates: Iterable[Sequence[float]]) -> Coordinate | None:
    """Return the arithmetic coordinate centroid, or ``None`` for no points."""

    points = tuple(tuple(float(value) for value in point) for point in coordinates)
    if not points:
        return None
    count = len(points)
    return tuple(sum(point[axis] for point in points) / count for axis in range(3))


def euclidean_distance(a: Sequence[float], b: Sequence[float]) -> float:
    """Return the three-dimensional Euclidean distance between two points."""

    return float(math.sqrt(sum((float(a[index]) - float(b[index])) ** 2 for index in range(3))))


def _atom_sort_key(atom: Atom) -> tuple[object, ...]:
    return (
        atom.chain,
        atom.residue_number,
        atom.insertion_code,
        atom.residue_name,
        atom.atom_name,
        atom.element,
        atom.coordinates,
    )


def build_interface_profile(
    atoms: Iterable[Atom],
    contact_cutoff_A: float = 5.0,
) -> InterfaceProfile:
    """Build an inter-chain interface from standard-residue heavy ``ATOM`` records."""

    cutoff = float(contact_cutoff_A)
    if cutoff < 0.0:
        raise ValueError("contact_cutoff_A must be nonnegative")

    protein_atoms = tuple(
        sorted(
            (
                atom
                for atom in atoms
                if atom.record == "ATOM"
                and atom.residue_name in STANDARD_AMINO_ACIDS
                and atom.element != "H"
            ),
            key=_atom_sort_key,
        )
    )
    by_chain: defaultdict[str, list[Atom]] = defaultdict(list)
    for atom in protein_atoms:
        by_chain[atom.chain].append(atom)
    chains = sorted(by_chain)

    interface_residues: set[str] = set()
    interface_coordinates: list[Coordinate] = []
    pair_contact_counts: dict[str, int] = {}
    for chain_index, chain_a in enumerate(chains):
        atoms_a = by_chain[chain_a]
        coordinates_a = [atom.coordinates for atom in atoms_a]
        for chain_b in chains[chain_index + 1 :]:
            atoms_b = by_chain[chain_b]
            coordinates_b = [atom.coordinates for atom in atoms_b]
            neighbors = cKDTree(coordinates_b).query_ball_point(coordinates_a, r=cutoff)
            touched_a = {index for index, hits in enumerate(neighbors) if hits}
            touched_b = {index for hits in neighbors for index in hits}
            contact_count = sum(len(hits) for hits in neighbors)
            if not contact_count:
                continue

            pair_name = f"{chain_a or '_'}-{chain_b or '_'}"
            pair_contact_counts[pair_name] = contact_count
            for index in sorted(touched_a):
                atom = atoms_a[index]
                interface_residues.add(atom.residue_id)
                interface_coordinates.append(atom.coordinates)
            for index in sorted(touched_b):
                atom = atoms_b[index]
                interface_residues.add(atom.residue_id)
                interface_coordinates.append(atom.coordinates)

    return InterfaceProfile(
        chain_count=len(chains),
        protein_residue_count=len({atom.residue_id for atom in protein_atoms}),
        interface_residues=tuple(interface_residues),
        interface_atom_coordinates=tuple(interface_coordinates),
        pair_contact_counts=pair_contact_counts,
    )
