"""Pure post-static evidence mapping for the frozen OIPS analysis."""
from __future__ import annotations
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from pathlib import Path
import re
from types import MappingProxyType
import numpy as np
import pandas as pd
from .clustering import parse_residues, simple_set
from .io import RANKING_COLUMNS
from .structure import STANDARD_AMINO_ACIDS, Atom, centroid, euclidean_distance, read_structure_atoms

def _columns(value: str) -> tuple[str, ...]:
    return tuple(value.split())
TOPK_QC_COLUMNS = _columns('pdb_id cluster_v2_id Within_PDB_rank OIPS-P_static\n mappability center_available residue_set_available cluster_diameter_A spatial_continuity\n tool_support_count center_dispersion_A core_envelope_ratio contributing_chain_count\n dominant_chain_fraction pocket_interface_overlap distance_to_interface_A possible_over_merging\n clear_exclusion_flag possible_split_subpocket possible_surface_noise_cluster\n nearest_missing_residue_distance_sequence_positions missing_residue_proximity_flag\n missing_residue_assessment_source possible_crystal_contact_risk boundary_sensitive QC_status QC_rule_version')
REFERENCE_COLUMNS = _columns('pdb_id cluster_v2_id Within_PDB_rank DCC_A\n reference_contact_residue_count reference_contact_overlap_count contact_precision contact_recall\n residue_IoU functional_residue_overlap functional_residue_status interface_overlap\n distance_to_interface_A chain_contribution_count R_auto_rule_pass R_auto_rule\n reference_ligand_annotation requested_reference_ligand_key selected_reference_ligand_key\n reference_ligand_atom_count reference_selection_status reference_selection_unresolved')
MD_COLUMNS = _columns('pdb_id MD_run Simulation_context static_top_cluster_v2_id static_top_rank\n persistent_MD_contact_residue_count best_MD_mapped_cluster_v2_id best_MD_cluster_Jaccard\n best_MD_cluster_precision best_MD_cluster_recall best_MD_cluster_center_distance_A\n static_top_MD_overlap_count static_top_MD_cluster_coverage static_top_MD_contact_coverage\n static_top_MD_Jaccard static_top_MD_center_distance_A D_dyn_run_score Concordance_call\n MD_mapping_rule_version source_persistent_contacts')
MD_INPUT_COLUMNS = _columns('pdb_id MD_run Simulation_context Persistent_MD_contact_residues MD_contact_center D_dyn_run_score')
LABEL_COLUMNS = _columns('pdb_id cluster_v2_id Within_PDB_rank automated_evidence_label\n label_reason independent_evidence_details R_auto_DCC_A R_auto_contact_recall R_auto_residue_IoU\n QC_status unresolved_flag label_scope interpretation_guardrail')
REDOCKING_COLUMNS = _columns('pdb_id reference_cluster_v2_id reference_cluster_static_rank\n Raw_ligand_RMSD_A GlideScore_kcal_per_mol RMSD_threshold_call Reference_pose_recovered\n Failure_or_warning_reason redocking_role')
UNRESOLVED_COLUMNS = _columns('issue_id level pdb_id cluster_v2_id issue_type status\n required_user_or_future_input current_handling')
CONVENIENCE_COLUMNS = _columns('automated_evidence_label label_reason unresolved_flag QC_status\n DCC_A contact_recall contact_precision residue_IoU R_auto_rule_pass')
VALID_MAPPABILITY = frozenset(_columns('center_and_residue_mappable center_only_mappable residue_only_mappable'))
IGNORE_HET = frozenset(_columns('HOH WAT DOD NA CL K CA MG MN ZN FE CU CO NI CD SO4 PO4 GOL EDO PEG ACT ACE FMT DMS DMSO'))
R_AUTO_RULE_TEXT = 'DCC<=6A+recall>=0.10 OR DCC<=10A+(recall>=0.20 or IoU>=0.15) OR DCC<=12A+recall>=0.50+IoU>=0.15'

@dataclass(frozen=True)
class QCRules:
    overmerge_diameter_gt_A: float = 10.5
    overmerge_core_ratio_lt: float = 0.25
    overmerge_secondary_units_min: int = 4
    overmerge_tool_support_max: int = 2
    overmerge_median_iou_lt: float = 0.15
    exclusion_diameter_gt_A: float = 11.5
    exclusion_core_ratio_lt: float = 0.1
    exclusion_tool_support_min: int = 2
    split_tool_support_min: int = 2
    split_dispersion_gt_A: float = 4.0
    split_median_iou_lt: float = 0.15
    split_core_ratio_lt: float = 0.2
    surface_noise_tool_support: int = 1
    surface_noise_geometry_lt: float = 45.0
    surface_noise_ligandability_lt: float = 45.0
    surface_noise_interface_fraction_lt: float = 0.1

@dataclass(frozen=True)
class ReferenceRules:
    contact_cutoff_A: float = 4.5
    near_dcc_max_A: float = 6.0
    near_recall_min: float = 0.1
    middle_dcc_max_A: float = 10.0
    middle_recall_min: float = 0.2
    middle_iou_min: float = 0.15
    far_dcc_max_A: float = 12.0
    far_recall_min: float = 0.5
    far_iou_min: float = 0.15

@dataclass(frozen=True)
class MDRules:
    concordant_iou_min: float = 0.15
    concordant_precision_min: float = 0.2
    concordant_center_max_A: float = 6.0
    partial_iou_min: float = 0.05
    alternative_iou_min: float = 0.05
    alternative_center_max_A: float = 6.0
    boundary_shift_center_max_A: float = 8.0

@dataclass(frozen=True)
class PosthocInputs:
    rankings: pd.DataFrame
    structures: Mapping[str, Path]
    references: pd.DataFrame
    md_runs: pd.DataFrame
    redocking: pd.DataFrame
    systems: pd.DataFrame

    def __post_init__(self) -> None:
        for name in ('rankings', 'references', 'md_runs', 'redocking', 'systems'):
            object.__setattr__(self, name, getattr(self, name).copy(deep=True))
        paths = {str(key).upper(): Path(value) for key, value in self.structures.items()}
        object.__setattr__(self, 'structures', MappingProxyType(paths))

@dataclass(frozen=True)
class PosthocResult:
    topk_qc: pd.DataFrame
    reference_mapping: pd.DataFrame
    reference_evaluable: pd.Series
    md_mapping: pd.DataFrame
    evidence_labels: pd.DataFrame
    redocking_mapping: pd.DataFrame
    unresolved_cases: pd.DataFrame
    convenience_master: pd.DataFrame
    representative_case_ids: tuple[str, ...]

def _require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f'{label} is missing required columns: {missing}')

def _number(value: object, name: str, *, missing_ok: bool=False) -> float:
    if value is None or (isinstance(value, str) and (not value.strip())) or pd.isna(value):
        if missing_ok:
            return math.nan
        raise ValueError(f'{name} may not be missing')
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be numeric') from exc
    if not math.isfinite(result):
        if missing_ok and math.isnan(result):
            return math.nan
        raise ValueError(f'{name} must be finite')
    return result

def _integer(value: object, name: str) -> int:
    number = _number(value, name)
    if not number.is_integer():
        raise ValueError(f'{name} must be an integer')
    return int(number)

def _boolean(value: object, name: str, *, default: bool | None=None) -> bool:
    if value is None or (isinstance(value, str) and (not value.strip())) or pd.isna(value):
        if default is not None:
            return default
        raise ValueError(f'{name} may not be missing')
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str) and value.strip().casefold() in {'true', 'false'}:
        return value.strip().casefold() == 'true'
    raise ValueError(f'{name} must be true or false')

def _center_from_row(row: Mapping[str, object] | pd.Series) -> tuple[float, float, float] | None:
    values = [_number(row.get(column), column, missing_ok=True) for column in ('medoid_center_x', 'medoid_center_y', 'medoid_center_z')]
    missing = [math.isnan(value) for value in values]
    if all(missing):
        return None
    if any(missing):
        raise ValueError('medoid center coordinates must be all present or all missing')
    return tuple(values)

def _parse_center(value: object) -> tuple[float, float, float] | None:
    if value is None or (isinstance(value, str) and (not value.strip())) or pd.isna(value):
        return None
    numbers = re.findall('[-+]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][-+]?\\d+)?', str(value))
    if len(numbers) < 3:
        raise ValueError('MD_contact_center must contain three numeric coordinates')
    center = tuple((float(number) for number in numbers[:3]))
    if not all((math.isfinite(number) for number in center)):
        raise ValueError('MD_contact_center must contain finite coordinates')
    return center

def _overlap(cluster: set[str], evidence: set[str]) -> tuple[int, float, float, float]:
    left, right = (simple_set(cluster), simple_set(evidence))
    if not left or not right:
        return (0, math.nan, math.nan, math.nan)
    intersection = left & right
    return (len(intersection), len(intersection) / len(left), len(intersection) / len(right), len(intersection) / len(left | right))

def _unique_keys(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    if frame.duplicated(list(columns)).any():
        raise ValueError(f'{label} contains duplicate keys: {list(columns)}')

def _missing_residues(path: Path) -> tuple[dict[str, set[int]], str]:
    missing: defaultdict[str, set[int]] = defaultdict(set)
    if path.suffix.lower() in {'.pdb', '.ent'}:
        for line in path.read_text(encoding='utf-8', errors='ignore').splitlines():
            if line.startswith('REMARK 465'):
                match = re.search('REMARK 465\\s+[A-Z0-9]{3}\\s+([A-Za-z0-9])\\s+(-?\\d+)', line)
                if match:
                    missing[match.group(1)].add(int(match.group(2)))
    if missing:
        return (dict(missing), 'PDB_REMARK_465')
    observed: defaultdict[str, set[int]] = defaultdict(set)
    for atom in read_structure_atoms(path):
        match = re.match('-?\\d+', atom.residue_number)
        if atom.record == 'ATOM' and atom.residue_name in STANDARD_AMINO_ACIDS and match:
            observed[atom.chain].add(int(match.group()))
    for chain, numbers in observed.items():
        for left, right in zip(sorted(numbers), sorted(numbers)[1:]):
            if 2 <= right - left <= 25:
                missing[chain].update(range(left + 1, right))
    return (dict(missing), 'observed_numbering_gap_proxy' if missing else 'no_gap_detected')

def _nearest_missing(residues: set[str], missing: Mapping[str, set[int]]) -> float:
    distances: list[int] = []
    for residue in residues:
        parts = residue.split(':')
        match = re.match('-?\\d+', parts[2]) if len(parts) >= 3 else None
        if match and missing.get(parts[0]):
            distances.append(min((abs(int(match.group()) - number) for number in missing[parts[0]])))
    return float(min(distances)) if distances else math.nan

def run_topk_qc(rankings: pd.DataFrame, structures: Mapping[str, Path], top_k: int=3, rules: QCRules=QCRules()) -> pd.DataFrame:
    """Return deterministic automated structural QC for the static top-k."""
    required = ('pdb_id', 'cluster_v2_id', 'Within_PDB_rank', 'OIPS-P_static_recomputed', 'mappability', 'cluster_diameter_A', 'spatial_continuity', 'tool_support_count', 'center_dispersion_A', 'core_envelope_ratio', 'contributing_chain_count', 'dominant_chain_fraction', 'interface_fraction', 'distance_to_interface_A', 'same_tool_secondary_unit_count', 'G_geo', 'P_lig', 'center_available_representatives', 'residue_available_representatives', 'pairwise_residue_iou_median', 'boundary_sensitive', 'envelope_residues')
    _require_columns(rankings, required, 'rankings')
    _unique_keys(rankings, ('pdb_id', 'cluster_v2_id'), 'rankings')
    if int(top_k) < 1:
        raise ValueError('top_k must be positive')
    paths = {str(key).upper(): Path(value) for key, value in structures.items()}
    selected = rankings.loc[pd.to_numeric(rankings['Within_PDB_rank'], errors='raise').le(top_k)].copy()
    selected['pdb_id'] = selected['pdb_id'].astype(str).str.upper()
    cache: dict[str, tuple[dict[str, set[int]], str, str]] = {}
    rows: list[dict[str, object]] = []
    for row in selected.sort_values(['pdb_id', 'Within_PDB_rank', 'cluster_v2_id'], kind='mergesort').to_dict('records'):
        pdb_id = str(row['pdb_id'])
        if pdb_id not in paths or not paths[pdb_id].is_file():
            raise FileNotFoundError(f'prepared structure is missing for {pdb_id}')
        if pdb_id not in cache:
            path = paths[pdb_id]
            missing, source = _missing_residues(path)
            cryst = 'not_evaluable_no_crystal_symmetry_model'
            if path.suffix.lower() in {'.pdb', '.ent'} and any((line.startswith('CRYST1') for line in path.read_text(encoding='utf-8', errors='ignore').splitlines())):
                cryst = 'not_evaluable_no_generated_crystal_mates'
            cache[pdb_id] = (missing, source, cryst)
        diameter = _number(row['cluster_diameter_A'], 'cluster_diameter_A', missing_ok=True)
        ratio = _number(row['core_envelope_ratio'], 'core_envelope_ratio', missing_ok=True)
        dispersion = _number(row['center_dispersion_A'], 'center_dispersion_A', missing_ok=True)
        median_iou = _number(row['pairwise_residue_iou_median'], 'pairwise_residue_iou_median', missing_ok=True)
        support = _integer(row['tool_support_count'], 'tool_support_count')
        secondary = _integer(row['same_tool_secondary_unit_count'], 'same_tool_secondary_unit_count')
        spatial = _boolean(row['spatial_continuity'], 'spatial_continuity')
        overmerge = diameter > rules.overmerge_diameter_gt_A and ratio < rules.overmerge_core_ratio_lt or not spatial or (secondary >= rules.overmerge_secondary_units_min and support <= rules.overmerge_tool_support_max and (median_iou < rules.overmerge_median_iou_lt))
        exclusion = not spatial or (diameter > rules.exclusion_diameter_gt_A and ratio < rules.exclusion_core_ratio_lt and (support >= rules.exclusion_tool_support_min))
        split = not overmerge and support >= rules.split_tool_support_min and (dispersion > rules.split_dispersion_gt_A or median_iou < rules.split_median_iou_lt or ratio < rules.split_core_ratio_lt)
        geometry = _number(row['G_geo'], 'G_geo', missing_ok=True)
        ligandability = _number(row['P_lig'], 'P_lig', missing_ok=True)
        interface = _number(row['interface_fraction'], 'interface_fraction', missing_ok=True)
        surface = support == rules.surface_noise_tool_support and geometry < rules.surface_noise_geometry_lt and (math.isnan(ligandability) or ligandability < rules.surface_noise_ligandability_lt) and (math.isnan(interface) or interface < rules.surface_noise_interface_fraction_lt)
        mappability = str(row['mappability'])
        insufficient = mappability == 'center_only_mappable' or (support == 1 and math.isnan(ligandability)) or surface
        boundary = _boolean(row['boundary_sensitive'], 'boundary_sensitive')
        if mappability not in VALID_MAPPABILITY:
            status = 'QC_unmappable'
        elif overmerge:
            status = 'QC_possible_overmerge'
        elif split:
            status = 'QC_possible_split'
        elif boundary:
            status = 'QC_boundary_sensitive'
        elif insufficient:
            status = 'QC_insufficient_evidence'
        else:
            status = 'QC_pass'
        residues = parse_residues(row['envelope_residues'])
        distance = _nearest_missing(residues, cache[pdb_id][0])
        rows.append({'pdb_id': pdb_id, 'cluster_v2_id': row['cluster_v2_id'], 'Within_PDB_rank': _number(row['Within_PDB_rank'], 'Within_PDB_rank'), 'OIPS-P_static': _number(row['OIPS-P_static_recomputed'], 'OIPS-P-static'), 'mappability': mappability, 'center_available': _integer(row['center_available_representatives'], 'center_available_representatives') > 0, 'residue_set_available': _integer(row['residue_available_representatives'], 'residue_available_representatives') > 0, 'cluster_diameter_A': diameter, 'spatial_continuity': spatial, 'tool_support_count': support, 'center_dispersion_A': dispersion, 'core_envelope_ratio': ratio, 'contributing_chain_count': _integer(row['contributing_chain_count'], 'contributing_chain_count'), 'dominant_chain_fraction': _number(row['dominant_chain_fraction'], 'dominant_chain_fraction', missing_ok=True), 'pocket_interface_overlap': interface, 'distance_to_interface_A': _number(row['distance_to_interface_A'], 'distance_to_interface_A', missing_ok=True), 'possible_over_merging': overmerge, 'clear_exclusion_flag': exclusion, 'possible_split_subpocket': split, 'possible_surface_noise_cluster': surface, 'nearest_missing_residue_distance_sequence_positions': distance, 'missing_residue_proximity_flag': not math.isnan(distance) and distance <= 2, 'missing_residue_assessment_source': cache[pdb_id][1], 'possible_crystal_contact_risk': cache[pdb_id][2], 'boundary_sensitive': boundary, 'QC_status': status, 'QC_rule_version': 'automated_QC_v2_20260711'})
    return pd.DataFrame(rows, columns=TOPK_QC_COLUMNS).sort_values(['pdb_id', 'Within_PDB_rank', 'cluster_v2_id'], kind='mergesort').reset_index(drop=True)

def reference_rule_pass(dcc: float, recall: float, iou: float, rules: ReferenceRules=ReferenceRules()) -> bool:
    """Evaluate the inclusive frozen reference-association rule."""
    if math.isnan(float(dcc)) or math.isnan(float(recall)):
        return False
    return bool(dcc <= rules.near_dcc_max_A and recall >= rules.near_recall_min or (dcc <= rules.middle_dcc_max_A and (recall >= rules.middle_recall_min or (not math.isnan(float(iou)) and iou >= rules.middle_iou_min))) or (dcc <= rules.far_dcc_max_A and recall >= rules.far_recall_min and (not math.isnan(float(iou))) and (iou >= rules.far_iou_min)))

def _ligand_key(value: object) -> tuple[str, str, str, str] | None:
    parts = str(value).split(':')
    return (parts[0], parts[1], parts[2], parts[3] if len(parts) > 3 else '') if len(parts) >= 3 else None

def _select_ligand(atoms: Sequence[Atom], requested: tuple[str, str, str, str] | None):
    groups: defaultdict[tuple[str, str, str, str], list[Atom]] = defaultdict(list)
    for atom in atoms:
        if atom.record == 'HETATM' and atom.residue_name not in IGNORE_HET:
            groups[atom.residue_name, atom.chain, atom.residue_number, atom.insertion_code].append(atom)
    if requested is None:
        return (None, [], 'invalid_project_reference_key')
    if requested in groups:
        return (requested, groups[requested], 'exact_project_reference_key')
    same = [(key, value) for key, value in groups.items() if key[0] == requested[0]]
    if len(same) == 1:
        return (same[0][0], same[0][1], 'unique_resname_fallback')
    if same:
        key, values = sorted(same, key=lambda item: (-len(item[1]), item[0]))[0]
        return (key, values, 'ambiguous_resname_largest_group_fallback')
    return (None, [], 'reference_ligand_not_found')

def map_reference_evidence(rankings: pd.DataFrame, structures: Mapping[str, Path], references: pd.DataFrame, rules: ReferenceRules=ReferenceRules()) -> pd.DataFrame:
    """Map curated reference ligands after the static ranking is frozen."""
    ranking_required = ('pdb_id', 'cluster_v2_id', 'Within_PDB_rank', 'medoid_center_x', 'medoid_center_y', 'medoid_center_z', 'envelope_residues', 'interface_fraction', 'distance_to_interface_A', 'contributing_chain_count')
    reference_required = ('pdb_id', 'selected_ligand', 'ligand_status', 'reference_ligand_annotation', 'override_applied', 'decision_id')
    _require_columns(rankings, ranking_required, 'rankings')
    _require_columns(references, reference_required, 'references')
    _unique_keys(rankings, ('pdb_id', 'cluster_v2_id'), 'rankings')
    _unique_keys(references, ('pdb_id',), 'references')
    refs = references.copy()
    refs['pdb_id'] = refs['pdb_id'].astype(str).str.upper()
    refs['override_applied'] = [_boolean(value, 'override_applied') for value in refs['override_applied']]
    ref_lookup = refs.set_index('pdb_id')
    paths = {str(key).upper(): Path(value) for key, value in structures.items()}
    work = rankings.copy()
    work['pdb_id'] = work['pdb_id'].astype(str).str.upper()
    rows: list[dict[str, object]] = []
    for pdb_id, target in work.groupby('pdb_id', sort=True):
        if pdb_id not in ref_lookup.index:
            raise ValueError(f'reference annotation is missing for {pdb_id}')
        if pdb_id not in paths or not paths[pdb_id].is_file():
            raise FileNotFoundError(f'prepared structure is missing for {pdb_id}')
        annotation = ref_lookup.loc[pdb_id]
        atoms = read_structure_atoms(paths[pdb_id])
        requested_text = str(annotation['reference_ligand_annotation'] if annotation['override_applied'] else annotation['selected_ligand'])
        selected, ligand_atoms, selection_status = _select_ligand(atoms, _ligand_key(requested_text))
        ligand_center = centroid(atom.coordinates for atom in ligand_atoms)
        if ligand_center is not None and not all(math.isfinite(value) for value in ligand_center):
            raise ValueError(f'selected reference ligand centroid must be finite for {pdb_id}')
        contacts = {atom.residue_id for atom in atoms if atom.record == 'ATOM' and atom.residue_name in STANDARD_AMINO_ACIDS and any((euclidean_distance(atom.coordinates, ligand.coordinates) <= rules.contact_cutoff_A for ligand in ligand_atoms))}
        unresolved = selection_status in {'ambiguous_resname_largest_group_fallback', 'reference_ligand_not_found', 'invalid_project_reference_key'}
        selected_text = ':'.join(selected) if selected else ''
        ordered = target.sort_values(['Within_PDB_rank', 'cluster_v2_id'], kind='mergesort')
        for row in ordered.to_dict('records'):
            center = _center_from_row(row)
            dcc = euclidean_distance(center, ligand_center) if center is not None and ligand_center is not None else math.nan
            overlap, precision, recall, iou = _overlap(parse_residues(row['envelope_residues']), contacts)
            rows.append({
                'pdb_id': pdb_id, 'cluster_v2_id': row['cluster_v2_id'],
                'Within_PDB_rank': _number(row['Within_PDB_rank'], 'Within_PDB_rank'),
                'DCC_A': dcc, 'reference_contact_residue_count': len(simple_set(contacts)),
                'reference_contact_overlap_count': overlap, 'contact_precision': precision,
                'contact_recall': recall, 'residue_IoU': iou,
                'functional_residue_overlap': math.nan,
                'functional_residue_status': 'not_available_not_imputed',
                'interface_overlap': _number(row['interface_fraction'], 'interface_fraction', missing_ok=True),
                'distance_to_interface_A': _number(row['distance_to_interface_A'], 'distance_to_interface_A', missing_ok=True),
                'chain_contribution_count': _integer(row['contributing_chain_count'], 'contributing_chain_count'),
                'R_auto_rule_pass': reference_rule_pass(dcc, recall, iou, rules),
                'R_auto_rule': R_AUTO_RULE_TEXT,
                'reference_ligand_annotation': annotation['reference_ligand_annotation'],
                'requested_reference_ligand_key': requested_text,
                'selected_reference_ligand_key': selected_text,
                'reference_ligand_atom_count': len(ligand_atoms),
                'reference_selection_status': selection_status,
                'reference_selection_unresolved': unresolved,
            })
    return pd.DataFrame(rows, columns=REFERENCE_COLUMNS).sort_values(
        ['pdb_id', 'Within_PDB_rank', 'cluster_v2_id'], kind='mergesort'
    ).reset_index(drop=True)

def derive_reference_evaluable(reference_mapping: pd.DataFrame) -> pd.Series:
    """Return target availability independently of whether any candidate is R_auto."""
    required = ('pdb_id', 'selected_reference_ligand_key', 'reference_ligand_atom_count')
    _require_columns(reference_mapping, required, 'reference mapping')
    values: dict[str, bool] = {}
    grouped = reference_mapping.groupby(reference_mapping['pdb_id'].astype(str).str.upper(), sort=True)
    for pdb_id, rows in grouped:
        atom_counts = pd.to_numeric(rows['reference_ligand_atom_count'], errors='raise')
        usable = rows['selected_reference_ligand_key'].astype(str).str.len().gt(0) & atom_counts.gt(0)
        values[pdb_id] = bool(usable.any())
    return pd.Series(values, name='reference_evaluable', dtype=bool)

def map_md_evidence(rankings: pd.DataFrame, md_runs: pd.DataFrame, rules: MDRules=MDRules()) -> pd.DataFrame:
    """Map persistent MD contact regions to immutable static candidates."""
    ranking_required = ('pdb_id', 'cluster_v2_id', 'Within_PDB_rank', 'envelope_residues', 'medoid_center_x', 'medoid_center_y', 'medoid_center_z')
    _require_columns(rankings, ranking_required, 'rankings')
    _require_columns(md_runs, MD_INPUT_COLUMNS, 'MD evidence')
    _unique_keys(rankings, ('pdb_id', 'cluster_v2_id'), 'rankings')
    work, source = rankings.copy(), md_runs.copy()
    work['pdb_id'] = work['pdb_id'].astype(str).str.upper()
    source['pdb_id'] = source['pdb_id'].astype(str).str.upper()
    _unique_keys(source, ('pdb_id', 'Simulation_context', 'MD_run'), 'MD evidence')
    targets = set(work['pdb_id'])
    if not set(source['pdb_id']).issubset(targets):
        raise ValueError('MD evidence contains a target absent from rankings')
    rows: list[dict[str, object]] = []
    for md in source.sort_values(['pdb_id', 'Simulation_context', 'MD_run'], kind='mergesort').to_dict('records'):
        pdb_id = str(md['pdb_id'])
        target = work.loc[work['pdb_id'].eq(pdb_id)].copy()
        top = target.sort_values(['Within_PDB_rank', 'cluster_v2_id'], kind='mergesort').iloc[0]
        persistent, md_center = parse_residues(md['Persistent_MD_contact_residues']), _parse_center(md['MD_contact_center'])
        comparisons: list[dict[str, object]] = []
        for cluster in target.to_dict('records'):
            overlap, precision, recall, iou = _overlap(parse_residues(cluster['envelope_residues']), persistent)
            center = _center_from_row(cluster)
            distance = euclidean_distance(center, md_center) if center is not None and md_center is not None else math.nan
            comparisons.append({'cluster_v2_id': cluster['cluster_v2_id'], 'overlap': overlap, 'precision': precision, 'recall': recall, 'iou': iou, 'center_distance': distance})
        comparison = pd.DataFrame(comparisons)
        if persistent:
            ordered = comparison.assign(sort_iou=comparison['iou'].fillna(-1), sort_distance=comparison['center_distance'].fillna(np.inf))
            best = ordered.sort_values(['sort_iou', 'overlap', 'sort_distance', 'cluster_v2_id'], ascending=[False, False, True, True], kind='mergesort').iloc[0]
            static = comparison.loc[comparison['cluster_v2_id'].eq(top['cluster_v2_id'])].iloc[0]
            siou = _number(static['iou'], 'static MD IoU', missing_ok=True)
            spre = _number(static['precision'], 'static MD precision', missing_ok=True)
            sdist = _number(static['center_distance'], 'static MD center distance', missing_ok=True)
            biou = _number(best['iou'], 'best MD IoU', missing_ok=True)
            bdist = _number(best['center_distance'], 'best MD center distance', missing_ok=True)
            if siou >= rules.concordant_iou_min or spre >= rules.concordant_precision_min and sdist <= rules.concordant_center_max_A:
                call = 'concordant'
            elif siou >= rules.partial_iou_min:
                call = 'partially_concordant'
            elif str(best['cluster_v2_id']) != str(top['cluster_v2_id']) and (biou >= rules.alternative_iou_min or bdist <= rules.alternative_center_max_A):
                call = 'static_dynamic_conflict'
            elif sdist <= rules.boundary_shift_center_max_A:
                call = 'boundary_shift'
            else:
                call = 'static_dynamic_conflict'
        else:
            best, static = pd.Series(dtype=object), pd.Series(dtype=object)
            call = 'apo_only_context' if str(md['Simulation_context']).strip().casefold() == 'apo' else 'insufficient_MD_evidence'
        rows.append({
            'pdb_id': pdb_id, 'MD_run': md['MD_run'], 'Simulation_context': md['Simulation_context'],
            'static_top_cluster_v2_id': top['cluster_v2_id'], 'static_top_rank': _number(top['Within_PDB_rank'], 'static_top_rank'),
            'persistent_MD_contact_residue_count': len(simple_set(persistent)),
            'best_MD_mapped_cluster_v2_id': best.get('cluster_v2_id', ''),
            'best_MD_cluster_Jaccard': _number(best.get('iou'), 'best_MD_cluster_Jaccard', missing_ok=True),
            'best_MD_cluster_precision': _number(best.get('precision'), 'best_MD_cluster_precision', missing_ok=True),
            'best_MD_cluster_recall': _number(best.get('recall'), 'best_MD_cluster_recall', missing_ok=True),
            'best_MD_cluster_center_distance_A': _number(best.get('center_distance'), 'best_MD_cluster_center_distance_A', missing_ok=True),
            'static_top_MD_overlap_count': _number(static.get('overlap'), 'static_top_MD_overlap_count', missing_ok=True),
            'static_top_MD_cluster_coverage': _number(static.get('precision'), 'static_top_MD_cluster_coverage', missing_ok=True),
            'static_top_MD_contact_coverage': _number(static.get('recall'), 'static_top_MD_contact_coverage', missing_ok=True),
            'static_top_MD_Jaccard': _number(static.get('iou'), 'static_top_MD_Jaccard', missing_ok=True),
            'static_top_MD_center_distance_A': _number(static.get('center_distance'), 'static_top_MD_center_distance_A', missing_ok=True),
            'D_dyn_run_score': _number(md['D_dyn_run_score'], 'D_dyn_run_score', missing_ok=True),
            'Concordance_call': call, 'MD_mapping_rule_version': 'cluster_v2_posthoc_MD_mapping_20260711',
            'source_persistent_contacts': 'previously_chain_reconciled_MD_contact_region',
        })
    present = {str(row['pdb_id']) for row in rows}
    for pdb_id in sorted(targets - present):
        top = work.loc[work['pdb_id'].eq(pdb_id)].sort_values(['Within_PDB_rank', 'cluster_v2_id'], kind='mergesort').iloc[0]
        rows.append({
            'pdb_id': pdb_id, 'MD_run': '', 'Simulation_context': 'not_available',
            'static_top_cluster_v2_id': top['cluster_v2_id'], 'static_top_rank': top['Within_PDB_rank'],
            'persistent_MD_contact_residue_count': 0, 'best_MD_mapped_cluster_v2_id': '',
            'best_MD_cluster_Jaccard': math.nan, 'best_MD_cluster_precision': math.nan,
            'best_MD_cluster_recall': math.nan, 'best_MD_cluster_center_distance_A': math.nan,
            'static_top_MD_overlap_count': math.nan, 'static_top_MD_cluster_coverage': math.nan,
            'static_top_MD_contact_coverage': math.nan, 'static_top_MD_Jaccard': math.nan,
            'static_top_MD_center_distance_A': math.nan, 'D_dyn_run_score': math.nan,
            'Concordance_call': 'MD_not_available',
            'MD_mapping_rule_version': 'cluster_v2_posthoc_MD_mapping_20260711',
            'source_persistent_contacts': 'not_available',
        })
    return pd.DataFrame(rows, columns=MD_COLUMNS).sort_values(
        ['pdb_id', 'Simulation_context', 'MD_run'], kind='mergesort'
    ).reset_index(drop=True)

def assign_evidence_labels(rankings: pd.DataFrame, reference_mapping: pd.DataFrame, topk_qc: pd.DataFrame, md_mapping: pd.DataFrame, *, rules: MDRules=MDRules()) -> pd.DataFrame:
    """Assign R/A/U/X evidence states without altering static ranking."""
    required = ('pdb_id', 'cluster_v2_id', 'Within_PDB_rank', 'tool_support_count', 'interface_fraction', 'distance_to_interface_A', 'spatial_continuity', 'mappability')
    _require_columns(rankings, required, 'rankings')
    _require_columns(reference_mapping, REFERENCE_COLUMNS, 'reference mapping')
    _require_columns(topk_qc, TOPK_QC_COLUMNS, 'top-k QC')
    _require_columns(md_mapping, MD_COLUMNS, 'MD mapping')
    _unique_keys(reference_mapping, ('pdb_id', 'cluster_v2_id'), 'reference mapping')
    ref = reference_mapping.set_index(['pdb_id', 'cluster_v2_id'])
    qc = topk_qc.set_index(['pdb_id', 'cluster_v2_id'])
    dynamic: defaultdict[tuple[str, str], list[str]] = defaultdict(list)
    for row in md_mapping.to_dict('records'):
        cluster_id = str(row['best_MD_mapped_cluster_v2_id'])
        if cluster_id and cluster_id.casefold() != 'nan' and (_number(row['persistent_MD_contact_residue_count'], 'persistent_MD_contact_residue_count') > 0) and (_number(row['best_MD_cluster_Jaccard'], 'best_MD_cluster_Jaccard', missing_ok=True) >= rules.alternative_iou_min or _number(row['best_MD_cluster_center_distance_A'], 'best_MD_cluster_center_distance_A', missing_ok=True) <= rules.alternative_center_max_A):
            dynamic[str(row['pdb_id']), cluster_id].append(str(row['MD_run']))
    rows: list[dict[str, object]] = []
    for row in rankings.sort_values(['pdb_id', 'Within_PDB_rank', 'cluster_v2_id'], kind='mergesort').to_dict('records'):
        key = (str(row['pdb_id']), str(row['cluster_v2_id']))
        if key not in ref.index:
            raise ValueError(f'reference mapping is missing candidate {key[1]}')
        rrow = ref.loc[key]
        qrow = qc.loc[key] if key in qc.index else None
        clear = _boolean(qrow['clear_exclusion_flag'], 'clear_exclusion_flag') if qrow is not None else False
        qc_status = str(qrow['QC_status']) if qrow is not None else 'not_in_static_top3'
        unresolved_ref = _boolean(rrow['reference_selection_unresolved'], 'reference_selection_unresolved')
        r_auto = _boolean(rrow['R_auto_rule_pass'], 'R_auto_rule_pass') and (not clear) and (not unresolved_ref)
        spatial = _boolean(row['spatial_continuity'], 'spatial_continuity')
        interface = _number(row['interface_fraction'], 'interface_fraction', missing_ok=True)
        interface_distance = _number(row['distance_to_interface_A'], 'distance_to_interface_A', missing_ok=True)
        interface_support = not clear and _integer(row['tool_support_count'], 'tool_support_count') >= 2 and (interface >= 0.35) and (interface_distance <= 8.0) and spatial
        md_runs = sorted(dynamic.get(key, []))
        independent: list[str] = []
        if interface_support:
            independent.append(f'assembly_interface_support:overlap={interface:.3f},distance={interface_distance:.2f}A')
        if md_runs:
            independent.append('MD_contact_region_support:' + ';'.join(md_runs))
        if clear or str(row['mappability']) not in VALID_MAPPABILITY:
            label, reason = ('X_auto', 'automated_QC_clear_exclusion')
        elif r_auto:
            label = 'R_auto'
            reason = f"reference_rule_pass:DCC={_number(rrow['DCC_A'], 'DCC_A', missing_ok=True):.2f}A,contact_recall={_number(rrow['contact_recall'], 'contact_recall', missing_ok=True):.3f},IoU={_number(rrow['residue_IoU'], 'residue_IoU', missing_ok=True):.3f}"
        elif independent:
            label, reason = ('A_auto', ';'.join(independent))
        else:
            label, reason = ('U_auto', 'acceptable_static_candidate_without_sufficient_R_auto_or_A_auto_evidence')
        unresolved = label == 'U_auto' or unresolved_ref or qc_status in {'QC_boundary_sensitive', 'QC_possible_overmerge', 'QC_possible_split', 'QC_unmappable', 'QC_insufficient_evidence'}
        rank = _number(row['Within_PDB_rank'], 'Within_PDB_rank')
        rows.append({'pdb_id': key[0], 'cluster_v2_id': key[1], 'Within_PDB_rank': rank, 'automated_evidence_label': label, 'label_reason': reason, 'independent_evidence_details': ';'.join(independent), 'R_auto_DCC_A': rrow['DCC_A'], 'R_auto_contact_recall': rrow['contact_recall'], 'R_auto_residue_IoU': rrow['residue_IoU'], 'QC_status': qc_status, 'unresolved_flag': unresolved, 'label_scope': 'static_top3' if rank <= 3 else 'non_top3_posthoc_mapping', 'interpretation_guardrail': 'reproducible_evidence_state_not_manual_or_experimental_confirmation'})
    return pd.DataFrame(rows, columns=LABEL_COLUMNS).reset_index(drop=True)

def map_redocking_evidence(labels: pd.DataFrame, redocking: pd.DataFrame) -> pd.DataFrame:
    """Attach redocking evidence using categories recomputed from numeric RMSD."""
    _require_columns(labels, ('pdb_id', 'cluster_v2_id', 'Within_PDB_rank', 'automated_evidence_label'), 'labels')
    required = ('pdb_id', 'Docking_software', 'GlideScore_kcal_per_mol', 'Raw_ligand_RMSD_A', 'RMSD_threshold_call', 'Failure_or_warning_reason')
    _require_columns(redocking, required, 'redocking evidence')
    _unique_keys(redocking, ('pdb_id',), 'redocking evidence')
    first = labels.loc[labels['automated_evidence_label'].eq('R_auto')].sort_values(['pdb_id', 'Within_PDB_rank', 'cluster_v2_id'], kind='mergesort').groupby('pdb_id', sort=True).first()
    rows: list[dict[str, object]] = []
    for row in redocking.sort_values('pdb_id', kind='mergesort').to_dict('records'):
        pdb_id = str(row['pdb_id']).upper()
        rmsd = _number(row['Raw_ligand_RMSD_A'], 'Raw_ligand_RMSD_A')
        if rmsd < 0:
            raise ValueError('Raw_ligand_RMSD_A must be nonnegative')
        if rmsd <= 2.0:
            call, recovered = ('RMSD <= 2 A', True)
        elif rmsd <= 3.0:
            call, recovered = ('2 A < RMSD <= 3 A', False)
        else:
            call, recovered = ('RMSD > 3 A', False)
        rows.append({'pdb_id': pdb_id, 'reference_cluster_v2_id': first.loc[pdb_id, 'cluster_v2_id'] if pdb_id in first.index else '', 'reference_cluster_static_rank': first.loc[pdb_id, 'Within_PDB_rank'] if pdb_id in first.index else math.nan, 'Raw_ligand_RMSD_A': rmsd, 'GlideScore_kcal_per_mol': _number(row['GlideScore_kcal_per_mol'], 'GlideScore_kcal_per_mol', missing_ok=True), 'RMSD_threshold_call': call, 'Reference_pose_recovered': recovered, 'Failure_or_warning_reason': row['Failure_or_warning_reason'], 'redocking_role': 'posthoc_reference_region_chemical_validation_not_static_ranking'})
    return pd.DataFrame(rows, columns=REDOCKING_COLUMNS).sort_values('pdb_id', kind='mergesort').reset_index(drop=True)

def build_unresolved_cases(labels: pd.DataFrame, topk_qc: pd.DataFrame, md_mapping: pd.DataFrame, systems: pd.DataFrame) -> pd.DataFrame:
    """Build the explicit unresolved-evidence ledger without Task 6 metrics."""
    _require_columns(labels, LABEL_COLUMNS, 'labels')
    _require_columns(topk_qc, TOPK_QC_COLUMNS, 'top-k QC')
    _require_columns(md_mapping, MD_COLUMNS, 'MD mapping')
    _require_columns(systems, ('pdb_id', 'pocket_category', 'category_mapping_status'), 'systems')
    _unique_keys(systems, ('pdb_id',), 'systems')
    rows: list[dict[str, object]] = []

    def add(level: str, pdb_id: str, cluster_id: str, issue_type: str, status: str, required: str, handling: str) -> None:
        rows.append({'issue_id': f'U{len(rows) + 1:03d}', 'level': level, 'pdb_id': pdb_id, 'cluster_v2_id': cluster_id, 'issue_type': issue_type, 'status': status, 'required_user_or_future_input': required, 'current_handling': handling})
    top = labels.loc[labels['Within_PDB_rank'].le(3)].sort_values(['pdb_id', 'Within_PDB_rank', 'cluster_v2_id'], kind='mergesort')
    for row in top.to_dict('records'):
        if _boolean(row['unresolved_flag'], 'unresolved_flag'):
            add('candidate', row['pdb_id'], row['cluster_v2_id'], 'automated_label_or_boundary_unresolved', f"{row['automated_evidence_label']}|{row['QC_status']}", 'Optional targeted structural inspection or additional orthogonal evidence', 'Retained with unresolved_flag; not forced into a binary biological assignment')
    for pdb_id in sorted(set(md_mapping.loc[md_mapping['Concordance_call'].eq('MD_not_available'), 'pdb_id'])):
        add('target', pdb_id, '', 'MD_not_available', 'open', 'Ligand-containing or otherwise hypothesis-matched MD trajectory if dynamic refinement is desired', 'D_dyn omitted; static priority remains unchanged')
    counts = labels.loc[labels['automated_evidence_label'].eq('R_auto')].groupby('pdb_id').size()
    for pdb_id, count in counts.loc[counts.gt(1)].sort_index().items():
        add('target', pdb_id, '', 'reference_region_split_across_multiple_cluster_v2_boundaries', f'{int(count)}_R_auto_clusters', 'Optional boundary-focused inspection; no input is required for current metrics', 'Reference metric uses the first R_auto rank and preserves all fragments')
    system_work = systems.copy()
    system_work['pdb_id'] = system_work['pdb_id'].astype(str).str.upper()
    for pdb_id in sorted(system_work['pdb_id']):
        add('target', pdb_id, '', 'formal_interface_energetics_not_available', 'PISA_or_equivalent_fields_NA', 'PISA interface area, buried SASA, DeltaG, hydrogen bonds, salt bridges, and hotspot annotation', 'O_rel uses only the declared 5 A heavy-atom interface proxy; missing fields are not treated as zero')
        add('target', pdb_id, '', 'crystal_contact_risk_not_evaluable', 'no_generated_crystal_mates', 'Symmetry-expanded crystal-contact analysis if this question is important for the target', 'No positive or negative crystal-contact claim is made')
    provisional = system_work.loc[system_work['category_mapping_status'].astype(str).str.contains('provisional', case=False, na=False)].sort_values('pdb_id', kind='mergesort')
    for row in provisional.to_dict('records'):
        add('target', row['pdb_id'], '', 'pocket_archetype_assignment_provisional', str(row['pocket_category']), 'Author confirmation of the four-class pocket archetype', 'Category-stratified results remain descriptive')
    return pd.DataFrame(rows, columns=UNRESOLVED_COLUMNS)

def _config_data(config: object) -> Mapping[str, object]:
    data = getattr(config, 'data', config)
    if not isinstance(data, Mapping):
        raise TypeError('config must be a mapping or ManuscriptConfig')
    return data

def representative_case_ids(config: object) -> tuple[str, ...]:
    analysis = _config_data(config).get('analysis')
    if not isinstance(analysis, Mapping):
        raise ValueError('configuration analysis section is missing')
    values = analysis.get('representative_case_ids')
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise ValueError('representative_case_ids must be a sequence')
    result = tuple((str(value).strip().upper() for value in values))
    if not result or any((not value for value in result)) or len(set(result)) != len(result):
        raise ValueError('representative_case_ids must be nonempty unique target IDs')
    return result

def _rules_from_config(config: object) -> tuple[int, QCRules, ReferenceRules, MDRules]:
    posthoc = _config_data(config).get('posthoc')
    if not isinstance(posthoc, Mapping):
        raise ValueError('configuration posthoc section is missing')
    qc, reference, md = (posthoc.get(name) for name in ('top3_qc', 'reference', 'md'))
    if not all((isinstance(value, Mapping) for value in (qc, reference, md))):
        raise ValueError('configuration posthoc rule sections are missing')
    assert isinstance(qc, Mapping) and isinstance(reference, Mapping) and isinstance(md, Mapping)

    def values(section: Mapping[str, object], cls):
        return cls(**{name: _integer(section[name], name) if isinstance(default, int) else _number(section[name], name) for name, default in ((field, getattr(cls(), field)) for field in cls.__dataclass_fields__)})
    return (_integer(qc['top_k'], 'top_k'), values(qc, QCRules), values(reference, ReferenceRules), values(md, MDRules))

def _validate_inputs(inputs: PosthocInputs) -> None:
    rankings = inputs.rankings
    if rankings.columns.tolist() != list(RANKING_COLUMNS) or rankings.empty:
        raise ValueError('post-hoc analysis requires a nonempty exact static ranking schema')
    _unique_keys(rankings, ('pdb_id', 'cluster_v2_id'), 'rankings')
    target_ids = set(rankings['pdb_id'].astype(str))
    if not target_ids or any((not re.fullmatch(r'[0-9][A-Z0-9]{3}', value) for value in target_ids)):
        raise ValueError('static rankings must contain uppercase PDB target IDs')
    expected_order = rankings.sort_values(
        ['pdb_id', 'Within_PDB_rank', 'cluster_v2_id'], kind='mergesort'
    ).reset_index(drop=True)
    if not rankings.reset_index(drop=True)[['pdb_id', 'Within_PDB_rank', 'cluster_v2_id']].equals(
        expected_order[['pdb_id', 'Within_PDB_rank', 'cluster_v2_id']]
    ):
        raise ValueError('static rankings must be stably sorted')
    contracts = ((inputs.references, 'references'), (inputs.redocking, 'redocking'), (inputs.systems, 'systems'))
    for frame, label in contracts:
        _require_columns(frame, ('pdb_id',), label)
        values = frame['pdb_id'].astype(str)
        if len(frame) != len(target_ids) or values.nunique() != len(target_ids):
            raise ValueError(f'{label} must contain one unique row per ranked target')
        if set(values) != target_ids:
            raise ValueError(f'{label} target coverage does not match static rankings')
    if set(inputs.structures) != target_ids:
        raise ValueError('prepared structures must exactly cover ranked targets')
    if any((not path.is_file() or path.is_symlink() for path in inputs.structures.values())):
        raise FileNotFoundError('one or more prepared structures are missing')
    _require_columns(inputs.md_runs, MD_INPUT_COLUMNS, 'MD evidence')
    md_runs = inputs.md_runs.copy()
    original_md_ids = md_runs['pdb_id'].astype(str)
    md_runs['pdb_id'] = original_md_ids.str.upper()
    _unique_keys(md_runs, ('pdb_id', 'Simulation_context', 'MD_run'), 'MD evidence')
    if not original_md_ids.equals(md_runs['pdb_id']):
        raise ValueError('MD pdb_id values must be uppercase')
    if not set(md_runs['pdb_id']).issubset(target_ids):
        raise ValueError('MD target coverage is outside static rankings')

def _convenience_master(rankings: pd.DataFrame, labels: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    output = rankings.copy(deep=True)
    keys = pd.MultiIndex.from_frame(output[['pdb_id', 'cluster_v2_id']])
    label_lookup = labels.set_index(['pdb_id', 'cluster_v2_id'])
    ref_lookup = reference.set_index(['pdb_id', 'cluster_v2_id'])
    for column in CONVENIENCE_COLUMNS[:4]:
        output[column] = label_lookup.reindex(keys)[column].to_numpy()
    for column in CONVENIENCE_COLUMNS[4:]:
        output[column] = ref_lookup.reindex(keys)[column].to_numpy()
    return output.loc[:, [*RANKING_COLUMNS, *CONVENIENCE_COLUMNS]]

def run_posthoc(inputs: PosthocInputs, config: object) -> PosthocResult:
    """Run the complete deterministic post-hoc mapping without filesystem writes."""
    _validate_inputs(inputs)
    top_k, qc_rules, reference_rules, md_rules = _rules_from_config(config)
    topk_qc = run_topk_qc(inputs.rankings, inputs.structures, top_k=top_k, rules=qc_rules)
    reference = map_reference_evidence(inputs.rankings, inputs.structures, inputs.references, rules=reference_rules)
    md = map_md_evidence(inputs.rankings, inputs.md_runs, rules=md_rules)
    labels = assign_evidence_labels(inputs.rankings, reference, topk_qc, md, rules=md_rules)
    redocking = map_redocking_evidence(labels, inputs.redocking)
    unresolved = build_unresolved_cases(labels, topk_qc, md, inputs.systems)
    cases = representative_case_ids(config)
    ranked_targets = set(inputs.rankings['pdb_id'].astype(str))
    if not set(cases).issubset(ranked_targets):
        raise ValueError('representative case IDs must occur in static rankings')
    return PosthocResult(topk_qc=topk_qc, reference_mapping=reference, reference_evaluable=derive_reference_evaluable(reference), md_mapping=md, evidence_labels=labels, redocking_mapping=redocking, unresolved_cases=unresolved, convenience_master=_convenience_master(inputs.rankings, labels, reference), representative_case_ids=cases)
