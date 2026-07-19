# Release data dictionary

This dictionary covers every CSV in the canonical frozen release under `results/reference/`. Column order is normative and comes from `io.py`, `posthoc.py`, `statistics.py`, `figures.py`, and `validation.py`. Scientific definitions are in [oips-formula.md](oips-formula.md).

## Serialization and missingness

All release CSVs use UTF-8, LF line endings, `%.15g` float formatting, lowercase `true`/`false`, explicit column order, and stable sorting. An empty numeric cell means unavailable (`NA`); it does not mean zero. Text fields may be empty only when the producing contract permits it, for example an unavailable MD run identifier. Semicolon-delimited fields are display/audit lists inside one CSV cell, not additional rows.

Suffixes have consistent meanings: `_A` is an Ångström quantity; `_count` and `_N` are counts; `_fraction`, `_precision`, `_recall`, `_IoU`, proportions, MRR, and Top-k summary values are normally in `[0,1]`; `_flag`, `_evaluable`, `_changed`, and `_recovered` fields are booleans. Exceptions and exact missing behavior are stated below.

## Canonical CSV inventory

### Clustering: five CSVs

| Canonical path | Grain and key | Purpose | Ordered columns |
|---|---|---|---|
| `results/reference/clustering/cluster_v2_candidates.csv` | One row per candidate; key `pdb_id,cluster_v2_id`. | Reference-blind cluster geometry, consensus, mappability, and audit state before scoring. | `pdb_id, cluster_v2_id, medoid_unit_id, medoid_tool, medoid_pocket_id, medoid_center_x, medoid_center_y, medoid_center_z, cluster_diameter_A, center_dispersion_A, tool_support_count, supporting_tools, representative_pockets_per_tool, same_tool_units, same_tool_secondary_unit_count, raw_record_count, core_residue_count, envelope_residue_count, core_envelope_ratio, core_residues, envelope_residues, contributing_chains, contributing_chain_count, dominant_chain_fraction, cluster_chain_entropy, mappability, center_available_representatives, residue_available_representatives, pairwise_residue_iou_median, pairwise_residue_iou_min, spatial_continuity, boundary_sensitive` |
| `results/reference/clustering/cluster_v2_membership.csv` | One raw record assignment within a candidate; stable key fields `pdb_id,cluster_v2_id,tool,raw_row_id`. | Full raw-record membership and the one-formal-vote-per-tool audit. | `pdb_id, cluster_v2_id, tool, same_tool_unit_id, raw_row_id, raw_pocket_id, formal_tool_representative, representative_pocket_id, formal_vote_count` |
| `results/reference/clustering/tool_record_to_cluster_v2_mapping.csv` | One row per source `row_id`; stable order `pdb_id,tool,row_id`. | Exhaustive mapping of mapped and excluded source records. | `row_id, pdb_id, tool, pocket_id, center_x, center_y, center_z, residue_count, same_tool_unit_id, same_tool_group_size, representative_for_tool_unit, representative_pocket_id, cluster_v2_id, mapping_status, exclusion_reason` |
| `results/reference/clustering/excluded_unmappable_records.csv` | One excluded record; key `row_id`. | Audit of records lacking both a complete center and a residue set. | `row_id, pdb_id, tool, pocket_id, center_method, residue_count, exclusion_reason, retained_in_audit` |
| `results/reference/clustering/cluster_v2_boundary_audit.csv` | One row per candidate; key `pdb_id,cluster_v2_id`. | Compact boundary and spatial-continuity audit. | `pdb_id, cluster_v2_id, tool_support_count, same_tool_secondary_unit_count, cluster_diameter_A, center_dispersion_A, pairwise_residue_iou_median, pairwise_residue_iou_min, core_envelope_ratio, spatial_continuity, boundary_sensitive, boundary_reason` |

### Static scoring: two CSVs

| Canonical path | Grain and key | Purpose | Ordered columns |
|---|---|---|---|
| `results/reference/static/cluster_v2_master_table.csv` | One row per candidate; key `pdb_id,cluster_v2_id`. | Candidate fields plus the five static modules and `OIPS-P_static`; contains no reference, MD, redocking, or evidence-label fields. | `pdb_id, cluster_v2_id, medoid_unit_id, medoid_tool, medoid_pocket_id, medoid_center_x, medoid_center_y, medoid_center_z, cluster_diameter_A, center_dispersion_A, tool_support_count, supporting_tools, representative_pockets_per_tool, same_tool_units, same_tool_secondary_unit_count, raw_record_count, core_residue_count, envelope_residue_count, core_envelope_ratio, core_residues, envelope_residues, contributing_chains, contributing_chain_count, dominant_chain_fraction, cluster_chain_entropy, mappability, center_available_representatives, residue_available_representatives, pairwise_residue_iou_median, pairwise_residue_iou_min, spatial_continuity, C_cons, G_geo, P_lig, Q_evidence, O_rel_formal, interface_fraction, interface_recall, cluster_interface_residue_count, cluster_chain_count, distance_to_interface_A, interface_distance_score, boundary_sensitive, OIPS-P_static` |
| `results/reference/static/cluster_v2_static_rankings.csv` | One row per candidate; key `pdb_id,cluster_v2_id`, sorted by `pdb_id,Within_PDB_rank,cluster_v2_id`. | Immutable static master plus formula recomputation and average within-target rank. | All columns of `cluster_v2_master_table.csv`, followed by `OIPS-P_static_recomputed, Within_PDB_rank, tie_flag, tie_size`. |

### Post-hoc analysis: twenty CSVs

| Canonical path | Grain and key | Purpose | Ordered columns |
|---|---|---|---|
| `results/reference/analysis/final_top3_automated_QC.csv` | One static top-3 candidate; key `pdb_id,cluster_v2_id`. | Automated structural QC without reranking. | `pdb_id, cluster_v2_id, Within_PDB_rank, OIPS-P_static, mappability, center_available, residue_set_available, cluster_diameter_A, spatial_continuity, tool_support_count, center_dispersion_A, core_envelope_ratio, contributing_chain_count, dominant_chain_fraction, pocket_interface_overlap, distance_to_interface_A, possible_over_merging, clear_exclusion_flag, possible_split_subpocket, possible_surface_noise_cluster, nearest_missing_residue_distance_sequence_positions, missing_residue_proximity_flag, missing_residue_assessment_source, possible_crystal_contact_risk, boundary_sensitive, QC_status, QC_rule_version` |
| `results/reference/analysis/final_reference_mapping.csv` | One candidate-reference comparison; key `pdb_id,cluster_v2_id`. | Ligand selection, DCC/contact overlap, and inclusive `R_auto` predicate. | `pdb_id, cluster_v2_id, Within_PDB_rank, DCC_A, reference_contact_residue_count, reference_contact_overlap_count, contact_precision, contact_recall, residue_IoU, functional_residue_overlap, functional_residue_status, interface_overlap, distance_to_interface_A, chain_contribution_count, R_auto_rule_pass, R_auto_rule, reference_ligand_annotation, requested_reference_ligand_key, selected_reference_ligand_key, reference_ligand_atom_count, reference_selection_status, reference_selection_unresolved` |
| `results/reference/analysis/final_md_cluster_v2_mapping.csv` | One real MD run or one target-level unavailable row; key `pdb_id,Simulation_context,MD_run`. | Persistent MD region mapping to the immutable static candidates. | `pdb_id, MD_run, Simulation_context, static_top_cluster_v2_id, static_top_rank, persistent_MD_contact_residue_count, best_MD_mapped_cluster_v2_id, best_MD_cluster_Jaccard, best_MD_cluster_precision, best_MD_cluster_recall, best_MD_cluster_center_distance_A, static_top_MD_overlap_count, static_top_MD_cluster_coverage, static_top_MD_contact_coverage, static_top_MD_Jaccard, static_top_MD_center_distance_A, D_dyn_run_score, Concordance_call, MD_mapping_rule_version, source_persistent_contacts` |
| `results/reference/analysis/final_automated_evidence_labels.csv` | One candidate; key `pdb_id,cluster_v2_id`. | Reproducible `R_auto/A_auto/U_auto/X_auto` evidence states. | `pdb_id, cluster_v2_id, Within_PDB_rank, automated_evidence_label, label_reason, independent_evidence_details, R_auto_DCC_A, R_auto_contact_recall, R_auto_residue_IoU, QC_status, unresolved_flag, label_scope, interpretation_guardrail` |
| `results/reference/analysis/final_redocking_cluster_v2_mapping.csv` | One target; key `pdb_id`. | Numeric RMSD categories attached to the first reference-associated cluster. | `pdb_id, reference_cluster_v2_id, reference_cluster_static_rank, Raw_ligand_RMSD_A, GlideScore_kcal_per_mol, RMSD_threshold_call, Reference_pose_recovered, Failure_or_warning_reason, redocking_role` |
| `results/reference/analysis/unresolved_cases.csv` | One issue; key `issue_id`. | Explicit ledger of unresolved evidence, metadata, and interpretation needs. | `issue_id, level, pdb_id, cluster_v2_id, issue_type, status, required_user_or_future_input, current_handling` |
| `results/reference/analysis/final_cluster_v2_master_table.csv` | One candidate; key `pdb_id,cluster_v2_id`. | Convenience join only; it does not replace or mutate the static ranking. | All ordered columns of `cluster_v2_static_rankings.csv`, followed by `automated_evidence_label, label_reason, unresolved_flag, QC_status, DCC_A, contact_recall, contact_precision, residue_IoU, R_auto_rule_pass`. |
| `results/reference/analysis/target_level_candidate_prioritization.csv` | One target; key `pdb_id`. | Reference and first-supported target endpoints. | `pdb_id, Pocket_category, Protein_family, formal_cluster_v2_count, reference_evaluable, R_auto_cluster_count, reference_first_rank, reference_rank_percentile, reference_top1, reference_top3, reference_top5, reference_reciprocal_rank, first_supported_evaluable, first_supported_rank, first_supported_top1, first_supported_top3, first_supported_reciprocal_rank, static_top_cluster_v2_id, static_top_evidence_label` |
| `results/reference/analysis/final_candidate_prioritization_metrics.csv` | One overall row; key `analysis_level,group`. | Overall endpoint summary. | `analysis_level, group, reference_evaluable_N, Reference_Top1, Reference_Top3, Reference_Top5, Reference_MRR, Reference_median_rank, Reference_median_rank_percentile, first_supported_evaluable_N, First_supported_Top1, First_supported_Top3, First_supported_MRR` |
| `results/reference/analysis/final_category_metrics.csv` | One pocket category; key `analysis_level,group`. | The same endpoint summary stratified by pocket archetype. | Same ordered columns as `final_candidate_prioritization_metrics.csv`. |
| `results/reference/analysis/final_bootstrap_intervals.csv` | One bootstrap method and metric; key `bootstrap_method,metric`. | Target-resampled and family-clustered percentile intervals. | `bootstrap_method, metric, point_estimate, CI_2.5_percent, CI_97.5_percent, iterations, random_seed` |
| `results/reference/analysis/final_family_sensitivity.csv` | One excluded protein family; key `excluded_family`. | Leave-one-family-out endpoint sensitivity. | `excluded_family, excluded_target_N, remaining_target_N, reference_evaluable_N, Reference_Top1, Reference_Top3, Reference_Top5, Reference_MRR, Reference_median_rank, Reference_median_rank_percentile, first_supported_evaluable_N, First_supported_Top1, First_supported_Top3, First_supported_MRR` |
| `results/reference/analysis/final_orel_ablation_targets.csv` | One target; key `pdb_id`. | Per-target full-versus-without-`O_rel_formal` comparison. | `analysis_level, pdb_id, Pocket_category, Protein_family, Full_reference_rank, Without_O_rel_reference_rank, reference_rank_change_without_minus_full, Full_first_supported_rank, Without_O_rel_first_supported_rank, first_supported_rank_change_without_minus_full, Full_reference_Top3, Without_O_rel_reference_Top3, reference_Top3_inclusion_change, rank_Spearman_rho, Full_top_cluster_v2_id, Without_O_rel_top_cluster_v2_id, top_cluster_identity_changed` |
| `results/reference/analysis/final_orel_ablation_categories.csv` | One pocket category; key `Pocket_category`. | Category-level ablation endpoints. | `analysis_level, Pocket_category, target_N, Full_reference_Top1, Without_O_rel_reference_Top1, Full_reference_Top3, Without_O_rel_reference_Top3, Full_reference_MRR, Without_O_rel_reference_MRR, mean_rank_Spearman_rho, top_cluster_identity_change_fraction` |
| `results/reference/analysis/orel_ablation_cluster_rankings.csv` | One candidate; key `pdb_id,cluster_v2_id`. | Candidate-level ablated scores and ranks. | `pdb_id, cluster_v2_id, Within_PDB_rank, Without_O_rel_score, Without_O_rel_rank` |
| `results/reference/analysis/representative_case_results.csv` | One configured case; key `pdb_id`. | Compact source table for 5J89, 5TBM, and 4W9H case panels. | `pdb_id, static_top_cluster_v2_id, static_top_score, static_top_label, static_top_tool_support, static_top_chains, static_top_interface_overlap, reference_cluster_v2_id, reference_rank, reference_DCC_A, without_O_rel_reference_rank, without_O_rel_top_cluster_v2_id, MD_calls, redocking_RMSD_A, redocking_status` |
| `results/reference/analysis/weight_sensitivity_scenarios.csv` | One perturbed module and direction; key `scenario_id`. | Ten one-at-a-time ±20% weight scenarios and their rank-stability/end-point summaries. | `scenario_id, perturbed_module, direction, multiplier, target_N, baseline_top1_retained_N, mean_top3_jaccard, median_spearman_rho, Reference_Top1_N, Reference_Top3_N, First_supported_Top1_N, First_supported_Top3_N` |
| `results/reference/analysis/weight_sensitivity_targets.csv` | One target; key `pdb_id`. | Across-scenario rank ranges and Top-1 identity stability. | `pdb_id, baseline_top_cluster_v2_id, top1_retention_count, top1_changed_scenarios, baseline_reference_rank, minimum_reference_rank, maximum_reference_rank, baseline_first_supported_rank, minimum_first_supported_rank, maximum_first_supported_rank` |
| `results/reference/analysis/single_tool_target_ranks.csv` | One target and method; key `pdb_id,method`. | Native regional counts and first `R_auto` rank for five tools and OIPS-P; complete-case membership is explicit. | `pdb_id, method, output_status, mappable_region_count, native_top_cluster_v2_id, first_reference_associated_rank, complete_case` |
| `results/reference/analysis/single_tool_complete_case_metrics.csv` | One method; key `method`. | Common 14-target Top-1/Top-3/Top-5/MRR comparison. | `method, N, Top1_N, Top1, Top3_N, Top3, Top5_N, Top5, MRR` |

### Figure source data: four CSVs

Each file has one row per plotted source record, stable order `panel,record_id,series`, and the same ordered schema: `panel, record_id, series, group, x, y, x_label, y_label, lower, upper, annotation`.

- `results/reference/figure_source_data/repository_summary_figure_1_candidate_landscape_source_data.csv`
- `results/reference/figure_source_data/repository_summary_figure_2_qc_and_orel_ablation_source_data.csv`
- `results/reference/figure_source_data/repository_summary_figure_3_posthoc_evidence_source_data.csv`
- `results/reference/figure_source_data/repository_summary_figure_4_representative_cases_source_data.csv`

## Shared clustering and static columns

Definitions in this section apply everywhere the named column is repeated.

| Column(s) | Type/range and definition | Missing behavior |
|---|---|---|
| `pdb_id` | Uppercase four-character target identifier matching `[0-9][A-Z0-9]{3}`. | Required. |
| `cluster_v2_id` | Stable target-prefixed candidate identifier such as `1A9M_V2C001`. | Empty only for an excluded raw record that has no cluster assignment. |
| `row_id`, `raw_row_id` | Unique public feature-row integer; `raw_row_id` is the same identifier in membership. | Required in their row grains. |
| `tool`, `medoid_tool` | Source detector name; `medoid_tool` is the medoid unit's detector. | Required. |
| `pocket_id`, `raw_pocket_id`, `medoid_pocket_id` | Source pocket identifier; prefixes state whether it is an input, member, or medoid identifier. | Required for mapped/source rows. |
| `medoid_unit_id` | Selected same-tool unit used as candidate medoid. | Required for a candidate. |
| `medoid_center_x`, `medoid_center_y`, `medoid_center_z`; `center_x`, `center_y`, `center_z` | Cartesian center coordinates in Å. Medoid fields describe the candidate; unsuffixed fields in the mapping describe the source record. | A center is usable only when all three coordinates are present; all three may be empty. |
| `center_method` | Released source label for how a record center was obtained. | May be empty only if the curated source supplied an empty label. |
| `cluster_diameter_A` | Maximum Euclidean distance among available centers of all same-tool units in the candidate; nonnegative. | Defined as 0 for zero or one available center. |
| `center_dispersion_A` | Root-mean-square distance of available formal-representative centers from the medoid; nonnegative. | Empty when the medoid center is unavailable. |
| `tool_support_count` | Number of formal tool representatives, a nonnegative integer and 1–5 for released candidates. | Required. |
| `supporting_tools` | Semicolon-separated sorted formal-support tool names. | Required for a candidate. |
| `representative_pockets_per_tool` | Semicolon-separated `tool:pocket` audit mapping for formal representatives. | Required for a candidate. |
| `same_tool_unit_id`, `same_tool_units` | Stable same-tool consolidation identifier; plural field is the semicolon-separated candidate membership list. | Empty in raw mapping only when the record was excluded before unit construction. |
| `same_tool_group_size` | Number of raw records consolidated into the raw record's same-tool unit; positive integer when mapped. | Empty for excluded records. |
| `same_tool_secondary_unit_count` | Candidate units beyond the one formal representative per supported tool; nonnegative integer. | Required for a candidate. |
| `raw_record_count` | Number of source feature rows represented by a candidate; positive integer. | Required. |
| `representative_for_tool_unit`, `formal_tool_representative`, `retained_in_audit` | Booleans marking the unit representative, formal cluster/tool vote, or retained exclusion-audit row. | Required in the corresponding tables. |
| `formal_vote_count` | `1` only for the raw row that is the formal representative of its cluster/tool; otherwise `0`. | Required; cluster/tool sums cannot exceed 1. |
| `representative_pocket_id` | Pocket identifier selected for the same-tool unit or formal tool vote, depending on table. | Empty only for an excluded mapping row. |
| `residue_count` | Source residue count, nonnegative integer. | Required; zero denotes no released residue set. |
| `core_residue_count`, `envelope_residue_count` | Sizes of the formal-representative consensus core and union envelope; nonnegative integers. | Required. |
| `core_residues`, `envelope_residues` | Semicolon-separated canonical full residue identifiers in the core and envelope. | Empty represents an empty set. |
| `core_envelope_ratio` | `core_residue_count/envelope_residue_count`, range `[0,1]`. | Empty for an empty envelope. |
| `contributing_chains`, `contributing_chain_count` | Sorted semicolon list and count of chains represented in the envelope. | Empty list and count 0 for an empty envelope. |
| `dominant_chain_fraction` | Largest chain's fraction of envelope residues, range `(0,1]`. | Empty for an empty envelope. |
| `cluster_chain_entropy` | Normalized Shannon entropy of envelope-residue counts by chain, range `[0,1]`. | Empty in the clustering table for no contributing chains; the static scoring helper defines empty input as 0 internally. |
| `mappability` | `center_and_residue_mappable`, `center_only_mappable`, or `residue_only_mappable`, based on available formal representatives. | Required for a candidate. |
| `center_available_representatives`, `residue_available_representatives` | Counts of formal representatives with complete centers or nonempty residue sets. | Required, nonnegative integers. |
| `pairwise_residue_iou_median`, `pairwise_residue_iou_min` | Median and minimum of available formal-representative pairwise residue IoUs, range `[0,1]`. | Empty if no pair has two nonempty residue sets. |
| `spatial_continuity` | Boolean `cluster_diameter_A <= 12.0`. | Required. |
| `boundary_sensitive` | Boolean strict-threshold audit flag for multi-tool candidates. | Required; false is not proof of biological correctness. |
| `boundary_reason` | Semicolon list of triggered diameter, dispersion, core/envelope, and median-IoU audit rules. | Empty when `boundary_sensitive=false`. |
| `mapping_status` | `mapped_to_cluster_v2` or `excluded_unmappable`. | Required for every source record. |
| `exclusion_reason` | `unmappable_no_center_and_no_residue_set` for the released exclusion path. | Empty for mapped records. |

### Static score and rank columns

| Column(s) | Type/range and definition | Missing behavior |
|---|---|---|
| `C_cons` | Tool-weighted formal consensus score, `[0,100]`. | Required for valid candidates. |
| `G_geo` | Mean of up to the top three available formal-representative geometry values; finite real. | Empty when all geometry values are unavailable. |
| `P_lig` | Mean of up to the top two available formal-representative ligandability values; finite real. | Empty when all ligandability values are unavailable. |
| `Q_evidence` | Bounded static representation-completeness score, `[0,100]`. | Required for a valid nonempty candidate. |
| `O_rel_formal` | Static oligomer-interface relevance score, `[0,100]`. | The distance subterm may be omitted and renormalized; the module remains defined by the declared no-interface branches. |
| `interface_fraction` | Fraction of envelope residues that are static interface residues, `[0,1]`. | Empty for an empty envelope. |
| `interface_recall` | Fraction of static interface residues recovered by the envelope, `[0,1]`. | Empty when the prepared structure has no detected interface. |
| `cluster_interface_residue_count` | Size of the envelope/interface intersection; nonnegative integer. | Required. |
| `cluster_chain_count` | Number of chains in the full envelope residue identifiers used by `O_rel_formal`; nonnegative integer. | Zero for an empty envelope. |
| `distance_to_interface_A` | Minimum medoid-center distance to a static interface atom, nonnegative Å. | Empty if the center or interface coordinates are unavailable. |
| `interface_distance_score` | Configured piecewise transform of `distance_to_interface_A`, normally `[10,100]`. | Empty when distance is unavailable. |
| `OIPS-P_static` | Missing-aware weighted static score. | Empty only if all five modules are unavailable; that state is invalid in the frozen bundle. |
| `OIPS-P_static_recomputed` | Independent recomputation from serialized modules and configured weights. | Must equal `OIPS-P_static` within `1e-10`. |
| `Within_PDB_rank` | Descending average rank within target; range `[1,N]`, potentially fractional. | Required in the frozen ranking. |
| `tie_flag`, `tie_size` | Whether the exact recomputed score is shared and the size of that tie group. | Required; `tie_size>=1`. |

## Post-hoc QC, reference, MD, and label columns

### Top-3 QC

| Column(s) | Definition | Missing behavior or values |
|---|---|---|
| `center_available`, `residue_set_available` | Whether at least one formal representative has a center or residue set. | Required booleans. |
| `pocket_interface_overlap` | Copy of static `interface_fraction`. | Empty when static value is unavailable. |
| `possible_over_merging`, `clear_exclusion_flag`, `possible_split_subpocket`, `possible_surface_noise_cluster` | Deterministic QC booleans from the configured strict rules. | Missing numeric comparisons are false except where the rule explicitly treats missing ligandability/interface as low information. |
| `nearest_missing_residue_distance_sequence_positions` | Minimum absolute sequence-number difference between an envelope residue and an inferred missing residue. | Empty when no comparable missing position exists. |
| `missing_residue_proximity_flag` | True when the preceding distance is at most 2 positions. | False when distance is unavailable. |
| `missing_residue_assessment_source` | `PDB_REMARK_465`, `observed_numbering_gap_proxy`, or `no_gap_detected`. | Required. |
| `possible_crystal_contact_risk` | Explicit non-evaluability state: `not_evaluable_no_generated_crystal_mates` or `not_evaluable_no_crystal_symmetry_model`. | Never interpreted as a positive or negative crystal-contact finding. |
| `QC_status` | `QC_pass`, `QC_boundary_sensitive`, `QC_possible_split`, `QC_possible_overmerge`, `QC_insufficient_evidence`, or `QC_unmappable`. | Required. |
| `QC_rule_version` | Stable producing-rule label `automated_QC_v2_20260711`. | Required. |

### Reference mapping

| Column(s) | Definition | Missing behavior or values |
|---|---|---|
| `DCC_A`, `R_auto_DCC_A`, `reference_DCC_A` | Candidate medoid-to-selected-ligand-centroid distance in Å; prefixed variants are copies in label/case tables. | Empty if either center is unavailable. |
| `reference_contact_residue_count` | Number of standard-amino-acid residues within 4.5 Å of the selected ligand. | Zero when no contacts resolve. |
| `reference_contact_overlap_count` | Candidate-envelope/reference-contact intersection size. | Zero if either set is empty. |
| `contact_precision`, `R_auto_contact_recall`, `contact_recall`, `R_auto_residue_IoU`, `residue_IoU` | Precision, recall, and IoU overlap quantities defined in the formula contract; prefixed fields are copies. | Empty if either compared set is empty. |
| `functional_residue_overlap`, `functional_residue_status` | Reserved functional-residue result and its provenance state. | The released numeric field is empty and status is `not_available_not_imputed`. |
| `interface_overlap`, `chain_contribution_count` | Copies of static interface fraction and contributing-chain count for interpretation. | Interface copy may be empty. |
| `R_auto_rule_pass` | Boolean result of the inclusive DCC/recall/IoU predicate before unresolved-selection and QC label guards. | False when required quantities are unavailable. |
| `R_auto_rule` | Human-readable frozen predicate string. | Required. |
| `reference_ligand_annotation` | Reviewed project ligand key, including the corrected `6VN:D:503:` key for 5LGE. | Required. |
| `requested_reference_ligand_key` | Key actually requested after applying the reviewed override rule. | Required. |
| `selected_reference_ligand_key` | Structure ligand group selected by exact key or declared fallback. | Empty when no ligand is selected. |
| `reference_ligand_atom_count` | Number of atoms in the selected ligand group. | Zero when selection fails. |
| `reference_selection_status` | `exact_project_reference_key`, `unique_resname_fallback`, `ambiguous_resname_largest_group_fallback`, `invalid_project_reference_key`, or `reference_ligand_not_found`. | Required. |
| `reference_selection_unresolved` | True for ambiguous fallback, invalid key, or ligand not found. | Required boolean. |

### MD mapping

| Column(s) | Definition | Missing behavior or values |
|---|---|---|
| `MD_run`, `Simulation_context` | Released run label and context. | `MD_run` is empty and context is `not_available` for a target without MD. |
| `static_top_cluster_v2_id`, `static_top_rank` | Immutable static rank-1 cluster and its rank. | Required. |
| `persistent_MD_contact_residue_count` | Size of the released persistent contact set. | Zero for empty or unavailable contact sets. |
| `best_MD_mapped_cluster_v2_id` | Candidate chosen by IoU, overlap, center distance, and identifier. | Empty when no persistent contacts exist or MD is unavailable. |
| `best_MD_cluster_Jaccard`, `best_MD_cluster_precision`, `best_MD_cluster_recall` | Overlap values for the best MD-mapped candidate. | Empty without a persistent contact set. |
| `best_MD_cluster_center_distance_A` | Best candidate-center to MD-region-center distance, Å. | Empty if either center is unavailable. |
| `static_top_MD_overlap_count` | Static-top envelope/MD-contact intersection size. | Empty when MD contacts are unavailable. |
| `static_top_MD_cluster_coverage`, `static_top_MD_contact_coverage`, `static_top_MD_Jaccard` | Precision, recall, and IoU for the immutable static top candidate. | Empty when either set is empty. |
| `static_top_MD_center_distance_A` | Static-top center to MD-region center distance, Å. | Empty if either center is unavailable. |
| `D_dyn_run_score` | Released precomputed dynamic score for the run; finite real when present. | Empty when unavailable; never enters static scoring. |
| `Concordance_call` | `concordant`, `partially_concordant`, `boundary_shift`, `static_dynamic_conflict`, `apo_only_context`, `insufficient_MD_evidence`, or `MD_not_available`. | Required. |
| `MD_mapping_rule_version` | Stable label `cluster_v2_posthoc_MD_mapping_20260711`. | Required. |
| `source_persistent_contacts` | Provenance label for the released contact region. | `not_available` for unavailable MD. |

### Evidence labels, redocking, and unresolved ledger

| Column(s) | Definition | Missing behavior or values |
|---|---|---|
| `automated_evidence_label`, `static_top_evidence_label`, `static_top_label` | `R_auto`, `A_auto`, `U_auto`, or `X_auto`; prefixed fields copy the applicable candidate label. | Required for ranked candidates. |
| `label_reason` | Machine-generated reason for the automated label. | Required. |
| `independent_evidence_details` | Semicolon-delimited assembly-interface and/or MD-run support details. | Empty when no independent support is present. |
| `unresolved_flag` | True for unresolved evidence, ligand selection, or specified QC states. | Required. |
| `label_scope` | `static_top3` or `non_top3_posthoc_mapping`. | Required. |
| `interpretation_guardrail` | Fixed reminder that labels are reproducible evidence states, not manual/experimental confirmation. | Required. |
| `reference_cluster_v2_id`, `reference_cluster_static_rank` | First `R_auto` candidate and its immutable static rank for a target. | Empty when no `R_auto` candidate exists. |
| `Raw_ligand_RMSD_A`, `redocking_RMSD_A` | Nonnegative ligand RMSD in Å; case-table field is a copy. | Required in the released redocking summary. |
| `GlideScore_kcal_per_mol` | Released docking score in kcal/mol. | Empty when unavailable. |
| `RMSD_threshold_call`, `redocking_status` | `RMSD <= 2 A`, `2 A < RMSD <= 3 A`, or `RMSD > 3 A`; case field is a copy. | Recomputed from numeric RMSD. |
| `Reference_pose_recovered` | True only for `Raw_ligand_RMSD_A <= 2.0`. | Required. |
| `Failure_or_warning_reason` | Released warning/failure text from the minimal redocking summary. | May be empty if no warning was recorded. |
| `redocking_role` | Fixed statement that redocking is post-hoc chemical validation, not static ranking input. | Required. |
| `issue_id` | Stable sequential unresolved-case identifier such as `U001`. | Required and unique. |
| `level` | Scope of the issue, for example candidate or target. | Required. |
| `issue_type`, `status` | Machine-readable issue class and current state. | Required. |
| `required_user_or_future_input`, `current_handling` | What could resolve the issue and how the current analysis treats it. | Required narrative fields. |

## Endpoint, uncertainty, and ablation columns

| Column(s) | Definition | Missing behavior or range |
|---|---|---|
| `Pocket_category`, `Protein_family` | Reviewed target category and family from `data/metadata/systems.tsv`. | Required; three category mappings are explicitly provisional in the input metadata. |
| `formal_cluster_v2_count`, `target_N` | Candidate count for a target or target count in a category. | Nonnegative integer. |
| `reference_evaluable`, `reference_evaluable_N` | Whether a target has a resolved nonempty ligand and the count of such targets. | Required; denominators use only evaluable targets. |
| `R_auto_cluster_count` | Number of `R_auto` candidates for a target. | Nonnegative integer. |
| `reference_first_rank`, `reference_rank`, `Full_reference_rank` | Minimum static rank among `R_auto` candidates. | Empty when no `R_auto` candidate exists. |
| `reference_rank_percentile` | `1-(rank-1)/(candidate_count-1)`, with 1 for a one-candidate target. | Empty without a finite reference rank; range `[0,1]`. |
| `reference_top1`, `reference_top3`, `reference_top5`; `Full_reference_Top1`, `Full_reference_Top3` | Inclusive rank-threshold indicators or their group means. | Target indicators are false for an evaluable target without `R_auto`; group means lie in `[0,1]`. |
| `reference_reciprocal_rank`, `Reference_MRR`, `Full_reference_MRR` | Target reciprocal reference rank and its evaluable-target mean. | Target value is 0 for evaluable target without `R_auto`; summaries are `[0,1]` or empty for no denominator. |
| `first_supported_evaluable`, `first_supported_evaluable_N` | Whether/count of targets having at least one `R_auto` or `A_auto` candidate. | Required. |
| `first_supported_rank`, `Full_first_supported_rank` | Minimum static rank among `R_auto` or `A_auto`. | Empty when none exists. |
| `first_supported_top1`, `first_supported_top3`; `First_supported_Top1`, `First_supported_Top3` | Inclusive target indicators and their first-supported-evaluable means. | Ratios lie in `[0,1]`; target indicators are false when not evaluable. |
| `first_supported_reciprocal_rank`, `First_supported_MRR` | Target reciprocal first-supported rank and its mean. | Empty for a target without support; summary empty for no denominator. |
| `Reference_Top1`, `Reference_Top3`, `Reference_Top5` | Overall, category, family-sensitivity, or bootstrap reference Top-k proportions. | `[0,1]` or empty for no evaluable target. |
| `Reference_median_rank`, `Reference_median_rank_percentile` | Median finite reference rank and percentile among reference-evaluable rows with a mapped rank. | Empty when no finite values exist. |
| `analysis_level` | `overall`, `pocket_archetype`, or `target`, according to table. | Required. |
| `group` | `all_targets` or a pocket-category label. | Required. |
| `bootstrap_method` | `target_resampling` or `family_clustered`. | Required. |
| `metric` | One of the seven bootstrapped Top-k/MRR metric names. | Required. |
| `point_estimate` | Metric on the unresampled 21-target table. | `[0,1]` for the released metrics. |
| `CI_2.5_percent`, `CI_97.5_percent` | Nan-aware 2.5th and 97.5th bootstrap quantiles. | `[0,1]`; empty only if every replicate is unavailable. |
| `iterations`, `random_seed` | Bootstrap iteration count and RNG seed; frozen as 10,000 and `20260710`. | Required integers. |
| `excluded_family`, `excluded_target_N`, `remaining_target_N` | Family removed in sensitivity analysis and removed/retained target counts. | Required; counts are nonnegative integers. |
| `Without_O_rel_score`, `Without_O_rel_rank` | Missing-aware score and descending average rank after removing only `O_rel_formal`. | Score/rank empty only if remaining modules cannot score. |
| `Without_O_rel_reference_rank`, `without_O_rel_reference_rank` | Minimum ablated rank among `R_auto` candidates; lowercase-prefixed field is the representative-case copy. | Empty when no `R_auto` candidate exists. |
| `Without_O_rel_first_supported_rank` | Minimum ablated rank among `R_auto` or `A_auto`. | Empty when none exists. |
| `reference_rank_change_without_minus_full`, `first_supported_rank_change_without_minus_full` | Ablated rank minus full rank. | Empty when either endpoint rank is unavailable; positive means worse. |
| `Without_O_rel_reference_Top1`, `Without_O_rel_reference_Top3` | Ablated reference Top-k indicator or category proportion. | Boolean at target level; `[0,1]` at category level. |
| `reference_Top3_inclusion_change` | Ablated Top-3 indicator minus full Top-3 indicator. | One of `-1,0,1`. |
| `rank_Spearman_rho`, `mean_rank_Spearman_rho` | Per-target Spearman correlation of full/ablated ranks and category mean. | `[-1,1]`; empty with fewer than two distinct ranks. |
| `Full_top_cluster_v2_id`, `Without_O_rel_top_cluster_v2_id`, `without_O_rel_top_cluster_v2_id` | Full and ablated rank-1 identifiers; lowercase field is the case-table copy. | Required when ranks are available. |
| `top_cluster_identity_changed`, `top_cluster_identity_change_fraction` | Per-target identity-change boolean and category mean. | Boolean and `[0,1]`, respectively. |
| `static_top_score`, `static_top_tool_support`, `static_top_chains`, `static_top_interface_overlap` | Static rank-1 score, tool count, contributing-chain list, and interface fraction for a representative case. | Interface value may be empty; other fields are required for configured cases. |
| `MD_calls` | Sorted semicolon list of non-apo MD concordance calls for the case. | `MD_not_available` when no applicable call exists. |

## Weight-sensitivity and native single-tool columns

| Column(s) | Definition | Missing behavior or range |
|---|---|---|
| `scenario_id`, `perturbed_module`, `direction` | Stable scenario identifier, one of the five static modules, and `decrease` or `increase`. | Required; the five-by-two Cartesian set is complete. |
| `multiplier` | `0.8` or `1.2`; the altered weight remains subject to candidate-level missing-aware renormalization. | Positive finite value. |
| `baseline_top1_retained_N`, `top1_retention_count`, `target_N` | Scenario- or target-level Top-1 identity stability counts and target denominator. | Nonnegative integers; target retention is out of ten scenarios. |
| `mean_top3_jaccard` | Mean across targets of baseline/scenario rank-at-most-3 set Jaccard. | `[0,1]`. |
| `median_spearman_rho` | Median target-level baseline/scenario candidate-rank Spearman correlation. | `[-1,1]`; empty only if correlations cannot be computed. |
| `top1_changed_scenarios` | Semicolon-separated scenario identifiers. | Empty when all scenarios retain the baseline Top-1 candidate. |
| `baseline_reference_rank`, `minimum_reference_rank`, `maximum_reference_rank`; first-supported counterparts | Baseline and across-scenario rank range for the fixed evidence-labeled candidate sets. | Empty when the corresponding fixed set has no candidate. |
| `method` | One of the five detector names or `OIPS-P`. | Required. |
| `complete_case`, `N` | Common-target membership and the common method denominator. | Required; frozen main analysis has 14 targets. |
| `mappable_region_count` | Number of mapped, cluster-deduplicated same-tool units, or cluster-v2 candidates for OIPS-P. | Positive integer on complete cases. |
| `output_status` | `unavailable`, `unmappable`, `no_hit`, or `reference_hit`. | Required. |
| `first_reference_associated_rank` | Earliest native regional rank whose mapped cluster is `R_auto`. | Empty for no hit; its Top-k contributions are false and reciprocal-rank contribution is zero. |
| `native_top_cluster_v2_id` | Cluster-v2 region selected first by the native tool order or OIPS-P rank. | Required on complete-case rows. |
| `Top1_N`, `Top3_N`, `Top5_N`; `Top1`, `Top3`, `Top5`, `MRR` | Complete-case hit counts and proportions/mean reciprocal rank. | Counts range from 0 to `N`; proportions lie in `[0,1]`. |

## Figure source columns

| Column | Definition | Missing behavior |
|---|---|---|
| `panel` | Figure-contract panel identifier. | Required. |
| `record_id` | Stable source-record identifier within the panel. | Required. |
| `series` | Plotted series name. | Required. |
| `group` | Category, method, case, or state used for grouping/coloring. | May be empty when a panel has no grouping dimension. |
| `x`, `y` | Numeric plotted coordinates or values. | Either may be empty when the corresponding axis value is categorical or not used. |
| `x_label`, `y_label` | Human-readable axis/category labels carried with the source row. | May be empty when not used by that mark. |
| `lower`, `upper` | Lower and upper interval endpoints. | Empty for marks without intervals. |
| `annotation` | Text annotation rendered or available for the mark. | Empty when no annotation is required. |
