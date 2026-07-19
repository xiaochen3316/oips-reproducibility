# Methods

## Scope and analysis boundary

This repository accompanies the manuscript *Choosing among Competing Pockets in Oligomeric Proteins: An OIPS-Assisted, Traceable Multi-Evidence Analysis*. It reproduces the released processed analysis for 21 oligomeric protein targets. The frozen input contains 1,742 pocket-tool records from DoGSiteScorer, DoGSite3, CavityPlus, CASTpFold, and SiteMap. It does not rerun those services, commercial pocket calculations, docking, or molecular-dynamics (MD) simulations from first principles.

The workflow has two physically and logically separate evidence layers:

| Layer | Permitted inputs | Products | Prohibited influence |
|---|---|---|---|
| Static candidate construction and prioritization | The 14-column tool feature table, prepared structures, and declared configuration | cluster-v2 candidates, static module scores, `OIPS-P_static`, and within-target ranks | Reference ligands, MD, redocking, literature labels, and post-hoc QC cannot alter clustering, scores, or ties. |
| Post-hoc evaluation and interpretation | The immutable static ranking plus reference, MD, redocking, system-category, and structure inputs | QC, evidence mappings and labels, target metrics, uncertainty, sensitivity, ablation, and figures | Post-hoc fields are written only to analysis outputs; the two static CSVs remain unchanged. |

`Q_evidence` is a static input-completeness module. Its name does not denote reference, MD, redocking, or literature evidence. `D_dyn_run_score` is carried only in the post-hoc MD mapping and never enters `OIPS-P_static`.

## Curated inputs and structure handling

The public feature table retains exactly `row_id`, `pdb_id`, `tool`, `pocket_id`, `display_order`, `sitemap_rank`, `center_x`, `center_y`, `center_z`, `center_method`, `residue_count`, `residue_set_json`, `pocket_geometry_score`, and `pocket_ligandability_score`. Empty numeric cells mean unavailable; they are not converted to zero. Records lacking both a complete center and a nonempty residue set are retained in the exclusion audit rather than silently discarded.

Prepared PDB or mmCIF coordinate files are parsed from `ATOM` and `HETATM` records. Static oligomer interfaces are constructed from non-hydrogen `ATOM` records belonging to the 20 standard amino acids. A protein residue is an interface residue when one of its heavy atoms is within 5.0 Å, inclusively, of a heavy atom from another chain. The resulting interface profile records interface residues and coordinates, chain count, protein-residue count, and inter-chain atom-contact counts. This is a geometric proxy; no symmetry-expanded crystal mates or interface energetics are generated.

Residue comparisons use chain, residue number, and insertion code after dropping the residue-name component. Thus a released identifier such as `A:GLY:10:` is compared as `A:10:`. This reconciles otherwise equivalent residue labels without changing chain or sequence position.

## Same-tool consolidation

Only records from the same target and same tool are considered for same-tool consolidation. Pairwise duplicate compatibility uses residue intersection-over-union (IoU), containment, center distance, and the configured DoGSite hierarchy rules. All positive similarity and maximum-distance thresholds are inclusive. Missing residue sets make IoU and containment unavailable; missing centers make center distance unavailable. The exact predicates and thresholds are specified in [oips-formula.md](oips-formula.md).

Groups are built deterministically in native-rank and `row_id` order. A record can join a group only if it is compatible with every existing member and the maximum diameter among available centers remains at most 8.0 Å. The representative minimizes the declared average distance/residue dissimilarity cost, with penalties of 0.20 for a missing center and 0.20 for a missing residue set; native rank, hierarchy depth, and `row_id` provide deterministic tie breaks. This stage produced 1,417 same-tool units from 1,564 mappable records; 178 source records were unmappable and remained in the audit.

## Cross-tool cluster-v2 reconstruction

Same-tool units are linked across different tools using one of three inclusive conditions: center distance at most 6.0 Å; center distance at most 10.0 Å together with residue IoU at least 0.20; or residue IoU at least 0.35. Candidate addition and merging also require a medoid-constrained group: the maximum distance between available centers must be at most 12.0 Å, and every non-medoid unit must be compatible with the medoid. Deterministic similarity and identifier rules resolve competing assignments.

A cluster may retain secondary same-tool units for traceability, but exactly one formal representative per tool contributes a vote and static module values. The formal representative minimizes cross-tool dissimilarity, then missingness, residue-set size, native rank, and source `row_id`. The residue envelope is the union of formal-representative residues. For a one-representative cluster, all envelope residues form the core; otherwise a residue enters the core when supported by at least `max(2, ceil(n/2))` of the `n` formal representatives.

The frozen reconstruction contains 733 cluster-v2 candidates. All satisfy the 12.0 Å diameter cap; the maximum is 11.914302329553335 Å. Each cluster-tool pair contributes at most one formal vote. Boundary-sensitive status is an audit flag, not an exclusion: for clusters supported by at least two tools, it is raised when diameter is greater than 9.0 Å, center dispersion is greater than 4.0 Å, core/envelope ratio is below 0.20, or median pairwise residue IoU is below 0.15. These strict audit inequalities differ intentionally from the inclusive 12.0 Å validity cap. Ninety-six frozen candidates are boundary-sensitive.

## Static scoring and ranking

Five modules are calculated from formal representatives and the prepared static structure:

- `C_cons`: tool-weighted formal support, with a 0.5 factor for a single-tool cluster.
- `G_geo`: mean of the highest three available geometry scores, or all available scores when fewer than three exist.
- `P_lig`: mean of the highest two available ligandability scores, or all available scores when fewer than two exist.
- `Q_evidence`: a bounded score derived from representative count, fractions with centers and residue sets, and presence of SiteMap.
- `O_rel_formal`: oligomer-interface relevance derived from pocket/interface residue overlap, chain context, center-to-interface distance, and interface extent.

The configured weights are 0.22, 0.18, 0.24, 0.12, and 0.24 for `C_cons`, `G_geo`, `P_lig`, `Q_evidence`, and `O_rel_formal`, respectively. `OIPS-P_static` is the weighted mean over available modules only. Missing modules remain missing and their weights are removed from the denominator. If every module were unavailable, the score would be unavailable; a valid released cluster always retains enough static information to score.

Candidates are ranked within each target by descending recomputed `OIPS-P_static`. Exact score ties receive average ranks; neither post-hoc evidence nor identifiers break a scientific tie. `cluster_v2_id` is used only for stable serialization among equal rank values. The frozen 733-row ranking has no ties, but fractional ranks remain part of the contract.

### One-at-a-time module-weight sensitivity

Each of the five configured module weights is independently multiplied by 0.8 and 1.2, yielding ten perturbation scenarios. In every scenario, missing modules are still removed from both numerator and denominator, so the available weights are renormalized candidate by candidate. Candidate membership, module values, missingness, cluster identifiers, and all post-static evidence labels remain fixed; only the perturbed static score and its within-target rank are recomputed.

Each scenario reports baseline Top-1 identity retention, mean target-level Top-3 Jaccard overlap, median within-target Spearman rank correlation, first-reference Top-1 and Top-3 counts, and first-supported Top-1 and Top-3 counts. The target-level table also records the minimum and maximum reference and first-supported ranks over all ten scenarios. In the frozen analysis, 20 of 21 targets retain the same Top-1 candidate in all ten scenarios. The only change is for 7O2I when `P_lig` is decreased by 20%, where V2C003 replaces V2C001.

## Frozen static boundary and post-hoc mapping

The serialized static ranking is copied before post-hoc analysis and tested for equality afterward. Post-hoc analysis adds no columns to `results/reference/static/cluster_v2_master_table.csv` or `results/reference/static/cluster_v2_static_rankings.csv`. A convenience master containing copied static columns plus selected evidence fields is written separately as `results/reference/analysis/final_cluster_v2_master_table.csv`.

### Automated top-3 QC

The top three ranked candidates per target undergo deterministic structural QC for possible over-merging, clear exclusion, possible splitting, low-information surface noise, proximity to missing sequence positions, and boundary sensitivity. PDB `REMARK 465` records are preferred for missing-residue assessment; otherwise short internal numbering gaps of 2–25 positions are used as a proxy. Crystal-contact risk remains explicitly not evaluable because the workflow does not generate crystallographic symmetry mates. QC statuses preserve candidates and expose uncertainty; they do not rerank the static list.

### Reference-ligand mapping

For each target, the curated ligand key is selected from the prepared structure. The reviewed 5LGE key is `6VN:D:503:`. Reference-contact residues are standard-amino-acid atoms within 4.5 Å, inclusively, of any selected ligand atom. Every static candidate is compared with the ligand centroid and contact set. The association rule combines center-to-ligand distance (DCC), contact recall, and residue IoU using inclusive thresholds. A target is reference-evaluable when a nonempty selected ligand key resolves to at least one ligand atom, independently of whether any candidate passes the association rule.

### MD mapping

Released MD evidence consists of precomputed persistent contact regions, optional region centers, and `D_dyn_run_score`; full trajectories are not read. Each MD run is compared with every immutable candidate. The best-mapped candidate is selected by descending IoU, descending overlap count, ascending center distance, and then identifier. The static top candidate is classified as concordant, partially concordant, boundary-shifted, conflicting, apo-only context, insufficient evidence, or MD not available by the declared inclusive rules. There are 23 released MD runs across 15 targets; six targets receive explicit `MD_not_available` rows.

### Evidence labels and redocking

`R_auto` denotes a candidate that passes the reference rule without a clear QC exclusion or unresolved ligand selection. `A_auto` denotes independent assembly-interface or MD-region support when `R_auto` is absent. `U_auto` denotes an acceptable candidate without sufficient automated reference or alternative support, and `X_auto` is reserved for clear automated exclusion or invalid mappability. These are reproducible evidence states, not manual or experimental confirmation.

Redocking is attached to the first `R_auto` cluster per target after static ranking. Pose recovery is recomputed from numeric RMSD: at most 2.0 Å is recovered; greater than 2.0 and at most 3.0 Å is intermediate; greater than 3.0 Å is not recovered. Redocking does not change static scores or ranks.

## Endpoints, uncertainty, and sensitivity

For a reference-evaluable target, the reference rank is the smallest static rank among `R_auto` candidates. A target with a resolved reference but no `R_auto` candidate contributes failure to Top-k and zero to reciprocal rank. The first-supported rank is the smallest rank among `R_auto` or `A_auto` candidates and is summarized only over targets for which such a candidate exists. Top-1, Top-3, Top-5, mean reciprocal rank (MRR), median reference rank, and median rank percentile are calculated as defined in [oips-formula.md](oips-formula.md).

Uncertainty uses 10,000 bootstrap iterations and seed `20260710`. Target resampling draws 21 target rows with replacement. Family-clustered resampling draws the 13 protein-family labels with replacement and includes every target belonging to each sampled family. Each method initializes its own random-number generator with the same seed. Percentile limits are the 2.5th and 97.5th nan-aware quantiles.

Family sensitivity removes one protein family at a time and recomputes the endpoint table. Oligomer-relevance ablation removes only `O_rel_formal`, renormalizes the remaining available module weights, reranks by descending ablated score with average ties, and compares target ranks, Top-k inclusion, top identity, and Spearman rank correlation. It is a sensitivity analysis, not a retrained model or causal effect estimate.

### Native single-tool comparison

The comparator starts from the same released raw pocket records. Records are first consolidated into the same-tool units used by the main reconstruction. Units are ordered by retained native tool output order (`display_order`, then `sitemap_rank`, then `row_id`). A unit is mapped through the exhaustive source-to-cluster-v2 crosswalk; when more than one unit from a tool maps to the same cluster, only the earliest native occurrence is retained. A single-tool region is reference-associated when its mapped cluster passes the same final `R_auto` definition used for OIPS-P.

The complete-case comparison is restricted to the 14 targets having at least one mappable region from every one of the five tools. For each tool and OIPS-P, the analysis reports the rank of the first reference-associated region and aggregates Top-1, Top-3, Top-5, and MRR over the same 14 targets. A target without a reference-associated region contributes failure to Top-k and zero to MRR. This preserves native tool ordering and a common evaluation rule; it does not recalibrate or retrain any method.

## Determinism and frozen results

CSV outputs use UTF-8, LF line endings, `%.15g` float formatting, empty fields for missing numeric values, lowercase booleans, declared column order, and stable sort keys. Input, configuration, output, and environment hashes are recorded in manifests and reports. Validation recomputes static scores, checks cross-table identities, checks the diameter cap, rebuilds figure source data, and compares frozen numeric snapshots using named tolerances.

The frozen overall reference results are Top-1 12/21, Top-3 18/21, Top-5 19/21, and MRR 0.723015873015873. First-supported results are Top-1 20/21, Top-3 21/21, and MRR 0.9761904761904762. These are retrospective results for this released 21-target set and should be interpreted with the limitations in [limitations.md](limitations.md).
