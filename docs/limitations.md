# Limitations

## Intended interpretation

OIPS is a traceable prioritization framework for comparing competing pocket candidates in this released oligomeric-protein set. The frozen results are retrospective descriptions of 21 targets, not evidence that OIPS is a universal pocket predictor, a prospective success-rate estimate, or a benchmark proving superiority over the five contributing pocket tools. Tool outputs are inputs to the consensus; the repository does not rerun or independently benchmark those tools.

`OIPS-P_static` is a ranking score, not a calibrated probability of ligandability, druggability, binding, or biological relevance. A score difference has no validated probabilistic interpretation. `R_auto`, `A_auto`, `U_auto`, and `X_auto` are reproducible evidence states, not experimental truth labels or manual expert conclusions.

## Dataset size and dependence

The analysis contains 21 targets, 733 formal candidates, 13 protein-family labels, and four pocket categories. Several families contain a single target, while HIV-1 protease and PD-L1 contribute repeated related structures. Consequently, target rows are not uniformly independent biological systems and category/family estimates can be dominated by a few observations.

Target-resampled and family-clustered bootstrap intervals describe sampling variation within this curated set. They do not account for uncertainty in tool generation, structure preparation, ligand annotation, category assignment, or choice of targets, and they do not establish performance on unseen families. Leave-one-family-out results are a sensitivity analysis, not external validation. Three pocket-category mappings remain explicitly provisional; category-stratified results are descriptive.

## Processed-input rather than first-principles reproduction

The repository reproduces clustering, static scoring, post-hoc mapping, statistics, and figures from committed processed inputs. It does not reproduce the following upstream computations from first principles:

- DoGSiteScorer, DoGSite3, CavityPlus, CASTpFold, or commercial SiteMap production runs;
- commercial redocking runs or licensed docking binaries;
- full MD setup, equilibration, production trajectories, or persistent-contact extraction; or
- preparation of every paired structure from an unmodified archive entry.

Raw service packages, job/session identifiers, result URLs, commercial logs, licensed binary formats, complete MD trajectories, inputs, checkpoints, and representative solvent snapshots are not Git payloads. Their availability and redistribution state are represented only in the external manifest when known. The current external inventory contains explicit not-available/not-archived states and no assigned external-archive DOI.

## Incomplete upstream metadata

The committed prepared structures are sufficient for this released processed analysis, and their file hashes are fixed. However, `data/metadata/systems.tsv` marks every target `author_confirmation_required` and records only approved category fields plus structure provenance. The discoverable source materials did not verify many upstream details, including structure version/date, biological assembly ID, chain mapping, retained/deleted hetero groups, missing-atom treatment, protonation tool and pH, force field, water model, MD engine/version, timestep, ensemble, random seed, trajectory processing, and exact pocket-tool versions.

These fields remain `not_available_in_source_materials`. They must not be interpreted as default settings, and the missing information prevents an independent reader from reconstructing every upstream preparation and simulation choice.

## Static candidate construction

Same-tool consolidation and cross-tool clustering use fixed geometric/residue thresholds and a deterministic greedy/merge order. The thresholds were frozen for traceability; they were not learned in a held-out optimization. Alternative residue canonicalization, centers, linkage rules, or thresholds could change candidate boundaries. A same-tool unit with a missing center can still join through residue containment, and missing centers are omitted from diameter calculations.

The 12.0 Å center-diameter cap constrains only available centers. It does not prove that an irregular residue envelope is a single physical pocket. Conversely, the medoid constraint and one-formal-vote-per-tool rule may split biologically connected subregions. Boundary, possible-split, and possible-overmerge flags expose some of this ambiguity but do not quantify all clustering uncertainty. Flagged candidates are retained rather than manually corrected.

Residue comparison drops the residue-name component and retains chain, residue number, and insertion code. This reconciles naming discrepancies but cannot repair an incorrect chain or sequence-number mapping. Residue sets and centers are detector-derived summaries; the public workflow does not return to raw grids, surfaces, or cavity meshes.

## Static scoring

`C_cons` depends on declared tool weights and treats the five tools as complementary formal votes, although their algorithms and outputs are not statistically independent. `G_geo` and `P_lig` average only the highest available representative values; missingness changes the contributing set. Although `OIPS-P_static` renormalizes missing modules rather than zero-filling them, candidates scored from different available-module subsets are not based on identical information.

The public schema requires finite geometry and ligandability values but does not enforce a universal 0–100 meaning across providers. `Q_evidence` measures static representation completeness and SiteMap presence; it is not independent experimental evidence. The module weights are fixed analysis choices, not uncertainty estimates.

Average within-target ranking preserves exact ties, but the frozen snapshot happens to contain no static ties. The absence of ties in this dataset does not guarantee tie-free behavior in another dataset.

## Oligomer-interface proxy

`O_rel_formal` uses a 5.0 Å inter-chain heavy-atom contact proxy on each committed prepared structure. It does not use PISA or an equivalent calculation of buried surface area, solvation free energy, hydrogen bonds, salt bridges, hotspots, or assembly stability. It does not expand crystallographic symmetry mates; therefore crystal-contact risk is explicitly not evaluable. Interface geometry can depend on the prepared chain set and assembly choice, which are not fully documented upstream.

The piecewise distance transform, interface-fraction transform, chain-context scores, and weights are frozen heuristic components. The `O_rel_formal` ablation changes only that module and renormalizes the other available modules. Differences between full and ablated ranks show sensitivity to this declared term; they do not identify a causal contribution or validate the term independently.

## Reference-ligand evaluation

Reference association is geometric and retrospective. It uses the curated ligand key, a ligand centroid, residues within 4.5 Å, and inclusive DCC/recall/IoU rules. The reviewed 5LGE correction is data-driven, but other ligand selection errors or biologically irrelevant crystallographic ligands could affect results. Exact-key selection is validated in the frozen snapshot; fallback modes remain in the code for transparent handling of other inputs.

An evaluable target without a passing `R_auto` candidate contributes failure to reference Top-k and zero to reference MRR. In contrast, first-supported metrics exclude targets without any `R_auto` or `A_auto` candidate. These denominators answer different questions and must not be compared as if they were identical endpoints. Multiple `R_auto` fragments for a target are retained; summary metrics use the smallest static rank.

No functional-residue annotation set was available for the released calculation. `functional_residue_overlap` is therefore `NA` with status `not_available_not_imputed` and contributes to no label or metric.

## MD and redocking evidence

The release contains 23 summarized MD runs for 15 targets and explicit unavailable rows for six targets. It maps previously extracted persistent contact sets and optional centers; it does not inspect trajectories. The underlying simulation contexts are heterogeneous, and an apo run with no persistent contacts is labeled `apo_only_context`, not negative evidence. `D_dyn_run_score` is carried through as a source summary and never affects static ranking.

MD concordance categories depend on fixed overlap and center-distance thresholds. A best-mapped alternative region can reflect a genuine dynamic site, a boundary shift, contact-definition choices, or upstream alignment/mapping error. The category alone cannot distinguish those mechanisms.

Redocking is represented by one minimal summary row per target. RMSD categories are recomputed from the released numeric RMSD, but raw poses, grids, protocols, logs, and licensed software artifacts are absent. Redocking is post-hoc reference-region validation and neither validates non-reference pockets nor changes `OIPS-P_static`.

## Figures, numerical precision, and release state

Figures are generated from validated analysis tables and four source-data CSVs. They summarize the released calculations; they do not provide additional independent evidence. Tables use deterministic `%.15g` serialization, while scientific comparisons use named absolute tolerances to avoid treating binary floating-point representation as exact decimal identity.

The manuscript author order, affiliations, and correspondence have been verified for documentation. ORCIDs, public repository URL, code DOI, data DOI, and manuscript DOI are not yet available. Until those identifiers and the remaining rights/metadata checks are completed, availability statements are publication drafts and the package must not be described as a finalized archival release.
