# Journal-facing statements

## Status

These statements accompany the manuscript *Choosing among Competing Pockets in Oligomeric Proteins: An OIPS-Assisted, Traceable Multi-Evidence Analysis*. Author order, affiliations, and correspondence below are verified. The public repository URL, code DOI, data DOI, manuscript DOI, and author ORCIDs are not yet available. Bracketed release gates must be resolved before submission or publication; they are not placeholders that may be silently omitted.

## Verified manuscript metadata

### Title

Choosing among Competing Pockets in Oligomeric Proteins: An OIPS-Assisted, Traceable Multi-Evidence Analysis

### Authors

Xiao Chen<sup>1</sup>, Yifan Zhu<sup>1</sup>, Shilei Zhao<sup>1</sup>, Xin Zhang<sup>1</sup>, Haopeng Sun<sup>2,*</sup>, Yao Chen<sup>1,*</sup>, Xin Xue<sup>3,*</sup>

### Affiliations

1. Nanjing University of Chinese Medicine, Nanjing 210023, China.
2. School of Pharmacy, China Pharmaceutical University, Nanjing, 211198, China.
3. National and Local Collaborative Engineering Center of Chinese Medicinal Resources Industrialization and Formulae Innovative Medicine, Jiangsu Collaborative Innovation Center of Chinese Medicinal Resources Industrialization, Jiangsu Key Laboratory for High Technology Research of TCM Formulae, Nanjing University of Chinese Medicine, Nanjing 210023, China.

### Correspondence

Correspondence and requests for materials should be addressed to Haopeng Sun (`sunhaopeng@163.com`), Yao Chen (`300630@njucm.edu.cn`), or Xin Xue (`xuexin@njucm.edu.cn`).

## Data Availability Statement — pre-publication draft

The curated processed input data, prepared structure files approved for redistribution, frozen cluster-v2 and OIPS analysis tables, figure source data, and validation reports supporting this study are provided with the companion repository for *Choosing among Competing Pockets in Oligomeric Proteins: An OIPS-Assisted, Traceable Multi-Evidence Analysis* **[public repository URL: not yet assigned]**. The versioned data archive identifier is **[data DOI: not yet assigned]**.

The repository reproduces the released processed analysis. It does not reproduce upstream proprietary web-service or commercial pocket calculations, commercial redocking, or complete molecular-dynamics simulations from first principles. Full MD inputs and trajectories, restart/checkpoint files, representative solvent snapshots, and large raw pocket-output packages are not included in the Git payload. `data/external_archive_manifest.tsv` records their known availability, archive, completeness, and license states without exposing local paths, service sessions, job identifiers, or restricted formats. Assets marked unavailable or not archived must not be described as deposited.

RCSB PDB-derived prepared coordinates remain subject to their recorded source terms and attribution. Commercial-software binaries, logs, license details, and raw third-party service packages are not redistributed. The exact local scientific payload, source identifiers, byte counts, SHA-256 values, media types, access states, and license/terms statements are recorded in `data/manifest.tsv`, `data/SHA256SUMS`, and `data/metadata/asset-rights.tsv`.

## Code Availability Statement — pre-publication draft

The Python code, configuration, tests, deterministic toy example, and one-command validation/reproduction workflow are available in the companion repository **[public repository URL: not yet assigned]**. The archived software identifier is **[code DOI: not yet assigned]**. The exact release commit, tag, configuration hash, input hashes, environment hash, and key-result hashes will be recorded in `release/manifest.json` and must agree with the archived release before this statement is finalized.

The code license covers original repository code; it does not grant rights to third-party services, commercial software, or restricted source artifacts. Compatibility ranges are declared in `pyproject.toml`, and the release-tested environment is pinned separately. A standard reproduction writes to `results/reproduced/`; `results/reference/` is immutable and used only for comparison.

## Reproducibility Statement

Static candidate reconstruction and `OIPS-P_static` ranking use only the released 14-column pocket-feature table, prepared structures, and declared configuration. Reference-ligand, MD, redocking, and literature-derived fields are physically separated and are loaded only after the static ranking is frozen. Post-hoc analysis cannot alter cluster membership, static module values, within-target ranks, or ties.

The frozen workflow validates 1,742 source records across 21 targets, including 1,564 mappable records, 178 audited exclusions, 1,417 same-tool units, 733 formal cluster-v2 candidates, and 96 boundary-sensitive candidates. It checks one formal vote per tool, an inclusive 12.0 Å maximum diameter, missing-aware score recomputation, cross-table identities, bootstrap metadata, sensitivity results, figure-source reconstruction, manifests, checksums, and sensitive-content rules. The frozen overall reference endpoints are Top-1 12/21, Top-3 18/21, Top-5 19/21, and MRR 0.723015873015873; these are retrospective results for the released set.

## Citation statement — release gate

The accompanying manuscript title is the title shown above. The software title used by Citation File Format metadata is `OIPS Reproducibility Package`; these two titles are not interchangeable. A preferred paper citation must not be added until a final publication citation exists.

- Public repository URL: **not yet available**.
- Code DOI: **not yet available**.
- Data DOI: **not yet available**.
- Manuscript DOI: **not yet available**.
- Author ORCIDs: **not yet available**.

The final software citation must use the verified author order and spelling and the exact archived version. No generic team author, inferred ORCID, provisional URL, or invented DOI may be substituted.

## Third-party materials statement

Normalized factual measurements and minimal team-derived summaries are distributed only when the rights audit marks them `approved_for_repository`. RCSB PDB coordinates retain RCSB source attribution and terms. Raw CASTpFold, ProteinsPlus/DoGSite, CavityPlus, and Schrödinger-derived packages, as well as licensed software artifacts, are not relicensed by the authors and are not included unless their exact asset class is explicitly approved. The authoritative asset-by-asset boundary is `data/metadata/asset-rights.tsv`; narrative statements do not override it.

## Statements still requiring author confirmation

The following manuscript declarations were not supplied with the verified identity metadata and are therefore not drafted as facts:

- CRediT author-contribution roles;
- funding sources and grant numbers;
- competing-interest declaration;
- acknowledgements;
- any ethics statement, if applicable; and
- author ORCIDs.

These items must be provided and approved by the authors before journal submission. Their absence must not be converted into “none declared” or any other substantive statement without confirmation.
