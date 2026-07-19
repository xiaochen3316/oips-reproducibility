# Third-party materials and redistribution boundaries

This repository combines team-authored code and curation with selected
third-party-origin scientific inputs and normalized measurements. The BSD and
CC BY licenses in this repository apply only where the authors have authority
to grant them. They do not replace the terms of source databases, services, or
commercial software.

The record-level authority is `data/metadata/asset-rights.tsv`; the local asset
inventory is `data/manifest.tsv`; excluded or externally archived material is
listed by stable, non-sensitive labels in `data/external_archive_manifest.tsv`.

## Prepared structural coordinates

The paired structures under `data/structures/` are team-prepared derivatives
of Protein Data Bank entries obtained through RCSB PDB. The source PDB
identifiers and persistent identifiers are recorded in `data/manifest.tsv`.
RCSB states that PDB archive data are available under CC0 1.0; users should
still cite the original structure publications and the PDB entries.

- Source and policy: <https://www.rcsb.org/pages/usage-policy>
- Repository treatment: transformed coordinates are included; source identity,
  checksum, and transformation provenance are retained.
- License boundary: neither the root BSD license nor the team CC BY notice
  overrides the PDB source terms.

## Normalized pocket features

`data/static/tool_pocket_features.csv` contains normalized factual fields used
by the public analysis, including centers, residues, ranks, and selected scalar
measurements derived from CASTpFold, ProteinsPlus tools
(DoGSiteScorer/DoGSite3), CavityPlus, and SiteMap workflows.

- Raw downloads, service packages, logs, binary files, credentials, private
  URLs, and job identifiers are excluded.
- Provider terms and method citations remain applicable; inclusion of
  normalized factual measurements does not relicense a provider's software or
  raw output package.
- Exact service/software versions and access dates not established by the
  reviewed source records remain blocking items in
  `PRE_PUBLICATION_CHECKLIST.md`.

Relevant provider information retained by the rights audit:

- CASTpFold: <https://cfold.bme.uic.edu/castpfold/>
- ProteinsPlus: <https://proteins.plus/>
- DoGSite3 service help: <https://proteins.plus/help/dogsite3_rest>
- CavityPlus method record: <https://doi.org/10.1093/nar/gky380>
- Schrodinger terms: <https://www.schrodinger.com/eula/>

## Molecular-dynamics and redocking summaries

The CSV files under `data/posthoc/` contain only the reviewed, minimal derived
fields required for post-hoc mapping and evaluation. Full MD inputs,
trajectories, restart/checkpoint files, representative solvent-context
snapshots, commercial docking project files, software binaries, license
details, benchmark disclosures, and logs are not redistributed.

These summaries are not used to construct cluster-v2 candidates or calculate
OIPS-P_static. Their source and review status are recorded in the rights and
manual-decision tables. Any future external deposit must retain vendor and
structure terms and must be reviewed independently from this Git release.

## Python dependencies

Python dependencies are installed from their normal distribution channels and
are not vendored in this repository. Each dependency remains under its own
license. Compatibility ranges are declared in `pyproject.toml`; the exact
release-tested versions will be recorded in `environment/constraints.txt`.

## Team-authored content

Subject to the exceptions above and any record-level notice, the original
documentation, curated tabular data, and figures are offered under CC BY 4.0.
The source code is offered under BSD-3-Clause. When redistributing a mixed
artifact, preserve this file, the manifests, source attribution, and all
applicable notices.
