# OIPS Reproducibility Package

[简体中文](README.zh-CN.md) | [Reproduction guide](REPRODUCE.md)

This is the **public code and processed-data reproducibility repository** for
the OIPS study. It contains the processed data, deterministic code, equations,
and frozen reference outputs accompanying the manuscript:

> *Choosing among Competing Pockets in Oligomeric Proteins: An
> OIPS-Assisted, Traceable Multi-Evidence Analysis*

The package reconstructs cluster-v2 pocket candidates, calculates the static
OIPS-P score, maps post-hoc evidence, computes evaluation endpoints, and
regenerates the released figures and reports. Its purpose is to make the
reported processed analysis inspectable and repeatable.

## Quick start

Use Python 3.11 or 3.12 in a clean environment from the repository root:

```console
python -m pip install --editable ".[test]"
python -m oips_repro reproduce --config config/manuscript.yaml --output results/reproduced
python -m oips_repro verify --config config/manuscript.yaml --bundle results/reproduced --snapshot tests/scientific/data/expected_summary.json
```

The second command is the one-command reproduction workflow. It validates the
released inputs before writing a new bundle. The third command checks that
bundle against the frozen scientific snapshot. See [REPRODUCE.md](REPRODUCE.md)
for platform-neutral setup, individual pipeline stages, expected files, and
safe overwrite behavior.

## Frozen headline results

The released snapshot contains 21 targets and 1,742 source-tool records. Of
these, 1,564 records are mappable, 178 are retained in the exclusion audit,
and same-tool consolidation yields 1,417 units. Cross-tool clustering produces
733 cluster-v2 candidates. The maximum frozen cluster diameter is
11.914302329553335 A, below the inclusive 12 A cap, and no cluster contains
more than one formal vote from the same tool.

For the 21-target reference-ligand endpoint, the frozen static ranking gives
Top-1 = 0.5714285714, Top-3 = 0.8571428571, Top-5 = 0.9047619048, and
MRR = 0.7230158730. Machine-readable expected values and tolerances live in
`tests/scientific/data/expected_summary.json`; the files under
`results/reference/` are the immutable canonical outputs.

## Static ranking and post-hoc evidence

OIPS-P_static is computed only from released static pocket features,
cross-tool consensus, geometry, ligandability, evidence quality, and oligomer
relevance. Reference ligands, molecular-dynamics summaries, and redocking
summaries do not enter candidate construction or the static score. They are
mapped after the static ranking has been frozen and are used only for audit,
interpretation, and evaluation.

The exact equations, inclusivity rules, missing-value behavior, and ranking
conventions are defined in [docs/oips-formula.md](docs/oips-formula.md). The
narrative workflow is described in [docs/methods.md](docs/methods.md).

## Repository map

- `src/oips_repro/`: deterministic implementation and command-line interface.
- `config/`: frozen manuscript configuration, schema, and figure contract.
- `data/`: curated public input snapshot, checksums, rights audit, and external
  archive inventory.
- `results/reference/`: immutable validated tables, figure source data, and
  reports used as the comparison target.
- `results/reproduced/`: user-generated output; never used as the reference.
- `figures/manuscript/`: repository summary figures. Their names deliberately
  avoid asserting a final manuscript figure number; the manuscript mapping
  must be confirmed in `FIGURE_SOURCE_DATA_INDEX.tsv` before publication.
- `docs/`: methods, equations, data dictionary, provenance, limitations, and
  journal-facing statements.
- `tests/`: unit, integration, security, and frozen scientific checks.

The analysis bundle also includes the ten one-at-a-time ±20% module-weight
perturbations and the native single-tool comparison on the 14 five-tool
complete cases. Their scenario-, target-, and method-level CSV files are under
`results/reference/analysis/` and are rebuilt by the standard workflow.

For article review, [`MANUSCRIPT_RESULT_TRACEABILITY.tsv`](MANUSCRIPT_RESULT_TRACEABILITY.tsv)
maps Figures 1–7 and Tables 1–3 to the exact released source tables and marks
the final structural artwork or trajectory-level panels that remain
author-managed.

## Scope and limitations

This repository reproduces the released processed analysis. It does **not**
rerun proprietary web services, commercial docking software, or full molecular
dynamics simulations from first principles. Raw service packages, licensed
binaries and logs, trajectories, restart files, and checkpoints are excluded.
Their status and redistribution boundary are recorded in
[THIRD_PARTY.md](THIRD_PARTY.md), `data/metadata/asset-rights.tsv`, and
`data/external_archive_manifest.tsv`.

The released structures and normalized service-derived measurements retain
their source terms and attribution requirements. Missing upstream versions,
access dates, and archive identifiers are not guessed; they remain explicit
pre-publication checks.

## Citation and release status

Software-author metadata is recorded in [CITATION.cff](CITATION.cff). The
manuscript, code archive, data archive, and public repository do not yet have
verified public identifiers in this pre-publication package, so no DOI or
repository URL is asserted here. Add them only after reservation and landing-
page verification, following [PRE_PUBLICATION_CHECKLIST.md](PRE_PUBLICATION_CHECKLIST.md).

The planned software version is `1.0.0`, with intended tag
`v1.0.0-manuscript`. Until the checklist is complete and that exact commit has
been archived, cite the software metadata in `CITATION.cff` together with the
specific commit supplied with the materials.

## Licensing

Code is licensed under the BSD 3-Clause License in [LICENSE](LICENSE).
Original team-authored documentation, tabular data, and figures are licensed
under CC BY 4.0 in [LICENSES/CC-BY-4.0.txt](LICENSES/CC-BY-4.0.txt), except
where a file or record identifies third-party terms. Third-party material is
not relicensed by this repository; see [THIRD_PARTY.md](THIRD_PARTY.md).
