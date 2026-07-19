# Reproducing the released analysis

## Requirements

- Python 3.11 or 3.12.
- A clean checkout of this repository.
- No proprietary software or network service is required for the released
  processed-analysis workflow.

Python 3.11 is the manuscript baseline. The dependency ranges are declared in
`pyproject.toml`; the exact release-tested environment will be frozen in
`environment/constraints.txt` before public release.

## 1. Create an environment

Using `venv`:

```console
python -m venv .venv
```

Activate it with the command appropriate for your shell, then install the
package and test dependencies:

```console
python -m pip install --upgrade pip
python -m pip install --editable ".[test]"
```

Alternatively, use the reviewed Conda specification:

```console
conda env create --file environment/environment.yml
conda activate oips-repro
python -m pip install --editable .
```

## 2. Validate the released inputs

```console
python -m oips_repro validate-data --config config/manuscript.yaml
```

This checks the configuration, manifests, checksums, schemas, table relations,
rights metadata, frozen counts, and scientific invariants before computation.
A nonzero exit means the workflow must not proceed.

## 3. Run the complete workflow

```console
python -m oips_repro reproduce --config config/manuscript.yaml --output results/reproduced
```

The command runs clustering, static scoring, post-hoc analysis, figure-source
generation, SVG rendering, and report generation in the configured order. It
also records the input and output hashes needed for verification.

The workflow refuses to overwrite an existing managed bundle unless `--force`
is supplied. `--force` is intentionally restricted to recognized generated
files inside the selected output directory; it does not authorize deletion of
unknown files or modification of `results/reference/`.

## 4. Verify against the frozen snapshot

```console
python -m oips_repro verify --config config/manuscript.yaml --bundle results/reproduced --snapshot tests/scientific/data/expected_summary.json
```

Successful verification confirms the frozen row counts, distributions,
numeric metrics, representative cases, formulas, ranks, figure source values,
manifest context, and file hashes. Review the generated Markdown reports under
`results/reproduced/reports/` for human-readable evidence.

## Optional: run individual stages

The one-command workflow is preferred. For debugging, the same public CLI
exposes these stages:

```console
python -m oips_repro cluster --config config/manuscript.yaml --output results/staged/clustering
python -m oips_repro score --config config/manuscript.yaml --cluster-dir results/staged/clustering --output results/staged/static
python -m oips_repro analyze --config config/manuscript.yaml --cluster-dir results/staged/clustering --static-dir results/staged/static --posthoc-data data/posthoc --output results/staged/analysis
python -m oips_repro figures --config config/manuscript.yaml --analysis results/staged/analysis --output results/staged/figures
```

These commands are deterministic for the same validated inputs, configuration,
dependency environment, and random seed (`20260710`).

## Run the test suite

```console
python -m pytest -q
```

Before a public release, also run the release scanner with the completed
release manifest:

```console
python scripts/release_check.py --repository-root . --data-manifest data/manifest.tsv --release-manifest release/manifest.json
```

## Expected output layout

The complete bundle contains:

- `clustering/`: candidate, membership, mapping, boundary-audit, and exclusion
  tables.
- `static/`: static master table and rankings.
- `analysis/`: post-hoc mappings, endpoints, uncertainty, sensitivity,
  ablation, unresolved cases, and representative cases.
- Within `analysis/`, `weight_sensitivity_*.csv` contains the ten one-at-a-time
  ±20% module-weight perturbations and `single_tool_*.csv` contains the native
  detector comparison and 14-system complete-case metrics.
- `figure_source_data/`: one source CSV for each repository summary figure;
  final manuscript numbering is intentionally left for author confirmation.
- `figures/`: regenerated vector figures.
- `reports/`: rebuild, analysis, numeric cross-check, and validation reports.
- a machine-readable run manifest recording provenance and hashes.

The canonical comparison bundle is `results/reference/`. Never edit it to make
a failed reproduction pass; investigate the input, environment, or code change
instead.

## What this does not reproduce

This workflow starts from reviewed, normalized pocket features and selected
derived evidence. It does not reproduce proprietary service calls, commercial
docking, or full molecular-dynamics simulations from raw coordinates. See
`THIRD_PARTY.md`, `docs/provenance.md`, and `docs/limitations.md` for the exact
boundary.
