# OIPS toy reproduction fixture

This directory is a tiny, self-contained, CPU-only repository for exercising
the public OIPS pipeline without network access. It uses three synthetic
prepared structures labelled with the manuscript representative IDs. The
coordinates and scores are deliberately minimal test data, not biological
claims or replacements for the public snapshot.

The feature table contains three candidate regions per target from multiple
tools, plus one record with neither a center nor residues so the exclusion
path is exercised. One target has a small MD summary; the other two exercise
the explicit `MD_not_available` state. All input paths are repository-relative
and all locally declared data payloads are covered by `data/SHA256SUMS`.

Run from the project repository with, for example:

```text
python -m oips_repro reproduce --config examples/toy/config/manuscript.yaml --output <temporary-directory>
```

Do not commit generated result bundles beneath this fixture.
