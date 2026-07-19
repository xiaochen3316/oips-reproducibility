# Reproduction report

## Pipeline executed

1. `validate-data`
2. one-command `reproduce`
3. frozen-snapshot `verify`

The clean reproduction used the released `config/manuscript.yaml` and wrote to
an external ephemeral output directory so the frozen reference bundle remained
untouched.

## Results

- source-tool records: 1,742
- mappable records: 1,564
- excluded/unmappable records: 178
- same-tool consolidated units: 1,417
- cluster-v2 candidates: 733 across 21 targets
- boundary-sensitive candidates: 96/733
- rankings: 733; ties: 0
- maximum static-score recomputation difference: 0.0
- analysis tables: 20
- bootstrap rows: 14
- regenerated repository summary figures: 4 (12 figure-stage files including
  SVG, PNG, and source CSV)
- weight perturbation scenarios: 10
- native single-tool complete-case targets: 14
- reproduced bundle files verified: 44
- snapshot checks: 24 passed, 0 failed

The reproduced outputs matched the frozen scientific snapshot within declared
tolerances, including the manuscript's ±20% weight perturbation and 14-system
single-tool comparison. The workflow does not rerun upstream pocket services,
commercial docking, or full MD trajectories.

Final determination: **processed-analysis reproduction passed**.
