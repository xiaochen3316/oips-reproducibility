# Contributing

Thank you for helping improve the OIPS Reproducibility Package. Contributions
should preserve the traceability of the released scientific analysis.

## Before proposing a change

1. Open a focused issue or describe the scientific or technical reason in the
   pull request.
2. Work from a current branch and keep unrelated changes separate.
3. Do not add proprietary binaries, raw service packages, trajectories,
   checkpoints, credentials, private URLs, job identifiers, personal data, or
   absolute local paths.
4. Check `THIRD_PARTY.md` and `data/metadata/asset-rights.tsv` before adding any
   externally sourced material.

## Development setup

```console
python -m pip install --editable ".[test]"
python -m pytest -q
```

Use Python 3.11 or 3.12. Keep scientific functions deterministic and separate
pure calculation from filesystem and command-line concerns.

## Scientific changes

Any change to a formula, threshold, missing-value rule, ranking convention,
manual decision, input asset, or frozen result must include:

- a clear scientific rationale;
- a targeted test that fails before the implementation change;
- updated configuration or decision-audit records;
- regenerated outputs in a non-reference directory;
- an explicit review of whether the frozen reference bundle should change;
- updated methods, equations, data dictionary, limitations, and changelog when
  applicable.

Never modify `results/reference/` merely to make a failing test pass. Reference
changes require documented scientific approval and a new release decision.

## Documentation and data

- Use UTF-8, LF line endings, repository-relative paths, and stable sorting.
- Define every equation symbol and state ranges, threshold inclusivity, and
  missing-value behavior.
- Keep unavailable metadata explicit rather than guessing values.
- Update checksums and manifests through the reviewed preparation workflow.
- Preserve third-party attribution and redistribution boundaries.

## Pull-request checks

Run the complete test suite and, when release files are affected, the release
scanner. In the pull request, report the commands and their exact outcomes.
All contributors must follow `CODE_OF_CONDUCT.md`.
