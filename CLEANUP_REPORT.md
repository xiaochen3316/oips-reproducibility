# Cleanup report

The public candidate was created from a separate extraction of the supplied
ZIP. The original ZIP, 7z, and audit extraction were not modified.

Excluded from the public tree:

- `.superpowers/sdd/task-7-report.md` — internal development report.
- `docs/superpowers/plans/` — internal implementation plan.
- `docs/superpowers/specs/` — internal design record.
- Local `.venv`, Python bytecode, pytest cache, and editable-install metadata
  created during validation — reproducible local build artifacts.

The scientific definitions present in the excluded internal documents were
already represented in `docs/methods.md`, `docs/oips-formula.md`,
`docs/provenance.md`, `docs/data-dictionary.md`, configuration, tests, and
frozen outputs; no scientific rule was migrated from an internal file.

Four auto-generated plots and their source tables were renamed with the
`repository_summary_figure_` prefix because final manuscript numbering has not
been author-confirmed. No numeric content was changed.
