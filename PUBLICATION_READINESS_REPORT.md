# Publication readiness report

## Decision: READY_WITH_MANUAL_REVIEW

The released processed-analysis workflow, frozen outputs, tests, and security
scan pass locally. The repository is suitable for author review and eventual
GitHub upload after the manual items below are resolved.

### Completed technical checks

- The converted archive was compared file-by-file with the prior scientific
  payload: all 158 prior public files were byte-identical, and the 415 extra
  entries were development metadata, caches, worktrees, duplicate reproduced
  outputs, or tests rather than missing scientific inputs.
- The manuscript's ±20% weight sensitivity and 14-system single-tool
  comparison now have record-level code, frozen outputs, documentation, and
  tests. Their published values are reproduced exactly.

### Blocking issues for final tagged release

- Third-party transformed-structure and normalized-data redistribution must
  receive final author/rights sign-off.
- Repository URL, tag, final commit, and archive DOI values do not yet exist.

### Non-blocking issues

- Final manuscript molecular-rendering artwork remains author-managed; the
  released data sources and regeneration boundary are indexed in
  `MANUSCRIPT_RESULT_TRACEABILITY.tsv`.
- Exact upstream service/software versions and access dates remain unavailable
  where the reviewed source records did not establish them.

### Manual confirmations required

- author order, affiliations, and public corresponding-author emails;
- confirmation of final artwork/version mapping;
- third-party audit and PDB attribution;
- scientific lead, rights reviewer, and corresponding-author approval.

### Excluded files

Internal development plans/reports, local environments and caches, raw
archives, raw service/commercial-software bundles, trajectories, checkpoints,
logs, binaries, credentials, private job identifiers, and duplicate generated
outputs.

### Test summary

- Full processed reproduction, frozen snapshot verification, and the complete
  test suite pass; exact final counts are recorded in `TEST_REPORT.md`.

### Third-party uncertainties

The included prepared PDB-derived structures and normalized factual outputs
have documented provenance and provisional repository approval. Raw upstream
and commercial files remain excluded. Final record-level human confirmation is
still required before a public tag.
