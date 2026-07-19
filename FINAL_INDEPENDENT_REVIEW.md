# Final independent review

## PASS_WITH_ACTIONS

The repository was reread from its public candidate tree. README commands map
to real CLI entry points; input validation, full processed reproduction,
snapshot comparison, the complete test suite, and the tracked-file release scanner pass.
Static inputs and post-static evidence are physically separated, and the code
path constructs and ranks candidates before reference/MD/redocking mapping.
Post-static evidence does not modify frozen OIPS-P_static ranking.

Traceability was confirmed for cluster reconstruction, missing-aware scoring,
Top-3 QC, both evaluation endpoints, target and family bootstrap, LOFO, O_rel
ablation, redocking, MD concordance, representative cases, unresolved cases,
the ±20% weight sensitivity, native single-tool complete-case comparison,
and repository figure source data. Denominators are explicit: Top-3 automated
QC evaluates 63 candidates (including 12/63 boundary-sensitive), whereas the
all-candidate audit evaluates 733 candidates (96/733 boundary-sensitive).

Actions required before creating the public repository or final tag:

- obtain third-party and author-metadata confirmation;
- confirm the final manuscript artwork/version mapping.

The converted supplement comparison and the two formerly missing manuscript
analyses are complete and no longer appear in the action list.

Actions that may wait until a real GitHub/archive record exists:

- verified repository URL, release tag, final commit SHA, and DOI fields.

Never upload the original ZIP/7z, `.git` from this review copy, environments,
caches, internal plans, raw proprietary projects, trajectories, checkpoints,
credentials, or private identifiers.

Recommended GitHub repository name: `oips-reproducibility`.
