# Archive comparison report

## Sources

- Primary source: `OIPS-reproducibility-v1.0.0-prepublication-9c2d0b9(1).zip`.
- Supplemental source: `oips-reproducibility(1).7z`.

The original ZIP was readable and contained 203 archive entries (158 regular
files; approximately 12 MB after extraction). The initially supplied 7z could
not be read in this runtime, so the author supplied a direct ZIP conversion of
its unpacked contents. The converted archive contained the same 158-file public
source snapshot: every corresponding file had an identical SHA-256 digest.

The converted archive also contained 415 non-source extras in the inspected
release-check tree, dominated by embedded `.git` history and duplicated
`results/reproduced/` outputs. Other archive branches contained worktrees,
tamper-test outputs, local environments, caches, and repeated release-check
copies. These were excluded. No scientifically necessary file was present only
in the converted archive, and no file was copied from its Git history or
generated-output trees.

Status: **comparison complete; no missing public scientific payload found**.
