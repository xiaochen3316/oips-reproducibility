# Pre-publication release checklist

**Status: `pre_publication_incomplete`**

This is the authoritative release gate. Do not create the public GitHub
release, the `v1.0.0-manuscript` tag, or a final archive until every blocking
item below is checked and independently verified.

## Confirmed manuscript metadata

- [x] Author order and spelling: Xiao Chen; Yifan Zhu; Shilei Zhao; Xin Zhang;
  Haopeng Sun; Yao Chen; Xin Xue.
- [x] Affiliations mapped to each author and recorded in `CITATION.cff` and
  `docs/journal-statements.md`.
- [x] Corresponding authors and email addresses confirmed: Haopeng Sun
  (`sunhaopeng@163.com`), Yao Chen (`300630@njucm.edu.cn`), and Xin Xue
  (`xuexin@njucm.edu.cn`).
- [x] Manuscript title confirmed in manuscript-facing documentation.

## Blocking identifier checks

- [ ] Confirm each author's ORCID, or record explicitly that the author has no
  ORCID or does not wish to include it. Do not infer identifiers from names.
- [ ] Create and verify the public GitHub organization/repository URL, then add
  the exact URL to README, CFF, and the release manifest.
- [ ] Reserve the code archive DOI and verify its draft landing page.
- [ ] Reserve the data archive DOI and verify its draft landing page.
- [ ] Add the manuscript DOI only after it has been assigned and verified.
- [ ] Ensure unavailable identifiers remain null or omitted; never use dummy
  DOI or URL values.

## Blocking third-party checks

- [ ] Confirm the exact versions and access dates for CASTpFold,
  ProteinsPlus/DoGSiteScorer/DoGSite3, CavityPlus, SiteMap, docking, and MD
  workflows where source records permit.
- [ ] Complete a final record-level redistribution review against
  `data/metadata/asset-rights.tsv` and `data/manifest.tsv`.
- [ ] Confirm that every transformed PDB entry has source attribution and the
  intended RCSB/PDB terms.
- [ ] Confirm that raw service bundles, commercial binaries/logs, trajectories,
  restart files, checkpoints, credentials, and private job identifiers are not
  present in the Git payload.
- [ ] Confirm any external archive deposit, completeness, checksum, persistent
  identifier, and access condition in `data/external_archive_manifest.tsv`.

## Version and archive checks

- [ ] Confirm that `pyproject.toml`, package `__version__`, `CITATION.cff`, and
  `release/manifest.json` all state version `1.0.0`.
- [ ] Confirm intended tag `v1.0.0-manuscript` points to the exact verified
  commit; do not move the tag after archival deposit.
- [ ] Freeze the exact Python 3.11 release environment in
  `environment/constraints.txt` and record its hash.
- [ ] Pass the complete suite on Python 3.11 and 3.12.
- [ ] Validate `CITATION.cff` with the pinned `cffconvert --validate` command.
- [ ] Pass `scripts/release_check.py` with both data and release manifests.
- [ ] Reproduce and verify the full bundle from a clean clone.
- [ ] Confirm the clean clone contains no untracked file required for success.

## DOI and release sequence

1. Create draft archival records and reserve the code and data DOIs.
2. Add only verified identifiers to README, CFF, and the release manifest.
3. Re-run all validation from a clean clone and record the exact commit.
4. Create `v1.0.0-manuscript` on that commit and publish the GitHub release.
5. Archive the exact tagged release.
6. Verify archive checksums, metadata, and DOI landing pages.
7. Add the manuscript DOI later, if necessary, through a traceable metadata-only
   update.

## Final sign-off

- [ ] Scientific lead approval.
- [ ] Data/rights review approval.
- [ ] Corresponding-author approval from Haopeng Sun, Yao Chen, and Xin Xue.
- [ ] Release operator records the final commit, tag, archive records, and UTC
  release time in `release/manifest.json`.
