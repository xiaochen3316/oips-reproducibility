# Placeholder update checklist

After the public repository and archival records exist, update only verified
values:

- `README.md` and `README.zh-CN.md`: repository URL, release tag, commit, and
  archive identifiers where appropriate.
- `CITATION.cff`: add valid `repository-code`, release date, and DOI fields;
  omit fields until valid rather than inserting placeholder URLs.
- `release/manifest.json`: repository URL, code/data/manuscript DOI, final
  clean commit, release status, and verified tag.
- `docs/code-and-data-availability-template.md` and
  `docs/journal-statements.md`: replace `REPOSITORY_URL_PENDING`,
  `RELEASE_TAG_PENDING`, `COMMIT_SHA_PENDING`, `CODE_DOI_PENDING`, and
  `DATA_DOI_PENDING`.
- `PRE_PUBLICATION_CHECKLIST.md`: mark only independently verified checks.
- `FIGURE_SOURCE_DATA_INDEX.tsv`: replace `PENDING_AUTHOR_CONFIRMATION` with
  final manuscript figure mappings.

The archive label `9c2d0b9` is not treated as a verified Git commit.
