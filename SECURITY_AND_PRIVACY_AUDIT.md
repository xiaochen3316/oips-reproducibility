# Security and privacy audit

## Result

No live API key, GitHub token, password, private key, authorization header,
commercial license key, or machine-specific absolute path was found in the
public candidate. Pattern literals used by security tests and the release
scanner are test/rule definitions, not credentials.

The following intended public author information requires final human
confirmation before upload:

- author names and affiliations in `CITATION.cff`;
- corresponding-author email addresses in `CITATION.cff` and journal-facing
  documentation.

The public tree contains no `.env`, private-key file, virtual environment,
cache, raw archive, DOCX, or commercial software project format. The release
scanner performs an additional fail-closed scan of tracked files in CI.

Status: **no secret-related block to local Git preparation**.
