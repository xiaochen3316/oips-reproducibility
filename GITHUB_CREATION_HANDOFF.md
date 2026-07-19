# GitHub creation handoff

Recommended repository name: `oips-reproducibility`.

After the manual actions in `READY_FOR_GITHUB.md` are complete:

1. Create an empty public GitHub repository without auto-generating files.
2. Add the verified remote URL locally.
3. Replace identifier placeholders and regenerate the release manifest.
4. Run all tests, full reproduction, snapshot verification, and release scan
   from a clean clone on Python 3.11 and 3.12.
5. Commit the verified metadata update and record its full SHA.
6. Only then create the final manuscript tag and GitHub release; archive that
   exact tag and verify DOI landing pages.

This handoff intentionally contains no token, remote command, fabricated URL,
DOI, tag, or final commit.
