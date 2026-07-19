# Ready for GitHub review

The repository is technically reproducible locally, but public release remains
**READY_WITH_MANUAL_REVIEW** rather than unconditionally ready.

Before upload, authors must:

1. confirm the transformed-PDB and normalized third-party-data review;
2. confirm author metadata and final manuscript artwork/version mapping;
3. create and verify the repository/archive identifiers; and
4. run the complete suite on Python 3.11 and 3.12 through GitHub Actions.

The supplied converted archive has been compared with the earlier snapshot,
and the two manuscript analyses previously missing from the public package are
now implemented and frozen. These are no longer release blockers.

Never upload the supplied ZIP/7z, a virtual environment, cache, raw proprietary
outputs, trajectories, checkpoints, credentials, or internal planning files.
