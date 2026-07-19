# Test report

- Date: 2026-07-18
- Platform: Linux, x86_64
- Python: 3.12.13
- Isolated test environment: successful (`venv --copies` outside repository)
- Editable package installation: successful
- CLI help: successful
- Input validation: 6 checks passed, 0 failed
- Pytest: **216 passed, 0 failed, 0 skipped**
- Main frozen snapshot: **24 passed, 0 failed**
- Toy frozen snapshot: **22 passed, 0 failed**
- Release scanner: **183 tracked files checked, 0 failures**
- Public-release blocker from tests: no

An initial environment-only attempt failed because the host-created in-tree
venv contained unusable runtime symlinks; that environment and all generated
caches were removed. A clean copied-interpreter environment then installed the
declared dependencies and passed the complete suite. The full 21-target
reproduction generated 44 verified files, including 20 analysis tables; the
independent toy workflow generated the same stage structure.
