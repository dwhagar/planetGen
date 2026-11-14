**Headline:** A stray scratch script, `m.py`, is gone from the repository root, and a test now fails if another one lands there.

### Fixed
- `m.py` (a throwaway measuring helper that came in with the GEN.1309 A commit) is removed; nothing used it.
- New `test_repo_root_holds_no_scratch_files.py` fails when a tracked `.py` file other than `setup.py` sits in the repository root (TEST.745).

generation-output: unchanged
generation-affects: none
