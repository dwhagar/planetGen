**Headline:** The check that deploy scripts only name files that exist no longer mistakes a record key for a folder.

### Changed
- **Deploy script path check (TEST.566).** `test_scripts_only_name_repo_paths_that_exist` read `e.setdefault("src", "crossref")` in `scripts/settle_references.py` as the path `src/crossref`. It now reads a path only when it is built by a `join(...)`, as the check was meant to. Test only; no change to any script or output.
