**Headline:** Every release note now has to open with a Headline line, so no release goes out without one.

### Changed
- `scripts/bump_version.py` refuses a release note that does not start with a `**Headline:**` line (the rule has been optional since it was restored on 2026-10-14). `--check`, `--check-pr` and the stamp all apply it; `changes/README.md` says so.
