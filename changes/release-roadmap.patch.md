**Headline:** The plan set has a release roadmap from v8.0 to v10.0, proposed for Boss: a headline, the open items, what each release needs and a list of ideas not yet filed.

### Added

- `docs/plan/release-roadmap.html`, built by `scripts/build_todo_docs.py` from `docs/plan/release-roadmap.md` and `docs/TODO.md`. It groups every open item into 8.0, 8.1 (infinite zoom), 8.2 (accounts and email invitations), 8.3 to 8.8, 9.0 to 9.3 and 10.0, with open, ready and finished counts per release, the headline changes it proposes and a "Suggested, not filed" list for each release.
- Tests that every item is placed once, that no item sits in an earlier release than something it needs, and that every release has a headline.

### Changed

- This is a proposal: no `Headline:` line in `docs/TODO.md` moves, no `.major.md` note is added and no suggestion is filed as an item until Boss approves.

TODO items resolved: none / Added: none
