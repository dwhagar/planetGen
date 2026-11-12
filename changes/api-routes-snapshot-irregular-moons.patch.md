**Headline:** The route snapshot and the link audit now know the two irregular-moon routes that GEN.1086 added, which turns main green again.

### Changed

- `src/tests/fixtures/api_routes.json` lists `/api/objects/<ref>/irregular-moons` and `/api/objects/<ref>/irregular-moons/<rank>` (added without a version bump, since a new route is additive).
- `test_api_link_ids.py` names both in `NOT_AUDITED` with the reason (the planet is the ref in the path), as the belt pages are.

generation-output: unchanged
