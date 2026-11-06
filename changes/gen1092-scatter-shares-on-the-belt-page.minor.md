**Headline:** The information page of an outer disk now says where the disk's mass went: how much stayed cold or got stirred, sits in a resonance, was scattered, thrown out, parked beside the giant, put in the inner belt or caught as moons.

### Added
- **Where the disk's mass went (GEN.1092).** The page of an outer disk (`/object/belt:<ID>/belt`) has a table of the groups the outermost giant sorted the disk into, each with its share of the disk, adding up to the whole. `GET /api/objects/<ref>/belt-page` carries it as `belt.fates` (a list of `{population, share}`, empty for a belt that is not an outer disk), a new field on an existing route, so the API version does not change.

generation-output: unchanged

TODO items resolved: none (GEN.1092 stays open until the shares that left the disk are handed to the inner belt, the moons and the rogue builders).
