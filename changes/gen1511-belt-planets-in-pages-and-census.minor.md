**Headline:** An asteroid belt's page lists the planets it made, and the designation of a belt body that became a planet leads to that planet.

### Added
- **The belt's information page has a "Bodies made planets" section (GEN.1511).** It lists the heaviest round bodies of the belt that GEN.1510 made planets, each with its name (a link to its own page), its designation in the belt, class, mass, diameter and distance; `GET /api/objects/<ref>/belt-page` has the same list as `planets`. Their mass is not in the belt's totals.
- **The old designation still works.** `/object/belt:<ID>/member/<rank>` of a body made a planet goes to that planet's page; `GET /api/objects/<ref>/belt-members/<rank>` answers `promoted` with the planet's reference.
- **A planet made from a belt names it.** Its property sheet has a "Made from a belt" row with the belt and the designation the body had there.
- **A census of belt bodies made planets** (`tests/audit/belt_promotion_census.py`): planets by class and by star class, per belt, the generation time, and the rules of the promotion checked over a sample of systems.
generation-output: none
