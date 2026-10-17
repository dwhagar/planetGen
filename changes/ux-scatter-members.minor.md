**Headline:** A star or phenomenon the scatter drew with company now shows its members and its total mass before its sector is generated.

### Added

- New page `/systems/uncharted/<id>`: an uncharted star's data, its members (kind and mass) and the orbits that bind them (separation, period, eccentricity). A star already built into a system sends the visitor to the system.
- New `GET /api/uncharted-systems/<id>` with the star, its `system_id` and its `members` and `nodes`.
- A Members column in the Uncharted stars table and the Phenomena table ("3 members, 2.4 solar masses", a dash for a single object), a members and orbits list on the uncharted phenomenon page, and a label on the Galaxy Map and in an unfilled neighbour's panel for a star with company.
- `member_count` and `system_mass_solar` on the scatter rows the uncharted-systems, phenomena, bright-star and uncharted-cell routes return.

TODO items resolved: UX.114 / Added: none
