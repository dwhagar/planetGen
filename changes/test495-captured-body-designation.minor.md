**Headline:** A body that a planet captures keeps its belt designation and gains a link to the planet, and a body from outside the belts prints the planet's ID and a rank.

### Added
- **Designation of captured bodies (TEST.495).** `object_uid.captured_member(belt_id, rank, planet_id)` returns the belt member's unchanged designation `<belt ID>.<rank hex>` with a parent link to the capturing planet of the same system; `outside_reservoir_member(planet_id, rank)` prints `<planet ID>.<rank hex>` for a body from a reservoir that is no belt of the system; `parse_captured` reads either back through `parse_id`. All are pure; a capture across systems, by the belt itself or by a system's own ID is refused. Nothing is stored or shown yet.

No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.495.
