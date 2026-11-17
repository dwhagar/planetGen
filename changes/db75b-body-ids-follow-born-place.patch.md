**Headline:** A body's public ID now follows the place it was born in, so a redraw that trims a planet or a moon leaves every body it keeps with its ID.

### Changed
- **A planet's slot is its place in the star's plan, a belt's its rank, a moon's its place among the planet's born moons (DB.75b).** The slot an ID follows (`store._body_slots`, `body_numbers`) was the body's index in the list of its star's orbits or its planet's moons, so a draw that kept one planet fewer from the middle of the plan moved every planet after it up one place and each took the ID of the one before. A planet now carries `plan_slot` (the slot of `planet_architecture.Plan.seed` it was made in, which no trim moves), kept as its `body_history.born_index` (the column already held a moon's place in its born moons, so there is no migration) and given back when the system is loaded; the slot is `planet:<star>:orbit:N`. A planet promoted from a belt is `planet:<star>:belt:O:R` (the belt's number and the member's rank, `planets.belt_order` and `belt_rank`), a belt is `belt:<star>:rank:N` (its rank among the star's belts), a moon is `moon:<planet slot>:born:N` (`body_history.born_index`). A body stored with no such number (a galaxy made before it, a planet a pulsar has) is placed by its index as before. A fresh system numbers its bodies in the order it always did, so no ID changes in a galaxy made by this version. The key stays the star's role, not a seed: the slot-path form of GEN.556 will read the host chain.

### Fixed
- **A redraw dropped the `body_history` and `lost_moons` rows of the new bodies.** `replace_star_system_content` moved the body rows to the kept system and left the history rows under the stand-in system, which deleting it removed, so a system redrawn once had no birth numbers for its next draw or edit. They now move with the bodies, and the old bodies' rows are deleted.

generation-output: unchanged
generation-affects: none
