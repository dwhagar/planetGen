**Headline:** a free-floating brown dwarf no longer carries a one-in-five "has a captured moon" flag, and a free planet's chance of a moon is no longer a fixed table figure: the moons of a thrown-out planet are the ordinary moon routine's, inside the edge its ejection left.

### Changed

- `generation/phenomena/hosts.py` (GEN.558): the `planet-satellites` host table and `ROGUE_PLANET_MOON` (the fixed 20 percent) are removed with no shim; the planet host class names `Planet.add_moons` as what makes a planet's moons, a bound planet's outer limit being its Hill radius and a free planet's the edge `starless_host.build_rogue_body` sets from what the ejection left (reused: `ejection_events`, `rogue_moons`, `clearance_ledger.moon_hill_fraction`; nothing new is drawn for the physics rogues, whose moons were already `moons` rows).
- `generation/phenomena/rogue.py`, `describe/phenomena.py`, `describe/sheet.py`, `describe/stored_phenomena.py`, `web/system_pages.py`, `db/store.py`, `db/column_units.py`: the brown dwarf's `has_moons` attribute, its sentence ("a smaller companion body, likely captured after ejection"), its sheet row and its stored value are gone. A brown dwarf is 13 to 80 Jupiter masses, above the planet chain's ceiling (13), so it has no moons until the star-like chain of the unified design (free brown dwarfs with planet tables) makes its own; the draw it spent on the flag is not spent, so the later draws of a brown dwarf move.
- The chance that a thrown-out planet keeps a moon is now an output: in a run of 40 donors (G, K, M and F stars) 45 rogues were made and 4 kept at least one moon, about 9 percent (the 20 percent of the table was a design default; no source replaces it, so 20 percent is a loose upper test target only).
- Migration 0244 drops `rogue_planets.has_moons` (v207 left it NULL for the physics rogue); `db/models.py` and `schema.sql` follow, with the v243 old-schema fixture and a migration test.
- Tests: the host-table and fill-stage tests no longer read the removed table; the stand-in rows of the host table file lose the planet row.

generation-output: changed
generation-affects: fill
