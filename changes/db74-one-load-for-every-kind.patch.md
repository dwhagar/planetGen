**Headline:** The star, planet, moon, belt and comet loaders read their rows back with the same field lists that write them, and a loaded star keeps its stored metallicity.

### Changed
- **`object_rows.load_dict` reads a stored row back through the kind's field list (third slice of DB.74).** Each `Field` now names the `from_dict` keys it fills and inverts its own conversion (`scaled` divides kilometres back to AU, a NULL reflex offset comes back as 0.0, a 0/1 flag as a `bool`). `_star_row_to_dict`, `_black_hole_row_to_dict`, `_neutron_star_row_to_dict`, `_planet_or_moon_row_to_dict`, `_belt_row_to_dict` and `_comet_row_to_dict` are built from `STAR_FIELDS`, the planet and moon field lists (now `BODY_FIELDS` and `PLANET_ONLY_FIELDS`, replacing `_BODY_COLUMNS`, `_PLANET_ONLY_COLUMNS` and `body_row_values`), `BELT_FIELDS`, `COMET_FIELDS` and the remnants' own columns, so a column is spelled once for the insert, the in-place edit and the load. `db.edits` updates a body from the same list. The dicts loaded are the same as before for random rows of every table; the rows written are the same for seeded systems.
- **A loaded star keeps its stored `metallicity_fe_h`.** The old star loader never read the column the insert writes, so a star read back from the database carried none. The new loader returns it.
- `tests/test_object_rows.py` writes seeded stars, planets, moons, belts and comets and reads them back to the dict `to_dict` makes, and checks the helpers' inverses.

TODO items resolved: none (DB.74 stays open: the core and profile tables follow). Added: none.
