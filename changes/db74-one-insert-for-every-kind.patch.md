**Headline:** Ten hand-written row inserts became one insert driven by a list of fields per kind.

### Changed
- **`db/object_rows.py` writes any kind's row from a list of fields (second slice of DB.74).** A `Field` names the columns it fills and how to read them from the object and the call's context; `insert_row` makes the `INSERT`. `store.py` declares `STAR_FIELDS`, `PLANET_FIELDS`, `MOON_FIELDS`, `BELT_FIELDS`, `COMET_FIELDS`, `BLACK_HOLE_FIELDS`, `NEUTRON_STAR_FIELDS`, `ROGUE_PLANET_FIELDS`, `INTERSTELLAR_COMET_FIELDS` and `ASTEROID_FIELD_FIELDS`, and the ten `insert_*` functions keep only what is theirs: the name, the child rows, the moons of a planet. The spin, the galactic orbit, the placement and the movement vector of a free object are fields shared by the kinds that have them. The rows written are the same: a run over seeded stars, planets, moons, comets and phenomena gave the same columns and values before and after.
- `tests/test_object_rows.py` checks that each kind's fields name columns of its table once and fill every column the table requires, and the constructors and the anchored remnant (its motion and spin left to its star).

TODO items resolved: none (DB.74 stays open: the generic load and the core and profile tables follow). Added: none.
