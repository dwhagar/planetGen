**Headline:** The columns of the orbital event tables now declare their units.

### Changed
- The 36 numeric columns of `orbital_events` and `orbital_event_bodies` (DB.56) are declared in the column registry (`db/column_units.py`): kilometres, km/s, kilograms, parsecs, joules, J/kg (a new local unit with its `_j_per_kg` suffix), Jupiter masses (`_mj`), solar masses, ratios, counts and flags. `hydrogen_limit_sol` is a mass although its name has no `mass`, so it is a suffix exception. Nothing in the database changes.
