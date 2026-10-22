**Headline:** Every numeric column of the database now declares the unit it is stored in, in one list.

### Changed
- **The column registry (`db/column_units.py`).** The 974 numeric columns of the 63 tables are listed by table and unit (`TABLES`), each unit a key of the unit registry (`physics.unit_registry`), or of `LOCAL_UNITS` for the few it does not hold yet (volumes, rates, a dose rate, logarithmic scales): lengths in `km`, `au`, `pc`, `mpc` or `ly`, masses in `kg` or `Msun`, times in `h`, `yr`, `Gyr` or `unix_s`, and so on. A column with no unit of its own says what it is (`ratio`, `count`, `index`, `id`, `key`, `seed`, `flag`, `unit-vector`, `row-defined`); six columns holding four quantities whose code states no unit (`atm_density`, `atm_molar_density`, `formed_at`, `parabolic_mean_anomaly`) are declared `unspecified` and wait for a source. Nothing in the data moves: values and units stay as stored.
- `point_columns.PointLayout` no longer carries a position unit and a mass unit of its own; it reads them from the registry (`position_unit`, `mass_unit`), and `LENGTH_M` and `MASS_KG` are gone (the mass unit of the solar-mass columns is now `Msun`).
- `tests/test_column_units.py` fails when a numeric column has no declaration, when a declaration names no column, and when the unit in a column's name is not the unit it is declared in.

TODO items resolved: DB.76. Added: none.
