**Headline:** `physics.hydrosphere`, `bulk_draw` and the other modules on the same import chain now import on their own, and two item codes that showed in text are cut (GEN.1401).

### Changed

- `hydrosphere` and `icy_tide` read the melting point of water from `aquifers`, the one definition, so the chain `hydrosphere` to `planet_type` to `orbit_elements` to `interior_heat` to `icy_tide` back to `hydrosphere` no longer fails when `hydrosphere` or `bulk_draw` is the first import. `tests/test_physics_imports_alone.py` imports each of the chain's modules in a fresh interpreter.
- The item codes GEN.1047 (a field description in the sheet) and GEN.1360 (a relation note in the stellar abundances) no longer show in the text the code prints.

generation-output: unchanged
