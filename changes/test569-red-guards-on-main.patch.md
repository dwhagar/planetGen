**Headline:** The guards that were red on main are green: reproducible maths in five places, the unit registry for mass conversions, a sector that no longer carries a memory address, and the sector page's map text.

### Changed
- `detmath.hypot`, `dist` and `matmul` replace `math.hypot`, `math.dist` and `@` in `generation/rogue_landing.py`, `physics/ring_darkening.py`, `physics/system_bodies.py` and `galaxy/acceleration.py`, so the last bit does not depend on the machine. The pinned galaxy is unchanged by them (checked on the commit before DOC.86).
- 17 hand-written mass conversions (`api/routes.py`, `galaxy/potential.py`, `physics/belt_ledger.py`, `field_mass.py`, `giant_heat.py`, `giant_interior.py`, `hydrosphere.py`, `life_readings.py`, `surface_air.py`) go through `unit_registry.converter`, bit for bit; the stale `open_clusters.py` entry is off the kept list.
- A system keeps its stars' architectures in a list and not in a dictionary keyed by `id(star)` (`StarSystem.architectures`, GEN.1111): the address in the key made a sector differ between runs, so `test_a_sector_does_not_depend_on_hash_seed_or_locale` failed.
- `util/secret_key.py` may import `secrets` (the form-signing key never feeds a draw); the sector page test follows DOC.86's serial comma.

### Known
- `test_the_pinned_galaxy_is_the_pinned_one` stays red on main until train 10 rewrites the golden: DOC.86 (#4786) changed generated sentences, which are stored in the sectors.

generation-output: unchanged
TODO items resolved: TEST.569 (all but the pinned galaxy, which waits for train 10's golden)
