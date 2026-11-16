**Headline:** The belt mass module imports on its own again.

### Fixed
- GEN.1633: `physics/belt_mass.py` failed to import in a fresh interpreter (a circular import through `collision_outcome`, `planets`, `generation/irregular_moons`, `irregular_population` and `belt_ledger`). `physics/planets.py` now imports `irregular_moons` where it uses it, and `tests/test_physics_imports_alone.py` covers `belt_mass` and the four modules of the cycle. Nothing generated changes.

generation-output: none
