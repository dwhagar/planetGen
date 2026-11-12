**Headline:** Every stored field of a star, planet and moon is now marked as a birth value or as present state, with tests that a remake keeps the first and works out the second again; nothing reads the list yet.

### Changed

- GEN.582 (light): new `generation/body_state_fields.py` with `STAR_BIRTH` and `PLANET_BIRTH` (a moon is a Planet) and `birth_fields`, `present_fields`, `is_birth`. A field is present unless named birth. Tests: every stored field is exactly one of the two; birth names exist; a rocky planet moved 12 percent outward and remade (`remake_after_orbit_change`) keeps every birth value and moves the present fields that depend on its orbit. Classifies only; the Field table (GEN.565, 8.2) absorbs the list.

generation-output: unchanged
