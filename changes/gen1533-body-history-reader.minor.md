**Headline:** A planet or moon loaded from the database now carries the birth eccentricity, spin and capture draws that its history left, so the stored galaxy has the numbers a run has.

### Added

- GEN.1533: `store.load_star_system` reads `body_history` (GEN.1497, written since v243 and read by nothing until now) through `physics/body_history.attach`, which sets `_orbit_birth` and `_birth_spin` on each restored planet and moon, the two names a generated body keeps them under. The spin is attached whole or not at all, and a row written under another `FORMULA_VERSION` is not read. Reuse: the table, the column names and the body codes are GEN.1497's. This changes no stored value and no generated world: nothing that runs on a restored body reads the two names yet (`planets.draw_birth_spin` still draws, which keeps the edits as they are), so the passes that will (the air history of GEN.1739, and a moved body's spin history, a later slice that is a lock item because it changes what an edit produces) start from the real birth numbers instead of a fresh draw. The capture draws are a `FLOAT` column and come back to six digits.

generation-output: unchanged
generation-affects: none
