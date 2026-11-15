**Headline:** A planet that is moved to resolve an overlap keeps the spin it was born with, and its air is solved on the orbital period it has now.

### Changed

- GEN.1597: `generation/validation.reconcile_moved_planet` drew a new birth spin for a moved planet and its moons (`generate_spin` without `birth=`) and set the orbital period after the air was solved, so the air saw the old period. It now passes each body's own `_birth_spin` (the orbit moved, the planet did not), and sets the period before the air. Stored values change for a planet that was moved: its day length, tilt and what follows from them, and the draws after it.

generation-output: changed
generation-affects: fill
