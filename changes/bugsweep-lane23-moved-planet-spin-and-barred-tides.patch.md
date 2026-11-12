**Headline:** a planet that is moved to resolve an overlap keeps the spin it was born with, and a body a system drops for the bar on habitable worlds no longer pulls tides on the bodies that stay.

### Changed

- `generation/validation.reconcile_moved_planet` (GEN.1597): the planet and its moons start the spin history from the birth spin they already have (`birth=`) instead of drawing a new one, and the orbital period is set before the air is solved, not after. Stored values change for a planet that was moved (its day length and tilt, and the draws after it).
- `physics/planets.generate_moons` (GEN.1599): the disk-full rejections stop at the call's `rejection_limit`, as the too-large ones do. Nothing is stored differently (a disk exists only round a giant, whose limit is the same number).
- `generation/system._drop_barred_bodies` (GEN.1604) now runs before the tide pass, so the tides, the heating and the effect-ledger rows of the bodies that stay no longer count a dropped planet or moon (seed 433 of a bar-habitable system had 40 such rows). Stored ocean tides and the ledger of systems that drop a body change.
- `physics/ocean_tides.py` (GEN.1609): confirmed as designed, no value changes: a hycean ocean is deep under a hydrogen envelope, not the thin layer the tide formula holds for, so its tide columns stay empty; the reason is now in the module text and a test pins it.

generation-affects: fill
