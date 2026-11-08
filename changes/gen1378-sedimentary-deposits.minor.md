**Headline:** A rocky world with an ocean and a history of life now also carries its sedimentary and biogenic deposits (iron formations, phosphorite, evaporite, limestone, chert, coal and kerogen), in the same table as its magmatic ore.

### Added
- **Seven sedimentary and biogenic deposit types (GEN.1378, `physics/sedimentary_deposits.py`).** After the life pass each body is assessed for banded iron formation (a neutral ocean from the first oxygenic photosynthesis, for about 1.2 Gyr), phosphorite (an open ocean with life and at least 1 umol/L of phosphate from the phosphorus budget, the grade capped by the crust's P2O5), evaporite (an open ocean with land that dries at 270 K or more, graded by the dominant salt of the body's own brine), limestone, chert (from complex cells), coal (flora of 3 or more on a wet land, from multicellular life) and kerogen. Each is a row in `planet_deposit_types`, `moon_deposit_types` and `rogue_planet_deposit_types` with status, grade, endowment and the window it formed in. A body whose life never reached the stage a type needs has the type `absent_by_state`. Stand-ins SD1 to SD9 are in section 6.222 of the gap-filling plan.

### Changed
- `ore_deposits.poisson` and `ore_deposits.lognormal_sum` are public so the two modules share the count and the size law. The stored deposits of every rocky body gain seven rows and the next epoch is needed.

generation-output: changed
TODO items resolved: GEN.1378
