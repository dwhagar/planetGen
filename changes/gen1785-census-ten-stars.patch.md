**Headline:** The class census draws ten stars a cell, so one star's age no longer decides a cell's class and life shares.

### Changed
- **The class census spreads a cell over ten stars (GEN.1785).** A cell drew one star and put all 300 planets round it, so the whole cell had that star's age (0.159 Gyr for the G2V ecosphere cell on main: 50 of its 51 living bodies at abiogenesis). It now draws ten stars of 30 bodies each, each on its own seeded stream. The cell's moons built bare are made whole and each body gets the end of its making (`StarSystem._finish_planet`, `_finish_moon`), as a system's bodies do. A body with no class letter is in no class share. A cell the reference drew nothing in (no moon of an M dwarf's hot planets) expects none, and a body there is a flag.
- **The system census counts promoted belt bodies apart.** `planets` is the slot planets; `promoted_count` is the belt members made planets (GEN.1510), with its own range.
- **Both censuses are renewed once from the new draws** (commit and sample in each file's reference); no range was widened to cover a fault.
