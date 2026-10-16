**Headline:** Any sector in the grid's reach can be asked for and is generated from a profile of its place, in the disk or out in the halo.

### Changed

- New `planetgen.galaxy.profile`: the sector profile, a pure function of the galaxy shape and a cell's address (no database, no random number) that returns the expected systems, the mix of stellar populations with each one's age window, [Fe/H], velocity dispersion and rotation, the cold-gas fraction and ambient medium, and the dark matter in the cell. `FillContext` takes a profile.
- The reach (shell 65,535 by default) is the only refusal. `ensure_sector_generated`, `planetgen galaxy --shell --band --wedge`, `--shell --band` and named `--shells`/`--bands` spans generate and save any cell of the reach from its profile, empty or not, with no warning and no `--strict` stop. `outside_ok` is gone, and `ensure_sector_generated` reports `in_reach` next to `qualifies`. A bare `--shell I`, a span with no range, the random start and the maps' planned sectors still mean the outline.
- A neighbourhood around a halo sector takes every cell of the reach within the radius and says how many sectors that is, and how many systems they expect, before it runs; around an outline sector it stays inside the outline as before.
- The `--shell/--band/--wedge` form stops with a message, not a traceback, for a band its shell does not hold.

TODO items resolved: GEN.235. Added: none.
