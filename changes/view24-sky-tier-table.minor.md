**Headline:** The generator can now say how deep the star backfill must reach, and how much of it fits, to store the sky a system sees at a given limiting magnitude.

### Changed

- `generation/sky_tiers.py` (VIEW.24, the sky tier table): `floor_luminosity_sol` gives the dimmest star that still shows from a distance at a limiting magnitude (the star of the best colour, through the dust), `extinction_curve` samples the dust of the galaxy along 48 lines of sight from an observer, `reach_pc` gives how far a star of a given luminosity shows (the clearest line), and `build` gives the table: the luminosity bands from 50 Lsun up and the mass rings, each with the radius of its sphere and the cells in it, cut by a budget of cell visits (`CELL_VISIT_BUDGET`, 3e8: at a 7.5 limit it keeps the bands under about 1,000 Lsun). `floor_table` lists the floor against distance and which ledger draws it. Nothing reads it yet; the charting step (VIEW.25) will. Registry rows SKT1 and SKT2 (6.176).

generation-output: unchanged
