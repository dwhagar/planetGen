**Headline:** The gas a giant planet gathers now has the water and CO2 that drifting pebbles evaporated into it, so fewer giants come out carbon-rich.

### Changed

- The carbon-to-oxygen ratio of a gas giant's envelope is now set by the disc gas at the age the giant gathered it (GEN.1473, `physics/disk_drift.py`): pebbles drift inward faster than the gas and evaporate at each ice line, so the gas inside the water line is oxygen-rich while the pebbles last and the gas inside the CO2 line gets CO2 once the pebble front has reached it. The age comes from the star's mass, metals and stored snow line and the giant's formation temperature; there is no new random draw.
- In a sample of 2,787 generated giants the share with C/O over 0.8 falls from 18 to 12 percent (Sun-like stars) and from 18 to 11 percent (K stars), the median stays where it was, and the share under 0.35 rises by about 2 points; M dwarf giants change the most (median 0.62 to 0.44). The measure is in `reports/gen1473-drift-measure-2026-11-09.md`.
- Stand-ins DE1 to DE10 are in section 6.242 of the registry (the pebble flux and front are read from the paper, the end of the pebbles, the viscous gas flux, the lag and the gas-gathering age are registered stand-ins).
- Generated giants' atmospheres change, so this is a new epoch; planets, orbits, masses and everything else are unchanged.

generation-output: changed
generation-affects: fill
