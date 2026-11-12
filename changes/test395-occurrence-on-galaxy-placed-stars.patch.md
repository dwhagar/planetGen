**Headline:** The planet occurrence check now also runs on stars the galaxy places, so its host mix is the sky's (about 80 percent M dwarfs) and it reads the final classes and stored life.

### Changed

- New census `tests/audit/prevalence_galaxy.py` with the live test `tests/test_prevalence_galaxy_census.py` (TEST.395): 100 seeded sectors at the reference solar position, filled the way a galaxy sector is (`census.fill_region_sector`), give about 300 systems; each is filed under its host group by mass (`planet_occurrence.host_group`) and held to the same observed bands as the system draw's occurrence census. A group of 25 stars or more is a cell of its own with its range widened for its own count; the `galaxy` cell pools the groups that have a band for the measure and averages their ranges by their stars, so a failure names the cell and the measure that left its band.
- A flag carries the known-gap item of its group in the stored occurrence file (GEN.728, GEN.1537). On main the sample shows 0 new flags and 6 known gaps: planets of 1 to 4 Earth radii per M dwarf 0.97 against 1.9 to 3.1 observed, planets per K dwarf 0.91 against 2.2 to 4.4, and no giant at 2 to 8 AU round the 54 K and G stars.
- The report (`occurrence-galaxy-latest.md`) lists the final class mix of the sample and the bodies with life per star (1.2), report only. Classes O (one body in 1,318), P and Q are not given a target rate here.
- Tests and notes only: no generator code changes.

generation-output: unchanged
