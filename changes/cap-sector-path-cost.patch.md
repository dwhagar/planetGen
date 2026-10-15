**Headline:** Working out the path of a body held by a heavy black hole no longer takes minutes in a core sector.

### Changed
- A body held by a heavy mass (the thousands of rogue planets that orbit a core black hole) used to run all 20,000 integration steps of its sector path, about half a second each, so a core sector with a 10 to 90 million solar-mass hole took many minutes to settle and `test_galaxy_gen` timed out about one run in three (GEN.232). A path is now cut off once its direction has turned through one full revolution (`sector_path.MAX_TURN_RAD`), the most the spline's 48 knots can follow anyway, and the hard step ceiling is 2,000. A held body now costs about 5 ms, whatever the hole's mass, and every body's cost is bounded the same way. Paths that leave the sector are unchanged to the last bit. The integrator also visits the point masses once per step instead of twice (`sector_path._pull`).

TODO items resolved: GEN.232. Added: none.
