**Headline:** The habitability levels page now takes every number from the limits the score uses, so it cannot say one thing while the rating does another.

### Changed
- `web/lib/habitability_ref.py` builds the Pressure, Temperature, Air, Water and Radiation tables from `habitability.HUMAN` (the `protection_limits` values) and the gas, water and chaotropicity constants. A test fails when a number on the page differs from its constant, and when a changed limit does not move the page.
- Oxygen moves from Pressure to Air, where the score reads it, and is split into too little (a mask works) and too much (50 to 160 kPa). Carbon monoxide, hydrogen sulfide and sulfur dioxide each get their own row with their own numbers, water gets a chaotropicity row, and the equipment rule names the 140 C heat-shielded-suit limit.
- Docstrings that called the overall rating PHI-4 say PHI-6; the stored key `"PHI-4"` and the `phi4` columns stay.

TODO items resolved: GEN.219 / Added: none
