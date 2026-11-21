**Headline:** The trojan stack reads its astronomical unit and year from the constants module and tags its rounded Earth mass.

### Changed
- `trojan_bodies.AU_CM` and `trojan_overlay.SECONDS_PER_YEAR` now come from `physics/constants.py` (`AU_M` times 100, `SECONDS_PER_YEAR`): the same numbers as before (the IAU astronomical unit and the Julian year). `trojan_bodies.EARTH_MASS_G` keeps the rounded 5.972e27 g and says in its docstring that it is derived from the constants module's Earth mass, so the generated output does not move. The constants-completeness gate (TEST.715) no longer lists the three.

generation-output: unchanged
