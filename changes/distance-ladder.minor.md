### Changed
- **Every distance is shown in its most meaningful unit.** One helper,
  `stellarObjects.utils.format_distance_m` (with `_km`, `_au`, `_ly` and
  `_pc` wrappers, re-exported by `html/lib/fmt.py`), and its browser copy
  `html/static/distance.js` pick the largest of km < AU < mpc < cpc < ly
  < pc < kpc < Mpc < Gpc the value is at least 1 of. Parsec values add
  ly in parentheses ("4.2 pc (13.7 ly)"), or AU below 0.01 ly ("2.4 mpc
  (495 AU)"). The system, sector, galaxy and phenomenon pages, the
  System, Sector and phenomenon map readouts, the wiki sector page and
  the generated text (planet, star, binary, belt and heliosphere
  distances) all use it.
- **Planet, moon and star radii are always km in scientific notation**
  (`utils.format_body_radius_km`), including a rogue planet's.
- **The distance constants are exact:** AU = 149,597,870,700 m,
  lightyear = 9,460,730,472,580,800 m, parsec = 3.085677581491367e16 m,
  with every conversion derived from them. Stored values shift by about
  one part in 70,000.

### Removed
- The unused `LY_THRESHOLD`, `HELIOSPHERE_DISPLAY_THRESHOLD_LY` and
  `ROUND_HABITABLE_ZONE_AU(_SMALL)` display constants.
