### Changed
- **Speeds and time periods are always shown in a meaningful unit**
  (UX.13, UX.14). Speeds go through one shared ladder, km/h, km/s, Mm/s,
  then multiples of c from a tenth of light speed ("29.8 km/s",
  "4.5 Mm/s", "0.25 c"); periods and durations through another, µs, ms,
  s, minutes, hours, days, years, ky, My and Gy ("27.3 days",
  "1.88 years", "236 My"), each to three significant figures. Orbital
  periods and speeds of planets, moons, comets, binaries, wide binaries
  and facilities, galactic orbits, pulsar spin periods, NAV travel
  times, the facility form's live readout and the admin pages' uptimes
  all use them, in Python (`stellarObjects.utils.format_speed_kms`,
  `format_duration_seconds`, `format_period_years`) and in the browser
  (`static/speed.js`, `static/period.js`). The old "x years y days z
  hours" period text (`years_to_time_string`) is gone.
