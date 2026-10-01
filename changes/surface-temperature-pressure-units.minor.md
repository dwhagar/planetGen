### Changed
- **Surface temperatures and pressures show customary units alongside
  metric.** A planet's or moon's surface temperature now reads in K with
  °C and °F ("288 K (15 °C, 59 °F)") in its description and on the
  System Map's info panel, and atmospheric, internal and core pressures
  read on a Pa, kPa, MPa, GPa ladder with atm and psi alongside
  ("101 kPa (1 atm, 14.7 psi)"), through the new
  `stellarObjects.utils.format_temperature_k` and `format_pressure_pa`.
  Star temperatures stay in K alone.
