**Headline:** The last six columns without a declared unit now have one, read from the code that uses them.

### Changed
- `atm_density` on planets and moons is kg/m³, `atm_molar_density` kg/mol, `asteroid_fields.formed_at` a Unix time and `comets.parabolic_mean_anomaly` a plain number (DB.86); the `unspecified` marker is removed from the column registry and a test fails if a column goes without a unit. No stored value changes.
