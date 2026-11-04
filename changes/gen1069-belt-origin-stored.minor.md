**Headline:** A belt now keeps where it came from: what formed it, what the giants cleared from it, the sentences of what shaped it, and its largest body, stored with the belt.

### Changed
- **Belts store their origin and their largest body (GEN.1069, schema v204, migration 0204).** `asteroid_belts` gain `origin_cause`, `co_causes`, `formation_snow_line_au`, `primordial_mass_earth`, `cleared_mass_earth`, `scattered_mass_earth`, `stirred`, `origin_history`, `origin_sentence`, `scatter_fates` and `largest_body_km`. A belt read back from the database says the same things as one just generated, so the belt page can show where a belt came from and how heavy its largest body is. The largest body is the one settled with the belt's mass (GEN.1000), and a belt cut back to the orbits its star holds takes the largest body of the shorter belt with the same draw. Nothing about which belts exist or their masses changes; a belt made before v201 has these empty until the galaxy is regenerated.

generation-output: unchanged

TODO items resolved: none.
