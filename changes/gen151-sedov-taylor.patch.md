**Headline:** Supernova remnants are sized by the real Sedov-Taylor blast-wave law, so a remnant in thin gas grows larger than one in dense gas.

### Changed
- **Remnant radii ignored the surrounding gas.** Every remnant used one fixed coefficient, which is a blast in gas of about 220 atoms per cubic centimetre, although the remnant classes list densities from 0.1 to 100. The radius is now `1.03 ly (E51 / n)^(1/5) (t / yr)^(2/5)` for an explosion of 1e51 erg in the remnant's own density `n`, held at the snowplow onset radius `45.7 ly * n^(-3/7)` where radiative cooling begins. A 1,000-year remnant in 1 atom per cubic centimetre is 5.0 pc.
- The new `sedov_taylor_radius_ly` function in `generation/phenomena/supernova_remnant.py` is the single formula; the plausibility check uses it too.

TODO items resolved: GEN.151 / Added: none
