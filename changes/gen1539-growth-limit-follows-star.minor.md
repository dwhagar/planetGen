**Headline:** The distance beyond which a star's disk has built less, and where an outer-disk belt can start, now follows the star's mass (11.4 au for the Sun, mass to the half power) instead of 11.4 au for every star, so small stars no longer get an outer-disk belt at the Sun's distance.

### Changed

- `tuning.PLANET_OUTER_GROWTH_AU` (11.36 au, a stand-in tuned on the Sun) is now the Sun's value; `physics/disk_size.growth_limit_au(M)` scales it as the gas disk's median edge is (GEN.1520), so the same 11.4 percent of the edge: 7.0 au for 0.38 Msun, 9.5 au for 0.7 Msun, 12.9 au for 1.28 Msun. It is read where the core taper starts (`planet_mass.oligarch_core_kg`) and where the oligarch chain stops and an outer-disk belt can start (`belt_birth` to `belt_origin`). The Sun is unchanged; stand-in DS5 (section 6.274 of the gap plan), tier educated guess from physics, tweak L1.
- Found in the GEN.1520 after-sweep: every spectral type showed a median belt at exactly 11.4 AU because the start was fixed. Planet and belt counts, the cores of outer planets and the places of belts change for stars that are not Sun-like; the before and after sweep is in `reports/gen1520-disk-size-by-mass-2026-11-11.md`.

generation-output: changed
generation-affects: plan, scatter, fill, population
