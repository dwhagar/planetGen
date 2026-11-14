**Headline:** A planet's belt-dust influx on its sheet now counts the many times a drifting grain crosses its orbit (the Earth gets about a million kilograms a year from the asteroid belt alone, not 45), and the Venus-cloud species V8 reads that influx from the system's own belts.

### Changed

- `physics/influx.py` (GEN.1543): a grain of belt dust does not cross a planet's orbit once. Poynting-Robertson drag takes it inward over the course of tens of thousands of orbits, and it crosses the planet's orbit radius twice an orbit while its semi-major axis is within its eccentricity of the planet's. The planet's share of each pass was already `R^2 / (4 a^2)`; the crossings (`crossings`, `e c sqrt(a) / (pi beta sqrt(G M))`, about 3e4 at 1 AU) multiply it, saturating at the whole supply. A comet nucleus still crosses once.
- Check: a Solar-System asteroid belt (2.4e21 kg, 4.6 Gyr) gives the Earth about 7.6e5 kg a year, against 45 before and about 1e7 measured; the Jupiter-family comets, which are the larger share of the real dust, are not in that belt. Stand-ins IP1 (grain beta 0.01) and IP2 (grain eccentricity 0.1), IF3 reworded.
- V8's dust supply now reads the system's own belts and comets (`generation/dust_supply.apply`, one call before the late accretion pass, which keeps `_dust_influx_kg_yr` on each rocky planet and moon; a moon takes its planet's distance with its own radius). The gate is an absolute floor of `dust_min_influx_kg_m2_yr` (2.4e-10, a tenth of Venus's by this model with a Solar-System asteroid belt; stand-in CV14). The fixed-Sun `cosmic_dust` law of GEN.1542 and its stand-ins CV11 to CV13 are gone.
- V8 is still off by default, so no stored value changes with the defaults. Only the derived row on a planet's sheet changes (the star's luminosity now sets the grains' drift); nothing stored moves, and the late-accretion pass, which reads `influx.reach` and `shield` only, is unchanged.

generation-output: unchanged
generation-affects: none
