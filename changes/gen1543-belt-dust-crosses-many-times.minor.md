**Headline:** A planet's belt-dust influx on its sheet now counts the many times a drifting grain crosses its orbit, so the Earth gets about a million kilograms a year from the asteroid belt alone instead of 45.

### Changed

- `physics/influx.py` (GEN.1543): a grain of belt dust does not cross a planet's orbit once. Poynting-Robertson drag takes it inward over the course of tens of thousands of orbits, and it crosses the planet's orbit radius twice an orbit while its semi-major axis is within its eccentricity of the planet's. The planet's share of each pass was already `R^2 / (4 a^2)`; the crossings (`crossings`, `e c sqrt(a) / (pi beta sqrt(G M))`, about 3e4 at 1 AU) multiply it, saturating at the whole supply. A comet nucleus still crosses once.
- Check: a Solar-System asteroid belt (2.4e21 kg, 4.6 Gyr) gives the Earth about 7.6e5 kg a year, against 45 before and about 1e7 measured; the Jupiter-family comets, which are the larger share of the real dust, are not in that belt. Stand-ins IP1 (grain beta 0.01) and IP2 (grain eccentricity 0.1), IF3 reworded.
- Only the derived row on a planet's sheet changes (the star's luminosity now sets the grains' drift); nothing stored moves, and the late-accretion pass, which reads `influx.reach` and `shield` only, is unchanged.

generation-output: unchanged
generation-affects: none
