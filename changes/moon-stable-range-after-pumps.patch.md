**Headline:** A moon that a Kozai pump, a resonance or the replay of the moons' history left outside its planet's stable range is now lost, so a generated system no longer fails its own lunar check.

### Fixed

- The settle pass of a system drops the moons whose orbit lies outside the stable range of the eccentricities the pumps gave (`planets.drop_unstable_moons`, after the tide has worn them as far as it will). Nothing checked them before: the planet's pump and the moons' pump raise the eccentricity of a planet and its moons, the Domingos limit falls with both and can reach 0, and the moon stayed (found on the Train 61 gate as `test_wide_binaries_never_overlap`: a moon at 1.3e4 km in a stable range of 6765 to 0 km).
- The replay of the moons' history (`moon_history._settle`) holds each moon to the same limit as the lunar check, the semi-major axis the fit allows for the moon's own eccentricity, so a moon the drop took cannot come back when a planet's end stages are made again (`refinish_rocky`).

generation-output: changed
generation-affects: fill
