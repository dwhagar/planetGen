**Headline:** the tidal heating of a planet or moon now counts every body that pulls on it, and the spin torques of the other bodies use the exact mean of their distance to the sixth power.

### Changed

- `tidal_heating.moon_tides_of` (GEN.1503, GEN.1446 part 3) adds to a planet's moons, and to a moon, the tide of every companion star and sibling planet (and for a moon the host star at its planet's place and the other moons of its planet), each against the body's spin as the moons' tides already were. The effect ledger notes every source's heating power (`tide_sources.heating_by_source`, at a nominal `k2/Q` so the shares are what is kept).
- `tide_sources` takes the torque and heat of such a source from the exact mean of `1 / d^6` of two circles instead of the square of the mean of `1 / d^3` (registry row TS6 rewritten); that also corrects the spin torques of GEN.1502 for close neighbours.
- Census over 30 systems (98 bodies): tidal flux moves on 39 bodies (most from values near zero; two close planets 0.013 AU apart go from 5e-27 to 5e-5 W/m2) and internal heat flux on 22, by under 0.1 percent; temperatures move on one body, by 1e-6 K.

generation-affects: fill
