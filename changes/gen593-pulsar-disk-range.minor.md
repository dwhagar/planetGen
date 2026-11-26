**Headline:** A millisecond pulsar's disk planets can now be as heavy as the range says (up to 10 Earth masses, PSR B1257+12's heaviest is 4.3), and the two pulsar-planet origins that are not built say what they wait for.

### Changed

- The disk origin drew its mass range cut at 3.3 Earth masses, where a fixed 10,000 km radius ceiling at a mean density ended it; the cut is gone (GEN.593, `PULSAR_DISK_MAX_RADIUS_KM` removed, `pulsar_planets.disk_mass_range_earth` removed). A disk planet is drawn from 0.01 to 10 Earth masses whole and its radius follows from its mass by the terrestrial law (a 10 Earth mass one is about 11,800 km, class S, carbon-rich).
- The captured-giant and ablated-companion origins stay unbuilt, now named with what they wait for (`tuning.PULSAR_PLANET_PREREQUISITES`: GEN.163 and GEN.555), with the stand-in written in the module: a bare pulsar, logged with the origin and its prerequisite.

generation-output: changed
generation-affects: fill
