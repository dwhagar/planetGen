**Headline:** class M now reads the sourced human wet-bulb limit (31 C) as its upper line instead of a 28 C margin under it, so an Earth-like world, whose hottest-day wet bulb is 28 C at its median, no longer sits on the edge of the class.

### Changed

- `planet_class_fit.json`: the `wet_bulb_c` line of class M moves from 28 C to 31 C, the value O and R already use and the uncompensable limit of `protection_limits` ([R4521]). Choice: keep the reading (the wet bulb of the daytime peak of the hottest band, which the limit is a limit on) and fix the line, rather than read the band mean, because the sourced limit is a limit on the peak exposure; the daytime peak is about 15 K over the band mean, so a mean reading would have needed a lower, unsourced line. DM2 in the stand-in registry is updated (the 28 C line is gone).
- New `tests/test_wet_bulb_edge.py` keeps the lines of M, O and R equal to the sourced limit.
- Changes the class of bodies whose hottest-day wet bulb is between 28 and 34 C, so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
