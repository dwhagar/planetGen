**Headline:** New functions read from the temperature of every place how much of a world is livable, how hot and cold it gets, and where a base would go (no stored value changes yet).

### Changed

- `physics/temperature_coverage.py` (GEN.1264) reads a `TemperatureSeries` (GEN.1263): the shares of area and time in the Blue, Green, Yellow and Red bands (cold and hot counted separately), the one livable share (the Green band or better, area times time) with its two pieces (the share of the world livable most of the time, and the time a typical place of it is livable), and the areas livable all the time, nearly all the time and at some time.
- The coldest and hottest hour anywhere, by latitude band, by longitude, by quarter of the orbit, and, for a world that turns once an orbit, for the day side, the terminator belt and the night side (with the livable share of each, so a world livable only in its belt can be told).
- The best site (the cell livable the largest share of the time, a tie going to the more moderate) with its own coldest and hottest hour and the longest unbroken spell out of the livable band there (how long a person would stay indoors).
- The wet bulb rule is applied to each sample. Three definitions are registered in GEN.725 (rows SS12 to SS14 of section 6.91). Nothing calls it yet (the score, columns and migration are GEN.1265 and GEN.1266), so no generated output changes.
