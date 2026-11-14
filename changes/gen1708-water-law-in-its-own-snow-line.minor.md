**Headline:** Planets in the habitable zone are no longer born so wet: the water law is read in the snow line it was fitted in.

### Changed

- GEN.1708: the Bern water law (the dry, low, high and rich classes by formation distance and mass) takes its distance in snow lines of 2.69 AU at one solar luminosity ([R0600]), but the generator read it against the star's stored snow line, the disk-era line of GEN.1611 (1.6 AU round a Sun-like star). A planet at 1 AU stood at 0.62 snow lines where the law's own unit makes it 0.37, so habitable-zone planets drew the rich class far too often. `delivered_water.law_snow_ratio` now converts: the distance over the stored line times the Sun-like disk line (1.59 AU) over 2.69 AU (about 0.59), so the star's own line scatter and the line's change with the star's mass stay in the ratio. The stored line is unchanged and still decides whether a body is icy or rocky. The air's water-class rows (b1 to b4) read the same ratio. In 80 Sun-like systems the planets of 0.4 to 6 Earth masses between 0.7 and 1.6 AU were 6 rich, 8 high and 7 dry of 21 before and 1 rich, 8 high, 11 dry and 1 low after (all planets of that mass: dry 64 to 73, rich 134 to 124 of about 297). Registry row DW4.

generation-output: changed
generation-affects: fill
