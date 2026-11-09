**Headline:** A planet the obliquity tide has puffed up now loses its envelope to the star's X-rays at the radius it really has.

### Changed

- The loss of a hydrogen-helium envelope to the star's XUV (`envelope_escape.envelope_left_earth`) was worked out on the radius of the uninflated law, and the tide's factor of up to 2 was put on the radius afterwards. The rate goes as about R^4 to R^6 on the XUV branch and R^10 on the boil-off one, so an inflated planet lost far less than it should. The loss now asks for the tide's factor at every step (`radius_ratio`), so a tilted, hot planet keeps less envelope than the same planet upright. A planet with no tilt, an upright one, a cold one and every caller outside the settle pass work out the same numbers as before.
- On 300 seeded systems (1,217 planets with a captured envelope) 2 end tide-inflated and one of them is old enough to lose envelope: 0.052 to 0.037 Earth masses left. The output moves for well under 1 percent of the envelope planets today and for more if the tilt of GEN.1447 makes big tilts common.
- Whether a planet was inflated (`_tidal_radius_km`) is read from the tide's factor that `envelope_capture.kept_with_tide` returns, not from a second, uninflated envelope: with the loss on the tide, that second envelope is bigger than the tide's and made the comparison wrong.

generation-output: changed
generation-affects: fill
