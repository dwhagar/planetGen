**Headline:** A star's colour can be worked out from its temperature, the first part of the photometry library.

### Changed
- **Blackbody colour of a star (GEN.416).** The new `physics/photometry.py` takes a spectrum to a colour: the CIE 1931 curves (the Wyman, Sloan and Shirley fit), linear sRGB (D65) with the brightest channel 1, the sRGB curve, the 8-bit colour and the hex. `star_colour(temperature_k)` gives a star's colour: a Planck spectrum from 4500 K up, and for cooler stars the BT-Settl offset from the blackbody's chromaticity with a seeded spread the caller draws. It matches the table of the gas-optics design to 0.001 in chromaticity from 3000 K up. Nothing is stored and no generated value changes; the atmosphere colour work (GEN.418, GEN.419) and the view from a planet (VIEW.11) build on it.

generation-output: unchanged
