**Headline:** A planet outside a companion star's orbit is no longer lit thousands of times too brightly by that companion.

### Fixed
- **GEN.1790: the light, XUV and far-UV of a companion star were wrong for any body outside the companion's orbit.** `mean_inverse_square` averaged `1/(separation^2 - d^2)` and floored it, which is only right for a source outside the body's orbit. It now uses the exact coplanar mean `1/abs(R^2 - d^2)` for a source inside or outside, keeping a floor only very near `d = R` (stand-in RF9, a numerical guard). The Zeta Ophiuchi example, which put bodies past the ionization ceiling, now has none in 40 seeds. A 160-seed default census found no generated body outside a companion's orbit, so the change moves only systems that put planets there.
- generation-output: changed
- generation-affects: fill
