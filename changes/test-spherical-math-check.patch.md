**Headline:** The math check now tests the shell, band and wedge grid the program uses.

### Changed
- The three math checks that exercised the old cylinder grid (ring cells filling an annulus, the address round trip, the ring count rule) are replaced by five on the spherical grid: shells tile their spherical shell, shells are 4 pc deep around a 2 pc core ball with an odd count of bands from pole to pole, every band holds 6*m wedges (m only 2, 3, 5, 7) and the caps are one cell, every cell's center is found again in its own cell, and seeded random points lie inside the shell, band and wedge their address names. New tests show the gate fails when the wedge rule, the lookup or the cell volume is broken. The run still takes about a second.
