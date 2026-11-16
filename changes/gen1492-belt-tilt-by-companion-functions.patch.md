**Headline:** The functions for a companion star's tilt of a belt are built, not wired yet.

### Changed
- GEN.1492 (second part, slice 0): `physics/belt_forcing.py` gains `planet_nodal_rate`, `forced_normal` and `tilted_orbit`: the axis a belt member's orbit plane precesses about under a tilted companion and the star's planets (a heavy giant holds the belt flat, a bare belt follows the companion), and the member's inclination and node on the cone about it. Registry rows FE3 (now the tilt theory) and FE5 (the 0.5 rad limit), wiring register lines and a design section. Nothing reads them yet, so nothing generated changes; slice 1 (lock item) wires them into the belt's members.

generation-output: none
