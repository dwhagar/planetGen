**Headline:** The 3D system view draws a trojan region as an arc about its giant's L4 or L5 that follows the giant and librates, not as a full ring.

### Added
- MAP.433: the scene's belt row of a trojan region carries `region` (the giant, the point, the libration period, amplitude and phase, the arc's half width, the mean inclination), and `systemview3d.js` draws it as an arc turned each update to the giant's longitude plus 60 (L4) or minus 60 (L5) degrees plus the libration, counted from Unix time zero as the stored phase is.

### Fixed
- A trojan region is no belt in the count a stored irregular moon or a promoted planet names its belt by (`irregular_moons.belts_of`, the readers in `query`, the sheet and the moon pages, the scene), and it has its own body-number slot under its giant (`trojan:L4`), so a region ahead of a belt no longer shifts that belt's number.
