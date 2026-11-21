**Headline:** The trojan page text reads `trojans.pair_text`, and the unread Python Lagrange geometry is removed (the gravity layer has its own).

### Changed
- GEN.1820: a trojan region's page row `Lagrange point` and a trojan planet's sheet row `Shares an orbit with` are worded by `trojans.pair_text` (`L4 of Sun and Jupiter`), where they were the bare point or a hand-built sentence. A region whose giant or star is not known keeps the bare point.

### Removed
- `physics/lagrange_geometry.py` and its test: nothing read it. The System Map's gravity layer (MAP.169) computes the five points in `systemgravity.js` (`lagrangePoints`, tested in `systemgravity.test.mjs` against the closed forms) and the map icon is labelled by `trojans.point_label`, so the Python copy was redundant. The design notes and the wiring register say so.

generation-output: unchanged
