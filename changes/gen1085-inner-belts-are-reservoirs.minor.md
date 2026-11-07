**Headline:** A planet can now catch bodies from a belt near its orbit, and a belt body that is caught keeps its belt designation, leaves the belt's mass and shows as a moon.

### Added
- **Inner belts give planets irregular moons (GEN.1085, part 2).** A belt that is not an outer disk, whose nearest edge is within about 3 AU of a giant's orbit (scaled by the giant's Hill radius), is a reservoir for it in the encounter window, with the outer disk's size law and yield scaled by the belt's primordial mass. A body at least as heavy as the belt's lightest stored member is that member (the nearest in mass that no earlier planet of the star caught): it keeps its rank and its designation `<belt ID>.<rank hex>`, its mass and makeup as the belt shows them, and the belt's row stays with the status `moon`. A lighter body of the belt is composed from the belt's stretch and has no rank. The belt's mass falls by what its planets hold. Schema 224 adds `belt_order` and `belt_rank` to `irregular_moons`. Stand-in rows CD2 and CD3 of section 6.193.

### Changed
- A planet with no outer disk to catch from can now still hold bodies of an inner belt in reach, and no binary exchange (a Triton) is drawn where there is no outer disk.

generation-output: changed

TODO items resolved: GEN.1085 (with part 1).
