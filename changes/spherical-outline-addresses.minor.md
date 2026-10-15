### Added
- The galaxy outline on spherical shells (`skeleton.build_shell_extents`, `ShellBounds`): one band range per shell, built by bisection from the density bound; the default galaxy is 3,765 shells and about 24.20 billion candidate cells (the cylinder's 23.92 billion plus 1.2% in shells 928 to 1020, where the bound dips between the equator and the pole). Nothing stores or reads it yet.
- Sector designations on the shell grid (`spherical.provisional_sector_designation`): shell, band and wedge in 12, 12 and 16 bits, the same 40 bits as the birth-sector field of an object ID, so a designation is the first group of every ID born in the sector. `skeleton.outline_object_id_layout` is the plan step's width check.
- Sphere enumeration, the sectors along a segment, and cell description for the shell grid, for the generation paths and pages that move to it.

### Changed
- `object_uid` names its fields shell, band and wedge (the widths and the packing are unchanged).

TODO items resolved: GEN.225, GEN.226 (the designation, ID fields and width check; the sector ID and stored rows switch with DB.32) / Added: none
- The wedge count of a band is Boss's uniform rule (2026-10-15 07:17Z): `6 * m`, `m` a number with no prime factor above 7 nearest (in ratio) to the 4 pc arc, in every band of every shell from 4 out (caps stay 6 at the rim, shells 1 to 3 stay the hand table). Shell 4 now has 206 cells, the default outline 24.20 billion.
