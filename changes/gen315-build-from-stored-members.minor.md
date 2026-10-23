**Headline:** A star the galaxy scatter placed with one companion star now builds as that pair, on the orbit the scatter stored.

### Changed
- **A pre-placed star is built from the tree its scatter row stored (GEN.315).** A star with one companion star becomes that pair: the second star is built from its stored seed and birth mass, and the pair's separation, eccentricity, plane, phase and period are the stored node's (`BinaryStarProxy` and `WideBinaryPair` take it as `orbit=`), written to `orbit_nodes` as it was stored. A scattered star no longer draws a companion of its own, with or without a stored tree. A tree the pair classes cannot hold (three or more members, a neutron star, black hole or brown dwarf among them) builds the placed star alone and keeps its stored rows for the N-star build; a built tree's `scatter_members` and `scatter_nodes` are deleted when the sector is saved, and by the uncharted-system job.

generation-output: changed
