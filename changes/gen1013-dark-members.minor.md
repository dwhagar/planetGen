**Headline:** A brown dwarf, neutron star or black hole in a multiple system now stays in its orbit tree, so the orbits, the stability and the gravity are the whole tree's, not a tree with the dark member cut out.

### Changed
- **Dark members stay in the tree (GEN.1013).** The orbit draw already made a brown dwarf, neutron star or black hole a member of the tree with its real mass and radius, and the build then cut it out and closed the orbit above it, so the nodes above had a different mass from the one the stability screens were run with. The member now hangs in the built system as a dark member: its node, slot, kind, mass and radius are kept (`StarSystem.dark_members`), every node keeps the mass of everything below it, the planet windows and the star offsets count it, and it gives no light and has no planets. `multistar.without_members` is removed.
- **The pair is two built stars (GEN.1013).** Where no node holds two built stars (a star beside a sub-pair of a star and a brown dwarf) the pair is the two heaviest built stars on the lowest node above both. A tree with one built star is a single star with the dark members about it.
- **Stored (GEN.1013).** `dark_members` (schema v188, migration 0188) holds one row for each dark member; the clearance verifier reads them as part of the tree.
- generation-output: changed
