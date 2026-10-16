**Headline:** Stars, black holes, neutron stars and planetary nebulae in the scatter now come with their companions.

### Added
- **The scatter writes member trees (GEN.262).** The mass pass, the luminosity pass, their backfills and the phenomenon pass give each placed object that has company its stars, brown dwarfs, white dwarfs, neutron stars and black holes in `scatter_members` and `scatter_nodes`, with the row's `member_count` and `system_mass_solar` set (the row's point is the system's barycentre). A star row always draws its tree. A stellar black hole has a companion 1.25% of the time (a star 1.0%, a black hole 0.2%, a neutron star 0.05%), an intermediate-mass black hole a star 1%, a neutron star and a planetary nebula's central white dwarf at their tables. The tree is a function of the row's seed alone, and its nodes' phases are stamped with the orbit epoch (a fixed origin before any orbit update), so two scatters with the same seed store the same members. `planetgen plan` says how many of the stars and phenomena it placed have companions.

### Changed
- **The black-hole and neutron-star rates are rates of rows.** A row is a system keyed by its placed object, so a class's object rate is divided by one plus the same-class companions per row (1.002 for a stellar black hole), which keeps the density of the class. Seeded scatters therefore differ from before by about 0.2% of their black holes.
- **Deleting scatter rows deletes their members.** Clearing a scatter, one star pass, an unfinished band or the unbuilt stars of a cell takes the rows' members and nodes with them.
- **A pre-main-sequence member's `state` is `PMS`**, which fits the 16-character column (the first tree to be written exposed that `pre-main-sequence` did not).

TODO items resolved: GEN.262. Added: none.
