**Headline:** The node-planet stream test compares what settling leaves alone (masses) and two regenerations, since settling can respace a planet's orbit.

### Changed
- TEST.787: `test_a_nodes_planets_are_drawn_from_a_stream_of_its_own` compared the regenerated planets' distances with the stored ones, which have been through the settle pass; settling spaces the orbits again once the eccentricities are final and, with GEN.1709's wider inner systems, moves one planet of the sample from 71.8 to 122.3 AU. It now compares the masses with the stored planets' and two regenerations' distances with each other.
- Tests only; no generation output changes.
