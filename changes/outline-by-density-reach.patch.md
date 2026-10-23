**Headline:** The plan's outline takes clusters and dwarf galaxies in by the radius where their density reaches one system per sector, and prints how many sectors really expect a system.

### Changed

- The outline took a cluster or satellite in by its tidal radius, so a diffuse dwarf galaxy whose tidal sphere covers the host (the old Sagittarius row, 15.6 kpc) added billions of candidate sectors though its density never reaches one expected system in a sector. A component now reaches as far as its density does (the radius where its density times a sector's volume falls to one system, grown by the sector's half-diagonal, never past the tidal radius), and one that holds under one system in a sector adds nothing. The plan line now prints the candidate sectors (still an upper bound, since the outline keeps the best longitude and the whole run of bands up to a component's) and about how many of them expect at least one system, integrated from the density. For the 2139DD58 galaxy the reaching count is about 1.0 billion. A random start is now redrawn until its sector expects at least one system (it was drawn uniformly over the candidates, and most of them were empty); a start given an explicit density is unchanged.

generation-output: changed

TODO items resolved / Added: GEN.627 / Added: none
