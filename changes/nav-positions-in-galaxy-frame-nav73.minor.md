**Headline:** A route over several sectors now measures true distances: a routed system sits where the database stores it.

### Changed

- The route search added a system's sector-local offset to its sector's centre without turning it by the sector's own axes (`frame.sector_orientation`: local +X points out from the galactic axis), so a routed system sat up to a few light-years from its stored `galaxy_x_pc`, `galaxy_y_pc` and `galaxy_z_pc` (2.7 to 4.5 ly measured at a longitude of about 58 degrees; none only on the +X axis). `corridor._systems_of` and `query._load_nav_endpoint` now use `frame.local_to_galaxy_pc`, and a route between two systems of one sector that is searched in the galaxy frame comes back in that sector's own axes with `frame.galaxy_to_local_pc` (it was a plain shift). The hop lengths and courses of every cross-sector route change to the true ones; stored values are unchanged.
- The graph round a route (NAV.65) no longer needs its workaround: it reads the route's own positions, which are the index's, and the hops keep the route's lengths because those are now the true ones.
- Tests: a sector at a longitude of 58 degrees where every routed system equals its stored galaxy position, the hop lengths of a cross-sector route equal the distances between the stored positions, and a route inside one sector comes back in the sector's axes; the NAV.65 graph tests now run at that longitude too.

generation-output: unchanged
