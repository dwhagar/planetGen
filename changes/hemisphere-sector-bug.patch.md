### Fixed
- **The 3D Galaxy Map drew a generated neighborhood as a lopsided
  half-sphere.** Each map tile returns at most 250 placed sectors, and past
  that it kept the lowest ids. Neighborhoods are generated shell by shell
  outward from the core, so the lowest ids are the core-facing shells: a
  100 ly sphere of ~2,770 sectors showed only its inner ~1,000, as a bowl
  with a flat face where the cap ran out. A full tile now returns every Nth
  sector by id, so the sample covers the whole neighborhood.
