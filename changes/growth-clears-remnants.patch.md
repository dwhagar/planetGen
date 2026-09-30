### Fixed
- **Sector growth could place star systems inside a black hole's or
  neutron star's Hill sphere.** `SpaceSector._fine_tune_position` now
  checks every massive neighbor (systems and placed compact remnants),
  not just other systems (TODO item 45).
