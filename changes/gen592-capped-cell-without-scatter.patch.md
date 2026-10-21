**Headline:** A test now holds that a dense-core cell filled without the star scatter makes the systems its resolution gives.

### Changed
- A cell over the cell star ceiling that is filled with no star scatter drawn makes the number of systems `cell_resolution` gives (GEN.592). The reported case, a capped cell making none, does not reproduce on the current generator: the probe seed's first dense-core cell resolves 3,973 systems and makes 3,989. The test keeps it so.

### Fixed
- The black-hole stream test (GEN.563) scattered a shell of ten million rows twice and held 11 GB, so it was killed in a full run; it now uses a shell of a few hundred thousand rows and takes seconds.
