### Fixed
- The planet-offset tests no longer assume the primary star is the system frame's origin: a close pair's origin is its barycenter and a wide pair's is its primary star, so a random system used to fail them about one run in four (TEST.108).
- The Sector Map nebula click test waits until the pointer, with the real mouse on the spot, still reads the nebula before it clicks, so a pick that settles after load no longer fails it (TEST.109).
