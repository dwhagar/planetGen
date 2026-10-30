**Headline:** Listing what is near a place reads each system's name only for the page shown, so a wide search in a dense area is about three times faster.

### Changed

- PERF.217: `planetgen.db.near` lists the systems of a sphere from the covering cell index (id and position) and reads the name and sector of a system only for the systems on the page (`_system_labels`, a read by id), instead of reading the table row of every system in the sphere. At 50 pc around the Sun's density (62,000 systems) the list took 1,170 ms with the name and 360 ms for the ids alone in the profile of PERF.216. The rows, their order, names, parents and counts are the same; a test compares every page of a sphere with the names read directly and checks that the system list asks for no name. No stored value changes.
