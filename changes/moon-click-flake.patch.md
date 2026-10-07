### Fixed

- The System Map browser test no longer fails now and then: it clicks each body on its dot instead of the middle of its marker, which can be empty map when the body's label sits to one side, and it measures from the star to the planet with moons, so a system with only one planet still works.
- The Galaxy Map drill-down browser test no longer clicks the block it is already in when, under load, the map's tooltip still names it.
